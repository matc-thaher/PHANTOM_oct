function [M_new, R_new, c_new] = change_mass_definition(M, c, z, Delta_in, type_in, Delta_out, type_out, cosmo)
% change_mass_definition   Convert spherical overdensity mass definitions
%
%   [M_new, R_new, c_new] = change_mass_definition(M, c, z, Delta_in, type_in, Delta_out, type_out, cosmo)
%
%   Assumes a fixed NFW profile in physical coordinates and finds the new
%   spherical overdensity radius where the mean enclosed density equals the
%   new threshold. Replicates colossus.halo.mass_defs.changeMassDefinition
%   with profile='nfw'.
%
%   INPUTS
%   M        : halo mass [Msun/h], scalar or vector
%   c        : concentration R_in / r_s, scalar or vector
%   z        : redshift (scalar)
%   Delta_in : overdensity number for input  mdef  (e.g. 200)
%   type_in  : reference density for input   mdef  'c' (critical) or 'm' (mean)
%   Delta_out: overdensity number for output mdef  (e.g. 200)
%   type_out : reference density for output  mdef  'c' or 'm'
%              use 'vir' for type_in/type_out to get Bryan & Norman virial
%   cosmo    : PHANTOM cosmology struct with fields:
%                .rho_crit0   critical density at z=0 [Msun/h/(Mpc/h)^3]
%                .E(z)        dimensionless Hubble parameter
%                .Omega_m     matter density parameter
%                .Omega_m_z   Omega_m(z) function handle (for vir)
%
%   OUTPUTS
%   M_new    : new mass   [Msun/h]
%   R_new    : new radius [Mpc/h]
%   c_new    : new concentration R_new / r_s
%
%   USAGE EXAMPLES
%   % M200c -> M200m
%   [M200m, R200m, c200m] = change_mass_definition(M200c, c200c, z, 200,'c', 200,'m', cosmo);
%
%   % M200c -> Mvir
%   [Mvir, Rvir, cvir] = change_mass_definition(M200c, c200c, z, 200,'c', 0,'vir', cosmo);
%
%   REFERENCE
%   Colossus: Diemer 2018, ApJS 239, 35
%             mass_defs.changeMassDefinition (profile='nfw')

    % ---- Density thresholds at redshift z --------------------------------
    rho_c_z = cosmo.rho_crit0 .* cosmo.E(z).^2;
    rho_m_z = cosmo.Omega_m   .* rho_c_z;          % rho_mean(z)

    rho_in  = density_threshold(Delta_in,  type_in,  rho_c_z, rho_m_z, cosmo, z);
    rho_out = density_threshold(Delta_out, type_out, rho_c_z, rho_m_z, cosmo, z);

    % ---- NFW native parameters from input (M, c, mdef) ------------------
    R_in  = (3 .* M ./ (4 .* pi .* rho_in)).^(1/3);   % R_Delta_in [Mpc/h]
    rs    = R_in ./ c;                                  % scale radius [Mpc/h]

    % rhos from NFW mass normalisation:
    %   M = 4*pi*rhos*rs^3 * f(c),  f(c) = ln(1+c) - c/(1+c)
    fc    = log(1 + c) - c ./ (1 + c);
    rhos  = M ./ (4 .* pi .* rs.^3 .* fc);

    % ---- Find new concentration c_new by solving NFW xDelta --------------
    % Mean density within r equals rho_out:
    %   3*rhos/x^3 * f(x) = rho_out   where x = r/rs
    %
    % Rearranged: g(x) = f(x)/x^3 - rho_out/(3*rhos) = 0
    % f(x) = ln(1+x) - x/(1+x)
    % g is monotonically decreasing -> safe for bisection

    n      = numel(M);
    c_new  = zeros(size(M));
    R_new  = zeros(size(M));

    target = rho_out ./ (3 .* rhos);   % scalar or vector

    for i = 1:n
        c_new(i) = bisect_nfw(@(x) nfw_mean_density_ratio(x) - target(i), 1e-4, 1e3);
    end

    R_new = rs .* c_new;
    M_new = (4/3) .* pi .* R_new.^3 .* rho_out;
end

% ---- Local helpers -------------------------------------------------------

function val = nfw_mean_density_ratio(x)
    % f(x)/x^3  where f(x) = ln(1+x) - x/(1+x)
    val = (log(1 + x) - x ./ (1 + x)) ./ x.^3;
end

function rho = density_threshold(Delta, type, rho_c_z, rho_m_z, cosmo, z)
    if strcmpi(type, 'c')
        rho = Delta .* rho_c_z;
    elseif strcmpi(type, 'm')
        rho = Delta .* rho_m_z;
    elseif strcmpi(type, 'vir')
        % Bryan & Norman 1998 virial overdensity w.r.t. critical
        Om_z      = cosmo.Omega_m_z(z);
        Delta_vir = 18*pi^2 + 82*(Om_z - 1) - 39*(Om_z - 1)^2;
        rho       = Delta_vir .* rho_c_z;
    else
        error('change_mass_definition: unknown type "%s". Use ''c'', ''m'', or ''vir''.', type);
    end
end

function x = bisect_nfw(f, a, b)
    % Simple bisection for monotone f on [a,b]
    tol = 1e-8;
    fa  = f(a);
    for iter = 1:100
        mid = 0.5*(a + b);
        fm  = f(mid);
        if abs(fm) < tol || (b - a) < tol*mid
            x = mid;  return;
        end
        if sign(fm) == sign(fa)
            a = mid;  fa = fm;
        else
            b = mid;
        end
    end
    x = 0.5*(a + b);
end
