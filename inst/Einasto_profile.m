function [rho, rhos, rs, fc] = Einasto_profile(r, M, c, z, cosmo, Delta)
% Einasto_profile  Einasto (1965) density profile with alpha_e from Gao+2008
%
%   [rho, rhos, rs] = Einasto_profile(r, M, c, z, cosmo, Delta)
%
%   The Einasto shape parameter alpha_e is computed internally from the
%   dimensionless peak-height nu(M,z) using the empirical fitting formula
%   of Gao et al. (2008, MNRAS 387, 536), eq. (5):
%
%       alpha_e = 0.155 + 0.0095 * nu^2
%
%   where nu = delta_crit(z) / sigma(M,z) is the peak-height parameter,
%   i.e. the ratio of the linear collapse threshold to the rms density
%   fluctuation within spheres of mass M at redshift z.
%   nu < 1 : low-mass / common haloes    -> alpha_e ~ 0.16
%   nu ~ 3 : rare, massive clusters      -> alpha_e ~ 0.24-0.30
%
%   NOTE: alpha_e is capped at 0.3 following Benson (2012), since Gao+2008
%   did not probe nu > ~3.5 and extrapolation beyond that is not justified.
%
%   INPUTS
%   r       : radii [Mpc/h], scalar or vector
%   M       : halo mass [Msun/h]
%   c       : concentration  R_Delta / r_s
%   z       : redshift
%   cosmo   : cosmology struct with fields:
%               .rhocrit0  critical density at z=0 [Msun/h / (Mpc/h)^3]
%               .E(z)      dimensionless Hubble parameter E(z) = H(z)/H0
%               .nu(M,z)   peak-height nu(M,z) = delta_crit(z)/sigma(M,z)
%   Delta   : overdensity w.r.t. critical density (e.g. 200)
%
%   OUTPUTS
%   rho     : density profile at r [Msun/h / (Mpc/h)^3]
%   rhos    : density at scale radius rs [same units]
%   rs      : scale radius [Mpc/h]
%   fc     : concentration-dependent factor
%
%
%   Profile form:
%     rho(r) = rhos * exp( -(2/alpha_e) * [ (r/rs)^alpha_e - 1 ] )
%
%   References:
%     Einasto 1965
%     Merritt et al. 2006, AJ, 132, 2685
%     Gao et al. 2008, MNRAS, 387, 536  [alpha_e(nu) fitting formula, eq. 5]

    % ---- Convert M to Mvir for Gao+2008 nu_vir (Colossus convention) --------
    % Colossus always uses nu_vir for alpha_e regardless of input mdef.
    % We convert M_Delta -> Mvir using change_mass_definition with NFW profile.
    [Mvir, ~, ~] = change_mass_definition(M, c, z, Delta, 'c', 0, 'vir', cosmo);

    % ---- Peak-height nu using Mvir (matches Colossus convention) ----
    nu      = cosmo.nu(Mvir, z);
    
    % ---- Gao+2008 shape parameter (eq. 5) ---------------------------
    % alpha_e = 0.155 + 0.0095 * nu^2
    % Capped at 0.3: Gao+2008 did not constrain nu >> 3.5, and
    % extrapolating to higher peak heights is not justified.
    alpha_e = min(0.155 + 0.0095 .* nu.^2, 0.3);

    % ---- Critical density at redshift z -----------------------------
    rho_c   = cosmo.rho_crit0 .* cosmo.E(z).^2;

    % ---- Virial / overdensity radius and scale radius ---------------
    R_Delta = (3 .* M ./ (4 .* pi .* Delta .* rho_c)).^(1/3);
    rs      = R_Delta ./ c;

    % ---- Normalise rhos so that M(<R_Delta) = M ---------------------
    % Enclosed mass of Einasto profile:
    %   M(<r) = 4*pi*rhos*rs^3 * exp(2/alpha_e) / alpha_e
    %           * (alpha_e/2)^(3/alpha_e) * gamma(3/alpha_e, (2/alpha_e)*x^alpha_e)
    % where gamma(a,x) is the lower incomplete gamma function.
    % At r = R_Delta, x = c:
    hn      = 2 ./ alpha_e;                         % = 2/alpha_e
    x_edge  = hn .* c.^alpha_e;                     % argument at R_Delta
    Ig      = gammainc(x_edge, 3./alpha_e) ...
              .* gamma(3./alpha_e);                  % lower incomplete gamma * Gamma
    fc      = Ig;
    rhos    = M ./ (4.*pi .* rs.^3 ...
              .* exp(hn) ./ alpha_e ...
              .* hn.^(-3./alpha_e) .* Ig);

    % ---- Einasto density profile ------------------------------------
    % rho(r) = rhos * exp( -(2/alpha_e) * [ (r/rs)^alpha_e - 1 ] )
    x   = r(:) ./ rs;
    rho = rhos .* exp(-hn .* (x.^alpha_e - 1));

end