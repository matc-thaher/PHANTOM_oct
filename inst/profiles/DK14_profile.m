function rho = DK14_profile(r, M, c, z, cosmo, Delta, selected_by, Gamma)
% DK14_profile  Diemer & Kravtsov (2014) density profile
%
%   rho = DK14_profile(r, M, c, z, cosmo, Delta, selected_by, Gamma)
%
%   The DK14 profile multiplies an Einasto inner profile by a truncation
%   (splashback) function and adds a power-law outer (infalling) term:
%
%     rho(r) = rho_inner(r) * f_trans(r) + rho_outer(r)
%
%   where
%     rho_inner(r) : Einasto profile via Einasto_profile.m
%                    (alpha_e from Gao+2008 inside that function)
%     f_trans(r)   = [ 1 + (r/rt)^beta ]^(-gamma_t/beta)
%     rho_outer(r) = rho_m * b_e * (r/r_ref)^(-s_e)
%
%   PARAMETER DEFAULTS  (following Colossus / DK14 deriveParameters logic)
%   -----------------------------------------------------------------------
%   selected_by = 'M'  (mass-selected sample)
%       beta    = 4
%       gamma_t = 8
%       rt      = R200m * (1.9 - 0.18*nu200m)          [DK14 eq. 6]
%
%   selected_by = 'Gamma'  (mass + accretion-rate selected sample)
%       beta    = 6
%       gamma_t = 4
%       rt      = R200m * 0.54 * (1 + 0.53*exp(-Gamma)) [DK14 Gamma form]
%       Both z and Gamma must be provided.
%
%   Outer term defaults:
%       s_e   = 1.5   (outer power-law slope)
%       rho_m = cosmo.rho_m0 * (1+z)^3  (mean matter density at z)
%
%   INPUTS
%   r           : radii [Mpc/h], scalar or vector
%   M           : halo mass [Msun/h]  (as M_200m)
%   c           : concentration c_200m = R_200m / r_s
%   z           : redshift
%   cosmo       : cosmology struct with fields:
%                   .rhocrit0   critical density at z=0 [Msun/h/(Mpc/h)^3]
%                   .E(z)       E(z) = H(z)/H0 function handle
%                   .nu(M,z)    peak-height nu(M,z) function handle
%                   .Omega_m    matter density parameter (for rho_m0)
%                   .rho_m0     mean matter density at z=0 [Msun/h/(Mpc/h)^3]
%   Delta       : overdensity w.r.t. critical density (typically 200)
%   selected_by : (optional) 'M' [default] or 'Gamma'
%                   'M'     -> sample selected by mass only
%                   'Gamma' -> sample selected by mass AND accretion rate
%   Gamma       : (optional) mass accretion rate as in DK14
%                   Required (and used) only when selected_by = 'Gamma'
%
%   OUTPUT
%   rho         : total DK14 density at each r [Msun/h / (Mpc/h)^3]
%
%   DEPENDENCIES
%   Einasto_profile.m  (must reside in the same folder or on the MATLAB path)
%
%   REFERENCES
%   Diemer & Kravtsov 2014, ApJ 789, 1                    (DK14 profile)
%   Diemer 2022, ApJ 925, 182                             (updated form)
%   Gao et al. 2008, MNRAS 387, 536                       (alpha_e via nu)
%   Colossus deriveParameters() documentation
%   https://bdiemer.bitbucket.io/colossus/halo_profile_dk14.html

    % ---- Input defaults -------------------------------------------------
    if nargin < 7 || isempty(selected_by)
        selected_by = 'M';
    end
    if nargin < 8
        Gamma = [];
    end

    % If selected_by = 'Gamma' but no Gamma value supplied,
    % fall back silently to mass-selected defaults instead of erroring.
    if strcmp(selected_by, 'Gamma') && isempty(Gamma)
        warning('DK14_profile: selected_by = ''Gamma'' but Gamma is empty. Falling back to mass-selected defaults (beta=4, gamma_t=8).');
        selected_by = 'M';
    end

    r = r(:);

    % ---- Mean matter density at z ---------------------------------------
    rho_m_z  = cosmo.rhom(z);          % mean matter density at z

    % ---- R_200m and nu_200m  --------------------------------------------
    % Needed to calibrate rt (DK14 eq. 6 was calibrated for nu_200m).
    % R_200m is computed as the overdensity radius at Delta = 200 w.r.t.
    % critical, consistent with M being M_200m.
    % rho_c    = cosmo.rhocrit(z); %cosmo.rho_crit0 .* cosmo.E(z).^2;
    % R200m    = (3 .* M ./ (4 .* pi .* Delta .* rho_c)).^(1/3);
    % R200m    = (3 .* M ./ (4 .* pi .* 200 .* rho_m_z)).^(1/3);   % true R200m

    % M200m approximation for nu calibration (DK14 eq.6 uses nu_200m)
    % M200m_approx = M .* (rho_m_z ./ rho_c);
    % Exact M200c -> M200m conversion via NFW profile (no Colossus needed)
    [M200m, R200m, ~] = change_mass_definition(M, c, z, Delta, 'c', 200, 'm', cosmo);
    nu200m            = cosmo.nu(M200m, z);    % peak height at M_200m

    % ---- Derive beta, gamma_t, rt  (Colossus deriveParameters logic) ----
    %
    % Case 1: selected_by = 'M'
    %   (beta, gamma_t) = (4, 8)                      [DK14 Table 1]
    %   rt = R200m * (1.9 - 0.18 * nu200m)            [DK14 eq. 6]
    %
    % Case 2: selected_by = 'Gamma'
    %   (beta, gamma_t) = (6, 4)                      [DK14 Table 1]
    %   rt = R200m * 0.54 * (1 + 0.53 * exp(-Gamma))  [DK14 Gamma-based rt]
    %
    % In both cases rt is clipped to a physically reasonable lower bound
    % to avoid rt collapsing to zero for very high-nu halos.
    if strcmp(selected_by, 'M')
        beta    = 4;
        gamma_t = 8;
        rt      = R200m .* (1.9 - 0.18 .* nu200m);   % DK14 eq. 6
    else   % 'Gamma'
        beta    = 6;
        gamma_t = 4;
        rt      = R200m .* 0.54 .* (1 + 0.53 .* exp(-Gamma));
    end

    % Numerical safety: rt must be positive
    rt = max(rt, 0.01 .* R200m);

    % ---- Outer power-law slope default ----------------------------------
    s_e = 1.5;

    % ---- Inner Einasto profile  (delegates to Einasto_profile.m) --------
    % alpha_e is computed inside via Gao+2008: alpha_e = 0.155+0.0095*nu^2
    [rho_inner, ~, ~] = Einasto_profile(r, M, c, z, cosmo, Delta);
    % Inside DK14_profile.m, replace the Einasto call with:
    % [rho_inner, ~, ~] = Einasto_profile(r, M200m, c_new, z, cosmo, 200, '200m');

    % ---- Truncation (splashback) function --------------------------------
    % f_trans(r) = [ 1 + (r/rt)^beta ]^(-gamma_t/beta)
    %   r << rt : f_trans -> 1       (inner Einasto unaffected)
    %   r ~  rt : f_trans ~ 0.5      (splashback transition)
    %   r >> rt : f_trans -> 0       (sharp suppression)
    f_trans = (1 + (r ./ rt).^beta).^(-gamma_t ./ beta);

    % % ---- Outer power-law (infalling) term --------------------------------
    % % rho_outer = rho_m * (r)^(-s_e)
    % % b_e * r_ref^s_e are absorbed into rho_m by convention (b_e ~ 1,
    % % r_ref = 1 Mpc/h so the factor is unity in these units).
    % rho_outer = rho_m .* (r).^(-s_e);

    % ---- Total DK14 profile ---------------------------------------------
    rho = rho_inner .* f_trans; % + rho_outer;

end