function b = halo_bias_seljak04(sigma, delta_c, cosmo)
% Seljak & Warren (2004, MNRAS 355, 129) halo bias.
% Numerically calibrated fit; Eq. 5 (base model).
% Optionally adds Eq. 6 cosmology-dependent correction if cosmo is supplied.
%
% b(x) = 0.53 + 0.39 x^0.45 + 0.13/(40x + 1) + 5e-4 x^1.5
%         + (log10(x)/10) * [0.4(Om_0 - 0.3 + n_s - 1.0)]   [Eq.6, optional]
% where x = nu = delta_c / sigma.
%
% Inputs:
%   sigma    : rms variance sigma(M,z), vector
%   delta_c  : collapse threshold (default EdS value)
%   cosmo    : cosmology struct (optional); uses cosmo.Om0 and cosmo.ns
%              for the Eq.6 cosmology correction. Omit for Eq.5 only.

    if nargin < 2 || isempty(delta_c)
        delta_c = collapse_overdensity();
    end

    nu = delta_c ./ sigma;

    % Equation 5 — base fit
    b0 = 0.53 + (0.39 .* nu.^0.45) ...
             + (0.13 ./ ((40 .* nu) + 1)) ...
             + (5e-4 .* nu.^1.5);

    % --- Eq. 6 cosmology correction ---
    if nargin < 3 || isempty(cosmo)
        b = b0;
        return
    end

    % Extract cosmological parameters
    % m = Omega_m (Seljak & Warren use present-day value; z-dependence
    % can be added via cosmo.Omega_m_z(z) if needed)
    m        = cosmo.Omega_m;
    ns       = cosmo.ns;
    sigma8   = cosmo.sigma8;
    h        = cosmo.h;

    % alpha_s: running of the spectral index d(ns)/d(ln k)
    % Not stored in cosmo struct — defaults to 0 (scale-free spectrum)
    if isfield(cosmo, 'alpha_s')
        alpha_s = cosmo.alpha_s;
    else
        alpha_s = 0;
    end

    % Eq. 6 correction — three bracketed groups with coefficients 0.4, 0.3, 0.8
    correction = log10(nu) .* ( 0.4 * (m    - 0.3  + ns - 1  ) ...
                              + 0.3 * (sigma8 - 0.9 + h  - 0.7) ...
                              + 0.8 *  alpha_s );

    b = b0 + correction;
end
