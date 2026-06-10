function b = halo_bias_jing98(sigma, delta_c, cosmo)
% Jing (1998, ApJ 503, L9) halo bias.
% Calibrated from scale-free simulations; n_eff is the local slope
% of the linear power spectrum at the non-linear scale k_nl.
%
% b(nu) = (0.5/nu^4 + 1)^(0.06 - 0.02*n_eff) * (1 + (nu^2 - 1)/delta_c)
%
% Inputs:
%   sigma    : rms variance sigma(M,z), vector
%   delta_c  : collapse threshold (default from collapse_overdensity)
%   cosmo    : cosmology struct with field cosmo.n_eff (local spectral slope)
%              If cosmo.n_eff is absent, n_eff = -2 is used as a fallback.

    if nargin < 3 || isempty(delta_c)
        delta_c = collapse_overdensity();
    end
    if nargin < 2 || ~isfield(cosmo, 'n_eff') || isempty(cosmo.n_eff)
        n_eff = -2.0;   % rough CDM value near the non-linear scale
    else
        n_eff = cosmo.n_eff;
    end

    nu  = delta_c ./ sigma;
    nu2 = nu.^2;

    b = (0.5 ./ nu2.^2 + 1).^(0.06 - (0.02 * n_eff)) .* ...
        (1 + (nu2 - 1) / delta_c);
end
