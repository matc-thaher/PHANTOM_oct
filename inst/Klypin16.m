function [c, valid] = Klypin16(M, z, cosmo, cosmo_name, mdef, formula)
% Klypin16_concentration  Klypin et al. (2016) concentration model
%
%   [c, valid] = Klypin16(M, z, cosmo, cosmo_name, mdef, formula)
%
%   Implements BOTH fitting formulae from Appendix A of Klypin+2016.
%   The user selects which formula to use via the 'formula' argument.
%
%   FORMULA 'cM'  — eq. (A1), mass-based:
%     c(M,z) = C0(z) * (M/1e12)^(-gamma(z)) * [1 + (M/M0(z))^0.4]
%     Parameters C0, gamma, M0 are linearly interpolated from Tables A1-A4.
%
%   FORMULA 'cnu' — eq. (A5), peak-height-based:
%     c(sigma) = b0(z) * [1 + 7.37*(sigma/a0(z))^0.75]
%                       * [1 + 0.14*(sigma/a0(z))^(-2)]
%     where sigma = sigma(M,z) from the linear power spectrum.
%     Parameters a0, b0 interpolated from Tables A5-A8.
%
%   INPUTS
%   M          : halo mass [Msun/h], scalar or vector
%   z          : redshift (scalar)
%   cosmo      : cosmology struct; must have field:
%                  cosmo.sigmaM(M, z)  — only needed for formula='cnu'
%   cosmo_name : 'planck13' | 'bolshoi'
%   mdef       : '200c'     | 'vir'
%   formula    : 'cM'       | 'cnu'
%
%   OUTPUTS
%   c          : concentration, same shape as M
%   valid      : logical mask, same shape as M
%                  true where inputs are within calibrated range:
%                  M > 1e10 Msun/h  AND  z <= max(z_bins)
%
%   Reference: Klypin, Yepes, Gottlober, Prada & Hess 2016, MNRAS 457, 4340

% ---- load parameters ---------------------------------------------------
P = Klypin16_Table(formula, cosmo_name, mdef);

% ---- validity mask -----------------------------------------------------
valid = (M(:) > 1e10) & (z <= P.z_bins(end));

% ---- interpolate parameters to redshift z -----------------------------
switch lower(formula)

    % ====================================================================
    case 'cm'
    % eq. (A1):  c = C0 * (M/1e12)^(-gamma) * [1 + (M/M0)^0.4]
    % ====================================================================
        C0_z    = interp1(P.z_bins, P.C0,    z, 'linear', 'extrap');
        gamma_z = interp1(P.z_bins, P.gamma, z, 'linear', 'extrap');
        M0_z    = interp1(P.z_bins, P.M0,    z, 'linear', 'extrap');

        c = C0_z .* (M(:) ./ 1e12).^(-gamma_z) .* (1 + (M(:) ./ M0_z).^0.4);

    % ====================================================================
    case 'cnu'
    % eq. (A5):  c = b0 * [1 + 7.37*(sigma/a0)^0.75] * [1 + 0.14*(sigma/a0)^-2]
    % ====================================================================
        a0_z = interp1(P.z_bins, P.a0, z, 'linear', 'extrap');
        b0_z = interp1(P.z_bins, P.b0, z, 'linear', 'extrap');

        sigma   = cosmo.sigmaM(M(:), z);
        x       = sigma ./ a0_z;                    % sigma / a0

        c = b0_z .* (1 + 7.37 .* x.^0.75) .* (1 + 0.14 .* x.^(-2));

    otherwise
        error('Klypin16_concentration: unknown formula "%s".', formula);
end

% restore input shape
c     = reshape(c,     size(M));
valid = reshape(valid, size(M));
end