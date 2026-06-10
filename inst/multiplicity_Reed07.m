% function f = multiplicity_Reed07(sigma, z, cosmo, delta_c)
% % Reed et al. (2007), MNRAS 374, 2, Eq. 11  ← recommended Reed model
% % Incorporates the effective spectral index n_eff via cosmo.neff.
% % This is the more physically motivated of the two Reed models.
% %
% % INPUT:
% %   sigma   : sigma(M,z)
% %   z       : redshift
% %   cosmo   : PHANTOM cosmo struct (needed for neff)
% %   delta_c : (optional) collapse overdensity
% 
%     if nargin < 4 || isempty(delta_c)
%         delta_c = collapse_overdensity('corrections', true, 'z', z);
%     end
% 
% 
%     % Effective spectral index at kappa = 1 (Lagrangian scale)
%     M       = cosmo.R_of_M([]);     % placeholder; neff needs M, not sigma
%     % Invert sigma -> M via cosmo to get neff: use nu and sigmaM inversion
%     % Practical approach: compute neff at the Lagrangian radius of each sigma
%     % using dlnsigma_dlnM and the chain rule neff = -3 - 6*dlns_dlnM
%     % (d ln sigma / d ln R = 3 * d ln sigma / d ln M, so neff = -3 - 2*(3*dlns/dlnM))
%     % Here we use the approximate formula from Reed07 Eq. 14 as fallback
%     % if cosmo is not passed, but prefer numerical neff via cosmo.neff.
%     if nargin >= 3 && ~isempty(cosmo)
%         % Recover approximate M from sigma at this z using sigmaM inversion
%         % via a simple root-find on a precomputed grid
%         M_grid    = logspace(6, 17, 500);
%         sig_grid  = cosmo.sigmaM(M_grid, z);
%         M_of_sig  = interp1(log(sig_grid(end:-1:1)), log(M_grid(end:-1:1)), ...
%                             log(sigma(:).'), 'linear', 'extrap');
%         M_arr     = exp(M_of_sig);
%         n_eff     = cosmo.neff(M_arr, z, 1);   % kappa = 1
%     else
%         % Approximate n_eff from Reed07 Eq. 14
%         mz = 0.55 - 0.32 .* (1 - 1./(1+z)).^5;
%         rz = -1.74 - 0.8 .* abs(log(1./(1+z))).^0.8;
%         n_eff = mz .* log(1./sigma) + rz;
%     end
% 
%     A   = 0.3222;
%     ca  = 0.764;
%     c   = 1.08;
%     a   = ca / c;     % = 0.707
%     p   = 0.3;
% 
%     lsi = log(1./sigma);     % ln(sigma^-1)
%     nu  = delta_c ./ sigma;
% 
%     G1 = exp(-(lsi - 0.4).^2 ./ (2.*0.6^2));
%     G2 = exp(-(lsi - 0.75).^2 ./ (2.*0.2^2));
% 
%     f = A .* sqrt(2.*a./pi) ...
%          .* (1 + (a.*nu.^2).^-p + 0.6.*G1 + 0.4.*G2) ...
%          .* nu .* exp(-(ca.*nu.^2./2) - (0.03.*nu.^0.6 ./ (n_eff + 3).^2));
% 
% end

function f = multiplicity_Reed07(sigma, z, cosmo, delta_c)
% Reed et al. (2007), MNRAS 374, 2, Eq. 11

    sigma_in = sigma;
    sigma    = sigma(:);   % force column vector

    if nargin < 4 || isempty(delta_c)
        delta_c = collapse_overdensity('corrections', true, 'z', z, 'cosmo', cosmo);
    end

    if nargin >= 3 && ~isempty(cosmo)
        M_grid   = logspace(6, 17, 500);
        sig_grid = cosmo.sigmaM(M_grid, z);

        % Make both interpolation arrays monotonic and column-like
        xgrid = log(sig_grid(end:-1:1));
        ygrid = log(M_grid(end:-1:1));
        xq    = log(sigma);

        M_of_sig = interp1(xgrid, ygrid, xq, 'linear', 'extrap');
        M_arr    = exp(M_of_sig);

        n_eff = cosmo.neff(M_arr, z, 1);
        n_eff = n_eff(:);   % force column
    else
        mz = 0.55 - 0.32 .* (1 - 1./(1+z)).^5;
        rz = -1.74 - 0.8 .* abs(log(1./(1+z))).^0.8;
        n_eff = mz .* log(1./sigma) + rz;
        n_eff = n_eff(:);
    end

    A   = 0.3222;
    ca  = 0.764;
    c   = 1.08;
    a   = ca / c;
    p   = 0.3;

    lsi = log(1./sigma);
    nu  = delta_c ./ sigma;

    G1 = exp(-(lsi - 0.4).^2 ./ (2 .* 0.6^2));
    G2 = exp(-(lsi - 0.75).^2 ./ (2 .* 0.2^2));

    f = A .* sqrt(2 .* a ./ pi) ...
        .* (1 + (a .* nu.^2).^-p + 0.6 .* G1 + 0.4 .* G2) ...
        .* nu ...
        .* exp(-(ca .* nu.^2 ./ 2) - (0.03 .* nu.^0.6 ./ (n_eff + 3).^2));

    f = reshape(f, size(sigma_in));   % restore original shape
end