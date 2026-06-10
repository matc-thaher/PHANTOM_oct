function [c, mask] = Diemer19_general(M, z, cosmo, params, profile_name)
% DIEMER19_GENERAL  Diemer & Joyce (2019) / Ishiyama+21 concentration engine.
%
%   [c, mask] = Diemer19_general(M, z, cosmo, params)
%   [c, mask] = Diemer19_general(M, z, cosmo, params, profile_name)
%
%   profile_name : (optional) string passed to build_lhs_table.
%                  Default: 'nfw'
%
%   Compatible with MATLAB R2020a+ and GNU Octave.

    if nargin < 5 || isempty(profile_name)
        profile_name = 'nfw';
    end

    % --- build (or retrieve cached) inversion table ---
    [conc_interp, lhs_min_interp, lhs_max_interp] = build_lhs_table(profile_name);

    % --- cosmological quantities ---
    sigma     = cosmo.sigmaM(M, z);
    delta_c   = 1.686;
    nu        = delta_c ./ sigma;
    n_eff     = cosmo.neff(M, z, params.kappa);
    alpha_eff = cosmo.alphaEff(z);

    % --- flatten to column vectors ---
    nu        = nu(:);
    n_eff     = n_eff(:);
    alpha_eff = alpha_eff .* ones(size(nu));

    % --- Diemer19 / Ishiyama21 parameter combinations ---
    A_n     = params.a0 .* (1 + params.a1 .* (n_eff + 3));
    B_n     = params.b0 .* (1 + params.b1 .* (n_eff + 3));
    C_alpha = 1 - params.c_alpha .* (1 - alpha_eff);

    rhs = log10(A_n ./ nu .* (1 + nu.^2 ./ B_n));

    % --- validity mask ---
    lhs_min_n = lhs_min_interp(n_eff);
    lhs_max_n = lhs_max_interp(n_eff);
    mask      = (rhs >= lhs_min_n) & (rhs <= lhs_max_n);

    % --- invert G(c,n) -> log10(c) ---
    c      = nan(size(nu));
    log10c = conc_interp(rhs, n_eff);
    c(mask)  = 10.^log10c(mask) .* C_alpha(mask);
    c(~mask) = NaN;

    c    = reshape(c, size(M));
    mask = reshape(mask, size(M));
end