function c = Diemer19_zero_general(M, z, cosmo, params, varargin)
% Diemer19_zero_general  Direct root-finding version of the Ishiyama et al. (2021) concentration model
%
%   c = Diemer19_zero_general(M, z, cosmo, mode)
%   c = Diemer19_zero_general(M, z, cosmo, mode, profile_name)
%
%   Implements the Ishiyama et al. (2021) concentration-mass relation using
%   direct root finding with fzero, while allowing the left-hand-side G(c)
%   inversion to depend on the profile through profile_mu(c, profile_name).
%
%   This is intended to mirror the same profile-dependent definition used by
%   build_lhs_table(profile_name) / Diemer19_general, but solved per halo
%   directly instead of through a precomputed interpolation table.
%
%   INPUTS
%   M            : halo mass [Msun/h], scalar or vector
%   z            : redshift (scalar)
%   cosmo        : cosmology struct built by cosmology() + attach_linear_components()
%   mode         : '200c_all' | '200c_relaxed' | 'vir_all' | 'vir_relaxed' |
%                  '500_all' | '500_relaxed'
%   profile_name : (optional) string, profile used in profile_mu().
%                  Default: 'nfw'
%
%   OUTPUT
%   c            : concentration, same shape as M
%
%   Notes
%   -----
%   - If profile_name = 'nfw', the result should closely track the table-based
%     Ishiyama21 implementation, but may not be bitwise identical because
%     this routine uses direct fzero solves rather than interpolation.
%   - For non-NFW profiles, the inversion is done using mu(c) from
%     profile_mu(c, profile_name), matching the definition used by
%     build_lhs_table(profile_name).
%
%   Reference: Ishiyama et al. 2021, MNRAS 506, 4210
%              Diemer & Joyce 2019, ApJ 871, 168

    if nargin < 5 || isempty(varargin{1})
        profile_name = 'nfw';
    else
        profile_name = varargin{1};
    end

    % P      = Ishiyama21_Table(mode);
    % kappa  = P.kappa;
    % a0     = P.a0;  a1 = P.a1;
    % b0     = P.b0;  b1 = P.b1;
    % cAlpha = P.cAlpha;

    sigma   = cosmo.sigmaM(M, z);
    delta_c = 1.686;
    nu      = delta_c ./ sigma;

    neff  = cosmo.neff(M, z, params.kappa);
    alpha = cosmo.alphaEff(z);

    A_eff = params.a0 .* (1 + params.a1 .* (neff + 3));
    B_eff = params.b0 .* (1 + params.b1 .* (neff + 3));
    C_eff = 1  - params.c_alpha .* (1 - alpha);

    rhs = log10(A_eff ./ nu .* (1 + nu.^2 ./ B_eff));

    rhs      = rhs(:);
    neff_col = neff(:);
    Ceff_col = C_eff(:);

    c_unnorm = nan(size(rhs));

    for i = 1:numel(rhs)
        rhs_i = rhs(i);
        n_i   = neff_col(i);
        expo  = (5 + n_i) / 6;

        obj = @(cv) lhs_profile(cv, expo, profile_name) - rhs_i;

        c_lo = 1e-1;
        c_hi = 1e3;
        f_lo = obj(c_lo);
        f_hi = obj(c_hi);
        expanded = false;

        for k = 1:30
            if isfinite(f_lo) && isfinite(f_hi) && (f_lo * f_hi < 0)
                expanded = true;
                break;
            end
            c_lo = c_lo / 1.5;
            c_hi = c_hi * 1.5;
            f_lo = obj(c_lo);
            f_hi = obj(c_hi);
        end

        if ~expanded
            % warning('Ishiyama21_zero:noBracket', ...
            %     ['No valid root bracket found for halo %d (n_eff=%.3f, rhs=%.4e, ' ...
            %      'profile=%s). Returning NaN.'], ...
            %     i, n_i, rhs_i, profile_name);
            c_unnorm(i) = NaN;
            continue;
        end

        c_unnorm(i) = fzero(obj, [c_lo, c_hi]);
    end

    c = reshape(Ceff_col .* c_unnorm, size(M));
end

