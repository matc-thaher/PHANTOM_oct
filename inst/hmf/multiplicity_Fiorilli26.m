function f = multiplicity_Fiorilli26(sigma, z, mdef, cosmo, include_unbound)
% multiplicity_Fiorilli26   Evolution Mapping halo multiplicity function
%
%   f = multiplicity_Fiorilli26(sigma, z, mdef, cosmo)
%   f = multiplicity_Fiorilli26(sigma, z, mdef, cosmo, true)
%
%   Implements the HMF model of Fiorilli, Ruiz, Sanchez & Esposito (2026),
%   arXiv:2511.16730. The non-universality is captured by two quantities:
%
%     x_tilde  — a Gaussian-kernel integral over the growth factor history,
%                encoding the recent rate of structure formation
%     n_eff    — the local logarithmic slope of P(k) at the scale of
%                collapse (nu = 1), encoding the P(k) shape
%
%   The functional form is:
%
%     f(nu) = A0 * nu * (A*nu^a + B*nu^b) * exp(-C * nu^2)
%
%   where A, a, B, b are linear functions of n_eff and x_tilde, and all
%   parameters depend on the overdensity threshold Delta.
%
% INPUTS:
%   sigma          : rms linear fluctuation (array)
%   z              : redshift (scalar)
%   mdef           : mass definition string, e.g. '200m', '300m', '1000m'
%                    Must be a spherical overdensity wrt mean density.
%                    Supported range: 150m -- 1600m.
%   cosmo          : PHANTOM cosmo struct with fields:
%                      cosmo.growthFactor(z)            — D(z)
%                      cosmo.growthFactor(z,'deriv',1)  — dD/dz
%                      cosmo.Om(z)                      — Omega_m(z)
%                      cosmo.sigmaR(R, z)               — sigma(R, z)
%   include_unbound: (optional, default false) include unbound particles
%                    in the fit. Activates the SO parameter set.
%
% OUTPUT:
%   f    : multiplicity function f(sigma), same size as sigma
%
% Reference:
%   Fiorilli, Ruiz, Sanchez & Esposito (2026), arXiv:2511.16730

    if nargin < 5 || isempty(include_unbound)
        include_unbound = false;
    end

    % --- Parse Delta from mdef -------------------------------------------
    % Expects format like '200m', '1000m'. Extracts the integer.
    tok = regexp(mdef, '^(\d+)m$', 'tokens');
    if isempty(tok)
        error('multiplicity_Fiorilli26: mdef must be of the form ''NNNm'' (e.g. ''200m''). Got: %s', mdef);
    end
    Delta_m = str2double(tok{1}{1});

    if Delta_m < 150 || Delta_m > 1600
        error('multiplicity_Fiorilli26: overdensity Delta=%g is outside the calibrated range [150, 1600].', Delta_m);
    end

    % --- Tabulated parameters --------------------------------------------
    if include_unbound
        pars_dict = fiorilli26_params_unbound();
    else
        pars_dict = fiorilli26_params_bound();
    end

    % --- Interpolate or look up parameters for Delta_m -------------------
    key = sprintf('d%d', Delta_m);
    if isfield(pars_dict, key)
        params = pars_dict.(key);
    else
        params = fiorilli26_interpolate(Delta_m, include_unbound);
    end

    % --- Collapse overdensity and peak height ----------------------------
    delta_c = 1.686;   % EdS approximation; replace with cosmo.delta_c(z) if available
    nu      = delta_c ./ sigma(:);

    % --- Effective power spectrum slope at nu = 1 ------------------------
    % k_eff at nu=1: sigma(R,z)=delta_c at the pivot scale.
    % Use a fixed pivot R corresponding to sigma(R,z)~delta_c at z.
    % Here we adopt a numerical derivative of ln(sigma) wrt ln(R) at the
    % pivot radius, then convert: n_eff = -2*d(ln sigma)/d(ln R) - 3
    R_pivot = fzero(@(R) cosmo.sigmaR(R, z) - delta_c, [0.01, 100], ...
                optimset('TolX', 1e-6, 'Display', 'off'));   % R where sigma(R,z)=delta_c
    eps_R    = 1e-3;
    sig_hi   = cosmo.sigmaR(R_pivot * exp( eps_R), z);
    sig_lo   = cosmo.sigmaR(R_pivot * exp(-eps_R), z);
    dlns_dlnR = (log(sig_hi) - log(sig_lo)) ./ (2 * eps_R);
    n_eff    = -2 .* dlns_dlnR - 3;

    % --- x_tilde: formation history integral -----------------------------
    x_tilde = compute_x_tilde(z, params.memory, cosmo);

    % --- Assemble parameters ---------------------------------------------
    B0 = -0.6;
    if include_unbound && ~strcmp(mdef, '150m')
        B0 = -0.55;
    end

    A_par = 1          + (params.mA_neff * n_eff);
    a_par =              params.ma_neff * n_eff;
    B_par = B0         + (params.mB_int  * x_tilde) + (params.mB_neff * n_eff);
    b_par =              params.mb_neff * n_eff;

    % --- Multiplicity function -------------------------------------------
    f = params.A0 .* nu .* ((A_par .* nu.^a_par) + (B_par .* nu.^b_par)) ...
        .* exp(-params.C .* nu.^2);

    f = reshape(f, size(sigma));
end


% =========================================================================
% Private helpers
% =========================================================================

function x_tilde = compute_x_tilde(z_final, memory, cosmo)
% Gauss-Legendre quadrature of the x_tilde formation history integral.
% Eq. (X) of Fiorilli+2026: a Gaussian-kernel-weighted integral over
% Omega_m / f_{growth} as a function of past redshift, in log(1+z) space.
%
%   x_tilde = (sqrt(2)/(memory*sqrt(pi))) *
%             integral_{z_final}^{z_max} kernel(z') * Om(z')/f(z') dz'
%
% where kernel = exp(-0.5*(ln D(z') - ln D(z_final))^2 / memory^2)
% and f = -((1+z)/D) * dD/dz is the growth rate.

    N        = 256;
    logz_max = log(1 + z_final);
    logz_min = max(2.0 + logz_max, log(1 + 20.0));

    [xi, w]  = gauss_legendre(N);   % points and weights on [-1, 1]
    logz     = 0.5 * (logz_max - logz_min) .* xi + 0.5 * (logz_max + logz_min);
    zprime   = exp(logz) - 1;

    D_final  = cosmo.D(z_final);
    D        = arrayfun(@(zp) cosmo.D(zp), zprime);
    % central finite difference for dD/dz since cosmo.D has no derivative flag
    eps_z    = 1e-4;
    dDdz     = arrayfun(@(zp) (cosmo.D(zp + eps_z) - cosmo.D(zp - eps_z)) / (2 * eps_z), zprime);

    f_growth    = -(1 + zprime) ./ D .* dDdz;
    Om = arrayfun(@(zp) cosmo.Omega_m_z(zp), zprime);
    delta_logD  = log(D) - log(D_final);

    ker      = exp(-0.5 .* delta_logD.^2 ./ memory^2);
    integrand   = ker .* Om ./ f_growth .* (-1);
    prefactor   = sqrt(2) / (memory * sqrt(pi));
    x_tilde     = prefactor .* sum(integrand .* w) .* 0.5 .* (logz_max - logz_min);
end


function [x, w] = gauss_legendre(N)
% Gauss-Legendre quadrature nodes and weights on [-1, 1].
% Uses the eigenvalue method (standard, O(N^2)).
    i   = (1:N-1)';
    b   = i ./ sqrt(4*i.^2 - 1);
    J   = diag(b, 1) + diag(b, -1);
    [V, D] = eig(J);
    x   = diag(D);
    w   = 2 .* V(1,:)'.^2;
    [x, idx] = sort(x);
    w   = w(idx);
end


function params = fiorilli26_interpolate(Delta_m, include_unbound)
% Interpolate Fiorilli+2026 parameters for arbitrary Delta_m using
% the fitting functions from their Table B1 / Appendix B.
% Each parameter has a distinct functional form fitted over log10(Delta).

    log10D = log10(Delta_m);

    if include_unbound
        % SO (unbound) coefficient set
        mB_int  = -0.67728012 * log10D^2 + 3.75566407 * log10D - 4.94960408 ...
                  + 0.25495863 / (log10D - 2);
        ma_neff = -0.03773026 * log10D + 0.30843905;
        mb_neff = -0.9344137  * (log10D - 2.12605188)^0.17719982;
        C       =  0.03821661 * log10D^2 + 0.2379606;
        memory  = power_delta_interp(Delta_m, -1.39884619, 0.5369116);
        mB_neff =  0.31280277 * log10D - 0.1894617;
        A0      =  0.44974526 * log10D^(-0.81108109);
        mA_neff = -0.43998955 * log10D + 0.45321855;
    else
        % Bound-only coefficient set
        mB_int  = power_interp( log10D, 8.11141425, -2.48399173);
        A0      = linear_interp(log10D, 0.49928071, -0.66872544);
        mA_neff = power_offset_interp(log10D, 0.29173886, 0.67490284, 2.1300181);
        ma_neff = linear_interp(log10D, -0.01925892, 0.27439948);
        mb_neff = power_offset_interp(log10D, -0.86082424, 0.19605053, 2.04425832);
        C       = quadratic_interp(log10D, 0.06849132, -0.169926, 0.47536706);
        memory  = power_delta_interp(Delta_m, -1.40294073, 0.65348985);
        mB_neff = linear_interp(log10D, 0.28805721, -0.19086668);
        mA_neff_from_mB = 0.29173886 * (mB_int - 2.1300181)^0.67490284;
        mA_neff = mA_neff_from_mB;   % derived from mB_int per colossus logic
    end

    params = struct('A0', A0, 'mA_neff', mA_neff, 'ma_neff', ma_neff, ...
                    'mB_int', mB_int, 'mB_neff', mB_neff, 'mb_neff', mb_neff, ...
                    'C', C, 'memory', memory);
end


function y = linear_interp(x, a, b),       y = a * x + b;                       end
function y = power_interp(x, a, b),        y = a * x^b;                          end
function y = quadratic_interp(x, a, b, c), y = a * x^2 + b * x + c;             end
function y = power_offset_interp(x, a, b, c), y = a * (x - c)^b;                end
function y = power_delta_interp(delta, a, b)
    if delta < 1400
        y = 10^(a * (delta/1000)^b);
    else
        y = 0;
    end
end

function pars_dict = fiorilli26_params_unbound()
% Tabulated Fiorilli+2026 parameters for SO/unbound set (include_unbound=true).
    pars_dict = struct();
    pars_dict.d150  = struct('A0', 0.29198269, 'mA_neff', -0.17174035, ...
        'ma_neff', 0.24837154, 'mB_int', 1.58460952, 'mB_neff', 0.25799066, ...
        'mb_neff', -0.55326303, 'C', 0.41977053, 'memory', 10.^-0.52598578);
    pars_dict.d200  = struct('A0', 0.45957271, 'mA_neff', 0.03585533, ...
        'ma_neff', 0.20953602, 'mB_int', 0.94150542, 'mB_neff', 0.09651153, ...
        'mb_neff', -0.65669468, 'C', 0.43820068, 'memory', 10.^-0.61715944);
    pars_dict.d350  = struct('A0', 0.62105717, 'mA_neff', 0.15018411, ...
        'ma_neff', 0.20751836, 'mB_int', 0.69400279, 'mB_neff', 0.02885085, ...
        'mb_neff', -0.78662850, 'C', 0.48113462, 'memory', 10.^-0.76773492);
    pars_dict.d500  = struct('A0', 0.67199492, 'mA_neff', 0.17670884, ...
        'ma_neff', 0.20250074, 'mB_int', 0.62098254, 'mB_neff', 0.00735698, ...
        'mb_neff', -0.86698663, 'C', 0.52001207, 'memory', 10.^-0.94672069);
    pars_dict.d650  = struct('A0', 0.69723014, 'mA_neff', 0.19426423, ...
        'ma_neff', 0.21079469, 'mB_int', 0.58418105, 'mB_neff', -0.00417243, ...
        'mb_neff', -0.90109142, 'C', 0.54730761, 'memory', 10.^-1.08658836);
    pars_dict.d800  = struct('A0', 0.75664790, 'mA_neff', 0.22187978, ...
        'ma_neff', 0.20903328, 'mB_int', 0.52466727, 'mB_neff', -0.02521139, ...
        'mb_neff', -0.91235026, 'C', 0.56254499, 'memory', 10.^-1.22404472);
    pars_dict.d1000 = struct('A0', 0.82292715, 'mA_neff', 0.24610851, ...
        'ma_neff', 0.19017363, 'mB_int', 0.46059674, 'mB_neff', -0.04668131, ...
        'mb_neff', -0.89273329, 'C', 0.57568765, 'memory', 10.^-1.58110205);
    pars_dict.d1200 = struct('A0', 0.88388947, 'mA_neff', 0.26732326, ...
        'ma_neff', 0.18801687, 'mB_int', 0.42620312, 'mB_neff', -0.05892036, ...
        'mb_neff', -0.88487159, 'C', 0.59130843, 'memory', 10.^-1.89703952);
    pars_dict.d1600 = struct('A0', 0.99982238, 'mA_neff', 0.29953702, ...
        'ma_neff', 0.18831457, 'mB_int', 0.37450430, 'mB_neff', -0.08052848, ...
        'mb_neff', -0.92539003, 'C', 0.63060366, 'memory', 0.0);
end

function pars_dict = fiorilli26_params_bound()
% Tabulated Fiorilli+2026 parameters for bound-only set (include_unbound=false).
    pars_dict = struct();
    pars_dict.d150  = struct('A0', 0.420179, 'mA_neff', 0.039496, ...
        'ma_neff', 0.231785, 'mB_int', 1.209711, 'mB_neff', 0.147693, ...
        'mb_neff', -0.587214, 'C', 0.430895, 'memory', 10.^-0.365159);
    pars_dict.d200  = struct('A0', 0.473317, 'mA_neff', 0.082095, ...
        'ma_neff', 0.233696, 'mB_int', 1.031511, 'mB_neff', 0.102342, ...
        'mb_neff', -0.633970, 'C', 0.445163, 'memory', 10.^-0.535402);
    pars_dict.d350  = struct('A0', 0.616211, 'mA_neff', 0.164419, ...
        'ma_neff', 0.218697, 'mB_int', 0.771759, 'mB_neff', 0.035216, ...
        'mb_neff', -0.764383, 'C', 0.486745, 'memory', 10.^-0.736517);
    pars_dict.d500  = struct('A0', 0.680618, 'mA_neff', 0.200122, ...
        'ma_neff', 0.220196, 'mB_int', 0.683583, 'mB_neff', 0.008838, ...
        'mb_neff', -0.803405, 'C', 0.516894, 'memory', 10.^-0.924594);
    pars_dict.d650  = struct('A0', 0.713565, 'mA_neff', 0.221220, ...
        'ma_neff', 0.227979, 'mB_int', 0.643937, 'mB_neff', -0.003431, ...
        'mb_neff', -0.818353, 'C', 0.542734, 'memory', 10.^-1.075271);
    pars_dict.d800  = struct('A0', 0.790642, 'mA_neff', 0.247776, ...
        'ma_neff', 0.219600, 'mB_int', 0.574126, 'mB_neff', -0.026488, ...
        'mb_neff', -0.862872, 'C', 0.561443, 'memory', 10.^-1.126259);
    pars_dict.d1000 = struct('A0', 0.862865, 'mA_neff', 0.272441, ...
        'ma_neff', 0.211811, 'mB_int', 0.509796, 'mB_neff', -0.049462, ...
        'mb_neff', -0.830703, 'C', 0.575597, 'memory', 10.^-1.403230);
    pars_dict.d1200 = struct('A0', 0.846278, 'mA_neff', 0.273011, ...
        'ma_neff', 0.215391, 'mB_int', 0.509352, 'mB_neff', -0.046756, ...
        'mb_neff', -0.843784, 'C', 0.598936, 'memory', 10.^-1.651781);
    pars_dict.d1400 = struct('A0', 0.874217, 'mA_neff', 0.285563, ...
        'ma_neff', 0.215301, 'mB_int', 0.461444, 'mB_neff', -0.064397, ...
        'mb_neff', -0.793947, 'C', 0.597651, 'memory', 0.0);
    pars_dict.d1600 = struct('A0', 0.891903, 'mA_neff', 0.294428, ...
        'ma_neff', 0.228580, 'mB_int', 0.470713, 'mB_neff', -0.060580, ...
        'mb_neff', -0.827234, 'C', 0.626810, 'memory', 0.0);
end