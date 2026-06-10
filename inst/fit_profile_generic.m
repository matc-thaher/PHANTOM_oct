function result = fit_profile_generic(r, rho, profile_name, r_cut_low, r_cut_high, varargin)
% FIT_PROFILE_GENERIC  MATLAB/Octave-compatible generic profile fitter.
%
% Fits supported density profiles to simulation data in LINEAR density
% space (matching the approach of lsqcurvefit on raw density values),
% then computes all goodness-of-fit statistics in log10 space.
%
% This approach is numerically more stable than fitting in log-space
% directly, and is consistent with standard practice in halo profile
% fitting (e.g. Schive et al. 2014, Navarro et al. 1997).
%
% BASIC USAGE
%   result = fit_profile_generic(r, rho, 'nfw',           r_lo, r_hi)
%   result = fit_profile_generic(r, rho, 'soliton',       r_lo, r_hi)
%   result = fit_profile_generic(r, rho, 'einasto_fit',   r_lo, r_hi)
%   result = fit_profile_generic(r, rho, 'hernquist_fit', r_lo, r_hi)
%   result = fit_profile_generic(r, rho, 'dk14_fit',      r_lo, r_hi)
%
% OPTIONAL NAME-VALUE ARGUMENTS
%   'alpha' : Einasto/DK14 shape parameter alpha_e.
%             Default = 0.18  (Gao et al. 2008, MNRAS 387, median nu~1.5)
%             Physical range: 0.05 -- 0.35
%             Used only for 'einasto_fit' and 'dk14_fit'.
%             Pass this when you have a prior from the halo peak height:
%               alpha_e = 0.155 + 0.0095 * nu^2   (Gao+2008, eq. 5)
%  'r2_in' : the first radius value from where user want to build the fit
%   'Verbose'         : true (default) | false
%
% EXAMPLES
%   % default alpha
%   res = fit_profile_generic(r, rho, 'einasto_fit', r_c, r_vir);
%
%   % user-set alpha
%   res = fit_profile_generic(r, rho, 'einasto_fit', r_c, r_vir, 'alpha', 0.22);
%   res = fit_profile_generic(r, rho, 'dk14_fit',    r_c, r_vir, 'alpha', 0.20);
%
% INPUTS
%   r             : radial bins [kpc], column vector
%   rho           : density [Msun/kpc^3], column vector
%   profile_name  : string — 'nfw' | 'soliton' | 'einasto_fit' |
%                             'hernquist_fit' | 'dk14_fit'
%   r_cut_low     : lower radius cutoff [kpc]
%                   Soliton: ignored (fit always starts from r > 0)
%                   NFW / others: fit uses r > r_cut_low
%   r_cut_high    : upper radius cutoff [kpc]
%                   Soliton: fit uses r <= r_cut_high
%                   NFW / others: fit uses r <= r_cut_high
%
% RADIAL MASK CONVENTION (FDM split)
%   soliton              : 0  < r <= r_cut_high   (inner core region)
%   nfw / einasto / etc. : r_cut_low < r <= r_cut_high  (outer halo)
%
% OUTPUTS
%   result.pfit        : best-fit parameter vector [p1, p2] or [p1,p2,p3]
%   result.rho_fit     : model density evaluated at ALL input r
%   result.r_eval      : radius grid used to evaluate rho_fit
%                        equals r if 'r2' not given, or extended
%                        grid starting from r2 if r2 < min(r)
%   result.profile     : profile name string
%   result.r_used      : radii actually used in the fit
%   result.rho_used    : data densities used in the fit
%   result.n_pts       : number of data points used  (N)
%   result.n_params    : number of free parameters   (k)
%   result.dof         : degrees of freedom  (N - k)
%
%   Goodness-of-fit statistics (computed in log10 space after fitting)
%   result.log_rss     : sum of squared log10 residuals
%   result.rmse        : root mean square error [dex]
%                        < 0.05 dex  (~12% error) — good fit
%                        > 0.15 dex  (~40% error) — review fit
%   result.chi2_red    : reduced chi-squared  (log_rss / dof)
%                        ~ 1  statistically consistent
%                        >> 1 systematic misfit
%   result.R2          : coefficient of determination in log10 space
%                        > 0.99 recommended for clean profile fits
%   result.max_resid   : maximum absolute residual [dex]
%   result.AIC         : Akaike Information Criterion
%   result.BIC         : Bayesian Information Criterion
%                        Lower AIC/BIC = preferred model.
%                        |delta_AIC| > 10 is strong evidence.
%
% OCTAVE NOTE
%   Install the optim package once:
%       pkg install -forge optim
%   Load at the start of each Octave session:
%       pkg load optim
%   MATLAB uses lsqcurvefit. Octave uses leasqr from the optim package.

  % --- optional input: alpha -----------------------------------------
  opt.alpha   = 0.18;
  opt.r_ext   = [];
  opt.Verbose = false;
  for i = 1:2:length(varargin)
    opt.(varargin{i}) = varargin{i+1};
  end

  % --- environment detection -----------------------------------------
  is_octave = exist('OCTAVE_VERSION', 'builtin') ~= 0;

  % --- enforce column vectors ----------------------------------------
  r   = r(:);
  rho = rho(:);

  % --- radial mask (FDM split) ---------------------------------------
  switch lower(profile_name)
    case 'soliton'
      mask = (r > 0) & (r <= r_cut_low) & (rho > 0);
    case {'nfw', 'einasto_fit', 'hernquist_fit', 'dk14_fit'}
      mask = (r > r_cut_low) & (r <= r_cut_high) & (rho > 0);
    otherwise
      error('fit_profile_generic: unknown profile "%s".\nValid options: nfw, soliton, einasto_fit, hernquist_fit, dk14_fit.', ...
            profile_name);
  end

  r_fit   = r(mask);
  rho_fit = rho(mask);
  N       = numel(r_fit);

  if N < 5
  warning('fit_profile_generic: fewer than 5 data points after masking for "%s".', ...
          profile_name);

    result              = struct();
    result.pfit         = NaN;
    result.rho_fit      = NaN;
    result.r_eval       = NaN;
    result.profile      = profile_name;
    result.r_used       = r_fit;
    result.rho_used     = rho_fit;
    result.n_pts        = N;
    result.n_params     = NaN;
    result.dof          = NaN;
    result.log_rss      = Inf;
    result.rmse         = Inf;
    result.chi2_red     = Inf;
    result.R2           = NaN;
    result.max_resid    = NaN;
    result.AIC          = NaN;
    result.BIC          = NaN;
    return;
  end

  % --- profile selection, initial guess, parameter count -------------
  %
  % The fitting function signature MUST be fun(xdata, p) for lsqcurvefit.
  % Profile functions are wrapped to enforce this and to guarantee that
  % the output is always a column vector regardless of input shape —
  % this prevents lsqcurvefit from throwing size mismatch errors during
  % its internal Jacobian finite-difference steps.
  %
  % Fitting is done in LINEAR density space. Log-space statistics are
  % computed from the fitted model after convergence.

  switch lower(profile_name)
    case 'soliton'
      fun   = @(p, xdata) reshape(Soliton_profile(xdata(:), p(1), p(2)), [], 1);
      p0    = [max(rho_fit), mean(r_fit)];
      n_par = 2;

    case 'nfw'
      fun   = @(p, xdata) reshape(NFW_profile(xdata(:), p(1), p(2)), [], 1);
      p0    = [max(rho_fit), mean(r_fit)];
      n_par = 2;

    case 'einasto'
      fun   = @(p, xdata) reshape(Einasto_fit_wrapper(xdata(:), p(1), p(2), p(3)), [], 1);
      p0    = [max(rho_fit), mean(r_fit), alpha_in];
      n_par = 3;

    case 'hernquist'
      fun   = @(p, xdata) reshape(Hernquist_fit_wrapper(xdata(:), p(1), p(2)), [], 1);
      p0    = [max(rho_fit), mean(r_fit)];
      n_par = 2;

    case 'dk14'
      fun   = @(p, xdata) reshape(DK14_fit_wrapper(xdata(:), p(1), p(2), p(3)), [], 1);
      p0    = [max(rho_fit), mean(r_fit), alpha_in];
      n_par = 3;
  end

  dof = N - n_par;
  if dof < 1
    warning('fit_profile_generic: dof = %d for "%s". Fit is underdetermined.', ...
            dof, profile_name);
  end

  % --- solve: MATLAB vs Octave ---------------------------------------
  if ~is_octave
    % MATLAB — lsqcurvefit in linear density space
    opts = optimoptions('lsqcurvefit',      ...
                        'Display',           'off', ...
                        'MaxIterations',     2000,  ...
                        'FunctionTolerance', 1e-10, ...
                        'StepTolerance',     1e-10);
    try
      pfit = lsqcurvefit(fun, p0, r_fit, rho_fit, [], [], opts);
    catch ME
        warning('fit_profile_generic: lsqcurvefit failed for "%s".\n  MATLAB error: %s', ...
          profile_name, ME.message);

        result              = struct();
        result.pfit         = NaN;
        result.rho_fit      = NaN;
        result.r_eval       = NaN;
        result.profile      = profile_name;
        result.r_used       = r_fit;
        result.rho_used     = rho_fit;
        result.n_pts        = N;
        result.n_params     = n_par;
        result.dof          = dof;
        result.log_rss      = Inf;
        result.rmse         = Inf;
        result.chi2_red     = Inf;
        result.R2           = NaN;
        result.max_resid    = NaN;
        result.AIC          = NaN;
        result.BIC          = NaN;
        return;
    end

  else
    % Octave — leasqr from the optim package
    % leasqr signature: leasqr(x, y, p0, fun)
    % where fun must be fun(x, p) — same convention as lsqcurvefit
    %
    % Required:
    %   pkg install -forge optim   (once)
    %   pkg load optim             (each session)
    try
      [~, pfit] = leasqr(r_fit, rho_fit, p0, fun);
      pfit = pfit(:)';
    catch ME
        warning(['fit_profile_generic: leasqr failed for "%s".\n' ...
           '  Octave error: %s\n' ...
           '  Make sure the optim package is installed and loaded:\n' ...
           '    pkg install -forge optim\n' ...
           '    pkg load optim'], profile_name, ME.message);
        result          = struct();
        result.pfit     = NaN;
        result.rho_fit  = NaN;
        result.r_eval   = NaN;
        result.profile  = NaN;
        result.r_used   = NaN;
        result.rho_used = NaN;
        result.n_pts    = NaN;
        result.n_params = NaN;
        result.dof      = NaN;
        result.log_rss  = Inf;
        result.rmse     = Inf;
        result.chi2_red = Inf;
        result.profile  = profile_name;
        result.n_pts    = N;
        result.dof      = dof;
        return;
    end
  end

  % --- goodness-of-fit statistics in log10 space ---------------------
  rho_model = fun(pfit, r_fit);
  log_model  = log10(rho_model);
  log_data   = log10(rho_fit);
  residuals  = log_model - log_data;

  log_rss   = sum(residuals.^2);
  rmse      = sqrt(log_rss / N);
  chi2_red  = log_rss / max(dof, 1);
  log_mean  = mean(log_data);
  ss_tot    = sum((log_data - log_mean).^2);
  R2        = 1 - log_rss / ss_tot;
  max_resid = max(abs(residuals));

  % AIC and BIC under Gaussian residual assumption in log-space
  log_likelihood_proxy = -0.5 * N * log(log_rss / N);
  AIC = -2 * log_likelihood_proxy + 2 * n_par;
  BIC = -2 * log_likelihood_proxy + n_par * log(N);

  % --- build evaluation grid -----------------------------------------
  % If r2 is given and r2 < min(r), prepend a linspace from r2 to
  % min(r) with 20 points, then append the original r array.
  % The fit parameters do not change — only the grid is extended.
  if ~isempty(r2_in) && r2_in < min(r)
    r_inner   = linspace(r2_in, min(r), 20)';
    r_inner   = r_inner(1:end-1);   % drop last point to avoid duplicate at min(r)
    r_eval    = [r_inner; r];
  else
    r_eval    = r;
  end

  % --- package results -----------------------------------------------
  result.pfit      = pfit;
  result.rho_fit   = reshape(fun(pfit, r_eval), [], 1);
  result.r_eval    = r_eval;
  result.profile   = profile_name;
  result.r_used    = r_fit;
  result.rho_used  = rho_fit;
  result.n_pts     = N;
  result.n_params  = n_par;
  result.dof       = dof;
  result.log_rss   = log_rss;
  result.rmse      = rmse;
  result.chi2_red  = chi2_red;
  result.R2        = R2;
  result.max_resid = max_resid;
  result.AIC       = AIC;
  result.BIC       = BIC;

  % --- print summary -------------------------------------------------
  if opt.Verbose
    fprintf('\n--- Fit summary: %s ---\n',    profile_name);
    fprintf('  N = %d pts  |  dof = %d\n',   N, dof);
    fprintf('  log-RSS     = %.4e\n',         log_rss);
    fprintf('  RMSE        = %.4f dex\n',     rmse);
    fprintf('  chi2_red    = %.4f\n',         chi2_red);
    fprintf('  R2          = %.6f\n',         R2);
    fprintf('  max|resid|  = %.4f dex\n',     max_resid);
    fprintf('  AIC = %.2f  |  BIC = %.2f\n',  AIC, BIC);
  end
end
