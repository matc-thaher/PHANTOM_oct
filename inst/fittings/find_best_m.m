function [best_m, best_score, metric_list, stats] = find_best_m(profile_name, rc, r, rho, m_range, varargin)
% FIND_BEST_M  Scan cutoff multiplier m to find optimal r_cut = m * rc.
%
% By default, the optimal cutoff is chosen using reduced chi-squared, but
% the user may select any supported fit statistic through 'SelectionMetric'.
%
% INPUTS
%   profile_name : 'nfw' | 'soliton' | 'einasto_fit' | 'hernquist_fit' | 'dk14_fit'
%   rc           : reference radius (e.g., soliton core radius) [kpc]
%   r            : full radius array [kpc], column vector
%   rho          : full density array [Msun/kpc^3], column vector
%   m_range      : 1x2 vector [m_min, m_max], e.g. [2, 10]
% NAME-VALUE OPTIONS
%   'SelectionMetric' : 'chi2' (default) | 'rmse' | 'rss' | 'aic' | ...
%                       'bic' | 'maxresid' | 'r2'
%   'Verbose'         : true (default) | false
%
% OUTPUTS
%   best_m       : optimal cutoff multiplier
%   best_score   : best value of the chosen metric
%   metric_list  : values of the chosen metric across all m tested
%   stats        : struct containing all scanned statistics
%
% WHY chi2_red INSTEAD OF raw log-RSS?
%   Raw log-RSS decreases whenever data points are removed, so a fitter
%   with a very high r_cut (few points) always appears to win. chi2_red
%   = log_RSS / dof penalises models that fit only a handful of points,
%   giving a fair comparison across different cutoffs.
%
% NOTE ON DEGREES OF FREEDOM
%   dof = N_pts - N_params. If dof < 1 the fit is underdetermined and
%   that m value is skipped (set to Inf in chi2_list).

  opt.SelectionMetric = 'chi2';
  opt.Verbose         = true;
  for i = 1:2:length(varargin)
    opt.(varargin{i}) = varargin{i+1};
  end
  
  selection_metric = lower(opt.SelectionMetric);

  m_values = linspace(m_range(1), m_range(2), 100);
  r_max    = max(r);

  n = numel(m_values);
  stats = struct();
  stats.m_values   = m_values;
  stats.rss        = Inf(1, n);
  stats.chi2       = Inf(1, n);
  stats.rmse       = Inf(1, n);
  stats.R2         = -Inf(1, n);
  stats.maxresid   = Inf(1, n);
  stats.AIC        = Inf(1, n);
  stats.BIC        = Inf(1, n);

  for i = 1:n
    r_cut = m_values(i) * rc;

    try
      res = fit_profile_generic(r, rho, profile_name, r_cut, r_max);

      if res.dof < 1 || any(~isfinite([res.log_rss, res.rmse, res.chi2_red]))
        continue
      end

      stats.rss(i)      = res.log_rss;
      stats.chi2(i)     = res.chi2_red;
      stats.rmse(i)     = res.rmse;
      stats.R2(i)       = res.R2;
      stats.maxresid(i) = res.max_resid;
      stats.AIC(i)      = res.AIC;
      stats.BIC(i)      = res.BIC;

    catch
      continue
    end
  end

  switch selection_metric
    case "chi2"
      metric_list = stats.chi2;
      [best_score, idx] = min(metric_list);
      metric_label = 'chi2_red';

    case "rmse"
      metric_list = stats.rmse;
      [best_score, idx] = min(metric_list);
      metric_label = 'RMSE';

    case "rss"
      metric_list = stats.rss;
      [best_score, idx] = min(metric_list);
      metric_label = 'log-RSS';

    case "aic"
      metric_list = stats.AIC;
      [best_score, idx] = min(metric_list);
      metric_label = 'AIC';

    case "bic"
      metric_list = stats.BIC;
      [best_score, idx] = min(metric_list);
      metric_label = 'BIC';

    case "maxresid"
      metric_list = stats.maxresid;
      [best_score, idx] = min(metric_list);
      metric_label = 'max residual';

    case "r2"
      metric_list = stats.R2;
      [best_score, idx] = max(metric_list);
      metric_label = 'R^2';

    otherwise
      error(['Unknown SelectionMetric: %s. Use ''chi2'', ''rmse'', ''rss'', ', ...
             '''aic'', ''bic'', ''maxresid'', or ''r2''.'], selection_metric);
  end

  if ~isfinite(best_score)
    best_m = NaN;
    warning('find_best_m: no valid fit found for profile %s in the supplied m_range.', profile_name);
    return
  end

  best_m = m_values(idx);

  if opt.Verbose
    fprintf('\n--- find_best_m: %s ---\n', profile_name);
    fprintf('  Selection   = %s\n', metric_label);
    fprintf('  Best m      = %.2f\n', best_m);
    fprintf('  r_cut       = %.4f kpc\n', best_m * rc);
    fprintf('  chi2_red    = %.4f\n', stats.chi2(idx));
    fprintf('  log-RSS     = %.4e\n', stats.rss(idx));
    fprintf('  RMSE        = %.4f dex\n', stats.rmse(idx));
    fprintf('  R^2         = %.4f\n', stats.R2(idx));
    fprintf('  max_resid   = %.4f dex\n', stats.maxresid(idx));
    fprintf('  AIC         = %.4f\n', stats.AIC(idx));
    fprintf('  BIC         = %.4f\n', stats.BIC(idx));
  end
end
