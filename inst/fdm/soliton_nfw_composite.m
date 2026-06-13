function comp = soliton_nfw_composite(r, rho, r_c_guess, virRad, varargin)
% SOLITON_NFW_COMPOSITE  Composite soliton+NFW profile from simulation data.
%
% Procedure:
%   1. Fit soliton profile to the inner region (r < ~3.5 * r_c_guess).
%   2. Fit NFW profile to the outer region (r > ~3.5 * r_c_guess).
%   3. Find the intersection radius r_x where soliton(r_x) = NFW(r_x).
%   4. Build composite: soliton for r <= r_x, NFW for r > r_x.
%
% INPUTS
%   r          : radii [kpc], column vector
%   rho        : simulation density [Msun/kpc^3], column vector
%   r_c_guess  : initial estimate of soliton core radius [kpc]
%   virRad     : virial radius [kpc] (upper bound for NFW fit)
%
% NAME-VALUE OPTIONS
%   'mRangeSol'            : [m_min, m_max] for soliton cutoff scan
%                            default [2.0, 3.5]
%   'mRangeNFW'            : [m_min, m_max] for NFW cutoff scan
%                            default [3.5, 10.0]
%   'SelectionMetricSol'   : metric used to choose best m for soliton
%                            'chi2' (default) | 'rmse' | 'rss' | 'aic' | 'bic'
%   'SelectionMetricNFW'   : metric used to choose best m for NFW
%                            'chi2' (default) | 'rmse' | 'rss' | 'aic' | 'bic'
%   'RetryOnNoIntersection': retry with adjusted m-range if no crossing found
%                            false (default) | true
%   'RetryTarget'          : which range to adjust on each retry
%                            'soliton' (default) | 'nfw' | 'both'
%   'RetryField'           : which bound of the range to shift
%                            'min' (default) | 'max' | 'both'
%   'RetryStep'            : increment applied to the chosen bound per retry
%                            default 0.2
%   'RetryMax'             : maximum number of retry attempts
%                            default 10
%   'Verbose'              : print composite-level diagnostics
%                            true (default) | false
%   'VerboseFindBestM'     : print per-m scan diagnostics inside find_best_m
%                            true (default) | false
%
% NOTE
%   The full scan statistics are still computed inside find_best_m. These
%   options control which statistic is used to select best_m, whether the
%   code retries after a failed soliton-NFW intersection, and whether
%   diagnostic text is printed.
%
% OUTPUTS (struct)
%   comp.r_x                : intersection radius [kpc]
%   comp.rho_x              : density at intersection [Msun/kpc^3]
%   comp.rho_composite      : stitched density profile at all r [Msun/kpc^3]
%   comp.rho_soliton        : soliton model evaluated at all r
%   comp.rho_nfw            : NFW model evaluated at all r
%   comp.pfit_sol           : soliton fit parameters [rho0, rc]
%   comp.pfit_nfw           : NFW fit parameters [rhos, rs]
%   comp.r                  : input radius array
%   comp.r_cut_sol          : soliton fitting cutoff used [kpc]
%   comp.r_cut_nfw          : NFW fitting cutoff used [kpc]
%   comp.best_m_sol         : optimal cutoff multiplier for soliton
%   comp.best_m_nfw         : optimal cutoff multiplier for NFW
%   comp.best_score_sol     : best metric value for soliton scan
%   comp.best_score_nfw     : best metric value for NFW scan
%   comp.metric_list_sol    : metric values across soliton m scan
%   comp.metric_list_nfw    : metric values across NFW m scan
%   comp.stats_sol          : full statistics struct from soliton scan
%   comp.stats_nfw          : full statistics struct from NFW scan
%   comp.intersection_found : true if a genuine crossing was located
%   comp.intersection_method: 'polyxpoly' | 'intersectPolylines' |
%                             'sign-change/fzero' | 'fallback-r_cut_nfw'
%   comp.retry_used         : true if at least one retry was performed
%   comp.retry_count        : number of retries actually used
%   comp.mRangeSol_initial  : mRangeSol as originally supplied
%   comp.mRangeNFW_initial  : mRangeNFW as originally supplied
%   comp.mRangeSol_final    : mRangeSol used on the successful (or last) attempt
%   comp.mRangeNFW_final    : mRangeNFW used on the successful (or last) attempt
%   comp.selection_metric_sol : metric name used for soliton selection
%   comp.selection_metric_nfw : metric name used for NFW selection

  % -----------------------------------------------------------------------
  % Input parsing
  % -----------------------------------------------------------------------
  opt.mRangeSol             = [2.0, 3.5];
  opt.mRangeNFW             = [3.5, 10.0];
  opt.SelectionMetric       = 'chi2';
  opt.RetryOnNoIntersection = false;
  opt.RetryTarget           = 'soliton';
  opt.RetryField            = 'min';
  opt.RetryStep             = 0.2;
  opt.RetryMax              = 10;
  opt.Verbose               = true;
  opt.VerboseFindBestM      = true;

  for i = 1:2:length(varargin)
     opt.(varargin{i}) = varargin{i+1};
  end

  mRangeSol          = opt.mRangeSol;
  mRangeNFW          = opt.mRangeNFW;
  selectionMetricSol = lower(opt.SelectionMetric);
  selectionMetricNFW = lower(opt.SelectionMetric);
  verboseComposite   = logical(opt.Verbose);
  verboseFindBestM   = logical(opt.VerboseFindBestM);

  % -----------------------------------------------------------------------
  % Initialise output struct and retry bookkeeping
  % -----------------------------------------------------------------------
  r   = r(:);
  rho = rho(:);

  comp = struct();
  comp.intersection_found   = false;
  comp.intersection_method  = 'none';
  comp.retry_used           = false;
  comp.retry_count          = 0;
  comp.mRangeSol_initial    = mRangeSol;
  comp.mRangeNFW_initial    = mRangeNFW;
  comp.selection_metric = char(selectionMetricSol);

  max_attempts = 1;
  if opt.RetryOnNoIntersection
      max_attempts = opt.RetryMax + 1;
  end

  % -----------------------------------------------------------------------
  % Retry loop — runs once by default
  % -----------------------------------------------------------------------
  for attempt = 1:max_attempts

    % --- Step 1: find best cutoff for soliton ---
    [m_sol, best_score_sol, metric_list_sol, stats_sol] = ...
        find_best_m('soliton', r_c_guess, r, rho, mRangeSol, ...
                    'SelectionMetric', selectionMetricSol, ...
                    'Verbose', verboseFindBestM);
    r_cut_sol = m_sol * r_c_guess;
    res_sol   = fit_profile_generic(r, rho, 'soliton', r_cut_sol, r_cut_sol);
    p_sol     = res_sol.pfit;   % [rho0, rc]

    % --- Step 2: find best cutoff for NFW ---
    [m_nfw, best_score_nfw, metric_list_nfw, stats_nfw] = ...
        find_best_m('nfw', r_c_guess, r, rho, mRangeNFW, ...
                    'SelectionMetric', selectionMetricNFW, ...
                    'Verbose', verboseFindBestM);
    r_cut_nfw = m_nfw * r_c_guess;
    res_nfw   = fit_profile_generic(r, rho, 'nfw', r_cut_nfw, virRad);
    p_nfw     = res_nfw.pfit;   % [rhos, rs]

    % --- Step 3: find intersection radius ---
    % Build curves on a dense log-density grid and search for a crossing.
    r_search_min = min(r(r > 0));   % smallest physical radius in data
    r_dense      = logspace(log10(r_search_min), log10(virRad), 5000)';
    rho_s   = Soliton_profile(r_dense, p_sol(1), p_sol(2));
    rho_n   = NFW_profile(r_dense, p_nfw(1), p_nfw(2));
    y_sol   = log10(rho_s);
    y_nfw   = log10(rho_n);

    is_octave = exist('OCTAVE_VERSION', 'builtin') ~= 0;

    r_x   = NaN;
    rho_x = NaN;
    found_flag = false;
    used_intersection_method = '';

    % Option A: MATLAB Mapping Toolbox -> polyxpoly
    has_polyxpoly = ~is_octave && ...
                    exist('polyxpoly', 'file') == 2 && ...
                    license('test', 'MAP_Toolbox');

    % Option B: matgeom -> intersectPolylines
    has_intersectPolylines = exist('intersectPolylines', 'file') == 2;

    if has_polyxpoly
      [xi, yi] = polyxpoly(r_dense, y_sol, r_dense, y_nfw);
      if ~isempty(xi)
        [r_x, idx]              = max(xi);
        rho_x                   = 10.^yi(idx);
        used_intersection_method = 'polyxpoly';
        found_flag              = true;
      end

    elseif has_intersectPolylines
      P = intersectPolylines([r_dense, y_sol], [r_dense, y_nfw]);
      if ~isempty(P)
        [r_x, idx]              = max(P(:,1));
        rho_x                   = 10.^P(idx,2);
        used_intersection_method = 'intersectPolylines';
        found_flag              = true;
      end

    else
      diff_arr     = y_sol - y_nfw;
      sign_changes = find(diff(sign(diff_arr)) ~= 0);

      if ~isempty(sign_changes)
        i_cross = sign_changes(end);
        r_a     = r_dense(i_cross);
        r_b     = r_dense(i_cross + 1);
        diff_fn = @(rv) log10(Soliton_profile(rv, p_sol(1), p_sol(2))) - ...
                        log10(NFW_profile(rv, p_nfw(1), p_nfw(2)));
        try
          r_x = fzero(diff_fn, [r_a, r_b]);
        catch
          r_x = 0.5 * (r_a + r_b);
        end
        rho_x                   = Soliton_profile(r_x, p_sol(1), p_sol(2));
        used_intersection_method = 'sign-change/fzero';
        found_flag              = true;
      end
    end

    % --- Step 3 outcome ---
    if found_flag
      comp.intersection_found  = true;
      comp.intersection_method = used_intersection_method;
      comp.r_x                 = r_x;
      comp.rho_x               = rho_x;
      break
    end

    % No intersection on this attempt
    if attempt == max_attempts
      warning(['soliton_nfw_composite: no intersection found between soliton and NFW ', ...
               'after %d attempt(s). Using r_cut_nfw as fallback r_x.'], attempt);
      r_x   = r_cut_nfw;
      rho_x = Soliton_profile(r_x, p_sol(1), p_sol(2));
      used_intersection_method = 'fallback-r_cut_nfw';

      comp.intersection_found  = false;
      comp.intersection_method = used_intersection_method;
      comp.r_x                 = r_x;
      comp.rho_x               = rho_x;
      break
    end

    % Prepare next attempt
    comp.retry_used  = true;
    comp.retry_count = attempt;

    [mRangeSol, mRangeNFW] = snc_update_m_ranges(mRangeSol, mRangeNFW, ...
        opt.RetryTarget, opt.RetryField, opt.RetryStep);

    if verboseComposite
      fprintf(['No intersection on attempt %d. Retrying with ', ...
               'mRangeSol = [%.2f, %.2f], mRangeNFW = [%.2f, %.2f]\n'], ...
               attempt, mRangeSol(1), mRangeSol(2), mRangeNFW(1), mRangeNFW(2));
    end

  end   % retry loop

  % -----------------------------------------------------------------------
  % Unpack intersection result for stitching
  % -----------------------------------------------------------------------
  r_x   = comp.r_x;
  rho_x = comp.rho_x;

  if verboseComposite
    fprintf('Intersection: r_x = %.4f kpc | rho_x = %.4e Msun/kpc^3 | method = %s\n', ...
            r_x, rho_x, comp.intersection_method);
  end

  % -----------------------------------------------------------------------
  % Step 4: stitch composite profile
  % -----------------------------------------------------------------------
  rho_sol_full = Soliton_profile(r, p_sol(1), p_sol(2));
  rho_nfw_full = NFW_profile(r, p_nfw(1), p_nfw(2));

  rho_composite           = zeros(size(r));
  rho_composite(r <= r_x) = rho_sol_full(r <= r_x);
  rho_composite(r >  r_x) = rho_nfw_full(r >  r_x);

  % -----------------------------------------------------------------------
  % Package all outputs
  % -----------------------------------------------------------------------
  comp.r_cut_sol       = r_cut_sol;
  comp.r_cut_nfw       = r_cut_nfw;
  comp.best_m_sol      = m_sol;
  comp.best_m_nfw      = m_nfw;
  comp.best_score_sol  = best_score_sol;
  comp.best_score_nfw  = best_score_nfw;
  comp.metric_list_sol = metric_list_sol;
  comp.metric_list_nfw = metric_list_nfw;
  comp.stats_sol       = stats_sol;
  comp.stats_nfw       = stats_nfw;
  comp.mRangeSol_final = mRangeSol;
  comp.mRangeNFW_final = mRangeNFW;
  comp.rho_composite   = rho_composite;
  comp.rho_soliton     = rho_sol_full;
  comp.rho_nfw         = rho_nfw_full;
  comp.pfit_sol        = p_sol;
  comp.pfit_nfw        = p_nfw;
  comp.r               = r;

  % -----------------------------------------------------------------------
  % Plot
  % -----------------------------------------------------------------------
  figure;
  loglog(r, rho,           'k.',  'MarkerSize', 8,   'DisplayName', 'Simulation data');
  hold on;
  loglog(r, rho_sol_full,  'b--', 'LineWidth',  1.5, 'DisplayName', 'Soliton fit');
  loglog(r, rho_nfw_full,  'r--', 'LineWidth',  1.5, 'DisplayName', 'NFW fit');
  loglog(r, rho_composite, 'g-',  'LineWidth',  2.5, 'DisplayName', 'Composite');
  plot(r_x, rho_x, 'ko', 'MarkerSize', 10, 'MarkerFaceColor', 'y', ...
       'DisplayName', sprintf('r_x = %.3f kpc', r_x));
  xlabel('Radius [kpc]',             'FontSize', 14);
  ylabel('Density [M_{sun}/kpc^3]',  'FontSize', 14);
  legend('Location', 'southwest');
  title('Soliton + NFW Composite Profile');
  grid on;

end 