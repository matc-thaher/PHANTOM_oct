function [conc_interp, lhs_min_interp, lhs_max_interp] = build_lhs_table(profile_name)
% BUILD_LHS_TABLE  Build the G-function inversion table (renamed from get_Gc_table).
%
%   [conc_interp, lhs_min_interp, lhs_max_interp] = build_lhs_table(profile_name)
%
%   profile_name : 'nfw' (default) | 'hernquist' | 'einasto'
%
%   Compatible with MATLAB R2020a+ and GNU Octave.

    persistent cache

    cache_key = lower(strrep(profile_name, '-', '_'));

    if ~isempty(cache) && isfield(cache, cache_key)
        s              = cache.(cache_key);
        conc_interp    = s.conc_interp;
        lhs_min_interp = s.lhs_min_interp;
        lhs_max_interp = s.lhs_max_interp;
        return;
    end

    % --- add PHANTOM profiles to path (Octave-compatible) ---
    this_dir = fileparts(mfilename('fullpath'));
    prof_dir = fullfile(this_dir, '..', 'profiles');
    addpath(prof_dir);

    % --- grid sizes ---
    n_G = 80;  n_n = 40;  n_c = 80;
    n     = linspace(-4.0, 0.0, n_n);      % 1 x n_n
    c_log = linspace(-1.0, 3.0, n_c);      % 1 x n_c  log10(c)
    c     = 10.^c_log;                      % 1 x n_c

    % --- mu(c) called directly from PHANTOM profile functions ---
    mu = profile_mu(c, profile_name);       % 1 x n_c

    % --- vectorized lhs via implicit expansion, no loop over n ---
    exponents = (5 + n) / 6;                               % 1 x n_n
    lhs = log10( c(:) ./ (mu(:).^exponents) );             % n_c x n_n

    % --- enforce monotonic ascending region in c for each n ---
    mask_ascending = true(size(lhs));
    mask_ascending(1:end-1, :) = diff(lhs, 1, 1) > 0;

    % --- global LHS range ---
    lhs_min  = min(lhs(:));
    lhs_max  = max(lhs(:));
    lhs_grid = linspace(lhs_min, lhs_max, n_G);            % 1 x n_G

    log10c_table = zeros(n_G, n_n);
    mins = zeros(1, n_n);
    maxs = zeros(1, n_n);

    for j = 1:n_n
        m     = mask_ascending(:, j);                       % n_c x 1 logical
        lhs_j = lhs(m, j);                                  % ascending column
        c_j   = c_log(m)';                                  % log10(c), column

        % lhs_j is ascending after masking -> first/last = min/max
        mins(j) = lhs_j(1);
        maxs(j) = lhs_j(end);

        maskG     = (lhs_grid >= mins(j)) & (lhs_grid <= maxs(j));
        lhs_valid = lhs_grid(maskG);
        c_valid   = interp1(lhs_j, c_j, lhs_valid, 'linear', 'extrap');

        % pre-fill column then overwrite — removes mask_low/mask_high
        log10c_table(:, j)                  = c_j(1);
        log10c_table(maskG, j)              = c_valid(:);
        log10c_table(lhs_grid > maxs(j), j) = c_j(end);
    end

    % --- 2D interpolant: (lhs_grid, n) -> log10(c) ---
    [lhs_ndgrid, n_grid] = ndgrid(lhs_grid, n);
    conc_interp    = griddedInterpolant(lhs_ndgrid, n_grid, log10c_table, 'linear');

    % --- 1D interpolants: n -> min/max valid lhs ---
    lhs_min_interp = griddedInterpolant(n, mins, 'linear');
    lhs_max_interp = griddedInterpolant(n, maxs, 'linear');

    % --- cache ---
    s.conc_interp    = conc_interp;
    s.lhs_min_interp = lhs_min_interp;
    s.lhs_max_interp = lhs_max_interp;
    if isempty(cache), cache = struct(); end
    cache.(cache_key) = s;
end

