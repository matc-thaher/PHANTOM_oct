function comp = soliton_nfw_analytic(rho0_sol, rc_sol, Mvir, Rvir, c_nfw, r, plot_fig)
% SOLITON_NFW_ANALYTIC  Analytic soliton + NFW composite FDM profile.
%
%   comp = soliton_nfw_analytic(rho0_sol, rc_sol, Mvir, Rvir, c_nfw, r)
%
%   Constructs the composite density profile from analytic parameters:
%     - Soliton  : central density rho0 and core radius rc  (Schive+2014)
%     - NFW halo : derived analytically from Mvir, Rvir, concentration c
%
%   The transition radius r_x (soliton = NFW) is found with a three-tier
%   method that is compatible with both Octave and MATLAB:
%     Tier 1  polyxpoly        (MATLAB Mapping Toolbox, if available)
%     Tier 2  intersectPolylines (matgeom, if on path)
%     Tier 3  sign-change + fzero / midpoint fallback  (always available)
%
% --------------------------------------------------------------------------
% INPUTS
%   rho0_sol : soliton central density         [M_sun / kpc^3]
%   rc_sol   : soliton core radius             [kpc]
%   Mvir     : virial mass                     [M_sun]
%   Rvir     : virial radius                   [kpc]
%   c_nfw    : NFW concentration               [dimensionless]
%   r        : radii for output profile        [kpc], column or row vector
%
% options (optional):
%   plot_fig : True or false   ,   if the user want a plot
%
% OUTPUTS  (struct)
%   comp.r_x                 : transition radius              [kpc]
%   comp.rho_x               : density at r_x                 [M_sun/kpc^3]
%   comp.rho_composite       : composite profile at r         [M_sun/kpc^3]
%   comp.rho_soliton         : soliton profile at r           [M_sun/kpc^3]
%   comp.rho_nfw             : NFW profile at r               [M_sun/kpc^3]
%   comp.rs_nfw              : NFW scale radius               [kpc]
%   comp.rhos_nfw            : NFW scale density              [M_sun/kpc^3]
%   comp.r                   : input radius array             [kpc]
%   comp.intersection_found  : logical
%   comp.intersection_method : string identifying tier used
%
% --------------------------------------------------------------------------
% REFERENCES
%   Schive, Chiueh & Broadhurst 2014, Nat. Phys. 10, 496  (soliton profile)
%   Navarro, Frenk & White 1997, ApJ 490, 493              (NFW profile)
%   Schive+2014 PRL 113 261302, Eq.(3)                     (composite rule)
% --------------------------------------------------------------------------
    
    % ---- Optional plot flag ------------------------------------------
    if nargin < 7 || isempty(plot_fig)
        plot_fig = false;
    end

    r = r(:);   % force column vector

    % ------------------------------------------------------------------
    % Environment detection
    % ------------------------------------------------------------------
    is_octave = exist('OCTAVE_VERSION', 'builtin') ~= 0;

    % ------------------------------------------------------------------
    % Soliton profile on output grid
    % rho(r) = rho0 * [1 + 0.091*(r/rc)^2]^-8
    % Schive+2014 Nat.Phys., Eq.(S4) / Mocz+2017 Eq.(14)
    % ------------------------------------------------------------------
    rho_sol = Soliton_profile(r, rho0_sol, rc_sol);

    % ------------------------------------------------------------------
    % NFW analytic parameters from Mvir, Rvir, c
    % rs    = Rvir / c
    % rho_s = Mvir / (4*pi*rs^3 * f(c)),  f(c) = ln(1+c) - c/(1+c)
    % NFW 1997, ApJ 490, 493, Eq.(7)
    % ------------------------------------------------------------------
    nfw_struct = NFW_analytcl_Profile(Mvir, Rvir, c_nfw, r);
    rs_nfw     = nfw_struct.rs;
    rhos_nfw   = nfw_struct.rho_s;

    % Evaluate NFW on the full r grid (NFW_analytcl_Profile trims r > Rvir)
    rho_nfw = NFW_profile(r, rhos_nfw, rs_nfw);

    % ------------------------------------------------------------------
    % Find intersection r_x  —  three-tier, Octave-safe
    %
    % Dense log-spaced grid from r_min to Rvir.
    % Working in log10 space for numerical stability.
    % ------------------------------------------------------------------
    r_dense  = logspace(log10(min(r)), log10(Rvir), 3000)';
    y_sol    = log10(Soliton_profile(r_dense, rho0_sol, rc_sol));
    y_nfw    = log10(NFW_profile(r_dense, rhos_nfw, rs_nfw));

    r_x   = NaN;
    rho_x = NaN;
    found_flag               = false;
    used_intersection_method = '';

    % -- Tier 1: MATLAB Mapping Toolbox polyxpoly ----------------------
    has_polyxpoly = ~is_octave && ...
                    exist('polyxpoly', 'file') == 2 && ...
                    license('test', 'MAP_Toolbox');

    % -- Tier 2: matgeom intersectPolylines ----------------------------
    has_intersectPolylines = exist('intersectPolylines', 'file') == 2;

    if has_polyxpoly
        [xi, yi] = polyxpoly(r_dense, y_sol, r_dense, y_nfw);
        if ~isempty(xi)
            % take the outermost crossing (largest r) — soliton < NFW
            % beyond this point
            [r_x, idx]               = max(xi);
            rho_x                    = 10.^yi(idx);
            used_intersection_method = 'polyxpoly';
            found_flag               = true;
        end

    elseif has_intersectPolylines
        P = intersectPolylines([r_dense, y_sol], [r_dense, y_nfw]);
        if ~isempty(P)
            [r_x, idx]               = max(P(:,1));
            rho_x                    = 10.^P(idx, 2);
            used_intersection_method = 'intersectPolylines';
            found_flag               = true;
        end

    else
        % -- Tier 3: sign-change + fzero (Octave + MATLAB, no toolbox) --
        diff_arr     = y_sol - y_nfw;
        sign_changes = find(diff(sign(diff_arr)) ~= 0);

        if ~isempty(sign_changes)
            % outermost crossing
            i_cross = sign_changes(end);
            r_a     = r_dense(i_cross);
            r_b     = r_dense(i_cross + 1);

            diff_fn = @(rv) log10(Soliton_profile(rv, rho0_sol, rc_sol)) - ...
                            log10(NFW_profile(rv, rhos_nfw, rs_nfw));
            try
                r_x = fzero(diff_fn, [r_a, r_b]);
            catch
                r_x = 0.5 * (r_a + r_b);
            end
            rho_x                    = Soliton_profile(r_x, rho0_sol, rc_sol);
            used_intersection_method = 'sign-change/fzero';
            found_flag               = true;
        end
    end

    % -- Fallback if all tiers fail -------------------------------------
    if ~found_flag
        warning(['soliton_nfw_analytic: profiles do not cross on ' ...
                 '[r_min, Rvir]. Using r_x = 3*rc_sol as fallback. '  ...
                 'Check rho0_sol and NFW amplitude are consistent.']);
        r_x                      = 3 * rc_sol;
        rho_x                    = Soliton_profile(r_x, rho0_sol, rc_sol);
        used_intersection_method = 'fallback:3rc';
    end

    fprintf('Analytic intersection: r_x = %.4f kpc | rho_x = %.4e M_sun/kpc^3 [%s]\n', ...
            r_x, rho_x, used_intersection_method);

    % ------------------------------------------------------------------
    % Build composite profile
    % Soliton for r <= r_x,  NFW for r > r_x
    % Schive+2014 PRL 113 261302, Eq.(3)
    % ------------------------------------------------------------------
    rho_composite           = zeros(size(r));
    rho_composite(r <= r_x) = rho_sol(r <= r_x);
    rho_composite(r >  r_x) = rho_nfw(r >  r_x);

    % ------------------------------------------------------------------
    % Package output
    % 
    % ------------------------------------------------------------------
    comp.r_x                 = r_x;
    comp.rho_x               = rho_x;
    comp.rho_composite       = rho_composite;
    comp.rho_soliton         = rho_sol;
    comp.rho_nfw             = rho_nfw;
    comp.rs_nfw              = rs_nfw;
    comp.rhos_nfw            = rhos_nfw;
    comp.r                   = r;
    comp.intersection_found  = found_flag;
    comp.intersection_method = used_intersection_method;

    if plot_fig
        % ------------------------------------------------------------------
        % Plot
        % ------------------------------------------------------------------
        figure;
        loglog(r, rho_sol,       'b--', 'LineWidth', 1.5, ...
           'DisplayName', 'Soliton (analytic)');
        hold on;
        loglog(r, rho_nfw,       'r--', 'LineWidth', 1.5, ...
           'DisplayName', 'NFW (analytic)');
        loglog(r, rho_composite, 'g-',  'LineWidth', 2.5, ...
           'DisplayName', 'Composite');
        plot(r_x, rho_x, 'ko', 'MarkerSize', 10, 'MarkerFaceColor', 'y', ...
         'DisplayName', sprintf('r_x = %.3f kpc', r_x));
        xlabel('Radius [kpc]', 'FontSize', 14);
        ylabel('Density [M_{sun}/kpc^3]', 'FontSize', 14);
        legend('Location', 'southwest');
        title(sprintf('Analytic Composite: c = %.1f,  r_c = %.3f kpc', ...
                  c_nfw, rc_sol));
        % grid on;
    end
end
