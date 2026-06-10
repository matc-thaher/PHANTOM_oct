function prof = profile_obs(r, rho, quantity_out, cosmo_G)
% PROFILE_OBS  Derive radial profiles from any input density profile rho(r).
%
%   prof = profile_obs(r, rho, quantity_out)
%   prof = profile_obs(r, rho, quantity_out, G)
%
%   Computes enclosed mass, circular velocity, 3-D velocity dispersion
%   (isotropic Jeans), and projected surface density from an arbitrary
%   spherically symmetric density profile rho(r).
%
%   The function is profile-agnostic: pass NFW, soliton, composite FDM,
%   Einasto, Hernquist, or any numerically sampled rho(r).
%
%   A log-spaced radial grid is strongly recommended (use logspace) so
%   that both the inner core and the outer envelope are well resolved.
%
% --------------------------------------------------------------------------
% INPUTS
%   r            : radial grid [kpc], column or row vector, N points
%                  Must be monotonically increasing, no zeros.
%   rho          : density at each r [M_sun / kpc^3], same size as r
%   quantity_out : string or cell array of strings (see list below)
%                  Use 'all' to return every supported quantity.
%   G            : (optional) gravitational constant
%                  [kpc (km/s)^2 / M_sun], default = 4.3009e-6
%                  Pass your hc.G from halo_cosmo if available.
%
% --------------------------------------------------------------------------
% SUPPORTED quantity_out STRINGS
%   'M_enc'    Enclosed mass M(<r)                      [M_sun]
%   'V_circ'   Circular velocity V_c(r) = sqrt(GM/r)   [km/s]
%   'sigma_r'  3-D radial velocity dispersion, isotropic Jeans equation
%              sigma_r^2(r) = (1/rho) * int_r^inf rho * GM/r'^2 dr'
%                                                       [km/s]
%   'sigma_los' Line-of-sight (projected) velocity dispersion
%              sigma_los^2(R) = (2/Sigma) * int_R^inf rho*sigma_r^2 *
%                               r/sqrt(r^2-R^2) dr      [km/s]
%   'Sigma'    Projected surface density
%              Sigma(R) = 2 * int_R^inf rho*r/sqrt(r^2-R^2) dr
%                                                       [M_sun / kpc^2]
%   'all'      Returns all of the above
%
% --------------------------------------------------------------------------
% OUTPUT  struct with fields matching quantity_out strings, plus:
%   prof.r     : input radial grid [kpc]
%   prof.rho   : input density     [M_sun/kpc^3]
%
% --------------------------------------------------------------------------
% EQUATIONS
%   Enclosed mass (spherical shell integration):
%     M(<r) = 4*pi * int_0^r rho(r') r'^2 dr'
%     Computed via cumulative trapezoidal rule on the (log r, r^3*rho) grid
%     for accuracy across many decades.
%     Reference: standard; see e.g. Binney & Tremaine 2008, Eq.(2.22)
%
%   Circular velocity:
%     V_c(r) = sqrt(G * M(<r) / r)
%     Binney & Tremaine 2008, Eq.(2.30)
%
%   Isotropic Jeans equation (no anisotropy, beta=0):
%     d(rho*sigma_r^2)/dr = -rho * G*M(<r)/r^2
%     Integrated outward from r to infinity with boundary condition
%     sigma_r -> 0 as r -> inf:
%       sigma_r^2(r) = (1/rho(r)) * int_r^inf rho(r') * G*M(<r')/r'^2 dr'
%     Binney & Tremaine 2008, Eq.(4.215)
%
%   Projected surface density (Abel transform):
%     Sigma(R) = 2 * int_R^inf rho(r) * r / sqrt(r^2 - R^2) dr
%     Binney & Tremaine 2008, Eq.(2.135)
%     Computed by direct numerical quadrature; the integrand diverges as
%     r -> R^+ so the lower limit is offset to R*(1 + epsilon).
%
%   Line-of-sight velocity dispersion:
%     sigma_los^2(R) = (2/Sigma(R)) * int_R^inf rho(r)*sigma_r^2(r) *
%                      r / sqrt(r^2 - R^2) dr
%     Binney & Tremaine 2008, Eq.(4.225)
%     Same Abel-type integral as Sigma, weighted by rho*sigma_r^2.
%
% --------------------------------------------------------------------------
% NOTES ON NUMERICAL ACCURACY
%   All integrals use the trapezoidal rule on a log-spaced grid.
%   Transforming to log(r) converts the power-law integrands to smooth
%   functions that trapz handles accurately without large point counts.
%   For the Abel integrals (Sigma, sigma_los), each projected radius R_i
%   uses only the sub-array r >= R_i, so the singularity at r = R is
%   naturally avoided by the finite grid spacing.
%
%   Recommended grid: logspace(log10(r_min), log10(r_max), 500)
%   r_min should be well inside the core (e.g. rc/100 for FDM solitons).
%   r_max should reach at least the virial radius.
%
% --------------------------------------------------------------------------
% REFERENCES
%   Binney & Tremaine (2008) Galactic Dynamics, 2nd ed., Princeton
%   Navarro, Frenk & White 1997, ApJ 490, 493
%   Schive, Chiueh & Broadhurst 2014, Nat. Phys. 10, 496
%
% --------------------------------------------------------------------------

    % ------------------------------------------------------------------
    % 0.  Defaults and input preparation
    % ------------------------------------------------------------------
    if nargin < 4 || isempty(cosmo_G)
        G = 4.3009e-6;   % kpc (km/s)^2 / M_sun
    else
        G = cosmo_G;
    end

    r   = r(:);
    rho = rho(:);

    if numel(r) ~= numel(rho)
        error('profile_obs: r and rho must have the same length.');
    end
    if any(r <= 0)
        error('profile_obs: all radii must be positive (no zeros).');
    end

    N = numel(r);

    % Resolve 'all'
    all_quantities = {'M_enc','V_circ','sigma_r','sigma_los','Sigma'};
    if ischar(quantity_out)
        if strcmpi(quantity_out, 'all')
            quantity_out = all_quantities;
        else
            quantity_out = {quantity_out};
        end
    end
    scalar_out = numel(quantity_out) == 1;

    % Determine which computations are needed
    need_M      = any(strcmpi(quantity_out,'M_enc'))   || ...
                  any(strcmpi(quantity_out,'V_circ'))   || ...
                  any(strcmpi(quantity_out,'sigma_r'))  || ...
                  any(strcmpi(quantity_out,'sigma_los'));
    need_sigma_r = any(strcmpi(quantity_out,'sigma_r')) || ...
                   any(strcmpi(quantity_out,'sigma_los'));
    need_Sigma   = any(strcmpi(quantity_out,'Sigma'))   || ...
                   any(strcmpi(quantity_out,'sigma_los'));

    % ------------------------------------------------------------------
    % 1.  Enclosed mass  M(<r)
    %     M(<r) = 4*pi * int_0^r rho(r') r'^2 dr'
    %
    %     Work in log(r): dr' = r' d(ln r')
    %     => integrand = 4*pi * rho(r') * r'^3  d(ln r')
    %     Accurate for power-law rho(r) across many decades.
    %     Binney & Tremaine 2008, Eq.(2.22)
    % ------------------------------------------------------------------
    if need_M
        ln_r      = log(r);
        integrand = 4*pi * rho .* r.^3;         % d(ln r) integrand
        M_enc     = cumtrapz(ln_r, integrand);   % M(<r[i])

        % Handle the innermost shell: assume rho ~ const near r=0
        % (solid sphere approximation for the first point)
        M_enc(1)  = (4*pi/3) * rho(1) * r(1)^3;
    end

    % ------------------------------------------------------------------
    % 2.  Circular velocity  V_c(r) = sqrt(G * M(<r) / r)
    %     Binney & Tremaine 2008, Eq.(2.30)
    % ------------------------------------------------------------------
    if any(strcmpi(quantity_out,'V_circ'))
        V_circ = sqrt(G * M_enc ./ r);
    end

    % ------------------------------------------------------------------
    % 3.  Isotropic Jeans: 3-D radial velocity dispersion sigma_r(r)
    %
    %     sigma_r^2(r) = (1/rho(r)) * int_r^inf rho(r') * G*M(r')/r'^2 dr'
    %
    %     Boundary condition: sigma_r -> 0 as r -> r_max (outermost point).
    %     Integrate from outside inward using cumtrapz on reversed array.
    %     In log(r): integrand = rho(r') * G*M(r') / r'  d(ln r')
    %     Binney & Tremaine 2008, Eq.(4.215)
    % ------------------------------------------------------------------

    if need_sigma_r
        integrand_J = rho .* G .* M_enc ./ r;   % [M_sun (km/s)^2 / kpc^3]

        % Incremental trapezoid contributions on log(r) grid
        dln_r  = diff(ln_r);                          % N-1 values
        f_mid  = 0.5 * (integrand_J(1:end-1) + ...
                    integrand_J(2:end));           % N-1 values
        pieces = f_mid .* dln_r;                      % N-1 interval contributions

        % Tail sum: cumJ(i) = int_{r(i)}^{r_max} integrand d(ln r')
        % Last point = 0 (boundary condition sigma_r -> 0 at r_max)
        cumJ = [flipud(cumsum(flipud(pieces))); 0];   % N values

        sigma_r2 = cumJ ./ max(rho, 1e-300);
        sigma_r2 = max(sigma_r2, 0);
        sigma_r  = sqrt(sigma_r2);                    % [km/s]
    end

    % ------------------------------------------------------------------
    % 4.  Abel integrals: Sigma(R) and sigma_los(R)
    %
    %     Sigma(R)     = 2 * int_R^inf rho(r) * r/sqrt(r^2-R^2) dr
    %     sigma_los^2  = (2/Sigma) * int_R^inf rho*sigma_r^2 * r/sqrt(r^2-R^2) dr
    %
    %     For each projected radius R = r(i), integrate over r(j) >= r(i).
    %     The integrand diverges as r -> R^+; the finite grid spacing
    %     provides a natural regularisation (first sub-interval is skipped).
    %
    %     Substitution u = sqrt(r^2 - R^2) removes the singularity:
    %       int_R^inf f(r)*r/sqrt(r^2-R^2) dr = int_0^inf f(sqrt(u^2+R^2)) du
    %     We use this form for each R_i.
    %
    %     Binney & Tremaine 2008, Eq.(2.135) and Eq.(4.225)
    % ------------------------------------------------------------------
    if need_Sigma
        Sigma     = zeros(N, 1);
        if need_sigma_r
            sigma_los2 = zeros(N, 1);
        end

        for i = 1:N
            R_i = r(i);

            % Sub-array r >= R_i (skip r < R_i to avoid singularity)
            idx = r >= R_i;
            r_sub   = r(idx);
            rho_sub = rho(idx);

            if numel(r_sub) < 2
                Sigma(i) = 0;
                if need_sigma_r
                    sigma_los2(i) = 0;
                end
                continue
            end

            % u-substitution: u = sqrt(r^2 - R_i^2), r = sqrt(u^2 + R_i^2)
            % du/dr = r/sqrt(r^2-R_i^2)  =>  integrand in u-space = 2*rho(r)
            u_sub = sqrt(r_sub.^2 - R_i^2);

            % Remove duplicate u=0 point if it appears
            dup = u_sub == 0;
            if sum(dup) > 1
                keep = ~dup;
                keep(find(dup,1)) = true;
                u_sub   = u_sub(keep);
                rho_sub = rho_sub(keep);
                if need_sigma_r
                    sr2_sub = sigma_r2(idx);
                    sr2_sub = sr2_sub(keep);
                end
            else
                if need_sigma_r
                    sr2_sub = sigma_r2(idx);
                end
            end

            % Sigma: Sigma(R) = 2 * int_0^inf rho(sqrt(u^2+R^2)) du
            Sigma(i) = 2 * trapz(u_sub, rho_sub);

            % sigma_los: weighted by rho*sigma_r^2
            if need_sigma_r
                sigma_los2(i) = 2 * trapz(u_sub, rho_sub .* sr2_sub);
            end
        end

        % sigma_los^2 = integral / Sigma
        if need_sigma_r
            sigma_los = sqrt(max(sigma_los2 ./ max(Sigma, 1e-300), 0));
        end
    end

    % ------------------------------------------------------------------
    % 5.  Pack output
    % ------------------------------------------------------------------
    prof.r   = r;
    prof.rho = rho;

    for k = 1:numel(quantity_out)
        switch lower(quantity_out{k})
            case 'm_enc',     prof.M_enc     = M_enc;
            case 'v_circ',    prof.V_circ    = V_circ;
            case 'sigma_r',   prof.sigma_r   = sigma_r;
            case 'sigma_los', prof.sigma_los = sigma_los;
            case 'sigma',     prof.Sigma     = Sigma;
            otherwise
                error('profile_obs: unknown quantity ''%s''.', quantity_out{k});
        end
    end

    % Return Sigma if it was computed (may have been requested as 'Sigma'
    % or triggered by 'sigma_los')
    if need_Sigma && ~isfield(prof, 'Sigma')
        prof.Sigma = Sigma;
    end
end
