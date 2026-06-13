function P = Diemer15_Table(statistic)
% Diemer15_Table  Parameters for the Diemer & Kravtsov (2015) concentration model
%
%   P = Diemer15_Table(statistic)
%
%   Returns best-fit parameters from Table 3 of Diemer & Kravtsov (2015),
%   for the double power-law c-nu model (Eqs. 9-10).
%
%   STATISTIC options (case-insensitive)
%   -------------------------------------
%   'median'   Median concentration (default)
%   'mean'     Mean concentration
%
%   Fields of output struct P
%   -------------------------
%   P.kappa    kappa  — scale factor for Lagrangian radius in n_eff evaluation
%   P.phi0     phi_0  — c_min intercept
%   P.phi1     phi_1  — c_min slope w.r.t. n_eff
%   P.eta0     eta_0  — nu_min intercept
%   P.eta1     eta_1  — nu_min slope w.r.t. n_eff
%   P.alpha    alpha  — low-nu (falling) power-law slope
%   P.beta     beta   — high-nu (rising) power-law slope
%
%   Model form (Eqs. 9-10 of DK15):
%     c_min  = phi0 + phi1 * n_eff        (concentration floor)
%     nu_min = eta0 + eta1 * n_eff        (peak height at floor)
%     c200c  = (c_min/2) * [ (nu/nu_min)^(-alpha) + (nu/nu_min)^beta ]
%
%   Reference: Diemer & Kravtsov 2015, ApJ 799, 108, Table 3
%   Corrected: Diemer & Joyce 2019,  ApJ 871, 168, Table 3            

% Initialize parameters based on the selected statistic
P = struct();

switch lower(statistic)

    case 'median'
        % Table 3, median c-nu relation
        P.kappa = 1.00;
        P.phi0  = 6.58;
        P.phi1  = 1.27;
        P.eta0  = 7.28;
        P.eta1  = 1.56;
        P.alpha = 1.08;
        P.beta  = 1.77;

    case 'mean'
        % Table 3, mean c-nu relation
        P.kappa = 1.00;
        P.phi0  = 6.66;
        P.phi1  = 1.37;
        P.eta0  = 5.41;
        P.eta1  = 1.06;
        P.alpha = 1.22;
        P.beta  = 1.22;

    otherwise
        error('Diemer15_Table: unknown statistic "%s". Valid: ''median'', ''mean''.', statistic);
end
end