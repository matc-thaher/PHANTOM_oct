function f = multiplicity_FernandezGarcia26(sigma, M, z, mdef, cosmo)
% multiplicity_FernandezGarcia26   GPS+ halo multiplicity function
%
%   f = multiplicity_FernandezGarcia26(sigma, M, z, mdef)
%
%   Implements the generalised Press-Schechter + triaxial collapse (GPS+)
%   model of Fernandez-Garcia et al. (2026), arXiv:2512.05847.
%
%   The model predicts d(ln n)/d(ln M) via a numerical integral over
%   triaxial collapse geometry (the V integral, their Eq. X), then
%   converts to f(sigma) using a finite difference in M.
%
%   Unlike most multiplicity functions, this model is internally mass-
%   based, not purely sigma-based. Both sigma and M are required.
%
% INPUTS:
%   sigma  : rms linear fluctuation at mass M and redshift z (array)
%   M      : halo mass [M_sun/h], same size as sigma
%   z      : redshift (scalar)
%   mdef   : mass definition string — '200m' or 'vir' only
%            (GPS+ calibrated on Uchuu suite, M_200m and M_vir)
%   cosmo  : PHANTOM cosmo struct with fields:
%              cosmo.Om0, cosmo.rho_crit0
%              cosmo.sigmaR(R, z)   — sigma at Lagrangian radius R
%
% OUTPUT:
%   f      : multiplicity function f(sigma), same size as sigma
%
% DEPENDENCIES:
%   cosmo struct must expose:
%     cosmo.sigmaR(R, z)     — sigma(R, z) for Lagrangian radius R
%     cosmo.rho_m0           — mean matter density today [M_sun/h / (Mpc/h)^3]
%
% Reference:
%   Fernandez-Garcia, Betancort-Rijo, Prada et al. (2026)
%   arXiv:2512.05847

    if ~ismember(mdef, {'200m', 'vir'})
        error('multiplicity_FernandezGarcia26: mdef must be ''200m'' or ''vir''. Got: %s', mdef);
    end

    sigma = sigma(:);
    M     = M(:);

    % --- Model constants -------------------------------------------------
    aa = 1.089;
    bb = 0.652;
    DD = 1.0;
    EE = 0.17;
    FF = 0.087;

    % --- Mass-dependent b(M) and c(M) via log-polynomial fit -------------
    % Calibration nodes from Fernandez-Garcia+2026 (Table 1 / Fig. 2)
    M_nodes_b = [1e16, 1e15, 1e14, 6.5e10, 1e10, 1e9, 1e8, 1e7, 1e6];
    b_nodes   = [0.5259, 0.415, 0.328, 0.1764, 0.1552, 0.1308, 0.1179, 0.1045, 0.094];
    coeffs_b  = polyfit(log10(M_nodes_b), log10(b_nodes), 4);

    M_nodes_c = [3e15, 3e14, 3e13, 3e12, 3e11, 3e10, 3e9, 3e8, 3e7, 3e6, ...
                 1e10, 1e9, 1e8, 1e7, 1e6];
    c_nodes   = [0.613, 0.474, 0.373, 0.301, 0.249, 0.209, 0.1794, 0.1560, 0.1355, 0.1223, ...
                 0.1942, 0.168, 0.1466, 0.1298, 0.1161];
    coeffs_c  = polyfit(log10(M_nodes_c), log10(c_nodes), 4);

    b_M = 10.^polyval(coeffs_b, log10(M));
    c_M = 10.^polyval(coeffs_c, log10(M));

    % --- Effective delta_c with sigma-dependent correction ---------------
    % x = sigma / delta_c_EdS;  delta_c_EdS = 1.676
    x      = sigma ./ 1.676;
    U2     = (-0.00221 .* x.^3 + 0.03835 .* x.^2 + 0.17810 .* x - 0.01507).^2;
    sig_mod = sqrt(sigma.^2 + U2);

    term1   = (1 + 0.845 .* x - 0.04 .* x.^2 + 0.0025 .* x.^3).^bb;
    term2   = aa .* 1.365 .* (1 + EE .* b_M - FF .* b_M.^2).^DD;
    delta_c = term1 .* term2;
    cte     = delta_c ./ (sqrt(2) .* sig_mod);

    % --- Triaxial collapse integral V(c_M, cte) --------------------------
    % V = 3 * integral_0^1 erfc(cte * sqrt((1-exp(-c*xi^2))/(1+exp(-c*xi^2)))) * xi^2 dxi
    xi      = linspace(0, 1, 1000);

    % Vectorised: rows = mass points, columns = xi quadrature points
    cte_mat = cte(:);          % [N x 1]
    c_mat   = c_M(:);          % [N x 1]
    xi_mat  = xi(:)';          % [1 x 1000]

    arg      = (1 - exp(-c_mat .* xi_mat.^2)) ./ (1 + exp(-c_mat .* xi_mat.^2));
    integrand = erfc(cte_mat .* sqrt(arg)) .* xi_mat.^2;
    V = 3 .* trapz(xi, integrand, 2);   % [N x 1]

    % --- F(M) = erfc(0.98 * cte) / V ------------------------------------
    F_M = erfc(0.98 .* cte) ./ V;

    % --- Finite difference: dF/dM -> dn/dlnM -> f(sigma) ----------------
    % d(ln n)/d(ln M) evaluated via forward difference at M*(1+s)
    % This mirrors the colossus implementation exactly.
    s     = 0.01;
    M_p   = M .* (1 + s);

    % Lagrangian radius R = (3M / 4*pi*rho_m0)^{1/3}
    % sigma at M_p: caller must supply a cosmo struct or this uses a handle
    % passed at construction. Here we expose cosmo as a persistent or passed arg.
    % PHANTOM convention: cosmo.sigmaM(M, z) already converts M->R internally.
    rho_m0  = cosmo.rho_m0;                    % [M_sun/h / (Mpc/h)^3]
    R_p     = (3 .* M_p ./ (4 * pi * rho_m0)).^(1/3);
    sig_p   = cosmo.sigmaR(R_p, z);
    sig_p   = sig_p(:);

    x_p     = sig_p ./ 1.676;
    U2_p    = (-0.00221 .* x_p.^3 + 0.03835 .* x_p.^2 + 0.17810 .* x_p - 0.01507).^2;
    sig_mod_p = sqrt(sig_p.^2 + U2_p);
    b_p     = 10.^polyval(coeffs_b, log10(M_p));
    c_p     = 10.^polyval(coeffs_c, log10(M_p));
    term1_p = (1 + 0.845 .* x_p - 0.04 .* x_p.^2 + 0.0025 .* x_p.^3).^bb;
    term2_p = aa .* 1.365 .* (1 + EE .* b_p - FF .* b_p.^2).^DD;
    dc_p    = term1_p .* term2_p;
    cte_p   = dc_p ./ (sqrt(2) .* sig_mod_p);

    cte_p_mat = cte_p(:);
    c_p_mat   = c_p(:);
    arg_p     = (1 - exp(-c_p_mat .* xi_mat.^2)) ./ (1 + exp(-c_p_mat .* xi_mat.^2));
    integ_p   = erfc(cte_p_mat .* sqrt(arg_p)) .* xi_mat.^2;
    V_p       = 3 .* trapz(xi, integ_p, 2);
    F_p       = erfc(0.98 .* cte_p) ./ V_p;

    dF_dM    = (F_M - F_p) ./ (s .* M);
    dn_dlnM  = dF_dM .* rho_m0 ./ (1 + s/2);

    % --- Convert dn/dlnM -> f(sigma) -------------------------------------
    % f(sigma) = (M / rho_m0) * |d(ln sigma^-1)/d(ln M)|^{-1} * dn/dlnM
    % Use central finite difference for d(ln sigma)/d(ln M)
    dlns_dlnM = dlnsigma_dlnM(M, z, cosmo);   % PHANTOM utility
    f = dn_dlnM .* M ./ rho_m0 ./ abs(dlns_dlnM');
    f = f(:);
    f_smooth = smoothdata(f, 'sgolay', 15);

    % --- Size diagnostic (remove after fixing) ---
    fprintf('sigma:      %d x %d\n', size(sigma));
    fprintf('M:          %d x %d\n', size(M));
    fprintf('M_p:        %d x %d\n', size(M_p));
    fprintf('sig_p:      %d x %d\n', size(sig_p));
    fprintf('cte_p:      %d x %d\n', size(cte_p));
    fprintf('V_p:        %d x %d\n', size(V_p));
    fprintf('F_p:        %d x %d\n', size(F_p));
    fprintf('dF_dM:      %d x %d\n', size(dF_dM));
    fprintf('dn_dlnM:    %d x %d\n', size(dn_dlnM));
    fprintf('dlns_dlnM:  %d x %d\n', size(dlns_dlnM));
    fprintf('f:          %d x %d\n', size(f));
    % ---------------------------------------------

    f = reshape(f_smooth, size(sigma));
end

% function f = multiplicity_FernandezGarcia26(sigma, M, z, mdef, cosmo)
% % multiplicity_FernandezGarcia26   GPS+ halo multiplicity function
% %
% %   Implements Fernandez-Garcia et al. (2026), arXiv:2512.05847.
% %
% % INPUTS:
% %   sigma  : rms linear fluctuation at mass M (column array)
% %   M      : halo mass [M_sun/h], same size as sigma
% %   z      : redshift (scalar)
% %   mdef   : mass definition — '200m' or 'vir'
% %   cosmo  : PHANTOM cosmo struct
% %
% % OUTPUT:
% %   f      : multiplicity function f(sigma)
% 
%     if ~ismember(mdef, {'200m', 'vir'})
%         error('multiplicity_FernandezGarcia26: mdef must be ''200m'' or ''vir''. Got: %s', mdef);
%     end
% 
%     % Force column vectors
%     sigma = sigma(:);
%     M     = M(:);
% 
%     % --- Model constants -------------------------------------------------
%     aa = 1.089;
%     bb = 0.652;
%     DD = 1.0;
%     EE = 0.17;
%     FF = 0.087;
% 
%     % --- Polynomial fit coefficients for b(M) and c(M) ------------------
%     M_nodes_b = [1e16, 1e15, 1e14, 6.5e10, 1e10, 1e9, 1e8, 1e7, 1e6];
%     b_nodes   = [0.5259, 0.415, 0.328, 0.1764, 0.1552, 0.1308, 0.1179, 0.1045, 0.094];
%     coeffs_b  = polyfit(log10(M_nodes_b), log10(b_nodes), 4);
% 
%     M_nodes_c = [3e15, 3e14, 3e13, 3e12, 3e11, 3e10, 3e9, 3e8, 3e7, 3e6, ...
%                  1e10, 1e9, 1e8, 1e7, 1e6];
%     c_nodes   = [0.613, 0.474, 0.373, 0.301, 0.249, 0.209, 0.1794, 0.1560, ...
%                  0.1355, 0.1223, 0.1942, 0.168, 0.1466, 0.1298, 0.1161];
%     coeffs_c  = polyfit(log10(M_nodes_c), log10(c_nodes), 4);
% 
%     % --- Quadrature grid (shared by F_M and F_plus) ----------------------
%     xi     = linspace(0, 1, 1000);
%     xi_row = reshape(xi, 1, []);   % [1 x 1000]
% 
%     % --- Mean matter density ---------------------------------------------
%     rho_m0 = cosmo.rho_m0;   % [M_sun/h / (Mpc/h)^3]
% 
%     % =====================================================================
%     % Nested helper: computes F(M, sigma_in) given mass and sigma arrays
%     % =====================================================================
%     function F = compute_F(M_in, sig_in)
%         M_in   = M_in(:);
%         sig_in = sig_in(:);
% 
%         b_in = 10.^polyval(coeffs_b, log10(M_in));
%         c_in = 10.^polyval(coeffs_c, log10(M_in));
% 
%         x        = sig_in ./ 1.676;
%         U2       = (-0.00221.*x.^3 + 0.03835.*x.^2 + 0.17810.*x - 0.01507).^2;
%         sig_mod  = sqrt(sig_in.^2 + U2);
%         term1    = (1 + 0.845.*x - 0.04.*x.^2 + 0.0025.*x.^3).^bb;
%         term2    = aa .* 1.365 .* (1 + EE.*b_in - FF.*b_in.^2).^DD;
%         delta_c  = term1 .* term2;
%         cte      = delta_c ./ (sqrt(2) .* sig_mod);
% 
%         % Vectorised integral [N x 1000]
%         cte_col = reshape(cte,  [], 1);
%         c_col   = reshape(c_in, [], 1);
% 
%         arg       = (1 - exp(-c_col .* xi_row.^2)) ./ ...
%                     (1 + exp(-c_col .* xi_row.^2));
%         integrand = erfc(cte_col .* sqrt(arg)) .* xi_row.^2;
%         V         = 3 .* trapz(xi, integrand, 2);   % [N x 1]
% 
%         F = erfc(0.98 .* cte) ./ V;
%         F = F(:);
%     end
% 
%     % =====================================================================
%     % Main computation
%     % =====================================================================
% 
%     % --- F at M ----------------------------------------------------------
      % F_M = compute_F(M, sigma);
%     % --- F at M*(1+s) ----------------------------------------------------
%     s      = 0.01;
%     M_plus = M .* (1 + s);
%     R_plus = (3 .* M_plus ./ (4 * pi * rho_m0)).^(1/3);
%     R_plus = R_plus(:);
% 
%     sig_plus = cosmo.sigmaR(R_plus, z);
%     sig_plus = sig_plus(:);
% 
%     F_plus = compute_F(M_plus, sig_plus);
% 
%     % --- dF/dM and dn/dlnM ----------------------------------------------
%     dF_dM   = (F_M - F_plus) ./ (s .* M);
%     dn_dlnM = dF_dM .* rho_m0 ./ (1 + s/2);
% 
%     % --- Convert dn/dlnM -> f(sigma) ------------------------------------
%     % Derivative from same sig_plus — consistent with F_plus evaluation
%     dlns_dlnM = log(sig_plus ./ sigma) ./ log(1 + s);
% 
%     f = dn_dlnM .* M ./ rho_m0 ./ abs(dlns_dlnM);
%     f = f(:);
% 
% end
% 
% 
