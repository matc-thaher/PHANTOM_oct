function P = Klypin16_Table(formula, cosmo_name, mdef)
% Klypin16_Table  Parameters for the Klypin et al. (2016) concentration models
%
%   P = Klypin16_Table(formula, cosmo_name, mdef)
%
%   Two fitting formulae from Appendix A of Klypin+2016, MNRAS 457, 4340:
%
%   formula='cM'  — eq. (A1), Tables A1-A4:
%     c(M,z) = C0(z) * (M/1e12)^(-gamma(z)) * [1 + (M/M0(z))^0.4]
%     Available for: planck13+200c, planck13+vir, bolshoi+200c, bolshoi+vir
%
%   formula='cnu' — eq. (A5), Tables A5-A6:
%     c(sigma) = b0(z) * [1 + 7.37*(sigma/a0(z))^0.75]
%                       * [1 + 0.14*(sigma/a0(z))^(-2)]
%     Available for: planck13+200c, planck13+vir ONLY
%     (Bolshoi cnu tables are NOT provided in the paper)
%
%   INPUTS
%   formula    : 'cM' or 'cnu'
%   cosmo_name : 'planck13' or 'bolshoi'
%                NOTE: 'bolshoi' is only valid for formula='cM'
%   mdef       : '200c' or 'vir'
%
%   OUTPUT struct P
%   For 'cM':   P.z_bins, P.C0, P.gamma, P.M0  (M0 in Msun/h)
%   For 'cnu':  P.z_bins, P.a0, P.b0
%
%   Reference: Klypin et al. 2016, MNRAS 457, 4340, Appendix A Tables A1-A6

formula    = lower(formula);
cosmo_name = lower(cosmo_name);
mdef       = lower(mdef);

if ~ismember(formula, {'cm','cnu'})
    error('Klypin16_Table: unknown formula "%s". Valid: ''cM'', ''cnu''.', formula);
end
if ~ismember(cosmo_name, {'planck13','bolshoi'})
    error('Klypin16_Table: unknown cosmology "%s". Valid: ''planck13'', ''bolshoi''.', cosmo_name);
end
if ~ismember(mdef, {'200c','vir'})
    error('Klypin16_Table: unknown mdef "%s". Valid: ''200c'', ''vir''.', mdef);
end
if strcmp(formula,'cnu') && strcmp(cosmo_name,'bolshoi')
    error(['Klypin16_Table: the c-nu formula (eq. A5) is only fit to the planck13 ' ...
           'cosmology in Klypin+2016. Bolshoi c-nu tables do not exist in the paper. ' ...
           'Use formula=''cM'' for bolshoi, or switch to cosmology ''planck13''.']);
end

% ========================================================================
% c-M relation  eq.(A1):  c = C0*(M/1e12)^(-gamma) * [1+(M/M0)^0.4]
% M0 stored in Msun/h  (paper gives M0 in units of 1e12 Msun/h)
% ========================================================================
if strcmp(formula, 'cm')

    if strcmp(cosmo_name,'planck13') && strcmp(mdef,'200c')
        % Table A1 — Planck13, M200c, all halos
        P.z_bins = [0.00, 0.35, 0.50, 1.00, 1.44, 2.15, 2.50, 2.90, 4.10, 5.40];
        P.C0     = [7.40, 6.25, 5.65, 4.30, 3.53, 2.70, 2.42, 2.20, 1.92, 1.65];
        P.gamma  = [0.120, 0.117, 0.115, 0.110, 0.095, 0.085, 0.080, 0.080, 0.080, 0.080];
        P.M0     = [5.5e5, 1.0e5, 2.0e4, 9.0e2, 3.0e2, 4.2e1, 1.7e1, 8.5e0, 2.0e0, 3.0e-1] .* 1e12;

    elseif strcmp(cosmo_name,'planck13') && strcmp(mdef,'vir')
        % Table A3 — Planck13, Mvir, all halos
        P.z_bins = [0.00, 0.35, 0.50, 1.00, 1.44, 2.15, 2.50, 2.90, 4.10, 5.40];
        P.C0     = [9.75, 7.25, 6.50, 4.75, 3.80, 3.00, 2.65, 2.42, 2.10, 1.86];
        P.gamma  = [0.110, 0.107, 0.105, 0.100, 0.095, 0.085, 0.080, 0.080, 0.080, 0.080];
        P.M0     = [5.0e5, 2.2e4, 1.0e4, 1.0e3, 2.1e2, 4.3e1, 1.8e1, 9.0e0, 1.9e0, 4.2e-1] .* 1e12;

    elseif strcmp(cosmo_name,'bolshoi') && strcmp(mdef,'200c')
        % Table A5 — Bolshoi, M200c, all halos
        P.z_bins = [0.00, 0.50, 1.00, 1.44, 2.15, 2.50, 2.90, 4.10];
        P.C0     = [6.60, 5.25, 3.85, 3.00, 2.10, 1.80, 1.60, 1.40];
        P.gamma  = [0.110, 0.105, 0.103, 0.097, 0.095, 0.095, 0.095, 0.095];
        P.M0     = [2.0e6, 6.0e4, 8.0e2, 1.1e2, 1.3e1, 6.0e0, 3.0e0, 1.0e0] .* 1e12;

    elseif strcmp(cosmo_name,'bolshoi') && strcmp(mdef,'vir')
        % Table A6 — Bolshoi, Mvir, all halos
        P.z_bins = [0.00, 0.50, 1.00, 1.44, 2.15, 2.50, 2.90, 4.10];
        P.C0     = [9.00, 6.00, 4.30, 3.30, 2.30, 2.10, 1.85, 1.70];
        P.gamma  = [0.100, 0.100, 0.100, 0.100, 0.095, 0.095, 0.095, 0.095];
        P.M0     = [2.0e6, 7.0e3, 5.5e2, 9.0e1, 1.1e1, 6.0e0, 2.5e0, 1.0e0] .* 1e12;
    end

% ========================================================================
% c-nu relation  eq.(A5) — planck13 ONLY
%   c = b0 * [1 + 7.37*(sigma/a0)^0.75] * [1 + 0.14*(sigma/a0)^(-2)]
% ========================================================================
elseif strcmp(formula, 'cnu')

    if strcmp(mdef,'200c')
        % Table A7 — Planck13, M200c
        P.z_bins = [0.00, 0.38, 0.50, 1.00, 1.44, 2.50, 2.89, 5.41];
        P.a0     = [0.40, 0.65, 0.82, 1.08, 1.23, 1.60, 1.68, 1.70];
        P.b0     = [0.278, 0.375, 0.411, 0.436, 0.426, 0.375, 0.360, 0.351];

    elseif strcmp(mdef,'vir')
        % Table A8 — Planck13, Mvir
        P.z_bins = [0.00, 0.38, 0.50, 1.00, 1.44, 2.50, 5.50];
        P.a0     = [0.75, 0.90, 0.97, 1.12, 1.28, 1.52, 1.62];
        P.b0     = [0.567, 0.541, 0.529, 0.496, 0.474, 0.421, 0.393];
    end

end
end