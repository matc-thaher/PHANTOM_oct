function c_sup = suppression_factor(M, M_half, model, varargin)
% SUPPRESSION_FACTOR  Concentration suppression relative to CDM for FDM/WDM halos
%
% Usage:
%   c_sup = suppression_factor(M, M_half)                   % default: Laroche (2022)
%   c_sup = suppression_factor(M, M_half, 'laroche22')    % Laroche et al. (2022)
%   c_sup = suppression_factor(M, M_half, 'dentler22')    % Dentler et al. (2022)
%   c_sup = suppression_factor(M, M_half, 'schneider12')  % Schneider et al. (2012)
%   c_sup = suppression_factor(M, [], 'mah_wdm', z_obs, m_WDM_keV, cosmo)
%
% Inputs:
%   M       - Halo mass [M_sun], scalar or array
%   M_half  - Half-mode mass scale [M_sun]
%             FDM: M_half = 3.8e10 * m22^(-4/3) [M_sun]
%             WDM: M_half = M_hm from Eq. 9 of Schneider+2012
%   model   - (optional) suppression model: 'laroche2022' (default),
%             'dentler2022', or 'schneider2012'
%
% Outputs:
%   c_sup   - Dimensionless suppression factor (0 < c_sup <= 1)
%
% References:
%   Laroche et al. (2022) - Eq. for concentration suppression:
%       F(x) = (1 + a*x^b)^c,  x = M/M_half
%       parameters a=5.496, b=-1.648, c=-0.417
%
%   Dentler et al. (2022) - Two-factor suppression (Kawai et al. 2024, Eqs. 29-30):
%       Delta^FDM = (1 + M0/(f_coll*M))^(-gamma0) * (1 + gamma1*M0/M)^(-gamma2)
%       where M0 = 1.6e10 * m22^(-4/3)  [NOTE: M0 != M_half]
%              M0 = (1.6/3.8) * M_half  = M_half / 2.375
%       f_coll = 0.01, gamma0 = d ln(c_vir^B)/d ln(M_h) |_{Mh=4*M0} ~ 0.06,
%       gamma1 = 15, gamma2 = 0.3
%
%   Schneider et al. (2012), MNRAS 424, 684 - Eq. 39:
%       c_WDM/c_CDM = (1 + gamma1*(M_half/M))^(-gamma2)
%       gamma1 = 15,  gamma2 = 0.3
%       Calibrated on WDM N-body simulations, mWDM = 0.25-1.25 keV.
%       M_half is the WDM half-mode mass scale (their M_hm).

if nargin < 3 || isempty(model)
    model = 'laroche22';
end

switch lower(model)
    case 'laroche22'
        % Laroche et al. (2022) power-law suppression
        % c_vir(FDM)/c_vir(CDM) = F(M/M_half) = (1 + a*(M/M_half)^b)^c
        a    =  5.496;
        b    = -1.648;
        cnst = -0.417;
        x    = M ./ M_half;
        c_sup = (1 + a .* x.^b).^cnst;

    case 'dentler22'
        % Dentler et al. (2022) two-factor FDM suppression
        % As expressed in Kawai et al. (2024), Eqs. 29-30:
        %
        %   Delta^FDM = (1 + M0/(f_coll*M))^(-gamma0)
        %             * (1 + gamma1*(M0/M))^(-gamma2)
        %
        % M0 = 1.6e10 * m22^(-4/3)  [NOT the half-mode mass]
        %    = (1.6/3.8) * M_half
        %
        % gamma0 = d ln(c_vir^B)/d ln(M_h) evaluated at M_h = 4*M0
        %        ~ 0.06  (CDM low-mass slope from Kawai Eqs. 26-27)
        % gamma1 = 15,  gamma2 = 0.3,  f_coll = 0.01
        % M0 is NOT M_half; it differs by the ratio of prefactors (1.6 vs 3.8)
        M0     = M_half .* (1.6 / 3.8);   % = 1.6e10 * m22^(-4/3)
        f_coll = 0.01;
        gamma1 = 15.0;
        gamma2 = 0.3;
        % gamma0 = local log-slope of CDM c_vir-M relation at M_h = 4*M0
        % Pivot mass: 10^11 h^-1 Msun ~ 1.43e11 Msun  (with h = 0.7)
        M_pivot = 1e11 / 0.7;   % [M_sun]
        if 4 * M0 < M_pivot
            gamma0 = 0.06;      % low-mass CDM slope  (Kawai Eq. 26)
        else
            gamma0 = 0.12;      % high-mass CDM slope (Kawai Eq. 27)
        end
        % Two-factor suppression (Dentler Eq. 33 × Eq. 35 combined)
        factor1 = (1 + M0 ./ (f_coll .* M)).^(-gamma0);
        factor2 = (1 + gamma1 .* M0 ./ M).^(-gamma2);
        c_sup   = factor1 .* factor2;

    case 'schneider12'
        % Schneider, Smith, Maccio & Moore (2012), MNRAS 424, 684
        % Concentration suppression for WDM haloes, Eq. 39:
        %
        %   c_WDM/c_CDM = (1 + gamma1*(M_hm/M))^(-gamma2)
        %
        % Parameters from least-squares fit to N-body simulations:
        %   gamma1 = 15,  gamma2 = 0.3
        % Valid for mWDM = 0.25-1.25 keV, accurate to ~10 per cent.
        % M_half here is the WDM half-mode mass M_hm (their Eq. 9).
        gamma1 = 15.0;
        gamma2 = 0.3;
        c_sup  = (1 + gamma1 .* (M_half ./ M)).^(-gamma2);

   case 'mah_wdm'
        % Ludlow et al. (2016) physically motivated WDM concentration
        % suppression. Calls Ludlow16.m twice — once with CDM sigmaM and
        % once with a WDM-modified sigmaM — and returns c_WDM / c_CDM.
        %
        % Transfer function delegated entirely to transfer_function_WDM.m.
        % Default model is 'viel2005'. Pass 'bode2001' and any optional
        % name-value pairs (g_eff, nu) via tf_args if needed.
        %
        % CALL SIGNATURES:
        %   suppression_factor(M, [], 'ludlow2016_wdm', z_obs, m_WDM_keV, cosmo)
        %   suppression_factor(M, [], 'ludlow2016_wdm', z_obs, m_WDM_keV, cosmo, 'bode2001')
        %   suppression_factor(M, [], 'ludlow2016_wdm', z_obs, m_WDM_keV, cosmo, 'bode2001', 'g_eff', 100)
        %   suppression_factor(M, [], 'ludlow2016_wdm', z_obs, m_WDM_keV, cosmo, 'bode2001', 'nu', 1.2)
        %   suppression_factor(M, [], 'ludlow2016_wdm', z_obs, m_WDM_keV, cosmo, 'bode2001', 'g_eff', 100, 'nu', 1.2)
        %
        % NOTE: M_half (second argument) is unused here. Pass [] or any value.

        % --- Validate minimum required inputs ----------------------------
        if numel(varargin) < 3
            error(['suppression_factor:ludlow2016_wdm — requires at least ' ...
                   'z_obs, m_WDM_keV, cosmo after the model string.']);
        end

        z_obs     = varargin{1};   % observation redshift
        m_WDM_keV = varargin{2};   % WDM particle mass [keV]
        cosmo_cdm = varargin{3};   % PHANTOM cosmo struct (CDM)

        % --- Remaining varargin: tf_model + any name-value pairs ---------
        % Everything after cosmo is forwarded to transfer_function_WDM.m.
        % e.g. {'bode2001'} or {'bode2001','g_eff',100,'nu',1.2}
        tf_args = varargin(4:end);

        % --- Build WDM cosmo: only sigmaM is modified --------------------
        rho_m0 = cosmo_cdm.Omega_m * cosmo_cdm.rho_crit0;   % [M_sun/h / (Mpc/h)^3]

        cosmo_wdm        = cosmo_cdm;
        cosmo_wdm.sigmaM = @(M_in, z_in) ...
            sigmaM_WDM(M_in, z_in, cosmo_cdm, m_WDM_keV, tf_args{:});

        % --- Solve c for CDM and WDM -------------------------------------
        c_cdm = Ludlow16(M, z_obs, cosmo_cdm);
        c_wdm = Ludlow16(M, z_obs, cosmo_wdm);

        % Guard: halos below the free-streaming scale may not converge
        c_wdm(isnan(c_wdm)) = 0;
        c_cdm(c_cdm == 0)   = NaN;

        c_sup = c_wdm ./ c_cdm;
        c_sup = min(c_sup, 1.0);   % suppression only; ratio never exceeds 1

    otherwise
        error('suppression_factor: unknown model "%s". Choose ''laroche'', ''dentler'', ''schneider'', or ''ludlow''.', model);
end
end