function M_hm = halfmode_mass(model, varargin)
% halfmode_mass   Half-mode mass scale for WDM or FDM.
%
% USAGE:
%   M_hm = halfmode_mass('WDM', m_WDM, cosmo)
%   M_hm = halfmode_mass('FDM', m22)
%
% WDM: Schneider et al. (2012), MNRAS 424, 684, Eqs. 5, 8, 9.
%      Uses the Viel et al. (2005) transfer function with mu = 1.12.
%      INPUT:  m_WDM : thermal relic particle mass [keV]
%              cosmo : PHANTOM cosmo struct
%      OUTPUT: M_hm  : half-mode mass [h^-1 M_sun]
%
% FDM: Schive et al. (2016), Eq. (21).
%      INPUT:  m22  : boson mass in units of 1e-22 eV/c^2
%      OUTPUT: M_hm : half-mode mass [M_sun]

    switch upper(model)

        case 'WDM'
            m_WDM = varargin{1};
            cosmo  = varargin{2};
            Om_WDM = cosmo.Omega_m - cosmo.Omega_b;
            h      = cosmo.h;
            % Eq. 5: effective free-streaming length [Mpc/h]
            alpha_fs  = 0.049 * (m_WDM)^(-1.11) ...
                               * (Om_WDM / 0.25)^(0.11) ...
                               * (h / 0.7)^(1.22);
            % Eq. 8: half-mode length [Mpc/h], mu = 1.12
            mu        = 1.12;
            lambda_hm = 2*pi * alpha_fs * (2^(mu/5) - 1)^(-1/(2*mu));
            % Eq. 9: half-mode mass [h^-1 M_sun]
            rho_m = cosmo.rho_m0;
            M_hm  = (4*pi/3) * rho_m * (lambda_hm/2)^3;

        case 'FDM'
            m22  = varargin{1};
            M_hm = 3.8e10 .* m22.^(-4/3);

        otherwise
            error('halfmode_mass: unknown model ''%s''. Use ''WDM'' or ''FDM''.', model);

    end
end