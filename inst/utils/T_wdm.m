function T = T_wdm(k, m_WDM_keV, cosmo, model, varargin)
% transfer_function_WDM   WDM linear matter transfer function T(k)
%
%   T = transfer_function_WDM(k, m_WDM_keV, cosmo)
%   T = transfer_function_WDM(k, m_WDM_keV, cosmo, 'viel2005')
%   T = transfer_function_WDM(k, m_WDM_keV, cosmo, 'bode2001')
%   T = transfer_function_WDM(k, m_WDM_keV, cosmo, 'bode2001', 'g_eff', 100, 'nu', 1.0)
%
%   Returns T(k) such that P_WDM(k) = T(k)^2 * P_CDM(k).
%
% INPUTS:
%   k           : wavenumber array [h/Mpc]
%   m_WDM_keV   : WDM particle mass [keV], thermal relic
%   cosmo       : PHANTOM cosmo struct (needs .Om0, .h fields)
%   model       : (optional) transfer function model, default 'viel2005'
%   varargin    : optional name-value pairs (bode2001 only):
%                   'g_eff'  — effective relativistic d.o.f. at decoupling
%                              default: 1.5  (e.g. 100 for gravitino/gluino)
%                   'nu'     — shape fitting parameter
%                              default: 1.12 (Boltzmann best-fit, Viel+2005 Eq.6)
%                              Viel+2005 use nu=1.2 in their MCMC fits (Sec.IVA)
%
% OUTPUT:
%   T           : transfer function, same size as k, values in (0,1]
%
% Models:
%   'viel2005'  — Viel, Lesgourgues, Haehnelt, Matarrese & Riotto (2005)
%                 Phys. Rev. D 71, 063534, arXiv:astro-ph/0501562, Eqs. 6-7.
%                 Adopted by Schneider+2012 (their Eqs. 4-5).
%                 Standard choice for thermal relic WDM. nu and alpha
%                 are fixed from Boltzmann simulations; not user-adjustable
%                 in this model (see 'bode2001' if you need free nu).
%
%   'bode2001'  — Bode, Ostriker & Turok (2001), ApJ 556, 93
%                 arXiv:astro-ph/0010389, Eq. 16.
%                 Supports user-supplied g_eff and nu.

    if nargin < 4 || isempty(model)
        model = 'viel';
    end

    % --- Parse optional name-value pairs (used by bode2001) --------------
    g_eff = 1.5;   % default: thermal relic
    nu    = 1.2;   % default: Boltzmann best-fit (Viel+2005 Eq.6)

    for i = 1:2:length(varargin)
        switch lower(varargin{i})
            case 'g_eff'
                g_eff = varargin{i+1};
            case 'nu'
                nu = varargin{i+1};
            otherwise
                error('Unknown parameter: %s', varargin{i});
        end
    end

    switch lower(model)

        case {'viel', 'viel05'}
            % Viel, Lesgourgues, Haehnelt, Matarrese & Riotto (2005)
            % Phys. Rev. D 71, 063534, arXiv:astro-ph/0501562
            % Equations 6 and 7 (Section II, "Pure warm dark matter models")
            %
            % T(k) = [1 + (alpha*k)^(2*mu)]^(-5/mu),   mu = 1.12
            %
            % alpha = 0.049 * (m_WDM/keV)^{-1.11}
            %               * (Omega_WDM/0.25)^{0.11}
            %               * (h/0.7)^{1.22}   [h^{-1} Mpc]
            %
            % nu = 1.12 is the best fit from Boltzmann code for k < 5 h/Mpc.
            % Note: Viel+2005 adopt nu = 1.2 in their MCMC fits (their Sec.IVA);
            % nu is fixed here per the published formula. Use 'bode2001' with
            % the 'nu' option if you need a free shape parameter.
            mu    = 1.12;
            alpha = 0.049 * (m_WDM_keV)^(-1.11) ...
                          * (cosmo.Omega_m / 0.25)^(0.11) ...
                          * (cosmo.h   / 0.7 )^(1.22);   % [h^{-1} Mpc]
            T = (1 + (alpha .* k).^(2*mu)).^(-5/mu);

        case {'bode', 'bode01'}
            % Bode, Ostriker & Turok (2001), ApJ 556, 93
            % arXiv:astro-ph/0010389, Eq. 16
            %
            % T(k) = [1 + (alpha*k)^(2*nu)]^(-5/nu)
            %
            % alpha = 0.048 * (Omega_WDM/0.4)^{0.15}
            %               * (h/0.65)^{1.3}
            %               * (keV/m_WDM)^{1.15}
            %               * (1.5/g_eff)^{0.29}   [h^{-1} Mpc]
            %
            % g_eff : effective relativistic d.o.f. at decoupling
            %   1.5   — default thermal relic (Bode+2001 fiducial)
            %   ~100  — gravitino or particles decoupling near QCD transition
            %   1.5 to 106.75 depending on species and decoupling temperature
            %   Pass via: transfer_function_WDM(..., 'g_eff', 100)
            %
            % nu : shape fitting parameter, default 1.12
            %   Bode+2001 use nu=1.12 in their Eq.16.
            %   Viel+2005 use nu=1.2 in MCMC analysis (their Sec.IVA).
            %   Pass via: transfer_function_WDM(..., 'nu', 1.2)
            alpha = 0.048 * (cosmo.Omega_m / 0.4 )^(0.15) ...
                          * (cosmo.h   / 0.65)^(1.3)  ...
                          * (1.0 / m_WDM_keV)^(1.15)  ...
                          * (1.5 / g_eff)^(0.29);        % [h^{-1} Mpc]
            T = (1 + (alpha .* k).^(2*nu)).^(-5/nu);

        otherwise
            error('transfer_function_WDM: unknown model "%s". Choose ''viel2005'' or ''bode2001''.', model);
    end
end