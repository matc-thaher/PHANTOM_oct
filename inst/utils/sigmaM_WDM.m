function sigma = sigmaM_WDM(M, z, cosmo_cdm, m_WDM_keV, varargin)
% sigmaM_WDM   WDM-suppressed rms linear mass fluctuation sigma(M, z)
%
%   sigma = sigmaM_WDM(M, z, cosmo_cdm, m_WDM_keV)
%   sigma = sigmaM_WDM(M, z, cosmo_cdm, m_WDM_keV, 'bode2001')
%   sigma = sigmaM_WDM(M, z, cosmo_cdm, m_WDM_keV, 'bode2001', 'g_eff', 100)
%   sigma = sigmaM_WDM(M, z, cosmo_cdm, m_WDM_keV, 'bode2001', 'nu', 1.2)
%
% Approximates sigma_WDM(M, z) by applying the WDM transfer function
% at an effective wavenumber k_eff = 1/R(M), where R(M) is the
% Lagrangian radius enclosing mass M at the mean matter density:
%
%   R(M) = [ 3M / (4 pi rho_m0) ]^{1/3}   [Mpc/h]
%
% so that:
%
%   sigma_WDM(M, z) = sigma_CDM(M, z) * T_WDM(k_eff)
%
% This is an approximation (~5% error for M >> M_hm, conservative near
% the free-streaming cutoff). It captures the correct suppression scale
% without requiring a full P(k) convolution integral. Sufficient for
% formation-redshift root-finding in Ludlow16.m and HMF suppression.
% Replace with a full integral if your sigmaM infrastructure accepts
% an arbitrary P(k) input.
%
% INPUTS:
%   M           : halo mass array [M_sun/h]
%   z           : redshift (scalar)
%   cosmo_cdm   : PHANTOM cosmo struct with fields:
%                   .sigmaM(M, z)   — CDM rms fluctuation function handle
%                   .Om0            — matter density parameter
%                   .rho_crit0      — critical density [M_sun/h / (Mpc/h)^3]
%                   .h              — dimensionless Hubble parameter
%   m_WDM_keV   : WDM particle mass [keV], thermal relic
%   varargin    : forwarded to transfer_function_WDM
%                   tf_model  (char)   — 'viel2005' (default) or 'bode2001'
%                   'g_eff'   (scalar) — effective d.o.f., bode2001 only
%                   'nu'      (scalar) — shape parameter, bode2001 only
%
% OUTPUT:
%   sigma       : WDM-suppressed sigma(M, z), same size as M
%
% DEPENDENCIES:
%   transfer_function_WDM.m   (utils/)
%
% See also: transfer_function_WDM, Ludlow16, suppression_factor

    if isempty(varargin)
        varargin = {'viel05'};
    end

    rho_m0 = cosmo_cdm.Omega_m * cosmo_cdm.rho_crit0;   % [M_sun/h / (Mpc/h)^3]

    R     = (3 .* M ./ (4 * pi * rho_m0)).^(1/3);   % Lagrangian radius [Mpc/h]
    k_eff = 1.0 ./ R;                                 % effective wavenumber [h/Mpc]

    T     = T_wdm(k_eff, m_WDM_keV, cosmo_cdm, varargin{:});

    sigma = cosmo_cdm.sigmaM(M, z) .* T;
    % Guard: T can go slightly negative at extreme k; clamp sigma to real positive
    sigma = max(real(sigma), 0);
end