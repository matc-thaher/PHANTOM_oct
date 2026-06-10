function [T_FDM, P_FDM] = T_Schive25(k, m22, P_CDM)
% fdm_transfer   FDM linear matter power spectrum and transfer function.
%
% Implements Schive (2025, Living Reviews in Comp. Astrophysics),
% Eqs. (19) and (20). Original fitting formula from Hu et al. (2000).
%
% Eq. (19):  P_FDM(k) = T_FDM^2(k) * P_CDM(k)
%
% Eq. (20):  T_FDM(k) = cos(x^3) / (1 + x^8)
%            x = 1.61 * m22^(1/18) * k / k_Jeq
%            k_Jeq = 9 * m22^(1/2)   [Mpc^-1]
%
% The transfer function is redshift-independent: suppression is fixed
% during the radiation-dominated era (Hu et al. 2000).
%
% INPUT:
%   k      : wavenumber array [Mpc^-1]
%   m22    : FDM boson mass in units of 1e-22 eV/c^2  (dimensionless)
%   P_CDM  : CDM linear power spectrum at same k values (optional)
%            If not provided, only T_FDM is returned.
%
% OUTPUT:
%   P_FDM  : FDM linear power spectrum (same units as P_CDM)
%            Returns empty [] if P_CDM not supplied.
%   T_FDM  : FDM transfer function T(k), dimensionless, range [0,1]
%
% EXAMPLE:
%   k      = logspace(-2, 2, 500);       % [Mpc^-1]
%   T      = fdm_transfer(k, 1.0);       % m22 = 1, CDM not needed
%   [P, T] = fdm_transfer(k, 1.0, P_CDM);

    % Eq. 20: Jeans wavenumber at matter-radiation equality [Mpc^-1]
    k_Jeq = 9.0 * m22^(0.5);

    % Eq. 20: dimensionless wavenumber x
    x = 1.61 * m22^(1/18) .* (k ./ k_Jeq);

    % Eq. 20: FDM transfer function
    T_FDM = cos(x.^3) ./ (1 + x.^8);

    % Eq. 19: FDM power spectrum (only if P_CDM supplied)
    if nargin < 3 || isempty(P_CDM)
        P_FDM = [];
    else
        P_FDM = T_FDM.^2 .* P_CDM;
    end
end