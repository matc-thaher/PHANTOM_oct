function f = multiplicity_Lovell14(sigma, M, M_hm, z, cosmo, delta_c)
% Lovell et al. (2014), MNRAS 439, 300
% WDM suppression applied as a ratio to the CDM Sheth-Tormen f(sigma).
%
% f_WDM(sigma) = f_CDM(sigma) * (1 + gamma * M_hm/M)^beta
%
% where alpha = 0.6 (Lovell+2014, Eq. 5) and M_hm is the half-mode mass.
%
% INPUT:
%   sigma : sigma(M,z)
%   z     : redshift
%   M     : halo mass array [M_sun/h], same size as sigma
%   M_hm  : half-mode mass [M_sun/h]; compute from cosmo.transfer_model
%   delta_c : (optional)
%
% NOTE: M_hm = (4*pi/3)*rho_m * R_hm^3 where R_hm is where T_WDM^2 = 0.5

    if nargin < 6 || isempty(delta_c)
        delta_c = collapse_overdensity('corrections', true, 'z', z, 'cosmo', cosmo);
    end

    beta   = 0.99; gamma = 2.7;
    f_CDM   = multiplicity_ST(sigma, delta_c);
    f       = f_CDM .* (1 + (gamma .* M_hm ./ M)).^beta;
end