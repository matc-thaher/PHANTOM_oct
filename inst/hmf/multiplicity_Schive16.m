function f = multiplicity_Schive16(sigma, M, M_hm, z, cosmo, delta_c, base_model)
% Schive et al. (2016), ApJ 818, 89, & Schive 2025
% FDM (Fuzzy Dark Matter / Ultra-light axion) halo mass function.
%
% The FDM HMF is suppressed relative to CDM below the quantum Jeans mass.
% Schive+2016 parameterize the suppression as:
%
%   f_FDM = f_CDM * [1 + (M_h/(0.42M_hm))^-1.1]^-2.2    (Eq. 9)
%
% where M_hm is the FDM half-mode mass:
%   M_hm = 1.6e10 * (m_a / 1e-22 eV)^-1.5 * (Omega_m/0.3)^0.5 * h [M_sun]

% Validated by May & Springel (2023), MNRAS 524, 4256, who found broad
% agreement between this fitting function and the first self-consistent
% FDM halo mass function measured from Schrodinger-Poisson wave simulations.
% Du+2017 and Kulkarni+2022 both overpredict suppression relative to
% wave simulation measurements.
%
% INPUT:
%   sigma : sigma(M,z) — computed with CDM power spectrum
%   z     : redshift
%   M     : halo mass array [M_sun/h]
%   M_hm  : FDM half-mode mass [M_sun/h]
%   delta_c : (optional)
%   base_model : (optional) CDM base — default 'ST'
%                any model name accepted by multiplicity_dispatcher
%
% USAGE:
%   m_a  = 1e-22;   % axion mass [eV]
%   M_hm = 1.6e10 * (m_a/1e-22)^-1.5 * (cosmo.Omega_m/0.3)^0.5 * cosmo.h;
%   f    = multiplicity_FDM_Schive16(sigma, z, M, M_hm);

    
    if nargin < 6 || isempty(delta_c)
        delta_c = collapse_overdensity('corrections', true, 'z', z, 'cosmo', cosmo);
    end
    if nargin < 7 || isempty(base_model)
        base_model = 'ST';
    end

    f_CDM = multiplicity_dispatcher(base_model, sigma, delta_c);
    f     = f_CDM .* (1 + (M./(0.42.*M_hm)).^(-1.1)).^-2.2;
end