function c = Bhattacharya13(M, z, mdef, cosmo)
% Bhattacharya13_concentration  Bhattacharya et al. (2013) model
%
%   c = Bhattacharya13_concentration(M, z, mdef, cosmo)
%
%   Power-law in c-nu space:
%       c = A * (1+z)^B * nu^D
%   where nu = delta_c / sigma(M,z).
%   Calibrated for WMAP7 cosmology.
%   Valid range: 2e12 < M < 2e15 Msun/h,  0 < z < 2.
%
%   INPUTS
%   M     : halo mass [Msun/h]
%   z     : redshift
%   mdef  : '200c' | 'vir' | '200m'
%   cosmo : cosmology struct with cosmo.sigmaM
%
%   Reference: Bhattacharya et al. 2013, ApJ 766, 32, Table 2

delta_c = 1.686;
sigma   = cosmo.sigmaM(M, z);
nu      = delta_c ./ sigma;

switch lower(mdef)
    case '200c'
        % All haloes, NFW fit
        A =  5.9;  B = -0.35;  D = -0.44;
    case 'vir'
        A =  6.6;  B = -0.26;  D = -0.45;
    case '200m'
        A =  9.0;  B = -0.47;  D = -0.54;
    otherwise
        error('Bhattacharya13: unknown mdef "%s". Use 200c, vir, or 200m.', mdef);
end

c = A .* (1 + z).^B .* nu.^D;
end