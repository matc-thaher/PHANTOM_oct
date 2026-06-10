function f = multiplicity_Comparat17(sigma, z, cosmo, delta_c)
% multiplicity_Comparat17  Halo multiplicity f(sigma) for Comparat et al. (2017)
%
%   f = multiplicity_Comparat17(sigma, z, cosmo, delta_c)
%
%   This implementation follows the Comparat+2017 calibration of the
%   Bhattacharya+2011 functional form, valid at z=0 for virial SO masses.
%   Outside that regime it assumes universality.
%
% INPUTS:
%   sigma : rms linear fluctuation (array)
%   z     : redshift (scalar) – model calibrated at z~0
%   mdef  : mass definition string; should be 'vir' or equivalent
%   cosmo : PHANTOM cosmo struct with cosmo.collapseOverdensity(z)
%
% OUTPUT:
%   f     : multiplicity function f(sigma), same size as sigma

    if nargin < 4 || isempty(delta_c)
        delta_c = collapse_overdensity('corrections', true, 'z', z, 'cosmo', cosmo);
    end

    nu  = delta_c ./ sigma;
    nu2 = nu.^2;

    % Comparat+2017 parameters (updated wrt original Bhattacharya+2011)
    A = 0.324;
    a = 0.897;
    p = 0.624;
    q = 1.589;

    % Bhattacharya-type functional form
    % f(nu) = A * sqrt(2/pi) * exp(-a nu^2 / 2) * [1 + (a nu^2)^(-p)] * (nu sqrt(a))^q
    pref = A * sqrt(2.0/pi);
    core = exp(-0.5 * a .* nu2);
    corr = 1.0 + (a .* nu2).^(-p);
    pow  = (nu .* sqrt(a)).^q;

    f_vec = pref .* core .* corr .* pow;

    f = reshape(f_vec, size(sigma));
end