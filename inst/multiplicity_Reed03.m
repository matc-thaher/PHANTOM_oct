function f = multiplicity_Reed03(sigma, z, cosmo, delta_c)
% Reed et al. (2003), MNRAS 346, 565, Eq. 9
% ST correction at high mass. FOF.
    if nargin < 4 || isempty(delta_c), delta_c = collapse_overdensity('corrections', true, 'z', z, 'cosmo', cosmo); end
    f_ST = multiplicity_ST(sigma, delta_c);
    f    = f_ST .* exp(-0.7 ./ (sigma .* cosh(2.*sigma).^5));
end