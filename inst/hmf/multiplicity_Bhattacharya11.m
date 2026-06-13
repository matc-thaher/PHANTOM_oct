function f = multiplicity_Bhattacharya11(sigma, z, cosmo, delta_c)
% Bhattacharya et al. (2011), ApJ 732, 122, Table 4
% FOF, calibrated z = 0-2. Explicit z-dependence; delta_c without correction.
% mass range 6e11-3e15 solar mass
    if nargin < 3 || isempty(delta_c)
        delta_c = collapse_overdensity('corrections', true, 'z', z, 'cosmo', cosmo);   % paper uses uncorrected value
    end
    nu   = delta_c ./ sigma;
    nu2  = nu.^2;
    zp1  = 1 + z;
    A = 0.333 * zp1^-0.11;
    a = 0.788 * zp1^-0.01;
    p = 0.807;
    q = 1.795;
    f = A .* sqrt(2/pi) .* exp(-a.*nu2./2) .* (1 + (a.*nu2).^-p) .* (nu.*sqrt(a)).^q;
end