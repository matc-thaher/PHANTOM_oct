function f = multiplicity_Crocce10(sigma, z)
% Crocce et al. (2010), MNRAS 403, 1353, Eq. 22
% FOF calibration, z = 0-1. Explicit redshift evolution in all parameters.
    zp = 1 + z;
    A = 0.58  * zp^-0.13;
    a = 1.37  * zp^-0.15;
    b = 0.30  * zp^-0.084;
    c = 1.036 * zp^-0.024;
    sigma = sigma(:);
    f = A .* (sigma.^-a + b) .* exp(-c ./ sigma.^2);
end