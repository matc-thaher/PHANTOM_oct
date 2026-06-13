function f = multiplicity_RodriguezPuebla16(sigma, z)
% Rodriguez-Puebla et al. (2016), MNRAS 462, 893, Eq 32 & 25
% Virial SO definition, Planck cosmology, z = 0-7.
    A = 0.144 - 0.011.*z + 0.003.*z.^2;
    a = 1.351 + 0.068.*z + 0.006.*z.^2;
    b = 3.113 - 0.077.*z - 0.013.*z.^2;
    c = 1.187 + 0.009.*z;
    f = A .* ((sigma./b).^-a + 1) .* exp(-c./sigma.^2);
end