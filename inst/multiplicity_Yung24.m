function f = multiplicity_Yung24(sigma, z)
% Yung et al. (2024), MNRAS 530, 4868, Table A1 & Eq A2
% Virial SO, Planck cosmology. Calibrated at high-z using GUREFT simulations.
% Valid range: z ~ 0-12.
    A = 0.11416632 - 0.01486746.*z + 0.00137191.*z.^2;
    a = 1.05274399 + 0.02803087.*z - 0.00306126.*z.^2;
    b = 8.62813020 + 0.00384969.*z - 0.02349983.*z.^2;
    c = 1.13138924 + 0.01713172.*z - 0.00113630.*z.^2;
    f = A .* ((sigma./b).^-a + 1) .* exp(-c./sigma.^2);
end