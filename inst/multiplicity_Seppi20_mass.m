function f = multiplicity_Seppi20_mass(sigma, z, cosmo)
% 1D Seppi+2020 multiplicity marginalized over xoff and spin
    f = multiplicity_Seppi20(sigma, z, cosmo, [], [], false, true, true);
end