function Delta_vir = delta_vir_bn98(z, cosmo)
    Om_z = cosmo.Omega_m_z(z);
    x = Om_z - 1.0;
    Delta_vir = 18*pi^2 + 82*x - 39*x.^2;
end