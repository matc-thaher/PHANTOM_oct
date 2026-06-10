  function R = radius_from_mass( M,cosmo)
    % Mass in Msun/h → Radius in Mpc/h

        rho_m = cosmo.rho_m0;
        R = (3*M./(4*pi*rho_m)).^(1/3);
    end