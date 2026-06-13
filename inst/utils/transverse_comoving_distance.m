% =====================================================================
% transverse comoving distance in Mpc/h
% =====================================================================
function d = transverse_comoving_distance(z, cosmo)

    dc = comoving_distance(z, cosmo);

    if abs(cosmo.Omega_k) < 1e-12
        d = dc;
        return
    end

    c_km_s = 299792.458;
    H0 = cosmo.H0;
    Dh = (c_km_s / H0) * cosmo.h;   % Mpc/h

    if cosmo.Omega_k > 0
        sqrtOk = sqrt(cosmo.Omega_k);
        d = Dh / sqrtOk .* sinh(sqrtOk .* dc ./ Dh);
    else
        sqrtOk = sqrt(-cosmo.Omega_k);
        d = Dh / sqrtOk .* sin(sqrtOk .* dc ./ Dh);
    end
end