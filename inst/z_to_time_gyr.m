% =====================================================================
% lookback time / age
% =====================================================================
function out = z_to_time_gyr(z, cosmo)

    Mpc_km = 3.0856775814913673e19;
    H0_SI = cosmo.H0 / Mpc_km;
    sec_per_Gyr = 365.25 * 24 * 3600 * 1e9;

    integrand = @(zp) 1 ./ ((1 + zp) .* cosmo.E(zp));

    z = z(:);
    tL = zeros(size(z));

    for i = 1:length(z)
        tL(i) = (1 / H0_SI) * quadgk(integrand, 0, z(i));
    end

    t0 = (1 / H0_SI) * quadgk(integrand, 0, cosmo.zmax);

    out.lookback_Gyr = tL / sec_per_Gyr;
    out.t0_Gyr       = t0 / sec_per_Gyr;
    out.age_Gyr      = (t0 - tL) / sec_per_Gyr;
end