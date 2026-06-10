% =====================================================================
% line-of-sight comoving distance in Mpc/h
% =====================================================================
function d = comoving_distance(z, cosmo)

    c_km_s = 299792.458;
    H0 = cosmo.H0;

    z = z(:);
    d = zeros(size(z));

    for i = 1:length(z)
        d(i) = (c_km_s / H0) * quadgk(@(zp) 1 ./ cosmo.E(zp), 0, z(i));
    end

    % Convert Mpc to Mpc/h
    d = d * cosmo.h;
end