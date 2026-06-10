function Mstar = Child18_Mstar(z, delta_c, cosmo)
% Bisection on log10(M) to find sigma(Mstar, z) = delta_c
% sigma decreases with increasing M, so we bracket and bisect.

logM_lo = 6.0;
logM_hi = 17.0;

s_lo = cosmo.sigmaM(10^logM_lo, z);
s_hi = cosmo.sigmaM(10^logM_hi, z);

if s_lo < delta_c
    Mstar = 10^logM_lo;
    return
end
if s_hi > delta_c
    Mstar = 10^logM_hi;
    return
end

for iter = 1:200
    logM_mid = 0.5 * (logM_lo + logM_hi);
    s_mid    = cosmo.sigmaM(10^logM_mid, z);
    if abs(s_mid - delta_c) < 1e-6 || (logM_hi - logM_lo) < 1e-7
        break
    end
    if s_mid > delta_c
        logM_lo = logM_mid;
    else
        logM_hi = logM_mid;
    end
end

Mstar = 10^(0.5 * (logM_lo + logM_hi));
end