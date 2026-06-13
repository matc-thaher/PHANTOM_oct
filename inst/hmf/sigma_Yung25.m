function [sigma, dlns_dlnM] = sigma_Yung25(M_phys)
% sigma(M) fitting function from Yung+2025, Eq. 3.
% INPUT:  M_phys in physical M_sun (NO little-h).
% OUTPUT: sigma(M), dimensionless.
%
% WARNING: Do not pass M in h^-1 M_sun — see Yung+2025 Section 2.2.1.

    M12     = M_phys ./ 1e12;       % M_{h,12} = M / 10^12 M_sun
    y       = 1 ./ M12;
    sigma   = (26.8 .* y.^0.41) ./ ...
              (1 + 6.18.*y.^0.23 + 4.64.*y.^0.37);
    dlns_dlnM = gradient(log(sigma), log(M_phys));
end