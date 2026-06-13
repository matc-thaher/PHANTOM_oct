function rho = Hernquist_fit_wrapper(r, rhos, rs)
% HERNQUIST_FIT_WRAPPER  Hernquist profile for direct fitting.
%
% rho(r) = rhos / [ (r/rs) * (1 + r/rs)^3 ]
%
% Free parameters: rhos [Msun/kpc^3], rs [kpc]
%
% Reference: Hernquist 1990, ApJ 356, 359

  r   = r(:);
  x   = r ./ rs;
  rho = rhos ./ (x .* (1 + x).^3);
end