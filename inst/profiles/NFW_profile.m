function rho = NFW_profile(r, rhos, rs)
% NFW_profile  Navarro-Frenk-White (1997) density profile
%
%   [rho, rhos, rs] = NFW_profile(r, M, c, z, cosmo, Delta)
%
%   INPUTS
%   r      : radii [Mpc/h], scalar or vector
%   rhos   : density from simulation [Msun/kpc^3]
%   rs     : scale radius
%
%   OUTPUTS
%   rho    : density [Msun / kpc^3]
%
%   Reference: Navarro, Frenk & White 1997, ApJ 490, 493

  x = r ./ rs;
  rho = rhos ./ (x .* (1 + x).^2);
end