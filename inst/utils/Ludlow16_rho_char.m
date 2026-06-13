function rho_char = Ludlow16_rho_char(c, z_obs, cosmo)
% Ludlow16_rho_char  Characteristic density of an NFW halo
%
%   rho_char = Ludlow16_rho_char(c, z_obs, cosmo)
%
%   The mean enclosed density within r_{-2} = r200/c for an NFW halo
%   of concentration c observed at redshift z_obs:
%
%     <rho>(r_{-2}) = 200 * rho_crit(z_obs) * c^3 / (3 * f(c))
%
%   where  f(c) = ln(1+c) - c/(1+c)
%
%   This is the quantity that the Ludlow+2016 model equates to the
%   critical density of the Universe at the halo's formation redshift.
%
%   INPUTS
%   c      : concentration (scalar or array)
%   z_obs  : observation redshift
%   cosmo  : struct with cosmo.rhocrit(z) [Msun/h / (Mpc/h)^3]
%
%   OUTPUT
%   rho_char : mean enclosed density [same units as rhocrit]
%
%   Reference: Ludlow et al. 2016, MNRAS 460, 1214, Section 4.3

fc       = log(1 + c) - c ./ (1 + c);
rho_crit = cosmo.rhocrit(z_obs);
rho_char = (200/3) .* rho_crit .* c.^3 ./ fc;
end