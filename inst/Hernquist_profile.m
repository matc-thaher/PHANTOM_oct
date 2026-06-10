function [rho, rhos, rs, fc] = Hernquist_profile(r, M, c, z, cosmo, Delta)
% Hernquist_profile  Hernquist (1990) density profile
%
%   [rho, rhos, rs, fc] = Hernquist_profile(r, M, c, z, cosmo, Delta)
%
%   INPUTS
%   r      : radii [Mpc/h], scalar or vector
%   M      : halo mass [Msun/h]
%   c      : concentration  R_Delta / r_s
%   z      : redshift
%   cosmo  : cosmology struct (needs .rhocrit0, .E)
%   Delta  : overdensity w.r.t. critical density (e.g. 200)
%
%   OUTPUTS
%   rho    : density [Msun/h / (Mpc/h)^3]
%   rhos   : characteristic density [same units]
%   rs     : scale radius [Mpc/h]
%   fc     : concentration-dependent factor
%
%   Profile form:
%     rho(r) = rhos / [ (r/rs) * (1 + r/rs)^3 ]
%
%   Reference: Hernquist 1990, ApJ 356, 359

rho_c   = cosmo.rho_crit0 * cosmo.E(z)^2;
R_Delta = (3*M / (4*pi * Delta * rho_c))^(1/3);
rs      = R_Delta / c;

% Enclosed mass: M(<r) = M_tot * (r/rs)^2 / (1 + r/rs)^2
% => M = M_tot * c^2/(1+c)^2  at R_Delta, so:
fc   = c.^2 ./ (1 + c).^2;
Mtot = M ./ fc;                    % total (infinite) mass
rhos = Mtot ./ (2*pi .* rs.^3);   % normalisation from Hernquist 1990

x   = r(:) ./ rs;
rho = rhos ./ (x .* (1 + x).^3);
end