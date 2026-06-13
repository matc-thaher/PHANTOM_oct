function c = Ludlow16_fit(M, z, cosmo)
% Ludlow16_concentration_fit  Appendix C fitting formula, Ludlow+2016
%
%   c = Ludlow16_fit(M, z, cosmo)
%
%   Broken power law in peak height nu = delta_c / sigma(M,z).
%   Calibrated for Planck cosmology (Table 2).
%   Valid for -8 <= log10(M/[h-1 Msun]) <= 16.5, 0 <= z <= 1 (1+z <= 2).
%
%   INPUTS
%   M     : halo mass [Msun/h], scalar or vector
%   z     : redshift (scalar)
%   cosmo : struct with cosmo.sigmaM(M, z) handle
%
%   OUTPUT
%   c     : concentration (same shape as M)
%
%   Reference: Ludlow et al. 2016, MNRAS 460, 1214, Appendix C, Eqs. C1-C6

a       = 1.0 / (1.0 + z);           % scale factor

% --- Eq. C2-C6: redshift-dependent parameters -------------------------
c0      = 3.395 * (1 + z)^(-0.215);                          % C2
nu0     = 0.307 * (1 + z)^( 0.540);                          % C3
gamma1  = 0.628 * (1 + z)^(-0.047);                          % C4
gamma2  = 0.317 * (1 + z)^(-0.893);                          % C5
mu      = (4.135 - 0.564*a^(-1) - 0.210*a^(-2) ...
               + 0.0557*a^(-3) - 0.00348*a^(-4)) ...
               / cosmo.D(z);                                  % C6

% --- Peak height -------------------------------------------------------
delta_c = 1.686;
sigma   = cosmo.sigmaM(M(:), z);
nu      = delta_c ./ sigma;

% --- Eq. C1: broken power law ------------------------------------------
x       = nu ./ nu0;
c       = c0 .* x.^(-gamma1) ./ (1 + x.^(1/mu)).^(mu*(gamma1 - gamma2));

c = reshape(c, size(M));
end