function [c, z_form] = Ludlow16(M, z_obs, cosmo, f)
% Ludlow16_concentration  Ludlow et al. (2016) analytic concentration model
%
%   [c, z_form] = Ludlow16(M, z_obs, cosmo)
%   [c, z_form] = Ludlow16(M, z_obs, cosmo, f)
%
%   Implements the analytic c(M,z) model from Section 4.3 of Ludlow+2016.
%   The model works as follows:
%
%   1. Compute the Collapsed Mass History (CMH) from EPS theory (eq. 3):
%        Mcoll(z) = M0 * erfc( [delta_sc(z)-delta_sc(z0)] /
%                               sqrt(2*(sigma^2(fM0)-sigma^2(M0))) )
%
%   2. Find the formation redshift z_f at which Mcoll(z_f) = f * M0.
%
%   3. Set the characteristic density of the halo equal to the critical
%      density of the Universe at z_f:
%        <rho>(r_{-2}) = C * rho_crit(z_f)
%      where C is a calibration constant (C=550 from Ludlow+2016 eq.12).
%
%   4. Invert <rho>(r_{-2}) = 200*rho_crit(z_obs)*c^3 / (3*f(c))
%      numerically to find c.
%
%   INPUTS
%   M      : halo mass [Msun/h], scalar or vector
%   z_obs  : observation redshift (scalar)
%   cosmo  : struct with fields:
%              cosmo.sigmaM(M, z)    — rms linear fluctuation
%              cosmo.growthFactor(z) — linear growth factor D(z), D(0)=1
%              cosmo.rhocrit(z)      — critical density [Msun/h/(Mpc/h)^3]
%   f      : CMH progenitor fraction (default 0.02, paper eq.3)
%
%   OUTPUTS
%   c      : concentration, same shape as M
%   z_form : formation redshift, same shape as M
%
%   Reference: Ludlow, Bose, Angulo et al. 2016, MNRAS 460, 1214

if nargin < 4 || isempty(f)
    f = 0.02;
end

M      = M(:);
c      = zeros(size(M));
z_form = zeros(size(M));

% calibration constant C from eq.(12) of Ludlow+2016
C_cal = 650.0;

% for i = 1:numel(M)
%     M0 = M(i);
% 
%     % Step 1+2: find formation redshift from CMH
%     zf = Ludlow16_formation_z(M0, f, cosmo);
%     z_form(i) = zf;
% 
%     % Step 3: target characteristic density
%     % rho_target = C_cal * cosmo.rhocrit(zf);
%     rho_target = C_cal * cosmo.rho_crit0 * cosmo.E(zf)^2;
% 
%     % % Step 4: solve for c from NFW mean enclosed density relation
%     % %   rho_char(c) = (200/3) * rhocrit(z_obs) * c^3 / f(c) = rho_target
%     % % rho_obs = cosmo.rhocrit(z_obs);
%     % rho_obs = cosmo.rho_crit0 * cosmo.E(z_obs)^2;
%     % residual_c = @(cv) (200/3) * rho_obs * cv^3 / (log(1+cv) - cv/(1+cv)) ...
%     %                    - rho_target;
%     % 
%     % % bracket: c is between 1 and 1000
%     % c_lo = 1.0;
%     % c_hi = 1000.0;
%     % 
%     % if sign(residual_c(c_lo)) == sign(residual_c(c_hi))
%     %     c(i) = NaN;
%     % else
%     %     c(i) = fzero(residual_c, [c_lo, c_hi], ...
%     %                  optimset('TolX', 1e-6, 'Display', 'off'));
%     % end
%     % Step 4: solve for c using Einasto profile (alpha=0.18, fixed per Ludlow+16)
%     % M(<r_{-2}) / M(<R_200c) = einasto_mass_ratio(c, 0.18)
%     alpha_ein = 0.18;
%     c_arr     = logspace(0, 2, 500);
%     M_ratio   = einasto_mass_ratio(c_arr, alpha_ein);
% 
%     % Ratio of critical densities (eq. 12 rearranged)
%     rho_ratio = cosmo.E(zf)^2 / cosmo.E(z_obs)^2;
% 
%     % Solve: C_cal * rho_ratio = 200 * c^3 * M_ratio(c) / 3
%     % i.e.:  c^3 * M_ratio(c) = 3 * C_cal * rho_ratio / 200
%     lhs      =  c_arr.^3 .* M_ratio;
%     rhs_val  = C_cal * rho_ratio / 200;
% 
%     residual = lhs - rhs_val;
% 
%     if all(residual > 0) || all(residual < 0)
%         c(i) = NaN;
%     else
%         c(i) = interp1(residual, c_arr, 0.0, 'pchip');
%     end
%     end
alpha_ein = 0.18;
for i = 1:numel(M)
    M0 = M(i);
    [zf, c_val] = Ludlow16_formation_z(M0, f, cosmo, z_obs, [], C_cal, alpha_ein);
    z_form(i) = zf;
    c(i)      = c_val;
end

c      = reshape(c,      size(M));
z_form = reshape(z_form, size(M));
end