function [zf, c_out] = Ludlow16_formation_z(M0, f, cosmo, z_obs, c_arr, C_cal, alpha_ein)
% Ludlow16_formation_z  Formation redshift and concentration from coupled
%                        density + EPS condition (Ludlow+2016 eqs. 6-7)
%
%   [zf, c_out] = Ludlow16_formation_z(M0, f, cosmo, z_obs, c_arr, C_cal, alpha_ein)
%
%   Sweeps over a trial concentration array c_arr. For each c:
%     1. Finds zf from the density condition (eq.6):
%           E(zf)^2 = 200 * c^3 * M_ratio(c) / C_cal * E(z_obs)^2
%     2. Evaluates CMH(zf) = erfc(...) via Ludlow16_CMH
%   Then finds the c (and zf) where CMH(zf) = f (eq.7).
%
%   INPUTS
%   M0        : halo mass [Msun/h]
%   f         : progenitor fraction (default 0.02)
%   cosmo     : struct with cosmo.D(z), cosmo.sigmaM(M,z), cosmo.E(z)
%   z_obs     : observation redshift (default 0)
%   c_arr     : concentration trial array (default logspace(0,2,200))
%   C_cal     : calibration constant (default 650)
%   alpha_ein : Einasto shape parameter (default 0.18)
%
%   OUTPUTS
%   zf        : formation redshift
%   c_out     : concentration

if nargin < 4 || isempty(z_obs),     z_obs     = 0;                    end
if nargin < 5 || isempty(c_arr),     c_arr     = logspace(0, 2, 200);  end
if nargin < 6 || isempty(C_cal),     C_cal     = 650.0;                end
if nargin < 7 || isempty(alpha_ein), alpha_ein = 0.18;                 end

% M_ratio(c) = M(<r_{-2}) / M(<R_200c) — decreasing with c
M_ratio = einasto_mass_ratio(c_arr, alpha_ein);

% Step 1: for each c, find zf from density condition (eq.6)
rho_f_rho_c = 200.0 .* c_arr.^3 .* M_ratio ./ C_cal;
E2_obs      = cosmo.E(z_obs)^2;
zf_arr      = zeros(size(c_arr));

for j = 1:numel(c_arr)
    target = rho_f_rho_c(j) * E2_obs;
    obj    = @(zz) cosmo.E(zz)^2 - target;
    if obj(0) >= 0
        zf_arr(j) = 0;
    elseif obj(50) <= 0
        zf_arr(j) = 50;
    else
        zf_arr(j) = fzero(obj, [0,50], optimset('TolX',1e-5,'Display','off'));
    end
end

% Step 2: evaluate CMH at each zf — this is erfc(...) = Mcoll/M0
cmh_arr = Ludlow16_CMH(zf_arr, M0, f, cosmo);

% Step 3: find c where CMH(zf(c)) = f  (eq.7: M_ratio = erfc)
residual = M_ratio - cmh_arr;

% Find the LAST zero crossing (highest c root — physically correct)
idx = find(diff(sign(residual)) ~= 0, 1, 'last');

if isempty(idx)
    zf    = NaN;
    c_out = NaN;
else
    % Interpolate in the interval [idx, idx+1] only
    c_out = interp1(residual(idx:idx+1), c_arr(idx:idx+1),  0.0, 'linear');
    zf    = interp1(residual(idx:idx+1), zf_arr(idx:idx+1), 0.0, 'linear');
end
end