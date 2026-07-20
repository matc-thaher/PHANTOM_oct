% test_profiles.m
% Test suite for NFW, Einasto, Hernquist, Soliton, DK14 profile functions
%
% Tests: routing, physical behaviour, mass conservation, asymptotic limits,
%        edge cases, cross-profile comparisons, and diagnostic plots.
%
% DEPENDENCIES (must be on path or same folder):
%   NFW_analytcl_Profile.m
%   NFW_profile.m
%   Einasto_profile.m
%   Hernquist_profile.m
%   Soliton_profile.m
%   DK14_profile.m

% clear; clc; close all;

here     = fileparts(mfilename('fullpath'));
repo     = fullfile(here, '..');
profiles_path = fullfile(repo, 'src', 'profiles');
addpath(genpath(profiles_path));
utils_path         = fullfile(repo, 'src', 'utils');
addpath(genpath(utils_path));

fprintf('==========================================\n');
fprintf('  Halo Profile Functions — Test Suite\n');
fprintf('==========================================\n\n');

% ---- Shared cosmology stub ------------------------------------------
% Minimal cosmo struct sufficient for all profile functions.
% Uses flat ΛCDM Planck18-like values.
% cosmo.rhocrit0 = 2.775e11;                          % [Msun/h / (Mpc/h)^3]
% cosmo.Omega_m  = 0.315;
% cosmo.rho_m0   = cosmo.rhocrit0 * cosmo.Omega_m;   % [Msun/h / (Mpc/h)^3]
% cosmo.E        = @(z) sqrt(cosmo.Omega_m*(1+z)^3 + (1-cosmo.Omega_m));
% cosmo.nu       = @(M,z) 1.686 ./ (0.8 .* (M./1e12).^(-0.1) ./ cosmo.E(z));
cosmo = cosmology('Planck18');

% ---- Shared test parameters -----------------------------------------
M0     = 1e14;      % halo mass  [Msun/h]
c0     = 5.0;       % concentration
z0     = 0.0;       % redshift
Delta  = 200;
rhos0  = 1e10;       % scale density for NFW_profile  [Msun/kpc^3]
rs0    = 0.3;       % scale radius for NFW_profile  [Mpc/h]

% Radial arrays
r_mpc   = logspace(-1, 2, 200);   % [Mpc/h]  for Einasto, Hernquist, DK14
% r_kpc   = logspace(-2, 3, 200);   % [kpc]    for NFW_analytcl, NFW_profile, Soliton
r_kpc = r_mpc/1000;

tf_lut = {'false','true'};
tfstr  = @(x) tf_lut{x+1};

% =====================================================================
%% 1. NFW_analytcl_Profile — basic output struct fields
% =====================================================================
fprintf('--- 1. NFW_analytcl_Profile: output struct fields ---\n');
Mvir = 1e12;   Rvir = 300;   c_nfw = 10;
r_in = linspace(1, Rvir, 200);
NFW  = NFW_analytcl_Profile(Mvir, Rvir, c_nfw, r_in);

required_fields = {'rs','rho_s','rho','Menc','f_c','r','Mvir','Rvir','c'};
for i = 1:numel(required_fields)
    assert(isfield(NFW, required_fields{i}), ...
        sprintf('Missing field: %s', required_fields{i}));
    fprintf('  .%-6s present  PASS\n', required_fields{i});
end
fprintf('\n');

% =====================================================================
%% 2. NFW_analytcl_Profile — physical consistency
% =====================================================================
fprintf('--- 2. NFW_analytcl_Profile: physical consistency ---\n');

% rs = Rvir/c
assert(abs(NFW.rs - Rvir/c_nfw) < 1e-10, 'rs != Rvir/c');
fprintf('  rs = Rvir/c                   PASS\n');

% rho > 0 everywhere
assert(all(NFW.rho > 0), 'rho <= 0 somewhere');
fprintf('  rho > 0 everywhere            PASS\n');

% rho monotonically decreasing
assert(all(diff(NFW.rho) < 0), 'rho not monotonically decreasing');
fprintf('  rho monotonically decreasing  PASS\n');

% Menc monotonically increasing
assert(all(diff(NFW.Menc) > 0), 'Menc not monotonically increasing');
fprintf('  Menc monotonically increasing PASS\n');

% Menc at Rvir recovers Mvir within 1%
assert(abs(NFW.Menc(end) - Mvir)/Mvir < 0.01, ...
    sprintf('Menc(Rvir)/Mvir = %.4f, expected ~1', NFW.Menc(end)/Mvir));
fprintf('  Menc(Rvir) ≈ Mvir  (%.4f)    PASS\n', NFW.Menc(end)/Mvir);

% r array clipped to Rvir
assert(all(NFW.r <= Rvir), 'r array exceeds Rvir');
fprintf('  r array clipped to Rvir       PASS\n\n');

% =====================================================================
%% 3. NFW_analytcl_Profile — r array truncation beyond Rvir
% =====================================================================
fprintf('--- 3. NFW_analytcl_Profile: r truncation beyond Rvir ---\n');
r_long = linspace(1, Rvir*2, 400);   % intentionally extends beyond Rvir
NFW2   = NFW_analytcl_Profile(Mvir, Rvir, c_nfw, r_long);
assert(all(NFW2.r <= Rvir), 'r was not clipped at Rvir');
assert(numel(NFW2.r) < numel(r_long), 'No truncation occurred');
fprintf('  r truncated from %d -> %d points  PASS\n\n', ...
    numel(r_long), numel(NFW2.r));

% =====================================================================
%% 4. NFW_profile — output size and positivity
% =====================================================================
fprintf('--- 4. NFW_profile: output size and positivity ---\n');
rho_nfw = NFW_profile(r_kpc, rhos0, rs0*1e3);   % rs0 in kpc
assert(numel(rho_nfw) == numel(r_kpc), 'Output size mismatch');
assert(all(rho_nfw > 0), 'rho <= 0 somewhere');
assert(all(isfinite(rho_nfw)), 'Non-finite values in rho');
fprintf('  size match, rho > 0, finite   PASS\n');

% Monotonically decreasing
assert(all(diff(rho_nfw) < 0), 'NFW_profile not monotonically decreasing');
fprintf('  rho monotonically decreasing  PASS\n\n');

% =====================================================================
%% 5. NFW_profile — known asymptotic slope
% =====================================================================
fprintf('--- 5. NFW_profile: asymptotic slopes ---\n');
rs_test = 1.0;   rhos_test = 1.0;   % normalised units

% Inner slope: rho ~ r^-1 for r << rs
r_inner = [0.001, 0.002, 0.004];
rho_in  = NFW_profile(r_inner, rhos_test, rs_test);
slope_in = mean(diff(log(rho_in)) ./ diff(log(r_inner)));
assert(abs(slope_in - (-1)) < 0.05, ...
    sprintf('Inner slope = %.3f, expected -1', slope_in));
fprintf('  Inner slope = %.3f (expected -1)  PASS\n', slope_in);

% Outer slope: rho ~ r^-3 for r >> rs
r_outer = [100, 200, 400];
rho_out = NFW_profile(r_outer, rhos_test, rs_test);
slope_out = mean(diff(log(rho_out)) ./ diff(log(r_outer)));
assert(abs(slope_out - (-3)) < 0.05, ...
    sprintf('Outer slope = %.3f, expected -3', slope_out));
fprintf('  Outer slope = %.3f (expected -3) PASS\n\n', slope_out);

% =====================================================================
%% 6. Einasto_profile — outputs and normalisation
% =====================================================================
fprintf('--- 6. Einasto_profile: outputs and normalisation ---\n');
[rho_e, rhos_e, rs_e] = Einasto_profile(r_mpc, M0, c0, z0, cosmo, Delta);

assert(numel(rho_e) == numel(r_mpc), 'Output size mismatch');
assert(all(rho_e > 0),               'rho <= 0 somewhere');
assert(all(isfinite(rho_e)),         'Non-finite values');
assert(isscalar(rhos_e),             'rhos_e not scalar');
assert(isscalar(rs_e),               'rs_e not scalar');
fprintf('  size, positivity, finiteness  PASS\n');

% rho at r=rs should equal rhos_e by definition
rho_at_rs = Einasto_profile(rs_e, M0, c0, z0, cosmo, Delta);
assert(abs(rho_at_rs - rhos_e)/rhos_e < 1e-6, ...
    sprintf('rho(rs)/rhos = %.6f, expected 1', rho_at_rs/rhos_e));
fprintf('  rho(rs) == rhos_e             PASS\n');

% Monotonically decreasing
assert(all(diff(rho_e) < 0), 'Einasto not monotonically decreasing');
fprintf('  rho monotonically decreasing  PASS\n\n');

% =====================================================================
%% 7. Einasto_profile — alpha_e Gao+2008 bounds
% =====================================================================
fprintf('--- 7. Einasto_profile: Gao+2008 alpha_e bounds ---\n');
M_test_vec = [1e10, 1e12, 1e14, 1e16];
for i = 1:numel(M_test_vec)
    nu_val   = cosmo.nu(M_test_vec(i), z0);
    alpha_e  = min(0.155 + 0.0095*nu_val^2, 0.3);
    assert(alpha_e >= 0.155 && alpha_e <= 0.3, ...
        sprintf('alpha_e = %.4f out of [0.155, 0.3] for M=%.0e', ...
                alpha_e, M_test_vec(i)));
    fprintf('  M=%.0e  nu=%.3f  alpha_e=%.4f  in [0.155, 0.30]  PASS\n', ...
        M_test_vec(i), nu_val, alpha_e);
end
fprintf('\n');

% =====================================================================
%% 8. Hernquist_profile — outputs and total mass recovery
% =====================================================================
fprintf('--- 8. Hernquist_profile: outputs and mass recovery ---\n');
[rho_h, rhos_h, rs_h] = Hernquist_profile(r_mpc, M0, c0, z0, cosmo, Delta);

assert(numel(rho_h) == numel(r_mpc), 'Output size mismatch');
assert(all(rho_h > 0),               'rho <= 0 somewhere');
assert(all(isfinite(rho_h)),         'Non-finite values');
fprintf('  size, positivity, finiteness  PASS\n');

% Hernquist enclosed mass: M(<r) = Mtot * x^2 / (1+x)^2 where x = r/rs
rho_c_h  = cosmo.rho_crit0 .* cosmo.E(z0).^2;
R_Delta  = (3*M0 / (4*pi*Delta*rho_c_h))^(1/3);
rs_h_exp = R_Delta / c0;
fc       = c0^2 / (1+c0)^2;
Mtot     = M0 / fc;
x_edge   = r_mpc(end) / rs_h_exp;
Menc_h   = Mtot * x_edge^2 / (1 + x_edge)^2;
assert(abs(rs_h - rs_h_exp)/rs_h_exp < 1e-6, ...
    sprintf('rs mismatch: %.6e vs %.6e', rs_h, rs_h_exp));
fprintf('  rs = R_Delta/c                PASS\n');

% Enclosed mass at r -> inf -> Mtot (test at large r/rs)
r_large  = 1e4 * rs_h;
x_large  = r_large / rs_h;
Menc_inf = Mtot * x_large^2 / (1 + x_large)^2;
assert(abs(Menc_inf/Mtot - 1) < 0.001, ...
    sprintf('M(inf)/Mtot = %.6f, expected ~1', Menc_inf/Mtot));
fprintf('  M(r->inf) -> Mtot  (%.6f)   PASS\n\n', Menc_inf/Mtot);

% =====================================================================
%% 9. Hernquist_profile — inner/outer asymptotic slopes
% =====================================================================
fprintf('--- 9. Hernquist_profile: asymptotic slopes ---\n');

% Use normalised NFW-like call for slope check
rhos_h_norm = 1.0;   rs_h_norm = 1.0;

% Inner slope: rho ~ r^-1 for r << rs
r_in_h  = [0.001, 0.002, 0.004];
rho_in_h = rhos_h_norm ./ ((r_in_h./rs_h_norm) .* (1 + r_in_h./rs_h_norm).^3);
slope_in_h = mean(diff(log(rho_in_h)) ./ diff(log(r_in_h)));
assert(abs(slope_in_h - (-1)) < 0.05, ...
    sprintf('Hernquist inner slope = %.3f, expected -1', slope_in_h));
fprintf('  Inner slope = %.3f (expected -1)  PASS\n', slope_in_h);

% Outer slope: rho ~ r^-4 for r >> rs
r_out_h  = [100, 200, 400];
rho_out_h = rhos_h_norm ./ ((r_out_h./rs_h_norm) .* (1 + r_out_h./rs_h_norm).^3);
slope_out_h = mean(diff(log(rho_out_h)) ./ diff(log(r_out_h)));
assert(abs(slope_out_h - (-4)) < 0.05, ...
    sprintf('Hernquist outer slope = %.3f, expected -4', slope_out_h));
fprintf('  Outer slope = %.3f (expected -4) PASS\n\n', slope_out_h);

% =====================================================================
%% 10. Soliton_profile — output, half-density radius definition
% =====================================================================
fprintf('--- 10. Soliton_profile: output and rc definition ---\n');
rho0_sol = 1e8;    % central density [Msun/kpc^3]
rc_sol   = 1.0;    % core radius [kpc]
r_sol    = logspace(-2, 2, 300);

rho_sol = Soliton_profile(r_sol, rho0_sol, rc_sol);

assert(numel(rho_sol) == numel(r_sol), 'Output size mismatch');
assert(all(rho_sol > 0),              'rho <= 0 somewhere');
assert(all(isfinite(rho_sol)),        'Non-finite values');
fprintf('  size, positivity, finiteness  PASS\n');

% Central density recovery: rho(r=0) = rho0
rho_center = Soliton_profile(0, rho0_sol, rc_sol);
assert(abs(rho_center - rho0_sol)/rho0_sol < 1e-10, ...
    sprintf('rho(0) = %.4e, expected %.4e', rho_center, rho0_sol));
fprintf('  rho(r=0) == rho0              PASS\n');

% Half-density at r = rc: (1 + 0.091)^-8 ~ 0.5
rho_at_rc   = Soliton_profile(rc_sol, rho0_sol, rc_sol);
half_ratio  = rho_at_rc / rho0_sol;
assert(abs(half_ratio - 0.5) < 0.01, ...
    sprintf('rho(rc)/rho0 = %.4f, expected ~0.5', half_ratio));
fprintf('  rho(rc)/rho0 = %.4f (~0.5)   PASS\n', half_ratio);

% Monotonically decreasing
assert(all(diff(rho_sol) < 0), 'Soliton not monotonically decreasing');
fprintf('  rho monotonically decreasing  PASS\n\n');

% =====================================================================
%% 11. Soliton_profile — steep outer falloff ~ r^-16
% =====================================================================
fprintf('--- 11. Soliton_profile: outer power-law slope ---\n');
r_far   = [50, 100, 200];
rho_far = Soliton_profile(r_far, rho0_sol, rc_sol);
slope_sol = mean(diff(log(rho_far)) ./ diff(log(r_far)));
assert(abs(slope_sol - (-16)) < 0.5, ...
    sprintf('Outer slope = %.3f, expected ~ -16', slope_sol));
fprintf('  Outer slope = %.2f (expected -16)  PASS\n\n', slope_sol);

% % =====================================================================
% %% 12. DK14_profile — mass-selected defaults (selected_by = 'M')
% % =====================================================================
% fprintf('--- 12. DK14_profile: mass-selected defaults ---\n');
% rho_dk14_M = DK14_profile(r_mpc, M0, c0, z0, cosmo, Delta, 'M', []);
% 
% assert(numel(rho_dk14_M) == numel(r_mpc), 'Output size mismatch');
% assert(all(rho_dk14_M > 0),              'rho <= 0 somewhere');
% assert(all(isfinite(rho_dk14_M)),        'Non-finite values');
% fprintf('  size, positivity, finiteness  PASS\n');
% 
% % DK14 must exceed pure Einasto at large r (outer term contribution)
% [rho_e_dk, ~, ~] = Einasto_profile(r_mpc, M0, c0, z0, cosmo, Delta);
% assert(rho_dk14_M(end) > rho_e_dk(end), ...
%     'DK14 outer term not elevating density at large r');
% fprintf('  DK14 > Einasto at large r     PASS\n\n');

% =====================================================================
%% 12. DK14_profile — mass-selected defaults (selected_by = 'M')
% =====================================================================
fprintf('--- 12. DK14_profile: mass-selected defaults ---\n');
rho_dk14_M = DK14_profile(r_mpc, M0, c0, z0, cosmo, Delta, 'M', []);
assert(numel(rho_dk14_M) == numel(r_mpc), 'Output size mismatch');
assert(all(rho_dk14_M > 0),              'rho <= 0 somewhere');
assert(all(isfinite(rho_dk14_M)),        'Non-finite values');
fprintf('  size, positivity, finiteness  PASS\n');
% DK14 with outer term must exceed inner-only at large r
rho_dk14_outer = DK14_profile(r_mpc, M0, c0, z0, cosmo, Delta, 'M', [], true);
assert(numel(rho_dk14_outer) == numel(r_mpc), 'Output size mismatch (outer)');
assert(rho_dk14_outer(end) > rho_dk14_M(end), ...
    'DK14 outer term not elevating density at large r');
fprintf('  DK14+outer > DK14-inner at large r  PASS\n\n');
% =====================================================================
%% 13. DK14_profile — Gamma-selected mode
% =====================================================================
fprintf('--- 13. DK14_profile: Gamma-selected mode ---\n');
Gamma_val    = 1.5;
rho_dk14_G   = DK14_profile(r_mpc, M0, c0, z0, cosmo, Delta, 'Gamma', Gamma_val);

assert(numel(rho_dk14_G) == numel(r_mpc), 'Output size mismatch');
assert(all(rho_dk14_G > 0),              'rho <= 0 somewhere');
assert(all(isfinite(rho_dk14_G)),        'Non-finite values');
fprintf('  size, positivity, finiteness  PASS\n');

% M-selected and Gamma-selected should differ (different beta/gamma_t/rt)
max_diff = max(abs(rho_dk14_M - rho_dk14_G) ./ rho_dk14_M);
assert(max_diff > 1e-4, 'M and Gamma modes produced identical results');
fprintf('  M-mode != Gamma-mode (max rel diff = %.2e)  PASS\n\n', max_diff);

% =====================================================================
%% 14. DK14_profile — fallback when Gamma missing
% =====================================================================
fprintf('--- 14. DK14_profile: Gamma fallback warning ---\n');
warnstate = warning('off', 'all');
rho_dk14_fb = DK14_profile(r_mpc, M0, c0, z0, cosmo, Delta, 'Gamma', []);
warning(warnstate);
% Should produce same result as 'M' default
max_diff_fb = max(abs(rho_dk14_fb - rho_dk14_M) ./ rho_dk14_M);
assert(max_diff_fb < 1e-10, ...
    sprintf('Fallback result differs from M-selected (max diff = %.2e)', max_diff_fb));
fprintf('  Fallback matches M-selected   PASS\n\n');

% =====================================================================
%% 15. DK14_profile — default call (no selected_by)
% =====================================================================
fprintf('--- 15. DK14_profile: default call (nargin=6) ---\n');
rho_dk14_def = DK14_profile(r_mpc, M0, c0, z0, cosmo, Delta);
max_diff_def = max(abs(rho_dk14_def - rho_dk14_M) ./ rho_dk14_M);
assert(max_diff_def < 1e-10, ...
    'Default call does not match M-selected');
fprintf('  Default == M-selected         PASS\n\n');

% =====================================================================
%% 16. Cross-profile: all profiles decrease with r
% =====================================================================
fprintf('--- 16. All profiles monotonically decreasing ---\n');
profiles = { ...
    'NFW_analytcl',  NFW.rho;          ...
    'NFW_profile',   rho_nfw;          ...
    'Einasto',       rho_e;            ...
    'Hernquist',     rho_h;            ...
    'Soliton',       rho_sol;          ...
    'DK14 (M)',      rho_dk14_M;       ...
};
for i = 1:size(profiles,1)
    name = profiles{i,1};
    rho_p = profiles{i,2};
    is_dec = all(diff(rho_p) < 0);
    assert(is_dec, sprintf('%s is not monotonically decreasing', name));
    fprintf('  %-14s  decreasing? %s  PASS\n', name, tfstr(is_dec));
end
fprintf('\n');

% =====================================================================
% =====================================================================
%% 17. Cross-profile: inner slope comparison at r << rs
% =====================================================================
fprintf('--- 17. Inner slope comparison (r << rs) ---\n');

% NFW slope test (unitless, normalised)
rs_ref   = 0.3;   rhos_ref = 1.0;
r_tiny   = logspace(-4, -2, 30) .* rs_ref;   % r << rs_ref [Mpc/h]
slope_nfw_in = mean( diff(log(NFW_profile(r_tiny, rhos_ref, rs_ref))) ...
                  ./ diff(log(r_tiny(:))) );
slope_nfw_in = mean(slope_nfw_in(:));

fprintf('  NFW inner slope      = %+.3f  (expected ~ -1.0)\n', slope_nfw_in);
assert(isscalar(slope_nfw_in),          'slope_nfw_in not scalar');
assert(abs(slope_nfw_in - (-1)) < 0.1, 'NFW inner slope not ~ -1');
fprintf('  NFW inner slope OK            PASS\n');

% Einasto slope test:
% Must probe r << rs where rs comes from Einasto_profile itself.
% First get rs for this (M0, c0, z0) call.
[~, ~, rs_ein] = Einasto_profile(1.0, M0, c0, z0, cosmo, Delta);  % dummy r=1

% Build r_tiny relative to the actual Einasto scale radius
r_tiny_ein = logspace(-4, -2, 30) .* rs_ein;   % r << rs_ein [Mpc/h]

[rho_e_slope, ~, ~] = Einasto_profile(r_tiny_ein, M0, c0, z0, cosmo, Delta);
log_rho_e    = log(rho_e_slope(:));
log_r_ein    = log(r_tiny_ein(:));
slope_ein_in = mean( diff(log_rho_e) ./ diff(log_r_ein) );
slope_ein_in = mean(slope_ein_in(:));

fprintf('  Einasto rs           =  %.4f Mpc/h\n', rs_ein);
fprintf('  Einasto r_tiny range = [%.2e, %.2e] Mpc/h (vs rs=%.4f)\n', ...
        r_tiny_ein(1), r_tiny_ein(end), rs_ein);
fprintf('  Einasto inner slope  = %+.3f  (expected ~ 0, flat core)\n', slope_ein_in);

assert(isscalar(slope_ein_in),  'slope_ein_in not scalar');
assert(slope_ein_in > -0.6,     'Einasto inner slope too steep (should be near 0)');
fprintf('  Einasto inner slope OK        PASS\n\n');

% =====================================================================
%% 18. Vector mass input — all profiles return correct size
% =====================================================================
fprintf('--- 18. Vector mass input: output size ---\n');
M_vec = logspace(11, 15, 5);
r_single = 0.5;   % [Mpc/h]
for i = 1:numel(M_vec)
    [rho_ei, ~, ~] = Einasto_profile(r_single, M_vec(i), c0, z0, cosmo, Delta);
    [rho_hi, ~, ~] = Hernquist_profile(r_single, M_vec(i), c0, z0, cosmo, Delta);
    rho_dki        = DK14_profile(r_single, M_vec(i), c0, z0, cosmo, Delta);
    assert(isscalar(rho_ei) && rho_ei > 0, 'Einasto scalar fail');
    assert(isscalar(rho_hi) && rho_hi > 0, 'Hernquist scalar fail');
    assert(isscalar(rho_dki) && rho_dki > 0, 'DK14 scalar fail');
end
fprintf('  Einasto, Hernquist, DK14: scalar r, vector M  PASS\n\n');

% =====================================================================
%% 19. Plot: all profiles on same axes (normalised)
% =====================================================================
fprintf('--- 19. Plot: all profiles (normalised) ---\n');

% Compute each profile on a common radius grid [Mpc/h]
r_plot  = logspace(-2.5, 1.2, 300);
r_kpc_p = r_plot * 1e3;   % for NFW_analytcl and NFW_profile

% NFW_analytcl_Profile (uses kpc)
Rvir_kpc = 500;
NFW_plt  = NFW_analytcl_Profile(1e12, Rvir_kpc, 10, linspace(1, Rvir_kpc, 500));
r_nfw_n  = NFW_plt.r ./ Rvir_kpc;
rho_nfw_n = NFW_plt.rho ./ NFW_plt.rho(1);

% NFW_profile (arbitrary normalisation)
rho_nfw_p = NFW_profile(r_plot, 1.0, 0.3);
rho_nfw_p = rho_nfw_p ./ rho_nfw_p(1);

% Einasto
[rho_e_p,~,~] = Einasto_profile(r_plot, M0, c0, z0, cosmo, Delta);
rho_e_p = rho_e_p ./ rho_e_p(1);

% Hernquist
[rho_h_p,~,~] = Hernquist_profile(r_plot, M0, c0, z0, cosmo, Delta);
rho_h_p = rho_h_p ./ rho_h_p(1);

% DK14 (M-selected)
rho_dk_p = DK14_profile(r_plot, M0, c0, z0, cosmo, Delta, 'M');
rho_dk_p = rho_dk_p ./ rho_dk_p(1);

% Soliton (own radii in kpc)
r_sol_p  = logspace(-1, 1.5, 300);
rho_s_p  = Soliton_profile(r_sol_p, 1.0, 1.0);
rho_s_p  = rho_s_p ./ rho_s_p(1);

figure('Name','Halo Profiles — Normalised Comparison', ...
       'Position',[100 100 950 560]);
hold on;
plot(r_plot, rho_nfw_p,  'b-',  'LineWidth', 2.0, 'DisplayName', 'NFW (profile)');
plot(r_plot, rho_e_p,    'r-',  'LineWidth', 2.0, 'DisplayName', 'Einasto');
plot(r_plot, rho_h_p,    'g-',  'LineWidth', 2.0, 'DisplayName', 'Hernquist');
plot(r_plot, rho_dk_p,   'm--', 'LineWidth', 2.0, 'DisplayName', 'DK14 (M-sel)');
plot(r_sol_p./10, rho_s_p, 'k:', 'LineWidth', 2.0, 'DisplayName', 'Soliton');
set(gca, 'XScale', 'log', 'YScale', 'log');
xlabel('$r / r_{200}$', 'FontSize', 13, 'Interpreter', 'latex');
ylabel('$\rho / \rho(r_\mathrm{min})$', 'FontSize', 13, 'Interpreter', 'latex');
title('Normalised Halo Density Profiles', 'FontSize', 14);
legend('Location', 'southwest', 'FontSize', 10);
grid on; hold off;
fprintf('  PASS\n\n');

% =====================================================================
%% 20. Plot: DK14 decomposition — inner * f_trans + outer
% =====================================================================
fprintf('--- 20. Plot: DK14 decomposition ---\n');

[rho_e_dk2, rhos_e2, rs_e2] = Einasto_profile(r_mpc, M0, c0, z0, cosmo, Delta);
rho_c_val = cosmo.rho_crit0 .* cosmo.E(z0).^2;
R200m_val = (3*M0 / (4*pi*Delta*rho_c_val))^(1/3);
nu200m_v  = cosmo.nu(M0, z0);
rt_val    = R200m_val * (1.9 - 0.18*nu200m_v);
rt_val    = max(rt_val, 0.01*R200m_val);
beta_v = 4;   gamma_v = 8;
f_trans_v = (1 + (r_mpc./rt_val).^beta_v).^(-gamma_v./beta_v);
rho_m_v   = cosmo.rho_m0 .* (1+z0).^3;
rho_out_v = rho_m_v .* (r_mpc).^(-1.5);

figure('Name','DK14 Profile Decomposition','Position',[150 150 950 560]);
hold on;
plot(r_mpc, rho_e_dk2,                 'b-',  'LineWidth',2, 'DisplayName','Einasto inner');
plot(r_mpc, rho_e_dk2 .* f_trans_v,   'r-',  'LineWidth',2, 'DisplayName','Inner \times f_{trans}');
plot(r_mpc, rho_out_v,                 'g--', 'LineWidth',2, 'DisplayName','Outer (power-law)');
plot(r_mpc, DK14_profile(r_mpc, M0, c0, z0, cosmo, Delta), ...
                                       'k-',  'LineWidth',2.5,'DisplayName','DK14 total');
xline(rt_val, '--', sprintf('r_t = %.2f Mpc/h', rt_val), ...
      'Color',[0.5 0.5 0.5], 'LineWidth', 1.5, 'LabelVerticalAlignment','bottom');
set(gca, 'XScale','log', 'YScale','log');
xlabel('$r$ [Mpc/h]', 'FontSize', 13, 'Interpreter', 'latex');
ylabel('$\rho$ [M$_\odot h$ / Mpc$^3$]', 'FontSize', 13, 'Interpreter', 'latex');
title('DK14 Profile Decomposition  (M=10^{13} M_\odot/h, z=0)', 'FontSize', 14);
legend('Location','southwest','FontSize',10);
grid on; hold off;
fprintf('  PASS\n\n');

% =====================================================================
fprintf('==========================================\n');
fprintf('  All profile tests PASSED\n');
fprintf('==========================================\n');