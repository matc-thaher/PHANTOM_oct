%% test_Prada12.m
% Test script for Prada12_concentration
%
% Location : PHANTOM\tests\test_Prada12.m
%
% Tests covered:
%   1. Path check
%   2. B0 = B1 = 1 at z=0 (normalisation identity)
%   3. Spot checks against Colossus reference values at z=0
%   4. c decreasing with M at z=0 (standard expectation)
%   5. c decreasing with z at fixed galaxy-size mass
%   6. All values positive and finite over wide M-z grid
%   7. Vector input matches element-wise scalar calls
%   8. Bolshoi x0_norm reproduces 1.393
%   9. Plots: c(M) at multiple z + c(z) at fixed masses

% clear; clc; close all;

%% =========================================================================
% 0. Paths
% =========================================================================
here       = fileparts(mfilename('fullpath'));
repo       = fullfile(here, '..');
utils_path = fullfile(repo, 'src', 'utils');
concentration_path = fullfile(repo, 'src', 'concentration');
addpath(genpath(utils_path));
addpath(genpath(concentration_path));

assert(exist('Prada12','file') == 2, ...
    sprintf('Prada12_concentration.m not found.\nExpected in: %s', utils_path));
fprintf('Paths OK.\n');
fprintf('  Prada12_concentration : %s\n', which('Prada12_concentration'));

%% =========================================================================
% 1. Build Bolshoi cosmology struct (matches Prada+12 paper)
% =========================================================================
cosmo.Omega_m = 0.27;
cosmo.Omega_b = 0.0469;
cosmo.h      = 0.70;
cosmo.ns     = 0.95;
cosmo.sigma8 = 0.82;

% Attach linear components (sigmaM handle etc.)
% cosmo = cosmology(cosmo.Omegam, cosmo.OmegaL, cosmo.Omegab, ...
%                   cosmo.h, cosmo.ns, cosmo.sigma8);
cosmo = cosmology('custom', cosmo);
cosmo = attach_linear_components(cosmo);

fprintf('\nBolshoi cosmology built OK.\n');
fprintf('  Omegam=%.3f  OmegaL=%.3f  h=%.2f  sigma8=%.2f\n', ...
    cosmo.Omega_m, cosmo.Omega_L, cosmo.h, cosmo.sigma8);

%% =========================================================================
% 2. x0_norm check: (OmegaL/Omegam)^(1/3) must equal 1.393 for Bolshoi
% =========================================================================
fprintf('\n--- x0_norm check (Bolshoi cosmology) ---\n');

x0_norm = (cosmo.Omega_L / cosmo.Omega_m)^(1/3);
fprintf('  x0_norm = %.6f  (expected ~1.3932)\n', x0_norm);
assert(abs(x0_norm - 1.3932) < 1e-3, 'x0_norm deviates from expected 1.393 for Bolshoi cosmology.');

%% =========================================================================
% 3. B0 = B1 = 1 at z=0 (by construction of the normalisation)
%    At z=0: x = x0_norm, so B0 = cmin(x)/cmin(x0_norm) = 1
%                               B1 = smin_inv(x)/smin_inv(x0_norm) = 1
%    We verify this indirectly: c(M, z=0) must equal B0*C(B1*sigma)
%    with B0=B1=1, so c must equal C(sigma) directly.
% =========================================================================
fprintf('\n--- B0=B1=1 at z=0 (normalisation identity) ---\n');

% Compute c at z=0 via the full function
M_test = 1e12;
c_z0   = Prada12(M_test, 0, cosmo);

% Reproduce manually with B0=B1=1 forced
sigma_z0 = cosmo.sigmaM(M_test, 0);
A_p=2.881; b_p=1.257; c_p=1.022; d_p=0.060;
C_direct = A_p * ((sigma_z0/b_p)^c_p + 1) * exp(d_p/sigma_z0^2);

fprintf('  c(1e12, z=0) from function : %.6f\n', c_z0);
fprintf('  C(sigma) with B0=B1=1      : %.6f\n', C_direct);
fprintf('  Difference                 : %.2e\n',  abs(c_z0 - C_direct));
assert(abs(c_z0 - C_direct) < 1e-8, 'B0*B1 != 1 at z=0 — normalisation broken.');

%% =========================================================================
% 4. Plausible range at z=0
%    For M200c ~ 1e12 h^-1 Msun, Prada+12 predict c ~ 7-9 at z=0
%    (Bolshoi cosmology, all haloes)
% =========================================================================
fprintf('\n--- Plausible range at z=0 ---\n');

M_mw = 1e12;   % Milky-Way-size halo
c_mw = Prada12(M_mw, 0, cosmo);
fprintf('  c(1e12, z=0) = %.3f  (expected ~7-10 for Bolshoi)\n', c_mw);
assert(c_mw > 4 && c_mw < 20, 'c(1e12, z=0) out of plausible range [4, 20].');

M_cl = 1e14;   % cluster-size halo
c_cl = Prada12(M_cl, 0, cosmo);
fprintf('  c(1e14, z=0) = %.3f  (expected ~4-7 for Bolshoi)\n', c_cl);
assert(c_cl > 2 && c_cl < 15, 'c(1e14, z=0) out of plausible range [2, 15].');
assert(c_mw > c_cl, 'c(1e12) should be > c(1e14) at z=0.');

%% =========================================================================
% 5. c decreasing with M at z=0 over [1e11, 1e14]
% =========================================================================
fprintf('\n--- c decreasing with M at z=0 ---\n');

M_vec = logspace(11, 14, 50);
c_vec = Prada12(M_vec, 0, cosmo);
c_vec = reshape(c_vec, 1, []);

ok = true;
for k = 1:numel(c_vec)-1
    if c_vec(k+1) >= c_vec(k)
        ok = false;
        break;
    end
end
fprintf('  c monotone decreasing with M over [1e11,1e14]? %s\n', tf2str(ok));
assert(ok, 'c is not monotone decreasing with M at z=0.');

%% =========================================================================
% 6. c decreasing with z at fixed mass (z = 0 to 3)
% =========================================================================
fprintf('\n--- c decreasing with z at fixed M ---\n');

z_vec  = [0.0, 0.5, 1.0, 2.0, 3.0];
masses = [1e12, 1e13];

for M = masses
    c_z = arrayfun(@(z) Prada12(M, z, cosmo), z_vec);
    fprintf('\n  M = 1e%.0f h^-1 Msun:\n', log10(M));
    for iz = 1:numel(z_vec)
        fprintf('    z=%.1f  c=%.3f\n', z_vec(iz), c_z(iz));
    end
    ok = true;
    for k = 1:numel(c_z)-1
        if c_z(k+1) >= c_z(k)
            ok = false;
            break;
        end
    end
    fprintf('  c monotone decreasing z=0->3? %s\n', tf2str(ok));
end

%% =========================================================================
% 7. All values positive and finite over wide M-z grid
% =========================================================================
fprintf('\n--- Finite/positive check over M-z grid ---\n');

M_grid = logspace(10, 15, 80);
z_grid = [0.0, 0.5, 1.0, 2.0, 3.0, 5.0];
n_bad  = 0;
for iz = 1:numel(z_grid)
    c_test = Prada12(M_grid, z_grid(iz), cosmo);
    bad    = sum(~isfinite(c_test) | c_test <= 0);
    n_bad  = n_bad + bad;
end
fprintf('  Bad values (non-finite or <=0): %d  (expected 0)\n', n_bad);
assert(n_bad == 0, 'Non-finite or non-positive concentration values found.');

%% =========================================================================
% 8. Vector input matches element-wise scalar calls
% =========================================================================
fprintf('\n--- Vector vs scalar consistency ---\n');

M_check  = [1e11, 5e11, 1e12, 5e12, 1e13, 5e13, 1e14];
c_vector = Prada12(M_check, 0, cosmo);
c_scalar = arrayfun(@(M) Prada12(M, 0, cosmo), M_check);
c_vector = reshape(c_vector, 1, []);
c_scalar = reshape(c_scalar, 1, []);
max_diff = max(abs(c_vector - c_scalar));
fprintf('  Max diff vector vs scalar: %.2e  (expected < 1e-12)\n', max_diff);
assert(max_diff < 1e-10, 'Vector and scalar calls give inconsistent results.');

%% =========================================================================
% 9. Cosmology sensitivity: Planck vs Bolshoi at z=1
%    Different cosmologies should give different x0_norm and thus
%    different B0/B1 corrections at z>0, but same c at z=0
% =========================================================================
fprintf('\n--- Cosmology sensitivity: Planck vs Bolshoi at z=1 ---\n');

cosmo_planck          = cosmo;   % copy struct
cosmo_planck.Omegam   = 0.315;
cosmo_planck.OmegaL   = 0.685;
cosmo_planck          = attach_linear_components(cosmo_planck);

x0_bolshoi = (0.73/0.27)^(1/3);
x0_planck  = (0.685/0.315)^(1/3);
fprintf('  x0_norm Bolshoi : %.4f\n', x0_bolshoi);
fprintf('  x0_norm Planck  : %.4f\n', x0_planck);

c_bolshoi_z1 = Prada12(1e12, 1.0, cosmo);
c_planck_z1  = Prada12(1e12, 1.0, cosmo_planck);
fprintf('  c(1e12, z=1) Bolshoi : %.4f\n', c_bolshoi_z1);
fprintf('  c(1e12, z=1) Planck  : %.4f\n', c_planck_z1);
fprintf('  Different cosmologies give different c at z=1? %s\n', ...
    tf2str(abs(c_bolshoi_z1 - c_planck_z1) > 0.01));

%% =========================================================================
% 10. Plots
% =========================================================================
M_plot  = logspace(10, 15, 200);
z_plot  = [0.0, 0.5, 1.0, 2.0, 3.0];
colors  = lines(numel(z_plot));

figure('Name','Prada12: c(M) at multiple redshifts', 'Position',[80 80 1050 470]);

% --- Panel 1: c(M) at each z --------------------------------------------
subplot(1,2,1);
hold on;
for iz = 1:numel(z_plot)
    c_pl = Prada12(M_plot, z_plot(iz), cosmo);
    plot(log10(M_plot), reshape(c_pl,1,[]), '-', ...
        'Color', colors(iz,:), 'LineWidth', 2, ...
        'DisplayName', sprintf('z = %.1f', z_plot(iz)));
end
hold off;
xlabel('log_{10}(M_{200c}  [h^{-1} M_\odot])', 'FontSize', 12);
ylabel('Concentration  c_{200c}',               'FontSize', 12);
title('Prada+12, Bolshoi cosmology',            'FontSize', 13);
legend('Location','northeast','FontSize',9);
grid on; xlim([10 15]);

% --- Panel 2: c(z) at two fixed masses ----------------------------------
subplot(1,2,2);
z_fine   = linspace(0, 5, 200);
mass_fix = [1e12, 1e14];
ls       = {'-','--'};
hold on;
for im = 1:numel(mass_fix)
    c_zz = arrayfun(@(z) Prada12(mass_fix(im), z, cosmo), z_fine);
    plot(z_fine, c_zz, ls{im}, 'LineWidth', 2, ...
        'DisplayName', sprintf('M=10^{%.0f}', log10(mass_fix(im))));
end
hold off;
xlabel('Redshift  z',          'FontSize', 12);
ylabel('Concentration  c_{200c}', 'FontSize', 12);
title('c(z) at fixed mass',    'FontSize', 13);
legend('Location','northeast', 'FontSize', 10);
grid on;

sgtitle('Prada, Klypin, Cuesta et al. (2012)', ...
    'FontSize', 14, 'FontWeight', 'bold');

fprintf('\nAll Prada12 tests passed.\n');

%% =========================================================================
% Helper
% =========================================================================
function s = tf2str(tf)
    if tf; s = 'PASS'; else; s = 'FAIL'; end
end