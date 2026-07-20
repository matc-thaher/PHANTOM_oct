%% test_Bullock01_concentration.m
% Test script for Bullock01_concentration
%
% Builds a cosmology struct using the utils functions, then evaluates
% the Bullock01 concentration over a mass and redshift grid.
%
% Folder structure assumed:
%   src/utils/   — all utility functions (cosmology, sigma, neff, etc.)
%   src/concentration/  — concentration models (not needed for Bullock01)
%
% Run from the repo root:
%   >> cd /path/to/repo
%   >> test_Bullock01_concentration

%clear; clc; close all;

%% =========================================================================
% 0. Add paths
% =========================================================================
here     = fileparts(mfilename('fullpath'));        % F:\PHANTOM\tests
repo     = fullfile(here, '..');                   % F:\PHANTOM
utils_path         = fullfile(repo, 'src', 'utils');
concentration_path = fullfile(repo, 'src', 'concentration');

addpath(genpath(utils_path));
addpath(genpath(concentration_path));
%% =========================================================================
% 1. Build cosmology struct
%    cosmology() sets Omega_m, Omega_b, h, n_s, sigma8, etc.
%    attach_linear_components() attaches:
%      cosmo.E(z)          — dimensionless Hubble rate E(z) = H(z)/H0
%      cosmo.D(z)          — linear growth factor D(z), normalized D(0)=1
%      cosmo.T(k)          — Eisenstein & Hu 1998 transfer function
%      cosmo.Pk(k,z)       — linear matter power spectrum
%      cosmo.sigmaR(R,z)   — rms variance smoothed at radius R [Mpc/h]
%      cosmo.sigmaM(M,z)   — rms variance at Lagrangian radius of mass M
%      cosmo.neff(M,z,k)   — effective spectral index at kappa*R_L(M)
%      cosmo.alphaEff(z)   — effective growth exponent d ln D / d ln(1+z)
% =========================================================================
cosmo = cosmology('Planck18');   % Planck 2018 parameters

fprintf('Cosmology: Omega_m=%.3f, h=%.3f, sigma8=%.3f\n', ...
    cosmo.Omega_m, cosmo.h, cosmo.sigma8);

%% =========================================================================
% 2. Define mass and redshift grids
% =========================================================================
M_vec = logspace(11, 15, 40);   % halo masses [Msun/h],  10^11 to 10^15
z_vec = [0.0, 0.5, 1.0, 2.0];  % redshifts to compare

%% =========================================================================
% 3. Evaluate Bullock01 concentration
%    Default parameters: K = 2.9, F = 0.001
%    (Johnston et al. 2007 / Colossus recalibration)
% =========================================================================
fprintf('\nComputing Bullock01 concentrations...\n');

c_B01 = zeros(numel(z_vec), numel(M_vec));

for iz = 1:numel(z_vec)
    z = z_vec(iz);
    c_B01(iz, :) = Bullock01(M_vec, z, cosmo);
    fprintf('  z = %.1f  done  (c range: %.2f -- %.2f)\n', ...
        z, min(c_B01(iz,:)), max(c_B01(iz,:)));
end

%% =========================================================================
% 4. Sanity checks
% =========================================================================
fprintf('\n--- Sanity checks ---\n');

% 4a. Concentration should decrease with increasing mass (at fixed z)
for iz = 1:numel(z_vec)
    is_decreasing = all(diff(c_B01(iz,:)) < 0);
    fprintf('  z=%.1f: c decreasing with M? %s\n', ...
        z_vec(iz), mat2str(is_decreasing));
end

% 4b. Concentration should decrease with increasing redshift (at fixed M)
% Pick a mid-range mass index
i_M = round(numel(M_vec)/2);
fprintf('\n  c vs z at M = 10^%.1f Msun/h:\n', log10(M_vec(i_M)));
for iz = 1:numel(z_vec)
    fprintf('    z=%.1f  =>  c = %.3f\n', z_vec(iz), c_B01(iz, i_M));
end
is_decreasing_z = all(diff(c_B01(:, i_M)) < 0);
fprintf('  c decreasing with z? %s\n', mat2str(is_decreasing_z));

% 4c. Typical concentration values at z=0 should be in range ~4-20
c_z0 = c_B01(1, :);
in_range = all(c_z0 > 2 & c_z0 < 50);
fprintf('\n  z=0 concentrations in plausible range [2,50]? %s\n', ...
    mat2str(in_range));
fprintf('  z=0: c(1e12) = %.2f,  c(1e14) = %.2f\n', ...
    interp1(log10(M_vec), c_z0, 12), ...
    interp1(log10(M_vec), c_z0, 14));

%% =========================================================================
% 5. Parameter sensitivity: compare default vs original Bullock+01 values
% =========================================================================
fprintf('\n--- Parameter sensitivity (z=0) ---\n');

c_default  = Bullock01(M_vec, 0, cosmo);          % K=2.9, F=0.001
c_original = Bullock01(M_vec, 0, cosmo, 4.0, 0.01); % B01 original

fprintf('  At M=1e12 Msun/h:\n');
fprintf('    Default  (K=2.9, F=0.001): c = %.3f\n', ...
    interp1(log10(M_vec), c_default,  12));
fprintf('    Original (K=4.0, F=0.01) : c = %.3f\n', ...
    interp1(log10(M_vec), c_original, 12));

%% =========================================================================
% 6. Plot: c(M) at multiple redshifts
% =========================================================================
figure('Name', 'Bullock01 Concentration', 'Position', [100 100 800 500]);

colors = lines(numel(z_vec));
hold on;
for iz = 1:numel(z_vec)
    plot(log10(M_vec), c_B01(iz,:), ...
        'Color', colors(iz,:), 'LineWidth', 2, ...
        'DisplayName', sprintf('z = %.1f', z_vec(iz)));
end

% Overlay the original B01 parameters at z=0 for comparison
plot(log10(M_vec), c_original, 'k--', 'LineWidth', 1.5, ...
    'DisplayName', 'z=0, K=4.0 F=0.01 (B01 original)');

hold off;
xlabel('$log_{10}( M  [M_{\odot} / h] $)', 'FontSize', 13, 'Interpreter','latex');
ylabel('Concentration  c', 'FontSize', 13);
title('Bullock01 Concentration-Mass Relation', 'FontSize', 14);
legend('Location', 'northeast', 'FontSize', 11, 'Interpreter','latex');
grid on;
set(gca, 'FontSize', 12);

%% =========================================================================
% 7. Plot: c(z) at fixed masses
% =========================================================================
z_fine  = linspace(0, 3, 60);
M_fixed = [1e12, 1e13, 1e14];   % [Msun/h]

c_zfine = zeros(numel(M_fixed), numel(z_fine));
for iM = 1:numel(M_fixed)
    for iz = 1:numel(z_fine)
        c_zfine(iM, iz) = Bullock01(M_fixed(iM), z_fine(iz), cosmo);
    end
end

figure('Name', 'Bullock01 c(z)', 'Position', [150 150 800 500]);
colors2 = lines(numel(M_fixed));
hold on;
for iM = 1:numel(M_fixed)
    plot(z_fine, c_zfine(iM,:), ...
        'Color', colors2(iM,:), 'LineWidth', 2, ...
        'DisplayName', sprintf('$M = 10^{%.0f} {M_{\\odot}/h}$', log10(M_fixed(iM))));
end
hold off;
xlabel('Redshift  z', 'FontSize', 13);
ylabel('Concentration  c', 'FontSize', 13);
title('Bullock01: Concentration vs Redshift', 'FontSize', 14);
legend('Location', 'northeast', 'FontSize', 11, 'Interpreter','latex');
grid on;
set(gca, 'FontSize', 12);

fprintf('\nTest complete.\n');