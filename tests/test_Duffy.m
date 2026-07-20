%% test_Duffy08.m
% Test script for Duffy08_concentration and Duffy08_Table
%
% Location : PHANTOM\tests\test_Duffy08.m
% Run from : anywhere — paths are anchored to this file's location
%
% Tests covered:
%   1. All 24 table entries load without error
%   2. NFW and Einasto concentrations evaluated over mass-redshift grid
%   3. Sanity checks (c decreasing with M, c decreasing with z)
%   4. Profile comparison: NFW vs Einasto at z=0
%   5. Sample comparison: full vs relaxed at z=0
%   6. Redshift-range comparison: z0 vs z0_2 at z=0
%   7. Custom pivot mass warning and rescaling check
%   8. Plots

%clear; clc; close all;

%% =========================================================================
% 0. Paths
% =========================================================================
here     = fileparts(mfilename('fullpath'));   % PHANTOM\tests
repo     = fullfile(here, '..');              % PHANTOM\
utils_path = fullfile(repo, 'src', 'utils');
concentration_path = fullfile(repo, 'src', 'concentration');
addpath(genpath(utils_path));
addpath(genpath(concentration_path));

assert(exist('Duffy08', 'file') == 2, ...
    sprintf('Duffy08.m not found.\nExpected in: %s', utils_path));
assert(exist('Duffy08_Table', 'file') == 2, ...
    sprintf('Duffy08_Table.m not found.\nExpected in: %s', utils_path));

fprintf('Paths OK.\n');
fprintf('  Duffy08_concentration : %s\n', which('Duffy08_concentration'));
fprintf('  Duffy08_Table         : %s\n', which('Duffy08_Table'));

%% =========================================================================
% 1. Table integrity — every valid combination must load without error
% =========================================================================
fprintf('\n--- Table integrity check (all 24 entries) ---\n');

mdefs    = {'200c', 'vir', '200m'};
profiles = {'NFW', 'Einasto'};
samples  = {'full', 'relaxed'};
zranges  = {'z0', 'z0_2'};

n_ok  = 0;
n_err = 0;
for im = 1:numel(mdefs)
for ip = 1:numel(profiles)
for is = 1:numel(samples)
for iz = 1:numel(zranges)
    try
        P = Duffy08_Table(mdefs{im}, profiles{ip}, samples{is}, zranges{iz});
        assert(isfield(P,'A') && isfield(P,'B') && isfield(P,'C') && isfield(P,'M_pivot'), ...
            'Missing field in output struct.');
        n_ok = n_ok + 1;
    catch ME
        fprintf('  FAILED: %s | %s | %s | %s\n  -> %s\n', ...
            mdefs{im}, profiles{ip}, samples{is}, zranges{iz}, ME.message);
        n_err = n_err + 1;
    end
end; end; end; end

fprintf('  %d / 24 entries loaded OK,  %d failed.\n', n_ok, n_err);
assert(n_err == 0, 'One or more table entries failed to load.');

%% =========================================================================
% 2. Mass and redshift grids
% =========================================================================
M_vec = logspace(11, 15, 50);   % [Msun/h]
z_vec = [0.0, 0.5, 1.0, 2.0];

%% =========================================================================
% 3. Concentration over mass-redshift grid (NFW, full, z0_2, 200c)
% =========================================================================
fprintf('\n--- c(M,z) grid: NFW | full | z0_2 | 200c ---\n');

c_nfw = zeros(numel(z_vec), numel(M_vec));
for iz = 1:numel(z_vec)
    c_nfw(iz,:) = Duffy08(M_vec, z_vec(iz), '200c');
    fprintf('  z=%.1f  |  c range: [%.2f, %.2f]\n', ...
        z_vec(iz), min(c_nfw(iz,:)), max(c_nfw(iz,:)));
end

%% =========================================================================
% 4. Sanity checks
% =========================================================================
fprintf('\n--- Sanity checks ---\n');

% 4a. c decreasing with M at every z
for iz = 1:numel(z_vec)
    ok = all(diff(c_nfw(iz,:)) < 0);
    fprintf('  z=%.1f: c decreasing with M? %s\n', z_vec(iz), mat2str(ok));
end

% 4b. c decreasing with z at fixed M
i_mid = round(numel(M_vec)/2);
fprintf('\n  c vs z at M = 10^%.1f Msun/h:\n', log10(M_vec(i_mid)));
for iz = 1:numel(z_vec)
    fprintf('    z=%.1f  =>  c = %.3f\n', z_vec(iz), c_nfw(iz, i_mid));
end
ok_z = all(diff(c_nfw(:, i_mid)) < 0);
fprintf('  c decreasing with z? %s\n', mat2str(ok_z));

% 4c. All values positive and finite
ok_finite = all(isfinite(c_nfw(:))) && all(c_nfw(:) > 0);
fprintf('\n  All concentrations positive and finite? %s\n', mat2str(ok_finite));

% 4d. Plausible range at z=0: roughly 3-30 for 10^11-10^15 Msun/h
c_z0     = c_nfw(1,:);
in_range = all(c_z0 > 1 & c_z0 < 50);
fprintf('  z=0 values in plausible range [1,50]? %s\n', mat2str(in_range));

%% =========================================================================
% 5. Profile comparison: NFW vs Einasto (full, z0_2, 200c, z=0)
% =========================================================================
fprintf('\n--- Profile comparison: NFW vs Einasto (200c, full, z0_2, z=0) ---\n');

c_nfw_z0     = Duffy08(M_vec, 0, '200c', 'NFW',     'full', 'z0_2');
c_einasto_z0 = Duffy08(M_vec, 0, '200c', 'Einasto', 'full', 'z0_2');

for logM = [12, 13, 14]
    cn = interp1(log10(M_vec), c_nfw_z0,     logM);
    ce = interp1(log10(M_vec), c_einasto_z0, logM);
    fprintf('  M=1e%d:  NFW=%.2f  |  Einasto=%.2f  |  ratio=%.3f\n', ...
        logM, cn, ce, ce/cn);
end

%% =========================================================================
% 6. Sample comparison: full vs relaxed (NFW, z0_2, 200c, z=0)
% =========================================================================
fprintf('\n--- Sample comparison: full vs relaxed (NFW, 200c, z0_2, z=0) ---\n');

c_full    = Duffy08(M_vec, 0, '200c', 'NFW', 'full',    'z0_2');
c_relaxed = Duffy08(M_vec, 0, '200c', 'NFW', 'relaxed', 'z0_2');

% Relaxed haloes should be MORE concentrated than full sample
for logM = [12, 13, 14]
    cf = interp1(log10(M_vec), c_full,    logM);
    cr = interp1(log10(M_vec), c_relaxed, logM);
    ok = cr > cf;
    fprintf('  M=1e%d:  full=%.2f  |  relaxed=%.2f  |  relaxed>full? %s\n', ...
        logM, cf, cr, mat2str(ok));
end

%% =========================================================================
% 7. Redshift-range comparison: z0 vs z0_2 exactly at z=0
%    At z=0, (1+z)^C = 1 so z0_2 and z0 should give similar but not
%    identical results (different A values from separate fits)
% =========================================================================
fprintf('\n--- z-range comparison: z0 vs z0_2 at z=0 (NFW, 200c, full) ---\n');

c_z0fit  = Duffy08(M_vec, 0, '200c', 'NFW', 'full', 'z0');
c_z02fit = Duffy08(M_vec, 0, '200c', 'NFW', 'full', 'z0_2');

for logM = [12, 13, 14]
    cz0  = interp1(log10(M_vec), c_z0fit,  logM);
    cz02 = interp1(log10(M_vec), c_z02fit, logM);
    fprintf('  M=1e%d:  z0-fit=%.3f  |  z0_2-fit=%.3f  |  diff=%.2f%%\n', ...
        logM, cz0, cz02, 100*abs(cz0-cz02)/cz0);
end

%% =========================================================================
% 8. Custom pivot mass — check rescaling consistency
%    With a rescaled pivot, c(M=M_new_pivot) should equal A_new (by definition)
% =========================================================================
fprintf('\n--- Custom pivot mass check ---\n');

M_new = 1e14;   % alternative pivot

% Suppress the expected warning for this test
warning('off', 'all');
c_std    = Duffy08(M_new, 0, '200c', 'NFW', 'full', 'z0_2');
c_custom = Duffy08(M_new, 0, '200c', 'NFW', 'full', 'z0_2', M_new);
warning('on',  'all');

% Both calls evaluated at M = M_new_pivot; they should return the same c
% because A is rescaled to maintain the power law
fprintf('  c at M=1e14 with default pivot (2e12): %.4f\n', c_std);
fprintf('  c at M=1e14 with custom pivot (1e14) : %.4f\n', c_custom);
fprintf('  Difference: %.2e  (should be < 1e-10)\n', abs(c_std - c_custom));
assert(abs(c_std - c_custom) < 1e-8, ...
    'Custom pivot rescaling is inconsistent — check Duffy08_concentration.');

%% =========================================================================
% 9. Plot: c(M) — NFW vs Einasto at multiple redshifts
% =========================================================================
figure('Name','Duffy08: NFW vs Einasto c(M)','Position',[100 100 900 550]);

colors = lines(numel(z_vec));
ax1 = subplot(1,2,1);
hold on;
for iz = 1:numel(z_vec)
    cn = Duffy08(M_vec, z_vec(iz), '200c', 'NFW',     'full', 'z0_2');
    ce = Duffy08(M_vec, z_vec(iz), '200c', 'Einasto', 'full', 'z0_2');
    plot(log10(M_vec), cn, '-',  'Color', colors(iz,:), 'LineWidth', 2, ...
        'DisplayName', sprintf('NFW z=%.1f',     z_vec(iz)));
    plot(log10(M_vec), ce, '--', 'Color', colors(iz,:), 'LineWidth', 1.5, ...
        'DisplayName', sprintf('Einasto z=%.1f', z_vec(iz)));
end
hold off;
xlabel('log_{10}(M  [M_\odot/h])', 'FontSize', 12);
ylabel('Concentration  c', 'FontSize', 12);
title('NFW (—) vs Einasto (- -), M_{200c}', 'FontSize', 13);
legend('Location','northeast','FontSize',8);
grid on;

% =========================================================================
% 10. Plot: full vs relaxed, all three mass definitions at z=0
% =========================================================================
ax2 = subplot(1,2,2);
mdef_list  = {'200c','vir','200m'};
ls_full    = {'-', '--', ':'};
colors2    = lines(numel(mdef_list));
hold on;
for im = 1:numel(mdef_list)
    cf = Duffy08(M_vec, 0, mdef_list{im}, 'NFW', 'full',    'z0_2');
    cr = Duffy08(M_vec, 0, mdef_list{im}, 'NFW', 'relaxed', 'z0_2');
    plot(log10(M_vec), cf, ls_full{im}, 'Color', colors2(im,:), 'LineWidth', 2, ...
        'DisplayName', sprintf('%s full',    mdef_list{im}));
    plot(log10(M_vec), cr, ls_full{im}, 'Color', colors2(im,:), 'LineWidth', 1.2, ...
        'Marker', 'none', 'LineStyle', '-.', ...
        'DisplayName', sprintf('%s relaxed', mdef_list{im}));
end
hold off;
xlabel('$log_{10}(M  [M_{\odot}/h]$)', 'FontSize', 12, 'Interpreter','latex');
ylabel('Concentration  c', 'FontSize', 12);
title('Full (—) vs Relaxed (-.), NFW, z=0', 'FontSize', 13);
legend('Location','northeast','FontSize',8);
grid on;

sgtitle('Duffy et al. 2008 Concentration Model', 'FontSize', 14, 'FontWeight','bold');

fprintf('\nAll Duffy08 tests passed.\n');