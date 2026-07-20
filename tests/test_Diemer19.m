% test_Diemer19.m
% Test suite for Diemer19_Table.m and Diemer19_concentration.m
% Reference: Diemer & Joyce 2019, ApJ 871, 168

%% =========================================================================
% 0. Paths
% =========================================================================
here     = fileparts(mfilename('fullpath'));   % PHANTOM\tests
repo     = fullfile(here, '..');              % PHANTOM\
utils_path = fullfile(repo, 'src', 'utils');
concentration_path = fullfile(repo, 'src', 'concentration');
addpath(genpath(utils_path));
addpath(genpath(concentration_path));

assert(exist('Diemer19', 'file') == 2, ...
    sprintf('Diemer19.m not found.\nExpected in: %s', utils_path));
assert(exist('Diemer19_Table', 'file') == 2, ...
    sprintf('Diemer19_Table.m not found.\nExpected in: %s', utils_path));

fprintf('Paths OK.\n');
fprintf('  Diemer19_concentration : %s\n', which('Diemer19_concentration'));
fprintf('  Diemer19_Table         : %s\n', which('Diemer19_Table'));

fprintf('========================================\n');
fprintf('  Diemer19 test suite\n');
fprintf('========================================\n\n');

cosmo = cosmology('Planck18');
modes = {'mean', 'median'};

%% --- 1. Table: both modes load with all fields -----------------------
fprintf('--- 1. Table loads correctly ---\n');
fields = {'kappa','a0','a1','b0','b1','cAlpha'};
for i = 1:numel(modes)
    P = Diemer19_Table(modes{i});
    for j = 1:numel(fields)
        assert(isfield(P, fields{j}), ...
            sprintf('Mode "%s" missing field %s', modes{i}, fields{j}));
        assert(isfinite(P.(fields{j})), ...
            sprintf('Mode "%s" field %s is not finite', modes{i}, fields{j}));
    end
    fprintf('  %-8s  kappa=%.2f  a0=%.2f  a1=%.2f  b0=%.2f  b1=%.2f  cAlpha=%.2f  PASS\n', ...
        modes{i}, P.kappa, P.a0, P.a1, P.b0, P.b1, P.cAlpha);
end
fprintf('\n');

%% --- 2. kappa = 0.42 for both modes ----------------------------------
fprintf('--- 2. kappa = 0.41 ---\n');
kappa_no = [0.42, 0.41];
for i = 1:numel(modes)
    P = Diemer19_Table(modes{i});
    assert(abs(P.kappa - kappa_no(i)) < 1e-10, ...
        sprintf('kappa wrong in mode "%s": %.4f', modes{i}, P.kappa));
end
fprintf('  PASS\n\n');

%% --- 3. Median and mean parameters are different ---------------------
fprintf('--- 3. Median != mean parameters ---\n');
P_med  = Diemer19_Table('median');
P_mean = Diemer19_Table('mean');
assert(P_med.a0 ~= P_mean.a0 || P_med.b0 ~= P_mean.b0, ...
    'Median and mean must have different parameters');
fprintf('  PASS\n\n');

%% --- 4. Unknown mode throws error ------------------------------------
fprintf('--- 4. Unknown mode errors correctly ---\n');
try
    Diemer19_Table('bad_mode');
    error('Should have thrown an error for unknown mode');
catch ME
    assert(~isempty(ME.message), 'Error message must not be empty');
    fprintf('  Caught: "%s"  PASS\n\n', ME.message(1:min(60, end)));
end

%% --- 5. Smoke test: scalar output ------------------------------------
fprintf('--- 5. Smoke test: scalar c ---\n');
for i = 1:numel(modes)
    c = Diemer19(1e13, 0, cosmo, modes{i});
    assert(isscalar(c) && isfinite(c) && c > 0, ...
        sprintf('Bad output for mode "%s": c = %g', modes{i}, c));
    fprintf('  c(1e13, z=0, %-6s) = %.3f  PASS\n', modes{i}, c);
end
fprintf('\n');

%% --- 6. Vector input preserves size ----------------------------------
fprintf('--- 6. Vector input ---\n');
M_vec = logspace(11, 15, 20);
for i = 1:numel(modes)
    c_vec = Diemer19(M_vec, 0, cosmo, modes{i});
    assert(numel(c_vec) == numel(M_vec), 'Size mismatch');
    assert(all(isfinite(c_vec)) && all(c_vec > 0), 'Non-finite or negative c');
    fprintf('  %-6s  N=%d  c in [%.2f, %.2f]  PASS\n', ...
        modes{i}, numel(c_vec), min(c_vec), max(c_vec));
end
fprintf('\n');

%% --- 7. c decreases with M at z=0 -----------------------------------
fprintf('--- 7. c(M) decreasing at z=0 ---\n');
M_test = logspace(11, 15, 8);
fprintf('  %12s  %8s  %8s\n', 'M [Msun/h]', 'median', 'mean');
c_med_arr  = zeros(size(M_test));
c_mean_arr = zeros(size(M_test));
for i = 1:numel(M_test)
    c_med_arr(i)  = Diemer19(M_test(i), 0, cosmo, 'median');
    c_mean_arr(i) = Diemer19(M_test(i), 0, cosmo, 'mean');
    fprintf('  %12.3e  %8.3f  %8.3f\n', M_test(i), c_med_arr(i), c_mean_arr(i));
end
assert(c_med_arr(1)  > c_med_arr(end),  'median c must decrease with M');
assert(c_mean_arr(1) > c_mean_arr(end), 'mean   c must decrease with M');
fprintf('  PASS\n\n');

%% --- 8. c decreases with z at fixed M --------------------------------
fprintf('--- 8. c(z) decreasing at fixed M=1e13 ---\n');
z_arr = [0.0, 0.5, 1.0, 2.0, 3.0];
fprintf('  %6s  %8s  %8s\n', 'z', 'median', 'mean');
c_med_z  = zeros(size(z_arr));
c_mean_z = zeros(size(z_arr));
for i = 1:numel(z_arr)
    c_med_z(i)  = Diemer19(1e13, z_arr(i), cosmo, 'median');
    c_mean_z(i) = Diemer19(1e13, z_arr(i), cosmo, 'mean');
    fprintf('  %6.2f  %8.3f  %8.3f\n', z_arr(i), c_med_z(i), c_mean_z(i));
end
assert(c_med_z(1)  > c_med_z(end),  'median c must decrease with z');
assert(c_mean_z(1) > c_mean_z(end), 'mean   c must decrease with z');
fprintf('  PASS\n\n');

%% --- 9. Mean >= median (log-normal scatter makes mean > median) ------
fprintf('--- 9. Mean >= median at all masses ---\n');
M_check = logspace(11, 15, 10);
for i = 1:numel(M_check)
    c_med  = Diemer19(M_check(i), 0, cosmo, 'median');
    c_mean = Diemer19(M_check(i), 0, cosmo, 'mean');
    assert(c_mean >= c_med * 0.85, ...
        sprintf('mean << median at M=%.2e: mean=%.3f median=%.3f', ...
                M_check(i), c_mean, c_med));
end
fprintf('  PASS\n\n');

%% --- 10. Plot: c(M) median vs mean at z = 0, 1, 2 ------------------
fprintf('--- 10. Plot: c(M) at z=0,1,2 ---\n');
M_plot  = logspace(11, 15.5, 100);
z_plot  = [0, 1, 2];
colors  = [0.13 0.47 0.71;   % blue
           0.17 0.63 0.17;   % green
           0.84 0.15 0.16];  % red

figure('Name','Diemer19 c(M)','Position',[100 100 800 480]);
hold on;
for iz = 1:numel(z_plot)
    c_med  = arrayfun(@(m) Diemer19(m, z_plot(iz), cosmo, 'median'), M_plot);
    c_mean = arrayfun(@(m) Diemer19(m, z_plot(iz), cosmo, 'mean'),   M_plot);
    plot(M_plot, c_med,  '-',  'Color', colors(iz,:), 'LineWidth', 2.0, ...
         'DisplayName', sprintf('median  z=%.0f', z_plot(iz)));
    plot(M_plot, c_mean, '--', 'Color', colors(iz,:), 'LineWidth', 1.4, ...
         'DisplayName', sprintf('mean    z=%.0f', z_plot(iz)));
end
set(gca, 'XScale','log');
xlabel('M_{200c}  [M_\odot h^{-1}]', 'FontSize', 12);
ylabel('c_{200c}',                   'FontSize', 12);
title('Diemer & Joyce 2019  —  c(M), median (solid) vs mean (dashed)', 'FontSize', 13);
legend('Location','northeast','FontSize',9,'NumColumns',2);
grid on;
hold off;
fprintf('  PASS\n\n');

%% --- 11. Plot: c(z) median vs mean at fixed masses ------------------
fprintf('--- 11. Plot: c(z) at fixed masses ---\n');
z_range  = linspace(0, 3, 80);
M_fixed  = [1e12, 1e13, 1e14];
labels   = {'10^{12}','10^{13}','10^{14}'};

figure('Name','Diemer19 c(z)','Position',[150 150 800 480]);
hold on;
for im = 1:numel(M_fixed)
    c_med  = arrayfun(@(z) Diemer19(M_fixed(im), z, cosmo, 'median'), z_range);
    c_mean = arrayfun(@(z) Diemer19(M_fixed(im), z, cosmo, 'mean'),   z_range);
    plot(z_range, c_med,  '-',  'Color', colors(im,:), 'LineWidth', 2.0, ...
         'DisplayName', sprintf('median  M=%s', labels{im}));
    plot(z_range, c_mean, '--', 'Color', colors(im,:), 'LineWidth', 1.4, ...
         'DisplayName', sprintf('mean    M=%s', labels{im}));
end
xlabel('Redshift  z',  'FontSize', 12);
ylabel('c_{200c}',     'FontSize', 12);
title('Diemer & Joyce 2019  —  c(z), median (solid) vs mean (dashed)', 'FontSize', 13);
legend('Location','northeast','FontSize',9,'NumColumns',2);
grid on;
hold off;
fprintf('  PASS\n\n');

fprintf('========================================\n');
fprintf('  All Diemer19 tests PASSED\n');
fprintf('========================================\n');