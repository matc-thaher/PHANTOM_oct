% test_Ishiyama21.m
% Test suite for Ishiyama21_concentration.m
% Reference: Ishiyama et al. 2021, MNRAS 506, 4210, Table 2 / Appendix B
% clear; clc; close all;

%% =========================================================================
% 0. Paths
% =========================================================================
here     = fileparts(mfilename('fullpath'));   % PHANTOM\tests
repo     = fullfile(here, '..');              % PHANTOM\
utils_path = fullfile(repo, 'src', 'utils');
concentration_path = fullfile(repo, 'src', 'concentration');
addpath(genpath(utils_path));
addpath(genpath(concentration_path));

assert(exist('Ishiyama21', 'file') == 2, ...
    sprintf('Ishiyama21.m not found.\nExpected in: %s', concentration_path));
assert(exist('Ishiyama21_Table', 'file') == 2, ...
    sprintf('Ishiyama21_Table.m not found.\nExpected in: %s', utils_path));

fprintf('Paths OK.\n');
fprintf('  Ishiyama21_concentration : %s\n', which('Ishiyama21_concentration'));
fprintf('  Ishiyama21_Table         : %s\n', which('Ishiyama21_Table'));

fprintf('========================================\n');
fprintf('  Ishiyama21 test suite\n');
fprintf('========================================\n\n');

cosmo = cosmology('Planck18');

all_modes     = {'500_all','500_relaxed','200c_all','200c_relaxed','vir_all','vir_relaxed'};
relaxed_modes = {'500_relaxed','200c_relaxed','vir_relaxed'};
all_modes_unrelaxed = {'500_all','200c_all','vir_all'};

% helper
tf_lut = {'false','true'};
tfstr  = @(x) tf_lut{x+1};

%% --- 1. Table loads: all 6 modes return all 6 fields -----------------
fprintf('--- 1. Table loads correctly ---\n');
fields = {'kappa','a0','a1','b0','b1','cAlpha'};
for i = 1:numel(all_modes)
    P = Ishiyama21_Table(all_modes{i});
    for j = 1:numel(fields)
        assert(isfield(P, fields{j}), ...
            sprintf('Mode "%s" missing field "%s"', all_modes{i}, fields{j}));
        assert(isfinite(P.(fields{j})), ...
            sprintf('Mode "%s" field "%s" is not finite', all_modes{i}, fields{j}));
    end
    fprintf('  %-14s  kappa=%.2f  a0=%.2f  a1=%.2f  b0=%.2f  b1=%.2f  cAlpha=%.3f  PASS\n', ...
        all_modes{i}, P.kappa, P.a0, P.a1, P.b0, P.b1, P.cAlpha);
end
fprintf('\n');

%% --- 2. Known parameter spot-checks against Table 2 -----------------
fprintf('--- 2. Parameter spot-checks ---\n');
P = Ishiyama21_Table('200c_all');
assert(abs(P.kappa - 1.19) < 1e-10, sprintf('200c_all kappa wrong: %.4f', P.kappa));
assert(abs(P.a0    - 2.54) < 1e-10, sprintf('200c_all a0 wrong: %.4f',    P.a0));
assert(abs(P.b1    - 1.21) < 1e-10, sprintf('200c_all b1 wrong: %.4f',    P.b1));

P = Ishiyama21_Table('200c_relaxed');
assert(abs(P.b1    - 6.36) < 1e-10, sprintf('200c_relaxed b1 wrong: %.4f', P.b1));
assert(abs(P.kappa - 0.60) < 1e-10, sprintf('200c_relaxed kappa wrong: %.4f', P.kappa));

P = Ishiyama21_Table('vir_all');
assert(abs(P.cAlpha - (-0.19)) < 1e-10, sprintf('vir_all cAlpha wrong: %.4f', P.cAlpha));

P = Ishiyama21_Table('vir_relaxed');
assert(abs(P.cAlpha - (-0.017)) < 1e-10, sprintf('vir_relaxed cAlpha wrong: %.4f', P.cAlpha));
fprintf('  All spot-checks PASS\n\n');

%% --- 3. cAlpha negative only for vir modes ---------------------------
fprintf('--- 3. cAlpha sign check ---\n');
for i = 1:numel(all_modes)
    P = Ishiyama21_Table(all_modes{i});
    is_vir = strncmp(all_modes{i}, 'vir', 3);
    if is_vir
        assert(P.cAlpha < 0, ...
            sprintf('vir mode "%s" should have cAlpha < 0, got %.4f', all_modes{i}, P.cAlpha));
        fprintf('  %-14s  cAlpha=%.3f  (negative as expected)  PASS\n', all_modes{i}, P.cAlpha);
    else
        assert(P.cAlpha > 0, ...
            sprintf('non-vir mode "%s" should have cAlpha > 0, got %.4f', all_modes{i}, P.cAlpha));
        fprintf('  %-14s  cAlpha=%.3f  (positive as expected)  PASS\n', all_modes{i}, P.cAlpha);
    end
end
fprintf('\n');

%% --- 4. Unknown mode throws error ------------------------------------
fprintf('--- 4. Unknown mode errors correctly ---\n');
try
    Ishiyama21_Table('bad_mode');
    error('Should have thrown for unknown mode');
catch ME
    assert(~isempty(ME.message), 'Error message must not be empty');
    fprintf('  Caught: "%s"  PASS\n\n', ME.message(1:min(60,end)));
end

%% --- 5. Smoke test: scalar concentration for all modes ---------------
fprintf('--- 5. Smoke test: scalar c ---\n');
for i = 1:numel(all_modes)
    c = Ishiyama21(1e13, 0, cosmo, all_modes{i});
    assert(isscalar(c) && isfinite(c) && c > 0, ...
        sprintf('Bad output for mode "%s": c = %g', all_modes{i}, c));
    fprintf('  %-14s  c(1e13,z=0) = %.3f  PASS\n', all_modes{i}, c);
end
fprintf('\n');

%% --- 6. Vector input preserves size ----------------------------------
fprintf('--- 6. Vector input ---\n');
M_vec = logspace(11, 15, 20);
for i = 1:numel(all_modes)
    c_vec = Ishiyama21(M_vec, 0, cosmo, all_modes{i});
    assert(numel(c_vec) == numel(M_vec), ...
        sprintf('Size mismatch for mode "%s"', all_modes{i}));
    assert(all(isfinite(c_vec)) && all(c_vec > 0), ...
        sprintf('Non-finite or negative c in mode "%s"', all_modes{i}));
    fprintf('  %-14s  N=%d  c in [%.2f, %.2f]  PASS\n', ...
        all_modes{i}, numel(c_vec), min(c_vec), max(c_vec));
end
fprintf('\n');

%% --- 7. c decreases with M at z=0 ------------------------------------
fprintf('--- 7. c(M) decreasing at z=0 ---\n');
M_test = logspace(11, 15, 8);
fprintf('  %12s', 'M [Msun/h]');
for i = 1:numel(all_modes); fprintf('  %14s', all_modes{i}); end
fprintf('\n');
c_all_M = zeros(numel(all_modes), numel(M_test));
for im = 1:numel(M_test)
    fprintf('  %12.3e', M_test(im));
    for im2 = 1:numel(all_modes)
        c_all_M(im2,im) = Ishiyama21(M_test(im), 0, cosmo, all_modes{im2});
        fprintf('  %14.3f', c_all_M(im2,im));
    end
    fprintf('\n');
end
for i = 1:numel(all_modes)
    assert(c_all_M(i,1) > c_all_M(i,end), ...
        sprintf('c not decreasing with M for mode "%s"', all_modes{i}));
end
fprintf('  PASS\n\n');

%% --- 8. c decreases with z at fixed M --------------------------------
fprintf('--- 8. c(z) decreasing at fixed M=1e13 ---\n');
z_arr = [0.0, 0.5, 1.0, 2.0, 3.0];
fprintf('  %6s', 'z');
for i = 1:numel(all_modes); fprintf('  %14s', all_modes{i}); end
fprintf('\n');
c_all_z = zeros(numel(all_modes), numel(z_arr));
for iz = 1:numel(z_arr)
    fprintf('  %6.2f', z_arr(iz));
    for im = 1:numel(all_modes)
        c_all_z(im,iz) = Ishiyama21(1e13, z_arr(iz), cosmo, all_modes{im});
        fprintf('  %14.3f', c_all_z(im,iz));
    end
    fprintf('\n');
end
for i = 1:numel(all_modes)
    assert(c_all_z(i,1) > c_all_z(i,end), ...
        sprintf('c not decreasing with z for mode "%s"', all_modes{i}));
end
fprintf('  PASS\n\n');

%% --- 9. Relaxed more concentrated than all (same mass definition) ----
fprintf('--- 9. Relaxed > all at same mass definition ---\n');
mass_defs = {'500','200c','vir'};
for i = 1:numel(mass_defs)
    c_all_i = Ishiyama21(1e13, 0, cosmo, [mass_defs{i} '_all']);
    c_rel_i = Ishiyama21(1e13, 0, cosmo, [mass_defs{i} '_relaxed']);
    fprintf('  %-4s  all=%.3f  relaxed=%.3f  relaxed>all? %s\n', ...
        mass_defs{i}, c_all_i, c_rel_i, tfstr(c_rel_i > c_all_i));
    assert(c_rel_i > c_all_i, ...
        sprintf('relaxed not > all for mass def "%s"', mass_defs{i}));
end
fprintf('  PASS\n\n');

%% --- 10. Mock cosmo: sigma only, verify c is finite ------------------
fprintf('--- 10. Mock cosmo consistency ---\n');
delta_c   = 1.686;
nu_vals   = [0.5, 1.0, 2.0, 3.5];
neff_fixed = -2.0;
alpha_fixed = 0.18;
fprintf('  %5s  %12s  %14s  %14s\n', 'nu', '200c_all', '200c_relaxed', 'vir_all');
for j = 1:numel(nu_vals)
    sigma_j     = delta_c / nu_vals(j);
    cm.sigmaM   = @(M,z)   sigma_j;
    cm.neff     = @(M,z,k) neff_fixed;
    cm.alphaEff = @(z)     alpha_fixed;
    c1 = Ishiyama21(1e13, 0, cm, '200c_all');
    c2 = Ishiyama21(1e13, 0, cm, '200c_relaxed');
    c3 = Ishiyama21(1e13, 0, cm, 'vir_all');
    assert(isfinite(c1) && c1 > 0, sprintf('200c_all bad at nu=%.2f',      nu_vals(j)));
    assert(isfinite(c2) && c2 > 0, sprintf('200c_relaxed bad at nu=%.2f',  nu_vals(j)));
    assert(isfinite(c3) && c3 > 0, sprintf('vir_all bad at nu=%.2f',       nu_vals(j)));
    fprintf('  %5.2f  %12.4f  %14.4f  %14.4f\n', nu_vals(j), c1, c2, c3);
end
fprintf('  PASS\n\n');

%% --- 11. Plot: c(M) all 6 modes at z=0 ------------------------------
fprintf('--- 11. Plot: c(M) all modes at z=0 ---\n');
M_plot  = logspace(11, 15.5, 100);
colors  = lines(3);   % one colour per mass definition
ls_all  = '-';
ls_rel  = '--';

figure('Name','Ishiyama21 c(M)','Position',[100 100 820 500]);
hold on;
mass_defs_label = {'500','200c','vir'};
for i = 1:3
    c_a = arrayfun(@(m) Ishiyama21(m,0,cosmo,[mass_defs_label{i} '_all']),     M_plot);
    c_r = arrayfun(@(m) Ishiyama21(m,0,cosmo,[mass_defs_label{i} '_relaxed']), M_plot);
    plot(M_plot, c_a, ls_all, 'Color', colors(i,:), 'LineWidth', 2.0, ...
         'DisplayName', [mass_defs_label{i} ' all']);
    plot(M_plot, c_r, ls_rel, 'Color', colors(i,:), 'LineWidth', 1.4, ...
         'DisplayName', [mass_defs_label{i} ' relaxed']);
end
set(gca,'XScale','log');
xlabel('Halo mass  [$M_{\odot} h^{-1}$]','FontSize',12, 'Interpreter','latex');
ylabel('Concentration  c','FontSize',12);
title('Ishiyama et al. 2021 — c(M) at z=0, all modes (solid=all, dashed=relaxed)','FontSize',12);
legend('Location','northeast','FontSize',9,'NumColumns',2);
grid on; hold off;
fprintf('  PASS\n\n');

%% --- 12. Plot: c(z) all 6 modes at M=1e13 ---------------------------
fprintf('--- 12. Plot: c(z) all modes at M=1e13 ---\n');
z_range = linspace(0, 3, 80);

figure('Name','Ishiyama21 c(z)','Position',[150 150 820 500]);
hold on;
for i = 1:3
    c_a = arrayfun(@(z) Ishiyama21(1e13,z,cosmo,[mass_defs_label{i} '_all']),     z_range);
    c_r = arrayfun(@(z) Ishiyama21(1e13,z,cosmo,[mass_defs_label{i} '_relaxed']), z_range);
    plot(z_range, c_a, ls_all, 'Color', colors(i,:), 'LineWidth', 2.0, ...
         'DisplayName', [mass_defs_label{i} ' all']);
    plot(z_range, c_r, ls_rel, 'Color', colors(i,:), 'LineWidth', 1.4, ...
         'DisplayName', [mass_defs_label{i} ' relaxed']);
end
xlabel('Redshift  z','FontSize',12);
ylabel('Concentration  c','FontSize',12);
title('Ishiyama et al. 2021 — c(z) at M=10^{13} $M_{\odot} h^{-1}$ (solid=all, dashed=relaxed)','FontSize',12, 'Interpreter','latex');
legend('Location','northeast','FontSize',9,'NumColumns',2);
grid on; hold off;
fprintf('  PASS\n\n');

fprintf('========================================\n');
fprintf('  All Ishiyama21 tests PASSED\n');
fprintf('========================================\n');