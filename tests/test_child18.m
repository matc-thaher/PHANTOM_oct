% test_Child18.m
% Test suite for Child18_concentration.m
% Reference: Child et al. 2018, ApJ 859 55

%clear; clc; close all;

fprintf('========================================\n');
fprintf('  Child18 concentration test suite\n');
fprintf('========================================\n\n');



%% -- 0. Path check

here     = fileparts(mfilename('fullpath'));   % PHANTOM\tests
repo     = fullfile(here, '..');              % PHANTOM\
utils_path = fullfile(repo, 'src', 'utils');
concentration_path = fullfile(repo, 'src', 'concentration');
addpath(genpath(utils_path));
addpath(genpath(concentration_path));

assert(exist('Child18', 'file') == 2, ...
    sprintf('Child18.m not found.\nExpected in: %s', utils_path));
assert(exist('Child18_Table', 'file') == 2, ...
    sprintf('Child18_Table.m not found.\nExpected in: %s', utils_path));

fprintf('Paths OK.\n');
fprintf('  Child18_concentration : %s\n', which('Child18_concentration'));
fprintf('  Child18_Table         : %s\n', which('Child18_Table'));

cosmo = cosmology('Planck18');

%% --- 1. Smoke test: scalar input ---------------------------------------
fprintf('--- 1. Smoke test ---\n');
[c, Mstar] = Child18(1e13, 0, cosmo);
assert(isscalar(c),    'c must be scalar for scalar M input');
assert(isfinite(c),    'c must be finite');
assert(c > 0,          'c must be positive');
assert(isfinite(Mstar) && Mstar > 0, 'Mstar must be positive finite');
fprintf('  c(1e13, z=0) = %.3f,  Mstar = %.3e Msun/h  PASS\n\n', c, Mstar);

%% --- 2. Vector input preserves shape ----------------------------------
fprintf('--- 2. Vector input ---\n');
M_vec = logspace(11, 15, 20);
c_vec = Child18(M_vec, 0, cosmo);
assert(numel(c_vec) == numel(M_vec), 'output size must match input');
assert(all(isfinite(c_vec)),         'all c must be finite');
assert(all(c_vec > 0),               'all c must be positive');
fprintf('  N=%d, c in [%.2f, %.2f]  PASS\n\n', numel(c_vec), min(c_vec), max(c_vec));

%% --- 3. c decreases with M at z=0 -------------------------------------
fprintf('--- 3. c(M) decreasing at z=0 ---\n');
M_test = logspace(11, 15, 8);
c_test = Child18(M_test, 0, cosmo);
fprintf('  %12s  %6s\n', 'M [Msun/h]', 'c');
for i = 1:numel(M_test)
    fprintf('  %12.3e  %6.3f\n', M_test(i), c_test(i));
end
assert(c_test(1) > c_test(end), 'c must decrease from low to high M');
fprintf('  PASS\n\n');

%% --- 4. c decreases with z at fixed M ---------------------------------
fprintf('--- 4. c(z) decreasing at fixed M=1e13 ---\n');
z_arr = [0.0, 0.5, 1.0, 2.0, 3.0];
c_z   = zeros(size(z_arr));
for i = 1:numel(z_arr)
    c_z(i) = Child18(1e13, z_arr(i), cosmo);
end
fprintf('  z:  '); fprintf('%6.2f ', z_arr); fprintf('\n');
fprintf('  c:  '); fprintf('%6.3f ', c_z);   fprintf('\n');
assert(c_z(1) > c_z(end), 'c must decrease with z at fixed M');
fprintf('  PASS\n\n');

%% --- 5. High-mass plateau: c -> c0 ------------------------------------
fprintf('--- 5. High-mass concentration floor -> c0 ---\n');
% At M >> MT, c -> c0. Table 1 individual_all: c0 = 3.19
c0_expected = 3.19;
c_massive   = Child18(1e16, 3.0, cosmo, 'individual_all');
fprintf('  c(1e16, z=3) = %.3f  (expected near c0=%.2f)\n', c_massive, c0_expected);
assert(abs(c_massive - c0_expected) < 0.5, ...
    sprintf('c floor expected ~%.2f, got %.3f', c0_expected, c_massive));
fprintf('  PASS\n\n');

%% --- 6. Mstar decreases with z ----------------------------------------
fprintf('--- 6. Mstar decreasing with z ---\n');
z_arr2  = [0, 1, 2, 3];
Ms_arr  = zeros(size(z_arr2));
for i = 1:numel(z_arr2)
    [~, Ms_arr(i)] = Child18(1e13, z_arr2(i), cosmo);
end
fprintf('  z:            '); fprintf('%8.1f ', z_arr2);       fprintf('\n');
fprintf('  log10(Mstar): '); fprintf('%8.2f ', log10(Ms_arr)); fprintf('\n');
assert(Ms_arr(1) > Ms_arr(end), 'Mstar must decrease with z');
fprintf('  PASS\n\n');

%% --- 7. Relaxed sample more concentrated than all ---------------------
fprintf('--- 7. Relaxed > all halos ---\n');
c_all = Child18(1e13, 0, cosmo, 'individual_all');
c_rel = Child18(1e13, 0, cosmo, 'individual_relaxed');
fprintf('  c (all)     = %.3f\n', c_all);
fprintf('  c (relaxed) = %.3f\n', c_rel);
assert(c_rel > c_all, 'relaxed halos must be more concentrated');
fprintf('  PASS\n\n');

%% --- 8. Analytic check at M = M* (x=1) --------------------------------
fprintf('--- 8. Analytic value at M=Mstar ---\n');
% At M = Mstar: x=1, xT = 1/m
% c = A * 1^b / (1 + m^b) + c0
[~, Mstar0] = Child18(1e13, 0, cosmo, 'individual_all');
c_at_mstar  = Child18(Mstar0, 0, cosmo, 'individual_all');
m = -0.10; A = 3.44; b = 430.49; c0 = 3.19;
c_expected  = A * ((1/b)^m * (1 + (1/b))^(-m)-1) + c0;
fprintf('  c(Mstar, z=0) = %.4f  (analytic = %.4f)\n', c_at_mstar, c_expected);
assert(abs(c_at_mstar - c_expected) < 0.01, ...
    'c at M=Mstar does not match analytic Table 1 value');
fprintf('  PASS\n\n');

%% --- 9. All fit types run without error --------------------------------
fprintf('--- 9. All fit types ---\n');
types = {'individual_all','individual_relaxed','stack_nfw','stack_einasto'};
for i = 1:numel(types)
    cv = Child18(1e13, 0.5, cosmo, types{i});
    assert(isfinite(cv) && cv > 0, ...
        sprintf('fit_type "%s" returned bad value %.3f', types{i}, cv));
    fprintf('  %-22s  c = %.3f  PASS\n', types{i}, cv);
end
fprintf('\n');

fprintf('========================================\n');
fprintf('  All Child18 tests PASSED\n');
fprintf('========================================\n');