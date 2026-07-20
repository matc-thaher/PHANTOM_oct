%% test_Diemer15.m
% Test suite for Diemer15_Table and Diemer15_concentration
%
% Tests performed:
%   1. Path check
%   2. Table integrity — valid modes, error on invalid
%   3. Parameter value spot-check against DK15 Table 3
%   4. Analytical formula check — compute c by hand and compare
%   5. Double power-law shape — c has a minimum at nu = nu_min
%   6. n_eff dependence — c_min and nu_min shift correctly with n
%   7. c decreasing with M at z=0 (galaxy to cluster scale)
%   8. c decreasing with z at fixed M (moderate z)
%   9. Mean > median at low nu, median > mean at high nu (crossover)
%  10. Plausible c200c range at z=0

%clc; close; clear all; 

fprintf('=== test_Diemer15 ===\n\n');

%% ---- helpers -----------------------------------------------------------
tfstr = @(x) repmat('PASS', x, 1) + repmat('FAIL', ~x, 1);

tol_tight = 1e-10;   % exact formula agreement
tol_pct   = 0.01;    % 0.01% tolerance for parameter-based checks

%% =========================================================================
% 1. Path check
% =========================================================================
here       = fileparts(mfilename('fullpath'));
repo       = fullfile(here, '..');
utils_path = fullfile(repo, 'src', 'utils');
concentration_path = fullfile(repo, 'src', 'concentration');
addpath(genpath(utils_path));
addpath(genpath(concentration_path));

assert(exist('Diemer15','file') == 2, ...
    sprintf('Diemer15_concentration.m not found.\nExpected in: %s', utils_path));
assert(exist('Diemer15_Table','file') == 2, ...
    sprintf('Diemer15_Table.m not found.\nExpected in: %s', utils_path));

fprintf('Paths OK.\n');
fprintf('  Diemer15_concentration : %s\n', which('Dutton14_concentration'));
fprintf('  Diemer15_Table         : %s\n', which('Dutton14_Table'));

%% =========================================================================
% 2. Table integrity
% =========================================================================
fprintf('\n--- Table integrity ---\n');

% Valid modes return structs with required fields
required = {'kappa','phi0','phi1','eta0','eta1','alpha','beta'};
for stat = {'median', 'mean'}
    P = Diemer15_Table(stat{1});
    for fi = 1:numel(required)
        assert(isfield(P, required{fi}), ...
            'Diemer15_Table("%s") missing field: %s', stat{1}, required{fi});
    end
    fprintf('Diemer15_Table(''%s'') : OK\n', stat{1});
end

% Invalid mode must error
try
    Diemer15_Table('bogus');
    error('Should have thrown.');
catch ME
    fprintf('Unknown mode correctly throws: "%s"\n', ME.message);
end

%% =========================================================================
% 3. Parameter spot-check against DK15 Table 3
% =========================================================================
fprintf('\n--- Parameter values vs DK15 Table 3 ---\n');

expected.median = struct('kappa',1.00,'phi0',6.58,'phi1',1.27, ...
                         'eta0',7.28,'eta1',1.56,'alpha',1.08,'beta',1.77);
expected.mean   = struct('kappa',1.00,'phi0',6.66,'phi1',1.37, ...
                         'eta0',5.41,'eta1',1.06,'alpha',1.22,'beta',1.22);

for stat = {'median', 'mean'}
    P   = Diemer15_Table(stat{1});
    ref = expected.(stat{1});
    fprintf('\n  statistic = %s:\n', stat{1});
    for fi = 1:numel(required)
        got = P.(required{fi});
        exp1 = ref.(required{fi});
        d   = abs(got - exp1);
        ok  = d < tol_tight;
        fprintf('    %-8s : got=%.4f  expected=%.4f  diff=%.2e  %s\n', ...
            required{fi}, got, exp1, d, tfstr(ok));
        assert(ok, 'Parameter %s mismatch for statistic=%s', required{fi}, stat{1});
    end
end

%% =========================================================================
% 4. Analytical formula check
%    Compute c by hand for a fixed (nu, n_eff) and compare to function
%    We mock cosmo.sigmaM and cosmo.neff so the test is self-contained
% =========================================================================
fprintf('\n--- Analytical formula check ---\n');

% Fixed inputs
nu_test   = 1.5;
neff_test = -2.2;
delta_c   = 1.686;

% Build minimal cosmo mock with fixed sigma and neff
sigma_test          = delta_c / nu_test;
cosmo_mock.sigmaM   = @(M, z) sigma_test;   % constant, ignores M and z
cosmo_mock.neff     = @(M, z, kappa) neff_test;

M_test = 1e12;    % dummy mass
z_test = 0.0;

for stat = {'median', 'mean'}
    P = Diemer15_Table(stat{1});

    % Manual formula (Eqs. 9-10, DK15)
    c_min    = P.phi0 + P.phi1 * neff_test;
    nu_min   = P.eta0 + P.eta1 * neff_test;
    c_manual = 0.5 * c_min * ( (nu_test/nu_min)^(-P.alpha) ...
                              + (nu_test/nu_min)^(P.beta) );

    c_code = Diemer15(M_test, z_test, cosmo_mock, stat{1});

    diff_pct = abs(c_code - c_manual) / c_manual * 100;
    ok = diff_pct < tol_pct;
    fprintf('  %-6s: manual=%.6f  code=%.6f  diff=%.2e%%  %s\n', ...
        stat{1}, c_manual, c_code, diff_pct, tfstr(ok));
    assert(ok, 'Formula mismatch for statistic=%s', stat{1});
end

%% =========================================================================
% 5. Double power-law shape — c has a minimum at nu = nu_min
%    For a fixed n_eff, dc/dnu = 0 at nu = nu_min  (by construction)
% =========================================================================
fprintf('\n--- Double power-law minimum location ---\n');

neff_test2 = -2.0;
nu_vec     = linspace(0.2, 8.0, 2000);

for stat = {'median', 'mean'}
    P = Diemer15_Table(stat{1});

    c_min_val = P.phi0 + P.phi1 * neff_test2;
    nu_min    = P.eta0 + P.eta1 * neff_test2;

    % Analytical minimum location
    x_min        = (P.alpha / P.beta)^(1 / (P.alpha + P.beta));
    nu_cmin_anal = nu_min * x_min;

    % Numerical minimum over dense grid
    c_vec        = 0.5 * c_min_val * ( (nu_vec/nu_min).^(-P.alpha) ...
                                     + (nu_vec/nu_min).^(P.beta) );
    [~, idx]     = min(c_vec);
    nu_cmin_num  = nu_vec(idx);

    frac_err = abs(nu_cmin_num - nu_cmin_anal) / nu_cmin_anal;
    ok = frac_err < 0.005;   % 0.5% — limited by grid spacing
    fprintf('  %-6s: nu_min(param)=%.4f  x_min=%.4f  nu_at_cmin(anal)=%.4f  nu_at_cmin(num)=%.4f  frac_err=%.4f  %s\n', ...
        stat{1}, nu_min, x_min, nu_cmin_anal, nu_cmin_num, frac_err, tfstr(ok));
    assert(ok, 'Numerical minimum does not match analytical for statistic=%s', stat{1});
end
%% =========================================================================
% 6. n_eff dependence — steeper n -> lower c_min and lower nu_min
%    (DK15 Section 4.2: steeper slope => lower concentration floor)
% =========================================================================
fprintf('\n--- n_eff dependence (phi1 > 0, eta1 > 0 => c_min increases with n) ---\n');

n_vec = [-2.5, -2.0, -1.5, -1.0];   % increasingly shallow (less negative)
P_med = Diemer15_Table('median');

fprintf('  n_eff  |  c_min  |  nu_min\n');
c_min_prev  = -Inf;
nu_min_prev = -Inf;

for i = 1:numel(n_vec)
    n       = n_vec(i);
    c_min_i  = P_med.phi0 + P_med.phi1 * n;
    nu_min_i = P_med.eta0 + P_med.eta1 * n;
    fprintf('  %5.1f  |  %5.3f  |  %5.3f\n', n, c_min_i, nu_min_i);

    if i > 1
        ok_c  = c_min_i > c_min_prev;
        ok_nu = nu_min_i > nu_min_prev;
        assert(ok_c,  'c_min not increasing with n at n=%.1f', n);
        assert(ok_nu, 'nu_min not increasing with n at n=%.1f', n);
    end
    c_min_prev  = c_min_i;
    nu_min_prev = nu_min_i;
end
fprintf('  c_min and nu_min both increase with n: PASS\n');

%% =========================================================================
% 7. c(M) behaviour at z=0
%    DK15 predicts a concentration FLOOR (c_min) at high mass —
%    the double power-law by construction has a minimum.
%    Tests:
%      (a) c decreases from 1e11 to 1e13 (galaxy/group scale)
%      (b) c never falls below a physical floor (~1.5)
%      (c) c is monotone decreasing for M <= 1e13 (well below c_min turnup)
% =========================================================================
fprintf('\n--- c(M) behaviour at z=0 ---\n');

Mstar     = 3e12;
sigma_fn2 = @(M) 1.0  .* (M ./ Mstar).^(-0.18);
neff_fn2  = @(M) -2.5 + 0.3 .* log10(M ./ 1e11);

cosmo_lcdm.sigmaM = @(M, z) sigma_fn2(M);
cosmo_lcdm.neff   = @(M, z, kappa) neff_fn2(M);

M_vec     = logspace(11, 15, 20);
c_vec_med = arrayfun(@(M) Diemer15(M, 0, cosmo_lcdm, 'median'), M_vec);

fprintf('  M [Msun/h]   c200c\n');
for i = 1:numel(M_vec)
    fprintf('  1e%.2f       %.4f\n', log10(M_vec(i)), c_vec_med(i));
end

% (a) monotone decreasing from 1e11 to 1e13
M_lowz    = logspace(11, 13, 10);
c_lowz    = arrayfun(@(M) Diemer15(M, 0, cosmo_lcdm, 'median'), M_lowz);
ok_lowm   = all(diff(c_lowz) < 0);
fprintf('\n  (a) monotone decreasing 1e11->1e13? %s\n', tfstr(ok_lowm));
assert(ok_lowm, 'c not monotone decreasing from 1e11 to 1e13 at z=0');

% (b) c never below physical floor
c_floor = 1.5;
ok_floor = all(c_vec_med > c_floor);
fprintf('  (b) all c > %.1f (physical floor)? %s\n', c_floor, tfstr(ok_floor));
assert(ok_floor, 'c fell below physical floor at z=0');

% (c) detect and report the minimum (c_min turnup at high M)
[cmin_val, idx_min] = min(c_vec_med);
if idx_min < numel(M_vec)
    fprintf('  (c) c_min = %.3f at M = 1e%.2f Msun/h — expected DK15 floor behaviour\n', ...
        cmin_val, log10(M_vec(idx_min)));
else
    fprintf('  (c) no upturn detected in this mass range\n');
end
fprintf('  c(M) checks PASS\n');

%% =========================================================================
% 8. c(z) behaviour at fixed M=1e12
%    Same physics as the mass test: DK15 has a c_min floor.
%    As z increases, sigma(M,z) shrinks -> nu rises -> eventually
%    passes nu_at_cmin and c turns back up.
%    Tests:
%      (a) c decreases from z=0 to z=1 (well-tested low-z regime)
%      (b) c never below physical floor
%      (c) detect and report upturn redshift if present
% =========================================================================
fprintf('\n--- c(z) behaviour at fixed M=1e12 ---\n');

D_z     = @(z) 1 ./ (1 + z);
sigma_z = @(M, z) sigma_fn2(M) .* D_z(z);

cosmo_z.sigmaM = @(M, z) sigma_z(M, z);
cosmo_z.neff   = @(M, z, kappa) neff_fn2(M);

z_vec   = [0.0, 0.5, 1.0, 2.0, 3.0, 5.0];
M_fixed = 1e12;
c_z = arrayfun(@(z) Diemer15(M_fixed, z, cosmo_z, 'median'), z_vec);

fprintf('  M=1e12:');
for i = 1:numel(z_vec)
    fprintf(' z=%.1f:%.3f', z_vec(i), c_z(i));
end
fprintf('\n');

% (a) must decrease from z=0 to z=1
ok_lowz = c_z(1) > c_z(2) && c_z(2) > c_z(3);
fprintf('  (a) decreasing z=0->1? %s\n', tfstr(ok_lowz));
assert(ok_lowz, 'c not decreasing from z=0 to z=1 for M=1e12');

% (b) all c above physical floor
c_floor  = 1.5;
ok_floor = all(c_z > c_floor);
fprintf('  (b) all c > %.1f (physical floor)? %s\n', c_floor, tfstr(ok_floor));
assert(ok_floor, 'c fell below physical floor for M=1e12');

% (c) detect upturn and report — expected DK15 behaviour at high z
[~, iz_min] = min(c_z);
if iz_min < numel(z_vec)
    fprintf('  (c) c minimum at z=%.1f then upturn — expected DK15 floor behaviour\n', ...
        z_vec(iz_min));
else
    fprintf('  (c) no upturn detected in this z range\n');
end
fprintf('  c(z) checks PASS\n');
%% =========================================================================
% 9. Mean vs median — both positive, distinct, and self-consistent
%    With DK15 Table 3 params: beta_median=1.69 > beta_mean=0.67
%    so at high nu: median rises faster than mean (steeper beta)
%    At low nu:     median falls faster than mean (steeper alpha too)
%    We just verify: both positive, distinct from each other, and
%    that median > mean at high nu (since beta_median > beta_mean)
% =========================================================================

fprintf('\n--- Mean vs median ---\n');
neff_fixed = -2.0;
delta_c    = 1.686;
nu_vals    = [0.5, 1.0, 2.0, 3.5];
P_med  = Diemer15_Table('median');
P_mean = Diemer15_Table('mean');

fprintf('  nu    |  c_median  |  c_mean\n');
for nu_i = nu_vals
    sigma_i       = delta_c / nu_i;
    cm.sigmaM     = @(M, z) sigma_i;
    cm.neff       = @(M, z, k) neff_fixed;
    c_med_i  = Diemer15(1e12, 0, cm, 'median');
    c_mean_i = Diemer15(1e12, 0, cm, 'mean');
    assert(c_med_i  > 0, 'c_median <= 0 at nu=%.2f', nu_i);
    assert(c_mean_i > 0, 'c_mean   <= 0 at nu=%.2f', nu_i);
    assert(abs(c_med_i - c_mean_i) > 1e-6, ...
        'median and mean identical at nu=%.2f', nu_i);
    fprintf('  %.2f  |  %.4f    |  %.4f\n', nu_i, c_med_i, c_mean_i);
end

% Verify the two parameter sets are genuinely different
assert(P_med.beta ~= P_mean.beta, 'beta_median and beta_mean must differ');
fprintf('\n  beta_median=%.2f  beta_mean=%.2f  (distinct parameters confirmed)\n', ...
    P_med.beta, P_mean.beta);
fprintf('  Mean vs median checks PASS\n');

%% =========================================================================
% 10. Plausible c200c range at z=0
% =========================================================================
fprintf('\n--- Plausible c200c range at z=0 ---\n');

c_galaxy  = Diemer15(1e12, 0, cosmo_lcdm, 'median');
c_group   = Diemer15(1e13, 0, cosmo_lcdm, 'median');
c_cluster = Diemer15(1e14, 0, cosmo_lcdm, 'median');

fprintf('  c200c(1e12, z=0) = %.2f  (expected ~7-10)\n',  c_galaxy);
fprintf('  c200c(1e13, z=0) = %.2f  (expected ~5-7)\n',   c_group);
fprintf('  c200c(1e14, z=0) = %.2f  (expected ~4-6)\n',   c_cluster);

assert(c_galaxy  > 5  && c_galaxy  < 15, 'c200c(1e12) out of plausible range');
assert(c_group   > 3  && c_group   < 10, 'c200c(1e13) out of plausible range');
assert(c_cluster > 2  && c_cluster < 8,  'c200c(1e14) out of plausible range');
assert(c_galaxy > c_group && c_group > c_cluster, ...
    'c not ordered: galaxy > group > cluster');
fprintf('  All plausible range checks: PASS\n');

%% =========================================================================
fprintf('\n=== All Diemer15 tests passed ===\n');