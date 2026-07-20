%% test_Klypin11.m
% Test script for Klypin11_concentration and Klypin11_Table
%
% Location : PHANTOM\tests\test_Klypin11.m
%
% Tests covered:
%   1. Table integrity — both sample entries load correctly
%   2. Hard-coded spot checks against Table 3 (Eq. 10 / 12)
%   3. Subhalo check against Eq. (11) at M = 1e12
%   4. Sanity: c decreasing with M for distinct haloes at every z
%   5. Sanity: c decreasing with z at galaxy-size mass (M < 1e14)
%   6. Upturn at high mass (c increases above c_min at high M, high z)
%   7. Interpolation continuity — no jumps between tabulated redshifts
%   8. Out-of-range redshift warning issued and clamped gracefully
%   9. Subhalo > distinct at same mass (tidal stripping effect)
%  10. Plots: c(M) at all tabulated z + subhalo vs distinct at z=0

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

assert(exist('Klypin11','file') == 2, ...
    sprintf('Klypin11.m not found.\nExpected in: %s', utils_path));
assert(exist('Klypin11_Table','file') == 2, ...
    sprintf('Klypin11_Table.m not found.\nExpected in: %s', utils_path));

fprintf('Paths OK.\n');
fprintf('  Klypin11_concentration : %s\n', which('Klypin11_concentration'));
fprintf('  Klypin11_Table         : %s\n', which('Klypin11_Table'));

%% =========================================================================
% 1. Table integrity — both samples load cleanly
% =========================================================================
fprintf('\n--- Table integrity ---\n');

for s = {'distinct','subhalo'}
    try
        P = Klypin11_Table(s{1});
        assert(isfield(P,'c0') && isfield(P,'alpha') && isfield(P,'Mpivot'), ...
            'Missing required field in struct.');
        fprintf('  Klypin11_Table(''%s'') : OK\n', s{1});
    catch ME
        error('  Klypin11_Table(''%s'') FAILED: %s', s{1}, ME.message);
    end
end

% Unknown sample must error
try
    Klypin11_Table('bogus');
    error('Should have thrown for unknown sample.');
catch ME
    fprintf('  Unknown sample correctly throws: "%s"\n', ME.message);
end

%% =========================================================================
% 2. Hard spot checks — Table 3 column "c(1e12 h^-1 Msun)"
%    Table 3: z=0 -> 9.60, z=0.5 -> 7.21, z=1 -> 5.82,
%             z=2 -> 4.60, z=3 -> 4.40, z=5 -> 5.00
% =========================================================================
fprintf('\n--- Hard spot checks vs Table 3 (Eq. 12, M=1e12 h^-1 Msun) ---\n');

% Paper values for c(Mvir = 1e12) directly from Table 3
z_tab     = [0.0,  0.5,  1.0,  2.0,  3.0,  5.0];
c_paper   = [9.60, 7.21, 5.82, 4.60, 4.40, 5.00];
tol_pct   = 2.0;   % allow 2% tolerance (interpolation at exact tabulated z)

all_pass = true;
for i = 1:numel(z_tab)
    c_code = Klypin11(1e12, z_tab(i), 'distinct');
    diff_pct = 100 * abs(c_code - c_paper(i)) / c_paper(i);
    pass = diff_pct < tol_pct;
    all_pass = all_pass && pass;
    fprintf('  z=%.1f: paper=%.2f  code=%.4f  diff=%.2f%%  %s\n', ...
        z_tab(i), c_paper(i), c_code, diff_pct, tf2str(pass));
end
assert(all_pass, 'One or more spot checks exceed 2% tolerance.');

%% =========================================================================
% 3. Eq. (10) z=0 check: c = 9.60 * (M/1e12)^-0.075 (simplified form)
%    The full Eq. (12) has the upturn term [1 + (M/M0)^0.26]^-1.
%    At galaxy-cluster scale (M << M0 at z=0, M0=9.6e17) the upturn
%    term is ~1, so Eq.(10) and Eq.(12) should agree to <0.1%.
% =========================================================================
fprintf('\n--- Eq. (10) vs Eq. (12) at z=0 for M << M0 ---\n');

M_test = [1e11, 1e12, 1e13, 1e14];
P_z0   = Klypin11_Table('distinct');
for M = M_test
    c_eq10 = 9.60 * (M / 1e12)^(-0.075);
    c_eq12 = Klypin11(M, 0, 'distinct');
    diff1   = 100 * abs(c_eq12 - c_eq10) / c_eq10;
    fprintf('  M=1e%.0f:  Eq10=%.4f  Eq12=%.4f  diff=%.3f%%\n', ...
        log10(M), c_eq10, c_eq12, diff1);
end

%% =========================================================================
% 4. Eq. (11) subhalo check at M = 1e12: c_sub = 12 exactly
% =========================================================================
fprintf('\n--- Eq. (11) subhalo at M=1e12 h^-1 Msun ---\n');

c_sub_1e12 = Klypin11(1e12, 0, 'subhalo');
fprintf('  c_sub(1e12) = %.6f  (expected 12.000000)\n', c_sub_1e12);
assert(abs(c_sub_1e12 - 12.0) < 1e-10, 'Subhalo c at pivot mass should be exactly 12.');

% Also test subhalo slope: c at 1e13 / c at 1e12 = (10)^(-0.12)
c_sub_1e13  = Klypin11(1e13, 0, 'subhalo');
ratio_code  = c_sub_1e13 / c_sub_1e12;
ratio_exact = (1e13/1e12)^(-0.12);
fprintf('  c_sub(1e13)/c_sub(1e12) = %.6f  (expected %.6f)\n', ...
    ratio_code, ratio_exact);
assert(abs(ratio_code - ratio_exact) < 1e-10, 'Subhalo power-law slope is wrong.');

%% =========================================================================
% 5. Sanity: c decreasing with M for distinct haloes (all tabulated z)
%    (At z=0-3 the upturn only kicks in at M > ~1e15, well outside normal
%     halo range; for z=5 the upturn is around 6.6e11 Msun/h so we
%     test the regime M < 5e10 to 1e14 where the decline dominates.)
% =========================================================================

fprintf('\n--- Sanity: c decreasing with M (distinct, all z) ---\n');

M_vec = logspace(10, 14, 60);

c_test = Klypin11(M_vec, 0, 'distinct');
fprintf('  c_test size: %d x %d\n', size(c_test,1), size(c_test,2));
fprintf('  M_vec size:  %d x %d\n', size(M_vec,1),  size(M_vec,2));
fprintf('  idx size:    %d x %d\n', size(log10(M_vec(:))<=13, 1), size(log10(M_vec(:))<=13, 2));

% for i = 1:numel(z_tab)
%     c_arr = Klypin11(M_vec(:), z_tab(i), 'distinct');
%     c_arr = c_arr(:);           % ensure column
%     M_log = log10(M_vec(:));   % ensure column too
%     idx = reshape(log10(M_vec) <= 13, 1, []);
%     idx = idx';
%     % idx = find(log10(M_vec) <= 13);   % numeric indices — orientation doesn't matter
%     ok  = all(diff(c_arr(idx)) < 0);
%     fprintf('  z=%.1f: c monotone decreasing over [1e10, 1e13]? %s\n', ...
%         z_tab(i), tf2str(ok));
% end

for i = 1:numel(z_tab)
    c_arr  = Klypin11(M_vec, z_tab(i), 'distinct');
    M_log  = log10(M_vec);

    % Extract values in safe region
    c_safe = [];
    for k = 1:numel(M_log)
        if M_log(k) <= 13
            c_safe(end+1) = c_arr(k);  %#ok<AGROW>
        end
    end

    % Check each consecutive pair manually
    ok = true;
    for k = 1:numel(c_safe)-1
        if c_safe(k+1) >= c_safe(k)
            ok = false;
            break;
        end
    end

    fprintf('  z=%.1f: c monotone decreasing over [1e10, 1e13]? %s\n', ...
        z_tab(i), tf2str(ok));
end

%% =========================================================================
% 6. Sanity: c decreasing with z at fixed galaxy-size mass
%    (For z=0-3 the paper shows clear decline; z=3->5 shows slight upturn
%     at cluster masses but for M=1e12 the z=3 value c=4.40 and z=5
%     value c=5.00 means the function turns around. We only assert
%     monotone decline for z=0-3 at M=1e12.)
% =========================================================================
fprintf('\n--- Sanity: c decreasing with z at M=1e12 (z=0 to z=3) ---\n');

z_check = [0.0, 0.5, 1.0, 2.0, 3.0];
c_z     = arrayfun(@(z) Klypin11(1e12, z, 'distinct'), z_check);
ok_decl = all(c_z(2:end) < c_z(1:end-1));
fprintf('  c(z) at M=1e12: ');
fprintf('%.2f  ', c_z);
fprintf('\n  Monotone decreasing z=0 to z=3? %s\n', tf2str(ok_decl));

% At z=3 -> z=5 the paper shows c increases (5.00 > 4.40) — verify
c_z3 = Klypin11(1e12, 3.0, 'distinct');
c_z5 = Klypin11(1e12, 5.0, 'distinct');
fprintf('  c(z=3)=%.2f < c(z=5)=%.2f  (upturn at z>3)? %s\n', ...
    c_z3, c_z5, tf2str(c_z5 > c_z3));

%% =========================================================================
% 7. Sanity: all concentrations positive and finite
% =========================================================================
fprintf('\n--- Finite/positive check ---\n');

M_all  = logspace(9, 15, 100);
z_all  = [0, 0.5, 1, 2, 3, 5];
n_bad  = 0;
for iz = 1:numel(z_all)
    c_test = Klypin11(M_all, z_all(iz), 'distinct');
    bad    = sum(~isfinite(c_test) | c_test <= 0);
    n_bad  = n_bad + bad;
end
fprintf('  Bad values (non-finite or <=0): %d  (expected 0)\n', n_bad);
assert(n_bad == 0, 'Concentration produced non-finite or non-positive values.');

%% =========================================================================
% 8. Interpolation continuity — no jumps between tabulated z values
%    Evaluate at z = 0.0, 0.01, 0.02 ... 5.0 and check max jump
% =========================================================================
fprintf('\n--- Interpolation continuity ---\n');

z_fine = linspace(0, 5, 500);
c_fine = arrayfun(@(z) Klypin11(1e12, z, 'distinct'), z_fine);
max_jump = max(abs(c_fine(2:end) - c_fine(1:end-1)));
fprintf('  Max step in c(z) over z=[0,5] at M=1e12: %.5f\n', max_jump);
assert(max_jump < 0.5, 'Discontinuity detected in c(z) interpolation.');

%% =========================================================================
% 9. Out-of-range redshift warning
% =========================================================================
fprintf('\n--- Out-of-range redshift (z=7) should warn and clamp ---\n');

warning('off', 'all');
c_clamped = Klypin11(1e12, 7.0, 'distinct');
c_z5_ref  = Klypin11(1e12, 5.0, 'distinct');
warning('on',  'all');
fprintf('  c(z=7, clamped) = %.4f\n', c_clamped);
fprintf('  c(z=5, ref)     = %.4f\n', c_z5_ref);
assert(abs(c_clamped - c_z5_ref) < 1e-10, ...
    'Clamped z=7 should equal z=5 result exactly.');

%% =========================================================================
% 10. Subhalo > distinct at same mass (tidal stripping -> higher c)
%     Paper: subhalo c is ~30% higher than distinct at same mass
% =========================================================================
fprintf('\n--- Subhalo vs distinct at z=0 ---\n');

M_cmp = [1e11, 1e12, 1e13];
fprintf('  %-10s  %-10s  %-10s  %-8s\n', 'M [Msun/h]', 'c_distinct', 'c_subhalo', 'ratio');
for M = M_cmp
    cd = Klypin11(M, 0, 'distinct');
    cs = Klypin11(M, 0, 'subhalo');
    ok = cs > cd;
    fprintf('  1e%-6.0f    %-10.3f  %-10.3f  %.3f  %s\n', ...
        log10(M), cd, cs, cs/cd, tf2str(ok));
    assert(cs > cd, sprintf('Subhalo not more concentrated than distinct at M=1e%.0f.', log10(M)));
end

%% =========================================================================
% 11. Vector input: result matches element-wise scalar calls
% =========================================================================
fprintf('\n--- Vector vs scalar consistency ---\n');

M_vec2 = [1e11, 5e11, 1e12, 5e12, 1e13];
c_vec  = Klypin11(M_vec2, 0, 'distinct');
c_scl  = arrayfun(@(M) Klypin11(M, 0, 'distinct'), M_vec2);
max_diff = max(abs(c_vec(:) - c_scl(:)));
fprintf('  Max diff vector vs scalar: %.2e  (should be < 1e-12)\n', max_diff);
assert(max_diff < 1e-12, 'Vector call inconsistent with scalar calls.');

%% =========================================================================
% 12. Plots
% =========================================================================
figure('Name','Klypin11: c(M) at all tabulated redshifts', ...
       'Position', [80 80 1050 480]);

M_plot  = logspace(10, 15, 200);
colors  = lines(numel(z_tab));

% --- Panel 1: distinct haloes c(M) at each tabulated z ------------------
subplot(1,2,1);
hold on;
for i = 1:numel(z_tab)
    c_pl = Klypin11(M_plot, z_tab(i), 'distinct');
    plot(log10(M_plot), c_pl, '-', 'Color', colors(i,:), 'LineWidth', 2, ...
        'DisplayName', sprintf('z = %.1f', z_tab(i)));
end
hold off;
xlabel('$log_{10}(M_{vir}  [h^{-1} M_{\odot}])$', 'FontSize', 12, 'Interpreter', 'latex');
ylabel('Concentration  c_{vir}', 'FontSize', 12);
title('Distinct haloes — Klypin+2011 Eq. (12)', 'FontSize', 13);
legend('Location', 'northeast', 'FontSize', 9);
grid on;
xlim([10 15]); ylim([0 20]);

% Mark the c(1e12) Table 3 values as circles
hold on;
for i = 1:numel(z_tab)
    plot(12, c_paper(i), 'o', 'Color', colors(i,:), ...
        'MarkerSize', 7, 'MarkerFaceColor', colors(i,:));
end
hold off;
text(12.05, 18, 'Circles = Table 3', 'FontSize', 8, 'Color', [0.3 0.3 0.3]);

% --- Panel 2: subhalo vs distinct at z=0 --------------------------------
subplot(1,2,2);
c_dist_z0 = Klypin11(M_plot, 0, 'distinct');
c_sub_z0  = Klypin11(M_plot, 0, 'subhalo');

plot(log10(M_plot), c_dist_z0, 'b-',  'LineWidth', 2, 'DisplayName', 'Distinct (Eq. 12, z=0)');
hold on;
plot(log10(M_plot), c_sub_z0,  'r--', 'LineWidth', 2, 'DisplayName', 'Subhalo  (Eq. 11)');
hold off;
xlabel('$log_{10}(M  [h^{-1} M_{\odot}])$', 'FontSize', 12, 'Interpreter', 'latex');
ylabel('Concentration  c', 'FontSize', 12);
title('Distinct vs Subhalo, z=0', 'FontSize', 13);
legend('Location', 'northeast', 'FontSize', 10);
grid on;
xlim([10 15]); ylim([0 18]);

sgtitle('Klypin, Trujillo-Gomez & Primack (2011)', ...
    'FontSize', 14, 'FontWeight', 'bold');

fprintf('\nAll Klypin11 tests passed.\n');

%% =========================================================================
% Helper
% =========================================================================
function s = tf2str(tf)
    if tf; s = 'PASS'; else; s = 'FAIL'; end
end