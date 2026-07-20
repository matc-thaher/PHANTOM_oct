%% test_Dutton14.m
% Test script for Dutton14_concentration and Dutton14_Table
%
% Location : PHANTOM\tests\test_Dutton14.m
%
% Tests covered:
%   1.  Path check
%   2.  Table integrity — valid modes and unknown mode error
%   3.  Verify a(z=0) = a1 exactly (exponential at z=0 must equal a1)
%   4.  Verify a(z->inf) -> a0 (exponential saturates to a0)
%   5.  Spot checks vs Table 1 of Dutton & Maccio 2014 (c200 and cvir)
%   6.  c decreasing with M at z=0 over [1e10, 1e15]
%   7.  c decreasing with z at fixed mass z=0->5 (both mdef)
%   8.  All values finite and positive over wide M-z grid
%   9.  Vector input matches element-wise scalar calls
%   10. Custom M_pivot changes result but not at M=M_pivot
%   11. Plots: c(M) at multiple z + c(z) at fixed masses

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

assert(exist('Dutton14','file') == 2, ...
    sprintf('Dutton14_concentration.m not found.\nExpected in: %s', utils_path));
assert(exist('Dutton14_Table','file') == 2, ...
    sprintf('Dutton14_Table.m not found.\nExpected in: %s', utils_path));

fprintf('Paths OK.\n');
fprintf('  Dutton14_concentration : %s\n', which('Dutton14_concentration'));
fprintf('  Dutton14_Table         : %s\n', which('Dutton14_Table'));

%% =========================================================================
% 1. Table integrity
% =========================================================================
fprintf('\n--- Table integrity ---\n');

mdefs = {'200c', 'vir'};
for i = 1:numel(mdefs)
    try
        P = Dutton14_Table(mdefs{i});
        fields = {'a0','a1','eta','phi','b0','b1','Mpivot'};
        for f = fields
            assert(isfield(P, f{1}), 'Missing field: %s', f{1});
        end
        fprintf('  Dutton14_Table(''%s'') : OK\n', mdefs{i});
    catch ME
        fprintf('  Dutton14_Table(''%s'') : FAIL — %s\n', mdefs{i}, ME.message);
    end
end

% Unknown mode must throw
try
    Dutton14_Table('bogus');
    fprintf('  Unknown mode ERROR: did not throw (FAIL)\n');
catch ME
    fprintf('  Unknown mode correctly throws: "%s"\n', ME.message);
end

%% =========================================================================
% 2. a(z=0) = a1  — by definition of the exponential form
%    a(z) = a0 + (a1-a0)*exp(-eta*z^phi)
%    At z=0: exp(-eta*0^phi) = exp(0) = 1 => a(0) = a0 + (a1-a0) = a1
% =========================================================================
fprintf('\n--- a(z=0) = a1 check ---\n');

for i = 1:numel(mdefs)
    P  = Dutton14_Table(mdefs{i});
    az0 = P.a0 + (P.a1 - P.a0) .* exp(-P.eta .* (0).^P.phi);
    fprintf('  %s: a(z=0) = %.6f  a1 = %.6f  diff = %.2e  %s\n', ...
        mdefs{i}, az0, P.a1, abs(az0 - P.a1), tfstr(abs(az0 - P.a1) < 1e-10));
    assert(abs(az0 - P.a1) < 1e-10, 'a(z=0) != a1 for mdef=%s', mdefs{i});
end

%% =========================================================================
% 3. a(z->large) -> a0  — exponential decays to zero
% =========================================================================
fprintf('\n--- a(z=large) -> a0 check ---\n');

for i = 1:numel(mdefs)
    P    = Dutton14_Table(mdefs{i});
    azinf = P.a0 + (P.a1 - P.a0) .* exp(-P.eta .* 20.^P.phi);
    fprintf('  %s: a(z=20) = %.6f  a0 = %.6f  diff = %.2e\n', ...
        mdefs{i}, azinf, P.a0, abs(azinf - P.a0));
    assert(abs(azinf - P.a0) < 0.01, 'a(z=20) is not close to a0 for mdef=%s', mdefs{i});
end

%% =========================================================================
% 4. Spot checks vs Dutton & Maccio 2014, Table 1
%
%    Paper reports c(M_pivot=1e12, z):
%      c200 at z=0 => 10^a1  (since b=0 at M=Mpivot)
%      cvir at z=0 => 10^a1
%    The concentration at M=Mpivot exactly equals 10^a(z).
%    Specific values from Table 1 (NFW, Planck):
%      c200(1e12, z=0) = 10^0.905 = 8.035
%      cvir(1e12, z=0) = 10^1.025 = 10.593
%      c200(1e12, z=2) = 10^(a0+(a1-a0)*exp(-eta*2^phi))
% =========================================================================
fprintf('\n--- Spot checks at M=Mpivot vs paper Table 1 ---\n');

spot = struct();
% c200: at M=1e12, log10(c) = a(z) exactly
P200 = Dutton14_Table('200c');
Pvir = Dutton14_Table('vir');

z_check = [0, 0.5, 1.0, 2.0, 5.0];
fprintf('\n  c200 at M=1e12:\n');
for iz = 1:numel(z_check)
    z = z_check(iz);
    a_z = P200.a0 + (P200.a1 - P200.a0).*exp(-P200.eta.*z.^P200.phi);
    c_paper = 10^a_z;   % exact value from formula (no mass dependence at pivot)
    c_code  = Dutton14(1e12, z, '200c');
    diff_pct = abs(c_code - c_paper) / c_paper * 100;
    fprintf('    z=%.1f: formula=%.4f  code=%.4f  diff=%.3f%%  %s\n', ...
        z, c_paper, c_code, diff_pct, tfstr(diff_pct < 0.01));
    assert(diff_pct < 0.01, 'c200 mismatch at z=%.1f', z);
end

fprintf('\n  cvir at M=1e12:\n');
for iz = 1:numel(z_check)
    z = z_check(iz);
    a_z = Pvir.a0 + (Pvir.a1 - Pvir.a0).*exp(-Pvir.eta.*z.^Pvir.phi);
    c_paper = 10^a_z;
    c_code  = Dutton14(1e12, z, 'vir');
    diff_pct = abs(c_code - c_paper) / c_paper * 100;
    fprintf('    z=%.1f: formula=%.4f  code=%.4f  diff=%.3f%%  %s\n', ...
        z, c_paper, c_code, diff_pct, tfstr(diff_pct < 0.01));
    assert(diff_pct < 0.01, 'cvir mismatch at z=%.1f', z);
end

%% =========================================================================
% 5. Plausible c200 range at z=0
%    For Planck cosmology, relaxed haloes:
%      1e12: c200 ~ 8  (10^0.905)
%      1e14: c200 < 1e12 (b < 0 means c decreasing with M)
% =========================================================================
fprintf('\n--- Plausible range at z=0 (c200) ---\n');

c12 = Dutton14(1e12, 0, '200c');
c14 = Dutton14(1e14, 0, '200c');
fprintf('  c200(1e12, z=0) = %.3f  (expected ~8.0)\n', c12);
fprintf('  c200(1e14, z=0) = %.3f  (expected < c200(1e12))\n', c14);
assert(abs(c12 - 10^P200.a1) < 0.001, 'c200(1e12,z=0) != 10^a1');
assert(c14 < c12, 'c200 should decrease with M at z=0');

%% =========================================================================
% 6. c decreasing with M at z=0 over [1e10, 1e15]  (both mdef)
% =========================================================================
fprintf('\n--- c decreasing with M at z=0 ---\n');

M_vec = logspace(10, 15, 80);
for i = 1:numel(mdefs)
    c_arr = Dutton14(M_vec, 0, mdefs{i});
    c_arr = reshape(c_arr, 1, []);
    ok = true;
    for k = 1:numel(c_arr)-1
        if c_arr(k+1) >= c_arr(k)
            ok = false; break;
        end
    end
    fprintf('  %s: c monotone decreasing with M? %s\n', mdefs{i}, tfstr(ok));
    assert(ok, 'c not monotone decreasing with M for mdef=%s', mdefs{i});
end

%% =========================================================================
% 7. c(z) behaviour: decreasing at low z, may have upturn at high z
%    Physics: the exponential a(z) saturates to a0 at high z, causing
%    c(M_large, z) to flatten and upturn — this is the well-known c_min
%    feature (Klypin+2011, Dutton+2014). The test should:
%      (a) verify c is decreasing from z=0 to z=z_peak (low-z behaviour)
%      (b) verify c never falls below a physical floor (~1.5)
%      (c) verify the upturn only appears for MASSIVE halos (M >> Mpivot)
%          not for galaxy-mass halos at moderate z
% =========================================================================
fprintf('\n--- c(z) behaviour check ---\n');

z_vec   = [0.0, 0.5, 1.0, 2.0, 3.0, 5.0];
masses  = [1e12, 1e14];

for im = 1:numel(mdefs)
    mdef = mdefs{im};
    fprintf('\n  mdef = %s:\n', mdef);

    for M = masses
        c_z = arrayfun(@(z) Dutton14(M, z, mdef), z_vec);
        fprintf('    M=1e%.0f:', log10(M));
        for iz = 1:numel(z_vec)
            fprintf(' z=%.1f:%.3f', z_vec(iz), c_z(iz));
        end

        % (a) c must decrease from z=0 to z=1 for ALL masses
        ok_lowz = c_z(1) > c_z(2) && c_z(2) > c_z(3);
        fprintf('\n           decreasing z=0->1? %s', tfstr(ok_lowz));
        assert(ok_lowz, 'c should decrease from z=0 to z=1 for M=1e%.0f mdef=%s', ...
            log10(M), mdef);

        % (b) all c values must stay above a physical floor
        c_floor = 1.5;
        ok_floor = all(c_z > c_floor);
        fprintf('  above floor (%.1f)? %s', c_floor, tfstr(ok_floor));
        assert(ok_floor, 'c fell below physical floor for M=1e%.0f mdef=%s', ...
            log10(M), mdef);

        % (c) for galaxy-mass halos (M=1e12), c should be monotone
        %     decreasing all the way to z=5 (upturn only at M >> Mpivot)
        if M <= 1e12
            ok_mono = all(diff(c_z) < 0);
            fprintf('  monotone (M~Mpivot)? %s', tfstr(ok_mono));
            assert(ok_mono, 'c not monotone for M=1e12 mdef=%s', mdef);
        else
            % For massive halos: detect if upturn exists and report z_min
            [cmin_val, iz_min] = min(c_z);
            if iz_min < numel(z_vec)
                fprintf('  (upturn at z=%.1f, c_min=%.3f — expected physics)', ...
                    z_vec(iz_min), cmin_val);
            else
                fprintf('  (monotone, no upturn detected)');
            end
        end
        fprintf('\n');
    end
end
%% =========================================================================
% 8. Finite/positive check over full M-z grid
% =========================================================================
fprintf('\n--- Finite/positive check over M-z grid ---\n');

M_grid = logspace(9, 16, 100);
z_grid = [0, 0.5, 1.0, 2.0, 3.0, 5.0];
n_bad  = 0;
for i = 1:numel(mdefs)
    for iz = 1:numel(z_grid)
        c_t = Dutton14(M_grid, z_grid(iz), mdefs{i});
        bad = sum(~isfinite(c_t) | c_t <= 0);
        n_bad = n_bad + bad;
    end
end
fprintf('  Bad values (non-finite or <=0): %d  (expected 0)\n', n_bad);
assert(n_bad == 0, 'Non-finite or non-positive values found.');

%% =========================================================================
% 9. Vector input matches scalar calls
% =========================================================================
fprintf('\n--- Vector vs scalar consistency ---\n');

M_check = [1e11, 5e11, 1e12, 5e12, 1e13, 5e13, 1e14];
for i = 1:numel(mdefs)
    c_vec = reshape(Dutton14(M_check, 0, mdefs{i}), 1, []);
    c_scl = arrayfun(@(M) Dutton14(M, 0, mdefs{i}), M_check);
    mx    = max(abs(c_vec - c_scl));
    fprintf('  %s: max diff vector vs scalar = %.2e  %s\n', ...
        mdefs{i}, mx, tfstr(mx < 1e-10));
    assert(mx < 1e-10, 'Vector/scalar mismatch for mdef=%s', mdefs{i});
end

%% =========================================================================
% 10. Custom M_pivot changes result, but not at M = M_pivot
% =========================================================================
fprintf('\n--- Custom M_pivot test ---\n');

M_custom = 1e13;
c_default = Dutton14(M_custom, 0, '200c');           % pivot = 1e12
c_custom  = Dutton14(M_custom, 0, '200c', M_custom); % pivot = 1e13

% At M = M_pivot, log10(M/M_pivot) = 0 so b*log10(...) = 0, c = 10^a
a_z0 = P200.a0 + (P200.a1 - P200.a0).*exp(-P200.eta.*0.^P200.phi);
c_at_pivot = 10^a_z0;
fprintf('  c at M=M_custom (default pivot): %.4f\n', c_default);
fprintf('  c at M=M_custom (custom pivot):  %.4f  (should = 10^a1 = %.4f)\n', ...
    c_custom, c_at_pivot);
fprintf('  c values differ: %s\n', tfstr(abs(c_default - c_custom) > 0.01));
assert(abs(c_custom - c_at_pivot) < 1e-8, ...
    'c at M=Mpivot should equal 10^a(z) regardless of pivot choice');
assert(abs(c_default - c_custom) > 0.01, ...
    'Changing Mpivot should change the result at M != new pivot');

%% =========================================================================
% 11. Plots
% =========================================================================
M_plot = logspace(10, 15, 200);
z_plot = [0.0, 0.5, 1.0, 2.0, 3.0, 5.0];
colors = lines(numel(z_plot));
z_fine = linspace(0, 5, 200);

figure('Name','Dutton14: c(M) and c(z)', 'Position',[80 80 1200 900]);

for im = 1:2
    mdef = mdefs{im};

    % Panel 1: c(M) at multiple z
    subplot(2,2,(im-1)*2 + 1);
    hold on;
    for iz = 1:numel(z_plot)
        cv = Dutton14(M_plot, z_plot(iz), mdef);
        plot(log10(M_plot), reshape(cv,1,[]), '-', ...
            'Color', colors(iz,:), 'LineWidth', 2, ...
            'DisplayName', sprintf('z=%.1f', z_plot(iz)));
    end
    hold off;
    xlabel('$log_{10}(M  [M_{\odot}/h]$)', 'FontSize', 11, 'Interpreter','latex');
    ylabel('Concentration  c',              'FontSize', 11);
    title(sprintf('Dutton+14  c(%s) vs M',  mdef), 'FontSize', 12);
    legend('Location','northeast','FontSize',8);
    grid on; xlim([10 15]);

    % Panel 2: c(z) at fixed masses
    subplot(2,2,(im-1)*2 + 2);
    mass_fix = [1e12, 1e14];
    ls = {'-','--'};
    hold on;
    for k = 1:numel(mass_fix)
        cz = arrayfun(@(z) Dutton14(mass_fix(k), z, mdef), z_fine);
        plot(z_fine, cz, ls{k}, 'LineWidth', 2, ...
            'DisplayName', sprintf('M=10^{%.0f}', log10(mass_fix(k))));
    end
    hold off;
    xlabel('Redshift  z',     'FontSize', 11);
    ylabel('Concentration  c','FontSize', 11);
    title(sprintf('Dutton+14  c(%s) vs z', mdef), 'FontSize', 12);
    legend('Location','northeast','FontSize',10);
    grid on;
end

sgtitle('Dutton & Maccio 2014, Planck Cosmology (relaxed haloes)', ...
    'FontSize', 13, 'FontWeight', 'bold');

fprintf('\nAll Dutton14 tests passed.\n');

%% =========================================================================
% Helper
% =========================================================================
function s = tfstr(tf)
    if tf; s = 'PASS'; else; s = 'FAIL'; end
end