%% test_CHMR.m
% Test suite for schive_CHMR and thaher_CHMR
% Run from: F:\PHANTOM\tests\

% clear; clc;

%% ── Path setup ────────────────────────────────────────────────────────────
here      = fileparts(mfilename('fullpath'));
repo      = fullfile(here, '..');
chmr_path = fullfile(repo, 'src', 'chmr');
addpath(genpath(chmr_path));

%% ── Shared parameters ─────────────────────────────────────────────────────
H        = 70.2;
Omega_M0 = 0.2720;
m_eV     = 0.8e-22;
zeta_0   = 357.6746;
M_min0   = 4.4e7 * (m_eV/1e-22)^(-1.5) * (Omega_M0/0.27)^(-0.25) * (H/70)^0.5;

pass = 0;
fail = 0;

fprintf('====================================================\n');
fprintf('  CHMR Test Suite\n');
fprintf('====================================================\n\n');

%% ── Helper ────────────────────────────────────────────────────────────────
function assert_test(name, condition, tol_info)
    if nargin < 3, tol_info = ''; end
    if condition
        fprintf('  [PASS]  %s  %s\n', name, tol_info);
    else
        fprintf('  [FAIL]  %s  %s\n', name, tol_info);
    end
end

%% ════════════════════════════════════════════════════════════════════════
fprintf('── TEST 1: Scalar input at z=0 ──────────────────────────\n');
%  At z=0: a=1, zeta=zeta_0  →  Mc = 0.25*(Mh/M_min0)^(1/3)*M_min0
%  ════════════════════════════════════════════════════════════════════════
a     = 1;
zeta  = zeta_0;
Mh    = 1e9;

Mc     = schive_CHMR(Mh, a, zeta, zeta_0, M_min0);
Mc_new = thaher_CHMR(Mh, a, zeta, zeta_0, M_min0, 1.5, 2);

Mc_analytical = 0.25 * (Mh/M_min0)^(1/3) * M_min0;

assert_test('schive scalar z=0 vs analytical', abs(Mc - Mc_analytical)/Mc_analytical < 1e-10, ...
    sprintf('(Mc=%.4e, expected=%.4e)', Mc, Mc_analytical));
assert_test('thaher scalar z=0 <= schive (suppression active)', Mc_new <= Mc, ...
    sprintf('(Mc_new=%.4e, Mc=%.4e)', Mc_new, Mc));

fprintf('\n');

%% ════════════════════════════════════════════════════════════════════════
fprintf('── TEST 2: Thaher converges to Schive for Mh >> M_min0 ─────\n');
%  ════════════════════════════════════════════════════════════════════════
Mh_large = 1e13;    % far above M_min0
Mc_s  = schive_CHMR(Mh_large, a, zeta, zeta_0, M_min0);
Mc_t  = thaher_CHMR(Mh_large, a, zeta, zeta_0, M_min0, 1.5, 2);
rel_diff = abs(Mc_t - Mc_s) / Mc_s;

assert_test('thaher → schive for Mh >> M_min0 (rel diff < 1%)', rel_diff < 0.01, ...
    sprintf('(rel diff = %.2e)', rel_diff));

fprintf('\n');

%% ════════════════════════════════════════════════════════════════════════
fprintf('── TEST 3: Thaher suppression is stronger for Mh << M_min0 ─\n');
%  ════════════════════════════════════════════════════════════════════════
Mh_small = 0.01 * M_min0;
Mc_s_sm  = schive_CHMR(Mh_small, a, zeta, zeta_0, M_min0);
Mc_t_sm  = thaher_CHMR(Mh_small, a, zeta, zeta_0, M_min0, 1.5, 2);

assert_test('thaher < schive for Mh << M_min0', Mc_t_sm < Mc_s_sm, ...
    sprintf('(Mc_t=%.4e, Mc_s=%.4e)', Mc_t_sm, Mc_s_sm));

fprintf('\n');

%% ════════════════════════════════════════════════════════════════════════
fprintf('── TEST 4: Scale factor dependence (a^-0.5 scaling) ────────\n');
%  Doubling a should scale Mc by 2^(-0.5)
%  ════════════════════════════════════════════════════════════════════════
Mh   = 1e9;
a1   = 0.5;  a2 = 1.0;
Mc1  = schive_CHMR(Mh, a1, zeta, zeta_0, M_min0);
Mc2  = schive_CHMR(Mh, a2, zeta, zeta_0, M_min0);
expected_ratio = (a1/a2)^(-0.5);          % = sqrt(2)
actual_ratio   = Mc1 / Mc2;

assert_test('schive a^-0.5 scaling', abs(actual_ratio - expected_ratio)/expected_ratio < 1e-10, ...
    sprintf('(ratio=%.6f, expected=%.6f)', actual_ratio, expected_ratio));

Mc1t = thaher_CHMR(Mh, a1, zeta, zeta_0, M_min0, 1.5, 2);
Mc2t = thaher_CHMR(Mh, a2, zeta, zeta_0, M_min0, 1.5, 2);

assert_test('thaher a^-0.5 scaling', abs(Mc1t/Mc2t - expected_ratio)/expected_ratio < 1e-10, ...
    sprintf('(ratio=%.6f, expected=%.6f)', Mc1t/Mc2t, expected_ratio));

fprintf('\n');

%% ════════════════════════════════════════════════════════════════════════
fprintf('── TEST 5: Zeta dependence (zeta/zeta_0)^(1/6) scaling ─────\n');
%  ════════════════════════════════════════════════════════════════════════
zeta_test    = 2 * zeta_0;
Mc_zeta1     = schive_CHMR(Mh, a, zeta_0,    zeta_0, M_min0);
Mc_zeta2     = schive_CHMR(Mh, a, zeta_test, zeta_0, M_min0);
expected_zeta_ratio = 2^(1/6);
actual_zeta_ratio   = Mc_zeta2 / Mc_zeta1;

assert_test('schive (zeta/zeta_0)^(1/6) scaling', ...
    abs(actual_zeta_ratio - expected_zeta_ratio)/expected_zeta_ratio < 1e-10, ...
    sprintf('(ratio=%.6f, expected=%.6f)', actual_zeta_ratio, expected_zeta_ratio));

fprintf('\n');

%% ════════════════════════════════════════════════════════════════════════
fprintf('── TEST 6: Array input (vector Mh) ──────────────────────────\n');
%  ════════════════════════════════════════════════════════════════════════
Mh_arr = logspace(6, 11, 50);
Mc_arr     = schive_CHMR(Mh_arr, a, zeta, zeta_0, M_min0);
Mc_new_arr = thaher_CHMR(Mh_arr, a, zeta, zeta_0, M_min0, 1.5, 2);

assert_test('schive vector output size matches input', isequal(size(Mc_arr'), size(Mh_arr)), ...
    sprintf('(size=%s)', mat2str(size(Mc_arr))));
assert_test('thaher vector output size matches input', isequal(size(Mc_new_arr'), size(Mh_arr)), ...
    sprintf('(size=%s)', mat2str(size(Mc_new_arr))));
assert_test('all schive outputs positive', all(Mc_arr > 0));
assert_test('all thaher outputs positive', all(Mc_new_arr > 0));
assert_test('thaher <= schive element-wise for all Mh', all(Mc_new_arr <= Mc_arr + eps));

fprintf('\n');

%% ════════════════════════════════════════════════════════════════════════
fprintf('── TEST 7: thaher default parameters (eta=1.5, gamma=2) ─────\n');
%  Calling with and without explicit defaults should give identical results
%  ════════════════════════════════════════════════════════════════════════
Mc_explicit = thaher_CHMR(Mh, a, zeta, zeta_0, M_min0, 1.5, 2);
Mc_default  = thaher_CHMR(Mh, a, zeta, zeta_0, M_min0);

assert_test('thaher explicit defaults == nargin defaults', Mc_explicit == Mc_default, ...
    sprintf('(explicit=%.6e, default=%.6e)', Mc_explicit, Mc_default));

fprintf('\n');

%% ════════════════════════════════════════════════════════════════════════
fprintf('── TEST 8: gamma sensitivity (larger gamma = sharper bend) ──\n');
%  ════════════════════════════════════════════════════════════════════════
Mh_below = 0.5 * M_min0;     % below M_t, where suppression acts
Mc_g2 = thaher_CHMR(Mh_below, a, zeta, zeta_0, M_min0, 1.5, 2);
Mc_g4 = thaher_CHMR(Mh_below, a, zeta, zeta_0, M_min0, 1.5, 4);

assert_test('larger gamma gives stronger suppression below M_t', Mc_g4 < Mc_g2, ...
    sprintf('(gamma=2: %.4e, gamma=4: %.4e)', Mc_g2, Mc_g4));

fprintf('\n');


fprintf('\n');
fprintf('====================================================\n');
fprintf('  Done.\n');
fprintf('====================================================\n');