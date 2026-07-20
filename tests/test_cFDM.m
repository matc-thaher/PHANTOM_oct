% test_c_FDM.m
% Test suite for c_FDM dispatcher and suppression_factor
% Tests routing, defaults, M_half computation, suppression models,
% physical behaviour, error handling, and plots
%
% Mirrors the structure of test_c_CDM.m

here               = fileparts(mfilename('fullpath'));
repo               = fullfile(here, '..');
utils_path         = fullfile(repo, 'src', 'utils');
concentration_path = fullfile(repo, 'src', 'concentration');
suppression_path = fullfile(repo, 'src', 'suppression');
addpath(genpath(utils_path));
addpath(genpath(concentration_path));
addpath(genpath(suppression_path));

fprintf('==========================================\n');
fprintf('  c_FDM / suppression_factor test suite\n');
fprintf('==========================================\n\n');

cosmo  = cosmology('Planck18');
M0     = 1e13;       % reference mass      [M_sun]
m_fid  = 1e-22;      % fiducial boson mass [eV]
tf_lut = {'false','true'};
tfstr  = @(x) tf_lut{x+1};

% Ground-truth M_half for fiducial boson mass
M_half_fid = 3.8e10 * (m_fid / 1e-22)^(-4/3);   % = 3.8e10 M_sun

%% --- 1. Defaults: ishiyama21 CDM + laroche22 suppression -----------
fprintf('--- 1. Defaults: ishiyama21 + laroche22 ---\n');

[c_def,  Mh_def]  = c_FDM(M0, 0, m_fid, cosmo);
[c_expl, Mh_expl] = c_FDM(M0, 0, m_fid, 'laroche22', 'ishiyama21', cosmo);

assert(abs(c_def - c_expl) < 1e-10, ...
    sprintf('Default != ishiyama21+laroche22: %.6f vs %.6f', c_def, c_expl));
assert(abs(Mh_def - M_half_fid) / M_half_fid < 1e-10, ...
    sprintf('M_half wrong: %.4e  expected %.4e', Mh_def, M_half_fid));

fprintf('  c_FDM(default)            = %.4f\n', c_def);
fprintf('  c_FDM(ish21 + laroche22)  = %.4f  PASS\n', c_expl);
fprintf('  M_half(m=1e-22 eV)        = %.4e M_sun  PASS\n\n', Mh_def);

%% --- 2. M_half scales as 3.8e10 * m22^(-4/3) ------------------------
fprintf('--- 2. M_half scaling with boson mass ---\n');

m_test = [0.1e-22, 0.5e-22, 1.0e-22, 1.5e-22, 2.0e-22, 2.5e-22];
for i = 1:numel(m_test)
    m_ax          = m_test(i);
    m22           = m_ax / 1e-22;
    Mh_expected   = 3.8e10 * m22^(-4/3);
    [~, Mh_out]   = c_FDM(M0, 0, m_ax, cosmo);
    relerr        = abs(Mh_out - Mh_expected) / Mh_expected;
    assert(relerr < 1e-10, ...
        sprintf('M_half wrong for m22=%.1f: got %.4e expected %.4e', ...
                m22, Mh_out, Mh_expected));
    fprintf('  m22=%-4.1f  M_half=%.4e M_sun  PASS\n', m22, Mh_out);
end
fprintf('\n');

%% --- 3. suppression_factor: laroche22 bounds -----------------------
fprintf('--- 3. suppression_factor laroche22 in (0,1] ---\n');

M_grid = logspace(7, 15, 50);
for i = 1:numel(m_test)
    Mh  = 3.8e10 * (m_test(i)/1e-22)^(-4/3);
    sup = arrayfun(@(m) suppression_factor(m, Mh, 'laroche22'), M_grid);
    assert(all(sup > 0) && all(sup <= 1 + 1e-10), ...
        sprintf('Laroche suppression out of (0,1] for m22=%.1f', m_test(i)/1e-22));
    fprintf('  m22=%-4.1f  sup in [%.4f, %.4f]  PASS\n', ...
        m_test(i)/1e-22, min(sup), max(sup));
end
fprintf('\n');

%% --- 4. suppression_factor: dentler22 bounds -----------------------
fprintf('--- 4. suppression_factor dentler22 in (0,1] ---\n');

for i = 1:numel(m_test)
    Mh  = 3.8e10 * (m_test(i)/1e-22)^(-4/3);
    sup = arrayfun(@(m) suppression_factor(m, Mh, 'dentler22'), M_grid);
    assert(all(sup > 0) && all(sup <= 1 + 1e-10), ...
        sprintf('Dentler suppression out of (0,1] for m22=%.1f', m_test(i)/1e-22));
    fprintf('  m22=%-4.1f  sup in [%.4f, %.4f]  PASS\n', ...
        m_test(i)/1e-22, min(sup), max(sup));
end
fprintf('\n');

%% --- 5. suppression -> 1 as M -> inf (both models) ------------------
fprintf('--- 5. Suppression -> 1 as M >> M_half ---\n');

Mh_fid   = M_half_fid;
M_large  = 1e20;   % M >> M_half: suppression must saturate to 1
sup_la   = suppression_factor(M_large, Mh_fid, 'laroche22');
sup_de = suppression_factor(M_large, M_half_fid, 'dentler22');

assert(abs(sup_la - 1) < 1e-3, ...
    sprintf('Laroche suppression not -> 1 at large M: %.6f', sup_la));
assert(abs(sup_de - 1) < 1e-3, ...
    sprintf('Dentler suppression not -> 1 at large M: %.6f', sup_de));

fprintf('  laroche22  sup(M=1e20) = %.6f  PASS\n', sup_la);
fprintf('  dentler22  sup(M=1e20) = %.6f  PASS\n\n', sup_de);

%% --- 6. suppression -> 0 as M -> 0 (both models) --------------------
fprintf('--- 6. Suppression -> 0 as M << M_half ---\n');

M_tiny = 1e3;   % M << M_half
sup_la = suppression_factor(M_tiny, Mh_fid, 'laroche22');
sup_de = suppression_factor(M_tiny, Mh_fid, 'dentler22');

assert(sup_la < 0.01, ...
    sprintf('Laroche suppression not -> 0 at tiny M: %.6f', sup_la));
assert(sup_de < 0.01, ...
    sprintf('Dentler suppression not -> 0 at tiny M: %.6f', sup_de));

fprintf('  laroche22  sup(M=1e3)  = %.6e  PASS\n', sup_la);
fprintf('  dentler22  sup(M=1e3)  = %.6e  PASS\n\n', sup_de);

%% --- 7. Suppression is monotonically increasing with M ---------------
fprintf('--- 7. Suppression increases monotonically with M ---\n');

M_mono = logspace(7, 15, 200);
for i = 1:numel(m_test)
    Mh   = 3.8e10 * (m_test(i)/1e-22)^(-4/3);
    s_la = arrayfun(@(m) suppression_factor(m, Mh, 'laroche22'), M_mono);
    s_de = arrayfun(@(m) suppression_factor(m, Mh, 'dentler22'), M_mono);
    assert(all(diff(s_la) >= -1e-12), ...
        sprintf('Laroche not monotone for m22=%.1f', m_test(i)/1e-22));
    assert(all(diff(s_de) >= -1e-12), ...
        sprintf('Dentler not monotone for m22=%.1f', m_test(i)/1e-22));
    fprintf('  m22=%-4.1f  laroche monotone? %s   dentler monotone? %s  PASS\n', ...
        m_test(i)/1e-22, tfstr(all(diff(s_la)>=0)), tfstr(all(diff(s_de)>=0)));
end
fprintf('\n');

%% --- 8. Heavier boson -> less suppression at fixed M -----------------
fprintf('--- 8. Heavier boson mass -> less suppression at fixed M ---\n');

M_ref = 1e10;   % well within suppressed regime for all boson masses
sup_la_prev = 0;
sup_de_prev = 0;
for i = 1:numel(m_test)
    Mh     = 3.8e10 * (m_test(i)/1e-22)^(-4/3);
    sup_la = suppression_factor(M_ref, Mh, 'laroche22');
    sup_de = suppression_factor(M_ref, Mh, 'dentler22');
    if i > 1
        assert(sup_la > sup_la_prev, ...
            sprintf('Laroche: heavier boson did not reduce suppression (m22=%.1f)', ...
                    m_test(i)/1e-22));
        assert(sup_de > sup_de_prev, ...
            sprintf('Dentler: heavier boson did not reduce suppression (m22=%.1f)', ...
                    m_test(i)/1e-22));
    end
    sup_la_prev = sup_la;
    sup_de_prev = sup_de;
    fprintf('  m22=%-4.1f  sup_la=%.4f  sup_de=%.4f  PASS\n', ...
        m_test(i)/1e-22, sup_la, sup_de);
end
fprintf('\n');

%% --- 9. c_FDM < c_CDM at low mass (suppression is active) -----------
fprintf('--- 9. c_FDM < c_CDM at low halo mass ---\n');

M_low = 1e8;   % well below M_half
sup_models = {'laroche22', 'dentler22'};
for s = 1:numel(sup_models)
    [c_fdm_low, ~] = c_FDM(M_low, 0, m_fid, sup_models{s}, 'ishiyama21', cosmo);
    c_cdm_low      = c_CDM(M_low, 0, 'ishiyama21', cosmo);
    assert(c_fdm_low < c_cdm_low, ...
        sprintf('%s: c_FDM not < c_CDM at low mass (%.4f vs %.4f)', ...
                sup_models{s}, c_fdm_low, c_cdm_low));
    fprintf('  %-12s  c_FDM=%.4f  c_CDM=%.4f  c_FDM<c_CDM? %s  PASS\n', ...
        sup_models{s}, c_fdm_low, c_cdm_low, tfstr(c_fdm_low < c_cdm_low));
end
fprintf('\n');

%% --- 10. c_FDM -> c_CDM at high mass (suppression saturates) --------
fprintf('--- 10. c_FDM -> c_CDM at high halo mass ---\n');

M_high = 1e15;   % M >> M_half: suppression ~ 1
tol_rel = 1e-3;
for s = 1:numel(sup_models)
    [c_fdm_hi, ~] = c_FDM(M_high, 0, m_fid, sup_models{s}, 'ishiyama21', cosmo);
    c_cdm_hi      = c_CDM(M_high, 0, 'ishiyama21', cosmo);
    relerr        = abs(c_fdm_hi - c_cdm_hi) / c_cdm_hi;
    assert(relerr < tol_rel, ...
        sprintf('%s: c_FDM not close to c_CDM at high mass (relerr=%.2e)', ...
                sup_models{s}, relerr));
    fprintf('  %-12s  c_FDM=%.4f  c_CDM=%.4f  rel_err=%.2e  PASS\n', ...
        sup_models{s}, c_fdm_hi, c_cdm_hi, relerr);
end
fprintf('\n');

%% --- 11. c_FDM decreases with M at z=0 (both suppression models) ----
fprintf('--- 11. c_FDM(M) decreasing at z=0 ---\n');

M_lo = 1e11;  M_hi_test = 1e15;
for s = 1:numel(sup_models)
    [c_lo, ~] = c_FDM(M_lo,      0, m_fid, sup_models{s}, 'ishiyama21', cosmo);
    [c_hi, ~] = c_FDM(M_hi_test, 0, m_fid, sup_models{s}, 'ishiyama21', cosmo);
    assert(c_lo > c_hi, ...
        sprintf('%s: c_FDM not decreasing with M (%.3f vs %.3f)', ...
                sup_models{s}, c_lo, c_hi));
    fprintf('  %-12s  c(1e11)=%.3f  c(1e15)=%.3f  decreasing? %s  PASS\n', ...
        sup_models{s}, c_lo, c_hi, tfstr(c_lo > c_hi));
end
fprintf('\n');

%% --- 12. c_FDM decreases with z at fixed M --------------------------
fprintf('--- 12. c_FDM(z) decreasing at fixed M=1e13 ---\n');

z_lo = 0;  z_hi = 2;
for s = 1:numel(sup_models)
    [c_z0, ~] = c_FDM(M0, z_lo, m_fid, sup_models{s}, 'ishiyama21', cosmo);
    [c_z2, ~] = c_FDM(M0, z_hi, m_fid, sup_models{s}, 'ishiyama21', cosmo);
    assert(c_z0 > c_z2, ...
        sprintf('%s: c_FDM not decreasing with z (z=0:%.3f z=2:%.3f)', ...
                sup_models{s}, c_z0, c_z2));
    fprintf('  %-12s  c(z=0)=%.3f  c(z=2)=%.3f  decreasing? %s  PASS\n', ...
        sup_models{s}, c_z0, c_z2, tfstr(c_z0 > c_z2));
end
fprintf('\n');

%% --- 13. Vector input: output size matches input ---------------------
fprintf('--- 13. Vector input preserves size ---\n');

M_vec = logspace(8, 14, 25);
for s = 1:numel(sup_models)
    [c_vec, ~] = c_FDM(M_vec, 0, m_fid, sup_models{s}, 'ishiyama21', cosmo);
    assert(numel(c_vec) == numel(M_vec), ...
        sprintf('Size mismatch for %s', sup_models{s}));
    assert(all(isfinite(c_vec)) && all(c_vec > 0), ...
        sprintf('Non-finite/negative values for %s', sup_models{s}));
    fprintf('  %-12s  N=%d  c in [%.2f, %.2f]  PASS\n', ...
        sup_models{s}, numel(c_vec), min(c_vec), max(c_vec));
end
fprintf('\n');

%% --- 14. All c_CDM backends run without error via c_FDM -------------
fprintf('--- 14. All CDM backends route correctly through c_FDM ---\n');

cdm_calls = { ...
    'bullock01',  {cosmo}                            ; ...
    'duffy08',    {}                            ; ...
    'klypin11',   {}                            ; ...
    'prada12',    {cosmo}                       ; ...
    'dutton14',   {}                            ; ...
    'diemer15',   {cosmo}                       ; ...
    'ludlow16',   {cosmo}                       ; ...
    'klypin16',   {cosmo, 'planck13', '200c', 'cM'}   ; ...
    'child18',    {cosmo}                       ; ...
    'diemer19',   {cosmo}                       ; ...
    'ishiyama21', {cosmo}                       ; ...
};
for i = 1:size(cdm_calls, 1)
    mdl      = cdm_calls{i, 1};
    extra    = cdm_calls{i, 2};
    [c_out, Mh_out] = c_FDM(M0, 0, m_fid, 'laroche22', mdl, extra{:});
    assert(isscalar(c_out),   sprintf('Non-scalar from %s', mdl));
    assert(isfinite(c_out),   sprintf('Non-finite from %s', mdl));
    assert(c_out > 0,         sprintf('Non-positive from %s', mdl));
    assert(isfinite(Mh_out),  sprintf('Non-finite M_half from %s', mdl));
    fprintf('  %-12s  c=%.4f  M_half=%.3e  PASS\n', mdl, c_out, Mh_out);
end
fprintf('\n');

%% --- 15. Unknown suppression model throws clean error ---------------
fprintf('--- 15. Unknown suppression model errors correctly ---\n');
try
    c_FDM(M0, 0, m_fid, 'bad_sup_model', 'ishiyama21', cosmo);
    error('Should have thrown for unknown suppression model');
catch ME
    assert(~isempty(ME.message), 'Error message must not be empty');
    fprintf('  Caught: "%s"  PASS\n\n', ME.message(1:min(60,end)));
end

%% --- 16. Unknown CDM model throws clean error -----------------------
fprintf('--- 16. Unknown CDM model errors correctly ---\n');
try
    c_FDM(M0, 0, m_fid, 'laroche22', 'bad_cdm_model', cosmo);
    error('Should have thrown for unknown CDM model');
catch ME
    assert(~isempty(ME.message), 'Error message must not be empty');
    fprintf('  Caught: "%s"  PASS\n\n', ME.message(1:min(60,end)));
end

%% --- 17. Missing cosmo for physics-based CDM model ------------------
fprintf('--- 17. Missing cosmo errors correctly ---\n');
cdm_needs_cosmo = {'prada12','diemer15','ludlow16','klypin16', ...
                   'child18','diemer19','ishiyama21'};
for i = 1:numel(cdm_needs_cosmo)
    try
        c_FDM(M0, 0, m_fid, 'laroche22', cdm_needs_cosmo{i});
        error('Should have thrown for missing cosmo: %s', cdm_needs_cosmo{i});
    catch ME
        assert(~isempty(ME.message), ...
            sprintf('Empty error message for model %s', cdm_needs_cosmo{i}));
        fprintf('  %-12s  caught missing-cosmo error  PASS\n', cdm_needs_cosmo{i});
    end
end
fprintf('\n');

%% --- 18. Cross-model / suppression comparison table -----------------
fprintf('--- 18. Cross-model comparison table at z=0, m22=1 ---\n');

M_compare  = [1e10, 1e12, 1e13, 1e14];
cdm_subset = {'bullock01','duffy08','klypin11','dutton14','ishiyama21'};

fprintf('  %-12s  %-12s  %8s  %8s  %8s  %8s\n', ...
    'CDM model','sup model','1e10','1e12','1e13','1e14');
fprintf('  %s\n', repmat('-',1,72));
for s = 1:numel(sup_models)
    for i = 1:numel(cdm_subset)
        mdl = cdm_subset{i};
        row = zeros(1,4);
        extra = {};
        if any(strcmp(mdl, {'bullock01','prada12','diemer15','ludlow16','klypin16', ...
                             'child18','diemer19','ishiyama21'}))
            extra = {cosmo};
        end
        for j = 1:4
            [row(j), ~] = c_FDM(M_compare(j), 0, m_fid, ...
                                 sup_models{s}, mdl, extra{:});
        end
        fprintf('  %-12s  %-12s  %8.3f  %8.3f  %8.3f  %8.3f\n', ...
            mdl, sup_models{s}, row(1), row(2), row(3), row(4));
    end
end
fprintf('\n');

%% --- 19. Plot: c_FDM(M) — Laroche suppression, Ishiyama21 CDM -------
fprintf('--- 19. Plot: c_FDM(M) with Laroche (2022) suppression ---\n');

m_plot_vals = [0.1e-22, 0.5e-22, 1.0e-22, 1.5e-22, 2.0e-22, 2.5e-22];
m_labels    = {'0.1','0.5','1.0','1.5','2.0','2.5'};
M_plot      = logspace(7, 14, 120);
clrs        = lines(numel(m_plot_vals));

figure('Name','c_FDM — Laroche suppression','Position',[100 100 900 520]);
hold on;

% CDM reference
c_cdm_plot = arrayfun(@(m) c_CDM(m, 0, 'ishiyama21', cosmo), M_plot);
plot(M_plot, c_cdm_plot, 'k--', 'LineWidth', 2.5, 'DisplayName', 'CDM (Ishiyama21)');

for i = 1:numel(m_plot_vals)
    c_fdm_plot = arrayfun(@(m) c_FDM(m, 0, m_plot_vals(i), ...
                          'laroche22', 'ishiyama21', cosmo), M_plot);
    plot(M_plot, c_fdm_plot, 'Color', clrs(i,:), 'LineWidth', 1.8, ...
         'DisplayName', sprintf('m = %s \\times10^{-22} eV', m_labels{i}));
end

set(gca, 'XScale', 'log', 'YScale', 'log');
xlabel('Halo mass  $M$  [$M_{\odot}$]', 'FontSize', 12, 'Interpreter', 'latex');
ylabel('Concentration  $c$',            'FontSize', 12, 'Interpreter', 'latex');
title('FDM Concentration — Laroche (2022) Suppression', 'FontSize', 13);
legend('Location', 'northeast', 'FontSize', 8, 'NumColumns', 2);
grid on; hold off;
fprintf('  PASS\n\n');

%% --- 20. Plot: c_FDM(M) — Dentler suppression, Ishiyama21 CDM -------
fprintf('--- 20. Plot: c_FDM(M) with Dentler (2022) suppression ---\n');

figure('Name','c_FDM — Dentler suppression','Position',[120 120 900 520]);
hold on;

plot(M_plot, c_cdm_plot, 'k--', 'LineWidth', 2.5, 'DisplayName', 'CDM (Ishiyama21)');

for i = 1:numel(m_plot_vals)
    c_fdm_plot = arrayfun(@(m) c_FDM(m, 0, m_plot_vals(i), ...
                          'dentler22', 'ishiyama21', cosmo), M_plot);
    plot(M_plot, c_fdm_plot, 'Color', clrs(i,:), 'LineWidth', 1.8, ...
         'DisplayName', sprintf('m = %s \\times10^{-22} eV', m_labels{i}));
end

set(gca, 'XScale', 'log', 'YScale', 'log');
xlabel('Halo mass  $M$  [$M_{\odot}$]', 'FontSize', 12, 'Interpreter', 'latex');
ylabel('Concentration  $c$',            'FontSize', 12, 'Interpreter', 'latex');
title('FDM Concentration — Dentler (2022) Suppression', 'FontSize', 13);
legend('Location', 'northeast', 'FontSize', 8, 'NumColumns', 2);
grid on; hold off;
fprintf('  PASS\n\n');

%% --- 21. Plot: suppression factors only — Laroche vs Dentler --------
fprintf('--- 21. Plot: suppression factor comparison ---\n');

figure('Name','Suppression factor — Laroche vs Dentler','Position',[140 140 900 520]);
hold on;

% CDM reference line at Delta = 1
plot(M_plot, ones(size(M_plot)), 'k--', 'LineWidth', 2.5, ...
     'DisplayName', 'CDM  (\Delta = 1)');

for i = 1:numel(m_plot_vals)
    Mh   = 3.8e10 * (m_plot_vals(i)/1e-22)^(-4/3);
    s_la = arrayfun(@(m) suppression_factor(m, Mh, 'laroche22'), M_plot);
    s_de = arrayfun(@(m) suppression_factor(m, Mh, 'dentler22'), M_plot);
    lbl  = m_labels{i};
    % Laroche: solid
    plot(M_plot, s_la, '-',  'Color', clrs(i,:), 'LineWidth', 2.0, ...
         'DisplayName', sprintf('Laroche  %s\\times10^{-22}', lbl));
    % Dentler: dotted, same colour
    plot(M_plot, s_de, ':',  'Color', clrs(i,:), 'LineWidth', 2.2, ...
         'DisplayName', sprintf('Dentler  %s\\times10^{-22}', lbl));
end

set(gca, 'XScale', 'log');
ylim([-0.05, 1.15]);
xlabel('Halo mass  $M$  [$M_{\odot}$]', 'FontSize', 12, 'Interpreter', 'latex');
ylabel('Suppression  $\Delta$',         'FontSize', 12, 'Interpreter', 'latex');
title('Suppression Factor: Laroche (2022) vs Dentler (2022)', 'FontSize', 13);
legend('Location', 'southeast', 'FontSize', 7, 'NumColumns', 2);
grid on; hold off;
fprintf('  PASS\n\n');

%% --- 22. Plot: c_FDM(z) — both suppression models at M=1e13 ---------
fprintf('--- 22. Plot: c_FDM(z) both suppression models at M=1e13 ---\n');

z_range = linspace(0, 3, 80);

figure('Name','c_FDM vs redshift','Position',[160 160 900 520]);
hold on;

% CDM reference
c_cdm_z = arrayfun(@(z) c_CDM(M0, z, 'ishiyama21', cosmo), z_range);
plot(z_range, c_cdm_z, 'k--', 'LineWidth', 2.5, 'DisplayName', 'CDM (Ishiyama21)');

ls = {'-', ':'};
for s = 1:numel(sup_models)
    for i = 1:numel(m_plot_vals)
        c_fdm_z = arrayfun(@(z) c_FDM(M0, z, m_plot_vals(i), ...
                            sup_models{s}, 'ishiyama21', cosmo), z_range);
        lbl = sprintf('%s  m=%s\\times10^{-22}', sup_models{s}, m_labels{i});
        plot(z_range, c_fdm_z, ls{s}, 'Color', clrs(i,:), 'LineWidth', 1.8, ...
             'DisplayName', lbl);
    end
end

xlabel('Redshift  $z$', 'FontSize', 12, 'Interpreter', 'latex');
ylabel('Concentration  $c$', 'FontSize', 12, 'Interpreter', 'latex');
title('$c_{\rm FDM}(z)$ at $M = 10^{13}\,M_{\odot}$ — solid: Laroche, dotted: Dentler', ...
      'FontSize', 12, 'Interpreter', 'latex');
legend('Location', 'northeast', 'FontSize', 7, 'NumColumns', 2);
grid on; hold off;
fprintf('  PASS\n\n');

fprintf('==========================================\n');
fprintf('  All c_FDM / suppression_factor tests PASSED\n');
fprintf('==========================================\n');