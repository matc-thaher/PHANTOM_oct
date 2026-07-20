% test_c_CDM.m
% Test suite for the c_CDM concentration dispatcher
% Tests routing, defaults, error handling, and basic physical behaviour

here     = fileparts(mfilename('fullpath'));   % PHANTOM\tests
repo     = fullfile(here, '..');              % PHANTOM\
utils_path = fullfile(repo, 'src', 'utils');
concentration_path = fullfile(repo, 'src', 'concentration');
addpath(genpath(utils_path));
addpath(genpath(concentration_path));

fprintf('==========================================\n');
fprintf('  c_CDM dispatcher test suite\n');
fprintf('==========================================\n\n');

cosmo  = cosmology('Planck18');
M0     = 1e13;   % reference mass  [Msun/h]
tf_lut = {'false','true'};
tfstr  = @(x) tf_lut{x+1};

%% --- 1. Default model is ishiyama21 ----------------------------------
fprintf('--- 1. Default model = ishiyama21 ---\n');
c_def    = c_CDM(M0, 0, [],          cosmo);
c_ish    = c_CDM(M0, 0, 'ishiyama21', cosmo);
assert(abs(c_def - c_ish) < 1e-10, ...
    sprintf('Default != ishiyama21: %.6f vs %.6f', c_def, c_ish));
fprintf('  c_CDM(M,z,[])  = %.4f\n', c_def);
fprintf('  c_CDM(ishiyama21) = %.4f  PASS\n\n', c_ish);

%% --- 2. All models route without error (scalar output) ---------------
fprintf('--- 2. All models return finite positive scalar ---\n');

model_calls = { ...
    'bullock01',  {cosmo},                     ; ...
    'duffy08',    {'200c_NFW_full_z0_2'},      ; ...
    'duffy08',    {'200c_NFW_relaxed_z0_2'},   ; ...
    'klypin11',   {'distinct'},                ; ...
    'klypin11',   {'subhalo'},                 ; ...
    'prada12',    {cosmo},                     ; ...
    'dutton14',   {'200c'},                    ; ...
    'dutton14',   {'vir'},                     ; ...
    'diemer15',   {cosmo},                     ; ...
    'diemer15',   {cosmo, 'mean'},              ; ...
    'ludlow16',   {cosmo},                     ; ...
    'klypin16',   {cosmo, 'planck13', '200c', 'cM'},  ; ...
    'klypin16',   {cosmo, 'planck13', 'vir',  'cnu'}, ; ...
    'child18',    {cosmo},                     ; ...
    'child18',    {cosmo, 'individual_relaxed'}; ...
    'child18',    {cosmo, 'stack_nfw'},         ; ...
    'child18',    {cosmo, 'stack_einasto'},     ; ...
    'diemer19',   {cosmo},                     ; ...
    'diemer19',   {cosmo, 'mean'},              ; ...
    'ishiyama21', {cosmo},                     ; ...
    'ishiyama21', {cosmo, '200c_relaxed'},      ; ...
    'ishiyama21', {cosmo, '500_all'},           ; ...
    'ishiyama21', {cosmo, 'vir_relaxed'},       ; ...
};

for i = 1:size(model_calls, 1)
    mdl  = model_calls{i,1};
    args = model_calls{i,2};
    c    = c_CDM(M0, 0, mdl, args{:});
    assert(isscalar(c),    sprintf('Non-scalar output: %s', mdl));
    assert(isfinite(c),    sprintf('Non-finite output: %s', mdl));
    assert(c > 0,          sprintf('Non-positive output: %s', mdl));
    if numel(args) >= 2
        fprintf('  %-12s  mode=%-20s  c=%.3f  PASS\n', mdl, args{2}, c);
    else
        fprintf('  %-12s  (default mode)               c=%.3f  PASS\n', mdl, c);
    end
end
fprintf('\n');

%% --- 3. Short aliases route to same result as full names -------------
fprintf('--- 3. Short aliases ---\n');
alias_pairs = { ...
    'bullock01', 'b01'  ; ...
    'duffy08',   'd08'  ; ...
    'klypin11',  'k11'  ; ...
    'prada12',   'p12'  ; ...
    'dutton14',  'dm14' ; ...
    'diemer15',  'dk15' ; ...
    'ludlow16',  'l16'  ; ...
    'klypin16',  'k16'  ; ...
    'child18',   'c18'  ; ...
    'diemer19',  'dj19' ; ...
    'ishiyama21','ish21'; ...
};
% Replace the entire alias loop body:
for i = 1:size(alias_pairs,1)
    full_name  = alias_pairs{i,1};
    short_name = alias_pairs{i,2};
    c_full  = c_CDM(M0, 0, full_name,  cosmo);
    c_short = c_CDM(M0, 0, short_name, cosmo);
    assert(abs(c_full - c_short) < 1e-10, ...
        sprintf('Alias mismatch: %s vs %s (%.6f vs %.6f)', ...
                full_name, short_name, c_full, c_short));
    fprintf('  %-12s == %-6s  c=%.4f  PASS\n', full_name, short_name, c_full);
end
fprintf('\n');

%% --- 4. Unknown model throws a clean error ---------------------------
fprintf('--- 4. Unknown model errors correctly ---\n');
try
    c_CDM(M0, 0, 'bad_model', cosmo);
    error('Should have thrown for unknown model');
catch ME
    assert(~isempty(ME.message), 'Error message must not be empty');
    fprintf('  Caught: "%s"  PASS\n\n', ME.message(1:min(60,end)));
end

%% --- 5. Missing cosmo throws a clean error ---------------------------
fprintf('--- 5. Missing cosmo errors correctly ---\n');
models_needing_cosmo = {'diemer15','ludlow16','klypin16','bullock01', ...
                        'diemer19','ishiyama21','child18', 'prada12'};
for i = 1:numel(models_needing_cosmo)
    try
        c_CDM(M0, 0, models_needing_cosmo{i});
        error('Should have thrown for missing cosmo in %s', models_needing_cosmo{i});
    catch ME
        assert(~isempty(ME.message), ...
            sprintf('Empty error message for model %s', models_needing_cosmo{i}));
        fprintf('  %-12s  caught missing-cosmo error  PASS\n', models_needing_cosmo{i});
    end
end
fprintf('\n');

%% --- 6. Vector input: output size matches input ----------------------
fprintf('--- 6. Vector input preserves size ---\n');
M_vec   = logspace(11, 15, 20);
models_vec = {'bullock01','diemer19','ishiyama21','child18','ludlow16'};
for i = 1:numel(models_vec)
    mdl = models_vec{i};
    c_vec = c_CDM(M_vec, 0, mdl, cosmo);
    assert(numel(c_vec) == numel(M_vec), ...
        sprintf('Size mismatch for %s', mdl));
    assert(all(isfinite(c_vec)) && all(c_vec > 0), ...
        sprintf('Non-finite/negative values for %s', mdl));
    fprintf('  %-12s  N=%d  c in [%.2f, %.2f]  PASS\n', ...
        mdl, numel(c_vec), min(c_vec), max(c_vec));
end
fprintf('\n');

%% --- 7. c decreases with M at z=0 (all models) ----------------------
fprintf('--- 7. c(M) decreasing at z=0 ---\n');
M_lo = 1e11;  M_hi = 1e15;
all_models = {'bullock01','duffy08','klypin11','prada12','dutton14', ...
              'diemer15','ludlow16','klypin16','child18','diemer19','ishiyama21'};
for i = 1:numel(all_models)
    mdl = all_models{i};
        c_lo = c_CDM(M_lo, 0, mdl, cosmo);
        c_hi = c_CDM(M_hi, 0, mdl, cosmo);
    assert(c_lo > c_hi, ...
        sprintf('c not decreasing with M for %s (c_lo=%.3f c_hi=%.3f)', ...
                mdl, c_lo, c_hi));
    fprintf('  %-12s  c(1e11)=%.3f  c(1e15)=%.3f  decreasing? %s  PASS\n', ...
        mdl, c_lo, c_hi, tfstr(c_lo > c_hi));
end
fprintf('\n');

%% --- 8. c decreases with z at fixed M (all models) ------------------
fprintf('--- 8. c(z) decreasing at fixed M=1e13 ---\n');
z_lo = 0;  z_hi = 1;
for i = 1:numel(all_models)
    mdl = all_models{i};
    c_z0 = c_CDM(M0, z_lo, mdl, cosmo);
    c_z2 = c_CDM(M0, z_hi, mdl, cosmo);
    assert(c_z0 > c_z2, ...
        sprintf('c not decreasing with z for %s (z=0:%.3f z=2:%.3f)', ...
                mdl, c_z0, c_z2));
    fprintf('  %-12s  c(z=0)=%.3f  c(z=2)=%.3f  decreasing? %s  PASS\n', ...
        mdl, c_z0, c_z2, tfstr(c_z0 > c_z2));
end
fprintf('\n');

%% --- 9. Cross-model comparison table at z=0 -------------------------
fprintf('--- 9. Cross-model comparison at z=0 ---\n');
M_compare = [1e12, 1e13, 1e14];
fprintf('  %-12s  %10s  %10s  %10s\n', 'model', '1e12','1e13','1e14');
fprintf('  %s\n', repmat('-',1,48));
for i = 1:numel(all_models)
    mdl = all_models{i};
    row = zeros(1,3);
    for j = 1:3
        row(j) = c_CDM(M_compare(j), 0, mdl, cosmo);
    end
    fprintf('  %-12s  %10.3f  %10.3f  %10.3f\n', mdl, row(1), row(2), row(3));
end
fprintf('\n');

%% --- 10. Plot: c(M) all models at z=0 -------------------------------
fprintf('--- 10. Plot: c(M) all models at z=0 ---\n');
M_plot  = logspace(9, 15.5, 100);
colors  = lines(numel(all_models));

figure('Name','c_CDM dispatcher — c(M) z=0','Position',[100 100 900 520]);
hold on;
for i = 1:numel(all_models)
    mdl = all_models{i};
    c_plot = arrayfun(@(m) c_CDM(m, 0, mdl, cosmo), M_plot);
    plot(M_plot, c_plot, 'Color', colors(i,:), 'LineWidth', 1.8, ...
         'DisplayName', mdl);
end
set(gca,'XScale','log');
xlabel('$M_{200c}  [M_{\odot} h^{-1}]$','FontSize',12,'Interpreter','latex');
ylabel('Concentration  c','FontSize',12);
title('c\_CDM dispatcher — all models at z=0','FontSize',13);
legend('Location','northeast','FontSize',8,'NumColumns',2);
grid on; hold off;
fprintf('  PASS\n\n');

%% --- 11. Plot: c(z) all models at M=1e13 ----------------------------
fprintf('--- 11. Plot: c(z) all models at M=1e13 ---\n');
z_range = linspace(0, 3, 80);

figure('Name','c_CDM dispatcher — c(z) M=1e13','Position',[150 150 900 520]);
hold on;
for i = 1:numel(all_models)
    mdl = all_models{i};
    c_plot = arrayfun(@(z) c_CDM(M0, z, mdl, cosmo), z_range);
    
    plot(z_range, c_plot, 'Color', colors(i,:), 'LineWidth', 1.8, ...
         'DisplayName', mdl);
end
xlabel('Redshift  z','FontSize',12);
ylabel('Concentration  c','FontSize',12);
title('c\_CDM dispatcher — all models at $M=10^{13} M_{\odot} h^{-1}$','FontSize',13, 'Interpreter','latex');
legend('Location','northeast','FontSize',8,'NumColumns',2);
grid on; hold off;
fprintf('  PASS\n\n');

fprintf('==========================================\n');
fprintf('  All c_CDM dispatcher tests PASSED\n');
fprintf('==========================================\n');