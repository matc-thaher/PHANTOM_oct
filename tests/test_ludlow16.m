%% test_Ludlow16.m
% Test suite for Ludlow16_concentration and helper functions
%
% Tests:
%   1.  Path check
%   2.  CMH shape — Mcoll(0)=M0 and Mcoll decreases with z
%   3.  CMH EPS formula exact check at one point
%   4.  rho_char formula exact check
%   5.  formation_z bracket check — z_form > z_obs
%   6.  c output is positive and finite
%   7.  c(M) decreasing at z=0 (CDM monotone)
%   8.  c(z) decreasing at fixed M (CDM monotone)
%   9.  f dependence — lower f => earlier z_form => higher c
%  10.  C_cal sensitivity — larger C => higher c (sanity)
%  11.  Plausible c range at z=0

fprintf('=== test_Ludlow16 ===\n\n');
tfstr = @(x) repmat('PASS', x, 1) + repmat('FAIL', ~x, 1);

% ---- cosmology ----------------------------------------------------------
cosmo = cosmology('Planck18');   % real struct with all handles attached
rhoFn = @(z) cosmo.rho_crit0 .* (0.3*(1+z)^3 + 0.7); % flat LCDM mock
cosmo.rhocrit = rhoFn;

%% =========================================================================
% 1. Path check
% =========================================================================
here     = fileparts(mfilename('fullpath'));   % PHANTOM\tests
repo     = fullfile(here, '..');              % PHANTOM\
utils_path = fullfile(repo, 'src', 'utils');
concentration_path = fullfile(repo, 'src', 'concentration');
addpath(genpath(utils_path));
addpath(genpath(concentration_path));

assert(exist('Ludlow16', 'file') == 2, ...
    sprintf('Ludlow16.m not found.\nExpected in: %s', utils_path));

fprintf('Paths OK.\n');
fprintf('  Ludlow16_concentration : %s\n', which('Ludlow16_concentration'));

%% =========================================================================
% 2. CMH shape
% =========================================================================
fprintf('\n--- 2. CMH shape ---\n');
M0   = 1e12;
f    = 0.02;
z_vec = [0.0, 0.5, 1.0, 2.0, 5.0];

Mcoll = zeros(1, numel(z_vec));
for i = 1:numel(z_vec)
    Mcoll(i) = Ludlow16_CMH(z_vec(i), M0, f, cosmo, true);
end

fprintf('  z:     ');  fprintf(' %.2f', z_vec);  fprintf('\n');
fprintf('  Mcoll: ');  fprintf(' %.3e', Mcoll);  fprintf('\n');

% At z=0: delta_sc(z0)-delta_sc(z0)=0, erfc(0)=1 => Mcoll=M0
ok_z0 = abs(Mcoll(1) - M0) / M0 < 1e-6;
fprintf('  Mcoll(z=0) = M0? %s\n', tfstr(ok_z0));
assert(ok_z0);

% CMH increases with z (more mass collapsed at earlier times)
ok_mono = true;
for i = 1:(numel(Mcoll)-1)
    if Mcoll(i) <= Mcoll(i+1)
        ok_mono = false;
        break;
    end
end
fprintf('  Mcoll decreases with z? %s\n', tfstr(ok_mono));
assert(ok_mono, 'CMH not monotone decreasing with z');

%% =========================================================================
% 3. CMH EPS formula exact check
% =========================================================================
fprintf('\n--- 3. CMH exact formula check ---\n');

z_test   = 1.0;
M0_test  = 1e12;
f_test   = 0.02;

D0   = cosmo.D(0);
Dz   = cosmo.D(z_test);
d0   = 1.686 / D0;
dz   = 1.686 / Dz;
sfl  = cosmo.sigmaM(f_test * M0_test, 0);
sM0  = cosmo.sigmaM(M0_test, 0);
denom = sqrt(2 * (sfl^2 - sM0^2));

Mcoll_manual = M0_test * erfc((dz - d0) / denom);
Mcoll_code   = Ludlow16_CMH(z_test, M0_test, f_test, cosmo, true);

d_pct = abs(Mcoll_code - Mcoll_manual) / Mcoll_manual * 100;
ok    = d_pct < 0.001;
fprintf('  manual=%.6e  code=%.6e  diff=%.2e%%  %s\n', ...
    Mcoll_manual, Mcoll_code, d_pct, tfstr(ok));
assert(ok);

%% =========================================================================
% 4. rho_char formula exact check
% =========================================================================
fprintf('\n--- 4. rho_char exact formula check ---\n');

c_test    = 8.0;
z_rho     = 0.0;
fc_manual = log(1 + c_test) - c_test / (1 + c_test);
rho_manual = (200/3) * cosmo.rhocrit(z_rho) * c_test^3 / fc_manual;
rho_code   = Ludlow16_rho_char(c_test, z_rho, cosmo);

d_pct2 = abs(rho_code - rho_manual) / rho_manual * 100;
ok2    = d_pct2 < 0.001;
fprintf('  c=%.1f: manual=%.4e  code=%.4e  diff=%.2e%%  %s\n', ...
    c_test, rho_manual, rho_code, d_pct2, tfstr(ok2));
assert(ok2);

%% =========================================================================
% 5. formation_z — z_form > z_obs always
% =========================================================================
fprintf('\n--- 5. Formation redshift bracket check ---\n');

M_test_vec = [1e10, 1e12, 1e14];
for mi = 1:numel(M_test_vec)
    zf = Ludlow16_formation_z(M_test_vec(mi), 0.02, cosmo);
    ok_bracket = zf > 0;
    fprintf('  M=1e%.0f: z_form=%.3f > 0? %s\n', ...
        log10(M_test_vec(mi)), zf, tfstr(ok_bracket));
    assert(ok_bracket, 'z_form not positive for M=1e%.0f', log10(M_test_vec(mi)));
end

%% =========================================================================
% 6. c is positive and finite for a range of masses and redshifts
% =========================================================================
fprintf('\n--- 6. c positive and finite ---\n');

M_range = logspace(10, 14, 8);
z_range = [0.0, 0.5, 1.0];

for zi = 1:numel(z_range)
    for mi = 1:numel(M_range)
        [c_i, ~] = Ludlow16(M_range(mi), z_range(zi), cosmo);
        ok_i = isfinite(c_i) && c_i > 0;
        if ~ok_i
            fprintf('  FAIL: M=1e%.1f z=%.1f c=%.3f\n', ...
                log10(M_range(mi)), z_range(zi), c_i);
        end
        assert(ok_i, 'c not positive/finite at M=1e%.1f z=%.1f', ...
            log10(M_range(mi)), z_range(zi));
    end
end
fprintf('  All (M,z) combinations positive and finite: PASS\n');

%% =========================================================================
% 7. c(M) decreasing at z=0 (CDM monotone)
% =========================================================================
fprintf('\n--- 7. c(M) decreasing at z=0 ---\n');

M_vec = logspace(10, 14, 10);
c_vec = zeros(1, numel(M_vec));
for mi = 1:numel(M_vec)
    [c_vec(mi), ~] = Ludlow16(M_vec(mi), 0, cosmo);
end

fprintf('  M [Msun/h]   c\n');
for mi = 1:numel(M_vec)
    fprintf('  1e%.2f       %.3f\n', log10(M_vec(mi)), c_vec(mi));
end

ok_mono = true;
for mi = 1:(numel(c_vec)-1)
    if c_vec(mi) <= c_vec(mi+1)
        ok_mono = false;
        break;
    end
end
fprintf('  c decreasing with M at z=0? %s\n', tfstr(ok_mono));
assert(ok_mono, 'c not monotone decreasing with M at z=0');

%% =========================================================================
% 8. c(z) decreasing at fixed M=1e12
% =========================================================================
fprintf('\n--- 8. c(z) decreasing at fixed M=1e12 ---\n');

z_vec2  = [0.0, 0.5, 1.0];
M_fixed = 1e12;
c_z = zeros(1, numel(z_vec2));
zf_z = zeros(1, numel(z_vec2));
for zi = 1:numel(z_vec2)
    [c_z(zi), zf_z(zi)] = Ludlow16(M_fixed, z_vec2(zi), cosmo);
end

fprintf('  z:      ');  fprintf(' %.2f', z_vec2);    fprintf('\n');
fprintf('  c:      ');  fprintf(' %.3f', c_z);        fprintf('\n');
fprintf('  z_form: ');  fprintf(' %.3f', zf_z);       fprintf('\n');

ok_z = true;
for zi = 1:(numel(c_z)-1)
    if c_z(zi) <= c_z(zi+1)
        ok_z = false;
        break;
    end
end
fprintf('  c decreasing with z? %s\n', tfstr(ok_z));
assert(ok_z, 'c not decreasing with z for M=1e12');

%% =========================================================================
% 9. f dependence — lower f => earlier z_form => higher c
% =========================================================================
fprintf('--- 9. f dependence ---\n');
fprintf('  %6s  %8s  %6s\n', 'f', 'z_form', 'c');

f_vals  = [0.10, 0.05, 0.02, 0.01, 0.005];   % DESCENDING — large f first
M0      = 1e12;
z0      = 0;
zf_prev = -1;                                  % sentinel: any real z_form > -1

for i = 1:numel(f_vals)
    f       = f_vals(i);
    [c_val, zf] = Ludlow16(M0, z0, cosmo, f);
    fprintf('  %6.3f  %8.3f  %6.3f\n', f, zf, c_val);

    % As f decreases, z_form must increase (earlier formation time)
    assert(zf > zf_prev, ...
        sprintf('z_form not increasing as f decreases (f=%.3f gave z_form=%.3f, prev=%.3f)', ...
                f, zf, zf_prev));
    zf_prev = zf;
end
fprintf('  PASS: z_form increases monotonically as f decreases\n');

%% =========================================================================
% 10. Plausible c range at z=0
% =========================================================================
fprintf('\n--- 10. Plausible c range at z=0 ---\n');

masses = [1e12, 1e13, 1e14];
labels = {'1e12','1e13','1e14'};
lo     = [4, 3, 2];
hi     = [20, 15, 12];

c_prev2 = Inf;
for mi = 1:numel(masses)
    [c_mi, zf_mi] = Ludlow16(masses(mi), 0, cosmo);
    in_range = c_mi > lo(mi) && c_mi < hi(mi);
    dec_ok   = c_mi < c_prev2;
    fprintf('  c(M=%s, z=0) = %.2f  z_form=%.2f  range=[%d,%d]? %s  dec? %s\n', ...
        labels{mi}, c_mi, zf_mi, lo(mi), hi(mi), tfstr(in_range), tfstr(dec_ok));
    assert(in_range, 'c out of plausible range for M=%s', labels{mi});
    assert(dec_ok,   'c not decreasing at M=%s', labels{mi});
    c_prev2 = c_mi;
end

%% =========================================================================
% 11. Appendix C fit — consistency with main model and plausible range
% =========================================================================
fprintf('\n--- 11. Appendix C polynomial fit ---\n');

assert(~isempty(which('Ludlow16_fit')), ...
    'Ludlow16_concentration_fit not found');

M_vec_fit = logspace(10, 14, 8);
fprintf('  M [Msun/h]   c_fit    c_analytic\n');
c_prev_fit = Inf;

for mi = 1:numel(M_vec_fit)
    c_fit = Ludlow16_fit(M_vec_fit(mi), 0, cosmo);
    [c_ana, ~] = Ludlow16_concentration(M_vec_fit(mi), 0, cosmo);
    fprintf('  1e%.2f       %.3f    %.3f\n', ...
        log10(M_vec_fit(mi)), c_fit, c_ana);

    assert(c_fit > 0 && isfinite(c_fit), ...
        'c_fit not positive/finite at M=1e%.1f', log10(M_vec_fit(mi)));

    % monotone decreasing
    if mi > 1
        assert(c_fit < c_prev_fit, ...
            'c_fit not decreasing with M at index %d', mi);
    end
    c_prev_fit = c_fit;
end
fprintf('  Appendix C fit positive, finite, monotone: PASS\n');

% Redshift dependence: c_fit(z=0) > c_fit(z=1) at fixed M
c_fit_z0 = Ludlow16_fit(1e12, 0.0, cosmo);
c_fit_z1 = Ludlow16_fit(1e12, 1.0, cosmo);
ok_z_fit  = c_fit_z0 > c_fit_z1;
fprintf('  c_fit(z=0)=%.3f > c_fit(z=1)=%.3f? %s\n', ...
    c_fit_z0, c_fit_z1, tfstr(ok_z_fit));
assert(ok_z_fit);

fprintf('\n=== All Ludlow16 tests passed ===\n');