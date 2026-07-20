%% test_Klypin16.m
% Test suite for Klypin16_Table and Klypin16_concentration
%
% Tests:
%   1.  Path check
%   2.  Table integrity — all 8 table combinations, required fields
%   3.  Error handling — invalid formula / cosmo / mdef
%   4.  cM  parameter spot-check vs Tables A1-A4
%   5.  cnu parameter spot-check vs Tables A5-A8
%   6.  cM  analytical formula check (exact hand-computation)
%   7.  cnu analytical formula check (exact hand-computation)
%   8.  c(M) behaviour at z=0: decreasing galaxy->group range
%   9.  c(z) behaviour at fixed M=1e12: decreasing z=0->1
%  10.  cM vs cnu consistency — same qualitative trend
%  11.  Validity mask — low M and out-of-range z flagged false
%  12.  Plausible c ranges at z=0

fprintf('=== test_Klypin16 ===\n\n');
tfstr = @(x) repmat('PASS', x, 1) + repmat('FAIL', ~x, 1);

%% =========================================================================
% 1. Path check
% =========================================================================
here       = fileparts(mfilename('fullpath'));
repo       = fullfile(here, '..');
utils_path = fullfile(repo, 'src', 'utils');
concentration_path = fullfile(repo, 'src', 'concentration');
addpath(genpath(utils_path));
addpath(genpath(concentration_path));

assert(exist('Klypin16','file') == 2, ...
    sprintf('Klypin16.m not found.\nExpected in: %s', utils_path));
assert(exist('Klypin16_Table','file') == 2, ...
    sprintf('Klypin16_Table.m not found.\nExpected in: %s', utils_path));

fprintf('Paths OK.\n');
fprintf('  Klypin16_concentration : %s\n', which('Klypin16_concentration'));
fprintf('  Klypin16_Table         : %s\n', which('Klypin16_Table'));

%% =========================================================================
% 2. Table integrity — all 8 combinations exist and have correct fields
% =========================================================================
fprintf('\n--- 2. Table integrity ---\n');

combos_cm  = {{'planck13','200c'}, {'planck13','vir'}, {'bolshoi','200c'}, {'bolshoi','vir'}};
combos_cnu = {{'planck13','200c'}, {'planck13','vir'}};

for i = 1:4
    P = Klypin16_Table('cm', combos_cm{i}{1}, combos_cm{i}{2});
    assert(isfield(P,'z_bins') && isfield(P,'C0') && isfield(P,'gamma') && isfield(P,'M0'), ...
        'cM table missing fields: %s %s', combos_cm{i}{1}, combos_cm{i}{2});
    assert(numel(P.z_bins)==numel(P.C0) && numel(P.C0)==numel(P.gamma) && numel(P.gamma)==numel(P.M0), ...
        'cM table length mismatch: %s %s', combos_cm{i}{1}, combos_cm{i}{2});
    fprintf('  cM  %-10s %-5s : %d redshift nodes  OK\n', combos_cm{i}{1}, combos_cm{i}{2}, numel(P.z_bins));
end
for i = 1:2
    P = Klypin16_Table('cnu', combos_cnu{i}{1}, combos_cnu{i}{2});
    assert(isfield(P,'z_bins') && isfield(P,'a0') && isfield(P,'b0'), ...
        'cnu table missing fields: %s %s', combos_cnu{i}{1}, combos_cnu{i}{2});
    assert(numel(P.z_bins)==numel(P.a0) && numel(P.a0)==numel(P.b0), ...
        'cnu table length mismatch: %s %s', combos_cnu{i}{1}, combos_cnu{i}{2});
    fprintf('  cnu %-10s %-5s : %d redshift nodes  OK\n', combos_cnu{i}{1}, combos_cnu{i}{2}, numel(P.z_bins));
end

%% =========================================================================
% 3. Error handling
% =========================================================================
fprintf('\n--- 3. Error handling ---\n');
errors_expected = {
    {'bogus','planck13','200c'}, ...
    {'cm','badcosmo','200c'}, ...
    {'cm','planck13','baddef'}, ...
    {'cnu', 'bolshoi',  '200c'}};
for i = 1:numel(errors_expected)
    args = errors_expected{i};
    try
        Klypin16_Table(args{1}, args{2}, args{3});
        error('Should have thrown for: %s %s %s', args{1}, args{2}, args{3});
    catch ME
        fprintf('  Invalid (%s,%s,%s) throws: "%s"\n', args{1},args{2},args{3}, ME.message);
    end
end

%% =========================================================================
% 4. cM parameter spot-check at z=0 vs Tables A1-A4
% =========================================================================
fprintf('\n--- 4. cM parameter spot-check at z=0 ---\n');

ref_cm = struct();
ref_cm.planck13_200c = struct('C0',7.40,'gamma',0.120,'M0',5.5e17);
ref_cm.planck13_vir  = struct('C0',9.75,'gamma',0.110,'M0',5.0e17);
ref_cm.bolshoi_200c  = struct('C0',6.60,'gamma',0.110,'M0',2.0e18);
ref_cm.bolshoi_vir   = struct('C0',9.00,'gamma',0.100,'M0',2.0e18);

tags_cm = {'planck13_200c','planck13_vir','bolshoi_200c','bolshoi_vir'};
for i = 1:4
    parts = strsplit(tags_cm{i},'_');
    cn = parts{1}; md = parts{2};
    % if strcmp(cn,'bolshoi'), md_str = md; else, md_str = [parts{2},'_',parts{3}]; end
    % handle 200c having underscore
    sp   = strsplit(tags_cm{i},'_');
    cname = sp{1};
    if numel(sp)==3, mname=[sp{2},'_',sp{3}]; else, mname=sp{2}; end

    P   = Klypin16_Table('cm', cname, mname);
    ref = ref_cm.(tags_cm{i});

    dC0    = abs(P.C0(1)    - ref.C0)    / ref.C0;
    dgamma = abs(P.gamma(1) - ref.gamma) / ref.gamma;
    dM0    = abs(P.M0(1)    - ref.M0)    / ref.M0;

    ok = dC0<1e-6 && dgamma<1e-6 && dM0<1e-6;
    fprintf('  %-22s: C0=%.2f  gamma=%.3f  M0=%.2e  %s\n', ...
        tags_cm{i}, P.C0(1), P.gamma(1), P.M0(1), tfstr(ok));
    assert(ok, 'cM parameter mismatch for %s', tags_cm{i});
end

%% =========================================================================
% 5. cnu parameter spot-check at z=0 vs Tables A5-A8
% =========================================================================
fprintf('\n--- 5. cnu parameter spot-check at z=0 ---\n');

ref_cnu = struct();
ref_cnu.planck13_200c = struct('a0',0.40,'b0',0.278);
ref_cnu.planck13_vir  = struct('a0',0.75,'b0',0.567);

tags_cnu = {'planck13_200c','planck13_vir'};
for i = 1:2
    sp    = strsplit(tags_cnu{i},'_');
    cname = sp{1};
    if numel(sp)==3, mname=[sp{2},'_',sp{3}]; else, mname=sp{2}; end

    P   = Klypin16_Table('cnu', cname, mname);
    ref = ref_cnu.(tags_cnu{i});

    da0 = abs(P.a0(1) - ref.a0) / ref.a0;
    db0 = abs(P.b0(1) - ref.b0) / ref.b0;

    ok = da0<1e-6 && db0<1e-6;
    fprintf('  %-22s: a0=%.3f  b0=%.3f  %s\n', ...
        tags_cnu{i}, P.a0(1), P.b0(1), tfstr(ok));
    assert(ok, 'cnu parameter mismatch for %s', tags_cnu{i});
end

%% =========================================================================
% 6. cM analytical formula check
% =========================================================================
fprintf('\n--- 6. cM analytical formula check ---\n');

M_test = 5e13;
z_test = 0.0;
P = Klypin16_Table('cm', 'planck13', '200c');
C0_z    = interp1(P.z_bins, P.C0,    z_test);
gamma_z = interp1(P.z_bins, P.gamma, z_test);
M0_z    = interp1(P.z_bins, P.M0,    z_test);

c_manual = C0_z * (M_test/1e12)^(-gamma_z) * (1 + (M_test/M0_z)^0.4);

cosmo_dummy.sigmaM = @(M,z) 1.0;   % not used for cM
[c_code, ~] = Klypin16(M_test, z_test, cosmo_dummy, 'planck13', '200c', 'cM');

d_pct = abs(c_code - c_manual) / c_manual * 100;
ok = d_pct < 0.001;
fprintf('  manual=%.6f  code=%.6f  diff=%.2e%%  %s\n', c_manual, c_code, d_pct, tfstr(ok));
assert(ok);

%% =========================================================================
% 7. cnu analytical formula check
% =========================================================================
fprintf('\n--- 7. cnu analytical formula check ---\n');

sigma_test = 0.6;
z_test2    = 0.0;
P = Klypin16_Table('cnu', 'planck13', '200c');
a0_z = interp1(P.z_bins, P.a0, z_test2);
b0_z = interp1(P.z_bins, P.b0, z_test2);
x    = sigma_test / a0_z;
c_manual2 = b0_z * (1 + 7.37*x^0.75) * (1 + 0.14*x^(-2));

cosmo_cnu.sigmaM = @(M, z) sigma_test;
[c_code2, ~] = Klypin16(1e12, z_test2, cosmo_cnu, 'planck13', '200c', 'cnu');

d_pct2 = abs(c_code2 - c_manual2) / c_manual2 * 100;
ok2 = d_pct2 < 0.001;
fprintf('  manual=%.6f  code=%.6f  diff=%.2e%%  %s\n', c_manual2, c_code2, d_pct2, tfstr(ok2));
assert(ok2);

%% =========================================================================
% 8. c(M) behaviour: decreasing in galaxy-to-group range at z=0
% =========================================================================
fprintf('\n--- 8. c(M) behaviour at z=0 ---\n');

cosmo_lcdm.sigmaM = @(M,z) 0.82 .* (M ./ 3e12).^(-0.17);

M_gal2grp = logspace(11, 13, 12);

for frm = {'cM', 'cnu'}
    c_gal2grp = zeros(1, numel(M_gal2grp));
    for mi = 1:numel(M_gal2grp)
        [c_tmp, ~] = Klypin16(M_gal2grp(mi), 0, cosmo_lcdm, 'planck13', '200c', frm{1});
        c_gal2grp(mi) = c_tmp;
    end

    % Check monotone decreasing without diff — compare each consecutive pair
    ok_dec = true;
    for mi = 1:(numel(c_gal2grp)-1)
        if c_gal2grp(mi) <= c_gal2grp(mi+1)
            ok_dec = false;
            break;
        end
    end

    ok_pos = all(c_gal2grp > 0);

    fprintf('  %-4s c values: ', frm{1});
    fprintf('%.3f ', c_gal2grp);
    fprintf('\n');
    fprintf('  %-4s c decreasing 1e11->1e13 at z=0? %s\n', frm{1}, tfstr(ok_dec));
    fprintf('  %-4s all c > 0?                      %s\n', frm{1}, tfstr(ok_pos));
    assert(ok_dec,  'c not decreasing 1e11->1e13 for formula=%s', frm{1});
    assert(ok_pos,  'c has non-positive values for formula=%s',    frm{1});
end

%% =========================================================================
% 9. c(z) behaviour at fixed M=1e12
% =========================================================================
fprintf('\n--- 9. c(z) behaviour at fixed M=1e12 ---\n');

sigma_fn = @(M) 0.82 .* (M ./ 3e12).^(-0.17);
cosmo_z.sigmaM = @(M,z) sigma_fn(M) ./ (1+z);

z_vec   = [0.0, 0.5, 1.0];
M_fixed = 1e12;

for frm = {'cM', 'cnu'}
    c_z = zeros(1, numel(z_vec));
    for zi = 1:numel(z_vec)
        [c_tmp, ~] = Klypin16(M_fixed, z_vec(zi), cosmo_z, 'planck13', '200c', frm{1});
        c_z(zi) = c_tmp;
    end

    fprintf('  %-4s M=1e12:', frm{1});
    for i = 1:numel(z_vec)
        fprintf(' z=%.1f:%.3f', z_vec(i), c_z(i));
    end
    fprintf('\n');

    % (a) c at z=0 must be greater than z=0.5 for both formulae
    ok_lowz = c_z(1) > c_z(2);
    fprintf('  %-4s (a) c(z=0) > c(z=0.5)?  %s\n', frm{1}, tfstr(ok_lowz));
    assert(ok_lowz, 'c(z=0) not > c(z=0.5) for formula=%s', frm{1});

    % (b) all values above physical floor
    ok_pos = all(c_z > 1.0);
    fprintf('  %-4s (b) all c > 1 (floor)?   %s\n', frm{1}, tfstr(ok_pos));
    assert(ok_pos, 'c below physical floor for formula=%s', frm{1});

    % (c) cM must be fully monotone; cnu may have upturn at high z
    if strcmp(frm{1}, 'cM')
        ok_mono = c_z(2) > c_z(3);
        fprintf('  %-4s (c) fully monotone z=0->1? %s\n', frm{1}, tfstr(ok_mono));
        assert(ok_mono, 'cM not monotone decreasing z=0->1');
    else
        % cnu has a concentration floor — upturn at high z is expected
        [~, iz_min] = min(c_z);
        if iz_min < numel(z_vec)
            fprintf('  %-4s (c) c minimum at z=%.1f then upturn — expected cnu floor behaviour\n', ...
                frm{1}, z_vec(iz_min));
        else
            fprintf('  %-4s (c) no upturn detected\n', frm{1});
        end
    end
end

%% =========================================================================
% 10. cM vs cnu consistency — same qualitative c-M ordering
% =========================================================================
fprintf('\n--- 10. cM vs cnu qualitative consistency ---\n');

M_lo = 1e12;
M_hi = 1e14;
c_cM_lo  = Klypin16(M_lo, 0, cosmo_lcdm, 'planck13', '200c', 'cM');
c_cM_hi  = Klypin16(M_hi, 0, cosmo_lcdm, 'planck13', '200c', 'cM');
c_cnu_lo = Klypin16(M_lo, 0, cosmo_lcdm, 'planck13', '200c', 'cnu');
c_cnu_hi = Klypin16(M_hi, 0, cosmo_lcdm, 'planck13', '200c', 'cnu');

ok_cM  = c_cM_lo  > c_cM_hi;
ok_cnu = c_cnu_lo > c_cnu_hi;
fprintf('  cM:  c(1e12)=%.3f > c(1e14)=%.3f? %s\n', c_cM_lo,  c_cM_hi,  tfstr(ok_cM));
fprintf('  cnu: c(1e12)=%.3f > c(1e14)=%.3f? %s\n', c_cnu_lo, c_cnu_hi, tfstr(ok_cnu));
assert(ok_cM  && ok_cnu);

%% =========================================================================
% 11. Validity mask
% =========================================================================
fprintf('\n--- 11. Validity mask ---\n');

% Below minimum mass
[~, v_low] = Klypin16(1e9, 0, cosmo_lcdm, 'planck13', '200c', 'cM');
fprintf('  M=1e9  (below 1e10): valid=%d (expect 0)  %s\n', v_low, tfstr(~v_low));
assert(~v_low);

% Good mass, z=0
[~, v_ok] = Klypin16(1e12, 0, cosmo_lcdm, 'planck13', '200c', 'cM');
fprintf('  M=1e12, z=0       : valid=%d (expect 1)  %s\n', v_ok, tfstr(v_ok));
assert(v_ok);

% Beyond redshift range for planck13 cM (max z=5.4)
[~, v_hi] = Klypin16(1e12, 7.0, cosmo_lcdm, 'planck13', '200c', 'cM');
fprintf('  z=7.0 (> zmax=5.4): valid=%d (expect 0)  %s\n', v_hi, tfstr(~v_hi));
assert(~v_hi);

%% =========================================================================
% 12. Plausible c ranges at z=0
% =========================================================================
fprintf('\n--- 12. Plausible c ranges at z=0 (planck13, 200c) ---\n');

masses  = [1e12, 1e13, 1e14];
labels  = {'1e12','1e13','1e14'};
lo_lims = [4, 3, 2];
hi_lims = [15, 10, 8];

for frm = {'cM','cnu'}
    fprintf('\n  formula = %s:\n', frm{1});
    c_prev = Inf;
    for mi = 1:3
        [c_i, ~] = Klypin16(masses(mi), 0, cosmo_lcdm, 'planck13', '200c', frm{1});
        in_range = c_i > lo_lims(mi) && c_i < hi_lims(mi);
        dec_ok   = c_i < c_prev;
        fprintf('    c200c(M=%s) = %.2f  range=[%d,%d]? %s  decreasing? %s\n', ...
            labels{mi}, c_i, lo_lims(mi), hi_lims(mi), tfstr(in_range), tfstr(dec_ok));
        assert(in_range, 'c200c(%s) out of range for %s', labels{mi}, frm{1});
        assert(dec_ok,   'c not decreasing at M=%s for %s', labels{mi}, frm{1});
        c_prev = c_i;
    end
end

fprintf('\n=== All Klypin16 tests passed ===\n');