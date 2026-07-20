% RUN_SOLITON_NFW_ANALYTIC
%   Driver script: given a soliton central density rho_c and boson mass m22,
%   derive all FDM soliton + NFW halo quantities analytically and build the
%   composite profile.
%
%   What this script does, step by step:
%     1. Set rho_c and m22  (the two physical inputs)
%     2. Call soliton_obs to get rc, Mc, and Mh
%     3. Call halo_obs to get Rvir, c_nfw (via c_CDM/Ishiyama21)
%     4. Build the radial grid from rc/100 to Rvir
%     5. Call soliton_nfw_analytic to construct the composite profile
%
% --------------------------------------------------------------------------
% UNITS THROUGHOUT
%   mass      : M_sun
%   length    : kpc
%   velocity  : km/s
%   density   : M_sun / kpc^3
% --------------------------------------------------------------------------
% REFERENCES
%   Schive, Chiueh & Broadhurst 2014, Nat. Phys. 10, 496   — soliton profile
%   Schive+2014, PRL 113 261302                             — CHMR, Eq.(7)
%   Robles, Bullock & Boylan-Kolchin 2018, MNRAS            — observables
%   Bryan & Norman 1998, ApJ 495, 80                        — Delta_vir
%   Ishiyama et al. 2021, MNRAS 506, 4210                   — concentration
% --------------------------------------------------------------------------

clear;

% ==========================================================================
% PHANTOM PATH  (edit phantom_root or leave '' if already on path)
% ==========================================================================
phantom_root = '../';   % e.g. '/home/you/PHANTOM'

if ~isempty(phantom_root)
    subdirs = {'src/concentration','src/utils','src/profiles', ...
               'src/chmr','src/fdm', 'src/halo'};
    for k = 1:numel(subdirs)
        p = fullfile(phantom_root, subdirs{k});
        if exist(p, 'dir') && isempty(strfind(path(), p))
            addpath(p);
        end
    end
end

% ==========================================================================
% USER INPUTS
% ==========================================================================

% Boson mass parameter  m_22 = m / (10^-22 eV)
m22 = 0.5;

% Soliton central density  [M_sun / kpc^3]
rho_c = 1.0e10;

% Redshift
z = 0.0;

% Cosmology  (PHANTOM struct)
cosmo = cosmology('uchuu');

% Overdensity convention for virial quantities
overdensity = 'bryan_norman';   % 'bryan_norman' | '200c' | '200m' | '500c'

% Concentration model name for c_CDM  (string or numeric fixed value)
%   '' or NaN  -> default 'ishiyama21'
%   'dutton14','diemer19', etc.  -> any c_CDM model
%   8.5        -> fixed numeric, bypass model
concentration_choice = '';

% ==========================================================================
% SOLITON OBSERVABLES FROM rho_c
%     Uses soliton_obs (this thread) with 'rho_c' as input.
%
%     Equations (Robles+2018, Schive+2014 Nat.Phys., Schive+2014 PRL):
%       rc    : rho_c = 1.93e7 * m22^2 * (rc/kpc)^-4  [Robles+2018 Eq.(3)]
%       Mc    : Mc    = 5.5e9  * m22^-2 * (rc/kpc)^-1  [Schive+2014 Nat.Phys. Eq.(S3)]
%       Mh    : rc    = 1.6 * m22^-1 * a^(1/2) * (Mh/1e9)^(-1/3)  [Schive+2014 PRL Eq.(7)]
% ==========================================================================

a = 1.0 / (1.0 + z);   % scale factor

rc_sol = soliton_obs('rho_c', rho_c, 'r_c',  m22);         % [kpc]
Mc_sol = soliton_obs('rho_c', rho_c, 'M_c',  m22);         % [M_sun]
Mh     = soliton_obs('rho_c', rho_c, 'M_h',  m22, a);      % [M_sun]
V_c    = soliton_obs('rho_c', rho_c, 'V_c',  m22);         % [km/s]
E_tot  = soliton_obs('rho_c', rho_c, 'E_tot',m22);         % [M_sun (km/s)^2]

fprintf('\n=== Soliton observables from rho_c = %.3e M_sun/kpc^3, m22 = %.2f ===\n', ...
        rho_c, m22);
fprintf('  rc_sol  = %.4f  kpc\n',          rc_sol);
fprintf('  Mc_sol  = %.4e  M_sun\n',        Mc_sol);
fprintf('  Mh      = %.4e  M_sun\n',        Mh);
fprintf('  V_c     = %.4f  km/s\n',         V_c);
fprintf('  E_tot   = %.4e  M_sun (km/s)^2\n', E_tot);

% ==========================================================================
% HALO OBSERVABLES FROM Mh
%     Uses halo_obs (this thread).
%
%     Equations:
%       Rvir  : Mh = (4*pi/3)*Delta*rho_crit*Rvir^3  [B&N 1998]
%       c_nfw : c_CDM with Ishiyama+2021 (default)
%       rs    : rs = Rvir / c
% ==========================================================================

opts_halo.overdensity   = overdensity;
opts_halo.concentration = concentration_choice;

% Rvir  = halo_obs('M_h', Mh, 'r_vir',        cosmo, z, opts_halo);    % [kpc]
% c_nfw = halo_obs('M_h', Mh, 'concentration', cosmo, z, opts_halo);   % [1]
% Vvir  = halo_obs('M_h', Mh, 'V_vir',         cosmo, z, opts_halo);   % [km/s]
% rs    = halo_obs('M_h', Mh, 'r_s',           cosmo, z, opts_halo);   % [kpc]

[halo, hc] = halo_obs('M_h', Mh, ...
    {'r_vir','concentration','V_vir','r_s'}, cosmo, z, opts_halo);

Rvir  = halo.r_vir;
c_nfw = halo.concentration;
Vvir  = halo.V_vir;
rs    = halo.r_s;

fprintf('\n=== Halo observables from Mh = %.3e M_sun ===\n', Mh);
fprintf('  Rvir    = %.2f   kpc\n',  Rvir);
fprintf('  c_nfw   = %.3f\n',        c_nfw);
fprintf('  Vvir    = %.3f   km/s\n', Vvir);
fprintf('  rs      = %.4f   kpc\n',  rs);

% ==========================================================================
% RADIAL GRID
%     Start at rc/100 (well inside the soliton core) up to Rvir.
%     Log-spaced so both the core and the outer NFW envelope are resolved.
% ==========================================================================

r_min = rc_sol / 15;    % [kpc] — deep inside soliton core
r_max = Rvir;            % [kpc] — virial radius

N_pts = 200;
r     = logspace(log10(r_min), log10(r_max), N_pts)';   % [kpc], column

fprintf('\n=== Radial grid: %.4f kpc  to  %.2f kpc  (%d points) ===\n', ...
        r_min, r_max, N_pts);

% ==========================================================================
% BUILD COMPOSITE PROFILE
%     soliton_nfw_analytic handles the three-tier intersection finder
%     (polyxpoly -> intersectPolylines -> sign-change/fzero).
% ==========================================================================

comp = soliton_nfw_analytic(rho_c, rc_sol, Mh, Rvir, c_nfw, r);

fprintf('\n=== Composite profile ===\n');
fprintf('  Intersection method : %s\n',       comp.intersection_method);
fprintf('  r_x                 : %.4f kpc\n', comp.r_x);
fprintf('  rho_x               : %.4e M_sun/kpc^3\n', comp.rho_x);
fprintf('  r_x / rc_sol        : %.2f\n',     comp.r_x / rc_sol);
fprintf('  r_x / Rvir          : %.4f\n',     comp.r_x / Rvir);

% ==========================================================================
% 6.  SUMMARY TABLE
% ==========================================================================
fprintf('\n%s\n', repmat('-',1,55));
fprintf('%-28s  %14s  %s\n', 'Quantity', 'Value', 'Units');
fprintf('%s\n', repmat('-',1,55));
fprintf('%-28s  %14.3e  %s\n', 'rho_c (input)',    rho_c,   'M_sun/kpc^3');
fprintf('%-28s  %14.3f  %s\n', 'm22 (input)',      m22,     '10^-22 eV');
fprintf('%-28s  %14.4f  %s\n', 'rc_sol',           rc_sol,  'kpc');
fprintf('%-28s  %14.3e  %s\n', 'Mc_sol',           Mc_sol,  'M_sun');
fprintf('%-28s  %14.3e  %s\n', 'Mh',               Mh,      'M_sun');
fprintf('%-28s  %14.3f  %s\n', 'Rvir',             Rvir,    'kpc');
fprintf('%-28s  %14.3f  %s\n', 'c_nfw',            c_nfw,   '');
fprintf('%-28s  %14.4f  %s\n', 'rs_nfw',           comp.rs_nfw, 'kpc');
fprintf('%-28s  %14.3f  %s\n', 'V_c (soliton)',    V_c,     'km/s');
fprintf('%-28s  %14.3f  %s\n', 'Vvir (halo)',      Vvir,    'km/s');
fprintf('%-28s  %14.4f  %s\n', 'r_x',              comp.r_x,'kpc');
fprintf('%-28s  %14.3e  %s\n', 'rho_x',            comp.rho_x,'M_sun/kpc^3');
fprintf('%s\n', repmat('-',1,55));


% ==========================================================================
% Observables from density profile
% ==========================================================================
% From composite FDM profile directly
prof = profile_obs(comp.r, comp.rho_composite, 'all', hc.G);

figure;

% Left axis — surface density
yyaxis left
loglog(prof.r, prof.Sigma, 'b-', 'LineWidth', 2);
ylabel('\Sigma  [M_{sun}/kpc^2]');

% Right axis — velocity dispersions
yyaxis right
loglog(prof.r, prof.sigma_r,   '-',  'LineWidth', 2, ...
       'DisplayName', '\sigma_r');
hold on;
loglog(prof.r, prof.sigma_los, '--', 'LineWidth', 2, ...
       'DisplayName', '\sigma_{los}');
ylabel('\sigma  [km/s]');

xlabel('Radius [kpc]');
legend('Location', 'southwest');
grid on;
title('Profile observables');


%% analysis

% After computing prof.sigma_r on a log-log grid:
log_r     = log10(prof.r);
log_sig   = log10(prof.sigma_r);

% Local slope d(log sigma_r)/d(log r)
local_slope = diff(log_sig) ./ diff(log_r);
r_mid       = 10.^(0.5*(log_r(1:end-1) + log_r(2:end)));

% Find where slope changes most rapidly (inflection = max of |d slope/d log r|)
d_slope   = abs(diff(local_slope) ./ diff(log_r(1:end-1)));
r_mid2    = r_mid(1:end-1);

[~, idx]  = max(d_slope);
r_x_jeans = r_mid2(idx);

fprintf('r_x from sigma_r inflection: %.4f kpc\n', r_x_jeans);

%%

r_x_kin = find_rx_from_sigma(comp.r, prof.sigma_r, prof.sigma_los);

function r_x_kin = find_rx_from_sigma(r, sigma_r, sigma_los)
% FIND_RX_FROM_SIGMA
%   Estimate the soliton-NFW transition radius from the crossing of
%   sigma_r and sigma_los profiles.
%
%   Method: find the radius where |sigma_r - sigma_los| is minimised
%   after the first local peak of sigma_los (to avoid the noisy inner core).
%
%   INPUT
%     r         : radial grid [kpc]
%     sigma_r   : 3-D isotropic Jeans dispersion [km/s]
%     sigma_los : projected line-of-sight dispersion [km/s]
%
%   OUTPUT
%     r_x_kin   : estimated transition radius [kpc]

    r       = r(:);
    sigma_r = sigma_r(:);
    sigma_los = sigma_los(:);

    % Smooth both profiles lightly in log-space before differencing
    % to suppress particle noise (simulation data)
    log_r   = log10(r);
    win     = max(3, round(numel(r) * 0.04));   % ~4% of points, odd window
    if mod(win,2)==0, win = win+1; end

    sr_sm   = smooth_log(log_r, sigma_r,   win);
    sl_sm   = smooth_log(log_r, sigma_los, win);

    % Signed difference: positive where sigma_los > sigma_r (NFW-dominated)
    %                    negative where sigma_r  > sigma_los (core-dominated)
    diff_sig = sl_sm - sr_sm;

    % Find zero crossing(s) — sign change in diff_sig
    sign_ch  = find(diff(sign(diff_sig)) ~= 0);

    if isempty(sign_ch)
        % No clean crossing — fall back to minimum of |diff|
        [~, idx] = min(abs(diff_sig));
        r_x_kin  = r(idx);
        method   = 'min|delta_sigma|';
    else
        % Take the outermost crossing (innermost ones can be noise)
        i_c     = sign_ch(end);
        % Linear interpolation between bracketing points
        r_a     = r(i_c);       d_a = diff_sig(i_c);
        r_b     = r(i_c+1);     d_b = diff_sig(i_c+1);
        r_x_kin = r_a - d_a * (r_b - r_a) / (d_b - d_a);
        method  = 'zero-crossing';
    end

    fprintf('r_x (kinematic, %s): %.4f kpc\n', method, r_x_kin);
end

% ------------------------------------------------------------------
function y_sm = smooth_log(log_r, y, win)
% Gaussian-weighted smoother in log(r) space
    half = (win-1)/2;
    weights = exp(-0.5*((-half:half)/half*2).^2);
    weights = weights / sum(weights);
    y_sm = conv(y, weights, 'same');
end