% =========================================================
% Science use case: weak-lensing convergence profile
% PHANTOM — NFW, Einasto, FDM composite
% Reference: Wright & Brainerd (2000, ApJ 534, 34)
%            Bartelmann (1996, A&A 313, 697)
% =========================================================
clear; clc; close all;
% --------------------------------------------------------------------------
% 0.  PATHS & SETTINGS
% --------------------------------------------------------------------------
addpath('F:\PHANTOM\src\utils',           '-begin');
addpath('F:\PHANTOM\src\concentration',   '-begin');
addpath('F:\PHANTOM\src\profiles',        '-begin');
addpath('F:\PHANTOM\src\halo',             '-begin');
addpath('F:\PHANTOM\src\fdm',             '-begin');


% =========================================================
% Science use case: weak-lensing convergence profile
% =========================================================
cosmo  = cosmology('Planck18');

% --- Halo and lens parameters ---
Mh    = 1e14;       % Msun/h  — cluster-mass halo
c     = 5.0;        % concentration
z_l   = 0.3;        % lens redshift
z_s   = 1.0;        % source redshift
Delta = 200;

% % --- Radial grid in Mpc/h (for NFW, Hernquist, Einasto) ---
% r_Mpc = logspace(-2, 0.5, 500)';     % Mpc/h

% --- Radial grid in kpc (for Soliton and composite) ---
% r_kpc = r_Mpc * 1e3 / cosmo.h;      % kpc
r_kpc = logspace(-2, 3, 500)';           % 0.01 kpc to 1000 kpc

% --- Radial grid in Mpc/h (for NFW, Hernquist, Einasto) ---
r_Mpc = r_kpc * cosmo.h / 1e3;          % convert to Mpc/h for profile functions

% =========================================================
% 1. NFW profile
% =========================================================
rho_c    = cosmo.rho_crit0 * cosmo.E(z_l)^2;
R_Delta  = (3*Mh / (4*pi * Delta * rho_c))^(1/3);   % Mpc/h
rs_nfw   = R_Delta / c;
fc_nfw   = log(1+c) - c/(1+c);
rhos_nfw = Mh / (4*pi * rs_nfw^3 * fc_nfw);

rho_nfw  = NFW_profile(r_Mpc, rhos_nfw, rs_nfw);    % Msun/h / (Mpc/h)^3

% profile_obs expects kpc and Msun/kpc^3 — convert
% 1 Msun/h / (Mpc/h)^3 = (1/h) / (1e3/h kpc)^3 * Msun
%                       = h^2 / 1e9  Msun/kpc^3
unit_conv = cosmo.h^2 / 1e9;
prof_nfw  = profile_obs(r_kpc, rho_nfw * unit_conv, 'Sigma');
Sigma_nfw = prof_nfw.Sigma;    % Msun/kpc^2

% =========================================================
% 2. Hernquist profile
% =========================================================
[rho_hern, ~, ~, ~] = Hernquist_profile(r_Mpc, Mh, c, z_l, cosmo, Delta);
prof_hern  = profile_obs(r_kpc, rho_hern * unit_conv, 'Sigma');
Sigma_hern = prof_hern.Sigma;

% =========================================================
% 3. Einasto profile
% =========================================================
[rho_ein, ~, ~, ~] = Einasto_profile(r_Mpc, Mh, c, z_l, cosmo, Delta);
prof_ein  = profile_obs(r_kpc, rho_ein * unit_conv, 'Sigma');
Sigma_ein = prof_ein.Sigma;

% =========================================================
% 4. FDM composite using soliton_obs dispatcher
% =========================================================
m22      = 1.0;                          % boson mass in units of 1e-22 eV
Mvir_Msun = Mh / cosmo.h;               % Msun (no h factor)
Rvir_kpc  = R_Delta * 1e3 / cosmo.h;    % kpc

% Get core radius from halo mass via Schive+2014 CHMR
rc_kpc  = soliton_obs('M_h', Mvir_Msun, 'r_c', m22);       % kpc

% Get central density from core radius
rho0_sol = soliton_obs('r_c', rc_kpc, 'rho_c', m22);        % Msun/kpc^3

% Build composite profile and get surface density
comp     = soliton_nfw_analytic(rho0_sol, rc_kpc, Mvir_Msun, Rvir_kpc, c, r_kpc, false);
prof_fdm = profile_obs(r_kpc, comp.rho_composite, 'Sigma');
Sigma_fdm = prof_fdm.Sigma;              % Msun/kpc^2

% =========================================================
% Critical surface density
% =========================================================
c_km_s  = 299792.458;
G_kpc   = 4.3009e-6;           % kpc (km/s)^2 / Msun

d_C_l   = cosmo.comovingDistance(z_l);     % Mpc/h
d_C_s   = cosmo.comovingDistance(z_s);     % Mpc/h
D_l     = cosmo.angularDiameterDistance(z_l);
D_s     = cosmo.angularDiameterDistance(z_s);
D_ls    = (d_C_s - d_C_l) / (1 + z_s);    % Mpc/h, flat cosmology

% Convert distances Mpc/h -> kpc
D_l_kpc  = D_l  * 1e3 / cosmo.h;
D_s_kpc  = D_s  * 1e3 / cosmo.h;
D_ls_kpc = D_ls * 1e3 / cosmo.h;

% Sigma_cr in Msun/kpc^2
Sigma_cr = (c_km_s^2 / (4*pi*G_kpc)) * (D_s_kpc / (D_l_kpc * D_ls_kpc));

% =========================================================
% Convergence kappa = Sigma / Sigma_cr
% =========================================================
kappa_nfw  = Sigma_nfw  / Sigma_cr;
kappa_hern = Sigma_hern / Sigma_cr;
kappa_ein  = Sigma_ein  / Sigma_cr;
kappa_fdm  = Sigma_fdm  / Sigma_cr;

%%
% =========================================================
% Plot
% =========================================================
% gap     = 0.0;
% bottom2 = 0.13;
% h2      = 0.20;
% bottom1 = bottom2 + h2 + gap;
% h1      = 0.97 - bottom1;
% 
% fig = figure();
% 
% ax1 = axes('Position',[0.13 bottom1 0.85 h1]);
% hold(ax1,'on'); box(ax1,'on');
% loglog(ax1, r_kpc, kappa_nfw,  'b-',  'LineWidth',2.5, 'DisplayName','NFW');
% loglog(ax1, r_kpc, kappa_hern, 'r--', 'LineWidth',2.5, 'DisplayName','Hernquist');
% loglog(ax1, r_kpc, kappa_ein,  'g-.', 'LineWidth',2.5, 'DisplayName','Einasto');
% loglog(ax1, r_kpc, kappa_fdm,  'k:',  'LineWidth',2.5, 'DisplayName','FDM (soliton+NFW)');
% ylabel(ax1, '$\kappa(R)$', 'Interpreter','latex','FontSize',28);
% set(ax1,'XTickLabel',[],'XScale','log','YScale','log', ...
%         'FontSize',18,'TickDir','in','TickLength',[0.015 0.015],'LineWidth',1.8);
% legend(ax1,'Location','southwest','FontSize',16);
% xlim(ax1,[min(r_kpc) max(r_kpc)]);
% ylim(ax1, [min(kappa_hern) max(kappa_fdm)])
% 
% ax2 = axes('Position',[0.13 bottom2 0.85 h2]);
% hold(ax2,'on'); box(ax2,'on');
% semilogx(ax2, r_kpc, kappa_hern./kappa_nfw, 'r--', 'LineWidth',2, 'DisplayName','Hernquist/NFW');
% semilogx(ax2, r_kpc, kappa_ein ./kappa_nfw, 'g-.', 'LineWidth',2, 'DisplayName','Einasto/NFW');
% semilogx(ax2, r_kpc, kappa_fdm ./kappa_nfw, 'k:',  'LineWidth',2, 'DisplayName','FDM/NFW');
% xlabel(ax2, '$R\ [\mathrm{kpc}]$','Interpreter','latex','FontSize',28);
% ylabel(ax2, '$\kappa/\kappa_\mathrm{NFW}$','Interpreter','latex','FontSize',24);
% xlim(ax2,[min(r_kpc) max(r_kpc)]);
% ylim(ax2,[-1.0 max(kappa_fdm ./kappa_nfw)]);
% set(ax2,'XScale','log','FontSize',18,'TickDir','in', ...
%         'TickLength',[0.015 0.015],'LineWidth',1.8);
% legend(ax2,'Location','northeast','FontSize',14);
% 
% linkaxes([ax1 ax2],'x');
% exportgraphics(fig,'lensing_convergence.png','ContentType','vector');

% =========================================================
% Figure 1 — CDM profiles
% =========================================================
gap  = 0.0;
bot2 = 0.13;
h2   = 0.22;
bot1 = bot2 + h2 + gap;
h1   = 0.95 - bot1;

fig1 = figure();

axL1 = axes('Position',[0.14 bot1 0.83 h1]);
hold(axL1,'on'); box(axL1,'on');
loglog(axL1, r_kpc, kappa_nfw,  'b-',  'LineWidth',2.5, 'DisplayName','NFW');
loglog(axL1, r_kpc, kappa_hern, 'r--', 'LineWidth',2.5, 'DisplayName','Hernquist');
loglog(axL1, r_kpc, kappa_ein,  'g-.', 'LineWidth',2.5, 'DisplayName','Einasto');
ylabel(axL1,'$\kappa(R)$','Interpreter','latex','FontSize',32);
set(axL1,'XTickLabel',[],'XScale','log','YScale','log', ...
    'FontSize',22,'TickDir','in','LineWidth',1.8);
legend(axL1,'Location','best','FontSize',20);
xlim(axL1,[min(r_kpc) max(r_kpc)]);
ylim(axL1,[1e-5 max(kappa_nfw)*3]);

axL2 = axes('Position',[0.14 bot2 0.83 h2]);
hold(axL2,'on'); box(axL2,'on');
semilogx(axL2, r_kpc, kappa_hern./kappa_nfw, 'r--', 'LineWidth',2, 'DisplayName','Hernquist/NFW');
semilogx(axL2, r_kpc, kappa_ein ./kappa_nfw, 'g-.', 'LineWidth',2, 'DisplayName','Einasto/NFW');
% yline(axL2, 1, 'k-', 'LineWidth',1.2);
xlabel(axL2,'$R\ [\mathrm{kpc}]$','Interpreter','latex','FontSize',32);
ylabel(axL2,'$\kappa/\kappa_\mathrm{NFW}$','Interpreter','latex','FontSize',32);
xlim(axL2,[min(r_kpc) max(r_kpc)]);
ylim(axL2,[0.0 2.8]);
set(axL2,'XScale','log','YScale','linear','FontSize',22,'TickDir','in', ...
    'LineWidth',1.8);
legend(axL2,'Location','best','FontSize',20);

linkaxes([axL1 axL2],'x');
exportgraphics(fig1,'lensing_CDM.png','ContentType','vector');

% =========================================================
% Figure 2 — FDM vs NFW
% =========================================================
fig2 = figure();

axR1 = axes('Position',[0.14 bot1 0.83 h1]);
hold(axR1,'on'); box(axR1,'on');
loglog(axR1, r_kpc,     kappa_nfw, 'b-', 'LineWidth',2.5, 'DisplayName','NFW');
loglog(axR1, r_kpc, kappa_fdm, 'k:', 'LineWidth',2.5, 'DisplayName','FDM (soliton+NFW)');
xline(axR1, rc_kpc, '--', 'Color',[0.5 0.5 0.5], 'LineWidth',1.5, ...
    'Label','$r_c$','Interpreter','latex','FontSize',20, ...
    'LabelVerticalAlignment','top', 'HandleVisibility', 'off');
ylabel(axR1,'$\kappa(R)$','Interpreter','latex','FontSize',32);
set(axR1,'XTickLabel',[],'XScale','log','YScale','log', ...
    'FontSize',22,'TickDir','in','LineWidth',1.8);
legend(axR1,'Location','best','FontSize',20);
xlim(axR1,[min(r_kpc) max(r_kpc)]);
ylim(axR1,[1e-5 max(kappa_fdm)*3]);

axR2 = axes('Position',[0.14 bot2 0.83 h2]);
hold(axR2,'on'); box(axR2,'on');
kappa_nfw_interp = interp1(r_kpc, kappa_nfw, r_kpc, 'linear','extrap');
semilogx(axR2, r_kpc, kappa_fdm./kappa_nfw_interp, 'k:', 'LineWidth',2, ...
    'DisplayName','FDM/NFW');
% yline(axR2, 1, 'b-', 'LineWidth',1.2);
xline(axR2, rc_kpc, '--', 'Color',[0.5 0.5 0.5], 'LineWidth',1.5, 'HandleVisibility','off');
xlabel(axR2,'$R\ [\mathrm{kpc}]$','Interpreter','latex','FontSize',32);
ylabel(axR2,'$\kappa/\kappa_\mathrm{NFW}$','Interpreter','latex','FontSize',32);
xlim(axR2,[min(r_kpc) max(r_kpc)]);
ylim(axR2,[0.5 max(kappa_fdm./kappa_nfw_interp)*2]);
set(axR2,'XScale','log','YScale','log','FontSize',20,'TickDir','in', ...
    'LineWidth',1.8);
legend(axR2,'Location','best','FontSize',22);

linkaxes([axR1 axR2],'x');
exportgraphics(fig2,'lensing_FDM.png','ContentType','vector');