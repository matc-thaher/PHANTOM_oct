% ==========================================================================
%  PHANTOM HMF Validation — dn/dlnM comparison
%  PHANTOM  vs  Colossus
%
%  Figure 1 : CDM models  (PS, ST, Tinker08, Crocce10, Watson13, Yung25,
%                          FernandezGarcia26)
%  Figure 2 : WDM/FDM models
% ==========================================================================
clear; clc; close all;

% --------------------------------------------------------------------------
% 0.  PATHS & SETTINGS
% --------------------------------------------------------------------------
addpath('F:\PHANTOM\src\utils',           '-begin');
addpath('F:\PHANTOM\src\concentration',   '-begin');
addpath('F:\PHANTOM\src\suppression',     '-begin');
addpath('F:\PHANTOM\src\profiles',        '-begin');
addpath('F:\PHANTOM\src\hmf',             '-begin');
addpath('F:\PHANTOM\src\halo',             '-begin');

python_exe   = 'F:\anaconda3\python.exe';
bridge_col   = 'F:\PHANTOM\src\utils\colossus_bridge.py';
cosmo_id     = 'planck18';
z            = 0.0;

% --------------------------------------------------------------------------
% 1.  PHANTOM COSMOLOGY
% --------------------------------------------------------------------------
cosmo                = cosmology(cosmo_id);
cosmo.transfer_model = 'eh98';
cosmo                = attach_linear_components(cosmo);

% --------------------------------------------------------------------------
% 2.  MASS GRID  [h^{-1} Msun]
% --------------------------------------------------------------------------
M = logspace(8, 14, 120)';

% --------------------------------------------------------------------------
% 3.  CDM PHANTOM HMF
% --------------------------------------------------------------------------
delta_c = collapse_overdensity('corrections', true, 'z', z, 'cosmo', cosmo);
Delta   = 200;
Delta_vir_c = halo_obs('M_h', 1e12, 'delta_vir', cosmo, z);
mdef    = '200m';

% models_phantom = {
%     'PS',               {delta_c};
%     'ST',               {delta_c};
%     'Tinker08',         {Delta};
%     'Crocce10',         {z};
%     'Watson13',         {z};
%     'Yung25',           {z};
%     'FernandezGarcia26',{M, z, mdef, cosmo};
%     'Fiorilli26', {z, mdef, cosmo};
%     };

models_phantom = {
    'PS',                 {delta_c};
    'ST',                 {delta_c};
    'Tinker08',           {Delta};
    'Watson13',           {z};
    'Crocce10',           {z};
    'Comparat17',         {z, cosmo, delta_c};
    'Yung25',             {z};
    'Fiorilli26',         {z, mdef, cosmo, false};
};



N_cdm = size(models_phantom, 1);

% labels_cdm = {'Press-Schechter'; 'Sheth-Tormen'; 'Tinker 2008'; ...
%               'Crocce 2010'; 'Watson 2013'; 'Yung25'; 'Fiorilli26'};

labels_cdm = {'Press-Schechter'; 'Sheth-Tormen'; ...
           'Tinker 2008'; 'Watson 2013'; ...
          'Crocce 2010'; 'Comparat 2017'; ...
          'Yung25'; 'Fiorilli26'};

fprintf('Computing PHANTOM CDM HMFs ...\n');
dndlnM_phantom = zeros(length(M), N_cdm);
for i = 1:N_cdm
    args = models_phantom{i, 2};
    [dndlnM_phantom(:,i), ~] = halo_mass_function(M, z, cosmo, ...
                                    models_phantom{i,1}, args{:});
    fprintf('  %s done.\n', labels_cdm{i});
end
%%
% --------------------------------------------------------------------------
% 4.  COLOSSUS CDM HMF
% --------------------------------------------------------------------------
models_col = {'press74'; 'sheth99';...
              'tinker08'; 'watson13'; 'crocce10'; 'comparat17'; ...
               'yung25'; 'fiorilli26'};


% mdefs_col = {'fof'; 'fof'; 'fof'; 'fof'; ...
%              '200m'; '200m'; 'fof'; 'fof'; ...
%              'fof'; 'fof'; '200c'; 'fof'; ...
%              'vir'; 'vir'; 'vir'};

fprintf('Querying Colossus CDM HMFs ...\n');
dndlnM_col = zeros(length(M), N_cdm);
for i = 1:N_cdm
    raw = colossus_query('hmf', M, cosmo_id, models_col{i}, z,...
                         python_exe, bridge_col);
    dndlnM_col(:,i) = raw;
    fprintf('  %s done.\n', labels_cdm{i});
end

res_col = dndlnM_phantom ./ dndlnM_col;
%%
% --------------------------------------------------------------------------
% 5.  WDM / FDM PHANTOM HMF
%     Example parameters — adjust m_wdm, M_hm, m22 as needed
% --------------------------------------------------------------------------
m_wdm  = 3.0;                    % WDM particle mass [keV]
m22    = 1.0;                    % FDM particle mass [10^-22 eV]

% Half-mode mass for WDM — compute from particle mass
% M_hm from Schneider+2012 / Lovell+2014 convention
% Placeholder: replace with your PHANTOM M_hm function if available
Omega_m = cosmo.Omega_m;
h       = cosmo.H0 / 100;
M_hm    = 1.0e10 * (m_wdm / 1.0)^(-3.33) * (Omega_m/0.25)^0.5 * (h/0.7)^2;  % approx

base_model = 'ST';   % CDM base for FDM suppression

models_wdm_fdm = {
    'Schneider12', {M, M_hm, z, cosmo, delta_c};
    'Lovell14',    {M, M_hm, z, cosmo, delta_c};
    'Schive16',    {M, M_hm, z, cosmo, delta_c, base_model};
    'Du17',        {M, z, m22, cosmo};
};
N_wfdm    = size(models_wdm_fdm, 1);
labels_wfdm = {'Schneider 2012 (WDM)'; 'Lovell 2014 (WDM)'; ...
               'Schive 2016 (FDM)'; 'Du 2017 (FDM)'};

fprintf('Computing PHANTOM WDM/FDM HMFs ...\n');
dndlnM_wfdm = zeros(length(M), N_wfdm);
for i = 1:N_wfdm
    args = models_wdm_fdm{i, 2};
    [dndlnM_wfdm(:,i), ~] = halo_mass_function(M, z, cosmo, ...
                                 models_wdm_fdm{i,1}, args{:});
    fprintf('  %s done.\n', labels_wfdm{i});
end

% CDM reference for WDM/FDM ratio (use ST)
[dndlnM_CDM_ref, ~] = halo_mass_function(M, z, cosmo, 'ST', delta_c);
res_wfdm = dndlnM_wfdm ./ dndlnM_CDM_ref;
%%
% ==========================================================================
% FIGURE 1 — CDM: dn/dlnM + PHANTOM/Colossus residual
% ==========================================================================
% fprintf('Plotting Figure 1 (CDM) ...\n');
% 
% cmap_cdm = turbo(N_cdm);
% lw   = 3;
% lw2  = 2.8;
% 
% fig1 = figure();
% 
% % --- Top panel ---
% ax1 = axes('Position',[0.12 0.30 0.86 0.67]);
% hold(ax1,'on'); box(ax1,'on');
% 
% for i = 1:N_cdm
%     loglog(ax1, M, dndlnM_phantom(:,i), '-',  ...
%            'Color', cmap_cdm(i,:), 'LineWidth', lw,  'DisplayName', labels_cdm{i});
%     loglog(ax1, M, dndlnM_col(:,i),     '--', ...
%            'Color', cmap_cdm(i,:), 'LineWidth', lw2, 'HandleVisibility','off');
% end
% loglog(ax1, NaN,NaN,'k-',  'LineWidth',lw,  'DisplayName','PHANTOM');
% loglog(ax1, NaN,NaN,'k--', 'LineWidth',lw2, 'DisplayName','Colossus');
% 
% ylabel(ax1, '$dn/d\ln M\ [h^3\,\mathrm{Mpc}^{-3}]$', ...
%        'Interpreter','latex','FontSize',32);
% set(ax1,'XTickLabel',[],'XScale','log','YScale','log', ...
%         'FontSize',22,'TickDir','in','TickLength',[0.015 0.015], 'LineWidth', 2);
% legend(ax1,'Location','best','NumColumns',2,'FontSize',22,'Interpreter','none');
% xlim(ax1,[min(M) max(M)]);
% 
% % --- Bottom panel — joined axes, no gap ---
% ax2 = axes('Position',[0.12 0.08 0.86 0.22]);
% hold(ax2,'on'); box(ax2,'on');
% 
% for i = 1:N_cdm
%     semilogx(ax2, M, res_col(:,i), '-', ...
%              'Color', cmap_cdm(i,:), 'LineWidth', lw);
% end
% % yline(ax2, 1, 'k--', 'LineWidth', 1.0);
% 
% xlabel(ax2, '$M\ [h^{-1}\,\mathrm{M}_\odot]$', ...
%        'Interpreter','latex','FontSize',32);
% ylabel(ax2, '$Ratio$', ...
%        'Interpreter','latex','FontSize',32);
% xlim(ax2,[min(M) max(M)]);
% ylim(ax2,[0.8 1.2]);
% set(ax2,'XScale','log','FontSize',22,'TickDir','in','TickLength',[0.015 0.015], 'LineWidth', 2);
% 
% % Link x-axes so zoom/pan stays in sync
% linkaxes([ax1 ax2],'x');
% 
% exportgraphics(fig1,'hmf_CDM_comparison.png','ContentType','vector');
% fprintf('Figure 1 saved.\n');

% ==========================================================================
% FIGURE 1 — CDM: dn/dlnM + PHANTOM/Colossus residual
% ==========================================================================
fprintf('Plotting Figure 1 (CDM) ...\n');
cmap_cdm = turbo(N_cdm);
lw   = 3;
lw2  = 2.8;

gap     = 0.0;   % vertical gap between the two panels
bottom2 = 0.15;    % bottom of lower panel  — leaves room for xlabel
h2      = 0.20;    % height of lower panel
bottom1 = bottom2 + h2 + gap;   % = 0.345
h1      = 0.97 - bottom1;       % = 0.625  (leaves 0.03 headroom at top)

fig1 = figure();

% --- Top panel ---
ax1 = axes('Position',[0.12 bottom1 0.86 h1]);
hold(ax1,'on'); box(ax1,'on');
for i = 1:N_cdm
    loglog(ax1, M, dndlnM_phantom(:,i), '-',  ...
           'Color', cmap_cdm(i,:), 'LineWidth', lw,  'DisplayName', labels_cdm{i});
    loglog(ax1, M, dndlnM_col(:,i),     '--', ...
           'Color', cmap_cdm(i,:), 'LineWidth', lw2, 'HandleVisibility','off');
end
loglog(ax1, NaN, NaN, 'k-',  'LineWidth', lw,  'DisplayName', 'PHANTOM');
loglog(ax1, NaN, NaN, 'k--', 'LineWidth', lw2, 'DisplayName', 'Colossus');
ylabel(ax1, '$dn/d\ln M\ [h^3\,\mathrm{Mpc}^{-3}]$', ...
       'Interpreter','latex','FontSize',32);
set(ax1, 'XTickLabel', [], 'XScale', 'log', 'YScale', 'log', ...
         'FontSize', 22, 'TickDir', 'in', 'TickLength', [0.015 0.015], ...
         'LineWidth', 2);
legend(ax1, 'Location', 'best', 'NumColumns', 2, 'FontSize', 22, ...
       'Interpreter', 'none');
xlim(ax1, [min(M) max(M)]);

% --- Bottom panel ---
ax2 = axes('Position',[0.12 bottom2 0.86 h2]);
hold(ax2,'on'); box(ax2,'on');
for i = 1:N_cdm
    semilogx(ax2, M, res_col(:,i), '-', ...
             'Color', cmap_cdm(i,:), 'LineWidth', lw);
end
xlabel(ax2, '$M\ [h^{-1}\,\mathrm{M}_\odot]$', ...
       'Interpreter', 'latex', 'FontSize', 32);
ylabel(ax2, '$Ratio$', ...
       'Interpreter', 'latex', 'FontSize', 32);
xlim(ax2, [min(M) max(M)]);
ylim(ax2, [0.8 1.2]);
set(ax2, 'XScale', 'log', 'FontSize', 22, 'TickDir', 'in', ...
         'TickLength', [0.015 0.015], 'LineWidth', 2);

linkaxes([ax1 ax2], 'x');
exportgraphics(fig1, 'hmf_CDM_comparison.png', 'ContentType', 'vector');
fprintf('Figure 1 saved.\n');

% ==========================================================================
% FIGURE 2 — WDM/FDM: dn/dlnM (top) + ratio to CDM (bottom)
% ==========================================================================
fprintf('Plotting Figure 2 (WDM/FDM) ...\n');

cmap_wfdm = lines(N_wfdm);
ls_wfdm   = {'-'; '--'; '-'; '--'};   % solid=WDM, dashed=FDM convention

fig2 = figure();

% --- Top panel ---
ax3 = axes; %('Position',[0.12 0.30 0.86 0.67]);
hold(ax3,'on'); box(ax3,'on');

% CDM reference line
loglog(ax3, M, dndlnM_CDM_ref, 'k-', 'LineWidth', lw, 'DisplayName','CDM (ST)');

for i = 1:N_wfdm
    loglog(ax3, M, dndlnM_wfdm(:,i), ls_wfdm{i}, ...
           'Color', cmap_wfdm(i,:), 'LineWidth', lw, 'DisplayName', labels_wfdm{i});
end

ylabel(ax3, '$dn/d\ln M\ [h^3\,\mathrm{Mpc}^{-3}]$', ...
       'Interpreter','latex','FontSize',32);
xlabel(ax3, '$M\ [h^{-1}\,\mathrm{M}_\odot]$', ...
       'Interpreter','latex','FontSize',32);
set(ax3,'XScale','log','YScale','log', ...
        'FontSize',22,'TickDir','in','TickLength',[0.015 0.015], 'LineWidth', 2);
legend(ax3,'Location','best','NumColumns',1,'FontSize',22,'Interpreter','none');
xlim(ax3,[min(M) max(M)]);
% title(ax3, sprintf('WDM: $m_{\\rm wdm}=%.1f$ keV,  FDM: $m_{22}=%.1f$', ...
%       m_wdm, m22), 'Interpreter','latex','FontSize',11);

% % --- Bottom panel ---
% ax4 = axes('Position',[0.12 0.08 0.86 0.22]);
% hold(ax4,'on'); box(ax4,'on');
% 
% for i = 1:N_wfdm
%     semilogx(ax4, M, res_wfdm(:,i), ls_wfdm{i}, ...
%              'Color', cmap_wfdm(i,:), 'LineWidth', lw, 'DisplayName', labels_wfdm{i});
% end
% yline(ax4, 1, 'k--', 'LineWidth', 1.0);
% 
% xlabel(ax4, '$M\ [h^{-1}\,\mathrm{M}_\odot]$', ...
%        'Interpreter','latex','FontSize',14);
% ylabel(ax4, '$n_{\rm model}/n_{\rm CDM}$', ...
%        'Interpreter','latex','FontSize',12);
% xlim(ax4,[min(M) max(M)]);
% ylim(ax4,[0.0 1.2]);
% set(ax4,'XScale','log','FontSize',11,'TickDir','in','TickLength',[0.015 0.015]);
% 
% linkaxes([ax3 ax4],'x');

exportgraphics(fig2,'hmf_WDMFDM_comparison.png','ContentType','vector');
fprintf('Figure 2 saved.\n');

fprintf('\nAll done.\n');