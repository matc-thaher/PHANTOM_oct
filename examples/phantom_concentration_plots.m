% ==========================================================================
%  PHANTOM Concentration Plots  —  refactored with loops + turbo colormap
%  Figure 1 : All CDM c(M) models from PHANTOM + Colossus comparison
%  Figure 2 : CDM Ishiyama21 + FDM & WDM suppressed concentrations
% ==========================================================================
clear; clc; close all;

% -----------------------------------------------------------------------
% 0.  PATHS & USER SETTINGS
% -----------------------------------------------------------------------
addpath('F:\PHANTOM\src\utils',         '-begin');
addpath('F:\PHANTOM\src\concentration', '-begin');
addpath('F:\PHANTOM\src\suppression',   '-begin');
addpath('F:\PHANTOM\src\profiles',      '-begin');

python_exe  = 'F:\anaconda3\python.exe';
bridge_path = 'F:\PHANTOM\src\utils\colossus_bridge.py';
cosmo_id    = 'planck18';
z           = 0.0;

m_ax_eV     = 1e-22;
m22         = m_ax_eV / 1e-22;
m_WDM_keV   = 3.0;

M_vec = logspace(8, 15, 80);

% -----------------------------------------------------------------------
% 1.  PHANTOM COSMOLOGY STRUCT
% -----------------------------------------------------------------------
cosmo                = cosmology(cosmo_id);
cosmo.transfer_model = 'eh98';
cosmo                = attach_linear_components(cosmo);

cosmo_k                = cosmology('planck13');
cosmo_k.transfer_model = 'eh98';
cosmo_k                = attach_linear_components(cosmo);

% -----------------------------------------------------------------------
% 2.  PHANTOM CDM CONCENTRATIONS  —  loop over model list
% -----------------------------------------------------------------------
% Each row: { model_name , needs_cosmo }
phantom_models = {
    'bullock01',  true;
    'dutton14',   true;
    'klypin16',   true;
    'ludlow16',   true;
    'ishiyama21', true;
};
N_cdm = size(phantom_models, 1);


fprintf('Computing PHANTOM concentrations ...\n');
% c_phantom = cell(N_cdm, 1);
% for i = 1:N_cdm
%     mname       = phantom_models{i,1};
%     c_phantom{i} = c_CDM(M_vec, z, mname, cosmo);
%     fprintf('  %s done\n', mname);
% end
% PHANTOM — individual, parameters matched to Colossus defaults
c_bullock01 = c_CDM(M_vec, z, 'bullock01', cosmo);
c_dutton14  = c_CDM(M_vec, z, 'dutton14',  cosmo, '200c');
% c_diemer15  = c_CDM(M_vec, z, 'diemer15',  cosmo, 'median');
c_klypin16  = c_CDM(M_vec, z, 'klypin16',  cosmo_k, '200c', 'cM');
c_ludlow16  = c_CDM(M_vec, z, 'ludlow16',  cosmo);
c_ishi21    = c_CDM(M_vec, z, 'ishiyama21',cosmo, '200c_all');

c_phantom = {c_bullock01; c_dutton14; c_klypin16; c_ludlow16; c_ishi21};
fprintf('PHANTOM concentrations done.\n');

% -----------------------------------------------------------------------
% 3.  COLOSSUS CONCENTRATIONS  —  same model list, loop
% -----------------------------------------------------------------------
% Colossus model names that match the PHANTOM list above
colossus_models = {
    'bullock01';
    'dutton14';
    'klypin16_m';
    'ludlow16';
    'ishiyama21';
};

fprintf('Querying Colossus concentrations ...\n');
c_colossus = cell(N_cdm, 1);
for i = 1:N_cdm
    if strcmp(colossus_models{i}, 'klypin16_m')
        cosmo_id_use = 'planck13';
    else
        cosmo_id_use = cosmo_id;
    end
    c_colossus{i} = colossus_query('concentration', M_vec, cosmo_id_use, ...
                                    colossus_models{i}, z, python_exe, bridge_path);
    fprintf('  %s done\n', colossus_models{i});
end
fprintf('Colossus concentrations done.\n');

% -----------------------------------------------------------------------
% 4.  SUPPRESSED CONCENTRATIONS  (FDM / WDM)
% -----------------------------------------------------------------------
M_half_fdm = halfmode_mass('fdm', m22,       cosmo);
M_half_wdm = halfmode_mass('WDM', m_WDM_keV, cosmo);

% Grab Ishiyama21 — last entry in the loop above
c_ishi21 = c_phantom{end};

% Suppression models: { label , sup_tag , multiplier_on_ishi21 , linestyle }
sup_models = {
    '$CDM_{Ishiyama21}$',                                            'none',       []        , '-' ;
    sprintf('$FDM_{Laroche22} (m_{22}=1$)'),     'laroche22',  M_half_fdm, '-' ;
    sprintf('$FDM_{Dentler22} (m_{22}=1$)'),     'dentler22',  M_half_fdm, '--';
    sprintf('$WDM_{Schneider12} (%.1f keV)$', m_WDM_keV), 'schneider12', M_half_wdm, '-' ;
    };
N_sup = size(sup_models, 1);

sup_factor = cell(N_sup, 1);
c_sup      = cell(N_sup, 1);
for i = 1:N_sup
    tag = sup_models{i,2};
    if strcmp(tag,'none')
        sup_factor{i} = ones(size(M_vec));
        c_sup{i}      = c_ishi21;
    else
        sup_factor{i} = suppression_factor(M_vec, sup_models{i,3}, tag);
        c_sup{i}      = c_ishi21 .* sup_factor{i};
    end
end

%%
% -----------------------------------------------------------------------
% 5.  FIGURE 1 :  CDM models  +  Colossus  —  turbo colormap
% -----------------------------------------------------------------------
fprintf('Plotting Figure 1 ...\n');

cmap1 = turbo(N_cdm);   % one colour per model, consistent across both panels
lw    = 2;
lwc   = 2;

fig1 = figure;

% ---- top panel ----
ax1 = subplot(3,1,[1 2]);
hold(ax1,'on'); box(ax1,'on');

for i = 1:N_cdm
    mname = phantom_models{i,1};
    % PHANTOM  —  solid
    semilogy(ax1, log10(M_vec), c_phantom{i},  '-',  ...
             'Color', cmap1(i,:), 'LineWidth', lw, ...
             'DisplayName', mname);
    % Colossus  —  dashed, same colour
    semilogy(ax1, log10(M_vec), c_colossus{i}, '--', ...
             'Color', cmap1(i,:), 'LineWidth', lwc, ...
             'DisplayName', colossus_models{i});
end

ylabel(ax1,'$c$','Interpreter','latex','FontSize',32);
% title(ax1, sprintf('CDM Concentration--Mass Relations ($z = %.1f$, Planck18)', z), ...
%       'Interpreter','latex','FontSize',13);
xlim(ax1,[8 15]);
set(ax1,'XTickLabel',[],'XScale','linear','YScale','log','FontSize',22);
% grid(ax1,'on');


lgd1  = legend(ax1, 'Location','best',...
               'NumColumns', 2, 'FontSize',15,'Interpreter','none');
% lgd1.Box = 'off';

% ---- bottom panel : residuals ----
ax2 = subplot(3,1,3);
hold(ax2,'on'); box(ax2,'on');
% yline(ax2, 1, 'k--', 'LineWidth', 1.2);

for i = 1:N_cdm
    ph  = c_phantom{i}(:);
    col = c_colossus{i}(:);
    res = ph ./col ; %100 * (ph ./ col - 1); % 
    plot(ax2, log10(M_vec), res, '-', ...
         'Color', cmap1(i,:), 'LineWidth', lw, ...
         'DisplayName', phantom_models{i,1});
end

xlabel(ax2,'$\log_{10}(M\ [\mathrm{M}_\odot\,h^{-1}])$','Interpreter','latex','FontSize',32);
ylabel(ax2,'$\Delta c/c_{\rm Col}$','Interpreter','latex','FontSize',32);
xlim(ax2,[8 15]);
ylim(ax2, [0.97 1.03])
set(ax2,'FontSize',22);
% grid(ax2,'on');


% % Remove bottom x-axis tick marks and labels from top panel
% set(ax1, 'XTickLabel', [], 'XTick', 8:1:15);

% Tight gap between panels — set positions so they share exactly one line
left   = 0.10;
right  = 0.88;   % width
bot2   = 0.15;   % bottom of lower panel
h2     = 0.20;   % height of lower panel
gap    = 0.00;   % zero gap
h1     = 0.60;   % height of upper panel
bot1   = bot2 + h2 + gap;

set(ax1, 'Position', [left, bot1, right, h1], 'LineWidth',2, 'FontSize', 20);
set(ax2, 'Position', [left, bot2, right, h2], 'LineWidth',2, 'FontSize', 20);

% Remove the top panel's bottom x-axis line and the lower panel's top line
% so they share one clean border
ax1.XAxis.Visible = 'on';
set(ax1, 'XTickLabel', [], 'TickDir', 'in');
set(ax2, 'TickDir', 'in');

% Link x-axes so zoom/pan stays synchronised
linkaxes([ax1, ax2], 'x');

exportgraphics(fig1,'concentration_CDM_comparison.png','ContentType','vector');
fprintf('Figure 1 saved.\n');

% -----------------------------------------------------------------------
% 6.  FIGURE 2 :  CDM + FDM + WDM suppressed concentrations
%                 Single panel — concentration lines only
% -----------------------------------------------------------------------
fprintf('Plotting Figure 2 ...\n');

cmap2 = turbo(N_sup);
lw2   = 2.0;

fig2 = figure;

hold on; box on;

xline(log10(M_half_fdm), ':', 'Color', cmap2(2,:), 'LineWidth',1.2, ...
      'Label','$M_{1/2}^{\rm FDM}$','Interpreter','latex',...
      'LabelVerticalAlignment', 'bottom','FontSize',12, 'HandleVisibility', 'off');
xline(log10(M_half_wdm), ':', 'Color', cmap2(4,:), 'LineWidth',1.2, ...
      'Label','$M_{1/2}^{\rm WDM}$','Interpreter','latex',...
      'LabelVerticalAlignment', 'bottom','FontSize',12, 'HandleVisibility', 'off');

for i = 1:N_sup
    semilogy(log10(M_vec), c_sup{i}, sup_models{i,4}, ...
             'Color', cmap2(i,:), 'LineWidth', lw2, ...
             'DisplayName', sup_models{i,1});
end

xlabel('$\log_{10}(M\ [\mathrm{M}_\odot\,h^{-1}])$','Interpreter','latex','FontSize',32);
ylabel('$c$','Interpreter','latex','FontSize',32);
% title(sprintf('Suppressed Concentration--Mass Relations ($z = %.1f$)', z), ...
      % 'Interpreter','latex','FontSize',13);
xlim([8 15]);
set(gca,'XScale','linear','YScale','log','FontSize',22, 'LineWidth', 2);
lgd2  = legend('Location','best',...
               'NumColumns', 1, 'FontSize',20,'Interpreter','latex');
% grid on;

% colormap(turbo);
% cb                      = colorbar('Location','eastoutside');
% clim([0.5, N_sup+0.5]);
% cb.Ticks                = 1:N_sup;
% cb.TickLabels           = sup_models(:,1);
% cb.Label.String         = 'Dark matter model';
% cb.Label.FontSize       = 10;
% cb.FontSize             = 9;
% cb.TickLabelInterpreter = 'latex';

exportgraphics(fig2,'concentration_suppression.png','ContentType','vector');
fprintf('Figure 2 saved.\n');
fprintf('\nAll done.\n');
