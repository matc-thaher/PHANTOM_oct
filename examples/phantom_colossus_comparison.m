%% PHANTOM vs Colossus — live comparison
clear; close all; clc;

addpath('F:\PHANTOM\src\utils', '-begin');
addpath('F:\PHANTOM\src\concentration', '-begin');
addpath('F:\PHANTOM\src\suppression', '-begin');
addpath('F:\PHANTOM\src\profiles', '-begin');

python_exe  = "F:\anaconda3\python.exe";
bridge_path = "F:\PHANTOM\src\utils\colossus_bridge.py";
cosmo_id    = 'planck15';
z           = 0.0;

%% Initialize PHANTOM
cosmo = cosmology(cosmo_id);
cosmo = derive_cosmo_params(cosmo);

%% --- Power spectrum ---
k = logspace(-3, 2, 300)';

% PHANTOM
cosmo_eh     = cosmo; cosmo_eh.transfer_model = 'eh98';
cosmo_eh     = attach_linear_components(cosmo_eh);
cosmo_ehfull = cosmo; cosmo_ehfull.transfer_model = 'eh98_full';
cosmo_ehfull = attach_linear_components(cosmo_ehfull);

Pk_ph_eh     = cosmo_eh.Pk(k, z);
Pk_ph_ehfull = cosmo_ehfull.Pk(k, z);

% Colossus — called live
Pk_co_eh     = colossus_query('power_spectrum', k, cosmo_id, 'eisenstein98',    z, python_exe, bridge_path);
Pk_co_ehfull = colossus_query('power_spectrum', k, cosmo_id, 'eisenstein98',    z, python_exe, bridge_path);
% Note: Colossus 'eisenstein98' is the full model; 'eisenstein98zb' is zero-baryon
Pk_co_ehzb   = colossus_query('power_spectrum', k, cosmo_id, 'eisenstein98_zb',  z, python_exe, bridge_path);

%% --- Variance ---
R = logspace(-2, 2, 200)';

sigma_ph_eh     = cosmo_eh.sigmaR(R, z);
sigma_ph_ehfull = cosmo_ehfull.sigmaR(R, z);

sigma_co_eh     = colossus_query('variance', R, cosmo_id, 'eisenstein98',   z, python_exe, bridge_path);
sigma_co_ehzb   = colossus_query('variance', R, cosmo_id, 'eisenstein98_zb', z, python_exe, bridge_path);

%% --- Correlation function ---
r = logspace(-1, 2.5, 150)';

xi_ph_eh     = cosmo_eh.correlationFunction(r, z);
xi_ph_ehfull = cosmo_ehfull.correlationFunction(r, z);

xi_co_eh     = colossus_query('correlation_function', r, cosmo_id, 'eisenstein98',   z, python_exe, bridge_path);
xi_co_ehzb   = colossus_query('correlation_function', r, cosmo_id, 'eisenstein98_zb', z, python_exe, bridge_path);

%% --- Ratios ---
ratio_Pk_eh_zb    = Pk_ph_eh     ./ Pk_co_ehzb;
ratio_pk_eh_full  = Pk_ph_ehfull ./ Pk_co_ehfull;
ratio_var_eh_zb   = sigma_ph_eh' ./ sigma_co_ehzb;
ratio_var_eh_full = sigma_ph_ehfull' ./ sigma_co_eh;
ratio_xi_eh_zb    = xi_ph_eh     ./ xi_co_ehzb;
ratio_xi_eh_full  = xi_ph_ehfull ./ xi_co_eh;

%% --- Plot (same layout as before) ---
%% Plot 1: Power Spectrum Comparison
figure; hold on; box on;
plot(k, Pk_co_ehzb, '-', 'LineWidth', 3, 'Color', [0.8, 0.2, 0.2], 'DisplayName', 'Colossus EH-zb');
plot(k, Pk_ph_eh, '--', 'LineWidth', 3, 'Color', [0.8, 0.2, 0.2], 'DisplayName', 'PHANTOM EH-zb');
plot(k, Pk_co_ehfull, '-', 'LineWidth', 3, 'Color', [0.2, 0.6, 0.8], 'DisplayName', 'Colossus EH-full');
plot(k, Pk_ph_ehfull, '--', 'LineWidth', 3, 'Color', [0.2, 0.6, 0.8], 'DisplayName', 'PHANTOM EH-full');
% loglog(k, pk_colossus.P_CAMB, '-', 'LineWidth', 2, 'Color', [0.2, 0.8, 0.2], 'DisplayName', 'Colossus CAMB');
% loglog(k, P_camb_phantom, '--', 'LineWidth', 2, 'Color', [0.2, 0.8, 0.2], 'DisplayName', 'PHANTOM CAMB');
xlabel('k [h/Mpc]', 'FontSize', 32);
ylabel('P(k) [(Mpc/h)^3]', 'FontSize', 32);
% title('Power Spectrum: PHANTOM vs Colossus', 'FontSize', 14);
legend('Location', 'best', 'FontSize', 22);
set(gca, 'XScale', 'log', 'YScale', 'log', 'FontSize', 22, 'LineWidth', 2);
ax=gca; yl=ax.YLim; yl(2) = yl(2)*3;
ylim(yl)
% axis tight;
saveas(gcf, 'power_spectrum_comparison.png');

figure; hold on; box on;
plot(k, ratio_Pk_eh_zb, '-', 'LineWidth', 3, 'Color', [0.8, 0.2, 0.2], 'DisplayName', 'EH-zb');
plot(k, ratio_pk_eh_full, '-', 'LineWidth', 3, 'Color', [0.2, 0.6, 0.8], 'DisplayName', 'EH-full');
% semilogx(k, ratio_pk_camb, '-', 'LineWidth', 2, 'Color', [0.2, 0.8, 0.2], 'DisplayName', 'CAMB');
yline(1, 'k--', 'LineWidth', 2);
xlabel('k [h/Mpc]', 'FontSize', 32);
ylabel('Ratio: PHANTOM / Colossus', 'FontSize', 32);
% title('Power Spectrum Ratio', 'FontSize', 14);
legend('Location', 'best', 'FontSize', 22);
ylim([0.95, 1.05]);
set(gca, 'XScale', 'log', 'YScale', 'linear', 'FontSize', 22, 'LineWidth', 2);
axis tight;

saveas(gcf, 'power_spectrum_residual.png');

%% Plot 2: Variance Comparison

figure; hold on; box on;
loglog(R, sigma_co_ehzb, '-', 'LineWidth', 3, 'Color', [0.8, 0.2, 0.2], 'DisplayName', 'Colossus EH-zb');
loglog(R, sigma_ph_eh, '--', 'LineWidth', 3, 'Color', [0.8, 0.2, 0.2], 'DisplayName', 'PHANTOM EH-zb');
loglog(R, sigma_co_eh, '-', 'LineWidth', 3, 'Color', [0.2, 0.6, 0.8], 'DisplayName', 'Colossus EH-full');
loglog(R, sigma_ph_ehfull, '--', 'LineWidth', 3, 'Color', [0.2, 0.6, 0.8], 'DisplayName', 'PHANTOM EH-full');
% loglog(R, var_colossus.sigma_CAMB, '-', 'LineWidth', 2, 'Color', [0.2, 0.8, 0.2], 'DisplayName', 'Colossus CAMB');
% loglog(R, sigma_camb_phantom, '--', 'LineWidth', 2, 'Color', [0.2, 0.8, 0.2], 'DisplayName', 'PHANTOM CAMB');
xlabel('R [Mpc/h]', 'FontSize', 32);
ylabel('\sigma(R)', 'FontSize', 32);
% title('Variance: PHANTOM vs Colossus', 'FontSize', 14);
legend('Location', 'best', 'FontSize', 22);
set(gca, 'XScale', 'log', 'YScale', 'log', 'FontSize', 22, 'LineWidth', 2);
axis tight;
saveas(gcf, 'variance_comparison.png');

figure; hold on; box on;
plot(R, ratio_var_eh_zb, '-', 'LineWidth', 3, 'Color', [0.8, 0.2, 0.2], 'DisplayName', 'EH-zb');
plot(R, ratio_var_eh_full, '-', 'LineWidth', 3, 'Color', [0.2, 0.6, 0.8], 'DisplayName', 'EH-full');
% semilogx(R, ratio_var_camb, '-', 'LineWidth', 2, 'Color', [0.2, 0.8, 0.2], 'DisplayName', 'CAMB');
yline(1, 'k--', 'LineWidth', 2);
xlabel('R [Mpc/h]', 'FontSize', 32);
ylabel('Ratio: PHANTOM / Colossus', 'FontSize', 32);
% title('Variance Ratio', 'FontSize', 14);
legend('Location', 'best', 'FontSize', 22);
ylim([0.95, 1.05]);
set(gca, 'XScale', 'log', 'YScale', 'linear', 'FontSize', 22, 'LineWidth', 2);
axis tight;

saveas(gcf, 'variance_residual.png');

%% Plot 3: Correlation Function Comparison

figure; hold on; box on;
plot(r, abs(xi_co_ehzb), '-', 'LineWidth', 3, 'Color', [0.8, 0.2, 0.2], 'DisplayName', 'Colossus EH-zb');
plot(r, abs(xi_ph_eh), '--', 'LineWidth', 3, 'Color', [0.8, 0.2, 0.2], 'DisplayName', 'PHANTOM EH-zb');
plot(r, abs(xi_co_eh), '-', 'LineWidth', 3, 'Color', [0.2, 0.6, 0.8], 'DisplayName', 'Colossus EH-full');
plot(r, abs(xi_ph_ehfull), '--', 'LineWidth', 3, 'Color', [0.2, 0.6, 0.8], 'DisplayName', 'PHANTOM EH-full');
% loglog(r, abs(xi_colossus.xi_CAMB), '-', 'LineWidth', 2, 'Color', [0.2, 0.8, 0.2], 'DisplayName', 'Colossus CAMB');
% loglog(r, abs(xi_camb_phantom), '--', 'LineWidth', 2, 'Color', [0.2, 0.8, 0.2], 'DisplayName', 'PHANTOM CAMB');
xlabel('r [Mpc/h]', 'FontSize', 32);
ylabel('|\xi(r)|', 'FontSize', 32);
%title('Correlation Function: PHANTOM vs Colossus', 'FontSize', 14);
legend('Location', 'best', 'FontSize', 22);
set(gca, 'XScale', 'log', 'YScale', 'log', 'FontSize', 22, 'LineWidth', 2);
axis tight;
saveas(gcf, 'correlation_function_comparison.png');

figure; hold on; box on;
plot(r, ratio_xi_eh_zb, '-', 'LineWidth', 3, 'Color', [0.8, 0.2, 0.2], 'DisplayName', 'EH-zb');
plot(r, ratio_xi_eh_full, '-', 'LineWidth', 3, 'Color', [0.2, 0.6, 0.8], 'DisplayName', 'EH-full');
% semilogx(r, ratio_xi_camb, '-', 'LineWidth', 2, 'Color', [0.2, 0.8, 0.2], 'DisplayName', 'CAMB');
yline(1, 'k--', 'LineWidth', 2);
xlabel('r [Mpc/h]', 'FontSize', 32);
ylabel('Ratio: PHANTOM / Colossus', 'FontSize', 32);
% title('Correlation Function Ratio', 'FontSize', 14);
legend('Location', 'best', 'FontSize', 22);
% ylim([0.95, 1.05]);
set(gca, 'XScale', 'log', 'YScale', 'linear', 'FontSize', 22, 'LineWidth', 2);
axis tight;

saveas(gcf, 'correlation_function_residual.png');


%% Plot 1c: Power Spectrum — WDM and FDM suppression

% Assume you have transfer_model keys for WDM / FDM
cosmo_wdm = cosmo;
cosmo_wdm.transfer_model = 'viel05';         % <-- use your actual key
cosmo_wdm = attach_linear_components(cosmo_wdm);

cosmo_fdm = cosmo;
cosmo_fdm.transfer_model = 'schive25';    % <-- or your FDM key
cosmo_fdm = attach_linear_components(cosmo_fdm);

Pk_cdm = Pk_ph_ehfull;                    % CDM reference (EH98 full)
Pk_wdm = cosmo_wdm.Pk(k, z);
Pk_fdm = cosmo_fdm.Pk(k, z);

figure; hold on; box on;
plot(k, Pk_cdm,  '-',  'LineWidth', 3, 'Color', [0.1, 0.4, 0.8], 'DisplayName', 'CDM (EH98 full)');
plot(k, Pk_wdm,  '--', 'LineWidth', 3, 'Color', [0.8, 0.2, 0.2], 'DisplayName', 'WDM');
plot(k, Pk_fdm,  '-.', 'LineWidth', 3, 'Color', [0.2, 0.7, 0.3], 'DisplayName', 'FDM');
xlabel('k [h/Mpc]', 'FontSize', 32);
ylabel('P(k) [(Mpc/h)^3]', 'FontSize', 32);
% title('Linear Power Spectrum: CDM vs WDM/FDM', 'FontSize', 14);
legend('Location', 'best', 'FontSize', 22);
set(gca, 'XScale', 'log', 'YScale', 'log', 'FontSize', 22, 'LineWidth', 2);
ax=gca; yl=ax.YLim; yl(2) = yl(2)+1;
ylim([yl(1), yl(2)])
% axis tight;
saveas(gcf, 'power_spectrum_wdm_fdm.png');

%% Plot 2c: Variance — filter dependence


% Reuse CDM EH98-full power as baseline
filters = {'tophat', 'gaussian', 'sharpk', 'smoothk', 'vsmk'};
colors  = [0.1 0.4 0.8;
           0.8 0.2 0.2;
           0.2 0.7 0.3;
           0.6 0.4 0.8;
           0.9 0.6 0.1];
figure; hold on; box on;
for i = 1:numel(filters)
    sig = cosmo_ehfull.sigmaR(R, z, filters{i});   % <-- use high-level wrapper
    % sig = sigmaR_ehfull(R, z, f);
    plot(R, sig, 'LineWidth', 2, 'Color', colors(i,:), 'DisplayName', filters{i});
end

xlabel('R [Mpc/h]', 'FontSize', 32);
ylabel('\sigma(R)', 'FontSize', 32);
% title('Variance: effect of filter', 'FontSize', 14);
legend('Location', 'best', 'FontSize', 22);
set(gca, 'XScale', 'log', 'YScale', 'log', 'FontSize', 22, 'LineWidth', 2);
axis tight;
saveas(gcf, 'variance_filters.png');

%% Plot 3c: Correlation Function — integral vs FFTLog
% Integral (default) method
cosmo_int = cosmo_ehfull;
cosmo_int.corr_method = 'integral';
xi_int = cosmo_int.correlationFunction(r, z);

% FFTLog method
cosmo_fft = cosmo_ehfull;
cosmo_fft.corr_method = 'fftlog';
xi_fft = cosmo_fft.correlationFunction(r, z);

% Absolute xi
figure; hold on; box on;
plot(r, abs(xi_int), '-',  'LineWidth', 3, 'Color', [0.1, 0.4, 0.8], 'DisplayName', 'Integral');
plot(r, abs(xi_fft), '--', 'LineWidth', 3, 'Color', [0.8, 0.2, 0.2], 'DisplayName', 'FFTLog');
xlabel('r [Mpc/h]', 'FontSize', 32);
ylabel('|\xi(r)|', 'FontSize', 32);
% title('Correlation: Integral vs FFTLog', 'FontSize', 14);
legend('Location', 'northeast', 'FontSize', 22);
set(gca, 'XScale', 'log', 'YScale', 'log', 'FontSize', 22, 'LineWidth', 2);
axis tight;

saveas(gcf, 'correlation_methods.png');