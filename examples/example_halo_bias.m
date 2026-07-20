
    clear; clc; close all;

    % Paths
    addpath('F:\PHANTOM\src\utils',         '-begin');
    addpath('F:\PHANTOM\src\concentration', '-begin');
    addpath('F:\PHANTOM\src\profiles',      '-begin');
    addpath('F:\PHANTOM\src\halo',          '-begin');
    addpath('F:\PHANTOM\src\fdm',           '-begin');

    % Cosmology
    cosmo = cosmology('Planck18');

    % Mass grid [M_sun/h]
    log10M = linspace(10, 15, 200);
    M      = 10.^log10M;

    % Redshifts to illustrate
    z_list = [0.0, 0.5, 1.0, 2.0];

    % Collapse threshold
    delta_c = collapse_overdensity();  % or cosmo.delta_c(0) if you prefer

    % Containers
    b_tinker = zeros(numel(z_list), numel(M));

    for iz = 1:numel(z_list)
        z = z_list(iz);

        % sigma(M,z) from your PHANTOM pipeline; adjust to your actual helper
        % e.g. sigmaM(M,z,cosmo) or cosmo.sigmaM(M,z)
        sigma = cosmo.sigmaM(M, z);

        % Tinker+2010 bias via dispatcher (default CDM model)
        % Note: Delta=200, pass redshift and cosmo
        b_tinker(iz, :) = halo_bias_dispatcher('tinker10', sigma, delta_c, 200, z, cosmo);
    end

    % Optional: second model for comparison at z=0
    sigma0 = cosmo.sigmaM(M, 0);
    b_ST   = halo_bias_dispatcher('st', sigma0, delta_c);

    % Plot
    figure;
    hold on;

    cmap = lines(numel(z_list));
    for iz = 1:numel(z_list)
        plot(log10M, b_tinker(iz,:), 'Color', cmap(iz,:), ...
             'LineWidth', 1.5, 'DisplayName', sprintf('Tinker10, z=%.1f', z_list(iz)));
    end

    % Overplot ST at z=0
    plot(log10M, b_ST, '--k', 'LineWidth', 1.5, 'DisplayName', 'Sheth-Tormen (z=0)');

    xlabel('log_{10}(M / M_\odot h^{-1})');
    ylabel('Linear halo bias b(M,z)');
    title('PHANTOM: halo bias vs mass');
    legend('Location','northwest');
    grid on;
    box on;

    hold off;
