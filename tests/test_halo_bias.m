
    % clear; clc; close all;

    % ---------------------------------------------------------------------
    % 0. PATHS & SETTINGS
    % ---------------------------------------------------------------------
    addpath('F:\PHANTOM\src\utils',         '-begin');
    addpath('F:\PHANTOM\src\concentration', '-begin');
    addpath('F:\PHANTOM\src\profiles',      '-begin');
    addpath('F:\PHANTOM\src\halo',          '-begin');
    addpath('F:\PHANTOM\src\fdm',           '-begin');
    addpath('F:\PHANTOM\src\hmf',           '-begin');

    % ---------------------------------------------------------------------
    % 1. BUILD COSMOLOGY
    % ---------------------------------------------------------------------
    cosmo = cosmology('Planck18');

    % ---------------------------------------------------------------------
    % 2. BASIC TEST INPUTS
    % ---------------------------------------------------------------------
    sigma   = [2.0, 1.5, 1.0, 0.7, 0.5];      % example sigma(M,z) values
    delta_c = collapse_overdensity();
    z       = 0;
    Delta   = 200;
    M       = logspace(10, 14, numel(sigma)); % [h^-1 Msun], dummy for WDM
    Mhm     = 1e11;                           % dummy half-mode mass

    % ---------------------------------------------------------------------
    % 3. FUNCTION VISIBILITY CHECK
    % ---------------------------------------------------------------------
    fprintf('=== Function visibility check ===\n');
    fnames = { ...
        'halo_bias', ...
        'halo_bias_PS', ...
        'halo_bias_ST', ...
        'halo_bias_SMT01', ...
        'halo_bias_jing98', ...
        'halo_bias_seljak04', ...
        'halo_bias_Tinker10', ...
        'halo_bias_bhattacharya11', ...
        'halo_bias_comparat17', ...
        'halo_bias_pillepich10', ...
        'halo_bias_schneider12', ...
        'collapse_overdensity'};

    for i = 1:numel(fnames)
        fpath = which(fnames{i});
        if isempty(fpath)
            fprintf('[MISSING] %s\n', fnames{i});
        else
            fprintf('[FOUND]   %s -> %s\n', fnames{i}, fpath);
        end
    end
    fprintf('\n');

    % ---------------------------------------------------------------------
    % 4. MODEL TESTS
    % ---------------------------------------------------------------------
    fprintf('=== halo_bias branch tests ===\n');

    model_list = { ...
        'cole89',         {sigma, delta_c}; ...
        'st',             {sigma, delta_c}; ...
        'smt01',          {sigma, delta_c}; ...
        'jing98',         {sigma, delta_c, cosmo}; ...
        'seljak04',       {sigma, delta_c}; ...
        'tinker10',       {sigma, delta_c, Delta, z, cosmo}; ...
        'bhattacharya11', {sigma, delta_c, z}; ...
        'comparat17',     {sigma, delta_c}; ...
        'pillepich10',    {sigma, delta_c, cosmo}; ...
        'wdm',            {sigma, delta_c, M, Mhm}; ...
        'fdm',            {sigma, delta_c, Delta, z, cosmo}; ...
        };

    results = struct('model', {}, 'status', {}, 'size', {}, 'message', {});
    bias_store = struct();  % will hold b-vectors per model for plotting

    for i = 1:size(model_list,1)
        model = model_list{i,1};
        args  = model_list{i,2};

        fprintf('Testing model: %-15s ... ', model);

        try
            b = halo_bias(model, args{:});

            ok_numeric = isnumeric(b);
            ok_real    = isreal(b);
            ok_size    = isequal(size(b), size(sigma)) || isscalar(b);

            if ok_numeric && ok_real && ok_size
                fprintf('PASS\n');
                results(end+1).model   = model; %#ok<AGROW>
                results(end).status    = 'PASS';
                results(end).size      = mat2str(size(b));
                results(end).message   = '';

                % store for plotting; force row vector
                bias_store.(model) = b(:).'; 
            else
                fprintf('WARN\n');
                results(end+1).model   = model; %#ok<AGROW>
                results(end).status    = 'WARN';
                results(end).size      = mat2str(size(b));
                results(end).message   = 'Output returned, but size/type check was unexpected.';
            end

        catch ME
            fprintf('FAIL\n');
            results(end+1).model   = model; %#ok<AGROW>
            results(end).status    = 'FAIL';
            results(end).size      = '-';
            results(end).message   = ME.message;
        end
    end

    fprintf('\n');

    % ---------------------------------------------------------------------
    % 5. SUMMARY TABLE
    % ---------------------------------------------------------------------
    fprintf('=== Summary ===\n');
    for i = 1:numel(results)
        fprintf('%-15s : %-5s  size=%-10s  %s\n', ...
            results(i).model, results(i).status, results(i).size, results(i).message);
    end
    fprintf('\n');

    % ---------------------------------------------------------------------
    % 6. NUMERICAL SPOT CHECK
    % ---------------------------------------------------------------------
    fprintf('=== Numerical spot check ===\n');
    for i = 1:size(model_list,1)
        model = model_list{i,1};
        args  = model_list{i,2};

        try
            b = halo_bias(model, args{:});
            fprintf('%-15s -> ', model);
            disp(b);
        catch
            fprintf('%-15s -> [failed]\n', model);
        end
    end

    % ---------------------------------------------------------------------
    % 7. INVALID MODEL TEST
    % ---------------------------------------------------------------------
    fprintf('\n=== Invalid model test ===\n');
    try
        halo_bias('not_a_model', sigma, delta_c);
        fprintf('FAIL: invalid-model test did not throw an error.\n');
    catch ME
        fprintf('PASS: invalid-model test threw error as expected.\n');
        fprintf('Message: %s\n', ME.message);
    end

    % ---------------------------------------------------------------------
    % 8. PLOT b(sigma) FOR ALL SUCCESSFUL MODELS
    % ---------------------------------------------------------------------
    figure;
    hold on;
    model_names = fieldnames(bias_store);

    for i = 1:numel(model_names)
        m = model_names{i};
        b = bias_store.(m);

        % If a scalar slipped through, expand it for plotting
        if isscalar(b)
            b = repmat(b, size(sigma));
        end

        plot(sigma, b, '-o', 'DisplayName', m);
    end

    set(gca, 'XScale', 'log');
    xlabel('\sigma(M,z)');
    ylabel('b(\sigma)');
    title('Halo bias models via halo\_bias\_dispatcher');
    grid on;
    legend('Location','best');
    hold off;

