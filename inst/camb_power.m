function [k_tab, P_tab] = camb_power(cosmo, varargin)
% Load linear matter power spectrum from CAMB through Python
%
% INPUT:
%   cosmo : cosmology struct with fields
%       Omega_m, Omega_b, h, ns
%
% OPTIONAL NAME-VALUE PAIRS:
%   'python_exe' : full path to python executable (default: 'python')
%   'minkh'      : minimum k in h/Mpc            (default: 1e-4)
%   'maxkh'      : maximum k in h/Mpc            (default: 100)
%   'npoints'    : number of k samples           (default: 2000)
%   'As'         : primordial amplitude          (optional)
%
% OUTPUT:
%   k_tab : k array in h/Mpc
%   P_tab : linear matter power spectrum in (Mpc/h)^3

    % --- Parse optional name-value pairs (Octave-compatible) ----------
    python_exe = 'python';
    minkh      = 1e-4;
    maxkh      = 100;
    npoints    = 2000;
    As         = [];

    for i = 1:2:length(varargin)
        switch lower(varargin{i})
            case 'python_exe',  python_exe = varargin{i+1};
            case 'minkh',       minkh      = varargin{i+1};
            case 'maxkh',       maxkh      = varargin{i+1};
            case 'npoints',     npoints    = varargin{i+1};
            case 'as',          As         = varargin{i+1};
            otherwise
                error('camb_power: Unknown parameter: %s', varargin{i});
        end
    end

    % --- Cosmological parameters --------------------------------------
    H0    = 100 * cosmo.h;
    ombh2 = cosmo.Omega_b * cosmo.h^2;
    omch2 = (cosmo.Omega_m - cosmo.Omega_b) * cosmo.h^2;
    ns    = cosmo.ns;

    % --- Write a temporary Python script ------------------------------
    py_script  = [tempname, '.py'];
    output_file = [tempname, '.txt'];

    fid = fopen(py_script, 'w');
    fprintf(fid, 'import camb\n');
    fprintf(fid, 'from camb import model\n');
    fprintf(fid, 'import numpy as np\n');
    fprintf(fid, 'pars = camb.CAMBparams()\n');
    fprintf(fid, 'pars.set_cosmology(H0=%.10g, ombh2=%.10g, omch2=%.10g)\n', H0, ombh2, omch2);

    if isempty(As)
        fprintf(fid, 'pars.InitPower.set_params(ns=%.10g)\n', ns);
    else
        fprintf(fid, 'pars.InitPower.set_params(ns=%.10g, As=%.10g)\n', ns, As);
    end

    fprintf(fid, 'pars.set_matter_power(redshifts=[0.0], kmax=%.10g)\n', maxkh);
    fprintf(fid, 'pars.NonLinear = model.NonLinear_none\n');
    fprintf(fid, 'results = camb.get_results(pars)\n');
    fprintf(fid, 'kh, z, pk = results.get_matter_power_spectrum(minkh=%.10g, maxkh=%.10g, npoints=%d)\n', ...
            minkh, maxkh, npoints);
    fprintf(fid, 'out = np.column_stack([kh, pk[0]])\n');
    fprintf(fid, 'np.savetxt(r"%s", out, header="k Pk")\n', output_file);
    fclose(fid);

    % --- Call Python via system() -------------------------------------
    cmd = sprintf('"%s" "%s"', python_exe, py_script);
    [status, msg] = system(cmd);
    if status ~= 0
        delete(py_script);
        error('camb_power: Python call failed.\n%s', msg);
    end

    % --- Read output --------------------------------------------------
    data  = dlmread(output_file, '', 1, 0);   % skip header line
    k_tab = data(:, 1);
    P_tab = data(:, 2);

    % --- Clean up -----------------------------------------------------
    delete(py_script);
    delete(output_file);

    % --- Sanitise -----------------------------------------------------
    valid = isfinite(k_tab) & isfinite(P_tab) & (k_tab > 0) & (P_tab > 0);
    k_tab = k_tab(valid);
    P_tab = P_tab(valid);

    [k_tab, idx] = sort(k_tab);
    P_tab = P_tab(idx);
end