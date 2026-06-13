function [k_tab, T_tab] = load_camb_transfer(cosmo)
% Load CAMB transfer function table from the file specified in cosmo.
% Supports .mat, .csv, .txt, and .dat formats.
% Returns k_tab [h/Mpc] and T_tab [dimensionless], both column vectors,
% sorted by k and stripped of non-finite entries.

if ~isfield(cosmo, 'camb_transfer_file') || isempty(cosmo.camb_transfer_file)
    error('load_camb_transfer: cosmo.camb_transfer_file must be specified.');
end

[~, ~, ext] = fileparts(cosmo.camb_transfer_file);

switch lower(ext)

    case '.mat'
        d = load(cosmo.camb_transfer_file);

        if isfield(d, 'k_h')
            k_tab = d.k_h(:);
        elseif isfield(d, 'k')
            k_tab = d.k(:);
        elseif isfield(d, 'k_h_Mpc')
            k_tab = d.k_h_Mpc(:);
        else
            error('load_camb_transfer: MAT file must contain k_h, k, or k_h_Mpc.');
        end

        if isfield(d, 'T_camb')
            T_tab = d.T_camb(:);
        elseif isfield(d, 'T')
            T_tab = d.T(:);
        else
            error('load_camb_transfer: MAT file must contain T_camb or T.');
        end

    case {'.csv', '.txt', '.dat'}
        raw      = dlmread(cosmo.camb_transfer_file, '', 0, 0);
        % skip comment rows if present
        raw      = raw(all(isfinite(raw), 2), :);
        k_tab    = raw(:, 1);
        T_tab    = raw(:, 2);

    otherwise
        error('load_camb_transfer: Unsupported file type "%s". Use .mat, .csv, .txt, or .dat.', ext);
end

% Strip bad rows and sort
valid         = isfinite(k_tab) & isfinite(T_tab) & (k_tab > 0);
k_tab         = k_tab(valid);
T_tab         = T_tab(valid);
[k_tab, idx]  = sort(k_tab);
T_tab         = T_tab(idx);
end