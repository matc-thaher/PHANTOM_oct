function delta_c = collapse_overdensity(varargin)
% collapse_overdensity   Linear overdensity threshold for spherical collapse
%
% In an Einstein-de Sitter universe: delta_c = (3/5)*(3*pi/2)^(2/3) = 1.6865.
% Small corrections apply in non-EdS cosmologies (< 3% for any realistic
% flat LCDM cosmology; Kitayama & Suto 1996, Eq. A6).
%
% INPUT (all optional, name-value pairs):
%   'corrections'  logical  Apply Omega_m(z) correction (default: false)
%   'z'            scalar   Redshift; required if corrections = true
%   'cosmo'        struct   PHANTOM cosmo struct; required if corrections = true
%
% OUTPUT:
%   delta_c : collapse overdensity (scalar)
%
% USAGE:
%   delta_c = collapse_overdensity()
%   delta_c = collapse_overdensity('corrections', true, 'z', 0.5, 'cosmo', cosmo)

    %% ---- Defaults -------------------------------------------------------
    corrections = false;
    z           = [];
    cosmo       = [];

    %% ---- Parse name-value pairs -----------------------------------------
    i = 1;
    while i <= numel(varargin)
        key = varargin{i};
        val = varargin{i+1};
        switch lower(key)
            case 'corrections', corrections = val;
            case 'z',           z           = val;
            case 'cosmo',       cosmo       = val;
            otherwise
                error('collapse_overdensity: unknown argument "%s".', key);
        end
        i = i + 2;
    end

    % ---- EdS value ------------------------------------------------------
    delta_c = (3/5) * (3*pi/2)^(2/3);   % = 1.68647...

    % ---- Non-EdS correction (Kitayama & Suto 1996 / Colossus convention) -
    if corrections
        if isempty(z) || isempty(cosmo)
            error('collapse_overdensity: ''z'' and ''cosmo'' required when corrections = true.');
        end
        Om_z = cosmo.Omega_m_z(z);      % Omega_m(z) from PHANTOM cosmo handle
        if cosmo.flat
            delta_c = delta_c * Om_z^0.0055;
        else
            % Open universe without dark energy (Eke+1996 / Colossus)
            delta_c = delta_c * Om_z^0.0185;
        end
    end
end