function c = c_CDM(M, z, model, varargin)
% c_CDM  Unified CDM concentration-mass relation dispatcher
%
%   c = c_CDM(M, z)
%   c = c_CDM(M, z, model, cosmo)
%   c = c_CDM(M, z, model, cosmo, mode)
%
%   Dispatches to one of 11 concentration-mass models. Default: ishiyama21.
%
%   INPUTS
%   M     : halo mass [Msun/h], scalar or vector
%   z     : redshift (scalar)
%   model : string identifying the model (see list below)
%   cosmo : cosmology struct from cosmology.m  (required for physics-based models)
%   mode  : model-specific sample/fit string   (optional, uses default if omitted)
%
%   MODELS & SIGNATURES
%   -------------------
%   Simple fits — no cosmo needed:
%     'duffy08'     Duffy   et al.  2008  — c = c_CDM(M, z, 'duffy08',  cosmo, mode)
%     'klypin11'    Klypin  et al.  2011  — c = c_CDM(M, z, 'klypin11', cosmo, mode)
%     'prada12'     Prada   et al.  2012  — c = c_CDM(M, z, 'prada12',  cosmo)
%
%   Physics-based fits — cosmo required:
%     'bullock01'   Bullock et al.  2001  — c = c_CDM(M, z, 'bullock01', cosmo)
%     'dutton14'    Dutton & Maccio 2014  — c = c_CDM(M, z, 'dutton14',   cosmo, mode)
%     'diemer15'    Diemer & Kravtsov 2015— c = c_CDM(M, z, 'diemer15',   cosmo, mode)
%     'ludlow16'    Ludlow  et al.  2016  — c = c_CDM(M, z, 'ludlow16',   cosmo)
%     'klypin16'    Klypin  et al.  2016  — c = c_CDM(M, z, 'klypin16',   cosmo, mode)
%     'child18'     Child   et al.  2018  — c = c_CDM(M, z, 'child18',    cosmo, mode)
%     'diemer19'    Diemer & Joyce  2019  — c = c_CDM(M, z, 'diemer19',   cosmo, mode)
%     'ishiyama21'  Ishiyama et al. 2021  — c = c_CDM(M, z, 'ishiyama21', cosmo, mode)
%
%   DEFAULT MODES (used when mode is omitted)
%   ------------------------------------------
%   duffy08     : 'full_all'
%   klypin11    : 'all'
%   dutton14    : 'all'
%   diemer15    : 'median'
%   klypin16    : 'all'
%   child18     : 'individual_all'
%   diemer19    : 'median'
%   ishiyama21  : '200c_all'
%
%   Examples
%   --------
%   c = c_CDM(M, z)
%   c = c_CDM(M, z, 'bullock01')
%   c = c_CDM(M, z, 'duffy08')
%   c = c_CDM(M, z, 'duffy08', '200c_NFW_full_z0_2')
%   c = c_CDM(M, z, 'duffy08', '200c_NFW_relaxed_z0')
%   c = c_CDM(M, z, 'duffy08', 'vir_NFW_full_z0_2')
%   c = c_CDM(M, z, 'klypin11')
%   c = c_CDM(M, z, 'klypin11', 'distinct')
%   c = c_CDM(M, z, 'klypin11', 'subhalo')
%   c = c_CDM(M, z, 'prada12',     cosmo)
%   c = c_CDM(M, z, 'dutton14')
%   c = c_CDM(M, z, 'dutton14', '200c')
%   c = c_CDM(M, z, 'dutton14', 'vir')
%   c = c_CDM(M, z, 'diemer15',    cosmo)
%   c = c_CDM(M, z, 'diemer15',    cosmo, 'mean')
%   c = c_CDM(M, z, 'ludlow16',    cosmo)
%   c = c_CDM(M, z, 'klypin16',  cosmo)
%   c = c_CDM(M, z, 'klypin16',  cosmo, 'planck13_200c_cM')
%   c = c_CDM(M, z, 'klypin16',  cosmo, 'bolshoi_200c_cM')
%   c = c_CDM(M, z, 'klypin16',  cosmo, 'planck13_vir_cnu')
%   c = c_CDM(M, z, 'child18',     cosmo)
%   c = c_CDM(M, z, 'child18',     cosmo, 'individual_relaxed')
%   c = c_CDM(M, z, 'child18',     cosmo, 'stack_nfw')
%   c = c_CDM(M, z, 'child18',     cosmo, 'stack_einasto')
%   c = c_CDM(M, z, 'diemer19',    cosmo)
%   c = c_CDM(M, z, 'diemer19',    cosmo, 'mean')
%   c = c_CDM(M, z, 'ishiyama21',  cosmo)
%   c = c_CDM(M, z, 'ishiyama21',  cosmo, '200c_relaxed')
%   c = c_CDM(M, z, 'ishiyama21',  cosmo, '500_all')
%   c = c_CDM(M, z, 'ishiyama21',  cosmo, 'vir_relaxed')
%
%   References
%   ----------
%   Bullock  et al. 2001, MNRAS 321, 559
%   Duffy    et al. 2008, MNRAS 390, L64
%   Klypin   et al. 2011, ApJ  740, 102
%   Prada    et al. 2012, MNRAS 423, 3018
%   Dutton & Maccio 2014, MNRAS 441, 3359
%   Diemer & Kravtsov 2015, ApJ 799, 108
%   Ludlow   et al. 2016, MNRAS 460, 1214
%   Klypin   et al. 2016, MNRAS 457, 4340
%   Child    et al. 2018, ApJ  859, 55
%   Diemer & Joyce  2019, ApJ  871, 168
%   Ishiyama et al. 2021, MNRAS 506, 4210

if nargin < 3 || isempty(model)
    model = 'ishiyama21';
end

switch lower(model)

    % ================================================================
    % Bullock et al. 2001
    % ================================================================
    case {'bullock01', 'b01'}
        [cosmo, ~] = parse_args(varargin, '', 'bullock01');
        c = Bullock01(M, z, cosmo);

    % ================================================================
    % Duffy et al. 2008  — no cosmo, mode = 'mdef_profile_sample_zrange'
    % e.g. '200c_NFW_full_z0_2' or '200c_NFW_relaxed_z0_2'
    % ================================================================
    case {'duffy08', 'd08'}
        % mode string: 'mdef_profile_sample_zrange'
        % defaults:     200c   NFW     full    z0_2
        if isempty(varargin) || isstruct(varargin{1})
            mdef  = '200c'; profile = 'NFW'; sample = 'full'; zrange = 'z0_2';
        else
            mode  = varargin{1};
            parts = strsplit(mode, '_');
            mdef    = parts{1};                                          % 200c | vir | 200m
            profile = parts{2};                                          % NFW | Einasto
            sample  = parts{3};                                          % full | relaxed
            zrange  = strjoin(parts(4:end), '_');                        % z0 | z0_2
        end
        c = Duffy08(M, z, mdef, profile, sample, zrange);

    % ================================================================
    % Klypin et al. 2011  — no cosmo, third arg is sample
    % sample: 'distinct' (default) | 'subhalo'
    % ================================================================
    case {'klypin11', 'k11'}
        if isempty(varargin) || isstruct(varargin{1})
            sample = 'distinct';
        else
            sample = varargin{1};
        end
        c = Klypin11(M, z, sample);

    % ================================================================
    % Prada et al. 2012
    % ================================================================
    case {'prada12', 'p12'}
        [cosmo, ~] = parse_args(varargin, '', 'prada12');
        c = Prada12(M, z, cosmo);

    % ================================================================
    % Dutton & Maccio 2014  — no cosmo, third arg is mdef
    % mdef: '200c' (default) | 'vir'
    % ================================================================
    case {'dutton14', 'dm14'}
        if isempty(varargin) || isstruct(varargin{1})
            mdef = '200c';
        else
            mdef = varargin{1};
        end
        c = Dutton14(M, z, mdef);

    % ================================================================
    % Diemer & Kravtsov 2015
    % ================================================================
    case {'diemer15', 'dk15'}
        [cosmo, mode] = parse_args(varargin, 'median', 'diemer15');
        c = Diemer15(M, z, cosmo, mode);

    % ================================================================
    % Ludlow et al. 2016  — no mode argument
    % ================================================================
    case {'ludlow16', 'l16'}
        [cosmo, ~] = parse_args(varargin, '', 'ludlow16');
        c = Ludlow16(M, z, cosmo);

    case {'ludlow16_fit', 'l16_fit'}
        [cosmo, ~] = parse_args(varargin, '', 'ludlow16_fit');
        c = Ludlow16_fit(M, z, cosmo);
    % ================================================================
    % Klypin et al. 2016
    % mode string: 'cosmo_name_mdef_formula'
    % defaults:     planck13    200c   cM
    % e.g. 'planck13_200c_cM' | 'bolshoi_vir_cnu'
    % ================================================================
    % case {'klypin16', 'k16'}
    %     % [cosmo, mode] = parse_args(varargin, 'planck13_200c_cM', 'klypin16');
    %     % parts       = strsplit(mode, '_');
    %     % cosmo_name  = parts{1};          % planck13 | bolshoi
    %     % mdef_k16    = parts{2};          % 200c | vir
    %     % formula     = parts{3};          % cM | cnu
    %     % [c, ~]      = Klypin16(M, z, cosmo, cosmo_name, mdef_k16, formula);
    %     [cosmo, cosmo_name, mdef_k16, formula] = parse_args(varargin, ...
    %                                          'planck13', '200c', 'cM', 'klypin16');
    %     [c, ~] = Klypin16(M, z, cosmo, cosmo_name, mdef_k16, formula);
    case {'klypin16', 'k16'}
        % Extract cosmo struct — parse_args only handles one mode string,
        % so we parse the remaining args manually here.
        [cosmo, ~] = parse_args(varargin, '', 'klypin16');

        % Defaults
        cosmo_name = 'planck13';
        mdef_k16   = '200c';
        formula    = 'cM';

        % Read remaining varargin entries (skip cosmo struct if present)
        extra = varargin;
        if ~isempty(extra) && isstruct(extra{1})
            extra = extra(2:end);   % remove cosmo struct from front
        end

        % Match each string to its parameter by content
        for i = 1:numel(extra)
            v = extra{i};
            if ~ischar(v), continue; end
            if any(strcmp(v, {'planck13','bolshoi'}))
                cosmo_name = v;
            elseif any(strcmp(v, {'200c','vir'}))
                mdef_k16 = v;
            elseif any(strcmp(v, {'cM','cnu'}))
                formula = v;
            else
                error('c_CDM klypin16: unrecognised argument "%s".', v);
            end
        end

        [c, ~] = Klypin16(M, z, cosmo, cosmo_name, mdef_k16, formula);

    % ================================================================
    % Child et al. 2018
    % ================================================================
    case {'child18', 'c18'}
        [cosmo, mode] = parse_args(varargin, 'individual_all', 'child18');
        c = Child18(M, z, cosmo, mode);

       % ================================================================
    % Diemer & Joyce 2019
    % ================================================================
    case {'diemer19', 'dj19'}
        [cosmo, mode, profile_name] = parse_args(varargin, 'median', 'diemer19');
        c = Diemer19(M, z, cosmo, mode, profile_name);
    % ================================================================
    % Ishiyama et al. 2021  (default)
    % ================================================================
    case {'ishiyama21', 'ish21'}
        [cosmo, mode, profile_name] = parse_args(varargin, 'vir_all', 'ishiyama21');
        c = Ishiyama21(M, z, cosmo, mode, profile_name);

    % ================================================================
    % Ishiyama et al. 2021 — direct fzero implementation
    % ================================================================
    case {'ishiyama21_zero', 'ish21z'}
        [cosmo, mode, ~] = parse_args(varargin, '200c_all', 'ishiyama21_zero');
        c = Ishiyama21_zero(M, z, cosmo, mode);

    % ================================================================
    otherwise
        error(['c_CDM: unknown model "%s".\n' ...
               'Valid: bullock01, duffy08, klypin11, prada12, dutton14,\n' ...
               '       diemer15, ludlow16, klypin16, child18, diemer19, ishiyama21, ishiyama21_zero.'], ...
               model);
end
end


% =========================================================================
function [cosmo, mode, profile_name] = parse_args(v, default_mode, model_name)
    if isempty(v)
        error('c_CDM: model "%s" requires cosmo as 4th argument.', model_name);
    end
    cosmo = v{1};
    if numel(v) >= 2 && ~isempty(v{2})
        mode = v{2};
    else
        mode = default_mode;
    end
    if numel(v) >= 3 && ~isempty(v{3})
        profile_name = v{3};
    else
        profile_name = 'nfw';
    end
end