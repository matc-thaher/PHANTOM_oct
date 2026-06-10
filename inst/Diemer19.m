function c = Diemer19(M, z, cosmo, mode, varargin)
% Diemer19_concentration  Diemer & Joyce (2019) concentration model
%
%   c = Diemer19_concentration(M, z, cosmo, mode)
%   c = Diemer19_concentration(M, z, cosmo, mode, profile_name)
%
%   Identical functional form and root-finding method as Ishiyama21.
%   Only the parameter table differs (Diemer19_Table vs Ishiyama21_Table2).
%
%   INPUTS
%   M            : halo mass [Msun/h], scalar or vector
%   z            : redshift (scalar)
%   cosmo        : cosmology struct from cosmology() + attach_linear_components()
%   mode         : '200c_all' | '200c_relaxed' | 'vir_all' | 'vir_relaxed'
%   profile_name : (optional) string, profile for lhs-table. Default: 'nfw'
%                  See profile_mu.m for supported profiles.
%   method       : (optional) string, solver type: 'zero' (default) or 'table'.
%
%   OUTPUT
%   c     : concentration, same shape as M
%
%   Reference: Diemer & Joyce 2019, ApJ 871, 168

    profile_name = 'nfw';
    method       = 'fzero';

    if numel(varargin) >= 1 && ~isempty(varargin{1})
        profile_name = varargin{1};
    end
    if numel(varargin) >= 2 && ~isempty(varargin{2})
        method = lower(varargin{2});
    end

    P = Diemer19_Table(mode);
    params.kappa   = P.kappa;
    params.a0      = P.a0;
    params.a1      = P.a1;
    params.b0      = P.b0;
    params.b1      = P.b1;
    params.c_alpha = P.cAlpha;

    switch method
        case 'fzero'
            c = Diemer19_zero_general(M, z, cosmo, params, profile_name);
        case 'table'
            c = Diemer19_general(M, z, cosmo, params, profile_name);
        otherwise
            error('Diemer19: unknown method \"%s\". Use ''fzero'' or ''table''.', method);
    end
end
