function deriv = dlnsigma_dlnM(M, z, cosmo, varargin)
% dlnsigma_dlnM   Logarithmic derivative d(ln sigma) / d(ln M)
%
% Uses symmetric central finite differences in log-M space.
%
% INPUT:
%   M        : mass array [M_sun/h], any shape
%   z        : redshift (scalar)
%   cosmo    : PHANTOM cosmology struct
%   varargin : optional arguments in any order —
%                filter name  (char)   e.g. 'smoothk'
%                filter_params (struct) e.g. struct('beta', 5)
%                eps_fd       (scalar) step size in ln M (default: 1e-4)
%
% OUTPUT:
%   deriv    : d(ln sigma) / d(ln M), same shape as M
%
% USAGE (all equivalent):
%   deriv = dlnsigma_dlnM(M, z, cosmo);
%   deriv = dlnsigma_dlnM(M, z, cosmo, 'smoothk');
%   deriv = dlnsigma_dlnM(M, z, cosmo, 'smoothk', 1e-5);
%   deriv = dlnsigma_dlnM(M, z, cosmo, 1e-5);              % tophat, custom eps

    % ---- Parse varargin -------------------------------------------------
    eps_fd       = 1e-4;    % default step in ln M
    filter_args  = {};      % forwarded to cosmo.sigmaM

    for i = 1:numel(varargin)
        v = varargin{i};
        if ischar(v)
            filter_args{end+1} = v;         %#ok<AGROW>  filter name
        elseif isstruct(v)
            filter_args{end+1} = v;         %#ok<AGROW>  filter_params struct
        elseif isnumeric(v) && isscalar(v)
            eps_fd = v;                     % user-supplied step size
        else
            error('dlnsigma_dlnM: unrecognised argument at position %d.', i+3);
        end
    end

    % ---- Central finite difference in ln M ------------------------------
    sig_hi = cosmo.sigmaM(M * exp( eps_fd), z, filter_args{:});
    sig_lo = cosmo.sigmaM(M * exp(-eps_fd), z, filter_args{:});

    deriv  = (log(sig_hi) - log(sig_lo)) ./ (2 * eps_fd);
end