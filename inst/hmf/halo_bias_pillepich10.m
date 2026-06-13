function [b, varargout] = halo_bias_pillepich10(sigma, cosmo, mode, varargin)
% Pillepich, Porciani & Hahn (2010, MNRAS 402, 191) halo bias.
%
% MODE = 'gaussian' [default]:
%   Linear Gaussian bias as a polynomial fit in sigma^{-1} (Eq. 14-15).
%   Valid for -1.1 < ln(sigma^{-1}) < 0.8, i.e. roughly
%   2.4e10 < M/[h^{-1}Msun] < 5e15 at z=0.
%
%   Usage:
%       b0 = halo_bias_pillepich10(sigma, cosmo)
%       b0 = halo_bias_pillepich10(sigma, cosmo, 'gaussian')
%
% MODE = 'nongaussian':
%   Total bias b(k) = b0 + Delta_b with simulation-calibrated beta
%   correction (Eq. 16-20). Valid for 0.01 < k < 0.05 h/Mpc
%   and -80 <= f_NL <= 750.
%
%   Usage:
%       [b, Delta_b] = halo_bias_pillepich10(sigma, cosmo, 'nongaussian', ...
%                          k, f_NL)
%       [b, Delta_b] = halo_bias_pillepich10(sigma, cosmo, 'nongaussian', ...
%                          k, f_NL, z)   % default z=0
%
%   Extra inputs (non-Gaussian mode only):
%       k    : wavenumber [h/Mpc], same size as sigma
%       f_NL : primordial non-Gaussianity parameter (scalar)
%       z    : redshift [optional, default 0]
%
%   Outputs:
%       b       : total bias b0 + Delta_b
%       Delta_b : non-Gaussian correction alone (optional second output)
%
% Cosmological quantities drawn automatically from cosmo struct:
%   cosmo.T(k)      -- matter transfer function
%   cosmo.D(z)      -- linear growth factor, normalized D(0)=1
%   cosmo.Omega_m   -- matter density parameter at z=0
%   cosmo.H0        -- Hubble constant [km/s/Mpc]
%
% Notes:
%   - The Gaussian fit is a direct polynomial calibration on Pillepich+2010
%     simulations, NOT a Tinker/SMT functional form.
%   - The beta correction (Eq. 19-20) is essential: the analytical
%     peak-background split (Dalal+2008) overestimates Delta_b by
%     20-70% at k > 0.01 h/Mpc without it.
%   - f_NL convention: defined at z=infinity (Section 2.2 of the paper).
%     The g(inf)/g(0) = 1.3064 factor for WMAP5 is absorbed into alpha.
%
% Reference: Pillepich, Porciani & Hahn 2010, MNRAS 402, 191
%            arXiv:0811.4203

    if nargin < 3 || isempty(mode)
        mode = 'gaussian';
    end
    mode = lower(mode);

    % ================================================================
    % GAUSSIAN BIAS: Eq. 14-15
    % b0(sigma) = B0 + B1*sigma^{-1} + B2*sigma^{-2}
    % ================================================================
    B0 =  0.647;
    B1 = -0.540;
    B2 =  1.614;

    inv_sigma = 1 ./ sigma;
    b0 = B0 + (B1 .* inv_sigma) + (B2 .* inv_sigma.^2);

    if strcmp(mode, 'gaussian')
        b = b0;
        return
    end

    % ================================================================
    % NON-GAUSSIAN CORRECTION: Eq. 16-20
    % Delta_b = beta(k, f_NL) * f_NL * (b0 - 1) * Gamma / alpha(k,z)
    % ================================================================
    if ~strcmp(mode, 'nongaussian')
        error('halo_bias_pillepich10: mode must be ''gaussian'' or ''nongaussian''.');
    end

    if numel(varargin) < 2
        error('halo_bias_pillepich10: non-Gaussian mode requires k and f_NL.');
    end

    k    = varargin{1};   % [h/Mpc]
    f_NL = varargin{2};   % scalar

    if numel(varargin) >= 3
        z = varargin{3};
    else
        z = 0;
    end

    % --- Draw cosmological quantities from cosmo struct ---
    Omega_m = cosmo.Omega_m;
    H0      = cosmo.H0;          % [km/s/Mpc]
    Tk      = cosmo.T(k);        % transfer function evaluated at k
    Dz      = cosmo.D(z);        % growth factor at redshift z

    % --- Physical constants ---
    c_light  = 299792.458;       % [km/s]
    c_over_H = c_light / H0;    % Hubble radius [Mpc/h]

    % --- Gamma (Eq. 17, Poisson kernel prefactor) ---
    delta_c = 1.686;
    Gamma   = 3 .* delta_c .* Omega_m ./ (c_over_H.^2);   % [(h/Mpc)^2]

    % --- g(inf)/g(0) normalization for f_NL defined at z=inf (Sec. 2.2) ---
    g_ratio = 1.3064;   % WMAP5 value; small sensitivity to cosmology

    % --- alpha(k, z): Eq. 17 ---
    alpha = k.^2 .* Tk .* Dz ./ g_ratio;                  % [(h/Mpc)^2]

    % --- beta(k, f_NL): simulation-calibrated correction, Eq. 19-20 ---
    % Corrects for the overestimate of the analytical peak-background split.
    beta0 = 1.029;
    beta1 = 4.25e-4;   % [f_NL^{-1}]
    beta2 = 14.8;      % [h^{-1} Mpc], k must be in [h/Mpc]

    beta  = beta0 .* (1 - (beta1 .* f_NL)) .* (1 - (beta2 .* k));

    % --- Scale-dependent non-Gaussian correction (Eq. 18) ---
    Delta_b = beta .* f_NL .* (b0 - 1) .* Gamma ./ alpha;

    b = b0 + Delta_b;

    if nargout > 1
        varargout{1} = Delta_b;
    end

end