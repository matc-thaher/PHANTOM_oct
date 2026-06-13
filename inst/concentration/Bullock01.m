function c = Bullock01(M, z, cosmo, K, F, delta_c, corrections, variant)
% Bullock01   Concentration model of Bullock et al. (2001), with the
%             Maccio et al. (2008) revision as the default.
%
%   c = Bullock01(M, z, cosmo)
%   c = Bullock01(M, z, cosmo, K, F)
%   c = Bullock01(M, z, cosmo, K, F, [], [], 'original')
%
%   The collapse redshift z_coll is defined as the redshift at which a
%   progenitor of mass F*M satisfies sigma(F*M, z_coll) = delta_c, i.e.:
%
%       D(z_coll) / D(0) = delta_c / sigma(F*M, 0)
%
%   Two variants are available for the concentration--collapse relation:
%
%   'maccio08'  [DEFAULT]  Maccio et al. (2008) revision:
%
%       c = K * [ H(z_coll) / H(z) ]^(2/3)
%
%       Physically motivated: halo density scales with the critical density
%       rho_crit ~ H^2, so concentration ~ rho_coll/rho_now ~ H^2_coll/H^2_now.
%       The exponent 2/3 follows from c ~ (rho_coll/rho_now)^(1/3) and the
%       virial relation M ~ rho * r^3. Calibrated to WMAP1/3/5 cosmologies.
%       Reproduces the c-M relation over 5 orders of magnitude in halo mass.
%
%       Reference: Maccio, Dutton & van den Bosch (2008), MNRAS 391, 1940
%                  https://ui.adsabs.harvard.edu/abs/2008MNRAS.391.1940M
%       Default parameters: K = 3.85, F = 0.01
%
%   'original'  Bullock et al. (2001) original formulation:
%
%       c = K * (1 + z_coll) / (1 + z)
%
%       Valid approximation in a matter-dominated universe where H(z) ~ (1+z)^(3/2).
%       Breaks down at low redshift when dark energy flattens H(z). Use only
%       when reproducing results from pre-2008 literature.
%
%       Reference: Bullock et al. (2001), MNRAS 321, 559
%                  https://ui.adsabs.harvard.edu/abs/2001MNRAS.321..559B
%       Default parameters: K = 4.0, F = 0.01  (original B01 calibration)
%
%   INPUTS
%   M          : halo mass [Msun/h], scalar or vector
%   z          : redshift, scalar
%   cosmo      : cosmology struct with cosmo.sigmaM(M,z), cosmo.D(z),
%                and cosmo.H(z) [km/s/Mpc]
%   K          : proportionality constant  (default depends on variant)
%   F          : collapse mass fraction    (default depends on variant)
%   delta_c    : collapse overdensity      (default: computed from cosmo)
%   corrections: logical, apply delta_c redshift corrections (default: false)
%   variant    : 'maccio08' (default) or 'original'
%
%   OUTPUT
%   c          : concentration, same shape as M

if nargin < 8 || isempty(variant), variant = 'maccio08'; end

switch variant
    case 'maccio08'
        if nargin < 4 || isempty(K), K = 3.85; end
        if nargin < 5 || isempty(F), F = 0.01;  end
    case 'original'
        if nargin < 4 || isempty(K), K = 4.0;  end
        if nargin < 5 || isempty(F), F = 0.01;  end
    otherwise
        error('Bullock01: unknown variant "%s". Use ''maccio08'' or ''original''.', variant);
end

if nargin < 7 || isempty(corrections), corrections = false;   end
if nargin < 6 || isempty(delta_c)
        delta_c = collapse_overdensity('corrections', corrections, 'z', z);
end

% delta_c = 1.686;
M       = M(:);
c       = zeros(size(M));

for i = 1:numel(M)

    % Mass of the collapsing progenitor
    M_prog = F * M(i);

    % sigma(M_prog) at z=0  (growth factor D(z)/D(0) scales sigma)
    sigma0 = cosmo.sigmaM(M_prog, 0);

    % z_coll satisfies: D(z_coll)/D(0) * sigma0 = delta_c
    % => D(z_coll)/D(0) = delta_c / sigma0
    D_ratio_target = delta_c / sigma0;

    if D_ratio_target >= 1
        % z_coll lies in the future. Search for a bracket automatically.
        obj = @(zc) cosmo.D(zc) / cosmo.D(0) - D_ratio_target;
        % Try increasingly distant future epochs until a sign change is found
        z_future_candidates = [-0.1, -0.3, -0.5, -0.7, -0.9, -0.95, -0.99];
        bracket_found = false;
        for jj = 1:numel(z_future_candidates)
            z_lo = z_future_candidates(jj);
            if obj(z_lo) * obj(0) < 0
                z_coll = fzero(obj, [z_lo, 0]);
                bracket_found = true;
                break;
            end
        end
        if ~bracket_found
            % D(z) never reaches D_ratio_target even in far future
            % Halo will never collapse — assign invalid concentration
            z_coll = z_future_candidates(end);
            c(i) = 0;
            continue;
        end
    else
        % Invert D(z_coll)/D(0) = D_ratio_target
        obj    = @(zc) cosmo.D(zc) / cosmo.D(0) - D_ratio_target;
        z_coll = fzero(obj, [0, 30]);
    end

    switch variant
        case 'maccio08'
            c(i) = K * (cosmo.Hz(z_coll) / cosmo.Hz(z))^(2/3);
        case 'original'
            c(i) = K * (1 + z_coll) / (1 + z);
    end
end

end