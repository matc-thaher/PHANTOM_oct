function [result, hc] = halo_obs(quantity_in, value_in, quantity_out, cosmo, varargin)
% HALO_OBS  Compute halo observable properties from a single input quantity.
%
%   result = halo_obs(quantity_in, value_in, quantity_out, cosmo)
%   result = halo_obs(quantity_in, value_in, quantity_out, cosmo, z)
%   result = halo_obs(quantity_in, value_in, quantity_out, cosmo, z, opts)
%
%   Dispatcher that converts any supported halo quantity to any other,
%   routing through M_h as the canonical intermediate.
%
%   The PHANTOM cosmo struct (from cosmology() + derive_cosmo_params() +
%   attach_linear_components()) is the sole source of cosmological
%   background. No parameters are duplicated here.
%
% --------------------------------------------------------------------------
% INPUTS
%   quantity_in  : string — name of the input quantity (see table below)
%   value_in     : scalar — value of the input in the units listed below
%   quantity_out : string or cell array of strings — desired output(s)
%   cosmo        : PHANTOM cosmology struct
%   z            : redshift (default: 0)
%   opts         : struct with optional fields
%       .overdensity   'bryan_norman' (default) | '200c' | '200m' | '500c'
%       .concentration numeric value  — use this fixed c (skips model)
%                      string name    — passed directly to c_CDM dispatcher
%                      '' or NaN      — use default model ('ishiyama21')
%
% SUPPORTED quantity_in / quantity_out STRINGS
%   'M_h'        Virial halo mass                               [M_sun]
%   'r_vir'      Virial radius                                  [kpc]
%   'V_vir'      Circular velocity at r_vir                     [km/s]
%   'V_max'      Peak circular velocity of NFW halo             [km/s]
%   'sigma_v'    1-D velocity dispersion                        [km/s]
%   'T_vir'      Virial temperature (mu=0.59)                   [K]
%   'E_kin'      Kinetic energy  (virial: E_kin = -E_pot/2)     [M_sun (km/s)^2]
%   'E_pot'      Gravitational potential energy                  [M_sun (km/s)^2]
%   'E_tot'      Total energy E_kin + E_pot                     [M_sun (km/s)^2]
%   't_dyn'      Dynamical (free-fall) time                     [Gyr]
%   'rho_mean'   Mean density within r_vir                      [M_sun/kpc^3]
%   'rho_s'      NFW scale density rho_s                        [M_sun/kpc^3]
%   'r_s'        NFW scale radius                               [kpc]
%   'concentration' NFW concentration c = r_vir/r_s             [dimensionless]
%   'r_200c'     Radius where mean density = 200*rho_crit        [kpc]
%   'M_200c'     Mass within r_200c                             [M_sun]
%   'nu'         Peak height nu = delta_c / sigma(M,z)          [dimensionless]
%   'delta_vir'  Bryan & Norman Delta_vir(z)                    [dimensionless]
%   'Omega_m_z'  Matter density parameter Omega_m(z)            [dimensionless]
%   'zeta_z'     zeta(z) = Delta_vir*Omega_m(z)/(18*pi^2)       [dimensionless]
%   'rho_c_z'    Critical density at z                          [M_sun/kpc^3]
%
% --------------------------------------------------------------------------
% EQUATION REFERENCE MAP
%   E(z):           Friedmann eq.; Diemer 2018 (COLOSSUS) Eq.(7)
%   Omega_m(z):     standard LCDM; Peebles 1993 Eq.(5.111)
%   rho_crit(z):    3*H(z)^2 / (8*pi*G)
%   Delta_vir:      Bryan & Norman 1998, ApJ 495, 80, Eq.(6); flat LCDM
%                     x = Omega_m(z)-1; Delta = 18*pi^2 + 82*x - 39*x^2
%   zeta(z):        Bryan & Norman 1998, Eq.(A2);
%                   Schive+2014 PRL after Eq.(4)
%   r_vir:          M = (4*pi/3)*Delta*rho_crit*r_vir^3
%                   Robles+2018 Sec.2.1; Schive+2014 PRL after Eq.(4)
%   V_vir:          V_vir = sqrt(G*M_h/r_vir)
%   V_max (NFW):    V_max = V_vir*sqrt(0.216*c/f(c))
%                   Bullock+2001, MNRAS 321, 559, Eq.(3)
%   sigma_v:        sigma_v = V_vir/sqrt(2)  (isotropic isothermal)
%   T_vir:          T_vir = mu*m_p*V_vir^2/(2*k_B); mu=0.59
%                   Voit 2005, Rev.Mod.Phys. 77, 207, Eq.(2)
%   E_pot,E_kin:    virial theorem; Padmanabhan 1993 / Voit 2005 Eq.(2)
%   t_dyn:          t_dyn = sqrt(3*pi/(16*G*rho_mean))
%                   Padmanabhan 1993; Bryan & Norman 1998, Eq.(2)
%   NFW rho_s:      NFW 1997, ApJ 490, 493, Eq.(7)
%   Concentration:  c_CDM dispatcher (PHANTOM src/concentration/)
%                   Default: Ishiyama et al. 2021, MNRAS 506, 4210
%   nu:             delta_c/sigma(M,z); sigma from cosmo.sigmaM(M,z)
%                   Press & Schechter 1974; Lacey & Cole 1994
%
% --------------------------------------------------------------------------

    % ---- Parse optional arguments ----------------------------------------
    z    = 0.0;
    opts = struct();
    if nargin >= 5, z    = varargin{1}; end
    if nargin >= 6, opts = varargin{2}; end

    % Default opts
    if ~isfield(opts, 'overdensity'),   opts.overdensity   = 'bryan_norman'; end
    if ~isfield(opts, 'concentration'), opts.concentration = NaN;            end
    if ~isfield(opts, 'phantom_root'),  opts.phantom_root  = '';             end

    % ---- Add PHANTOM src folders to path (Octave + MATLAB compatible) ----
    %
    %   By default halo_obs assumes PHANTOM src folders are already on the
    %   path (typical when the caller has run startup.m or called addpath).
    %
    %   To let halo_obs manage the path itself, set:
    %       opts.phantom_root = '/path/to/PHANTOM';
    %   The five required sub-folders are then added once, only if absent.
    %
    %   Sub-folders added:
    %     src/concentration  — c_CDM and all concentration model files
    %     src/utils          — cosmology, derive_cosmo_params, etc.
    %     src/profiles       — NFW_profile, Soliton_profile, etc.
    %     src/chmr           — schive_CHMR, thaher_CHMR
    %     src/fdm            — soliton_nfw_analytic, etc.
    if ~isempty(opts.phantom_root)
        phantom_subdirs = { ...
            'src/concentration', ...
            'src/utils',         ...
            'src/profiles',      ...
            'src/chmr',          ...
            'src/fdm'            ...
        };
        for kd = 1:numel(phantom_subdirs)
            p = fullfile(opts.phantom_root, phantom_subdirs{kd});
            % exist(...,'dir') and isempty(strfind(...)) work in both
            % Octave >= 4.0 and MATLAB — no pkg or toolbox required
            if exist(p, 'dir') && isempty(strfind(path(), p))
                addpath(p);
            end
        end
    end

    % ---- Build halo cosmology struct (Bryan-Norman + unit constants) ------
    hc = halo_cosmo(cosmo, z, opts.overdensity);

    % ---- Purely cosmological outputs (no M_h needed) ---------------------
    cosmo_only = {'delta_vir', 'omega_m_z', 'zeta_z', 'rho_c_z'};

    % ---- Step 1: convert input -> M_h [M_sun] ----------------------------
    if ismember(lower(quantity_in), cosmo_only)
        M_h = NaN;
    else
        M_h = to_M_h(quantity_in, value_in, hc, cosmo, z, opts);
    end

    % ---- Step 2: convert M_h -> requested output(s) ----------------------
    if ischar(quantity_out)
        quantity_out = {quantity_out};
        scalar_out   = true;
    else
        scalar_out   = false;
    end

    result = struct();
    for k = 1:numel(quantity_out)
        qo = quantity_out{k};
        if ismember(lower(qo), cosmo_only)
            result.(qo) = from_cosmo(qo, hc);
        else
            result.(qo) = from_M_h(qo, M_h, hc, cosmo, z, opts);
        end
    end

    if scalar_out
        fname  = fieldnames(result);
        result = result.(fname{1});
    end
end

% ==========================================================================
% HALO_COSMO  —  Bryan-Norman overdensity + unit constants
% ==========================================================================
function hc = halo_cosmo(cosmo, z, overdensity)
% halo_cosmo  Augment PHANTOM cosmo with halo-boundary quantities at redshift z.
%
%   hc = halo_cosmo(cosmo, z, overdensity)
%
%   INPUT
%     cosmo       – PHANTOM struct from cosmology() + derive_cosmo_params()
%     z           – redshift (scalar)
%     overdensity – 'bryan_norman'|'vir' | '200c' | '200m' | '500c'
%
%   UNIT CONVENTION
%     mass [M_sun], length [kpc], velocity [km/s]
%     G = 4.3009e-6  kpc (km/s)^2 M_sun^-1
%
%   CONSTANTS  (not in PHANTOM cosmo)
%     G_kpc  = 4.3009e-6   kpc (km/s)^2 / M_sun
%     Gyr_s  = 3.1557e16   s / Gyr
%     km_kpc = 3.0857e16   km / kpc
%
%   COSMOLOGICAL BACKGROUND  (delegated to PHANTOM handles)
%     E(z)        = cosmo.E(z)             Friedmann; Diemer 2018 Eq.(7)
%     H(z)        = cosmo.Hz(z) [km/s/Mpc] / 1e3  -> [km/s/kpc]
%     Omega_m(z)  = cosmo.Omega_m_z(z)
%     rho_crit(z) = 3*H(z)^2 / (8*pi*G)   [M_sun/kpc^3]  — recomputed
%                   from Hz to avoid h-unit ambiguity in cosmo.rhocrit(z)
%                   (PHANTOM stores rho_crit0 = 2.775e11 M_sun h^2 Mpc^-3)
%
%   OVERDENSITY THRESHOLD  (Bryan & Norman 1998, ApJ 495, 80)
%     x         = Omega_m(z) - 1
%     Delta_BN  = 18*pi^2 + 82*x - 39*x^2          Eq.(6), flat LCDM
%     zeta(z)   = Delta_BN * Omega_m(z) / (18*pi^2) Eq.(A2)
%                 [= 1 at z=0 in EdS]

    G_kpc  = 4.3009e-6;   % kpc (km/s)^2 / M_sun
    Gyr_s  = 3.1557e16;   % s / Gyr
    km_kpc = 3.0857e16;   % km / kpc

    % --- Delegate to PHANTOM function handles ---
    Ez     = cosmo.E(z);
    Hz_Mpc = cosmo.Hz(z);          % km/s/Mpc  (PHANTOM convention)
    Hz     = Hz_Mpc / 1e3;         % km/s/kpc
    Om_z   = cosmo.Omega_m_z(z);

    % rho_crit recomputed from Hz to avoid h-unit conversion on cosmo.rhocrit
    %   rho_crit [M_sun/kpc^3] = 3*H^2/(8*pi*G)
    %   H in [km/s/kpc], G in [kpc (km/s)^2 / M_sun]
    rho_crit = 3 * Hz^2 / (8 * pi * G_kpc);

    % --- Bryan & Norman 1998, ApJ 495, 80, Eq.(6), flat LCDM ---
    x        = Om_z - 1;
    Delta_BN = 18*pi^2 + 82*x - 39*x^2;

    % --- Select overdensity ---
    switch lower(overdensity)
        case {'bryan_norman', 'vir'}
            Delta = Delta_BN;
        case '200c'
            Delta = 200;
        case '200m'
            Delta = 200 * Om_z;   % effective Delta wrt rho_crit
        case '500c'
            Delta = 500;
        otherwise
            error('halo_cosmo: unknown overdensity ''%s''.', overdensity);
    end

    % --- Pass-through scalars from PHANTOM cosmo ---
    hc.H0    = cosmo.H0;        % km/s/Mpc
    hc.h     = cosmo.h;         % dimensionless
    hc.Om0   = cosmo.Omega_m;   % z=0 matter fraction
    hc.OL0   = cosmo.Omega_L;   % z=0 Lambda fraction

    % --- z-dependent ---
    hc.z        = z;
    hc.Ez       = Ez;
    hc.Hz       = Hz;           % km/s/kpc
    hc.Hz_Mpc   = Hz_Mpc;       % km/s/Mpc
    hc.Om_z     = Om_z;
    hc.rho_crit = rho_crit;     % M_sun/kpc^3

    % --- Bryan & Norman quantities ---
    hc.Delta_BN = Delta_BN;                          % B&N 1998 Eq.(6)
    hc.Delta    = Delta;                             % chosen overdensity
    hc.rho_vir  = Delta * rho_crit;                  % M_sun/kpc^3
    hc.zeta     = Delta_BN * Om_z / (18*pi^2);       % B&N 1998 Eq.(A2)

    % --- Physical constants ---
    hc.G       = G_kpc;
    hc.Gyr_s   = Gyr_s;
    hc.km_kpc  = km_kpc;

    % --- Hubble time [Gyr] ---
    %   t_H = 1/H(z);  Hz [km/s/kpc] * km_kpc [km/kpc] -> [s^-1]
    %   1/s / Gyr_s -> Gyr
    hc.t_H = 1.0 / (Hz * km_kpc) / Gyr_s;
end

% ==========================================================================
% CONCENTRATION  —  delegates to c_CDM (PHANTOM dispatcher)
% ==========================================================================
function c = get_concentration(M_h, hc, cosmo, z, opts)
% get_concentration  Resolve concentration from opts, calling c_CDM.
%
%   opts.concentration behaviour:
%     numeric (not NaN) — use that fixed value; skip all models
%     string            — pass as model name to c_CDM
%     '' | NaN | absent — default to 'ishiyama21'
%
%   c_CDM expects M in [Msun/h]; all other masses in PHANTOM are [Msun].
%   Concentration models available via c_CDM:
%     'bullock01','duffy08','klypin11','prada12','dutton14',
%     'diemer15','ludlow16','klypin16','child18','diemer19','ishiyama21'
%   Reference: c_CDM.m (PHANTOM src/concentration/c_CDM.m)

    oc = opts.concentration;

    % Fixed numeric value — bypass model entirely
    if isnumeric(oc) && isscalar(oc) && ~isnan(oc)
        c = oc;
        return
    end

    % Resolve model name
    if ischar(oc) && ~isempty(oc)
        model = oc;
    else
        model = 'ishiyama21';   % PHANTOM default
    end

    % c_CDM expects M in [Msun/h]  —  divide by h
    if ~isempty(opts.mode)
        c = c_CDM(M_h / hc.h, z, model, cosmo, opts.mode);
    else
        c = c_CDM(M_h / hc.h, z, model, cosmo);
    end
    c = max(c, 2.0);   % numerical floor
end

% ==========================================================================
% CONVERT INPUT -> M_h  [M_sun]
% ==========================================================================
function M_h = to_M_h(qin, val, hc, cosmo, z, opts)
% Routes input quantity to M_h [M_sun].
%
%   EQUATIONS
%   r_vir -> M_h :  M = (4*pi/3)*Delta*rho_crit*r_vir^3
%                   Bryan & Norman 1998; Robles+2018 Sec.2.1
%   V_vir -> M_h :  V_vir = sqrt(G*M/r_vir), r_vir = f*M^(1/3)
%                   => M = (V_vir^2 * f / G)^(3/2)
%                      f = (3/(4*pi*Delta*rho_crit))^(1/3)
%   V_max -> M_h :  V_max = V_vir*sqrt(0.216*c/f(c)); Bullock+2001 Eq.(3)
%   sigma_v->M_h :  sigma_v = V_vir/sqrt(2)
%   T_vir -> M_h :  T_vir = mu*m_p*V_vir^2/(2*k_B); Voit 2005 Eq.(2)

    G  = hc.G;
    D  = hc.Delta;
    rc = hc.rho_crit;
    f  = (3 / (4*pi * D * rc))^(1/3);   % r_vir = f * M^(1/3)

    switch lower(qin)
        case 'm_h'
            M_h = val;

        case 'r_vir'
            % M = (4*pi/3)*Delta*rho_crit*r_vir^3   [B&N 1998]
            M_h = (4*pi/3) * D * rc * val^3;

        case 'v_vir'
            % V_vir = sqrt(G*M/r_vir),  r_vir = f*M^(1/3)
            % => V^2 = G*M^(2/3)/f  =>  M = (V^2*f/G)^(3/2)
            M_h = (val^2 * f / G)^(3/2);

        case 'v_max'
            % V_max = V_vir * sqrt(0.216*c/fc)  Bullock+2001 Eq.(3)
            % Need a seed M_h to get c; iterate once (c is weak fn of M)
            M_seed = (val^2 * f / G)^(3/2);   % rough seed via V~V_vir
            c   = get_concentration(M_seed, hc, cosmo, z, opts);
            fc  = log(1+c) - c/(1+c);
            Vvir = val / sqrt(0.216 * c / fc);
            M_h  = (Vvir^2 * f / G)^(3/2);

        case 'sigma_v'
            % sigma_v = V_vir / sqrt(2)
            Vvir = val * sqrt(2);
            M_h  = (Vvir^2 * f / G)^(3/2);

        case 't_vir'
            % T_vir = mu*m_p*V_vir^2/(2*k_B)   Voit 2005 Eq.(2)
            % mu = 0.59 (fully ionised primordial H+He, X=0.76)
            mu  = 0.59;
            m_p = 1.6726e-27;   % kg
            k_B = 1.3806e-23;   % J/K
            % val in K; V_vir in km/s => convert to m/s (*1e3)
            Vvir = sqrt(2 * k_B * val / (mu * m_p)) / 1e3;   % km/s
            M_h  = (Vvir^2 * f / G)^(3/2);

        case 'rho_mean'
            error('halo_obs: rho_mean is the same for all halos at given z; provide M_h or r_vir.');

        case 'concentration'
            error('halo_obs: concentration alone does not determine M_h; pair it with M_h via opts.concentration.');

        otherwise
            error('halo_obs: unsupported input quantity ''%s''.', qin);
    end
end

% ==========================================================================
% CONVERT M_h -> OUTPUT QUANTITY
% ==========================================================================
function val = from_M_h(qout, M_h, hc, cosmo, z, opts)
% Computes any derived halo quantity from M_h [M_sun].
%
%   All equations cited inline.

    G  = hc.G;
    D  = hc.Delta;
    rc = hc.rho_crit;

    % ------------------------------------------------------------------
    % Primary geometric quantities
    % ------------------------------------------------------------------
    % r_vir [kpc]   Bryan & Norman 1998; Robles+2018 Sec.2.1
    %   M = (4*pi/3)*Delta*rho_crit*r_vir^3
    r_vir = (3 * M_h / (4*pi * D * rc))^(1/3);

    % V_vir [km/s]
    V_vir = sqrt(G * M_h / r_vir);

    % ------------------------------------------------------------------
    % Concentration  —  c_CDM dispatcher (PHANTOM src/concentration/)
    % ------------------------------------------------------------------
    c   = get_concentration(M_h, hc, cosmo, z, opts);
    r_s = r_vir / c;
    fc  = log(1+c) - c/(1+c);

    % NFW scale density  NFW 1997, ApJ 490, 493, Eq.(7)
    %   M = 4*pi*rho_s*r_s^3*f(c)
    rho_s = M_h / (4*pi * r_s^3 * fc);

    % ------------------------------------------------------------------
    % Velocity quantities
    % ------------------------------------------------------------------
    % V_max (NFW)   Bullock+2001, MNRAS 321, 559, Eq.(3)
    %   V_max = V_vir * sqrt(0.216*c/f(c))
    V_max   = V_vir * sqrt(0.216 * c / fc);

    % 1-D velocity dispersion   isotropic isothermal; sigma^2 = V_vir^2/2
    sigma_v = V_vir / sqrt(2);

    % ------------------------------------------------------------------
    % Energetics  (virial theorem for a virialized halo)
    %   E_pot = -G*M_h^2 / r_vir
    %   E_kin = -E_pot/2                     virial theorem
    %   E_tot = E_kin + E_pot = E_pot/2
    %   Padmanabhan 1993; Voit 2005, Rev.Mod.Phys. 77, 207, Eq.(2)
    % ------------------------------------------------------------------
    E_pot = -G * M_h^2 / r_vir;
    E_kin = -E_pot / 2;
    E_tot =  E_kin + E_pot;

    % ------------------------------------------------------------------
    % Virial temperature   Voit 2005, Rev.Mod.Phys. 77, 207, Eq.(2)
    %   T_vir = mu*m_p*V_vir^2 / (2*k_B),   mu=0.59
    % ------------------------------------------------------------------
    mu    = 0.59;
    m_p   = 1.6726e-27;   % kg
    k_B   = 1.3806e-23;   % J/K
    T_vir = mu * m_p * (V_vir * 1e3)^2 / (2 * k_B);   % V in m/s

    % ------------------------------------------------------------------
    % Dynamical time   Padmanabhan 1993; Bryan & Norman 1998, Eq.(2)
    %   t_dyn = sqrt(3*pi / (16*G*rho_mean))
    %   G [kpc (km/s)^2/M_sun], rho_mean [M_sun/kpc^3]
    %   => G*rho [(km/s)^2/kpc^2]
    %   => sqrt(1/(G*rho)) [kpc/(km/s)]
    %   conversion: 1 kpc/(km/s) = km_kpc [km/kpc] / (km/s) = km_kpc s
    %               / Gyr_s -> Gyr
    % ------------------------------------------------------------------
    rho_mean  = D * rc;
    t_dyn_raw = sqrt(3*pi / (16 * G * rho_mean));   % kpc / (km/s)
    t_dyn     = t_dyn_raw * hc.km_kpc / hc.Gyr_s;  % Gyr

    % ------------------------------------------------------------------
    % r_200c and M_200c  (from NFW profile, independent of Delta choice)
    %   Find r_200c: mean enclosed density = 200*rho_crit
    %   M_NFW(<r) = 4*pi*rho_s*r_s^3 * [ln(1+r/r_s) - (r/r_s)/(1+r/r_s)]
    %   NFW 1997, ApJ 490, 493, Eq.(7)
    % ------------------------------------------------------------------
    r200c_val = find_r_delta(rho_s, r_s, 200, rc);
    M200c_val = nfw_mass(rho_s, r_s, r200c_val);

    % ------------------------------------------------------------------
    % Peak height   nu = delta_c / sigma(M, z)
    %   delta_c ~ 1.686  (spherical collapse in EdS)
    %   sigma(M, z) from cosmo.sigmaM(M, z)  [PHANTOM handle]
    %   cosmo.sigmaM expects M in [Msun/h] — divide by h
    %   Press & Schechter 1974; Lacey & Cole 1994, MNRAS 271, 676
    % ------------------------------------------------------------------
    delta_c = 1.686;
    sigma_M = cosmo.sigmaM(M_h / hc.h, z);
    nu_val  = delta_c / sigma_M;

    % ------------------------------------------------------------------
    % Route to requested output
    % ------------------------------------------------------------------
    switch lower(qout)
        case 'm_h',           val = M_h;
        case 'r_vir',         val = r_vir;
        case 'v_vir',         val = V_vir;
        case 'v_max',         val = V_max;
        case 'sigma_v',       val = sigma_v;
        case 't_vir',         val = T_vir;
        case 'e_kin',         val = E_kin;
        case 'e_pot',         val = E_pot;
        case 'e_tot',         val = E_tot;
        case 't_dyn',         val = t_dyn;
        case 'rho_mean',      val = rho_mean;
        case 'rho_s',         val = rho_s;
        case 'r_s',           val = r_s;
        case 'concentration', val = c;
        case 'r_200c',        val = r200c_val;
        case 'm_200c',        val = M200c_val;
        case 'nu',            val = nu_val;
        otherwise
            error('halo_obs: unsupported output quantity ''%s''.', qout);
    end
end

% ==========================================================================
% PURELY COSMOLOGICAL QUANTITIES  (no M_h needed)
% ==========================================================================
function val = from_cosmo(qout, hc)
    switch lower(qout)
        case 'delta_vir',  val = hc.Delta_BN;
        case 'omega_m_z',  val = hc.Om_z;
        case 'zeta_z',     val = hc.zeta;
        case 'rho_c_z',    val = hc.rho_crit;
        otherwise
            error('halo_obs: ''%s'' is not a purely cosmological quantity.', qout);
    end
end

% ==========================================================================
% NFW HELPER FUNCTIONS
% ==========================================================================
function M = nfw_mass(rho_s, r_s, r)
% M_NFW(<r) = 4*pi*rho_s*r_s^3 * [ln(1+x) - x/(1+x)],  x = r/r_s
% NFW 1997, ApJ 490, 493, Eq.(7)
    x = r / r_s;
    M = 4*pi * rho_s * r_s^3 * (log(1+x) - x/(1+x));
end

function r_delta = find_r_delta(rho_s, r_s, Delta, rho_crit)
% Find r_delta: mean enclosed density = Delta*rho_crit.
% Bisection on f(r) = M_NFW(<r)/((4*pi/3)*r^3) - Delta*rho_crit = 0
    target = Delta * rho_crit;
    rlo    = 1e-3 * r_s;
    rhi    = 1e4  * r_s;
    for iter = 1:100
        rmid         = 0.5 * (rlo + rhi);
        rho_mean_mid = nfw_mass(rho_s, r_s, rmid) / ((4*pi/3) * rmid^3);
        if rho_mean_mid > target
            rlo = rmid;
        else
            rhi = rmid;
        end
        if (rhi - rlo) / rmid < 1e-8, break; end
    end
    r_delta = 0.5 * (rlo + rhi);
end
