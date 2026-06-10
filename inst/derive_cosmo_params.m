function cosmo = derive_cosmo_params(cosmo)
% DERIVE_COSMO_PARAMS  Compute derived cosmological parameters and attach
% background expansion, density, time, and distance functions.
%
%   COSMO = DERIVE_COSMO_PARAMS(COSMO) fills in derived z=0 quantities and
%   attaches function handles for background cosmology calculations in a
%   Friedman-Lemaitre-Robertson-Walker model, following the same general
%   structure used in the COLOSSUS cosmology module.
%
%   Required input fields
%   ---------------------
%   cosmo.h          - Dimensionless Hubble parameter, H0 / 100
%   cosmo.Omega_m    - Matter density parameter at z = 0
%   cosmo.Omega_b    - Baryon density parameter at z = 0
%   cosmo.rho_crit0  - Critical density at z = 0 in the code's internal units
%
%   Optional input fields
%   ---------------------
%   cosmo.Tcmb       - CMB temperature in K (default: 2.7255)
%   cosmo.Neff       - Effective number of neutrino species (default: 3.046)
%   cosmo.relspecies - Logical flag to include relativistic species
%                      (default: false)
%   cosmo.flat       - Logical flag for flat cosmology (default: true)
%   cosmo.Omega_L    - Dark-energy density parameter at z = 0
%                      If flat = true and Omega_L is not supplied, it is set to
%                      1 - Omega_m - Omega_r.
%                      If flat = false, Omega_L must be supplied.
%   cosmo.de_model   - Dark-energy model:
%                      'lambda' or 'lcdm' : cosmological constant
%                      'w0' or 'wcdm'     : constant w dark energy
%                      'w0wa' or 'cpl'    : CPL model w(a) = w0 + wa(1-a)
%                      (default: 'lambda')
%   cosmo.w0         - Constant or present-day dark-energy equation-of-state
%                      parameter (default: -1)
%   cosmo.wa         - CPL evolution parameter (default: 0)
%   cosmo.zmax       - Upper redshift bound used for age integrations
%                      (default: 1e4)
%
%   Derived scalar fields added to COSMO
%   ------------------------------------
%   H0, h2, Omh2, Ombh2, Omega_c, Omch2
%   Omega_gamma, Omega_nu, Omega_r
%   Omega_k, a_eq, z_eq
%   w0, wa
%   age0
%
%   Function handles added to COSMO
%   -------------------------------
%   Expansion and densities:
%     cosmo.fde(z)                  - Dark-energy evolution factor
%     cosmo.E(z)                    - Dimensionless Hubble parameter H(z)/H0
%     cosmo.Hz(z)                   - Hubble parameter in km/s/Mpc
%     cosmo.rhocrit(z)              - Critical density
%     cosmo.rhom(z)                 - Matter density
%     cosmo.rhob(z)                 - Baryon density
%     cosmo.rhocdm(z)               - Cold dark matter density
%     cosmo.rhor(z)                 - Radiation density
%     cosmo.rhoL(z)                 - Dark-energy density term
%
%   Redshift-dependent density parameters:
%     cosmo.Omega_m_z(z)
%     cosmo.Omega_b_z(z)
%     cosmo.Omega_c_z(z)
%     cosmo.Omega_r_z(z)
%     cosmo.Omega_L_z(z)
%     cosmo.Omega_gamma_z(z)        - Attached only if relativistic species are used
%     cosmo.Omega_nu_z(z)           - Attached only if relativistic species are used
%     cosmo.Omega_k_z(z)            - Attached only if Omega_k ~= 0
%
%   Time functions:
%     cosmo.time(z)                 - Struct with fields:
%                                     .lookback_Gyr
%                                     .t0_Gyr
%                                     .age_Gyr
%     cosmo.lookbackTime(z)         - Lookback time in Gyr
%     cosmo.age(z)                  - Age of the universe at redshift z in Gyr
%
%   Distance functions (all returned in comoving Mpc/h):
%     cosmo.comovingDistance(z)             - Line-of-sight comoving distance
%     cosmo.transverseComovingDistance(z)   - Transverse comoving distance
%     cosmo.angularDiameterDistance(z)      - Angular-diameter distance
%     cosmo.luminosityDistance(z)           - Luminosity distance
%
%   Notes
%   -----
%   1. In this implementation, Omega_L is the only public dark-energy
%      density parameter at z = 0. Even for non-lambda dark energy models,
%      the redshift evolution is handled through cosmo.fde(z), not through
%      a separate Omega_de field.
%
%   2. For de_model = 'lambda', the dark-energy term is constant and
%      cosmo.fde(z) = 1.
%
%   3. For de_model = 'w0', the dark-energy scaling is
%         fde(z) = (1 + z)^(3 * (1 + w0))
%
%   4. For de_model = 'w0wa' or 'cpl', the CPL form is used:
%         fde(z) = (1 + z)^(3 * (1 + w0 + wa)) * exp(-3 * wa * z / (1 + z))
%
%   5. Distances are returned in comoving Mpc/h, consistent with the
%      conventions commonly used in large-scale structure and halo studies.
%
%   6. Times and distances are evaluated by direct numerical integration
%      each time they are called. This is simple and robust, but slower than
%      interpolation-based implementations such as COLOSSUS.
%
%   Example
%   -------
%   cosmo = derive_cosmo_params(cosmo);
%   Ez  = cosmo.E(1.0);
%   t0  = cosmo.age(0);
%   tL  = cosmo.lookbackTime(2.0);
%   dc  = cosmo.comovingDistance(1.0);
%   da  = cosmo.angularDiameterDistance(1.0);
%   dl  = cosmo.luminosityDistance(1.0);
%   Omz = cosmo.Omega_m_z(3.0);
%
%   See also: ATTACH_LINEAR_COMPONENTS
    
    % ---------------------------------------------------------------
    % Derived z=0 parameters
    % ---------------------------------------------------------------
    cosmo.H0    = 100 * cosmo.h;
    cosmo.h2    = cosmo.h^2;
    cosmo.Omh2  = cosmo.Omega_m * cosmo.h2;
    cosmo.Ombh2 = cosmo.Omega_b * cosmo.h2;

    cosmo.Omega_c = cosmo.Omega_m - cosmo.Omega_b;
    cosmo.Omch2   = cosmo.Omega_c * cosmo.h2;

    % ---------------------------------------------------------------
    % Defaults
    % ---------------------------------------------------------------
    if ~isfield(cosmo, 'Tcmb') || isempty(cosmo.Tcmb)
        cosmo.Tcmb = 2.7255;
    end
    if ~isfield(cosmo, 'Neff') || isempty(cosmo.Neff)
        cosmo.Neff = 3.046;
    end
    if ~isfield(cosmo, 'relspecies') || isempty(cosmo.relspecies)
        cosmo.relspecies = false;
    end
    if ~isfield(cosmo, 'flat') || isempty(cosmo.flat)
        cosmo.flat = true;
    end
    if ~isfield(cosmo, 'zmax') || isempty(cosmo.zmax)
        cosmo.zmax = 1e4;
    end
    if ~isfield(cosmo, 'de_model') || isempty(cosmo.de_model)
        cosmo.de_model = 'lambda';
    end

    % ---------------------------------------------------------------
    % Radiation: photons + massless neutrinos
    % ---------------------------------------------------------------
    if cosmo.relspecies
        cosmo.Omega_gamma = 4.4814665e-7 * cosmo.Tcmb^4 / cosmo.h2;
        cosmo.Omega_nu    = 0.22710731766 * cosmo.Neff * cosmo.Omega_gamma;
        cosmo.Omega_r     = cosmo.Omega_gamma + cosmo.Omega_nu;

        cosmo.a_eq = cosmo.Omega_r / cosmo.Omega_m;
        cosmo.z_eq = 1.0 / cosmo.a_eq - 1.0;
    else
        cosmo.Omega_gamma = 0.0;
        cosmo.Omega_nu    = 0.0;
        cosmo.Omega_r     = 0.0;
    end

    % ---------------------------------------------------------------
    % Curvature and dark energy at z=0
    % Keep only Omega_L as the dark-energy density parameter today
    % ---------------------------------------------------------------
    if cosmo.flat
        cosmo.Omega_k = 0.0;

        if ~isfield(cosmo, 'Omega_L') || isempty(cosmo.Omega_L)
            cosmo.Omega_L = 1.0 - cosmo.Omega_m - cosmo.Omega_r;
        end
    else
        if ~isfield(cosmo, 'Omega_L') || isempty(cosmo.Omega_L)
            error('For non-flat cosmology, please provide cosmo.Omega_L.');
        end
        cosmo.Omega_k = 1.0 - cosmo.Omega_m - cosmo.Omega_L - cosmo.Omega_r;
    end

    % ---------------------------------------------------------------
    % Dark-energy evolution factor fde(z)
    % E(z)^2 = Om(1+z)^3 + Or(1+z)^4 + Ok(1+z)^2 + OL * fde(z)
    % ---------------------------------------------------------------
    switch lower(cosmo.de_model)
        case {'lambda', 'lcdm'}
            cosmo.w0  = -1.0;
            cosmo.wa  = 0.0;
            cosmo.fde = @(z) ones(size(z));

        case {'w0', 'wcdm'}
            if ~isfield(cosmo, 'w0') || isempty(cosmo.w0)
                cosmo.w0 = -1.0;
            end
            cosmo.wa  = 0.0;
            cosmo.fde = @(z) (1 + z).^(3 * (1 + cosmo.w0));

        case {'w0wa', 'cpl'}
            if ~isfield(cosmo, 'w0') || isempty(cosmo.w0)
                cosmo.w0 = -1.0;
            end
            if ~isfield(cosmo, 'wa') || isempty(cosmo.wa)
                cosmo.wa = 0.0;
            end
            cosmo.fde = @(z) (1 + z).^(3 * (1 + cosmo.w0 + cosmo.wa)) .* ...
                             exp(-3 * cosmo.wa .* z ./ (1 + z));

        otherwise
            error('Unknown dark energy model "%s".', cosmo.de_model);
    end

    % ---------------------------------------------------------------
    % E(z)^2 and H(z)
    % ---------------------------------------------------------------
    E2_terms = @(z) ...
        cosmo.Omega_m .* (1 + z).^3 + ...
        cosmo.Omega_r .* (1 + z).^4 + ...
        cosmo.Omega_k .* (1 + z).^2 + ...
        cosmo.Omega_L .* cosmo.fde(z);

    cosmo.E  = @(z) sqrt(E2_terms(z));
    cosmo.Hz = @(z) cosmo.H0 .* cosmo.E(z);

    % ---------------------------------------------------------------
    % Critical and component densities
    % rho_crit0 is assumed already defined elsewhere in your code
    % ---------------------------------------------------------------
    cosmo.rhocrit = @(z) cosmo.rho_crit0 .* E2_terms(z);
    cosmo.rhom    = @(z) cosmo.rho_crit0 .* cosmo.Omega_m .* (1 + z).^3;
    cosmo.rhob    = @(z) cosmo.rho_crit0 .* cosmo.Omega_b .* (1 + z).^3;
    cosmo.rhocdm  = @(z) cosmo.rho_crit0 .* cosmo.Omega_c .* (1 + z).^3;
    cosmo.rhor    = @(z) cosmo.rho_crit0 .* cosmo.Omega_r .* (1 + z).^4;
    cosmo.rhoL    = @(z) cosmo.rho_crit0 .* cosmo.Omega_L .* cosmo.fde(z);

    % ---------------------------------------------------------------
    % Redshift-dependent density parameters
    % ---------------------------------------------------------------
    cosmo.Omega_m_z = @(z) cosmo.Omega_m .* (1 + z).^3 ./ E2_terms(z);
    cosmo.Omega_b_z = @(z) cosmo.Omega_b .* (1 + z).^3 ./ E2_terms(z);
    cosmo.Omega_c_z = @(z) cosmo.Omega_c .* (1 + z).^3 ./ E2_terms(z);
    cosmo.Omega_r_z = @(z) cosmo.Omega_r .* (1 + z).^4 ./ E2_terms(z);
    cosmo.Omega_L_z = @(z) cosmo.Omega_L .* cosmo.fde(z) ./ E2_terms(z);

    if cosmo.Omega_gamma > 0
        cosmo.Omega_gamma_z = @(z) cosmo.Omega_gamma .* (1 + z).^4 ./ E2_terms(z);
        cosmo.Omega_nu_z    = @(z) cosmo.Omega_nu    .* (1 + z).^4 ./ E2_terms(z);
    end
    if abs(cosmo.Omega_k) > 0
        cosmo.Omega_k_z = @(z) cosmo.Omega_k .* (1 + z).^2 ./ E2_terms(z);
    end

    % ---------------------------------------------------------------
    % Time functions
    % ---------------------------------------------------------------
    cosmo.time = @(z) z_to_time_gyr(z, cosmo);
    cosmo.lookbackTime = @(z) z_to_time_gyr(z, cosmo).lookback_Gyr;
    cosmo.age = @(z) z_to_time_gyr(z, cosmo).age_Gyr;
    cosmo.age0 = cosmo.age(0);

    % ---------------------------------------------------------------
    % Distance functions (returned in Mpc/h)
    % ---------------------------------------------------------------
    cosmo.comovingDistance = @(z) comoving_distance(z, cosmo);
    cosmo.transverseComovingDistance = @(z) transverse_comoving_distance(z, cosmo);
    cosmo.angularDiameterDistance = @(z) transverse_comoving_distance(z, cosmo) ./ (1 + z);
    cosmo.luminosityDistance      = @(z) transverse_comoving_distance(z, cosmo) .* (1 + z);

end








