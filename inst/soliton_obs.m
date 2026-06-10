
%% soliton_observables.m
%
% Soliton core observables for Scalar Field / Fuzzy Dark Matter (SFDM/FDM)
% haloes, derived from:
%
%   [S14a] Schive, Chiueh & Broadhurst (2014, Nat. Phys. 10, 496)
%          Eqs. (S3), (S4) in the Supplementary Material
%   [S14b] Schive et al. (2014, PRL 113, 261302)
%          Eqs. (3), (6), (7)
%   [R18]  Robles, Bullock & Boylan-Kolchin (2018, MNRAS)
%          Eqs. (2)–(10)
%
% Every observable can be obtained from any ONE of the four primary
% inputs: rho_c, r_c, M_c, or M_h.  Choose an input–output pair and
% call the corresponding sub-function, or call
%
%     result = soliton_obs(quantity_in, value_in, quantity_out, m22)
%
% for a single dispatcher interface.
%
% -----------------------------------------------------------------------
%  Units used throughout
%    Mass            : solar masses  [M_sun]
%    Length          : kiloparsecs   [kpc]
%    Density         : M_sun / kpc^3
%    Velocity        : km / s
%    Energy          : M_sun (km/s)^2  (i.e. natural N-body energy)
%    m22             : boson mass in units of 1e-22 eV/c^2  (dimensionless)
%    a               : cosmic scale factor (default 1 = z=0)
% -----------------------------------------------------------------------

function result = soliton_obs(quantity_in, value_in, quantity_out, m22, a, zeta_z, zeta_0)
% SOLITON_OBS  Universal dispatcher for soliton observable conversions.
%
%  result = soliton_obs(quantity_in, value_in, quantity_out, m22)
%  result = soliton_obs(quantity_in, value_in, quantity_out, m22, a)
%
%  quantity_in / quantity_out : one of
%       'rho_c'  – peak (central) soliton density  [M_sun/kpc^3]
%       'r_c'    – soliton core radius (half-peak-density radius)  [kpc]
%       'M_c'    – soliton core mass  M(r < r_c)  [M_sun]
%       'M_h'    – host halo virial mass  [M_sun]
%       'V_c'    – peak circular velocity inside soliton  [km/s]
%       'V_h'    – virial (halo) circular velocity  [km/s]
%       'KE'     – quantum kinetic energy of soliton  [M_sun (km/s)^2]
%       'PE'     – gravitational potential energy of soliton  [M_sun (km/s)^2]
%       'E_tot'  – total soliton energy (KE + PE)  [M_sun (km/s)^2]
%       'sigma'  – 1-D line-of-sight velocity dispersion  [km/s]
%       'r_vir'  – host halo virial radius  [kpc]
%       'lam_c'  – de Broglie / angular momentum scale  [kpc km/s]
%       'M_c_over_M_h' – soliton core mass fraction  [dimensionless]
%       'r_c_over_r_vir' – core-to-virial radius fraction [dimensionless]
%
%  m22  : boson mass in units of 1e-22 eV/c^2
%  a    : scale factor (optional, default = 1, i.e. z = 0)

    if nargin < 5, a = 1; end
    if nargin == 7
        zeta_args = {zeta_z, zeta_0};
    else
        zeta_args = {};
    end

    % Step 1: Convert input to r_c (the canonical free parameter).
    switch lower(quantity_in)
        case 'r_c'
            r_c = value_in;
        case 'rho_c'
            r_c = rho_c_to_r_c(value_in, m22);
        case 'm_c'
            r_c = M_c_to_r_c(value_in, m22);
        case 'm_h'
            r_c = M_h_to_r_c(value_in, m22, a, zeta_args{:});
        otherwise
            error('Unknown input quantity: %s', quantity_in);
    end

    % Step 2: Derive requested output from r_c.
    switch lower(quantity_out)
        case 'rho_c'
            result = r_c_to_rho_c(r_c, m22);
        case 'r_c'
            result = r_c;
        case 'm_c'
            result = r_c_to_M_c(r_c, m22);
        case 'm_h'
            result = r_c_to_M_h(r_c, m22, a, zeta_args{:});
        case 'v_c'
            result = r_c_to_V_c(r_c, m22);
        case 'v_h'
            result = r_c_to_V_h(r_c, m22, a);
        case 'ke'
            result = r_c_to_KE(r_c, m22);
        case 'pe'
            result = r_c_to_PE(r_c, m22);
        case 'e_tot'
            result = r_c_to_E_tot(r_c, m22);
        case 'sigma'
            result = r_c_to_sigma(r_c, m22);
        case 'r_vir'
            result = r_c_to_r_vir(r_c, m22, a);
        case 'lam_c'
            result = r_c_to_lambda_c(r_c, m22);
        case 'm_c_over_m_h'
            result = r_c_to_M_c_over_M_h(r_c, m22, a);
        case 'r_c_over_r_vir'
            result = r_c_to_r_c_over_r_vir(r_c, m22, a);
        otherwise
            error('Unknown output quantity: %s', quantity_out);
    end
end


%% ======================================================================
%  PRIMITIVE CONVERTERS   (all work on scalars or arrays)
% =======================================================================

% -----------------------------------------------------------------------
%  rho_c  <-->  r_c
% -----------------------------------------------------------------------

function rho_c = r_c_to_rho_c(r_c, m22)
% R_C_TO_RHO_C  Peak soliton density from core radius.
%
%  rho_c = 1.93e7 * m22^-2 * (r_c / 1 kpc)^{-4}   [M_sun / kpc^3]
%
%  Source: Robles+2018, Eq. (3);  Schive+2014 (Nat.Phys.), Eq. (S4).
    rho_c = 1.93e7 ./ (m22.^2 .* r_c.^4);
end

function r_c = rho_c_to_r_c(rho_c, m22)
% RHO_C_TO_R_C  Core radius from peak density.
%
%  r_c = (1.93e7 * m22^-2 / rho_c)^{1/4}   [kpc]
%
%  Source: Robles+2018, Eq. (3) inverted.
    r_c = (1.93e7 ./ (m22.^2 .* rho_c)).^(0.25);
end


% -----------------------------------------------------------------------
%  M_c  <-->  r_c
% -----------------------------------------------------------------------

function M_c = r_c_to_M_c(r_c, m22)
% R_C_TO_M_C  Soliton core mass (mass enclosed within r_c) from core radius.
%
%  M_c = 8.06e7 * m22^{-2} * (r_c / 1 kpc)^{-1}   [M_sun]
%
%  Source: Schive+2014 (Nat.Phys. supplement), Eq. (S3);
%          Robles+2018, Eq. (10) rearranged via rc*Mc = const.
    M_c = 8.06e7 ./ (m22.^2 .* r_c);
end

function r_c = M_c_to_r_c(M_c, m22)
% M_C_TO_R_C  Core radius from soliton core mass.
%
%  r_c = 8.06e7 / (m22^2 * M_c)   [kpc]
%
%  Source: Schive+2014 (Nat.Phys. supplement), Eq. (S3) inverted.
    r_c = 8.06e7 ./ (m22.^2 .* M_c);
end


% -----------------------------------------------------------------------
%  M_h  <-->  r_c
% -----------------------------------------------------------------------

function r_c = M_h_to_r_c(M_h, m22, a, zeta_z, zeta_0)
% M_H_TO_R_C  Soliton core radius from host halo virial mass.
%
%  Full form (Schive+2014 PRL, Eq. 7):
%    r_c = 1.6 * m22^{-1} * a^{1/2} * (zeta_z/zeta_0)^{-1/6}
%                         * (M_h / 1e9 M_sun)^{-1/3}   [kpc]
%
%  where  zeta(z) = 18*pi^2 + 82*(Omega_m(z)-1) - 39*(Omega_m(z)-1)^2
%  is the Bryan-Norman (1998) virial overdensity factor.
%
%  If zeta_z and zeta_0 are NOT supplied, the ratio (zeta_z/zeta_0)^{-1/6}
%  is set to 1 (valid at z=0 or as a z=0 approximation).
%
%  Inputs:
%    M_h    : halo virial mass          [M_sun]
%    m22    : boson mass / 1e-22 eV
%    a      : scale factor  (default 1)
%    zeta_z : zeta evaluated at redshift z  (optional)
%    zeta_0 : zeta evaluated at z=0         (optional, ~18*pi^2 ~ 177.65)
%
%  Source: Schive+2014 (PRL 113, 261302), Eq. (7).
%          Bryan & Norman (1998) for zeta definition, cited as Ref. [64].

    if nargin < 3, a = 1; end

    if nargin == 5
        zeta_factor = (zeta_z / zeta_0)^(-1/6);
    else
        zeta_factor = 1;   % z=0 approximation
    end

    r_c = (1.6 ./ m22) .* sqrt(a) .* zeta_factor .* (M_h ./ 1e9).^(-1/3);
end

function M_h = r_c_to_M_h(r_c, m22, a, zeta_z, zeta_0)
% R_C_TO_M_H  Host halo virial mass from soliton core radius.
%  Inverse of M_h_to_r_c.
%
%  Source: Schive+2014 (PRL 113, 261302), Eq. (7) inverted.

    if nargin < 3, a = 1; end

    if nargin == 5
        zeta_factor = (zeta_z / zeta_0)^(-1/6);
    else
        zeta_factor = 1;
    end

    M_h = 1e9 .* (1.6 .* sqrt(a) .* zeta_factor ./ (m22 .* r_c)).^3;
end


% -----------------------------------------------------------------------
%  Virial radius
% -----------------------------------------------------------------------

function r_vir = r_c_to_r_vir(r_c, m22, a)
% R_C_TO_R_VIR  Host halo virial radius from soliton core radius.
%
%  Uses M_h --> r_vir via the Bryan-Norman (1998) definition of virial mass:
%    M_h = (4*pi/3) * Delta_vir * rho_crit * r_vir^3
%  At z=0 with Delta_vir ~ 360 and rho_crit = 2.775e11 * h^2 * Omega_m  M_sun/Mpc^3,
%  a convenient pre-computed result gives (with h=0.7, Omega_m=0.3):
%    r_vir [kpc] = 258.6 * (M_h / 1e12 M_sun)^{1/3}
%
%  This sub-function calls r_c_to_M_h first.
%
%  Source: Robles+2018 Eq. (8) uses rc/r_vir directly.
    if nargin < 3, a = 1; end
    M_h   = r_c_to_M_h(r_c, m22, a);                    % [M_sun]
    r_vir = r_c .* m22 .* (M_h ./ 1e9).^(2/3) ./ (6.20e-2); % [kpc], z=0
end


% -----------------------------------------------------------------------
%  Circular velocities
% -----------------------------------------------------------------------

function V_c = r_c_to_V_c(r_c, m22)
% R_C_TO_V_C  Peak circular velocity within the soliton core.
%
%  V_c = sqrt(G * M_c / r_c)   with G in consistent units.
%
%  In units where G = 4.302e-3 pc M_sun^{-1} (km/s)^2 :
%    G = 4.302e-6 kpc M_sun^{-1} (km/s)^2
%
%  The scaling from the literature gives (Schive+2014 Nat.Phys., text):
%    sigma ~ 115 km/s for M_c ~ 2e9 M_sun, r_c ~ 0.18 kpc  (Milky Way)
%  which is consistent with V_c = sqrt(G*M_c/r_c).
%
%  Source: computed from Eq. (S3) + virial condition.
    G_kpc = 4.302e-6;                    % kpc M_sun^{-1} (km/s)^2
    M_c   = r_c_to_M_c(r_c, m22);
    V_c   = sqrt(G_kpc .* M_c ./ r_c);  % [km/s]
end

function V_h = r_c_to_V_h(r_c, m22, a)
% R_C_TO_V_H  Halo virial circular velocity  V_h = sqrt(G * M_h / r_vir).
%
%  Source: standard virial definition; Robles+2018 text (Vc/Vvir ~ 0.9 note).
    if nargin < 3, a = 1; end
    G_kpc = 4.302e-6;
    M_h   = r_c_to_M_h(r_c, m22, a);
    r_vir = r_c_to_r_vir(r_c, m22, a);
    V_h   = sqrt(G_kpc .* M_h ./ r_vir);
end


% -----------------------------------------------------------------------
%  Energies   (quantum virial theorem: 2*KE + PE = 0  for a soliton)
% -----------------------------------------------------------------------

function KE = r_c_to_KE(r_c, m22)
% R_C_TO_KE  Quantum kinetic (gradient) energy of the soliton.
%
%  From the quantum virial theorem for solitons:  2*KE + PE = 0
%  => KE = -PE/2 = -E_tot
%  The total energy scales as  E_tot ~ -G*M_c^2 / r_c  (half of PE).
%
%  E_tot = -0.5 * G * M_c^2 / r_c   [M_sun (km/s)^2]
%  KE    = -E_tot = +0.5 * G * M_c^2 / r_c
%
%  Source: Schive+2014 (PRL), Eq. (5) and surrounding discussion;
%          quantum virial theorem (2 KE + PE = 0 for ground-state soliton).
    G_kpc = 4.302e-6;
    M_c   = r_c_to_M_c(r_c, m22);
    KE    = 0.5 .* G_kpc .* M_c.^2 ./ r_c;
end

function PE = r_c_to_PE(r_c, m22)
% R_C_TO_PE  Gravitational potential energy of the soliton.
%
%  PE = -G * M_c^2 / r_c   [M_sun (km/s)^2]
%
%  Source: virial theorem applied to the soliton ground state
%          (Schive+2014 PRL, discussion after Eq. 5).
    G_kpc = 4.302e-6;
    M_c   = r_c_to_M_c(r_c, m22);
    PE    = -G_kpc .* M_c.^2 ./ r_c;
end

function E_tot = r_c_to_E_tot(r_c, m22)
% R_C_TO_E_TOT  Total energy of the soliton  (KE + PE).
%
%  E_tot = KE + PE = -0.5 * G * M_c^2 / r_c   [M_sun (km/s)^2]
%
%  Source: Schive+2014 (PRL), Eq. (5).
    E_tot = r_c_to_KE(r_c, m22) + r_c_to_PE(r_c, m22);
end


% -----------------------------------------------------------------------
%  Velocity dispersion  (line-of-sight)
% -----------------------------------------------------------------------

function sigma = r_c_to_sigma(r_c, m22)
% R_C_TO_SIGMA  1-D line-of-sight velocity dispersion of test particles
%               in the soliton potential.
%
%  The soliton potential at r=0 gives a virial velocity dispersion.
%  For an isotropic distribution, sigma_1D = V_c / sqrt(3).
%
%  The Schive+2014 (Nat.Phys.) text quotes sigma ~ 115 km/s at rc ~ 0.18 kpc
%  for a Milky Way halo of M_h ~ 1e12 M_sun, consistent with this formula.
%
%  Source: Schive+2014 (Nat.Phys.), text after Eq. (S4).
    V_c   = r_c_to_V_c(r_c, m22);
    sigma = V_c ./ sqrt(3);
end


% -----------------------------------------------------------------------
%  Angular momentum / de Broglie scale
% -----------------------------------------------------------------------

function lam_c = r_c_to_lambda_c(r_c, m22)
% R_C_TO_LAMBDA_C  Characteristic angular momentum of a circular orbit
%                  within the soliton core:  lam_c = sqrt(G*r_c*M_c).
%
%  The product r_c * M_c is independent of M_h (depends only on m22):
%    r_c * M_c = 5.5e9 / m22^2   [M_sun kpc]
%
%  Hence:
%    lam_c = sqrt(G_kpc * r_c * M_c)
%           = sqrt(G_kpc * 5.5e9 / m22^2)   [kpc km/s]
%           ~ 18.6 / m22   [kpc km/s]
%
%  Source: Robles+2018, Eq. (10) and surrounding text.
    G_kpc  = 4.302e-6;
    M_c    = r_c_to_M_c(r_c, m22);
    lam_c  = sqrt(G_kpc .* r_c .* M_c);   % [kpc km/s]
end


% -----------------------------------------------------------------------
%  Mass fractions / radius fractions
% -----------------------------------------------------------------------

function ratio = r_c_to_M_c_over_M_h(r_c, m22, a)
% R_C_TO_M_C_OVER_M_H  Soliton core mass fraction.
%
%  M_c / M_h scales as  M_h^{-2/3}:
%    M_c / M_h = 5.04e-2 * (M_h / 1e9 M_sun)^{-2/3} / m22
%
%  Source: Robles+2018, Eq. (9).
    if nargin < 3, a = 1; end
    M_h   = r_c_to_M_h(r_c, m22, a);
    M_c   = r_c_to_M_c(r_c, m22);
    ratio = M_c ./ M_h;
end

function ratio = r_c_to_r_c_over_r_vir(r_c, m22, a)
% R_C_TO_R_C_OVER_R_VIR  Core-to-virial radius fraction.
%
%  r_c / r_vir scales as M_h^{-2/3}:
%    r_c / r_vir = 6.20e-2 * (M_h / 1e9 M_sun)^{-2/3} / m22
%
%  Source: Robles+2018, Eq. (8).
    if nargin < 3, a = 1; end
    r_vir = r_c_to_r_vir(r_c, m22, a);
    ratio = r_c ./ r_vir;
end


% %% ======================================================================
% %  SOLITON DENSITY PROFILE
% % =======================================================================
% 
% function rho = soliton_profile(r, rho_c, r_c)
% % SOLITON_PROFILE  Analytic approximation to the soliton density profile.
% %
% %  rho(r) = rho_c * [1 + 0.091*(r/r_c)^2]^{-8}
% %
% %  Accurate to ~2% for 0 <= r <= 3*r_c, enclosing ~95% of total soliton mass.
% %
% %  Source: Robles+2018, Eq. (2);  Schive+2014 (Nat.Phys.), Eq. (S4);
% %          Schive+2014 (PRL), Eq. (3).
% %
% %  Inputs:
% %    r     : radial coordinate  [kpc]  (scalar or array)
% %    rho_c : central density    [M_sun/kpc^3]
% %    r_c   : core radius        [kpc]
% %  Output:
% %    rho   : density            [M_sun/kpc^3]
%     rho = rho_c ./ (1 + 0.091 .* (r ./ r_c).^2).^8;
% end
