function c = Prada12(M200c, z, cosmo)
% Prada12_concentration  Prada et al. (2012) concentration model
%
%   c = Prada12(M200c, z, cosmo)
%
%   Implements Eqs. (12)-(22) of Prada et al. 2012, MNRAS 423, 3018.
%
%   Step-by-step (from paper Section 5):
%     1. Compute x and D(a)                         Eqs. (12-13)
%     2. Compute sigma(M, z) from cosmo             Eq.  (5) / cosmo.sigmaM
%     3. Compute cmin(x) and sigma_min_inv(x)       Eqs. (19-22)
%     4. Compute B0(x) and B1(x)                    Eq.  (18)
%     5. Rescale sigma: sigma' = B1(x) * sigma      Eq.  (15)
%     6. Evaluate C(sigma')                         Eqs. (16-17)
%     7. Final concentration: c = B0(x) * C(sigma') Eq.  (14)
%
%   INPUTS
%   M200c : halo mass M_200c [h^-1 Msun], scalar or vector
%   z     : redshift, scalar
%   cosmo : cosmology struct with fields:
%             cosmo.Omegam   matter density parameter (z=0)
%             cosmo.OmegaL   cosmological constant (z=0)
%             cosmo.sigmaM   function handle: sigma(M, z)
%
%   OUTPUT
%   c     : concentration c_200c, same shape as M200c
%
%   Reference: Prada, Klypin, Cuesta, Betancort-Rijo & Primack 2012,
%              MNRAS 423, 3018, Eqs. (12)-(22)

% ---- Eq. (13): x parameter ---------------------------------------------
a = 1 ./ (1 + z);
x = (cosmo.Omega_L ./ cosmo.Omega_m).^(1/3) .* a;

% ---- Eq. (19): cmin(x) -------------------------------------------------
%   c0=3.681, c1=5.033, alpha=6.948, x0=0.424
c0    = 3.681;  c1   = 5.033;
alpha = 6.948;  x0_c = 0.424;
cmin  = c0 + (c1 - c0) .* (atan(alpha .* (x - x0_c)) ./ pi + 0.5);

% ---- Eq. (20): sigma_min_inv(x) ----------------------------------------
%   sigma0_inv=1.047, sigma1_inv=1.646, beta=7.386, x1=0.526
s0inv = 1.047;  s1inv = 1.646;
beta  = 7.386;  x1    = 0.526;
smin_inv = s0inv + (s1inv - s0inv) .* (atan(beta .* (x - x1)) ./ pi + 0.5);

% ---- Eq. (18): B0 and B1 -----------------------------------------------
%   Normalisation point: x at z=0 for the given cosmology
x0_norm  = (cosmo.Omega_L ./ cosmo.Omega_m).^(1/3);   % x at a=1 (z=0)

%   cmin and sigma_min_inv at the normalisation point
cmin_norm    = c0 + (c1 - c0) .* (atan(alpha .* (x0_norm - x0_c)) ./ pi + 0.5);
smin_inv_norm= s0inv + (s1inv - s0inv) .* (atan(beta .* (x0_norm - x1)) ./ pi + 0.5);

%   Paper Eq. (18): B0 = cmin(x)/cmin(1.393)
%                   B1 = sigma_min_inv(x)/sigma_min_inv(1.393)
%   where 1.393 is the value of x at z=0 for Bolshoi/MultiDark cosmology
%   (Omegam=0.27, OmegaL=0.73).  For a general cosmology we normalise at
%   z=0 of the *input* cosmology (x0_norm), which recovers exactly B0=B1=1
%   at z=0 by construction.
B0 = cmin     ./ cmin_norm;
B1 = smin_inv ./ smin_inv_norm;

% ---- Eq. (15): rescaled sigma ------------------------------------------
sigma   = cosmo.sigmaM(M200c, z);
sigma_p = B1 .* sigma;

% ---- Eqs. (16-17): C(sigma') -------------------------------------------
%   A=2.881, b=1.257, c_p=1.022, d=0.060
A_p  = 2.881;
b_p  = 1.257;
c_p  = 1.022;
d_p  = 0.060;

C_sig = A_p .* ((sigma_p ./ b_p).^c_p + 1) .* exp(d_p ./ sigma_p.^2);

% ---- Eq. (14): final concentration -------------------------------------
c = B0 .* C_sig;

% Restore input shape
c = reshape(c, size(M200c));
end