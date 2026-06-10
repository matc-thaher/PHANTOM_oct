function c = Diemer15(M200c, z, cosmo, statistic)
% Diemer15_concentration  Diemer & Kravtsov (2015) concentration model
%
%   c = Diemer15(M200c, z, cosmo)
%   c = Diemer15(M200c, z, cosmo, statistic)
%
%   Universal c_200c-nu model depending on peak height and the local
%   slope of the power spectrum.  Valid for any cosmology, mass, redshift.
%
%   INPUTS
%   M200c     : halo mass M_200c [Msun/h], scalar or vector
%   z         : redshift (scalar)
%   cosmo     : cosmology struct with cosmo.sigmaM, cosmo.neff
%   statistic : 'median' (default) | 'mean'
%
%   OUTPUT
%   c         : concentration c_200c, same size as M200c
%
%   Reference: Diemer & Kravtsov 2015, ApJ 799, 108

if nargin < 4 || isempty(statistic)
    statistic = 'median';
end

% Load parameters from Table 3
P = Diemer15_Table(statistic);

kappa = P.kappa;
phi0  = P.phi0;   phi1 = P.phi1;
eta0  = P.eta0;   eta1 = P.eta1;
alpha = P.alpha;
beta  = P.beta;

% Peak height  nu = delta_c / sigma(M, z)
delta_c = 1.686;
sigma   = cosmo.sigmaM(M200c, z);
nu      = delta_c ./ sigma;
nu      = max(nu, 0.1);          % floor to avoid divergence at very low nu

% Local power-spectrum slope at kappa-scaled Lagrangian radius (Eq. 12)
n_eff = cosmo.neff(M200c, z, kappa);

% Eq. (10): concentration floor and its location in nu
c_min  = phi0 + phi1 .* n_eff;
nu_min = eta0 + eta1 .* n_eff;

% Eq. (9): double power-law in nu/nu_min
c = 0.5 .* c_min .* ( (nu ./ nu_min).^(-alpha) + (nu ./ nu_min).^beta );
end