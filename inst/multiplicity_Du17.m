function f = multiplicity_Du17(sigma, M, z, m22, cosmo, delta_c)
% Du, Behrens & Niemeyer (2017), MNRAS 465, 941
% FDM halo mass function via excursion set with mass-dependent barrier.
%
% The FDM collapse barrier is:
%   delta_fdm(M,z) = G(M) * delta_cdm(z)
%
% where G(M) is a fitting function from Marsh (2016), Eq. (11) of Du+2017:
%   G(M) = hF(x)*exp(a3*x^-a4) + [1-hF(x)]*exp(a5*x^-a6)
%   x    = M / M_J
%   hF(x) = 0.5 * {1 - tanh[M_J * (x - a2)]}
%   M_J  = 1e8 * a1 * (m22)^(-3/2) * (Om*h^2/0.14)^(1/4) [h^-1 M_sun]
%   {a1..a6} = {3.4, 1.0, 1.8, 0.5, 1.7, 0.9}
%
% f(sigma) is computed with the Sheth-Tormen formula using delta_fdm
% as the effective barrier (consistent with Du+2017 Section 2.1).
%
% INPUT:
%   sigma   : sigma(M,z) from FDM linear P(k)
%   M       : halo mass array [h^-1 M_sun]
%   z       : redshift
%   m22     : FDM boson mass [1e-22 eV/c^2]
%   cosmo   : PHANTOM cosmo struct
%   delta_c : (optional) CDM collapse threshold

    if nargin < 6 || isempty(delta_c)
        delta_c = collapse_overdensity('corrections', true, 'z', z, 'cosmo', cosmo);
    end

    % Fitting parameters (Du+2017, Eq. 11; Marsh 2016)
    a1 = 3.4;  a2 = 1.0;
    a3 = 1.8;  a4 = 0.5;
    a5 = 1.7;  a6 = 0.9;

    % Jeans mass M_J [h^-1 M_sun] (Du+2017, Eq. 14)
    M_J = 1e8 * a1 * m22^(-3/2) * (cosmo.Omega_m * cosmo.h^2 / 0.14)^(0.25);

    % Dimensionless mass variable x = M / M_J
    x = M ./ M_J;

    % Interpolation function hF(x) (Du+2017, Eq. 13)
    hF = 0.5 .* (1 - tanh(M_J .* (x - a2)));

    % Growth suppression G(M) (Du+2017, Eq. 11)
    G = hF .* exp(a3 .* x.^(-a4)) + (1 - hF) .* exp(a5 .* x.^(-a6));

    % Mass-dependent FDM collapse barrier
    delta_fdm = G .* delta_c;

    % ST multiplicity with FDM barrier (Du+2017 Section 2.1 approximation)
    p  = 0.3;  q = 0.707;  A_ST = 0.3222;
    nu = delta_fdm ./ sigma;
    S  = sigma.^2;
    f  = A_ST .* sqrt(q / (2*pi)) .* nu .* (1 + (q .* nu.^2).^(-p)) ...
              .* exp(-q .* nu.^2 / 2) ./ S;
end