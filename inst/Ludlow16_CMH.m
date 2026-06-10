function Mcoll = Ludlow16_CMH(z, M0, f, cosmo)
% Ludlow16_CMH  Collapsed Mass History from EPS theory
%
%   Mcoll = Ludlow16_CMH(z, M0, f, cosmo)
%
%   Computes the total mass in collapsed progenitors more massive than
%   f * M0 at redshift z, using the EPS analytic form (eq. 3):
%
%     Mcoll(z) = M0 * erfc( [delta_sc(z) - delta_sc(z0)] /
%                            sqrt(2 * [sigma^2(f*M0) - sigma^2(M0)]) )
%
%   where delta_sc(z) = 1.26 / D(z)  (Ludlow+2016 use 1.26, not 1.686,
%   for the collapse threshold — see footnote 4 of the paper).
%
%   INPUTS
%   z      : redshift array at which to evaluate Mcoll
%   M0     : final halo mass at z0 [Msun/h]
%   f      : progenitor mass fraction threshold (paper uses f=0.02)
%   cosmo  : struct with fields:
%              cosmo.sigmaM(M, z)   — rms linear fluctuation at z
%              cosmo.growthFactor(z)— linear growth factor D(z),
%                                     normalised so D(0)=1
%
%   OUTPUT
%   Mcoll  : collapsed mass [Msun/h], same size as z
%
%   Reference: Ludlow et al. 2016, MNRAS 460, 1214, eq.(3)

    delta_sc0 = 1.686;                        % collapse threshold (footnote 4)

    % D0 = growth_factor_D(0, cosmo);
    % Dz = growth_factor_D(z, cosmo);
    D0 = cosmo.D(0);
    Dz = cosmo.D(z);

    delta0    = delta_sc0 ./ D0;             % at z0=0
    deltaz    = delta_sc0 ./ Dz;            % at redshift z

    sigma_fM0 = cosmo.sigmaM(f * M0, 0);
    sigma_M0  = cosmo.sigmaM(M0, 0);

    denom     = sqrt(2 * (sigma_fM0^2 - sigma_M0^2));
    arg       = (deltaz - delta0) ./ denom;

    Mcoll     =  erfc(arg); % M0 .*
end