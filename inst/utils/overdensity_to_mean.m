function [Delta_m, Delta_vir_c] = overdensity_to_mean(Delta, mdef, z, cosmo)
% overdensity_to_mean   Convert overdensity threshold to mean-density units
%                       and return the virial overdensity w.r.t. rho_crit.
%
% INPUT:
%   Delta : overdensity value (e.g. 200, 500); ignored when mdef = 'vir'
%   mdef  : 'mean' | 'crit' | 'vir'
%   z     : redshift
%   cosmo : PHANTOM cosmo struct
%
% OUTPUT:
%   Delta_m     : overdensity w.r.t. mean matter density
%   Delta_vir_c : virial overdensity w.r.t. rho_crit (Bryan & Norman 1998)
%                 always returned regardless of mdef — useful for Despali16
%
% USAGE:
%   [Delta_m, Delta_vir_c] = overdensity_to_mean(200, 'crit', z, cosmo);
%   [Delta_m, Delta_vir_c] = overdensity_to_mean(200, 'mean', z, cosmo);
%   [Delta_m, Delta_vir_c] = overdensity_to_mean([], 'vir',  z, cosmo);

    % Bryan & Norman (1998) virial overdensity — always computed
    x_BN        = cosmo.Omega_m_z(z) - 1;
    Delta_vir_c = 18*pi^2 + 82*x_BN - 39*x_BN^2;   % w.r.t. rho_crit

    switch lower(mdef)
        case 'mean'
            Delta_m = Delta;
        case 'crit'
            Delta_m = Delta .* cosmo.rhocrit(z) ./ cosmo.rhom(z);
        case 'vir'
            % Convert virial (w.r.t. rho_crit) to mean-density units
            Delta_m = Delta_vir_c .* cosmo.rhocrit(z) ./ cosmo.rhom(z);
        otherwise
            error('overdensity_to_mean: mdef must be ''mean'', ''crit'', or ''vir''.');
    end
end