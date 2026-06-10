function Mc = schive_CHMR(Mh, scale_factor, zeta, zeta_0, M_min0)
% Schive et al. core-halo mass relation
%
% INPUTS:
%   Mh           - Halo virial mass array [M_sun]  (N x 1 or 1 x N)
%   scale_factor - Cosmological scale factor array (1 x T)
%   zeta         - Zeta parameter array            (1 x T)
%   zeta_0       - Reference zeta (scalar or 1 x T)
%   M_min0       - Minimum halo mass [M_sun]       (scalar)
%
% OUTPUT:
%   Mc           - Soliton core mass [M_sun]        (N x T)

    Mc = 0.25 ...
        .* scale_factor.^(-0.5) ...
        .* (zeta ./ zeta_0).^(1/6) ...
        .* (Mh' ./ M_min0).^(1/3) ...
        .* M_min0;
end