function Mc_new = thaher_CHMR(Mh, scale_factor, zeta, zeta_0, M_min0, eta, gamma)
% Thaher's modified core-halo mass relation with low-mass suppression
%
% INPUTS:
%   Mh           - Halo virial mass array [M_sun]  (N x 1 or 1 x N)
%   scale_factor - Cosmological scale factor array (1 x T)
%   zeta         - Zeta parameter array            (1 x T)
%   zeta_0       - Reference zeta (scalar or 1 x T)
%   M_min0       - Minimum halo mass [M_sun]       (scalar)
%   eta          - Transition mass factor: M_t = eta * M_min0 (default: 1.5)
%   gamma        - Sharpness of the suppression bend          (default: 2)
%
% OUTPUT:
%   Mc_new       - Soliton core mass [M_sun]                  (N x T)

    if nargin < 6, eta   = 3/2; end
    if nargin < 7, gamma = 2;   end

    suppression = (1 + (M_min0 ./ Mh').^gamma) .^ ((-1/2) / eta);

    Mc_new = 0.25 ...
        .* scale_factor.^(-0.5) ...
        .* (zeta ./ zeta_0).^(1/6) ...
        .* (Mh' ./ M_min0).^(1/3) ...
        .* suppression ...
        .* M_min0;
end