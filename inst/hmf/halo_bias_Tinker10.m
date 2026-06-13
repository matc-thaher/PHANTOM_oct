function b = halo_bias_Tinker10(sigma, delta_c, Delta, z, cosmo)
% Tinker et al. (2010, ApJ 724, 878) halo bias.
% This is the current standard bias used alongside Tinker+2008 HMF.
%
% b(nu) = 1 - A*nu^a/(nu^a + delta_c^a) + B*nu^b + C*nu^c
%
% Table 2 parameters for Delta=200:
    if nargin < 2 || isempty(delta_c)
        delta_c = collapse_overdensity('corrections', true, 'z', z, 'cosmo', cosmo);   % EdS value 1.6865
    end
    y = log10(Delta);

    % Eq. 6 in Tinker+2010, Delta=200 fiducial
    A_t = 1.0 + (0.24 * y * exp(-(4/y)^4));
    a_t = (0.44 * y) - 0.88;
    B_t = 0.183;
    b_t = 1.5;
    C_t = 0.019 + (0.107*y) + (0.19*exp(-(4/y)^4));
    c_t = 2.4;

    nu = delta_c ./ sigma;

    b = 1 - (A_t .* nu.^a_t ./ (nu.^a_t + delta_c^a_t)) ...
          + (B_t .* nu.^b_t) ...
          + (C_t .* nu.^c_t);
end