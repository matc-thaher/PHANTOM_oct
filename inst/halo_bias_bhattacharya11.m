function b = halo_bias_bhattacharya11(sigma, delta_c, z)
% Bhattacharya et al. (2011, ApJ 732, 122) halo bias.
% Peak-background split applied to the Bhattacharya+2011 mass function
% (ST-type form with redshift-dependent parameters).
%
% b(nu,z) = 1 + (a(z)*nu^2 - q) / delta_c
%             + 2*p(z)/delta_c / (1 + (a(z)*nu^2)^p(z))
%
% Redshift-dependent best-fit parameters (Table 2):
%   a(z) = 0.788 * (1+z)^{-0.01}
%   p(z) = 0.807
%   q    = 1.795
%
% Note: The authors warn this does not match simulation bias as well as
% directly calibrated models (e.g. Tinker+2010).
%
% Reference: Bhattacharya et al. 2011, ApJ 732, 122  arXiv:1005.2239

    if nargin < 2 || isempty(delta_c)
        delta_c = collapse_overdensity();
    end
    if nargin < 3 || isempty(z)
        z = 0;
    end

    a_z = 0.788 .* (1 + z).^(-0.01);
    p_z = 0.807;
    q   = 1.795;

    nu  = delta_c ./ sigma;
    nu2 = nu.^2;

    b = 1 + (((a_z .* nu2) - q) / delta_c) ...
          + (((2 * p_z) / delta_c) ./ (1 + (a_z .* nu2).^p_z));
end
