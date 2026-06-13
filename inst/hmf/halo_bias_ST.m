function b = halo_bias_ST(sigma, delta_c)
% Sheth & Tormen (1999) halo bias — peak-background split applied to the
% ST ellipsoidal collapse mass function.
%
% Often cited as Sheth, Mo & Tormen (2001) in the literature, but the bias
% formula with a=0.707, p=0.3 is derived in Sheth & Tormen (1999).
%
% b(nu) = 1 + (a*nu^2 - 1)/delta_c + 2p/delta_c / (1 + (a*nu^2)^p)
%
% Reference: Sheth & Tormen 1999, MNRAS 308, 119   arXiv:astro-ph/9901122
    if nargin < 2 || isempty(delta_c)
        delta_c = collapse_overdensity();   % EdS value 1.6865
    end
    a = 0.707;  p = 0.3;
    nu = delta_c ./ sigma;
    nu2 = nu.^2;
    b  = 1 + ((a*nu2 - 1)/delta_c) + ((2*p/delta_c) ./ (1 + (a*nu2).^p));
end