function f = multiplicity_ST(sigma, delta_c)
% Sheth & Tormen (1999) multiplicity function
% Parameters: A=0.3222, a=0.707, p=0.3
     if nargin < 2 || isempty(delta_c)
        delta_c = collapse_overdensity();   % EdS value 1.6865
    end
    A = 0.3222;
    a = 0.707;
    p = 0.3;
    nu = delta_c ./ sigma;
    f  = A .* sqrt(2*a/pi) .* nu .* (1 + (nu.^2 * a).^(-p)) .* exp(-a*nu.^2/2);
end