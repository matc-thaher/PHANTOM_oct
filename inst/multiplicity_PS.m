function f = multiplicity_PS(sigma, delta_c)
% Press & Schechter (1974) multiplicity function
% f(sigma) = sqrt(2/pi) * (delta_c/sigma) * exp(-delta_c^2 / 2sigma^2)
    if nargin < 2 || isempty(delta_c)
        delta_c = collapse_overdensity();   % EdS value 1.6865
    end
    nu = delta_c ./ sigma;
    f  = sqrt(2/pi) .* nu .* exp(-nu.^2 / 2);
end