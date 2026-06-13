function b = halo_bias_PS(sigma, delta_c)
% Cole & Kaiser 1989, Mo & White 1996
% Press-Schechter linear halo bias
% b(nu) = 1 + (nu^2 - 1) / delta_c
    if nargin < 2 || isempty(delta_c)
        delta_c = collapse_overdensity();   % EdS value 1.6865
    end
    nu = delta_c ./ sigma;
    b  = 1 + ((nu.^2 - 1) / delta_c);
end