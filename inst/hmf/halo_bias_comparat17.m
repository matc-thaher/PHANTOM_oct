function b = halo_bias_comparat17(sigma, delta_c)
% Comparat et al. (2017, MNRAS 469, 4157) halo bias.
% Same functional form as Bhattacharya+2011 (ST-type peak-background split),
% but parameters re-fit directly to avoid systematic offsets.
% No redshift dependence; parameters updated relative to the published
% version in the Colossus implementation.
%
% b(nu) = 1 + (a*nu^2 - q) / delta_c
%           + 2*p/delta_c / (1 + (a*nu^2)^p)
%
%
% Reference: Comparat et al. 2017, MNRAS 469, 4157  arXiv:1611.06386

    if nargin < 4 || isempty(delta_c)
        delta_c = collapse_overdensity();
    end

    a = 0.897;
    p = 0.624;
    q = 1.589;

    nu  = delta_c ./ sigma;
    nu2 = nu.^2;

    b = 1 + ((a .* nu2 - q) / delta_c) ...
          + ((2 * p / delta_c) ./ (1 + (a .* nu2).^p));
end
