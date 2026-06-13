function f = multiplicity_Courtin11(sigma)
% Courtin et al. (2011), MNRAS 410, 1911
% ST functional form with different parameters and fixed delta_c = 1.673.
% FOF, no explicit redshift dependence.
    delta_c = 1.673;    % fixed in this model
    A = 0.348;
    a = 0.695;
    p = 0.1;
    nu = delta_c ./ sigma;
    f = A .* sqrt(a .* 2/pi) .* nu .* exp(-0.5.*nu.^2 .* a) .* (1 + (nu.^2 .* a).^(-p));
end