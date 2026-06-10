function f = multiplicity_Angulo12(sigma)
% Angulo et al. (2012), MNRAS 426, 2046, Eq. 2
% FOF, no redshift dependence. Note: corrects typo in published equation.
    f = 0.201 .* ((2.08./sigma).^1.7 + 1) .* exp(-1.172./sigma.^2);
end