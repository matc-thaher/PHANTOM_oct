function rho = Einasto_fit_wrapper(r, rhos, rs, alpha_e)
% EINASTO_FIT_WRAPPER  Einasto profile for direct fitting.
%
% rho(r) = rhos * exp( -(2/alpha_e) * [ (r/rs)^alpha_e - 1 ] )
%
% Free parameters: rhos [Msun/kpc^3], rs [kpc], alpha_e (shape, ~0.16-0.30)
%
% Reference: Einasto 1965; Merritt et al. 2006, AJ 132, 2685

  r    = r(:);
  exponent = -(2 ./ alpha_e) .* ((r ./ rs).^alpha_e - 1);
  rho  = rhos .* exp(exponent);
end