function rho = DK14_fit_wrapper(r, rhos, rs, alpha_e)
% DK14_FIT_WRAPPER  Diemer & Kravtsov (2014) profile for direct fitting.
%
% Uses fixed splashback parameters (beta=4, gamma_t=8, mass-selected)
% and a fixed outer term slope s_e=1.5.
% Only the inner Einasto amplitude (rhos), scale radius (rs), and
% shape parameter (alpha_e) are treated as free fitting parameters.
%
% rho_inner(r) = rhos * exp( -(2/alpha_e)*[ (r/rs)^alpha_e - 1 ] )
% f_trans(r)   = [ 1 + (r/rt)^4 ]^(-2)           [beta=4, gamma_t=8]
% rho_outer(r) = rhos * 1e-4 * (r/rs)^(-1.5)      [approximate outer term]
%
% For a physically exact outer term you need rho_m(z); this wrapper
% uses a scaled approximation so the fitter stays self-contained.
%
% Reference: Diemer & Kravtsov 2014, ApJ 789, 1

  r  = r(:);
  rt = 1.6 .* rs;   % approximate: rt ~ 1.6 rs for typical nu~1 halos

  rho_inner = Einasto_fit_wrapper(r, rhos, rs, alpha_e);
  f_trans   = (1 + (r ./ rt).^4).^(-2);
  rho_outer = rhos .* 1e-4 .* (r ./ rs).^(-1.5);

  rho = rho_inner .* f_trans + rho_outer;
end