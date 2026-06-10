function [r_x, rho_x, found_flag] = find_intersection_from_fits(sol, nfw, r, r_fallback)

  found_flag = false;
  r_x = NaN;
  rho_x = NaN;

  r_test = logspace(log10(min(r(r>0))), log10(max(r)), 2000);

  rho_sol = profile_value(sol, r_test);
  rho_nfw = profile_value(nfw, r_test);

  diff_arr = rho_sol - rho_nfw;
  idx = find(diff_arr(1:end-1) .* diff_arr(2:end) <= 0, 1, 'first');

  if ~isempty(idx)
      r1 = r_test(idx);
      r2 = r_test(idx+1);

      f = @(x) profile_value(sol, x) - profile_value(nfw, x);
      try
          r_x = fzero(f, [r1, r2]);
      catch
          r_x = 0.5 * (r1 + r2);
      end

      rho_x = profile_value(sol, r_x);
      found_flag = true;
      return
  end

  r_x = r_fallback;
  rho_x = profile_value(nfw, r_fallback);
end