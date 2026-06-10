function ratio = einasto_mass_ratio(c_array, alpha)
% einasto_mass_ratio  Ratio of Einasto enclosed mass within r_{-2}
%                     to enclosed mass within r_{-2}*c  (= R_200c).
%
%   ratio = einasto_mass_ratio(c_array, alpha)
%
%   Uses the regularised lower incomplete gamma function:
%     M(<r) proportional to gammainc( (2/alpha)*(r/r_{-2})^alpha , 3/alpha )
%
%   At r = r_{-2}  : argument = 2/alpha
%   At r = R_200c  : argument = (2/alpha) * c^alpha
%
%   Fixed alpha = 0.18 per Ludlow et al. (2016).
%   Do NOT use the Gao+2008 mass-dependent alpha here.
%
%   Reference: Ludlow et al. 2016, MNRAS 460, 1214

x_s   = 2 / alpha;                          % argument at r_{-2}
x_c   = x_s .* c_array .^ alpha;            % argument at R_200c
a_gam = 3 / alpha;                          % shape parameter

% gammainc(x,a) in MATLAB is the regularised lower incomplete gamma.
% The unnormalised forms cancel in the ratio so regularised is fine.
ratio = gammainc(x_s, a_gam) ./ gammainc(x_c, a_gam);
end