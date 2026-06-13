function c = Duffy08(M, z, mdef, profile, sample, redshift_range, M_pivot)
% Duffy08_concentration  Duffy et al. (2008) concentration model
%
%   c = Duffy08(M, z, mdef)
%   c = Duffy08(M, z, mdef, profile)
%   c = Duffy08(M, z, mdef, profile, sample)
%   c = Duffy08(M, z, mdef, profile, sample, redshift_range)
%   c = Duffy08(M, z, mdef, profile, sample, redshift_range, M_pivot)
%
%   Power-law fit (Eq. 4 of Duffy+08):
%       c = A * (M / M_pivot)^B * (1+z)^C
%
%   Parameters are loaded from Duffy08_Table, which contains the complete
%   Table 1 of the paper for both NFW and Einasto profiles.
%
%   Calibrated for WMAP5 cosmology:
%     (Om, Ob, OL, h, sigma8, ns) = (0.258, 0.0441, 0.742, 0.719, 0.796, 0.963)
%   Valid range: 1e11 < M < 1e15  h^-1 Msun,   0 < z < 2.
%
% -------------------------------------------------------------------------
%   INPUTS
%   M              : halo mass [Msun/h], scalar or vector
%   z              : redshift (scalar)
%   mdef           : '200c' | 'vir' | '200m'
%                    (default: '200c')
%   profile        : 'NFW' | 'Einasto'
%                    (default: 'NFW')
%   sample         : 'full' | 'relaxed'
%                    (default: 'full')
%   redshift_range : 'z0'  — use the z=0 fit (C=0, more accurate at z=0)
%                    'z0_2'— use the z=0 to 2 fit (C free, use for z>0)
%                    (default: 'z0_2')
%   M_pivot        : pivot mass [Msun/h]
%                    (default: 2e12, the value used in the paper)
%                    Override this if you want to rescale the relation to
%                    a different pivot, e.g. 1e14 to match Neto+07.
%
%   OUTPUT
%   c              : concentration, same shape as M
%
% -------------------------------------------------------------------------
%   NOTES
%   - NFW concentration is defined as c = R_Delta / r_s
%   - Einasto concentration is defined as c = R_Delta / r_{-2}
%     where r_{-2} is the radius where the log-slope equals -2 (analogous
%     to r_s in NFW).  Values are therefore close to NFW concentrations.
%   - For the z=0-only fit ('z0'), C=0 and the (1+z)^C term equals 1,
%     so the result is independent of redshift.  Use 'z0_2' for z > 0.
%   - The Einasto alpha parameter is NOT output here; use the Gao+08
%     relation alpha(nu) if you need it for profile evaluation.
%
% -------------------------------------------------------------------------
%   Reference: Duffy et al. 2008, MNRAS 390, L64, Table 1
% -------------------------------------------------------------------------

% ---- Default arguments -------------------------------------------------
if nargin < 4 || isempty(profile),        profile        = 'NFW';   end
if nargin < 5 || isempty(sample),         sample         = 'full';  end
if nargin < 6 || isempty(redshift_range), redshift_range = 'z0_2';  end
if nargin < 7 || isempty(M_pivot),        M_pivot        = 2e12;    end

% ---- Load parameters from table ----------------------------------------
P = Duffy08_Table(mdef, profile, sample, redshift_range);

% ---- Warn if user pivot differs from the paper's pivot -----------------
if M_pivot ~= P.M_pivot
    warning(['Duffy08_concentration: using M_pivot = %.3e instead of ' ...
             'the paper value %.3e.\nThis rescales A; B and C are unchanged.'], ...
             M_pivot, P.M_pivot);
    % Rescale A to the new pivot so the relation is self-consistent:
    % A_new * (M/M_new)^B = A_old * (M/M_old)^B
    % => A_new = A_old * (M_new/M_old)^B
    P.A = P.A * (M_pivot / P.M_pivot)^P.B;
end

% ---- Evaluate Eq. 4 ----------------------------------------------------
c = P.A .* (M ./ M_pivot).^P.B .* (1 + z).^P.C;

end