function c = Klypin11(M, z, sample)
% Klypin11_concentration  Klypin et al. (2011) concentration model
%
%   c = Klypin11_concentration(M, z, sample)
%
%   Eq. (12): c(Mvir,z) = c0(z) * (Mvir/1e12)^alpha
%                         * [1 + (Mvir/M0(z))^beta]^-1
%
%   c0(z) is back-solved from the tabulated c(1e12) and M0 so that
%   the formula exactly reproduces Table 3 at the pivot mass 1e12.
%
%   Reference: Klypin, Trujillo-Gomez & Primack 2011, ApJ 740, 102

if nargin < 3 || isempty(sample)
    sample = 'distinct';
end

P = Klypin11_Table(sample);
% CHANGE 1: right after the input parsing, before M = M(:)
M_input = M;      % save original shape
M = M(:);         % force column for internal math

% ---- Subhalo branch (Eq. 11) -------------------------------------------
if strcmpi(sample, 'subhalo')
    c = P.c0_sub .* (M ./ P.Mpivot).^P.alpha_sub;
    c = reshape(c, size(M_input));   % restore caller's shape
    return
end

% ---- Distinct halo branch (Eq. 12) -------------------------------------
z_tab   = P.z;
c1e12   = P.c1e12;     % tabulated c at M=1e12  (what the paper calls c(1e12))
M0_tab  = P.M0;
alpha   = P.alpha;
beta    = P.beta;
Mpivot  = P.Mpivot;    % 1e12 h^-1 Msun

% Warn and clamp for out-of-range redshift
if z > z_tab(end)
    warning('Klypin11_concentration: z=%.2f exceeds tabulated range (max z=%.1f). Using z=%.1f parameters.', ...
        z, z_tab(end), z_tab(end));
    z = z_tab(end);
end

% Interpolate c(1e12) and M0 log-linearly between table nodes
if z <= z_tab(1)
    c12_z = c1e12(1);
    M0_z  = M0_tab(1);
elseif z >= z_tab(end)
    c12_z = c1e12(end);
    M0_z  = M0_tab(end);
else
    % log-linear interpolation (both quantities are smooth & monotonic)
    c12_z = exp(interp1(z_tab, log(c1e12), z, 'linear'));
    % M0 spans Inf at z=0 — interpolate only over finite entries
    finite_idx = isfinite(M0_tab);
    if z <= z_tab(find(finite_idx, 1, 'first'))
        M0_z = Inf;
    else
        M0_z = exp(interp1(z_tab(finite_idx), log(M0_tab(finite_idx)), z, 'linear'));
    end
end

% Back-solve c0 so that c(M=1e12) = c12_z exactly:
%   c12_z = c0 * 1^alpha * [1 + (1e12/M0)^beta]^-1
%   => c0  = c12_z * [1 + (Mpivot/M0)^beta]
if isinf(M0_z)
    upturn_pivot = 1.0;         % (1e12/Inf)^beta = 0
else
    upturn_pivot = 1 + (Mpivot / M0_z)^beta;
end
c0_eff = c12_z * upturn_pivot;

% Eq. (12)
if isinf(M0_z)
    upturn = ones(size(M));
else
    upturn = 1 + (M ./ M0_z).^beta;
end

c = c0_eff .* (M ./ Mpivot).^alpha ./ upturn;
c = reshape(c, size(M_input));   % restore caller's shape
end