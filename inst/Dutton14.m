function c = Dutton14(M, z, mdef, M_pivot)
% Dutton14_concentration  Dutton & Maccio (2014) concentration model
%
%   c = Dutton14(M, z, mdef)
%   c = Dutton14(M, z, mdef, M_pivot)
%
%   Power-law fit:
%     log10(c) = a(z) + b(z) * log10(M / M_pivot)        Eq. (7)
%
%   with redshift-dependent parameters:
%     b(z) = b0 + b1*z                                    Eq. (10) or (12)
%     a(z) = a0 + (a1-a0) * exp(-eta * z^phi)             Eq. (11) or (13)
%
%   Calibrated for relaxed haloes in Planck 2013 cosmology.
%   Valid range: M > 1e10 h^-1 Msun,  0 <= z <= 5.
%
%   INPUTS
%   M       : halo mass [h^-1 Msun], scalar or vector
%   z       : redshift, scalar
%   mdef    : '200c' | 'vir'
%   M_pivot : (optional) pivot mass [h^-1 Msun], default = 1e12
%             Override only if you need a custom pivot; the paper
%             fixes this at 10^12 h^-1 Msun (Eq. 7).
%
%   OUTPUT
%   c       : concentration, same shape as M
%
%   Reference: Dutton & Maccio 2014, MNRAS 441, 3359, Eqs. (7)-(13)

% ---- Save input shape, flatten for computation -------------------------
M_input = M;
M       = M(:);

% ---- Load parameters ---------------------------------------------------
P = Dutton14_Table(mdef);

% ---- Override pivot if user supplies one -------------------------------
if nargin < 4 || isempty(M_pivot)
    M_pivot = P.Mpivot;   % 1e12 h^-1 Msun (paper default)
end

% ---- Eq. (10) or (12): b(z) = b0 + b1*z --------------------------------
b = P.b0 + P.b1 .* z;

% ---- Eq. (11) or (13): a(z) = a0 + (a1-a0)*exp(-eta*z^phi) ------------
a = P.a0 + (P.a1 - P.a0) .* exp(-P.eta .* z.^P.phi);

% ---- Eq. (7): log10(c) = a + b*log10(M/M_pivot) -----------------------
log10c = a + b .* log10(M ./ M_pivot);
c      = 10.^log10c;

% ---- Restore input shape -----------------------------------------------
c = reshape(c, size(M_input));
end