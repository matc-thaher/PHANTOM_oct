function P = Dutton14_Table(mdef)
% Dutton14_Table  Parameters for the Dutton & Maccio (2014) concentration model
%
%   P = Dutton14_Table(mdef)
%
%   Returns the best-fit parameters from Dutton & Maccio 2014, MNRAS 441, 3359.
%   Calibrated on relaxed haloes in Planck 2013 cosmology, z=0 to z=5.
%
%   MDEF options (case-insensitive)
%   --------------------------------
%   '200c'   M_200c definition  — Eqs. (10-11)
%   'vir'    M_vir  definition  — Eqs. (12-13)
%
%   Fields of output struct P
%   -------------------------
%   P.a0       constant term in exponential fit to a(z)   Eq. (11) or (13)
%   P.a1       zero-point at z=0                          Eq. (11) or (13)
%   P.eta      exponential decay coefficient              Eq. (11) or (13)
%   P.phi      power of z in exponential                  Eq. (11) or (13)
%   P.b0       slope at z=0                               Eq. (10) or (12)
%   P.b1       linear slope evolution with z              Eq. (10) or (12)
%   P.Mpivot   pivot mass [h^-1 Msun]                     fixed at 1e12
%
%   Reference: Dutton & Maccio 2014, MNRAS 441, 3359, Eqs. (7)-(13)

switch lower(mdef)

    case '200c'
        % c200 vs M200  — Eqs. (10) and (11)
        % b(z) = b0 + b1*z
        P.b0    = -0.101;
        P.b1    =  0.026;
        % a(z) = a0 + (a1 - a0) * exp(-eta * z^phi)
        P.a0    =  0.520;
        P.a1    =  0.905;   % value at z=0  (a(0) = a1 by construction)
        P.eta   =  0.617;
        P.phi   =  1.21;

    case 'vir'
        % cvir vs Mvir  — Eqs. (12) and (13)
        P.b0    = -0.097;
        P.b1    =  0.024;
        P.a0    =  0.537;
        P.a1    =  1.025;
        P.eta   =  0.718;
        P.phi   =  1.08;

    otherwise
        error('Dutton14_Table: unknown mdef "%s". Valid: ''200c'', ''vir''.', mdef);
end

P.Mpivot = 1e12;   % h^-1 Msun — fixed pivot, see Eq. (7)
end