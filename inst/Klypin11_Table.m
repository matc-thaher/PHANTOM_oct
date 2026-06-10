function P = Klypin11_Table(sample)
% Klypin11_Table  Parameters for Klypin et al. (2011) concentration model
%
%   P = Klypin11_Table(sample)
%
%   Returns best-fit parameters from Table 3 of Klypin, Trujillo-Gomez &
%   Primack (2011), ApJ 740, 102.  All fits are for M_vir (virial mass
%   defined via the Bryan & Norman 1998 top-hat overdensity).
%
%   SAMPLE options (case-insensitive)
%   ----------------------------------
%   'distinct'     All distinct (field) haloes, Eq. (10)/(12)
%   'subhalo'      Subhaloes, Eq. (11)
%
%   The functional form for distinct haloes is (Eq. 12):
%
%     c(Mvir, z) = c0(z) * (Mvir / 1e12)^(-0.075)
%                  * [ 1 + (Mvir / M0(z))^0.26 ]^(-1)
%
%   where c0(z) and M0(z) are tabulated in Table 3.
%   For the z=0 formula only (Eq. 10):
%     c(Mvir) = 9.60 * (Mvir / 1e12)^(-0.075)
%
%   For subhaloes (Eq. 11):
%     c(Msub) = 12   * (Msub / 1e12)^(-0.12)
%
%   Fields of output struct P
%   -------------------------
%   P.sample     input string (for bookkeeping)
%   P.z          redshift grid [1 x 7]
%   P.c0         c0(z)  [1 x 7]
%   P.M0         M0(z)  [h^-1 Msun], [1 x 7]
%   P.alpha      power-law slope of (Mvir/1e12)
%   P.cmin       minimum concentration at each z [1 x 7]
%   P.c1e12      c at Mvir=1e12 h^-1 Msun [1 x 7]   (Table 3 column)
%
%   For 'subhalo' only P.c0_sub and P.alpha_sub are filled; redshift
%   evolution is not tabulated for subhaloes in the paper.
%
%   Reference: Klypin, Trujillo-Gomez & Primack 2011, ApJ 740, 102, Table 3

switch lower(strtrim(sample))

    case 'distinct'
    P.z      = [0.0,   0.5,    1.0,    2.0,    3.0,    5.0];
    P.c0     = [9.60,  7.08,   5.45,   3.67,   2.83,   2.34];
    P.M0     = [Inf,   1.5e17, 2.50e15, 6.80e13, 6.30e12, 6.60e11];  % Inf at z=0
    P.c1e12  = [9.60,  7.21,   5.82,   4.60,   4.40,   5.00];        % exact Table 3
    P.cmin   = [NaN,   5.2,    5.1,    4.6,    4.2,    4.0];
    P.alpha  = -0.075;
    P.beta   =  0.26;
    P.Mpivot = 1e12;
    P.sample = 'distinct';

    case 'subhalo'
        % Eq. (11) — no redshift evolution tabulated
        P.c0_sub    = 12.0;
        P.alpha_sub = -0.12;
        P.Mpivot    = 1e12;   % [h^-1 Msun]
        P.sample    = 'subhalo';
        % Redshift fields left empty — not provided in the paper
        P.z   = 0.0;
        P.c0  = 12.0;
        P.M0  = NaN;
        P.c1e12 = 12.0;
        P.cmin  = NaN;
        P.alpha = -0.12;
        P.beta  = NaN;
        P.Mpivot = 1e12;

    otherwise
        error('Klypin11_Table: unknown sample "%s".\n  Valid: ''distinct'', ''subhalo''.', ...
            sample);
end
end