function P = Duffy08_Table(mdef, profile, sample, redshift_range)
% Duffy08_Table  Best-fit parameters for the Duffy et al. (2008) model
%
%   P = Duffy08_Table(mdef, profile, sample, redshift_range)
%
%   Returns the A, B, C parameters from Table 1 of Duffy et al. 2008,
%   MNRAS 390, L64.  All fits use M_pivot = 2e12 h^-1 Msun.
%
%   The power-law form is (Eq. 4):
%       c = A * (M / M_pivot)^B * (1+z)^C
%   For the z=0-only fits, C is fixed to 0.
%
% -------------------------------------------------------------------------
%   INPUTS
%   mdef           : halo mass definition
%                    '200c'  — M_200c (overdensity 200 x critical)
%                    'vir'   — M_vir  (Bryan & Norman 1998 overdensity)
%                    '200m'  — M_200m (overdensity 200 x mean)
%
%   profile        : density profile used in the fit
%                    'NFW'    — Navarro, Frenk & White 1997
%                    'Einasto' — Einasto 1965
%
%   sample         : halo sample
%                    'full'    — all haloes passing resolution criteria
%                    'relaxed' — haloes with centre-of-mass offset < 0.07*rvir
%
%   redshift_range : redshift range used in fit
%                    'z0'   — z = 0 only  (C forced to 0)
%                    'z0_2' — z = 0 to 2  (C free)
%
%   OUTPUT
%   P.A, P.B, P.C  — best-fit parameters
%   P.M_pivot      — pivot mass [h^-1 Msun]  = 2e12 for all entries
%
% -------------------------------------------------------------------------
%   Reference: Duffy et al. 2008, MNRAS 390, L64, Table 1
%              WMAP5 cosmology: (Om,Ob,OL,h,s8,ns)
%                               (0.258,0.0441,0.742,0.719,0.796,0.963)
% -------------------------------------------------------------------------

P.M_pivot = 2e12;   % h^-1 Msun — same for all entries in Table 1

% Build a key string for the switch: 'mdef_profile_sample_zrange'
key = lower(sprintf('%s_%s_%s_%s', mdef, profile, sample, redshift_range));

switch key

    % =====================================================================
    % NFW — M_200c
    % =====================================================================
    case '200c_nfw_full_z0'
        P.A =  5.74;  P.B = -0.097;  P.C =  0.00;
    case '200c_nfw_full_z0_2'
        P.A =  5.71;  P.B = -0.084;  P.C = -0.47;
    case '200c_nfw_relaxed_z0'
        P.A =  6.67;  P.B = -0.092;  P.C =  0.00;
    case '200c_nfw_relaxed_z0_2'
        P.A =  6.71;  P.B = -0.091;  P.C = -0.44;

    % =====================================================================
    % NFW — M_vir
    % =====================================================================
    case 'vir_nfw_full_z0'
        P.A =  7.96;  P.B = -0.091;  P.C =  0.00;
    case 'vir_nfw_full_z0_2'
        P.A =  7.85;  P.B = -0.081;  P.C = -0.71;
    case 'vir_nfw_relaxed_z0'
        P.A =  9.23;  P.B = -0.089;  P.C =  0.00;
    case 'vir_nfw_relaxed_z0_2'
        P.A =  9.23;  P.B = -0.090;  P.C = -0.69;

    % =====================================================================
    % NFW — M_200m
    % =====================================================================
    case '200m_nfw_full_z0'
        P.A = 10.39;  P.B = -0.089;  P.C =  0.00;
    case '200m_nfw_full_z0_2'
        P.A = 10.14;  P.B = -0.081;  P.C = -1.01;
    case '200m_nfw_relaxed_z0'
        P.A = 12.00;  P.B = -0.087;  P.C =  0.00;
    case '200m_nfw_relaxed_z0_2'
        P.A = 11.93;  P.B = -0.090;  P.C = -0.99;

    % =====================================================================
    % Einasto — M_200c
    % =====================================================================
    case '200c_einasto_full_z0'
        P.A =  6.48;  P.B = -0.127;  P.C =  0.00;
    case '200c_einasto_full_z0_2'
        P.A =  6.40;  P.B = -0.108;  P.C = -0.62;
    case '200c_einasto_relaxed_z0'
        P.A =  7.70;  P.B = -0.127;  P.C =  0.00;
    case '200c_einasto_relaxed_z0_2'
        P.A =  7.74;  P.B = -0.123;  P.C = -0.60;

    % =====================================================================
    % Einasto — M_vir
    % =====================================================================
    case 'vir_einasto_full_z0'
        P.A =  9.03;  P.B = -0.122;  P.C =  0.00;
    case 'vir_einasto_full_z0_2'
        P.A =  8.82;  P.B = -0.106;  P.C = -0.87;
    case 'vir_einasto_relaxed_z0'
        P.A = 10.79;  P.B = -0.125;  P.C =  0.00;
    case 'vir_einasto_relaxed_z0_2'
        P.A = 10.77;  P.B = -0.124;  P.C = -0.87;

    % =====================================================================
    % Einasto — M_200m
    % =====================================================================
    case '200m_einasto_full_z0'
        P.A = 11.84;  P.B = -0.124;  P.C =  0.00;
    case '200m_einasto_full_z0_2'
        P.A = 11.39;  P.B = -0.107;  P.C = -1.16;
    case '200m_einasto_relaxed_z0'
        P.A = 14.03;  P.B = -0.116;  P.C =  0.00;
    case '200m_einasto_relaxed_z0_2'
        P.A = 13.96;  P.B = -0.119;  P.C = -1.17;

    otherwise
        error(['Duffy08_Table: unrecognised combination:\n' ...
               '  mdef="%s", profile="%s", sample="%s", z_range="%s"\n\n' ...
               'Valid mdef    : 200c | vir | 200m\n' ...
               'Valid profile : NFW  | Einasto\n' ...
               'Valid sample  : full | relaxed\n' ...
               'Valid z_range : z0   | z0_2'], ...
               mdef, profile, sample, redshift_range);
end
end