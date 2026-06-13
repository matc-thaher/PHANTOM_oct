function delta_c = collapse_overdensity_fdm(M, cosmo, varargin)
% collapse_overdensity_fdm   Mass-dependent collapse threshold for FDM
%
% For fuzzy dark matter, quantum pressure raises the effective collapse
% barrier at masses near and below the Jeans mass M_J.  The threshold is
%
%   delta_c^fdm(M,z) = G(M) * delta_c^cdm(z)
%
% where G(M) = D_cdm(z) / D_fdm(M,z) is the ratio of the CDM growth
% factor to the suppressed FDM growth factor at mass scale M, following
% Marsh & Silk (2014) and Du et al. (2016).  The fitting function for G(M)
% is taken from Marsh (2016a) as implemented in Du et al. (2016), Eq. 11.
%
% For M >> M_J  :  G(M) -> 1  and  delta_c^fdm -> delta_c^cdm  (CDM limit)
% For M << M_J  :  G(M) >> 1 (collapse is suppressed by quantum pressure)
%
% INPUT:
%   M      [N x 1]  Halo mass in h^-1 M_sun
%   cosmo  struct   PHANTOM cosmo struct; must have cosmo.m22 set
%
% OPTIONAL name-value pairs:
%   'z'            scalar   Redshift for delta_c^cdm correction (default: 0)
%   'corrections'  logical  Apply Omega_m(z) correction to cdm baseline
%                           (default: false)
%
% OUTPUT:
%   delta_c  [N x 1]  FDM collapse overdensity at each mass M
%
% REFERENCES:
%   Marsh & Silk (2014), MNRAS 437, 2652
%   Marsh & Pop  (2015), MNRAS 451, 2479
%   Marsh        (2016), Phys. Rep. 643, 1   [fitting function for G]
%   Du, Behrens & Niemeyer (2016), MNRAS 465, 941
%
% USAGE:
%   delta_c = collapse_overdensity_fdm(M, cosmo)
%   delta_c = collapse_overdensity_fdm(M, cosmo, 'z', 0.5, 'corrections', true)

    %% ---- Defaults -------------------------------------------------------
    z           = 0;
    corrections = false;

    %% ---- Parse name-value pairs -----------------------------------------
    i = 1;
    while i <= numel(varargin)
        key = lower(varargin{i});
        val = varargin{i+1};
        switch key
            case 'z',           z           = val;
            case 'corrections', corrections = val;
            otherwise
                error('collapse_overdensity_fdm: unknown argument "%s".', key);
        end
        i = i + 2;
    end

    %% ---- Check cosmo has m22 --------------------------------------------
    if ~isfield(cosmo, 'm22') || isempty(cosmo.m22)
        error(['collapse_overdensity_fdm: cosmo.m22 (boson mass in units of ' ...
               '10^{-22} eV) must be set before calling this function.']);
    end
    m22 = cosmo.m22;   % m_a / (10^{-22} eV)

    %% ---- CDM baseline delta_c -------------------------------------------
    % EdS value: delta_c = (3/5)*(3*pi/2)^(2/3)
    delta_c_cdm = (3/5) * (3*pi/2)^(2/3);   % = 1.68647...

    if corrections
        Om_z = cosmo.Omega_m_z(z);
        if cosmo.flat
            delta_c_cdm = delta_c_cdm * Om_z^0.0055;
        else
            delta_c_cdm = delta_c_cdm * Om_z^0.0185;
        end
    end

    %% ---- FDM Jeans mass  M_J  (Du+2016, Eq. 14) ------------------------
    % M_J in h^{-1} M_sun
    Om_h2 = cosmo.Omh2;
    a1    = 3.4;
    M_J   = 1e8 * a1 * m22^(-1.5) * (Om_h2 / 0.14)^0.25 / cosmo.h;
    % Units: the formula gives h^{-1} M_sun directly when dividing by h.

    %% ---- Growth-factor ratio  G(M)  (Marsh 2016a / Du+2016, Eq. 11-13) -
    % Best-fit parameters from Du et al. (2016):
    %   {a1, a2, a3, a4, a5, a6} = {3.4, 1.0, 1.8, 0.5, 1.7, 0.9}
    % a1 is absorbed into M_J above; remaining params:
    a2 = 1.0;   a3 = 1.8;   a4 = 0.5;
    a5 = 1.7;   a6 = 0.9;

    x   = M(:) / M_J;   % dimensionless mass ratio  [N x 1]

    % Smooth transition function h_F(x)
    h_F = 0.5 * (1 - tanh(M_J * (x - a2)));
    % Note: the argument of tanh uses M_J as a dimensionless scale following
    % Du et al. (2016) Eq. 13; h_F -> 1 for x << a2, -> 0 for x >> a2.

    G   = h_F .* exp(a3 .* x.^(-a4)) + (1 - h_F) .* exp(a5 .* x.^(-a6));

    % Enforce G >= 1: below the CDM baseline makes no physical sense.
    G   = max(G, 1);

    %% ---- FDM collapse threshold -----------------------------------------
    delta_c = G * delta_c_cdm;   % [N x 1], broadcast scalar delta_c_cdm

end