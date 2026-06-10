function f = multiplicity_Schneider12(sigma, M, M_hm, z, cosmo, delta_c, variant)
% Schneider et al. (2012), MNRAS 424, 684
% WDM halo mass function for FIELD haloes.
%
% Variant 'A' (Eq. 27) — preferred, rms error < 5%:
%   f_WDM = f_ST(sigma_WDM) * (1 + M_hm/M)^-0.6
%   sigma must be computed from WDM linear P(k)
%
% Variant 'B' (Eq. 28) — simpler, slightly worse fit:
%   f_WDM = f_ST(sigma_CDM) * (1 + M_hm/M)^-1.16
%   sigma must be computed from CDM linear P(k)
%
% INPUT:
%   sigma   : sigma(M,z) — WDM P(k) for 'A', CDM P(k) for 'B'
%   M       : halo mass array [h^-1 M_sun]
%   M_hm    : half-mode mass [h^-1 M_sun]; use halfmode_mass()
%   delta_c : (optional) linear collapse threshold; default 1.686
%   variant : (optional) 'A' or 'B'; default 'A'
%
% NOTE: For subhalo mass functions use multiplicity_Lovell14 instead.

    if nargin < 6 || isempty(delta_c)
        delta_c = collapse_overdensity('corrections', true, 'z', z, 'cosmo', cosmo);
    end
    if nargin < 7 || isempty(variant)
        variant = 'A';
    end

    switch upper(variant)
        case 'A'
            exponent = 0.6;
        case 'B'
            exponent = 1.16;
        otherwise
            error('multiplicity_Schneider12: variant must be ''A'' or ''B''.');
    end

    f = multiplicity_ST(sigma, delta_c) .* (1 + M_hm ./ M).^(-exponent);
end