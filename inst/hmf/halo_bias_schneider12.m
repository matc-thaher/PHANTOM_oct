function b = halo_bias_schneider12(sigma, delta_c, M, Mhm)
% Schneider, Smith, Maccio & Moore (2012, MNRAS 424, 684) WDM halo bias.
%
% Applies the ST bias (Eq. 30) to WDM sigma(M,z), which must be
% pre-computed from the suppressed WDM power spectrum. The WDM physics
% enters entirely through sigma(M,z) — the functional form is identical
% to halo_bias_ST.
%
% For M > Mhm: b_WDM = b_ST(nu_WDM), indistinguishable from CDM bias.
% For M < Mhm: the bias upturn seen in simulations is driven by spurious
%   artificial haloes from particle discreteness (Section 4.2). These
%   bins are flagged with NaN and a warning is issued.
%
% Compute Mhm using: Mhm = halfmode_mass('WDM', m_WDM, cosmo)
%
% Inputs:
%   sigma   : rms variance sigma(M,z) from WDM power spectrum [same size as M]
%   delta_c : linear collapse threshold (default: EdS value ~1.686)
%   M       : halo mass [h^{-1} Msun], same size as sigma
%   Mhm     : half-mode mass [h^{-1} Msun], from halfmode_mass('WDM',...)
%
% Output:
%   b       : WDM linear halo bias; NaN where M < Mhm (spurious regime)
%
% Reference: Schneider, Smith, Maccio & Moore 2012, MNRAS 424, 684
%            arXiv:1112.0330   Section 4.2, Eq. 30

    if nargin < 2 || isempty(delta_c)
        delta_c = collapse_overdensity();
    end

    % ST bias applied to WDM sigma(M,z) — Section 4.2, Eq. 30
    b = halo_bias_SMT01(sigma, delta_c);

    % Mask spurious regime: M < Mhm is contaminated by artificial clumping.
    % Schneider+2012 Section 4.2 explicitly cautions against trusting bias
    % measurements below this scale.
    if nargin >= 4 && ~isempty(Mhm)
        spurious = M < Mhm;
        if any(spurious)
            warning('PHANTOM:halo_bias_schneider12', ...
                ['%d mass bin(s) fall below the half-mode mass ', ...
                 'Mhm = %.3e h^{-1}Msun. Bias here is contaminated by ', ...
                 'spurious artificial haloes (Schneider+2012, Sec. 4.2). ', ...
                 'These bins are set to NaN.'], ...
                sum(spurious), Mhm);
            b(spurious) = NaN;
        end
    end

end