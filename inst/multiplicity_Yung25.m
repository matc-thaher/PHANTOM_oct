% function f = multiplicity_Yung25(sigma, z)
% % Yung, Somerville & Iyer (2025), MNRAS 543, 3802
% % Eq. 2 + Table 1 — coefficients revised from Y24b to include z > 19 GUREFT data.
% % Valid range: 6 <= z <= 30, 6 <= log(M_h/M_sun) <= 13
% % Mass definition: Bryan & Norman (1998) virial, w.r.t. rho_crit
% %
% % NOTE: Colossus yung25 may carry Y24b coefficients — use Table 1 values below.
% 
%     % Table 1 coefficients [chi0, chi1, chi2]
%     % A_c = [0.21307778,  -0.01042236,   0.00013897];
%     % a_c = [0.94192066,   0.04453040,  -0.00202483];
%     % b_c = [3.27712602,  -0.01313422,   0.01027465];
%     % c_c = [1.15214631,   0.012866285, -0.00065572];
% 
%     A_c = [2.97165630e-01, -2.76808434e-03, -1.27528336e-04];
%     a_c = [1.65590338,      -5.50399410e-02, -1.63819807e-06];
%     b_c = [1.69700438,      -8.62801200e-02,  1.08082400e-02];
%     c_c = [1.16098576,       4.83463488e-03, -3.76272478e-04];
% 
%     % Polynomial evaluation in z
%     A = A_c(1) + A_c(2)*z + A_c(3)*z^2;
%     a = a_c(1) + a_c(2)*z + a_c(3)*z^2;
%     b = b_c(1) + b_c(2)*z + b_c(3)*z^2;
%     c = c_c(1) + c_c(2)*z + c_c(3)*z^2;
% 
%     sigma = sigma(:);
% 
%     % Eq. 2
%     f = A .* ((sigma./b).^(-a) + 1) .* exp(-c ./ sigma.^2);
% end

function [f, sigma_out, dlnSig_dlnM_out] = multiplicity_Yung25(sigma, z, variant, M_phys)
% multiplicity_Yung25   Yung, Somerville & Iyer (2025) multiplicity function.
%
% USAGE:
%   f = multiplicity_Yung25(sigma, z)
%       Default Colossus-like coefficients; uses sigma from PHANTOM pipeline.
%
%   f = multiplicity_Yung25(sigma, z, 'paper', M_phys)
%       Paper Table 1 coefficients + paper sigma(M) fit.
%       sigma_out and dlnSig_dlnM_out are returned for use in halo_mass_function.
%
% INPUT:
%   sigma   : sigma(M) array from PHANTOM pipeline
%   z       : redshift (scalar)
%   variant : 'colossus' (default) or 'paper'
%   M_phys  : physical mass [Msun], required when variant = 'paper'
%
% OUTPUT:
%   f               : multiplicity function values
%   sigma_out       : overridden sigma (paper variant only; else = input sigma)
%   dlnSig_dlnM_out : overridden derivative (paper variant only; else = [])

    if nargin < 3 || isempty(variant)
        variant = 'standard';
    end

    sigma_out       = sigma;   % default: pass through unchanged
    dlnSig_dlnM_out = [];      % default: no override

    switch lower(variant)

        case {'standard','default'}
            A_c = [2.97165630e-01, -2.76808434e-03, -1.27528336e-04];
            a_c = [1.65590338,     -5.50399410e-02, -1.63819807e-06];
            b_c = [1.69700438,     -8.62801200e-02,  1.08082400e-02];
            c_c = [1.16098576,      4.83463488e-03, -3.76272478e-04];

        case 'paper'
            if nargin < 4 || isempty(M_phys)
                error('multiplicity_Yung25: M_phys required for paper variant.');
            end
            A_c = [0.21307778,   -0.01042236,   0.00013897];
            a_c = [0.94192066,    0.04453040,  -0.00202483];
            b_c = [3.27712602,   -0.01313422,   0.01027465];
            c_c = [1.15214631,    0.012866285, -0.00065572];

            % Override sigma and its derivative using the paper fit
            [sigma_out, dlnSig_dlnM_out] = sigma_Yung25(M_phys(:));
            sigma = sigma_out(:);

        otherwise
            error('multiplicity_Yung25: unknown variant ''%s''. Use ''colossus'' or ''paper''.', variant);
    end

    A = A_c(1) + A_c(2)*z + A_c(3)*z^2;
    a = a_c(1) + a_c(2)*z + a_c(3)*z^2;
    b = b_c(1) + b_c(2)*z + b_c(3)*z^2;
    c = c_c(1) + c_c(2)*z + c_c(3)*z^2;

    sigma = sigma(:);
    f = A .* ((sigma ./ b).^(-a) + 1) .* exp(-c ./ sigma.^2);
end