function f = multiplicity_Watson13(sigma, z, Delta_m, Om_z)
% Watson et al. (2013), MNRAS 433, 1230
% Supports FOF (Delta_m = []) and SO (Delta_m in units of mean density).
%
% INPUT:
%   sigma   : sigma(M,z)
%   z       : redshift
%   Delta_m : overdensity w.r.t. mean (pass [] or omit for FOF mode)
%   Om_z    : Omega_m(z); required for SO mode at 0 < z < 6

    sigma = sigma(:);

    if nargin < 3 || isempty(Delta_m)
        % ---- FOF fit, Table 1, Eq. 12 ------------------------------------
        A     = 0.282;
        alpha = 2.163;
        beta  = 1.406;
        gamma = 1.210;
        f = A .* ((beta./sigma).^alpha + 1) .* exp(-gamma./sigma.^2);
        return
    end

    % ---- SO fit ----------------------------------------------------------
    Delta_178 = Delta_m / 178;

    if z == 0
        A = 0.194;  alpha = 1.805;  beta = 2.267;  gamma = 1.287;
    elseif z > 6
        A = 0.563;  alpha = 3.810;  beta = 0.874;  gamma = 1.453;
    else
        if nargin < 4 || isempty(Om_z)
            error('multiplicity_Watson13: Om_z required for SO mode at 0 < z <= 6.');
        end
        A     = Om_z .* (1.097 .* (1+z).^-3.216 + 0.074);
        alpha = Om_z .* (5.907 .* (1+z).^-3.058 + 2.349);
        beta  = Om_z .* (3.136 .* (1+z).^-3.599 + 2.344);
        gamma = 1.318;
    end
    


    f_178 = A .* ((beta./sigma).^alpha + 1) .* exp(-gamma./sigma.^2);

    % Delta correction, Eqs. 17-19 (note: Om_z needed here too)
    if nargin < 4 || isempty(Om_z), Om_z = 1; end
    C     = exp(0.023 .* (Delta_178 - 1));
    d     = -0.456 .* Om_z - 0.139;
    Gamma = C .* Delta_178.^d .* exp(0.072 .* (1 - Delta_178) ./ sigma.^2.130);
    f     = f_178 .* Gamma;
end