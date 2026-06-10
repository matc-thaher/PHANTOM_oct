    function T = T_EH98( k,cosmo)
        % Eisenstein & Hu 1998 transfer function (no-wiggle)
        % using equation 26, 28, 29, 30, 31 from paper

        Om = cosmo.Omega_m;
        Ob = cosmo.Omega_b;
        h  = cosmo.h;

        theta = 2.7255/2.7;

        s = 44.5*log(9.83./(Om*h.^2)) ./ sqrt(1 + 10*(Ob*h.^2).^(3/4));

        alpha_gamma = 1 - 0.328*log(431*(Om*h.^2))*(Ob/Om) ...
                    + 0.38*log(22.3*(Om*h.^2))*(Ob/Om).^2;

        Gamma_eff = Om*h * (alpha_gamma + ...
                   (1 - alpha_gamma) ./ (1 + (0.43*k*h*s).^4));

        q = k .* theta.^2 ./ Gamma_eff;

        L0 = log(2*exp(1) + 1.8*q);
        C0 = 14.2 + 731./(1 + 62.5*q);

        T0 = L0 ./ (L0 + C0.*q.^2);

        T = T0;
        
    end