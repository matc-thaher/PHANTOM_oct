  function D = growth_factor_D( z,cosmo)
        % Linear growth factor D(z) normalized to D(0)=1
        % Flat ΛCDM
        % written from appendix of eisenstein and Hu 98

            a = 1 ./ (1 + z);
            Om = cosmo.Omega_m;
            OL = cosmo.Omega_L;

            E2 = @(a) Om./a.^3 + OL;

            Om_a = @(a) Om ./ (a.^3 .* E2(a));
            OL_a = @(a) OL ./ E2(a);

            g = @(a) 5.*Om_a(a)./(2*( Om_a(a).^(4/7) - OL_a(a) + ...
                  (1 + Om_a(a)/2).*(1 + OL_a(a)/70) ));

            g0 = g(1.0);
            g_a = g(a);
            D = a .* g_a ./ g0;

        
    end