function neff = neff( M,z,kappa, cosmo)
        % small step in ln R

        R_L = radius_from_mass( M,cosmo);
        R = kappa.*R_L;

        d = 5e-3;
        
        % two point difference method for d(ln sigma)/d(ln R)
        % Rp = R .* (1+d);
        % Rm = R ./ (1+d);
        % 
        % sig_p = cosmo.sigmaR(Rp, z);
        % sig_m = cosmo.sigmaR(Rm, z);
        % 
        % dlnsigma = log(sig_p ./ sig_m) ./ (2*log(1+d)); %-log(1-d));

        % Five-point central stencil for d(ln sigma)/d(ln R)
        sig1 = cosmo.sigmaR(R*(1+2*d), z);
        sig2 = cosmo.sigmaR(R*(1+ d),  z);
        sig3 = cosmo.sigmaR(R*(1- d),  z);
        sig4 = cosmo.sigmaR(R*(1-2*d), z);

        dlnsigma = (-log(sig1) + 8*log(sig2) - 8*log(sig3) + log(sig4)) ...
               ./ (12*d);
        neff = -3 - 2 * dlnsigma;

    end

