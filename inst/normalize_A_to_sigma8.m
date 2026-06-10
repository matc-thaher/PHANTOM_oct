 function A = normalize_A_to_sigma8( cosmo)
    % Normalization amplitude A so that sigma8 matches

        R8 = 8.0; % Mpc/h

        Pk0_unit = @(k) cosmo.Pk0_unnorm(k);

        s2_unit = sigma_R2_given_Pk(R8, 0, cosmo, Pk0_unit);
        sigma_unit = sqrt(s2_unit);

        A = (cosmo.sigma8 / sigma_unit)^2;
    
 end