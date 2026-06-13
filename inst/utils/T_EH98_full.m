function T = T_EH98_full(k, cosmo)
    % Eisenstein & Hu 1998 full transfer function (with baryons)
    % k in h/Mpc, returns dimensionless transfer function T(k)

    Om  = cosmo.Omega_m;
    Ob  = cosmo.Omega_b;
    h   = cosmo.h;
    k   = k*cosmo.h;
    
    Tcmb = 2.7255 / 2.7;   % CMB temperature ratio (theta_cmb)

    Omh2 = Om * h^2;
    Obh2 = Ob * h^2;
    fb   = Ob / Om;         % baryon fraction
    fc   = 1 - fb;          % CDM fraction

    %% --- Redshift of matter-radiation equality (Eq. 2 & 3) ---
    z_eq = 2.5e4 * Omh2 * Tcmb^(-4);
    k_eq = 7.46e-2 * Omh2 * Tcmb^(-2);   % [h/Mpc]

    %% --- Redshift of drag epoch (Eq. 4) ---
    b1   = 0.313 * Omh2^(-0.419) * (1 + 0.607 * Omh2^0.674);
    b2   = 0.238 * Omh2^0.223;
    z_d  = 1291 * Omh2^0.251 / (1 + 0.659 * Omh2^0.828) ...
           * (1 + (b1 * Obh2^b2));

    %% --- Sound horizon at drag epoch (Eq. 5 & 6) ---
    % R_eq = 31.5e3 * Obh2 * Tcmb^(-4) * (1000/z_eq);         % drag epoch
    % R_d  = 31.5e3 * Obh2 * Tcmb^(-4) * (1000/z_d);          % epoch of matter-radiation equality
    R_eq = 31.5 * Obh2 / Tcmb^4 / (z_eq / 1e3);   % matter-radiation equality
    R_d  = 31.5 * Obh2 / Tcmb^4 / (z_d / 1e3);    % drag epoch
    s    = 2/(3*k_eq) * sqrt(6/R_eq) * ...
           log((sqrt(1+R_d) + sqrt(R_d+R_eq)) / (1 + sqrt(R_eq)));

    %% --- Silk damping scale (Eq. 7) ---
    k_silk = 1.6 * Obh2^0.52 * Omh2^0.73 * (1 + (10.4 * Omh2)^(-0.95));  % [h/Mpc]

    %% --- CDM transfer function T_c (Eqs. 11, 12, 17, 18) ---
    a1    = (46.9*Omh2)^0.670 * (1 + (32.1*Omh2)^(-0.532));
    a2    = (12.0*Omh2)^0.424 * (1 + (45.0*Omh2)^(-0.582));
    alpha_c = a1^(-fb) * a2^(-fb^3);

    bb1   = 0.944 / (1 + (458*Omh2)^(-0.708));
    bb2   = (0.395*Omh2)^(-0.0266);
    beta_c = 1 / (1 + bb1*(fc^bb2 - 1));

    % T_tilde for CDM (Eq. 10, 17 & 20)
    T_tilde = @(k, alpha, beta) ...
        log(exp(1) + 1.8*beta.*k/(13.41*k_eq)) ./ ...
        (log(exp(1) + 1.8*beta.*k/(13.41*k_eq)) + ...
        ((14.2/alpha + 386./(1 + 69.9*(k/(13.41*k_eq)).^1.08)) .* (k/(13.41*k_eq)).^2));

    f  = 1 ./ (1 + (k*s/5.4).^4);
    Tc = f .* T_tilde(k, 1, beta_c) + (1-f) .* T_tilde(k, alpha_c, beta_c);

    %% --- Baryon transfer function T_b (Eqs. 14, 15, 21, 22, 23, 24) ---
    y    = (1 + z_eq) / (1 + z_d);
    G_y  = y * ((-6*sqrt(1+y)) + (2+3*y)*log((sqrt(1+y)+1)/(sqrt(1+y)-1)));
    alpha_b = 2.07 * k_eq * s * (1+R_d)^(-3/4) * G_y;

    beta_b   = 0.5 + fb + (3 - 2*fb) * sqrt((17.2*Omh2)^2 + 1);
    beta_node = 8.41 * Omh2^0.435;

    s_tilde  = s ./ (1 + (beta_node./(k*s)).^3).^(1/3);

    j0       = @(x) sin(x)./x;   % spherical Bessel j0
    Tb_term1 = T_tilde(k, 1, 1) ./ (1 + (k*s/5.2).^2);
    Tb_term2 = (alpha_b ./ (1 + (beta_b./(k*s)).^3)) .* exp(-(k/k_silk).^1.4);
    Tb = (Tb_term1 + Tb_term2) .* j0(k .* s_tilde);

    %% --- Total transfer function (Eq. 16) ---
    T = fb * Tb + fc * Tc;

end