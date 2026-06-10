function D = growth_with_radiation(z, cosmo)
    % EH98 integral at low-z, Gnedin+2011 Eq.5 at high-z, linear transition
    % in log(a) around z_switch = 10
    % Uses cosmo.E, cosmo.Omega_m, cosmo.a_eq from derive_cosmo_params

    % Note: Ez_D uses radiation as constant (not (1+z)^4) to avoid divergence
    % at high-z in the integral — same approach as Colossus
    Ez_D      = @(zp) sqrt(cosmo.Omega_m.*(1+zp).^3 + cosmo.Omega_L + (cosmo.Omega_r));
    integrand = @(zp) (1+zp) ./ Ez_D(zp).^3;

    % Gnedin+2011 Eq.5 — matter-radiation regime analytic approximation
    function Dg = gnedin(a_in)
        x  = a_in ./ cosmo.a_eq;
        t1 = sqrt(1 + x);
        t2 = 2*t1 + ((2/3) + x) .* log((t1-1)./(t1+1));
        Dg = a_in + (2/3)*cosmo.a_eq + ((cosmo.a_eq/(2*log(2)-3)) .* t2);
    end

    z = z(:);
    a = 1 ./ (1 + z);
    D = zeros(size(z));

    % Transition zone around z=10 in log-space: [zt2=5, zt1=20]
    zt1 = 20.0;  zt2 = 5.0;
    mask_low  = z <= zt2;
    mask_high = z >= zt1;
    mask_mid  = ~mask_low & ~mask_high;

    % --- Low-z: integral ---
    if any(mask_low)
        zl = z(mask_low);
        Dl = arrayfun(@(zi) integral(integrand, zi, Inf, ...
                     'RelTol',1e-4,'AbsTol',1e-10, 'ArrayValued', true), zl);
        D(mask_low) = 5/2 * cosmo.Omega_m .* Ez_D(zl) .* Dl;
    end

    % --- High-z: Gnedin ---
    if any(mask_high)
        D(mask_high) = gnedin(a(mask_high));
    end

    % --- Transition: linear interpolation in log(a) ---
    if any(mask_mid)
        zm  = z(mask_mid);
        Dm  = arrayfun(@(zi) integral(integrand, zi, Inf, ...
                      'RelTol',1e-4,'AbsTol', 0, 'ArrayValued', true), zm);
        D_int = 5/2 * cosmo.Omega_m .* Ez_D(zm) .* Dm;
        D_gne = gnedin(a(mask_mid));
        at1   = log(1/(zt1+1));
        at2   = log(1/(zt2+1));
        loga  = log(a(mask_mid));
        D(mask_mid) = (D_int.*(loga-at1) + D_gne.*(at2-loga)) ./ (at2-at1);
    end

    % --- Normalize to z=0 ---
    % 1+z_eq term is not needed as it will be normalized by D0 and it will
    % automatically cancel out here.
    D0 = 5/2 * cosmo.Omega_m * Ez_D(0) * ...
         quadgk(integrand, 0, Inf, 'RelTol',1e-4,'AbsTol', 0, 'ArrayValued', true);
    D  = D ./ D0;
end