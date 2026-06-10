function D = growth_ODE(z, cosmo)
    % Linder & Jenkins 2003 Eq.11 ODE for G = D/a
    % X(a) = Omega_m/a^3 / (delta_H2/H0^2)  [Eq. 9-10]
    % w(a) = w0 + wa*(1-a)                   [Eq. 4.1, simulation section]

    z      = z(:);
    a_eval = 1 ./ (1 + z);
    a_min  = min(min(a_eval)*0.99, 1e-3);
    a_max  = max(a_eval) * 1.01;

    % w(a) equation of state
    switch lower(cosmo.de_model)
        case 'w0'
            w_a = @(a) cosmo.w0 + zeros(size(a));
        case 'w0wa'
            w_a = @(a) cosmo.w0 + cosmo.wa .* (1 - a);
        otherwise
            error('growth_ODE: unsupported de_model "%s"', cosmo.de_model);
    end

    % dark energy density term: delta_H2/H0^2 = E^2(z) - Omega_m*(1+z)^3
    % X(a) = matter density / dark energy density  [LJ2003 Eq.10]
    de_density = @(a) cosmo.E(1/a - 1)^2 - cosmo.Omega_m / a^3;

    % Linder & Jenkins 2003 Eq.11: ODE for G = D/a
    function dydx = odes(a, y)
        wa  = w_a(a);
        Xa  = (cosmo.Omega_m / a^3) / de_density(a);
        t1  = (3.5 - 1.5*wa/(1+Xa)) / a;
        t2  =  1.5*(1-wa) / ((1+Xa) * a^2);
        dydx = [y(2); -t1*y(2) - t2*y(1)];
    end

    % Solve ODE — a must be sorted for ode45
    [a_sorted, idx] = sort(a_eval);
    opts     = odeset('RelTol', 1e-6, 'AbsTol', 1e-6);
    sol      = ode45(@odes, [a_min, a_max], [1.0; 0.0], opts);
    G        = deval(sol, a_sorted);
    D_sorted = G(1,:)' .* a_sorted;   % D = G * a

    % Restore original order
    D      = zeros(size(z));
    D(idx) = D_sorted;

    % Normalize to z=0
    G0 = deval(sol, 1.0);
    D  = D ./ (G0(1) * 1.0);
end