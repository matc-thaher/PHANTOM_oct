% =========================================================================
function D = growth_integral(z, cosmo)
    % Heath 1977/ Peebles1980 / EH99 Eq.8 integral — works for non-flat, no radiation
    % D(z) = 5/2 * Om * E(z) * integral_{z}^{inf} (1+z')/E(z')^3 dz'
    Om = cosmo.Omega_m;

    % (1+z_eq) prefactor — cancels in D/D0 normalization but kept for correctness
    if cosmo.relspecies
        one_plus_zeq = 1.0 + cosmo.z_eq;
    else
        one_plus_zeq = 1.0;
    end

    % Integrand uses cosmo.E directly — no need to rebuild Ez
    integrand = @(zp) (1 + zp) ./ cosmo.E(zp).^3;


    z   = z(:);
    D   = zeros(size(z));
    for i = 1:numel(z)
        D(i) = quadgk(integrand, z(i), Inf, 'RelTol', 1e-4, 'AbsTol', 0, 'ArrayValued', true);
    end

    D   = (5/2) .* Om .* one_plus_zeq .* cosmo.E(z) .* D;
    D0  = (5/2) .* Om .* one_plus_zeq .* cosmo.E(0) .* quadgk(integrand, 0, Inf, 'RelTol', 1e-4, 'AbsTol', 0, 'ArrayValued', true);
    D   = D ./ D0;
end
