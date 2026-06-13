function D = growth_factor_D_ext(z, cosmo)
    % Extended linear growth factor D(z), normalized so D(z)/D(0) gives growth
    % Supports: flat/non-flat LCDM, radiation correction, w0/w0wa dark energy
    % Falls back to EH98 approximation when possible (same as growth_factor_D)
    %
    % For non-LCDM dark energy (w0wa), solves ODE (Linder & Jenkins 2003, Eq 11)
    % For radiation, uses Gnedin+2011 Eq 5 at high-z with linear transition

    % --- Determine which regime to use ---
    has_radiation = isfield(cosmo, 'relspecies') && cosmo.relspecies && cosmo.Omega_r > 0;
    is_flat       = ~isfield(cosmo, 'flat') || cosmo.flat;
    is_lcdm       = ~isfield(cosmo, 'de_model') || strcmpi(cosmo.de_model, 'lambda');

    if is_lcdm && is_flat && ~has_radiation
        % --- Regime 1: flat LCDM, no radiation — use EH98 analytic (your original) ---
        D = growth_factor_D(z, cosmo);

    elseif is_lcdm && ~has_radiation
        % --- Regime 2: non-flat LCDM, no radiation — use Heath/EH98 integral ---
        D = growth_integral(z, cosmo);

    elseif is_lcdm && has_radiation
        % --- Regime 3: flat LCDM + radiation — EH98 at low-z, Gnedin at high-z ---
        D = growth_with_radiation(z, cosmo);

    else
        % --- Regime 4: non-LCDM dark energy — ODE (Linder & Jenkins 2003) ---
        D = growth_ODE(z, cosmo);
    end

end

