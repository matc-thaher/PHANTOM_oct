function rho_th = density_threshold(z, mdef, cosmo)
    [mdef_type, mdef_delta] = parse_mass_definition(mdef);

    rho_c = cosmo.rhocrit(z);   % make sure units are used consistently
    rho_m = cosmo.rhom(z);

    switch lower(mdef_type)
        case 'c'
            rho_th = mdef_delta .* rho_c;
        case 'm'
            rho_th = mdef_delta .* rho_m;
        case 'vir'
            rho_th = delta_vir_bn98(z, cosmo) .* rho_c;
        otherwise
            error('density_threshold: invalid mdef ''%s''.', mdef);
    end
end