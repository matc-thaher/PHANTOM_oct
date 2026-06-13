    function alpha = alphaEff( z, cosmo)
    % alphaEff_Ishiyama  Effective growth exponent α_eff(z)
    %   Implements Eq. (B4) of Ishiyama+21:
    %       alpha_eff(z) = - d ln D(z) / d ln(1+z)
    %
    %   INPUT:
    %       z      : scalar or vector redshift
    %       cosmo  : struct with handle cosmo.D(z) giving the linear
    %                growth factor D(z), normalized e.g. D(0) = 1.
    %
    %   OUTPUT:
    %       alpha  : α_eff(z) at each input redshift

        % Small step in ln(1+z)-space for numerical derivative
        dz = 1e-3;

        zp = z + dz;
        zm = max(z - dz, 0);

        Dp = cosmo.D(zp);
        Dm = cosmo.D(zm);

        dln1pz = log((1+zp)/(1+zm));

        alpha = - (log(Dp) - log(Dm)) / dln1pz;
    end