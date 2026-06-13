function f = multiplicity_Despali16(sigma, z, mdef, cosmo,  delta_c, ellipsoidal)
% Despali et al. (2016), MNRAS 456, 2486, Eq. 12
% Universal ST-like form rescaled to any SO definition via Delta/Delta_vir.
%
% INPUT:
%   Delta       : overdensity threshold (e.g. 200 for 200c)
%   Delta_vir   : virial overdensity at same redshift (from cosmo.Delta_vir(z))
%   delta_c     : (optional) collapse overdensity
%   ellipsoidal : (optional) use ellipsoidal halo finder params (default false)

    if nargin < 5 || isempty(delta_c), delta_c = collapse_overdensity('corrections', true, 'z', z, 'cosmo', cosmo); end
    if nargin < 6 || isempty(ellipsoidal), ellipsoidal = false; end

    Delta     = density_threshold(z, mdef, cosmo) ./ cosmo.rhocrit(z);
    Delta_vir = density_threshold(z, 'vir', cosmo) ./ cosmo.rhocrit(z);

    x = log10(Delta / Delta_vir);

    if ellipsoidal     % Appendix Eq.A1
        A =  -0.1768.*x + 0.3953;
        a =   0.3268.*x.^2 + 0.2125.*x + 0.7057;
        p =  -0.04570.*x.^2 + 0.1937.*x + 0.2206;
    else
        A =  -0.1362.*x + 0.3292;
        a =   0.4332.*x.^2 + 0.2263.*x + 0.7665;
        p =  -0.1151.*x.^2 + 0.2554.*x + 0.2488;
    end

    nu_p = a .* delta_c^2 ./ sigma.^2;
    f    = 2 .* A .* sqrt(nu_p ./ (2*pi)) .* exp(-0.5.*nu_p) .* (1 + nu_p.^-p);
end