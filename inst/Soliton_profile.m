function rho = Soliton_profile(r, rho0, rc)
    % SOLITON_PROFILE  Computes the soliton (core) density profile for
    %                  ultra-light/fuzzy dark matter (ULDM/FDM) halos.
    %
    % This profile follows the fitting formula from Schive et al. (2014).
    %
    % USAGE:
    %   rho = soliton_profile(r, rho0, rc)
    %
    % INPUT:
    %   r    : Radial distance from halo center [kpc], scalar or array
    %   rho0 : Central (peak) soliton density [Msun/kpc^3]
    %   rc   : Soliton core radius [kpc], defined as the radius where the
    %          density drops to half its central value (half-density radius)
    %
    % OUTPUT:
    %   rho  : Soliton density at each radius r [Msun/kpc^3], same size as r
    %
    % FORMULA:
    %   rho(r) = rho0 * [ 1 + 0.091 * (r/rc)^2 ]^(-8)
    %
    %   The exponent -8 and coefficient 0.091 are empirical fitting parameters
    %   from Schive et al. (2014), calibrated against numerical simulations of
    %   the Schrodinger-Poisson equations.
    %
    
    % REFERENCE:
    %   Schive, H.-Y., Chiueh, T., & Broadhurst, T. (2014)
    %   "Cosmic Structure as the Quantum Interference of a Coherent Dark Wave"
    %   Nature Physics, 10, 496-499.
    %   https://doi.org/10.1038/nphys2996
    %
    % SEE ALSO:
    %   NFW_analytclPrfl

    rho = rho0 .* (1 + 0.091 * (r ./ rc).^2) .^ -8;
end