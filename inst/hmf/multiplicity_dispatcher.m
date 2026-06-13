function f = multiplicity_dispatcher(model, sigma, varargin)
% multiplicity_dispatcher   Unified interface to all PHANTOM multiplicity functions.
%
% USAGE:
%   % --- CDM ---
%   f = multiplicity_dispatcher('PS',                sigma, delta_c)
%   f = multiplicity_dispatcher('ST',                sigma, delta_c)
%   f = multiplicity_dispatcher('Angulo12',          sigma, delta_c)
%   f = multiplicity_dispatcher('Bhattacharya11',    sigma, delta_c)
%   f = multiplicity_dispatcher('Tinker08',          sigma, Delta)
%   f = multiplicity_dispatcher('Watson13',          sigma, z)
%   f = multiplicity_dispatcher('Crocce10',          sigma, z)
%   f = multiplicity_dispatcher('Reed03',            sigma, delta_c)
%   f = multiplicity_dispatcher('Reed07',            sigma, delta_c)
%   f = multiplicity_dispatcher('Courtin11',         sigma, delta_c)
%   f = multiplicity_dispatcher('Despali16',         sigma, Delta_c, Delta_vir_c, delta_c)
%   f = multiplicity_dispatcher('Bocquet16',         sigma, varargin{:})
%   f = multiplicity_dispatcher('Comparat17',        sigma, varargin{:})
%   f = multiplicity_dispatcher('Diemer20',          sigma, varargin{:})
%   f = multiplicity_dispatcher('RodriguezPuebla16', sigma, varargin{:})
%   f = multiplicity_dispatcher('Seppi20',           sigma, varargin{:})
%   f = multiplicity_dispatcher('Seppi20M',          sigma, varargin{:})
%   f = multiplicity_dispatcher('Yung24',            sigma, varargin{:})
%   f = multiplicity_dispatcher('Yung25',            sigma, z)
%   f = multiplicity_dispatcher('FernandezGarcia26', sigma, M, z, mdef, cosmo)
%   f = multiplicity_dispatcher('Fiorilli26',        sigma, z, mdef, cosmo)
%   f = multiplicity_dispatcher('Fiorilli26',        sigma, z, mdef, cosmo, include_unbound)
%   % --- WDM ---
%   f = multiplicity_dispatcher('Schneider12',       sigma, M, M_hm, z, delta_c, variant)
%   f = multiplicity_dispatcher('Lovell14',          sigma, M, M_hm, z, delta_c)
%   % --- FDM ---
%   f = multiplicity_dispatcher('Schive16',          sigma, M, M_hm, z, delta_c, base_model)
%   f = multiplicity_dispatcher('Du17',              sigma, M, z, m22, cosmo)
%   f = multiplicity_dispatcher('Du17',              sigma, M, z, m22, cosmo, delta_c)
%
% INPUTS:
%   model    : string key for the desired multiplicity-function model (case-insensitive)
%   sigma    : rms linear fluctuation sigma(M, z); array of any size
%   varargin : model-specific arguments listed above
%
% OUTPUT:
%   f        : multiplicity function f(sigma), same size as sigma
%
% MODELS IMPLEMENTED:
%
%   CDM
%   -----------------------------------------------------------------------
%   PS              Press & Schechter (1974). Excursion-set spherical collapse.
%                   Alias: none. Extra arg: delta_c.
%
%   ST              Sheth & Tormen (1999). Ellipsoidal-collapse correction.
%                   Aliases: sheth99, sheth01. Extra arg: delta_c.
%
%   Angulo12        Angulo et al. (2012). Recalibrated ST form.
%                   Extra arg: delta_c.
%
%   Bhattacharya11  Bhattacharya et al. (2011). Redshift-dependent ST-type fit.
%                   Extra arg: delta_c.
%
%   Tinker08        Tinker et al. (2008). SO calibration; Delta wrt mean density.
%                   Extra arg: Delta.
%
%   Watson13        Watson et al. (2013). FOF and SO calibrations with explicit
%                   redshift dependence. Extra arg: z.
%
%   Crocce10        Crocce et al. (2010). FOF fit; power-law redshift evolution.
%                   Extra arg: z.
%
%   Reed03          Reed et al. (2003). First WDM-motivated modification to PS;
%                   retains delta_c input. Extra arg: delta_c.
%
%   Reed07          Reed et al. (2007). Improved CDM fit with non-Gaussian
%                   corrections. Extra arg: delta_c.
%
%   Courtin11       Courtin et al. (2011). Virial-overdensity calibration.
%                   Extra arg: delta_c.
%
%   Despali16       Despali et al. (2016). Spherical overdensity with virial-
%                   to-critical conversion. Extra args: Delta_c, Delta_vir_c, delta_c.
%
%   Bocquet16       Bocquet et al. (2016). Hydro-simulation fit including
%                   baryonic effects.
%
%   Comparat17      Comparat et al. (2017). MultiDark simulation calibration.
%
%   Diemer20        Diemer (2020). Peak-height and environment-dependent fit.
%
%   RodriguezPuebla16  Rodriguez-Puebla et al. (2016). Halo and subhalo fit
%                      with separate parametrizations.
%
%   Seppi20         Seppi et al. (2020). Full sample; SO definition.
%
%   Seppi20M        Seppi et al. (2020). Mass-binned variant.
%
%   Yung24          Yung et al. (2024). Predecessor calibration.
%
%   Yung25          Yung et al. (2025). Updated redshift-dependent fit.
%                   Extra arg: z.
%
%   FernandezGarcia26  Fernandez-Garcia et al. (2026). Extra args: M, z, mdef, cosmo.
%
%   Fiorilli26      Fiorilli, Ruiz, Sanchez & Esposito (2026). Evolution-mapping
%                   model encoding formation history (x_tilde) and local power-
%                   spectrum slope (n_eff). SO mass definitions 150m--1600m.
%                   Extra args: z, mdef, cosmo [, include_unbound].
%
%   WDM
%   -----------------------------------------------------------------------
%   Schneider12     Schneider et al. (2012). CDM baseline suppressed by a
%                   power-law factor below the half-mode mass M_hm.
%                   Extra args: M, M_hm, z, delta_c, variant.
%
%   Lovell14        Lovell et al. (2014). Ratio fit n_WDM/n_CDM calibrated
%                   to WDM N-body runs. Extra args: M, M_hm, z, delta_c.
%
%   FDM
%   -----------------------------------------------------------------------
%   Schive16        Schive et al. (2016). CDM baseline with a two-parameter
%                   suppression function F(M/M_0) from FDM simulations.
%                   Extra args: M, M_hm, z, delta_c, base_model.
%
%   Du17            Du, Behrens & Niemeyer (2017). Mass-dependent FDM
%                   collapse barrier via Sheth-Tormen with G(M)*delta_c.
%                   Extra args: M, z, m22, cosmo [, delta_c].
%
% NOTES:
%   - All model keys are case-insensitive.
%   - WDM and FDM models apply a suppression correction to a CDM baseline;
%     they do not define a universal fit independent of dark matter model.
%   - For Fiorilli26, mdef must follow the format 'NNNm' (e.g. '200m').
%     The calibrated range is 150m to 1600m.
%   - For Du17, delta_c defaults to 1.686 if not supplied.

    switch lower(model)

        % ----------------------------------------------------------
        % CDM models
        % ----------------------------------------------------------
        case 'ps'
            f = multiplicity_PS(sigma, varargin{:});
        case {'st', 'sheth99', 'sheth01'}
            f = multiplicity_ST(sigma, varargin{:});
        case 'angulo12'
            f = multiplicity_Angulo12(sigma, varargin{:});
        case 'bhattacharya11'
            f = multiplicity_Bhattacharya11(sigma, varargin{:});
        case 'tinker08'
            f = multiplicity_Tinker08(sigma, varargin{:});
        case 'watson13'
            f = multiplicity_Watson13(sigma, varargin{:});
        case 'crocce10'
            f = multiplicity_Crocce10(sigma, varargin{:});
        case 'reed03'
            f = multiplicity_Reed03(sigma, varargin{:});
        case 'reed07'
            f = multiplicity_Reed07(sigma, varargin{:});
        case 'courtin11'
            f = multiplicity_Courtin11(sigma, varargin{:});
        case 'despali16'
            f = multiplicity_Despali16(sigma, varargin{:});
        case 'bocquet16'
            f = multiplicity_Bocquet16(sigma, varargin{:});
        case 'comparat17'
            f = multiplicity_Comparat17(sigma, varargin{:});
        case 'diemer20'
            f = multiplicity_Diemer20(sigma, varargin{:});
        case {'rodriguezpuebla16','rodriguez-puebla16'}
            f = multiplicity_RodriguezPuebla16(sigma, varargin{:});
        case {'seppi20','seppi20full'}
            f = multiplicity_Seppi20(sigma, varargin{:});
        case 'seppi20m'
            f = multiplicity_Seppi20_mass(sigma, varargin{:});
        case 'yung24'
            f = multiplicity_Yung24(sigma, varargin{:});
        case 'yung25'
            f = multiplicity_Yung25(sigma, varargin{:});
        case 'fernandezgarcia26'
            f = multiplicity_FernandezGarcia26(sigma, varargin{:});
        case 'fiorilli26'
            f = multiplicity_Fiorilli26(sigma, varargin{:});

        % ----------------------------------------------------------
        % WDM models
        % ----------------------------------------------------------
        case 'schneider12'
            f = multiplicity_Schneider12(sigma, varargin{:});
        case 'lovell14'
            f = multiplicity_Lovell14(sigma, varargin{:});

        % ----------------------------------------------------------
        % FDM models
        % ----------------------------------------------------------
        case 'schive16'
            f = multiplicity_Schive16(sigma, varargin{:});
        case 'du17'
            % Du, Behrens & Niemeyer (2017), MNRAS 465, 941
            % Mass-dependent FDM barrier via Sheth-Tormen with G(M)*delta_c
            % varargin: {M, z, m22, cosmo} or {M, z, m22, cosmo, delta_c}
            f = multiplicity_Du17(sigma, varargin{:});

        otherwise
            error('multiplicity_dispatcher: unknown model ''%s''.', model);
    end
end