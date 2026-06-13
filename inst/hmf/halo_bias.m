function b = halo_bias(model, sigma, delta_c, varargin)
% HALO_BIAS_DISPATCHER  Unified halo bias interface for PHANTOM.
%
% Dispatches to the appropriate halo bias model given a string key.
% All models take sigma = sigma(M, z) and delta_c as primary inputs.
% Additional model-specific arguments are passed through varargin.
%
% Usage:
%   b = halo_bias_dispatcher(model, sigma)
%   b = halo_bias_dispatcher(model, sigma, delta_c)
%   b = halo_bias_dispatcher(model, sigma, delta_c, extra)
%
% Inputs:
%   model    : string key selecting the bias model (see table below)
%   sigma    : rms variance sigma(M, z), scalar or vector [dimensionless]
%   delta_c  : linear collapse threshold (default: EdS value ~1.686)
%   varargin : model-specific extra arguments (see each model below)
%
% Supported model keys:
%
%   CDM models
%   ----------
%   'cole89'        Press-Schechter / Cole-Kaiser (1989)
%                   b(nu) = 1 + (nu^2 - 1) / delta_c
%                   Extra args: none
%
%   'st'            Sheth & Tormen (1999) constant-barrier bias
%                   Extra args: none
%
%   'smt01'         Sheth, Mo & Tormen (2001) moving-barrier bias
%                   Extra args: none
%
%   'jing98'        Jing (1998, ApJ 503, L9)
%                   Extra args: cosmo struct (optional, for n_eff)
%
%   'seljak04'      Seljak & Warren (2004, MNRAS 355, 129)
%                   Extra args: cosmo struct (optional, for Om_m, ns, sigma8, h)
%
%   'tinker10'      Tinker et al. (2010, ApJ 724, 878)   [DEFAULT CDM]
%                   Extra args: Delta (overdensity, default 200),
%                               z     (redshift,    default 0),
%                               cosmo (struct,      default empty)
%                   Full call: dispatcher('tinker10', sigma, delta_c, Delta, z, cosmo)
%
%   'bhattacharya11' Bhattacharya et al. (2011, ApJ 732, 122)
%                   Extra args: z (redshift, default 0)
%
%   'comparat17'    Comparat et al. (2017, MNRAS 469, 4157)
%                   Extra args: none
%
%   'pillepich10'   Pillepich, Porciani & Hahn (2010, MNRAS 402, 191)
%                   Gaussian mode only via dispatcher.
%                   For non-Gaussian bias call halo_bias_pillepich10 directly:
%                     halo_bias_pillepich10(sigma, cosmo, 'nongaussian', k, f_NL, z)
%                   Extra args: cosmo struct (required for internal T(k) and D(z))
%
%   Beyond-CDM models (sigma must be pre-computed with suppressed transfer function)
%   ---------------------------------------------------------------------------------
%   'wdm'           WDM implicit bias (Schneider et al. 2012, MNRAS 424, 684)
%                   Applies Tinker+2010 to WDM sigma(M,z).
%                   Extra args: Delta (default 200), z (default 0), cosmo
%   Full call    : b = halo_bias_dispatcher('wdm', sigma, delta_c, M, Mhm)
%
%   'fdm'           FDM implicit bias (no dedicated formula in literature, 2026)
%                   Applies Tinker+2010 to FDM sigma(M,z) as an approximation.
%                   Results near the Jeans mass scale should be treated with caution.
%                   Extra args: Delta (default 200), z (default 0), cosmo
%
% Output:
%   b        : linear halo bias b(M, z), same size as sigma
%
% References:
%   Cole & Kaiser 1989, MNRAS 237, 1127
%   Sheth & Tormen 1999, MNRAS 308, 119            arXiv:astro-ph/9901122
%   Sheth, Mo & Tormen 2001, MNRAS 323, 1          arXiv:astro-ph/9907024
%   Jing 1998, ApJ 503, L9                         arXiv:astro-ph/9805202
%   Seljak & Warren 2004, MNRAS 355, 129           arXiv:astro-ph/0403698
%   Tinker et al. 2010, ApJ 724, 878               arXiv:1001.3162
%   Bhattacharya et al. 2011, ApJ 732, 122         arXiv:1005.2239
%   Comparat et al. 2017, MNRAS 469, 4157          arXiv:1611.06386
%   Pillepich, Porciani & Hahn 2010, MNRAS 402, 191  arXiv:0811.4203
%   Schneider et al. 2012, MNRAS 424, 684          arXiv:1112.0330

    if nargin < 3 || isempty(delta_c)
        delta_c = collapse_overdensity();
    end

    switch lower(strtrim(model))

        % ----------------------------------------------------------
        % Cole & Kaiser (1989) / Press-Schechter bias
        % ----------------------------------------------------------
        case 'cole89'
            b = halo_bias_PS(sigma, delta_c);

        % ----------------------------------------------------------
        % Sheth & Tormen (1999) — constant barrier
        % ----------------------------------------------------------
        case 'st'
            b = halo_bias_ST(sigma, delta_c);

        % ----------------------------------------------------------
        % Sheth, Mo & Tormen (2001) — moving barrier
        % ----------------------------------------------------------
        case 'smt01'
            b = halo_bias_SMT01(sigma, delta_c);

        % ----------------------------------------------------------
        % Jing (1998)
        % ----------------------------------------------------------
        case 'jing98'
            if ~isempty(varargin)
                cosmo = varargin{1};
            else
                cosmo = struct();
            end
            b = halo_bias_jing98(sigma, delta_c, cosmo);

        % ----------------------------------------------------------
        % Seljak & Warren (2004)
        % ----------------------------------------------------------
        case 'seljak04'
            if ~isempty(varargin)
                cosmo = varargin{1};
            else
                cosmo = [];
            end
            b = halo_bias_seljak04(sigma, delta_c, cosmo);

        % ----------------------------------------------------------
        % Tinker et al. (2010)   [default CDM]
        % ----------------------------------------------------------
        case 'tinker10'
            Delta = 200;
            z     = 0;
            cosmo = struct();
            if length(varargin) >= 1 && ~isempty(varargin{1}), Delta = varargin{1}; end
            if length(varargin) >= 2 && ~isempty(varargin{2}), z     = varargin{2}; end
            if length(varargin) >= 3 && ~isempty(varargin{3}), cosmo = varargin{3}; end
            b = halo_bias_Tinker10(sigma, delta_c, Delta, z, cosmo);

        % ----------------------------------------------------------
        % Bhattacharya et al. (2011)
        % ----------------------------------------------------------
        case 'bhattacharya11'
            z = 0;
            if ~isempty(varargin) && ~isempty(varargin{1}), z = varargin{1}; end
            b = halo_bias_bhattacharya11(sigma, delta_c, z);

        % ----------------------------------------------------------
        % Comparat et al. (2017)
        % ----------------------------------------------------------
        case 'comparat17'
            b = halo_bias_comparat17(sigma, delta_c);

        % ----------------------------------------------------------
        % Pillepich, Porciani & Hahn (2010) — Gaussian mode only.
        % Non-Gaussian mode requires k and f_NL; call the function
        % directly: halo_bias_pillepich10(sigma, cosmo, 'nongaussian', k, f_NL, z)
        % ----------------------------------------------------------
        case 'pillepich10'
            if ~isempty(varargin)
                cosmo = varargin{1};
            else
                error('PHANTOM:halo_bias_dispatcher', ...
                    ['pillepich10 requires a cosmo struct as the first extra ' ...
                     'argument: dispatcher(''pillepich10'', sigma, delta_c, cosmo).']);
            end
            b = halo_bias_pillepich10(sigma, cosmo, 'gaussian');

        % ----------------------------------------------------------
        % WDM bias — Schneider, Smith, Maccio & Moore (2012, MNRAS 424, 684).
        % ST bias (Eq. 30) applied to WDM sigma(M,z), pre-computed from
        % the suppressed WDM transfer function. M and Mhm are required to
        % mask the spurious halo regime below the half-mode mass scale.
        % Compute Mhm via: halfmode_mass('WDM', m_WDM, cosmo)
        % Extra args:
        %   M   : halo mass array [h^{-1} Msun], same size as sigma
        %   Mhm : half-mode mass [h^{-1} Msun]
        % ----------------------------------------------------------
        case 'wdm'
            M   = [];
            Mhm = [];
            if length(varargin) >= 1 && ~isempty(varargin{1}), M   = varargin{1}; end
            if length(varargin) >= 2 && ~isempty(varargin{2}), Mhm = varargin{2}; end
            b = halo_bias_schneider12(sigma, delta_c, M, Mhm);

        % ----------------------------------------------------------
        % FDM implicit bias — no published FDM-specific formula (2026).
        % Tinker+2010 applied to FDM sigma(M,z) as an approximation.
        % sigma must be pre-computed with the FDM transfer function.
        % Results near the Jeans mass scale should be treated with caution.
        % ----------------------------------------------------------
        case 'fdm'
            Delta = 200;
            z     = 0;
            cosmo = struct();
            if length(varargin) >= 1 && ~isempty(varargin{1}), Delta = varargin{1}; end
            if length(varargin) >= 2 && ~isempty(varargin{2}), z     = varargin{2}; end
            if length(varargin) >= 3 && ~isempty(varargin{3}), cosmo = varargin{3}; end

            % FDM collapse threshold: mass-dependent barrier from Du et al. (2016),
            % Eq. (11-13), with the fitting function for G(M) from Marsh (2016a).
            % delta_c^fdm(M,z) = G(M) * delta_c^cdm(z), where G -> 1 for M >> M_J
            % (recovering the CDM limit) and G >> 1 near the Jeans mass.
            % No dedicated FDM bias fitting formula exists as of 2026; Tinker+2010
            % is applied to FDM sigma(M,z) with delta_c^fdm as a working approximation.
            % Marsh & Pop (2015, MNRAS 451, 2479) discuss the modified barrier;
            % see also Du, Behrens & Niemeyer (2016, MNRAS 465, 941).
            % Issue a warning when sigma or M is near the Jeans scale, where the
            % barrier remapping of Sheth et al. (2001) has not been recalibrated
            % for FDM and results carry systematic uncertainty of order unity.
            if isfield(cosmo, 'm22') && ~isempty(cosmo.m22)
                % M must be available here — pass it in varargin or retrieve from sigma
                delta_c = collapse_overdensity_fdm(M, cosmo, 'z', z, 'corrections', true);
                Om_h2 = cosmo.Omh2;
                M_J   = 1e8 * 3.4 * cosmo.m22^(-1.5) * (Om_h2 / 0.14)^0.25 / cosmo.h;
                if any(M(:) < 10 * M_J)
                    warning('PHANTOM:halo_bias_fdm_jeans', ...
                         ['FDM halo bias: one or more mass values fall within a decade ' ...
                            'of the Jeans mass M_J = %.2e h^{-1} M_sun (m22 = %.2f). ' ...
                            'The Tinker+2010 bias with delta_c^fdm is an approximation; ' ...
                            'treat these bias values as order-of-magnitude estimates.'], ...
                            M_J, cosmo.m22);
                end
            else
                warning('PHANTOM:halo_bias_dispatcher', ...
                    ['FDM halo bias: cosmo.m22 not set. Using CDM delta_c = 1.686 ' ...
                    'as fallback. Set cosmo.m22 to activate the mass-dependent ' ...
                    'FDM collapse barrier (Du+2016).']);
                delta_c = collapse_overdensity();   % CDM fallback
            end

            b = halo_bias_SMT01(sigma, delta_c, Delta, z, cosmo);

        % ----------------------------------------------------------
        % Unknown model key
        % ----------------------------------------------------------
        otherwise
            error('PHANTOM:halo_bias_dispatcher', ...
                ['Unknown halo bias model: ''%s''.\n' ...
                 'Valid keys: cole89, st, smt01, jing98, seljak04, tinker10, ' ...
                 'bhattacharya11, comparat17, pillepich10, wdm, fdm.'], model);
    end

end
