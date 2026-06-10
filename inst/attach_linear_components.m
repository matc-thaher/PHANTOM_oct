    function cosmo = attach_linear_components( cosmo)

        % % Expansion rate
        % cosmo.E  = @(z) sqrt(cosmo.Omega_m*(1+z).^3 + cosmo.Omega_L);

        % % Growth factor
        % cosmo.D = @(z) growth_factor_D( z,cosmo);
        % Growth factor — use extended version only when needed
        if isfield(cosmo,'relspecies') && cosmo.relspecies || ...
            isfield(cosmo,'flat') && ~cosmo.flat || ...
            isfield(cosmo,'de_model') && ~strcmpi(cosmo.de_model,'lambda')
            cosmo.D = @(z) growth_factor_D_ext(z, cosmo);
        else
            cosmo.D = @(z) growth_factor_D(z, cosmo);   % your original, untouched
        end

        % Transfer function
        % cosmo.T = @(k) T_EH98( k,cosmo);
        switch lower(cosmo.transfer_model)
            case 'eh98_full'
                cosmo.T = @(k) T_EH98_full(k, cosmo);

            case 'sugiyama95'
                cosmo.T = @(k) T_Sugiyama95(k, cosmo);

            case 'camb'
                 if ~isfield(cosmo,'python_exe') || isempty(cosmo.python_exe)
                    error(['For transfer_model = ''camb'', you must provide cosmo.python_exe ' ...
                        'with the full path to the Python executable, e.g. ' ...
                        '''C:\Python311\python.exe''.']);
                end
                if ~isfield(cosmo,'camb_minkh') || isempty(cosmo.camb_minkh)
                    cosmo.camb_minkh = 1e-4;
                end
                if ~isfield(cosmo,'camb_maxkh') || isempty(cosmo.camb_maxkh)
                    cosmo.camb_maxkh = 100.0;
                end
                if ~isfield(cosmo,'camb_npoints') || isempty(cosmo.camb_npoints)
                    cosmo.camb_npoints = 2000;
                end

                if isfield(cosmo,'As') && ~isempty(cosmo.As)
                    [k_tab, P_tab] = camb_power(cosmo, ...
                        'python_exe', cosmo.python_exe, ...
                        'minkh',      cosmo.camb_minkh, ...
                        'maxkh',      cosmo.camb_maxkh, ...
                        'npoints',    cosmo.camb_npoints, ...
                        'As',         cosmo.As);
                else
                    [k_tab, P_tab] = camb_power(cosmo, ...
                        'python_exe', cosmo.python_exe, ...
                        'minkh',      cosmo.camb_minkh, ...
                        'maxkh',      cosmo.camb_maxkh, ...
                        'npoints',    cosmo.camb_npoints);
                end

                % Store tables for later use if needed
                cosmo.k_camb = k_tab;
                cosmo.Pk0_camb_tab = P_tab;

                % Interpolated z=0 linear matter power spectrum
                cosmo.Pk0 = @(k) interp1(log(k_tab), P_tab, ...
                    min(max(log(k(:)), log(k_tab(1))), log(k_tab(end))), ...
                    'pchip');

                % Optional effective transfer function, only for plotting/diagnostics
                cosmo.T = @(k) sqrt(cosmo.Pk0(k) ./ k.^cosmo.ns);

            case 'viel05'
                % Warm dark matter transfer function (Viel+2005)
                % Required: cosmo.m_wdm_keV (thermal relic mass in keV), but we set
                % a default if missing (see below).
                if ~isfield(cosmo,'m_wdm_keV') || isempty(cosmo.m_wdm_keV)
                    warning(['transfer_model = ''viel05'': cosmo.m_wdm_keV not provided. ', ...
                        'Using default m_wdm_keV = 3.0 keV.']);
                    cosmo.m_wdm_keV = 3.0;
                end

                % Base CDM transfer for WDM (choose EH98 for now)
                if ~isfield(cosmo,'viel05_base') || isempty(cosmo.viel05_base)
                    cosmo.viel05_base = 'eh98';
                end

                switch lower(cosmo.viel05_base)
                    case 'eh98'
                        base_T = @(k) T_EH98(k, cosmo);

                    case 'eh98_full'
                        base_T = @(k) T_EH98_full(k, cosmo);

                    case 'sugiyama95'
                        base_T = @(k) T_Sugiyama95(k, cosmo);

                    otherwise
                        error('Unknown viel05_base ''%s''. Use eh98, eh98_full, sugiyama95, viel05, or bode01.', ...
                        cosmo.viel05_base);
                end

                % Suppression transfer T_WDM(k)
                supp_T  = @(k) T_wdm(k, cosmo.m_wdm_keV, cosmo, 'viel');

                % Total transfer: T_tot = T_WDM * T_base
                cosmo.T = @(k) T_total_suppressed(k, base_T, supp_T);

            case 'bode01'
                % Warm dark matter transfer function (Bode, Ostriker & Turok 2001)
                if ~isfield(cosmo,'m_wdm_keV') || isempty(cosmo.m_wdm_keV)
                    warning(['transfer_model = ''bode01'': cosmo.m_wdm_keV not provided. ', ...
                        'Using default m_wdm_keV = 3.0 keV.']);
                    cosmo.m_wdm_keV = 3.0;
                end

                % Base CDM transfer for WDM (choose EH98 for now)
                if ~isfield(cosmo,'bode01_base') || isempty(cosmo.bode01_base)
                    cosmo.bode01_base = 'eh98';
                end

                switch lower(cosmo.bode01_base)
                    case 'eh98'
                        base_T = @(k) T_EH98(k, cosmo);

                    case 'eh98_full'
                        base_T = @(k) T_EH98_full(k, cosmo);

                    case 'sugiyama95'
                        base_T = @(k) T_Sugiyama95(k, cosmo);

                    otherwise
                        error('Unknown bode01_base ''%s''. Use eh98, eh98_full, sugiyama95, viel05, or bode01.', ...
                        cosmo.bode01_base);
                end

                % Suppression transfer T_WDM(k)
                supp_T  = @(k) T_wdm(k, cosmo.m_wdm_keV, cosmo, 'bode');

                % Total transfer: T_tot = T_WDM * T_base
                cosmo.T = @(k) T_total_suppressed(k, base_T, supp_T);
            
            case 'schive25'
                % Fuzzy DM transfer (Schive 2025 / Hu+2000 form).
                % m22 = m / (1e-22 eV); default to 1 if not provided.
                if ~isfield(cosmo,'m22') || isempty(cosmo.m22)
                    warning(['transfer_model = ''schive25'': cosmo.m22 not provided. ', ...
                        'Using default m22 = 1.0 (m = 1e-22 eV).']);
                    cosmo.m22 = 1.0;
                end

                % Base CDM transfer model for Schive25.
                % User can choose via cosmo.schive25_base, with options:
                %   'eh98'      -> T_EH98
                %   'eh98_full' -> T_EH98_full
                %   'sugiyama95'-> T_Sugiyama95
                %
                % Default: 'eh98'.
                if ~isfield(cosmo,'schive25_base') || isempty(cosmo.schive25_base)
                    cosmo.schive25_base = 'eh98';
                end

                switch lower(cosmo.schive25_base)
                    case 'eh98'
                        base_T = @(k) T_EH98(k, cosmo);

                    case 'eh98_full'
                        base_T = @(k) T_EH98_full(k, cosmo);

                    case 'sugiyama95'
                        base_T = @(k) T_Sugiyama95(k, cosmo);

                    otherwise
                        error('Unknown schive25_base ''%s''. Use eh98, eh98_full, sugiyama95, viel05, or bode01.', ...
                        cosmo.schive25_base);
                end
                
                % Suppression transfer T_WDM(k)
                supp_T  = @(k) T_Schive25(k, cosmo.m22);

                % Total transfer: T(k) = T_FDM(k) * T_base(k),
                % with T_FDM from Schive/ Hu+2000.
                cosmo.T = @(k) T_total_suppressed(k, base_T, supp_T);

           case 'axioncamb'
                % AxionCAMB z=0 matter power spectrum from a precomputed .dat file
                % User must provide cosmo.axioncamb_file with full path to the .dat
                if ~isfield(cosmo,'axioncamb_file') || isempty(cosmo.axioncamb_file)
                    error(['For transfer_model = ''axioncamb'', provide cosmo.axioncamb_file ', ...
                        'pointing to params_axion_..._matterpower.dat.']);
                end

                data = load(cosmo.axioncamb_file);
                if size(data,2) < 2
                    error('axioncamb_file must have at least two columns [k, P(k)].');
                end

                k_ax  = data(:,1);   % [h/Mpc] or [1/Mpc], consistent with your k convention
                P_ax  = data(:,2);   % linear matter P(k) at z=0

                % Store tables
                cosmo.k_axioncamb     = k_ax;
                cosmo.Pk0_axioncamb   = P_ax;

                % Interpolated z=0 power spectrum (no extra sigma8 renormalization)
                cosmo.Pk0 = @(k) interp1(log(k_ax), P_ax, ...
                    min(max(log(k(:)), log(k_ax(1))), log(k_ax(end))), ...
                    'pchip');

                % Define an effective transfer function for diagnostics; not used in Pk0
                cosmo.T = @(k) sqrt(cosmo.Pk0(k) ./ k.^cosmo.ns);


           otherwise  % default: EH98 zero-baryon
                cosmo.T = @(k) T_EH98(k, cosmo);
        end

         % Default variance filter
        if ~isfield(cosmo,'variance_filter') || isempty(cosmo.variance_filter)
            cosmo.variance_filter = 'tophat';
        end

        % For sharpk/smoothk/vsmk, M = (4*pi/3) * rho_m * (c*R)^3
        % so R = (1/c) * (3M / 4*pi*rho_m)^(1/3)
        %
        % c defaults to 1 (top-hat / Gaussian — no calibration needed)
        % c = 2.5-2.7 for sharpk (Leo_2018, Ruderman_2026)
        % c = 3.0-3.7 for smoothk (Leo_2018)
        % c = 3.6     for vsmk    (Ruderman_2026)

        if ~isfield(cosmo, 'filter_c') || isempty(cosmo.filter_c)
            cosmo.filter_c = 1.0;
        end

        % % Unnormalized power spectrum
        % cosmo.Pk0_unnorm = @(k) k.^cosmo.ns .* cosmo.T(k).^2;
        % 
        % % Normalize P(k) to sigma8
        % A = normalize_A_to_sigma8(cosmo);
        % cosmo.Pk0 = @(k) A * cosmo.Pk0_unnorm(k);
        if ~strcmpi(cosmo.transfer_model,'camb') && ~strcmpi(cosmo.transfer_model,'axioncamb')
            cosmo.Pk0_unnorm = @(k) k.^cosmo.ns .* cosmo.T(k).^2;
            A = normalize_A_to_sigma8(cosmo);
            cosmo.Pk0 = @(k) A * cosmo.Pk0_unnorm(k);
        end

        cosmo.Pk  = @(k,z) cosmo.Pk0(k) .* (cosmo.D(z)/cosmo.D(0)).^2;

        % σ(R,z) and σ(M,z)
        cosmo.sigmaR = @(R,z,varargin) sigma_R(R, z, cosmo, varargin{:});
        cosmo.R_of_M = @(M) (1/cosmo.filter_c) * (3*M./(4*pi*cosmo.rho_m0)).^(1/3);
        cosmo.sigmaM = @(M,z,varargin) cosmo.sigmaR(cosmo.R_of_M(M), z, varargin{:});

        % n_eff and α_eff for Ishiyama21
        cosmo.neff     = @(M,z,kappa) neff(M,z,kappa,cosmo);
        cosmo.alphaEff = @(z) alphaEff(z,cosmo);

        if ~isfield(cosmo,'corr_method') || isempty(cosmo.corr_method)
            cosmo.corr_method = 'integral';
        end

        cosmo.correlationFunction = @(R,z,varargin) correlation_function(R, z, cosmo, varargin{:});

    end