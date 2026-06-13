    function cosmo = cosmology( name, user)

        if nargin < 1
            name = 'Planck18';
        end

        switch lower(name)
            case 'uchuu'        % also, planck 2020
                cosmo.Omega_m = 0.3089;
                cosmo.Omega_b = 0.0486;
                cosmo.h       = 0.6774;
                cosmo.ns      = 0.9667;
                cosmo.sigma8  = 0.8159;

            case 'planck18'
                cosmo.Omega_m = 0.3111;
                cosmo.Omega_b = 0.0490;
                cosmo.h       = 0.6766;
                cosmo.ns      = 0.9665;
                cosmo.sigma8  = 0.8102;

            case 'planck18-only'
                cosmo.Omega_m = 0.3153;
                cosmo.Omega_b = 0.0493;
                cosmo.h       = 0.6736;
                cosmo.ns      = 0.9649;
                cosmo.sigma8  = 0.8111;

            case 'planck15'
                cosmo.Omega_m = 0.3089;
                cosmo.Omega_b = 0.0486;
                cosmo.h       = 0.6774;
                cosmo.ns      = 0.9667;
                cosmo.sigma8  = 0.8159;

            case 'planck15-only'
                cosmo.Omega_m = 0.3080;
                cosmo.Omega_b = 0.0484;
                cosmo.h       = 0.6781;
                cosmo.ns      = 0.9677;
                cosmo.sigma8  = 0.8149;

            case 'planck13'
                cosmo.Omega_m = 0.3071;
                cosmo.Omega_b = 0.0483;
                cosmo.h       = 0.6777;
                cosmo.ns      = 0.9611;
                cosmo.sigma8  = 0.8288;

            case 'planck13-only'
                cosmo.Omega_m = 0.3175;
                cosmo.Omega_b = 0.0490;
                cosmo.h       = 0.6711;
                cosmo.ns      = 0.9624;
                cosmo.sigma8  = 0.8344;

            case 'wmap9'        % wmap9+bao
                cosmo.Omega_m = 0.2865;
                cosmo.Omega_b = 0.0463;
                cosmo.h       = 0.6932;
                cosmo.ns      = 0.9608;
                cosmo.sigma8  = 0.8200;

            case 'wmap9-only'        % wmap9
                cosmo.Omega_m = 0.2814;
                cosmo.Omega_b = 0.0464;
                cosmo.h       = 0.6970;
                cosmo.ns      = 0.9710;
                cosmo.sigma8  = 0.8200;

            case 'wmap9-ml'        
                cosmo.Omega_m = 0.2821;
                cosmo.Omega_b = 0.0461;
                cosmo.h       = 0.6970;
                cosmo.ns      = 0.9646;
                cosmo.sigma8  = 0.8170;

            case 'wmap7'       
                cosmo.Omega_m = 0.2743;
                cosmo.Omega_b = 0.0458;
                cosmo.h       = 0.7020;
                cosmo.ns      = 0.9680;
                cosmo.sigma8  = 0.8160;

           case 'wmap7-only'        
                cosmo.Omega_m = 0.2711;
                cosmo.Omega_b = 0.0451;
                cosmo.h       = 0.7030;
                cosmo.ns      = 0.9660;
                cosmo.sigma8  = 0.8090;

           case 'wmap7-ml'        
                cosmo.Omega_m = 0.2715;
                cosmo.Omega_b = 0.0455;
                cosmo.h       = 0.7040;
                cosmo.ns      = 0.9670;
                cosmo.sigma8  = 0.8100;

           case 'wmap5'        
                cosmo.Omega_m = 0.2732;
                cosmo.Omega_b = 0.0456;
                cosmo.h       = 0.7050;
                cosmo.ns      = 0.9600;
                cosmo.sigma8  = 0.8120;

           case 'wmap5-only'        
                cosmo.Omega_m = 0.2495;
                cosmo.Omega_b = 0.0432;
                cosmo.h       = 0.7240;
                cosmo.ns      = 0.9610;
                cosmo.sigma8  = 0.7870;

           case 'wmap5-ml'        
                cosmo.Omega_m = 0.2769;
                cosmo.Omega_b = 0.0459;
                cosmo.h       = 0.7020;
                cosmo.ns      = 0.9620;
                cosmo.sigma8  = 0.8170;

           case 'wmap3'        
                cosmo.Omega_m = 0.2342;
                cosmo.Omega_b = 0.0413;
                cosmo.h       = 0.7350;
                cosmo.ns      = 0.9510;
                cosmo.sigma8  = 0.7420;

           case 'wmap3-ml'        
                cosmo.Omega_m = 0.2370;
                cosmo.Omega_b = 0.0414;
                cosmo.h       = 0.7320;
                cosmo.ns      = 0.9540;
                cosmo.sigma8  = 0.7560;

           case 'wmap1'        
                cosmo.Omega_m = 0.2700;
                cosmo.Omega_b = 0.0463;
                cosmo.h       = 0.7200;
                cosmo.ns      = 0.9900;
                cosmo.sigma8  = 0.9000;

           case 'wmap1-ml'        
                cosmo.Omega_m = 0.3136;
                cosmo.Omega_b = 0.0497;
                cosmo.h       = 0.6800;
                cosmo.ns      = 0.9700;
                cosmo.sigma8  = 0.9000;

           case 'illustris'        
                cosmo.Omega_m = 0.2726;
                cosmo.Omega_b = 0.0456;
                cosmo.h       = 0.7040;
                cosmo.ns      = 0.9630;
                cosmo.sigma8  = 0.8090;

           case 'bolshoi'        
                cosmo.Omega_m = 0.2700;
                cosmo.Omega_b = 0.0469;
                cosmo.h       = 0.7000;
                cosmo.ns      = 0.9500;
                cosmo.sigma8  = 0.8200;

           case 'planck-mdark'        
                cosmo.Omega_m = 0.3070;
                cosmo.Omega_b = 0.0480;
                cosmo.h       = 0.6780;
                cosmo.ns      = 0.9600;
                cosmo.sigma8  = 0.8290;

            case 'millenium'        
                cosmo.Omega_m = 0.2500;
                cosmo.Omega_b = 0.0450;
                cosmo.h       = 0.7300;
                cosmo.ns      = 1.0000;
                cosmo.sigma8  = 0.9000;

            case 'eds'        
                cosmo.Omega_m = 1.0000;
                cosmo.Omega_b = 0.0000;
                cosmo.h       = 0.7000;
                cosmo.ns      = 1.0000;
                cosmo.sigma8  = 0.8200;

            case 'custom'
                cosmo = user;

            otherwise
                error('Unknown cosmology "%s"', name);
        end

        % --- Handle H0 vs h input for custom cosmology ---
        if isfield(cosmo, 'H0') && ~isfield(cosmo, 'h')
            cosmo.h = cosmo.H0 / 100.0;
        elseif ~isfield(cosmo, 'h') && ~isfield(cosmo, 'H0')
            error('Custom cosmology must provide either cosmo.h or cosmo.H0.');
        end
        
        cosmo.Omega_L   = 1 - cosmo.Omega_m;
        cosmo.rho_crit0 = 2.775e11;
        cosmo.rho_m0    = cosmo.rho_crit0 * cosmo.Omega_m; 

        if ~isfield(cosmo, 'transfer_model')
            cosmo.transfer_model = 'eh98';
        end
        % CAMB requires a file path — error early if missing
        if strcmpi(cosmo.transfer_model, 'camb')
            if ~isfield(cosmo, 'camb_transfer_file') || isempty(cosmo.camb_transfer_file)
                error(['transfer_model is set to ''camb'' but cosmo.camb_transfer_file is not set.\n' ...
                    'Please provide the path to the exported CAMB transfer function file.\n' ...
                    'Example: cosmo.camb_transfer_file = ''camb_transfer.mat'';']);
            end
        end


        % Attach all linear components
        cosmo = attach_linear_components( cosmo);   

        % update rho critical with z
        cosmo.rhocrit   = @(z) cosmo.rho_crit0 .* cosmo.E(z).^2;


        % % derived parameters:
        cosmo = derive_cosmo_params(cosmo);

        %
        cosmo.delta_c = @(z) collapse_overdensity('corrections', true, 'z', z, 'cosmo', cosmo);
        cosmo.nu      = @(M, z) cosmo.delta_c(z) ./ cosmo.sigmaM(M, z);
             
    end