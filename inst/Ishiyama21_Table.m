  function P = Ishiyama21_Table( mode)

        switch lower(mode)

            case '500_all'
                P.kappa  = 1.83; P.a0 = 1.95; P.a1 = 1.17;
                P.b0     = 3.57; P.b1 = 0.91; P.cAlpha = 0.26;

            case '500_relaxed'
                P.kappa  = 0.38; P.a0 = 1.44; P.a1 = 3.41;
                P.b0     = 2.86; P.b1 = 2.99; P.cAlpha = 0.42;

            case '200c_all'
                P.kappa  = 1.19; P.a0 = 2.54; P.a1 = 1.33;
                P.b0     = 4.04; P.b1 = 1.21; P.cAlpha = 0.22;

            case '200c_relaxed'
                P.kappa  = 0.60; P.a0 = 2.14; P.a1 = 2.63;
                P.b0     = 1.69; P.b1 = 6.36; P.cAlpha = 0.37;

            case 'vir_all'
                P.kappa  = 1.64; P.a0 = 2.67; P.a1 = 1.23;
                P.b0     = 3.92; P.b1 = 1.30; P.cAlpha = -0.19;

            case 'vir_relaxed'
                P.kappa  = 1.22; P.a0 = 2.52; P.a1 = 1.87;
                P.b0     = 2.13; P.b1 = 4.19; P.cAlpha = -0.017;

            otherwise
                error('Unknown halo definition "%s".', mode);
        end

    end