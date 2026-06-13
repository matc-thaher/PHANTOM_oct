% =========================================================================
% LOCAL: mu(c) dispatcher — calls PHANTOM profile functions directly
% =========================================================================
function mu = profile_mu(c, profile_name)
% Calls PHANTOM src/profiles/ functions directly and grabs their fc output.
% A minimal dummy cosmo/inputs struct is used since only fc (which depends
% only on c) is needed — not the density or scale radius outputs.

    dummy_cosmo.rhocrit0 = 1.0;
    dummy_cosmo.E        = @(z) 1.0;
    dummy_cosmo.nu       = @(M, z) 1.0;   % nu=1 -> alpha_e=0.164, fixed shape

    mu = zeros(size(c));

    switch lower(profile_name)

        case 'nfw'
            for k = 1:numel(c)
                NFW   = NFW_analytcl_Profile(1.0, 1.0, c(k), [0.5]);
                mu(k) = NFW.f_c;
            end

        case 'hernquist'
            for k = 1:numel(c)
                [~, ~, ~, fc] = Hernquist_profile([0.5], 1.0, c(k), 0.0, ...
                                                   dummy_cosmo, 200);
                mu(k) = fc;
            end

        case 'einasto'
            for k = 1:numel(c)
                [~, ~, ~, fc] = Einasto_profile([0.5], 1.0, c(k), 0.0, ...
                                                  dummy_cosmo, 200);
                mu(k) = fc;
            end

        otherwise
            error('build_lhs_table:unknownProfile', ...
                  'Profile "%s" not supported. Add a case to profile_mu().', ...
                  profile_name);
    end
end