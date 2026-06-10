function f = multiplicity_Diemer20(sigma, z, mdef, cosmo, delta_c)
% multiplicity_Diemer20  Splashback halo multiplicity f(sigma) for Diemer (2020)
%
%   f = multiplicity_Diemer20(sigma, z, mdef, cosmo)
%   f = multiplicity_Diemer20(sigma, z, mdef, cosmo, deltac_args)
%
%   This implements the Diemer (2020) splashback mass function for
%   dynamically measured splashback masses (mean or percentiles).
%
% INPUTS:
%   sigma       : rms linear fluctuation (array)
%   z           : redshift (scalar)
%   mdef        : splashback mass definition, e.g. 'sp-apr-mn' or 'sp-apr-p75'
%   cosmo       : PHANTOM cosmo struct with cosmo.collapseOverdensity(z, args)
%   deltac_args : (optional) struct of options for collapseOverdensity
%
% OUTPUT:
%   f           : multiplicity function f(sigma), same size as sigma

    if nargin < 4 || isempty(delta_c)
        delta_c = collapse_overdensity('corrections', true, 'z', z, 'cosmo', cosmo);
    end

    nu = delta_c ./ sigma(:);

    % Parameter sets from Diemer 2020, Table 2
    if strcmp(mdef, 'sp-apr-mn')
        A = 0.124399;
        a = 1.191457;
        b = 0.337871;
        c = 0.431710;

    elseif strncmp(mdef, 'sp-apr-p', 8)
        A  = 0.091878;
        a0 = 1.088267;
        b  = 0.242074;
        c0 = 0.445337;
        ap = 0.167465;
        cp = -0.068969;
        alpha = 1.756618;

        % Parse percentile XX from 'sp-apr-pXX'
        p_int = str2double(mdef(end-1:end));
        p     = p_int / 100.0;

        if p < 0.5 || p > 0.9
            error('multiplicity_Diemer20: percentile %d out of range [50..90].', p_int);
        end
        
        % equation (15)
        pprime = p.^alpha;
        a = a0 + ap * pprime;
        c = c0 + cp * pprime;

    else
        error('multiplicity_Diemer20: mass definition %s not supported.', mdef);
    end

    % f(nu) = A * [ (nu/b)^a + 1 ] * exp(-c nu^2 )
    f_vec = A .* ( (nu ./ b).^a + 1.0 ) .* exp(-c .* nu.^2);

    f = reshape(f_vec, size(sigma));
end