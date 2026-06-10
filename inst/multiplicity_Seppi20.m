function out = multiplicity_Seppi20(sigma, z, cosmo, delta_c, xoff, spin, int_over_sigma, int_over_xoff, int_over_spin)
% multiplicity_Seppi20  Seppi et al. (2020) 3D halo multiplicity model
%
%   out = multiplicity_Seppi20(sigma, z, cosmo, delta_c)
%   out = multiplicity_Seppi20(sigma, z, cosmo, delta_c, xoff, spin, ...
%                              int_over_sigma, int_over_xoff, int_over_spin)
%
%   Octave-compatible implementation of the Seppi+2020 3D abundance model
%   h(sigma, xoff, lambda), with optional marginalization over sigma, xoff,
%   and spin.
%
% INPUTS:
%   sigma           : array of sigma values
%   z               : redshift
%   cosmo           : cosmology struct
%   delta_c
%   xoff            : optional xoff grid, default logspace(-3.5,-0.3,50)
%   spin            : optional spin grid, default logspace(-3.5,-0.3,50)
%   int_over_sigma  : optional logical, default false
%   int_over_xoff   : optional logical, default true
%   int_over_spin   : optional logical, default true
%
% OUTPUT:
%   Depending on the integration flags, returns the corresponding
%   3D, 2D, or 1D distribution.

    if nargin < 4 || isempty(delta_c)
        delta_c = collapse_overdensity('corrections', true, 'z', z, 'cosmo', cosmo);
    end

    if nargin < 5 || isempty(xoff)
        xoff = logspace(-3.5, -0.3, 50);
    end
    if nargin < 6 || isempty(spin)
        spin = logspace(-3.5, -0.3, 50);
    end
    if nargin < 7 || isempty(int_over_sigma)
        int_over_sigma = false;
    end
    if nargin < 8 || isempty(int_over_xoff)
        int_over_xoff = true;
    end
    if nargin < 9 || isempty(int_over_spin)
        int_over_spin = true;
    end

    sigma = sigma(:);
    xoff  = xoff(:);
    spin  = spin(:);
    
    % Table A.4 and A.5
    zp   = 1 + z;
    A     = -22.004 * zp^(-0.0441);
    a     =  0.886  * zp^(-0.1611);
    q     =  2.285  * zp^( 0.0409);
    mu    = -3.326  * zp^(-0.1286);
    alpha =  5.623  * zp^( 0.1081);
    beta  = -0.391  * zp^(-0.3114);
    gamma =  3.024  * zp^( 0.0902);
    delta =  1.209  * zp^(-0.0768);
    e     = -1.105  * zp^( 0.6123);

    

    [SIGMA, XOFF, SPIN] = ndgrid(sigma, xoff, spin);
    nu = delta_c ./ SIGMA;

    ln10 = log(10.0);
    t1   = XOFF ./ 10.^(1.83 * mu);

    h_log10 = A ...
        + log10(sqrt(2.0 / pi)) ...
        + q .* log10(sqrt(a) .* nu) ...
        - (a ./ (2.0 * ln10)) .* nu.^2 ...
        + alpha .* log10(t1) ...
        - (1.0 / ln10) .* t1.^(0.05 * alpha) ...
        + gamma .* log10(SPIN ./ 10.^mu) ...
        - (1.0 / ln10) .* (t1 ./ SIGMA.^e).^beta .* (SPIN ./ 10.^mu).^delta;

    h = 10.^h_log10;

    g_xoff_spin  = seppi20_integrate_sigma(h, sigma);
    g_sigma_spin = seppi20_integrate_logxoff(h, xoff);
    g_sigma_xoff = seppi20_integrate_logspin(h, spin);

    f_xoff  = seppi20_integrate_sigma(g_sigma_xoff, sigma);
    f_spin  = seppi20_integrate_sigma(g_sigma_spin, sigma);
    f_sigma = seppi20_integrate_logxoff(g_sigma_xoff, xoff);

    if (~int_over_sigma) && (~int_over_xoff) && (~int_over_spin)
        out = h;
    elseif int_over_sigma && (~int_over_xoff) && (~int_over_spin)
        out = g_xoff_spin;
    elseif (~int_over_sigma) && int_over_xoff && (~int_over_spin)
        out = g_sigma_spin;
    elseif (~int_over_sigma) && (~int_over_xoff) && int_over_spin
        out = g_sigma_xoff;
    elseif int_over_sigma && int_over_xoff && (~int_over_spin)
        out = f_spin;
    elseif int_over_sigma && (~int_over_xoff) && int_over_spin
        out = f_xoff;
    elseif (~int_over_sigma) && int_over_xoff && int_over_spin
        out = f_sigma;
    else
        error('multiplicity_Seppi20: invalid integration flag combination.');
    end
end


function out = seppi20_integrate_sigma(arr, sigma)
    x = 1.0 ./ sigma(:);
    dims = ndims(arr);
    sz   = size(arr);

    if dims == 3
        out = zeros(sz(2), sz(3));
        i = 1;
        while i <= sz(2)
            j = 1;
            while j <= sz(3)
                if sz(1) == 1
                    out(i,j) = arr(1,i,j);
                else
                    out(i,j) = trapz(x, squeeze(arr(:,i,j)));
                end
                j = j + 1;
            end
            i = i + 1;
        end
    elseif dims == 2
        out = zeros(sz(2), 1);
        i = 1;
        while i <= sz(2)
            if sz(1) == 1
                out(i) = arr(1,i);
            else
                out(i) = trapz(x, arr(:,i));
            end
            i = i + 1;
        end
    else
        error('seppi20_integrate_sigma: unsupported array rank.');
    end
end


function out = seppi20_integrate_logxoff(arr, xoff)
    x = log10(xoff(:));
    dims = ndims(arr);
    sz   = size(arr);

    if dims == 3
        out = zeros(sz(1), sz(3));
        i = 1;
        while i <= sz(1)
            j = 1;
            while j <= sz(3)
                if sz(2) == 1
                    out(i,j) = arr(i,1,j);
                else
                    out(i,j) = trapz(x, squeeze(arr(i,:,j)));
                end
                j = j + 1;
            end
            i = i + 1;
        end
    elseif dims == 2
        out = zeros(sz(1), 1);
        i = 1;
        while i <= sz(1)
            if sz(2) == 1
                out(i) = arr(i,1);
            else
                out(i) = trapz(x, arr(i,:));
            end
            i = i + 1;
        end
    else
        error('seppi20_integrate_logxoff: unsupported array rank.');
    end
end

function out = seppi20_integrate_logspin(arr, spin)
    x  = log10(spin(:));
    sz = size(arr);

    if ndims(arr) ~= 3
        error('seppi20_integrate_logspin: expected rank-3 array.');
    end

    out = zeros(sz(1), sz(2));
    i = 1;
    while i <= sz(1)
        j = 1;
        while j <= sz(2)
            if sz(3) == 1
                out(i,j) = arr(i,j,1);
            else
                out(i,j) = trapz(x, squeeze(arr(i,j,:)));
            end
            j = j + 1;
        end
        i = i + 1;
    end
end

