function xi = correlation_function_integral(R, z, cosmo, Pk_handle)
% CORRELATION_FUNCTION_INTEGRAL  Fast vectorized quadrature for xi(R,z).
%
%   xi = CORRELATION_FUNCTION_INTEGRAL(R, z, cosmo, Pk_handle)
%
%   Evaluates
%       xi(R,z) = (1 / 2*pi^2) * \int dlnk k^3 P(k,z) sin(kR)/(kR)
%
%   Uses pre-computed k-grid and vectorized operations for speed.

    R = R(:).';  % Row vector for broadcasting
    
    % Growth factor
    D2 = (cosmo.D(z) / cosmo.D(0))^2;
    
    % k-grid for integration (logarithmic spacing)
    % More points = better accuracy, but slower
    Nk = 4096;  % Adjust for speed vs accuracy trade-off
    lnk = linspace(log(1e-6), log(1e4), Nk)';  % Column vector
    k = exp(lnk);
    
    % Power spectrum at all k values (vectorized)
    Pk = Pk_handle(k) * D2;  % Column vector [Nk x 1]
    
    % Compute kR for all combinations: [Nk x NR]
    kR = k * R;  % Broadcasting: [Nk x 1] * [1 x NR] = [Nk x NR]
    
    % Spherical Bessel j0(x) = sin(x)/x, vectorized
    j0 = sin(kR) ./ kR;
    j0(kR < 1e-10) = 1.0;  % Limit as kR -> 0
    
    % --- HIGH-kR DAMPING  (key fix for large-r accuracy) ----------------
    % Suppress integrand where kR > kR_cut.  Modes with kR >> 1 oscillate
    % rapidly and their net contribution is zero; without damping they only
    % add trapezoidal noise that blows up near the zero-crossing.
    kR_cut = 350.0;
    w = exp(-(kR ./ kR_cut).^2);
    
    % Integrand: k^3 * P(k) * j0(kR)  [Nk x NR]
    % integrand = (k.^3 .* Pk) .* j0;  % Broadcasting
    integrand = (k.^3 .* Pk) .* (j0 .* w);
    
    % Integrate over ln(k) using trapezoidal rule
    % Each column is the integral for one R value
    xi = (1 / (2*pi^2)) * trapz(lnk, integrand, 1);  % Integrate along dim 1
    
    % Return as column vector
    xi = xi(:);
end