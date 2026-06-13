function atilde = fftlog_fht(a, mu, q, dln, kr)
% FFTLOG_FHT  Discrete biased Hankel transform of a log-spaced sequence.
%
%   atilde = FFTLOG_FHT(a, mu, q, dln, kr)
%
%   Implements the core FFTLog algorithm:
%     FFT -> multiply by u_m -> inverse FFT
%
%   where
%     u_m = (kr)^(-2*pi*i*m/(N*dln)) * U_mu(q + 2*pi*i*m/(N*dln))
%
%   with
%     U_mu(x) = 2^x * Gamma((mu+1+x)/2) / Gamma((mu+1-x)/2).
%
%   Reference
%   ---------
%   Hamilton 2000; pyfftlog documentation.

    a = a(:);
    N = numel(a);

    c = fft(a);

    m = fftlog_modes(N);
    eta = q + 2*pi*1i*m / (N*dln);
    u   = kr.^(-2*pi*1i*m/(N*dln)) .* U_mu(mu, eta);

    % Even-N Nyquist mode: use real part only, per Hamilton's prescription
    if mod(N,2) == 0
        ny = N/2 + 1;
        u(ny) = real(u(ny));
    end

    c = c .* u;
    atilde = ifft(c);
end