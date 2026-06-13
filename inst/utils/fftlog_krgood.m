function krgood = fftlog_krgood(mu, q, dlnr, kr)
% FFTLOG_KRGOOD  Nearest low-ringing value of kr.
%
%   krgood = FFTLOG_KRGOOD(mu, q, dlnr, kr)
%
%   Adjusts kr to the nearest low-ringing value such that
%
%       (kr)^(-i*pi/dlnr) * U_mu(q + i*pi/dlnr)
%
%   is real. This is the low-ringing condition discussed by Hamilton.
%
%   Inputs
%   ------
%   mu    : Hankel order
%   q     : Bias parameter
%   dlnr  : Logarithmic spacing
%   kr    : Input kr value
%
%   Output
%   ------
%   krgood : Nearest low-ringing kr
%
%   Notes
%   -----
%   ln(krgood) lies within dlnr/2 of ln(kr).
%
%   References
%   ----------
%   Hamilton 2000, Section 6.
%   pyfftlog.krgood documentation.

    if kr <= 0
        error('kr must be positive.');
    end

    xp = (mu + 1 + q) / 2;
    xm = (mu + 1 - q) / 2;

    y = pi / (2*dlnr);

    zp = complex(xp, y);
    zm = complex(xm, y);

    % Phase of U_mu(q + i*pi/dlnr)
    argU = log(2) * (pi/dlnr) + imag(log(gamma(zp)) - log(gamma(zm)));

    % Snap ln(kr) to nearest allowed branch
    logkr = log(kr);
    logkr_good = dlnr * round((logkr + argU/pi) / dlnr) - argU/pi;

    krgood = exp(logkr_good);
end