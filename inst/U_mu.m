function U = U_mu(mu, x)
% U_MU  FFTLog kernel factor U_mu(x).
%
%   U = U_MU(mu, x)
%
%   Computes
%       U_mu(x) = 2^x * Gamma((mu + 1 + x)/2) / Gamma((mu + 1 - x)/2)
%
%   for complex x.
%
%   Reference
%   ---------
%   Hamilton 2000.

    U = 2.^x .* gamma((mu + 1 + x)/2) ./ gamma((mu + 1 - x)/2);
end