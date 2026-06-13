function b = halo_bias_SMT01(sigma, delta_c)
% Sheth, Mo & Tormen (2001, MNRAS 323, 1) halo bias.
% Derived from the peak-background split applied to the MOVING BARRIER
% of ellipsoidal collapse (Eq. 8 of SMT 2001). This is physically distinct
% from the ST99 bias (halo_bias_ST), which used a constant barrier.
%
% The moving barrier is:
%   B(sigma) = sqrt(a)*delta_c * [1 + beta*(sigma^2/(a*delta_c^2))^gamma]
%
% The resulting bias (Eq. 8, SMT 2001):
%   b(nu) = 1 + (1/(sqrt(a)*delta_c)) * [sqrt(a)*(a*nu^2)
%             + sqrt(a)*beta*(a*nu^2)^(1-gamma)
%             - (a*nu^2)^gamma / ((a*nu^2)^gamma + beta*(1-gamma)*(1-gamma/2))]
%
% Best-fit parameters (SMT 2001, below Eq. 8):
%   a = 0.707,  beta = 0.485,  gamma = 0.615
%
% NOTE: The ST99 function (halo_bias_ST) uses the peak-background split
% of the ST mass function with a CONSTANT barrier, yielding a simpler
% formula with p=0.3. That is a different (earlier) approximation.
% Use this function when you want the bias consistent with the full
% ellipsoidal collapse moving-barrier formulation.
%
% Reference: Sheth, Mo & Tormen 2001, MNRAS 323, 1   arXiv:astro-ph/9907024

    if nargin < 2 || isempty(delta_c)
        delta_c = collapse_overdensity();   % EdS value ~1.6865
    end

    % SMT 2001 moving-barrier parameters
    a     = 0.707;
    b     = 0.5;
    c     = 0.6;

    nu    = delta_c ./ sigma;        % peak height
    anu2  = a .* nu.^2;              % a*nu^2, appears repeatedly

    % Equation 8 of Sheth, Mo & Tormen (2001)
    term1 = sqrt(a) .* anu2;
    term2 = sqrt(a) .* b .* anu2.^(1 - c);
    term3 = anu2.^c ./ (anu2.^c + (b*(1-c)*(1-c/2)));

    b = 1 + ((1 ./ (sqrt(a) .* delta_c)) .* (term1 + term2 - term3));

end