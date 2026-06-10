function f = multiplicity_Tinker08(sigma, Delta, z)
% multiplicity_Tinker08   Halo multiplicity function of Tinker et al. (2008).
%
% Calibrated for spherical overdensity (SO) masses over 200 <= Delta_m <= 3200.
% Includes explicit redshift dependence (Eqs. 3-8 of the paper).
%
% USAGE:
%   f = multiplicity_Tinker08(sigma, Delta)
%   f = multiplicity_Tinker08(sigma, Delta, z)
%
% INPUT:
%   sigma : rms variance sigma(M, z), column vector [dimensionless]
%   Delta : overdensity w.r.t. mean matter density (200 <= Delta <= 3200)
%           NOTE: Delta must be w.r.t. MEAN density, not critical.
%           If your halo finder uses Delta_c, convert first:
%             Delta_m = Delta_c * (Omega_m_z / 1.0)
%   z     : redshift (default 0)
%
% OUTPUT:
%   f     : multiplicity function f(sigma), same size as sigma
%
% Reference: Tinker et al. 2008, ApJ 688, 709   arXiv:0803.2706
%            Eqs. 2-8, Table 2

    if nargin < 2 || isempty(Delta), Delta = 200; end
    if nargin < 3 || isempty(z),     z     = 0.0; end

    % --- Table 2 calibration nodes ---
    Delta_ref = [200, 300, 400, 600, 800, 1200, 1600, 2400, 3200];
    A0_ref    = [0.186, 0.200, 0.212, 0.218, 0.248, 0.255, 0.260, 0.260, 0.260];
    a0_ref    = [1.47,  1.52,  1.56,  1.61,  1.87,  2.13,  2.30,  2.53,  2.66 ];
    b0_ref    = [2.57,  2.25,  2.05,  1.87,  1.59,  1.51,  1.46,  1.44,  1.41 ];
    c0_ref    = [1.19,  1.27,  1.34,  1.45,  1.58,  1.80,  1.97,  2.24,  2.44 ];

    % --- Range check ---
    if Delta < Delta_ref(1)
        error('multiplicity_Tinker08: Delta = %g too small, minimum is %d.', ...
              Delta, Delta_ref(1));
    end
    if Delta > Delta_ref(end)
        error('multiplicity_Tinker08: Delta = %g too large, maximum is %d.', ...
              Delta, Delta_ref(end));
    end

    % --- Interpolate parameters at requested Delta (log spacing, Eq. 3) ---
    A0 = interp1(log(Delta_ref), A0_ref, log(Delta), 'linear', 'extrap');
    a0 = interp1(log(Delta_ref), a0_ref, log(Delta), 'linear', 'extrap');
    b0 = interp1(log(Delta_ref), b0_ref, log(Delta), 'linear', 'extrap');
    c0 = interp1(log(Delta_ref), c0_ref, log(Delta), 'linear', 'extrap');

    % --- Redshift evolution (Eqs. 5-8) ---
    alpha = 10 .^ ( -(0.75 ./ log10(Delta / 75.0)).^1.2 );
    A     = A0 * (1.0 + z)^(-0.14);
    a     = a0 * (1.0 + z)^(-0.06);
    b     = b0 * (1.0 + z)^(-alpha);
    c     = c0;                          % no redshift dependence on c

    % --- Multiplicity function (Eq. 2) ---
    sigma = sigma(:);    % force column vector
    f     = A .* ( (sigma ./ b).^(-a) + 1.0 ) .* exp(-c ./ sigma.^2);

end
