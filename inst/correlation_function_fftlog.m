function xi = correlation_function_fftlog(R, z, cosmo, Pk_handle)
% CORRELATION_FUNCTION_FFTLOG  FFTLog-based xi(R,z).
%
%   Uses a Hamilton-style discrete logarithmic Hankel transform.
%
%   For the spherical Bessel j0 transform, the corresponding Hankel order is
%   mu = 1/2 because j0(x) = sqrt(pi/(2x)) J_{1/2}(x).
%
%   The transform is applied to
%       A(k) = k^(3/2) P(k,z)
%   and the result is converted back to
%       xi(r,z) = B(r) / (2*pi^2 * r^(3/2)).
%
%   Notes
%   -----
%   - This routine uses a low-ringing choice of kr = k0*r0.
%   - Accuracy depends on grid size, k-range, and how well aliasing is
%     controlled.
%   - Always validate against the direct integral method.
%
%   Reference
%   ---------
%   Hamilton 2000.

    Rreq = R(:);

    D2 = (cosmo.D(z) / cosmo.D(0)).^2;

    % FFTLog settings
    N    = 4096;
    kmin = 1e-6;
    kmax = 1e4;
    q    = 0.0;
    mu   = 0.5;

    % Logarithmic spacing
    L   = log(kmax / kmin);
    dln = L / N;

    % Logarithmic k-grid centered at k0
    jc  = (N + 1) / 2;
    j   = (1:N).';
    k0  = sqrt(kmin * kmax);
    k   = k0 .* exp((j - jc) * dln);

    % Choose central kr close to unity, then adjust to nearest low-ringing kr
    kr_guess = 1.0;
    kr       = fftlog_krgood(mu, q, dln, kr_guess);
    r0       = kr / k0;

    % Corresponding logarithmic r-grid
    r = r0 .* exp((j - jc) * dln);

    % Input sequence for Hankel transform
    Pk = Pk_handle(k) .* D2;

    % HIGH-kR DAMPING: suppress the high-k tail using the central output
    % scale r0 as representative.  This mirrors the Colossus exponential
    % suppression at kR > 1000 and reduces ringing at large output radii.
    kR_cut = 350.0;
    w = exp(-(k .* r0 ./ kR_cut).^2);
    A  = k.^(1.5) .* Pk .* w;

    % FFTLog transform
    B = fftlog_fht(A, mu, q, dln, kr);

    % Convert to xi(r)
    xi_grid = real(B) ./ (2*pi^2 .* r.^(1.5));

    % Interpolate to requested radii
    % xi = interp1(log(r), xi_grid, log(Rreq), 'pchip', 'extrap');
    % --- Interpolate to requested radii ----------------------------------
    % Guard: only interpolate within the reliable interior of the r-grid
    %   (5% margin from each edge to avoid pchip extrapolation artefacts)
    r_lo = r(max(1,  round(0.05*N)));
    r_hi = r(min(N,  round(0.95*N)));

    xi = zeros(size(Rreq));

    % Interior: pchip on log(r)
    in = (Rreq >= r_lo) & (Rreq <= r_hi);
    if any(in)
        xi(in) = interp1(log(r), xi_grid, log(Rreq(in)), 'pchip');
    end

    % Exterior: fall back to direct integral (avoids bad extrapolation)
    out = ~in;
    if any(out)
        xi(out) = correlation_function_integral(Rreq(out), z, cosmo, Pk_handle);
    end

    if isscalar(R)
        xi = xi(1);
    end
end