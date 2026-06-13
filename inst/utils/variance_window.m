function W = variance_window(kR, filter_name, filter_params)
    % Fourier-space window functions for variance integral sigma^2(R,z)
    %
    % kR           : k*R array  (= k/k_M since k_M = 1/R)
    % filter_name  : 'tophat' | 'gaussian' | 'sharpk' | 'smoothk' | 'vsmk'
    % filter_params: optional struct with fields:
    %   .beta   - slope for smoothk           (default: 4,   Leo_2018)
    %   .beta1  - low-k  slope for vsmk       (default: 4.8, Ruderman_2026)
    %   .beta2  - high-k slope for vsmk       (default: 3.6, Ruderman_2026)
    %   .mu     - transition scale for vsmk   (default: 2.1, Ruderman_2026)
    %   .delta  - transition sharpness vsmk   (default: 12,  Ruderman_2026)
    %
    % NOTE: the mass calibration parameter c from Ruderman_2026 does NOT
    % enter the window function W(k,R). It only enters the mass-radius
    % relation M = (4*pi/3)*rho_m*(c*R)^3 and should be applied when
    % converting R to M, not here.
    %
    % References:
    %   tophat   : standard, e.g. Diemer_2018
    %   gaussian : standard, e.g. Diemer_2018
    %   sharpk   : Heaviside Theta(1-kR), Diemer_2018; Bond_1991
    %   smoothk  : Leo_2018, JCAP 2018, Eq. 9
    %   vsmk     : Ruderman_2026, arXiv:2602.01320, Eq. 3.6-3.7

    if nargin < 2 || isempty(filter_name)
        filter_name = 'tophat';
    end
    if nargin < 3 || isempty(filter_params)
        filter_params = struct();
    end

    % Default parameters
    beta  = 4.0;    % smoothk default (Leo_2018)
    beta1 = 4.8;    % vsmk small-k slope (Ruderman_2026, calibrated on WDM)
    beta2 = 3.6;    % vsmk large-k slope (Ruderman_2026, calibrated on DAO)
    mu    = 2.1;    % vsmk transition scale (Ruderman_2026, fixed across models)
    delta = 12.0;   % vsmk transition sharpness (Ruderman_2026, fixed across models)

    if isfield(filter_params, 'beta'),  beta  = filter_params.beta;  end
    if isfield(filter_params, 'beta1'), beta1 = filter_params.beta1; end
    if isfield(filter_params, 'beta2'), beta2 = filter_params.beta2; end
    if isfield(filter_params, 'mu'),    mu    = filter_params.mu;    end
    if isfield(filter_params, 'delta'), delta = filter_params.delta; end

    x = kR;  % x = k*R = k/k_M

    switch lower(filter_name)

        case {'tophat', 'top-hat'}
            % Real-space top-hat, Fourier transform
            % W = 3(sin(x) - x*cos(x)) / x^3
            W = 3 * (sin(x) - x .* cos(x)) ./ (x.^3);
            W(x == 0) = 1.0;

        case {'gaussian', 'gauss'}
            % Gaussian filter: W = exp(-x^2 / 2)
            W = exp(-0.5 * x.^2);

        case {'sharpk', 'sharp-k'}
            % Sharp-k: Heaviside Theta(1 - kR)
            % W = 1 if kR <= 1, else 0
            % This is Theta(1 - kR) as in COLOSSUS / Bond_1991
            W = double(x <= 1.0);

        case {'smoothk', 'smooth-k'}
            % Smooth-k filter (Leo_2018, Eq. 9):
            % W = 1 / (1 + (k/k_M)^beta)
            % where x = k/k_M = k*R
            W = 1.0 ./ (1.0 + x.^beta);

        case {'vsmk', 'variable-smoothk'}
            % Variable-slope smooth-k filter (Ruderman_2026, Eq. 3.6-3.7):
            %
            %   W = [1 + (k/k_M)^f(k)]^{-1}
            %
            %   f(k) = beta2 - (beta2 - beta1) * [1 + (mu * k/k_M)^delta]^{-1}
            %
            % x = k/k_M = k*R throughout
            % beta1: asymptotic slope for k/k_M << 1 (controls small-mass HMF slope)
            % beta2: asymptotic slope for k/k_M >> 1 (controls intermediate DAO regime)
            % mu   : transition scale in units of k_M (fixed: 2.1)
            % delta: sharpness of the transition     (fixed: 12)

            f = beta2 - (beta2 - beta1) ./ (1.0 + (mu .* x).^delta);
            W = 1.0 ./ (1.0 + x.^f);

        otherwise
            error(['Unknown variance filter "%s".\n' ...
                   'Valid options: tophat, gaussian, sharpk, smoothk, vsmk.'], ...
                   filter_name);
    end
end