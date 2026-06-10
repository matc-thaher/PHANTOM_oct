function T = T_Sugiyama95(k, cosmo)
% Transfer function: Sugiyama (1995)
%
% INPUT:
%   k     : wavenumber in h/Mpc
%   cosmo : cosmology struct
%
% REQUIRED FIELDS:
%   cosmo.Omega_m
%   cosmo.Omega_b
%   cosmo.h
%
% OPTIONAL:
%   cosmo.Tcmb0   (default = 2.7255 K)

    Om = cosmo.Omega_m;
    Ob = cosmo.Omega_b;
    h  = cosmo.h;

    if isfield(cosmo, 'Tcmb0') && ~isempty(cosmo.Tcmb0)
        Tcmb = cosmo.Tcmb0;
    else
        Tcmb = 2.7255;
    end

    % Effective shape parameter with baryon correction
    Gamma = Om * h * exp(-Ob * (1 + sqrt(2*h) / Om));

    % BBKS/Sugiyama dimensionless scale
    q = (Tcmb / 2.7)^2 .* k ./ Gamma;

    % Transfer function
    L0 = log(1 + 2.34 .* q) ./ (2.34 .* q);
    C0 = (1 + 3.89 .* q + (16.1 .* q).^2 + (5.46 .* q).^3 + (6.71 .* q).^4).^(-0.25);

    T = L0 .* C0;

    % Small-q limit
    T(q < 1e-9) = 1.0;

end