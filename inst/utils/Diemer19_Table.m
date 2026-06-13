function P = Diemer19_Table(mode)
% Diemer19_Table  Best-fit parameters for the Diemer & Joyce (2019) model
%
%   P = Diemer19_Table(mode)
%
%   Returns parameters from Table 2 of Diemer & Joyce 2019, ApJ 871, 168.
%
%   MODE options (case-insensitive)
%   --------------------------------
%   'median'  :  median concentration statistic
%   'mean'    :  mean   concentration statistic
%
%   Fields returned: P.kappa, P.a0, P.a1, P.b0, P.b1, P.cAlpha
%
%   Reference: Diemer & Joyce 2019, ApJ 871, 168, Table 2
%              https://doi.org/10.3847/1538-4357/ab4947

switch lower(mode)

    case 'mean'
        P.kappa  = 0.42;
        P.a0     = 2.37;
        P.a1     = 1.74;
        P.b0     = 3.39;
        P.b1     = 1.82;
        P.cAlpha = 0.20;

    case 'median'
        P.kappa  = 0.41;
        P.a0     = 2.45;
        P.a1     = 1.82;
        P.b0     = 3.20;
        P.b1     = 2.30;
        P.cAlpha = 0.21;

    otherwise
        error('Diemer19_Table: unknown mode "%s". Valid: median | mean.', mode);
end
end