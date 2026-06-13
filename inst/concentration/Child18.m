function [c, Mstar] = Child18(M, z, cosmo, fit_type)
% Child18_concentration  Concentration-mass relation, Child et al. 2018
%
%   [c, Mstar] = Child18(M, z, cosmo)
%   [c, Mstar] = Child18(M, z, cosmo, fit_type)
%
%   Implements Eq. (18) of Child et al. 2018, ApJ 859 55:
%
%       c = A * (M/M*)^b / [1 + (M/M_T)^(-b)] + c0
%
%   where M_T = m * M* is the transition mass above which c -> c0,
%   and M*(z) is the nonlinear mass scale where sigma(M*,z) = delta_c.
%
%   Parameters from Table 1 (valid 0 <= z <= 4, M/M* fitted over 8 decades).
%
%   INPUTS
%   M        : halo mass M200c [Msun/h], scalar or vector
%   z        : redshift (scalar)
%   cosmo    : cosmology struct from cosmology.m  (needs .sigmaM)
%   fit_type : 'individual_all'     (default)
%            | 'individual_relaxed'
%            | 'stack_nfw'
%            | 'stack_einasto'
%
%   OUTPUTS
%   c        : concentration c200c, same shape as M
%   Mstar    : nonlinear mass M*(z) [Msun/h]
%
%   Simulations: Q Continuum + Outer Rim, WMAP-7 cosmology
%   Reference  : Child et al. 2018, ApJ, 859, 55
%                https://doi.org/10.3847/1538-4357/aabf95

if nargin < 4 || isempty(fit_type)
    fit_type = 'individual_all';
end

% ---- Table 1 parameters ------------------------------------------------
% Eq.(18): c = A*(M/M*)^b / [1 + (M/MT)^(-b)] + c0,  MT = m*M*
% Columns : [b,      A,      m,        c0   ]
p = Child18_table();
if ~isfield(p, fit_type)
    error('Child18_concentration: unknown fit_type "%s".\nChoose: individual_all | individual_relaxed | stack_nfw | stack_einasto', fit_type);
end
par = p.(fit_type);
m   = par(1);
A   = par(2);
b   = par(3);
c0  = par(4);

% ---- Nonlinear mass scale M*(z): sigma(M*, z) = delta_c = 1.686 --------
delta_c = 1.686;
Mstar   = Child18_Mstar(z, delta_c, cosmo);

% ---- Eq. (18) -----------------------------------------------------------
% c = A*(M/Mstar)^b / [1 + (M/Mstar)^(-b)] + c0
% Note: for M >> MT the bracket -> 1 and c -> A*(M/Mstar)^b + c0 ~ c0
%       for M << MT the bracket -> (M/MT)^b and c -> A*(MT/Mstar)^b + c0
%       giving a power-law at low mass transitioning to plateau c0 at high mass
M  = M(:);
x  = M ./ (Mstar*b);
c  = (A .* ((x.^m .* (1 + x).^(-m))-1)) + c0;
c  = reshape(c, size(M));
end
