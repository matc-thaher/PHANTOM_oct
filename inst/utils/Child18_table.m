function p = Child18_table()
% Table 1 of Child et al. 2018, ApJ 859 55
% Parameters for Eq.(18): c = A*(M/M*)^b / [1+(M/MT)^(-b)] + c0
%                         MT = m * M*
%
% Columns: [b,      A,      m,        c0   ]
%                                        (transition  (concentration
%                                         mass ratio)  floor)
%
% Intrinsic scatter: sigma_c ~ c/3 for all fits (noted in paper text)

p.individual_all      = [-0.10,  3.44,   430.49,   3.19];
p.individual_relaxed  = [-0.09,  2.88,  1644.53,   3.54];
p.stack_nfw           = [-0.07,  4.61,   638.65,   3.59];
p.stack_einasto       = [-0.01,  63.2,   431.48,   3.36];
end