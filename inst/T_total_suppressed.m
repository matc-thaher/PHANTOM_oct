function Ttot = T_total_suppressed(k_in, base_T_fun, supp_T_fun)
% T_total_suppressed  Combined transfer T_tot(k) = T_supp(k) * T_base(k).
%
% INPUTS:
%   k_in       : wavenumber array
%   base_T_fun : handle for base CDM transfer, T_base(k)
%   supp_T_fun : handle for suppression transfer, T_supp(k) (WDM or FDM)
%
% OUTPUT:
%   Ttot       : total transfer T_supp(k) * T_base(k)

    k_in   = k_in(:);
    T_base = base_T_fun(k_in);
    T_supp = supp_T_fun(k_in);
    Ttot   = T_supp(:) .* T_base(:);
    Ttot   = reshape(Ttot, size(k_in));
end