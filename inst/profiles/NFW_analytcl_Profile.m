function NFW = NFW_analytcl_Profile( Mvir, Rvir, c, r)
 % NFW_build_from_given_params
 % Construct full NFW halo profile using Mvir, Rvir, and concentration c
 %
 % INPUT:
 %   Mvir : virial mass [Msun]
 %   Rvir : virial radius [kpc]
 %   c    : concentration
 %   r    : array of radii where density is desired [kpc], e.g., linspace(0.1, Rvir, 200)
 %
 % OUTPUT (NFW struct):
 %   .rs    : scale radius [kpc]
 %   .rho_s : scale density [Msun/kpc^3]
 %   .rho   : NFW density profile at r [Msun/kpc^3]
 %   .Menc  : Mass enclosed at r [Msun]
 %   .f_c   : concentration-dependent factor
 %   .r     : radius array [kpc]
 %   .Mvir, .Rvir, .c : stored for reference

 % ---- Scale radius ----
 rs = Rvir ./ c;

 % ---- f(c) normalization (NFW1996,1997) ----
 f_c = log(1+c) - c./(1+c);

 % ---- Scale density ----
 rho_s = Mvir ./ (4*pi*rs.^3 .* f_c);

 % ---- NFW Density Profile ----
 if r(end)>Rvir
     ridx = find(r<=Rvir);
     r = r(ridx);
 end

 % rho(r) = rho_s / [ (r/rs)*(1+r/rs)^2 ]
 x = r ./ rs;
 rho = rho_s ./ ( x .* (1 + x).^2 );
 rho(r==0) = rho_s; % avoid division by zero numerically

 % ---- Enclosed Mass Profile ----
 % M(<r) = 4π ρ_s r_s^3 [ ln(1+x) - x/(1+x) ]
 Menc = 4*pi*rho_s*rs^3 .* ( log(1 + x) - x./(1 + x) );

 % ---- Output Struct ----
 NFW.Mvir  = Mvir;
 NFW.Rvir  = Rvir;
 NFW.c     = c;
 NFW.rs    = rs;
 NFW.rho_s = rho_s;
 NFW.rho   = rho;
 NFW.Menc  = Menc;
 NFW.f_c   = f_c;
 NFW.r     = r;

end