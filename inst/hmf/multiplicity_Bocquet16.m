function f = multiplicity_Bocquet16(sigma, z, mdef, hydro)
% Bocquet et al. (2016), MNRAS 456, 2361
% Calibrated for 200m, 200c, 500c. Separate fits for DMO and hydro.
% mdef: '200m' | '200c' | '500c'
% hydro: true (default) = hydrodynamical fit; false = DMO fit

    if nargin < 4 || isempty(hydro), hydro = false; end

    switch lower(mdef)
        case '200m'
            if hydro
                A0=0.228; a0=2.15; b0=1.69; c0=1.30;
                Az=0.285; az=-0.058; bz=-0.366; cz=-0.045;
            else
                A0=0.175; a0=1.53; b0=2.55; c0=1.19;
                Az=-0.012; az=-0.040; bz=-0.194; cz=-0.021;
            end
        case '200c'
            if hydro
                A0=0.202; a0=2.21; b0=2.00; c0=1.57;
                Az=1.147; az=0.375; bz=-1.074; cz=-0.196;
            else
                A0=0.222; a0=1.71; b0=2.24; c0=1.46;
                Az=0.269; az=0.321; bz=-0.621; cz=-0.153;
            end
        case '500c'
            if hydro
                A0=0.180; a0=2.29; b0=2.44; c0=1.97;
                Az=1.088; az=0.150; bz=-1.008; cz=-0.322;
            else
                A0=0.241; a0=2.18; b0=2.35; c0=2.02;
                Az=0.370; az=0.251; bz=-0.698; cz=-0.310;
            end
        otherwise
            error('multiplicity_Bocquet16: mdef must be ''200m'', ''200c'', or ''500c''.');
    end

    zp = 1 + z;
    A = A0 .* zp.^Az;
    a = a0 .* zp.^az;
    b = b0 .* zp.^bz;
    c = c0 .* zp.^cz;

    f = A .* ((sigma./b).^-a + 1) .* exp(-c./sigma.^2);
end