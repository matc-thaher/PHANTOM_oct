function xi = correlation_function(R, z, cosmo, method, Pk_handle)
% CORRELATION_FUNCTION  Linear matter correlation function xi_mm(R,z).
%
%   xi = CORRELATION_FUNCTION(R, z, cosmo)
%   xi = CORRELATION_FUNCTION(R, z, cosmo, method)
%   xi = CORRELATION_FUNCTION(R, z, cosmo, method, Pk_handle)
%
%   Computes the linear matter correlation function
%
%       xi(R,z) = (1 / 2*pi^2) * \int_0^\infty dk k^2 P(k,z) j0(kR),
%
%   where j0(x) = sin(x)/x.
%
%   Inputs
%   ------
%   R         : Separation in comoving Mpc/h. Scalar or vector.
%   z         : Redshift.
%   cosmo     : Cosmology struct.
%   method    : Numerical method:
%               'integral'  - direct integral (default)
%               'fftlog'    - FFTLog discrete Hankel transform
%   Pk_handle : Optional handle for P(k) at z=0. Default: cosmo.Pk0
%
%   Output
%   ------
%   xi        : Linear matter correlation function.
%
%   Notes
%   -----
%   - The direct-integral method is the recommended default because it is
%     most transparent.
%   - The FFTLog method follows the discrete logarithmic Hankel transform
%     algorithm of Hamilton (2000), which is well suited to smooth spectra
%     spanning many decades in k.
%
%   References
%   ----------
%   Diemer 2018, COLOSSUS, Section 2.7.
%   Hamilton 2000, FFTLog.
%
%   Example
%   -------
%   R  = logspace(-2,2,100);
%   xi = correlation_function(R, 0, cosmo);           % default
%   xf = correlation_function(R, 0, cosmo, 'fftlog');

    if nargin < 4 || isempty(method)
        method = 'integral';
    end
    if nargin < 5 || isempty(Pk_handle)
        Pk_handle = cosmo.Pk0;
    end

    switch lower(method)
        case {'integral','direct'}
            xi = correlation_function_integral(R, z, cosmo, Pk_handle);

        case {'fftlog'}
            xi = correlation_function_fftlog(R, z, cosmo, Pk_handle);

        otherwise
            error('Unknown method "%s". Use "integral" or "fftlog".', method);
    end
end