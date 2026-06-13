function s = sigma_R( R, z, cosmo, filter_name)
    % sigma_R  Smoothed variance sigma(R,z)
    % R can be scalar or vector (Mpc/h)
    % z is scalar
    % filter_name optional: 'tophat', 'gaussian', 'sharpk', 'smoothk', 'vsmk'

    if nargin < 4 || isempty(filter_name)
        if isfield(cosmo, 'variance_filter') && ~isempty(cosmo.variance_filter)
            filter_name = cosmo.variance_filter;
        else
            filter_name = 'tophat';
        end
    end

    s2 = sigma_R2_given_Pk(R, z, cosmo, cosmo.Pk0, filter_name);
    s  = sqrt(s2);

    end