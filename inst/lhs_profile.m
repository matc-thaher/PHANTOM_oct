function y = lhs_profile(c, expo, profile_name)
    c = max(c, realmin);
    mu = profile_mu(c, profile_name);
    y  = log10(c ./ (mu.^expo));
end