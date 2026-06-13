function range_new = snc_adjust_one_range(range_old, retryField, retryStep)
  range_new = range_old;

  switch retryField
    case "min"
      range_new(1) = range_new(1) + retryStep;
    case "max"
      range_new(2) = range_new(2) + retryStep;
    case "both"
      range_new(1) = range_new(1) + retryStep;
      range_new(2) = range_new(2) + retryStep;
    otherwise
      error('soliton_nfw_composite: unknown RetryField ''%s''.', retryField);
  end

  if range_new(1) >= range_new(2)
    error('soliton_nfw_composite: invalid range after retry update: [%.3f, %.3f].', ...
          range_new(1), range_new(2));
  end
end