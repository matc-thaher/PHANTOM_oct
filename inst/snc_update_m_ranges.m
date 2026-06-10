% =========================================================================
% Local helper: adjust m-ranges for retry
% =========================================================================
function [mRangeSol_new, mRangeNFW_new] = snc_update_m_ranges( ...
    mRangeSol, mRangeNFW, retryTarget, retryField, retryStep)

  mRangeSol_new = mRangeSol;
  mRangeNFW_new = mRangeNFW;
  retryTarget = lower(retryTarget);
  retryField  = lower(retryField);

  switch retryTarget
    case "soliton"
      mRangeSol_new = snc_adjust_one_range(mRangeSol_new, retryField, retryStep);
    case "nfw"
      mRangeNFW_new = snc_adjust_one_range(mRangeNFW_new, retryField, retryStep);
    case "both"
      mRangeSol_new = snc_adjust_one_range(mRangeSol_new, retryField, retryStep);
      mRangeNFW_new = snc_adjust_one_range(mRangeNFW_new, retryField, retryStep);
    otherwise
      error('soliton_nfw_composite: unknown RetryTarget ''%s''.', retryTarget);
  end
end