# Package-wide numeric constants.
#
# Small shared magic numbers that were previously repeated as bare literals across several
# files. Collected here so the policy is stated once and changes in one place.

#' Default maximum profile depth (cm).
#'
#' The depth below which horizons are excluded from infilling / simulation unless a caller
#' overrides it. Every `max_depth = 250` default across the infilling and acquisition
#' functions resolves to this.
#' @noRd
DEFAULT_MAX_DEPTH_CM <- 250

#' Fallback half-width for a missing low/high bound, in representative-value units.
#'
#' When only the representative (`_r`) value of a low/representative/high triplet is present,
#' the missing `_l`/`_h` bounds are set to `_r -/+ RANGE_FALLBACK_HALFWIDTH`. Used for horizon
#' depths (`hzdept_*`/`hzdepb_*`) and component percentage (`comppct_*`). Deliberately tight:
#' a "no spread information" placeholder, not an uncertainty estimate.
#' @noRd
RANGE_FALLBACK_HALFWIDTH <- 2

#' Relative half-width of the Saxton-Rawls water-retention uncertainty band.
#'
#' `calculate_saxton_rawls_single()` returns `_l`/`_h` as the point estimate -/+ this fraction.
#' A pedotransfer-function output has no native interval; this is a nominal +/-15% band.
#' @noRd
SAXTON_RAWLS_RANGE_FRAC <- 0.15
