#' @title KSSL Reference Correlation Matrices (Optional Fallback)
#' @description Static, pre-computed, genetic-horizon-keyed (O/A/E/B/C/Cr,
#'   plus R for the texture matrix) correlation matrices fit once from KSSL
#'   lab data by an older, unrelated
#'   codebase (`code_ref/brdf/property_simulation.R`, outside this package).
#'   These functions let `monte-carlo.R`'s correlation-structure estimation
#'   optionally fall back to this reference data (via
#'   `config$monte_carlo$correlation_fallback = "kssl_global"`) instead of a
#'   plain identity matrix when there isn't enough SSURGO data to estimate
#'   correlations empirically. The underlying matrices are stored internally
#'   as package sysdata (see `data-raw/build_kssl_reference_correlations.R`)
#'   - not exported/documented as a public dataset, since they're
#'   implementation detail of this fallback feature, not a general-purpose
#'   soilSIM dataset.
#' @name kssl_reference_correlations
NULL

#' KSSL property correlation matrices (internal accessor)
#' @return Named list, one 9x9 correlation matrix per genhz key (`O`,`A`,`E`,`B`,`C`,`Cr`).
#' @keywords internal
.kssl_property_matrices <- function() {
  kssl_property_matrices
}

#' KSSL texture (sand/silt/clay) correlation matrices (internal accessor)
#'
#' Unused by `build_kssl_fallback_matrix()` in this version - `ilr1`/`ilr2`
#' are sourced from `.kssl_property_matrices()`'s bundled columns instead,
#' since that matrix has solid eigenvalues (0.08-0.22) across every genhz
#' group, unlike this standalone 3x3 matrix (exactly singular for `E`/`Cr`/`R`,
#' marginally positive-definite for the rest). Kept available for any future
#' use.
#' @return Named list, one 3x3 correlation matrix per genhz key (`O`,`A`,`E`,`B`,`C`,`Cr`,`R`).
#' @keywords internal
.kssl_texture_matrices <- function() {
  kssl_texture_matrices
}

#' Map soilSIM property names to KSSL reference-matrix column names
#'
#' Confirmed against the exact lookup table that built the source data
#' (`code_ref/brdf/property_simulation.R:804-816`).
#'
#' `om` -> `soc`: **caveat** - the KSSL matrix's `soc` column is actually
#' populated from `om_r` (organic matter %) in the source codebase, not true
#' lab-measured soil organic carbon. Treat it as an organic-matter proxy, not
#' true SOC.
#'
#' `ilr1`/`ilr2` map directly (same names) - valid only because
#' `monte-carlo.R`'s `composition_groups$texture$members` default is
#' `(sandtotal, silttotal, claytotal)`, matching the sequential binary
#' partition `compositions::ilr()` used to build the KSSL matrix (sand vs
#' silt+clay, then silt vs clay). See `distributions.R`'s ILR section
#' header for the positional-role convention this depends on.
#' @keywords internal
.kssl_property_name_map <- c(
  dbovendry   = "db",
  wthirdbar   = "wr_3b",
  wfifteenbar = "wr_15b",
  ph1to1h2o   = "ph",
  cec7        = "cec",
  om          = "soc",
  rfv         = "rfv",
  ilr1        = "ilr1",
  ilr2        = "ilr2"
)

#' Classify a horizon name into a generalized master-horizon designation
#'
#' Strips a leading lithologic-discontinuity digit prefix (e.g. `"2Bt2"` ->
#' `"Bt2"`), then matches against the master-horizon letters the KSSL
#' reference correlation matrices are keyed by: `O`, `A`, `E`, `B`, `C`,
#' `Cr`, `R`. `Cr` is checked before `C` so `"Cr"`/`"Crt"` aren't
#' misclassified as `"C"`.
#'
#' @param hzname Character vector of raw SSURGO horizon-designation text
#'   (e.g. `"Bt2"`, `"2Bw"`, `"Cr"`, `"R"`). `NA` is preserved as `NA`.
#' @return Character vector, same length as `hzname`, of `O`/`A`/`E`/`B`/`C`/`Cr`/`R`
#'   or `NA_character_` for unrecognized text.
#' @export
classify_genhz <- function(hzname) {
  x <- toupper(trimws(as.character(hzname)))
  na_mask <- is.na(x)
  x[na_mask] <- ""  # grepl(pattern, NA) returns NA, not FALSE - neutralize first
  x <- sub("^[0-9]+", "", x)

  patterns <- c(Cr = "^CR", R = "^R", O = "^O", A = "^A", E = "^E", B = "^B", C = "^C")
  out <- rep(NA_character_, length(x))
  for (label in names(patterns)) {
    unmatched <- is.na(out)
    out[unmatched & grepl(patterns[[label]], x)] <- label
  }
  out[na_mask] <- NA_character_
  out
}

#' Build a KSSL-derived fallback correlation matrix for a set of properties
#'
#' Starts from an identity matrix dimnamed by `properties`, then overlays the
#' *entire* principal submatrix of the matched KSSL property matrix at the
#' positions of whichever `properties` have a mapping in
#' `.kssl_property_name_map()`. Properties without a mapping stay at their
#' identity values (1 on the diagonal, 0 cross-correlation with everything,
#' including other mapped properties) - a deliberate, documented partial
#' degrade rather than an all-or-nothing failure.
#'
#' The result is positive-definite by construction: principal submatrices of
#' a positive-definite matrix are positive-definite, and a block-diagonal
#' matrix (KSSL block + identity block, zero cross-terms) with
#' positive-definite blocks is itself positive-definite. It is still passed
#' through [ensure_positive_definite_matrix()] defensively (cheap, and
#' `configure_correlation_structure()` already does this again on the final
#' matrix regardless).
#'
#' @param properties Character vector of property names to build a matrix
#'   for (in the order the returned matrix should be dimnamed).
#' @param genhz Optional single genhz value (`"O"`/`"A"`/`"E"`/`"B"`/`"C"`/`"Cr"`
#'   for the property matrix - `"R"` has no property-matrix entry). `NULL`
#'   (default) pools an unweighted average across every available genhz key
#'   (positive-definite matrices form a convex cone, so the average stays
#'   positive-definite).
#' @return A `length(properties) x length(properties)` correlation matrix
#'   dimnamed by `properties`, or `NULL` if fewer than 2 of `properties` have
#'   a KSSL mapping, or if `genhz` is supplied but not among the available
#'   keys.
#' @export
build_kssl_fallback_matrix <- function(properties, genhz = NULL) {
  kssl_matrices <- .kssl_property_matrices()

  if (!is.null(genhz)) {
    if (!genhz %in% names(kssl_matrices)) {
      return(NULL)
    }
    kssl_source <- kssl_matrices[[genhz]]
  } else {
    kssl_source <- Reduce(`+`, kssl_matrices) / length(kssl_matrices)
  }

  mapped_kssl_names <- .kssl_property_name_map[properties]
  names(mapped_kssl_names) <- properties
  mapped_kssl_names <- mapped_kssl_names[!is.na(mapped_kssl_names) & mapped_kssl_names %in% rownames(kssl_source)]

  if (length(mapped_kssl_names) < 2) {
    return(NULL)
  }

  result <- diag(length(properties))
  rownames(result) <- colnames(result) <- properties

  mapped_properties <- names(mapped_kssl_names)
  result[mapped_properties, mapped_properties] <- kssl_source[mapped_kssl_names, mapped_kssl_names]

  ensure_positive_definite_matrix(result)
}
