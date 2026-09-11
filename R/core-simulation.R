#' @title Component-Composition and Correlated-Triangular Property Simulation
#'
#' @description Component-composition simulation (`simulate_component_composition()`), correlated
#'   triangular-distribution sampling (`simulate_correlated_triangular()`), and per-cokey
#'   flexible-property simulation (`simulate_cokey_generalized()`), plus two small standalone
#'   horizon-data helpers (`remove_organic_layer()`, `slice_and_aggregate_soil_data()`).
#'
#'   `simulate_correlated_triangular()` draws uncorrelated normals with `stats::rnorm()` (no
#'   `MASS` dependency), and `simulate_cokey_generalized()`'s texture step uses
#'   `R/core-distributions.R`'s `ilr_forward()`/`ilr_inverse()` (no `compositions` dependency).
#' @name property_simulation
NULL

#' Remove Organic Layers and Adjust Depths
#'
#' Removes rows with organic horizons (where `hzname` contains a capital "O") and recalculates
#' the remaining horizons' depths within each `cokey` group so they stay cumulative and
#' re-anchored to 0, preserving each horizon's thickness.
#'
#' @param df A data frame containing soil horizon data, with columns `cokey`, `hzname`,
#'   `hzdept_r`, `hzdepb_r`, and optionally `hzdept_l`, `hzdepb_l`, `hzdept_h`, `hzdepb_h`.
#' @return A data frame with organic layers removed and depths adjusted within each `cokey` group.
#' @keywords internal
remove_organic_layer <- function(df) {

  process_group <- function(group) {
    filtered_group <- group[!grepl("O", group$hzname), ]

    if (nrow(filtered_group) == 0) {
      return(filtered_group)
    }

    initial_hzdept_r <- filtered_group$hzdept_r[1]

    filtered_group <- filtered_group |>
      dplyr::mutate(dplyr::across(dplyr::starts_with("hzdep"), ~ . - initial_hzdept_r))

    for (i in seq_len(nrow(filtered_group))) {
      if (filtered_group$hzdept_r[i] == 0) {
        if ("hzdept_l" %in% names(filtered_group)) {
          if (!is.na(filtered_group$hzdept_l[i])) {
            filtered_group$hzdept_l[i] <- 0
          }
        }
        if ("hzdept_h" %in% names(filtered_group)) {
          if (!is.na(filtered_group$hzdept_h[i])) {
            filtered_group$hzdept_h[i] <- 0
          }
        }
      }
    }

    thickness_r <- filtered_group$hzdepb_r - filtered_group$hzdept_r
    filtered_group$hzdept_r <- cumsum(c(0, thickness_r[-length(thickness_r)]))
    filtered_group$hzdepb_r <- cumsum(thickness_r)

    if (all(c("hzdept_l", "hzdepb_l") %in% names(filtered_group))) {
      if (!all(is.na(filtered_group$hzdept_l))) {
        thickness_l <- filtered_group$hzdepb_l - filtered_group$hzdept_l
        filtered_group$hzdept_l <- cumsum(c(0, thickness_l[-length(thickness_l)]))
        filtered_group$hzdepb_l <- cumsum(thickness_l)
      }
    }

    if (all(c("hzdept_h", "hzdepb_h") %in% names(filtered_group))) {
      if (!all(is.na(filtered_group$hzdept_h))) {
        thickness_h <- filtered_group$hzdepb_h - filtered_group$hzdept_h
        filtered_group$hzdept_h <- cumsum(c(0, thickness_h[-length(thickness_h)]))
        filtered_group$hzdepb_h <- cumsum(thickness_h)
      }
    }

    return(filtered_group)
  }

  df |>
    dplyr::group_by(cokey) |>
    dplyr::group_modify(~ process_group(.x)) |>
    dplyr::ungroup()
}

#' Slice and Aggregate Soil Data at Specified Depth Intervals
#'
#' Slices a horizon data frame into 1 cm increments based on `hzdept_r`/`hzdepb_r`, then
#' averages every numeric column within each requested depth range.
#'
#' @param df A data frame with `hzdept_r`/`hzdepb_r` depth-range columns.
#' @param depth_ranges A list of length-2 `c(top, bottom)` vectors.
#' @return A data frame with one row per depth range, mean values of soil properties for each.
#' @keywords internal
slice_and_aggregate_soil_data <- function(df, depth_ranges = list(c(0, 30), c(30, 100))) {

  if (!all(c("hzdept_r", "hzdepb_r") %in% colnames(df))) {
    stop("Input data frame must contain 'hzdept_r' and 'hzdepb_r' columns")
  }

  data_columns <- df |>
    dplyr::select(dplyr::where(is.numeric)) |>
    dplyr::select(-hzdept_r, -hzdepb_r) |>
    colnames()

  # Vectorized: row_idx repeats each source row's index
  # once per depth it covers, and depths_vec computes each row's 1cm depth sequence via
  # sequence()'s standard "concatenated per-group seq_len()" trick, instead of building and
  # rbind()-ing one tiny data.frame per depth.
  n_expand <- pmax(0, df$hzdepb_r - df$hzdept_r)
  row_idx <- rep(seq_len(nrow(df)), times = n_expand)
  depths_vec <- df$hzdept_r[row_idx] + sequence(n_expand) - 1

  aggregated_data <- df[row_idx, data_columns, drop = FALSE]
  aggregated_data$Depth <- as.numeric(depths_vec)
  rownames(aggregated_data) <- NULL

  results <- list()
  for (i in seq_along(depth_ranges)) {
    range <- depth_ranges[[i]]
    top <- range[1]
    bottom <- range[2]

    subset <- aggregated_data |>
      dplyr::select(dplyr::where(is.numeric)) |>
      dplyr::filter(Depth >= top & Depth < bottom)

    if (nrow(subset) > 0) {
      actual_max_depth <- max(subset$Depth)
      mean_values <- colMeans(subset, na.rm = TRUE)
      mean_values["hzdept_r"] <- top

      if (actual_max_depth + 1 < bottom) {
        mean_values["hzdepb_r"] <- actual_max_depth + 1
      } else {
        mean_values["hzdepb_r"] <- bottom
      }

      if ("compname" %in% colnames(aggregated_data)) {
        result_row <- c(compname = aggregated_data$compname[1], mean_values)
      } else {
        result_row <- mean_values
      }

      results <- append(results, list(result_row))
    } else {
      if ("compname" %in% colnames(aggregated_data)) {
        na_row <- c(compname = aggregated_data$compname[1],
                    rep(NA, length(data_columns) + 1),
                    hzdept_r = top,
                    hzdepb_r = bottom)
        names(na_row) <- c("compname", data_columns, "Depth", "hzdept_r", "hzdepb_r")
      } else {
        na_row <- c(rep(NA, length(data_columns) + 1),
                    hzdept_r = top,
                    hzdepb_r = bottom)
        names(na_row) <- c(data_columns, "Depth", "hzdept_r", "hzdepb_r")
      }

      results <- append(results, list(na_row))
    }
  }

  if (length(results) > 0) {
    result_df <- do.call(rbind, results) |> as.data.frame()

    numeric_cols <- c(data_columns, "Depth", "hzdept_r", "hzdepb_r")
    for (col in numeric_cols) {
      if (col %in% colnames(result_df)) {
        result_df[[col]] <- as.numeric(result_df[[col]])
      }
    }
  } else {
    result_df <- data.frame()
  }

  return(result_df)
}

#' Simulate Soil Component Composition
#'
#' Simulates component-composition percentages from low/rep/high (`comppct_l/r/h`) values via
#' `tri_dist()`, one component (`cokey`) at a time, and derives a single `sim_comppct` value per
#' component.
#'
#' @section Grain mismatch with `R/core-simulation.R`:
#' This returns one row per **component** (`cokey`), but `simulate_and_perturb_soil_profiles()`/
#' `simulate_profile_depths_by_mukey()` need `sim_comppct` on every **horizon** row. Callers must
#' `dplyr::left_join()` this output onto horizon-level data by `cokey` before passing it to those
#' functions - this function alone does not satisfy their requirement.
#'
#' @section `sim_comppct`'s derivation:
#' `sim_comppct <- round(sum(<n_simulations> triangular draws of comppct) / 100)`, which is
#' `round(n_simulations * comppct_r / 100)` in expectation (e.g. `comppct_r = 30`,
#' `n_simulations = 1000` -> `sim_comppct` ~= 300). It is not an independently-meaningful
#' simulation count.
#'
#' @param data Data frame with `mukey`, `cokey`, `compname`, `comppct_l`, `comppct_r`, `comppct_h`
#'   (one or more rows per component; deduplicated internally).
#' @param n_simulations Integer, number of triangular draws per component (default 1000).
#' @return A data frame, one row per component, with an added `sim_comppct` column.
#' @export
simulate_component_composition <- function(data, n_simulations = 1000) {

  data <- data |>
    dplyr::select(mukey, cokey, compname, comppct_l, comppct_r, comppct_h) |>
    dplyr::distinct()

  data$comppct_l[is.na(data$comppct_l)] <- data$comppct_r[is.na(data$comppct_l)] - RANGE_FALLBACK_HALFWIDTH
  data$comppct_h[is.na(data$comppct_h)] <- data$comppct_r[is.na(data$comppct_h)] + RANGE_FALLBACK_HALFWIDTH

  n_components <- nrow(data)
  compositions <- matrix(NA, nrow = n_simulations, ncol = n_components)

  for (i in 1:n_components) {
    low  <- data$comppct_l[i]
    mode <- data$comppct_r[i]
    high <- data$comppct_h[i]

    compositions[, i] <- tri_dist(n_simulations, a = low, b = high, c = mode)
  }

  sim_compositions <- apply(compositions, 2, function(col) {
    round(sum(col) / 100, 0)
  })

  data$sim_comppct <- sim_compositions

  return(data)
}

#' Simulate Correlated Samples from Triangular Distributions
#'
#' Generates `n` correlated samples from `length(params)` triangular distributions via a
#' Cholesky decomposition of `correlation_matrix` applied to standard-normal draws, then
#' transforms each column to its own triangular distribution via the inverse CDF.
#'
#' @param n Integer, number of samples to generate.
#' @param params List of `c(a, b, c)` triples, one per distribution - **note the order here is
#'   (lower, mode, upper)**, unlike `tri_dist()`'s own `(a = lower, b = upper, c = mode)`
#'   convention. A real source of confusion if assumed to match `tri_dist()`'s argument order.
#' @param correlation_matrix A square, positive-semi-definite correlation matrix,
#'   `length(params)` x `length(params)`.
#' @param random_seed Optional integer seed for reproducibility.
#' @return A matrix of correlated samples, `n` rows x `length(params)` columns.
#' @keywords internal
simulate_correlated_triangular <- function(n, params, correlation_matrix, random_seed = NULL) {

  if (!is.null(random_seed)) {
    set.seed(random_seed)
  }

  k <- length(params)
  # Equivalent to MASS::mvrnorm(n, mu = rep(0, k), Sigma = diag(k)) - a multivariate normal with
  # an identity covariance is just k independent standard normal draws, so no MASS dependency is
  # needed.
  uncorrelated_normal <- matrix(stats::rnorm(n * k), nrow = n, ncol = k)

  L <- chol(correlation_matrix)
  correlated_normal <- uncorrelated_normal %*% L

  samples <- matrix(NA, nrow = n, ncol = k)

  for (i in seq_along(params)) {
    a <- params[[i]][1]
    b <- params[[i]][2]
    c <- params[[i]][3]

    normal_var <- correlated_normal[, i]
    u <- stats::pnorm(normal_var)

    if (c == a) {
      samples[, i] <- rep(a, n)
    } else {
      condition <- u <= (b - a) / (c - a)
      samples[condition, i] <- a + sqrt(u[condition] * (c - a) * (b - a))
      samples[!condition, i] <- c - sqrt((1 - u[!condition]) * (c - a) * (c - b))
    }
  }

  return(samples)
}

#' Compute the Mode of a Numeric Vector
#'
#' @param x A vector of values.
#' @return The most frequently-occurring value in `x`.
#' @section Behavior:
#' `tabulate(match(x, ux))` on pre-sorted unique values counts occurrences without the
#' factor-coercion overhead of `table(x)`/`factor(x)`. On ties it returns the smallest value
#' among them (unique values are sorted ascending and `which.max()` returns the first maximum).
#' @keywords internal
compute_mode <- function(x) {
  ux <- sort(unique(x))
  counts <- tabulate(match(x, ux))
  ux[which.max(counts)]
}

#' Van Bemmelen-factor conversion from SSURGO organic matter % to an estimated soil organic
#' carbon %
#'
#' SSURGO's `chorizon` table reports organic matter (`om_l/_r/_h`), not organic carbon - but
#' `simulate_cokey_generalized()`'s simulated `soc` column (and the KSSL reference correlation
#' matrix's `soc` row/column, see `R/core-correlations.R`'s `.kssl_property_name_map`
#' docs) both need an actual SOC-scale value to be fused against SOLUS100's real `soc` variable
#' (organic carbon). `SOC = OM / 1.724` is the standard Van Bemmelen-factor approximation (the
#' same factor `remarginalized_awc()`'s `soc_to_om` argument uses in the opposite direction,
#' `R/core-fusion.R`). Applying this conversion needs no change to the KSSL correlation
#' *matrix* itself - Pearson correlation is invariant under a positive linear rescale of one
#' variable, and triangular-distribution parameters (min/mode/max) transform linearly too - only
#' the marginal values entering the correlated draw need rescaling, which is exactly where this
#' constant is used.
#' @keywords internal
OM_TO_SOC_FACTOR <- 1 / 1.724

#' Extend a Correlation Matrix with Identity Rows/Columns for Missing Names
#'
#' Adds a 1-on-diagonal, 0-cross-correlation row/column for each name in `missing_names` not
#' already in `m`'s dimnames - the same "uncorrelated with everything, safe default" convention
#' `build_kssl_fallback_matrix()` uses for properties without KSSL reference data. Used by
#' `simulate_cokey_generalized()` as a defensive
#' fallback so an unrecognized-but-requested property degrades to independence rather than
#' erroring on `m[keep_cols, keep_cols]`. The result stays positive-definite: a block-diagonal
#' matrix formed from a positive-definite block (`m`) and an identity block has no negative
#' eigenvalues.
#' @param m A square, dimnamed correlation matrix.
#' @param missing_names Character vector of names to add (assumed disjoint from `rownames(m)`).
#' @return `m` extended to `nrow(m) + length(missing_names)` rows/columns.
#' @keywords internal
extend_corr_matrix_with_identity <- function(m, missing_names) {
  old_names <- rownames(m)
  new_names <- c(old_names, missing_names)
  n_old <- length(old_names)
  n_new <- length(new_names)
  extended <- diag(n_new)
  dimnames(extended) <- list(new_names, new_names)
  extended[seq_len(n_old), seq_len(n_old)] <- m
  extended
}

#' Simulate Soil Properties for a Specific Cokey (Generalized)
#'
#' Simulates a flexible set of soil properties (whichever of bulk density, water retention,
#' texture, rock fragment volume, pH, CEC, organic carbon have complete `_l/_r/_h` triplets
#' present) for each row of a single cokey's horizon data, via
#' `simulate_correlated_triangular()`, using a genhz-keyed local correlation matrix. If texture
#' columns are present, sand/silt/clay are simulated jointly and converted to ILR coordinates
#' (`ilr1`/`ilr2`) before being included in the main correlated simulation, then converted back
#' after.
#'
#' @param sim_cokey A data frame of horizon rows for one cokey, with a `genhz` column and a
#'   `sim_comppct` column (number of realizations to simulate for that row - see
#'   `simulate_component_composition()`).
#' @param correlation_matrices A list of correlation matrices keyed by `genhz`, with row/column
#'   names matching (a subset of) `c("db", "wr_3b", "wr_15b", "ilr1", "ilr2", "rfv", "ph", "cec",
#'   "soc", "caco3", "ec", "ecec", "gypsum", "sar")`. The 5 chemistry properties
#'   `caco3`/`ec`/`ecec`/`gypsum`/`sar` have no real KSSL-fit correlation data (see
#'   `build_kssl_fallback_matrix()`) - `simulate_ssurgo_mapunit_draws()` builds its
#'   `correlation_matrices` through that function so they're present (as identity/uncorrelated)
#'   rather than missing; a direct caller supplying its own matrix without them will simply never
#'   populate `param_list[["caco3"]]` etc. (see `get_param_set()`'s per-row gating below), not error.
#' @param txt_correlation_matrices A list of texture correlation matrices keyed by `genhz`
#'   (sand/silt/clay order, matching the `simulate_correlated_triangular()` call for texture).
#' @param requested_properties `NULL` (default) simulates every property with a complete `_l/_r/_h`
#'   triplet, unchanged. Otherwise a character vector of property identifiers - in any vocabulary
#'   `normalize_requested_properties()` accepts (caller ids, SSURGO stems, output column names, or
#'   the internal `param_order` names) - restricting the simulation to just those (plus the texture
#'   coupling: any texture member pulls in the full sand/silt/clay draw and both ILR axes). The
#'   Gaussian-copula marginals of a restricted simulation match the corresponding columns of a full
#'   one under the default `joint_copula` vertical-correlation method (modulo RNG-stream Monte Carlo
#'   noise); see `simulate_ssurgo_mapunit_draws()` for the `gp_quantile_retrofit` caveat.
#' @return A data frame of simulated property values across all rows/realizations, with
#'   `compname`, `mukey`, `cokey`, `hzdept_r`, `hzdepb_r`, `simulation_number`, `unique_id`, and
#'   (when `sim_cokey` itself has a `bound_sd` column - see `attach_osd_boundary_distinctness()`)
#'   `bound_sd`. The `soc` column is an SSURGO-organic-matter-**derived SOC estimate**
#'   (`om * OM_TO_SOC_FACTOR`, the inverse Van Bemmelen factor) - not a lab-measured soil organic
#'   carbon value. See `OM_TO_SOC_FACTOR`'s own docs.
#' @keywords internal
simulate_cokey_generalized <- function(sim_cokey, correlation_matrices, txt_correlation_matrices = NULL,
                                        requested_properties = NULL) {

  # The 5 chemistry properties (caco3/ec/ecec/gypsum/sar) have no real KSSL-fit correlation
  # entry (see correlation_matrices' own @param doc) -
  # they simulate uncorrelated with everything else by default until real KSSL lab data covering
  # them becomes available for a re-fit.
  param_order <- c("db", "wr_3b", "wr_15b", "ilr1", "ilr2", "rfv", "ph", "cec", "soc",
                   "caco3", "ec", "ecec", "gypsum", "sar")

  # NULL -> keep everything (unchanged). Otherwise restrict to the requested subset; normalize here
  # (idempotently) so a direct caller can pass ids like "clay"/"ph" rather than param_order names.
  keep_params <- if (is.null(requested_properties)) {
    param_order
  } else {
    intersect(param_order, normalize_requested_properties(requested_properties))
  }
  want_texture <- any(c("ilr1", "ilr2") %in% keep_params)

  texture_cols_required <- c("sandtotal_l", "sandtotal_r", "sandtotal_h",
                             "silttotal_l", "silttotal_r", "silttotal_h",
                             "claytotal_l", "claytotal_r", "claytotal_h")

  get_param_set <- function(row, prefix) {
    lcol <- paste0(prefix, "_l")
    rcol <- paste0(prefix, "_r")
    hcol <- paste0(prefix, "_h")
    if (all(c(lcol, rcol, hcol) %in% names(row))) {
      return(c(row[[lcol]], row[[rcol]], row[[hcol]]))
    }
    return(NULL)
  }

  sim_data_out <- list()

  # Pooled, genhz-agnostic fallback matrices - an unweighted average of every available
  # genhz key's matrix (positive-definite matrices form a convex cone, so the average
  # stays positive-definite; the same reasoning build_kssl_fallback_matrix() documents for
  # its own genhz = NULL pooling path). Used below whenever a row's genhz can't be matched
  # to a specific key, instead of crashing.
  pooled_property_corr <- Reduce(`+`, correlation_matrices) / length(correlation_matrices)
  pooled_txt_corr <- if (!is.null(txt_correlation_matrices)) {
    Reduce(`+`, txt_correlation_matrices) / length(txt_correlation_matrices)
  } else {
    NULL
  }

  # ensure_positive_definite_matrix(txt_corr) is a deterministic function of genhz_val alone
  # (txt_corr is looked up per genhz_val, with the same pooled fallback when missing), so
  # cache it per genhz_val rather than recomputing eigen() on every row.
  pd_txt_corr_cache <- list()

  for (i in seq_len(nrow(sim_cokey))) {
    row <- sim_cokey[i, ]
    genhz_val <- as.character(row$genhz)
    # classify_genhz() returns NA for a raw hzname that doesn't parse into a KSSL-covered
    # master-horizon category (e.g. the generic "H1"/"H2"/"H3" designations some SSURGO
    # components use). `as.character(NA)` is NA_character_, and an NA_character_ key is
    # write-only in a base list: `l[[NA_character_]] <- x` creates an element that
    # `l[[NA_character_]]` then reads back as NULL. That silently defeated the per-genhz
    # `pd_txt_corr_cache` below - `txt_corr` came back NULL, `simulate_correlated_triangular()`
    # got a NULL correlation matrix, and the texture draw failed for EVERY such row (contained
    # by the tryCatch, but the row then contributed no texture at all). Normalize to a real
    # string key so all the `[[genhz_val]]` lookups miss cleanly and degrade to the pooled,
    # genhz-agnostic matrices, which is the intended graceful-degrade behavior.
    if (is.na(genhz_val) || !nzchar(genhz_val)) genhz_val <- "__unclassified__"

    local_corr <- correlation_matrices[[genhz_val]]
    if (is.null(local_corr)) {
      # classify_genhz() returns NA for a raw hzname that doesn't parse into one of the 7
      # KSSL-covered master-horizon categories (blank/garbled/nonstandard hzname text).
      # Degrade to the pooled, genhz-agnostic matrix rather than aborting the whole cokey.
      local_corr <- pooled_property_corr
    }

    has_texture <- all(texture_cols_required %in% names(row))

    ilr1_lrh <- NULL
    ilr2_lrh <- NULL

    if (want_texture && has_texture && !is.null(txt_correlation_matrices)) {
      # .kssl_texture_matrices() is documented as exactly singular for genhz E/Cr/R even
      # when correctly matched (see R/core-correlations.R) - chol() inside
      # simulate_correlated_triangular() would otherwise error deterministically for those
      # groups. Nudge
      # to the nearest positive-definite matrix defensively; cheap, and the same pattern
      # already used elsewhere in this package (e.g. build_kssl_fallback_matrix()). Cached per
      # genhz_val (see pd_txt_corr_cache's own comment above) since it's deterministic given
      # genhz_val alone - avoids recomputing an identical eigen() decomposition on every row.
      if (is.null(pd_txt_corr_cache[[genhz_val]])) {
        txt_corr_raw <- txt_correlation_matrices[[genhz_val]]
        if (is.null(txt_corr_raw)) {
          txt_corr_raw <- pooled_txt_corr
        }
        pd_txt_corr_cache[[genhz_val]] <- ensure_positive_definite_matrix(txt_corr_raw)
      }
      txt_corr <- pd_txt_corr_cache[[genhz_val]]
      # Order here (sand, silt, clay) must match txt_corr's own row/column order.
      params_txt <- list(
        c(row$sandtotal_l, row$sandtotal_r, row$sandtotal_h),
        c(row$silttotal_l, row$silttotal_r, row$silttotal_h),
        c(row$claytotal_l, row$claytotal_r, row$claytotal_h)
      )

      # A texture triplet with a leftover NA (e.g. infilling couldn't recover a value for
      # this specific row) makes simulate_correlated_triangular()'s internal `if (c == a)`
      # check evaluate to NA ("missing value where TRUE/FALSE needed"). Contain the failure to
      # "no texture for this row" instead: ilr1_lrh/ilr2_lrh
      # stay NULL (as initialized above), and the has_texture check further below already
      # treats NULL ilr1_lrh/ilr2_lrh as "texture unavailable for this row" gracefully.
      texture_result <- tryCatch({
        sim_txt <- simulate_correlated_triangular(as.integer(row$sim_comppct), params_txt, txt_corr)

        # ilr_forward()'s validated-to-match-compositions::ilr() transform, fed the physical
        # sand/silt/clay values by role (not by column position) - see file header.
        sim_txt_ilr <- ilr_forward(clay = sim_txt[, 3], sand = sim_txt[, 1], silt = sim_txt[, 2])
        ilr1_vals <- sim_txt_ilr[, "z1"]
        ilr2_vals <- sim_txt_ilr[, "z2"]
        list(
          ilr1_lrh = c(min(ilr1_vals), compute_mode(ilr1_vals), max(ilr1_vals)),
          ilr2_lrh = c(min(ilr2_vals), compute_mode(ilr2_vals), max(ilr2_vals))
        )
      }, error = function(e) {
        message("simulate_cokey_generalized(): texture simulation failed for row ", i, ": ", e$message)
        NULL
      })

      if (!is.null(texture_result)) {
        ilr1_lrh <- texture_result$ilr1_lrh
        ilr2_lrh <- texture_result$ilr2_lrh
      }
    }

    param_list <- vector("list", length(param_order))
    names(param_list) <- param_order

    db_set <- get_param_set(row, "dbovendry")
    if ("db" %in% keep_params && !is.null(db_set)) param_list[["db"]] <- db_set

    wr3_set <- get_param_set(row, "wthirdbar")
    if ("wr_3b" %in% keep_params && !is.null(wr3_set)) param_list[["wr_3b"]] <- wr3_set

    wr15_set <- get_param_set(row, "wfifteenbar")
    if ("wr_15b" %in% keep_params && !is.null(wr15_set)) param_list[["wr_15b"]] <- wr15_set

    if (want_texture && has_texture && !is.null(ilr1_lrh) && !is.null(ilr2_lrh)) {
      param_list[["ilr1"]] <- ilr1_lrh
      param_list[["ilr2"]] <- ilr2_lrh
    }

    rfv_set <- get_param_set(row, "rfv")
    if ("rfv" %in% keep_params && !is.null(rfv_set)) param_list[["rfv"]] <- rfv_set

    ph_set <- get_param_set(row, "ph1to1h2o")
    if ("ph" %in% keep_params && !is.null(ph_set)) param_list[["ph"]] <- ph_set

    cec_set <- get_param_set(row, "cec7")
    if ("cec" %in% keep_params && !is.null(cec_set)) param_list[["cec"]] <- cec_set

    # om -> soc conversion: SSURGO reports organic
    # matter, not organic carbon - convert via the Van Bemmelen factor so this triplet is
    # SOC-scale before it enters the correlated draw. See OM_TO_SOC_FACTOR's own docs for why
    # the KSSL correlation matrix itself needs no corresponding change (scale invariance).
    om_set <- get_param_set(row, "om")
    if ("soc" %in% keep_params && !is.null(om_set)) param_list[["soc"]] <- om_set * OM_TO_SOC_FACTOR

    # 5 chemistry properties - same SSURGO-stem-triplet
    # gating pattern as every property above, no conversion needed (SSURGO's units already match
    # SOLUS100's: caco3/gypsum are % by weight, ec is dS/m, ecec is cmol(+)/kg, sar is unitless).
    caco3_set <- get_param_set(row, "caco3")
    if ("caco3" %in% keep_params && !is.null(caco3_set)) param_list[["caco3"]] <- caco3_set

    ec_set <- get_param_set(row, "ec")
    if ("ec" %in% keep_params && !is.null(ec_set)) param_list[["ec"]] <- ec_set

    ecec_set <- get_param_set(row, "ecec")
    if ("ecec" %in% keep_params && !is.null(ecec_set)) param_list[["ecec"]] <- ecec_set

    gypsum_set <- get_param_set(row, "gypsum")
    if ("gypsum" %in% keep_params && !is.null(gypsum_set)) param_list[["gypsum"]] <- gypsum_set

    sar_set <- get_param_set(row, "sar")
    if ("sar" %in% keep_params && !is.null(sar_set)) param_list[["sar"]] <- sar_set

    param_list <- param_list[!vapply(param_list, is.null, logical(1))]

    if (!length(param_list)) {
      message("No recognized properties for row ", i, ". Skipping.")
      next
    }

    params_for_sim <- unname(param_list)

    keep_cols <- names(param_list)

    # Defensive: a keep_cols name absent from local_corr's own dimnames (e.g. a direct caller
    # supplying a correlation_matrices list that predates task P2's 5 new chemistry properties,
    # or any future property added to param_order before every caller's matrices catch up) would
    # otherwise make `local_corr[keep_cols, keep_cols]` below throw "subscript out of bounds" and
    # abort this row (this indexing happens before the tryCatch() further down). Extend local_corr
    # with an identity row/col for any such name instead - uncorrelated with everything, same
    # graceful-degrade convention build_kssl_fallback_matrix() and this function's own
    # genhz-mismatch/singular-texture-matrix fallbacks already use elsewhere in this file. The
    # production path (simulate_ssurgo_mapunit_draws()) never hits this: it builds
    # correlation_matrices via build_kssl_fallback_matrix() sized to the full param_order already.
    missing_cols <- setdiff(keep_cols, rownames(local_corr))
    if (length(missing_cols)) {
      local_corr <- extend_corr_matrix_with_identity(local_corr, missing_cols)
    }
    local_corr_sub <- local_corr[keep_cols, keep_cols, drop = FALSE]

    tryCatch({
      n_sim <- as.integer(row$sim_comppct)
      sim_data <- simulate_correlated_triangular(
        n = n_sim,
        params = params_for_sim,
        correlation_matrix = local_corr_sub
      )
      sim_data <- as.data.frame(sim_data)
      colnames(sim_data) <- names(param_list)

      if ("wr_3b" %in% names(param_list)) {
        sim_data[["wr_3b"]] <- sim_data[["wr_3b"]] / 100
      }
      if ("wr_15b" %in% names(param_list)) {
        sim_data[["wr_15b"]] <- sim_data[["wr_15b"]] / 100
      }

      if (has_texture && all(c("ilr1", "ilr2") %in% names(param_list))) {
        # ilr_inverse()'s `total = 100` default already returns percentage-scale values summing
        # to 100 - unlike compositions::ilrInv() (fraction-scale, summing to 1), no separate
        # "* 100" rescale step is needed here.
        sim_txt <- ilr_inverse(sim_data[["ilr1"]], sim_data[["ilr2"]], total = 100)
        sim_txt <- as.data.frame(sim_txt)[, c("sand", "silt", "clay")]
        colnames(sim_txt) <- c("sand_total", "silt_total", "clay_total")

        # drop = FALSE: with requested_properties restricting the simulation, exactly one
        # non-texture column can remain here (e.g. requested_properties = c("clay", "ph")), and
        # `df[, "ph"]` without drop = FALSE would collapse to an unnamed vector and lose the
        # column on the following cbind(). Harmless for the full-property path (many columns).
        sim_data <- sim_data[, setdiff(names(sim_data), c("ilr1", "ilr2")), drop = FALSE]
        sim_data <- cbind(sim_data, sim_txt)
      }

      sim_data$compname <- row$compname
      sim_data$mukey    <- row$mukey
      sim_data$cokey    <- row$cokey
      sim_data$hzdept_r <- row$hzdept_r
      sim_data$hzdepb_r <- row$hzdepb_r

      # OSD-derived boundary distinctness - optional column, present only when the caller ran
      # attach_osd_boundary_distinctness() on its input first (e.g. simulate_ssurgo_mapunit_draws()).
      if ("bound_sd" %in% names(row)) {
        sim_data$bound_sd <- row$bound_sd
      }

      sim_data$simulation_number <- seq_len(nrow(sim_data))
      sim_data$unique_id <- paste0(row$cokey, "-", sprintf("%02d", sim_data$simulation_number))

      sim_data_out <- append(sim_data_out, list(sim_data))
    }, error = function(e) {
      message("Error for row ", i, ": ", e$message)
    })
  }

  # dplyr::bind_rows() (not base do.call(rbind, ...)) - since one row's texture simulation
  # can fail while a sibling row in the same cokey succeeds (see the texture tryCatch()
  # above), sim_data_out can legitimately mix rows that have sand_total/silt_total/
  # clay_total columns with rows that don't. Base rbind() requires identical columns
  # across all inputs and errors ("numbers of columns of arguments do not match") on that
  # mismatch. bind_rows() unions columns and
  # fills the missing ones with NA, which is the correct semantics here (that row simply
  # has no simulated texture).
  # Preserve the pre-existing "no surviving rows -> NULL" contract (relied on by callers
  # and by this function's own tests) - bind_rows(list()) would otherwise return a 0-row
  # tibble instead of NULL.
  final_df <- if (length(sim_data_out) == 0) NULL else dplyr::bind_rows(sim_data_out)
  return(final_df)
}


# ============================================================================
# merged from core-simulation.R (P1 file reorg)
# ============================================================================

#' @title Soil Profile Depth Simulation
#'
#' @description SSURGO-based soil-profile horizon-depth simulation: top-down/bottom-up
#'   depth simulation, thickness-variability estimation via repeated simulation,
#'   OSD-distinctness-derived boundary perturbation, and
#'   `aqp::SoilProfileCollection`-level orchestration (single profile, whole
#'   collection, parallelized collection, or by-mukey).
#'
#'   This is the first place `aqp`/`SoilProfileCollection` objects are used in
#'   `soilSIM`. The random triangular-distribution sampler these functions rely
#'   on, `tri_dist()`, lives in `R/core-distributions.R`.
#' @name depth_simulation
NULL

#' Query SSURGO Database for Soil Data by Mukey
#'
#' Queries the NRCS SSURGO database (`mapunit`/`component`/`chorizon`/`chfrags`)
#' for horizon morphology, texture, bulk density, water-retention, and
#' rock-fragment-volume data for the given map unit key(s).
#'
#' This deliberately does not reuse `execute_ssurgo_query_working()` (which
#' also joins `mapunit`/`component`/`chorizon`/`chfrags`): that function only
#' ever selects `chf.fragsize_r`, never `chf.fragvol_l/r/h`, so its rock
#' fragment aggregation path (`aggregate_rock_fragment_volume_working()`)
#' currently produces no data. This function's entire purpose is summing
#' `fragvol_l/r/h` per horizon, so delegating would silently drop rock
#' fragment volume - a real regression, not a style difference.
#'
#' @param mukeys A vector of string or numeric map unit keys (mukeys) representing the map units to retrieve data for.
#'
#' @return A data frame with soil horizon data for the given mukeys, including aggregated rock fragment data.
#'
#' @keywords internal
fetch_ssurgo_aws_data <- function(mukeys) {
  # Step 1: Format mukey(s) for SQL
  formatted_mukey <- paste0("(", paste0("'", mukeys, "'", collapse = ","), ")")

  # Step 2: Construct the SQL query
  query <- sprintf("
    SELECT
      -- contextual data
      co.mukey, co.cokey, compname, comppct_l, comppct_r, comppct_h,
      -- horizon morphology
      ch.chkey, hzname, hzdept_l, hzdept_r, hzdept_h, hzdepb_l, hzdepb_r, hzdepb_h, hzthk_l, hzthk_r, hzthk_h,
      -- soil texture parameters
      sandtotal_l, sandtotal_r, sandtotal_h,
      silttotal_l, silttotal_r, silttotal_h,
      claytotal_l, claytotal_r, claytotal_h,
      -- bulk density and water retention parameters
      dbovendry_l, dbovendry_r, dbovendry_h,
      wthirdbar_l, wthirdbar_r, wthirdbar_h,
      wfifteenbar_l, wfifteenbar_r, wfifteenbar_h,
      -- rock fragment volume from chfrags table
      chf.fragvol_l AS rfv_l, chf.fragvol_r AS rfv_r, chf.fragvol_h AS rfv_h
    FROM mapunit AS mu
      INNER JOIN component AS co ON mu.mukey = co.mukey
      INNER JOIN chorizon AS ch ON co.cokey = ch.cokey
      LEFT JOIN chfrags AS chf ON ch.chkey = chf.chkey
    WHERE mu.mukey IN %s
    ORDER BY co.cokey, ch.hzdept_r ASC;", formatted_mukey)

  # Step 3: Submit the query and fetch results
  result <- soilDB::SDA_query(query)

  # Step 4: Group by `chkey` and sum rock fragment values (rfv_l, rfv_r, rfv_h)
  rfv_sum <- result |>
    dplyr::group_by(chkey) |>
    dplyr::summarise(
      rfv_l = sum(rfv_l, na.rm = TRUE),
      rfv_r = sum(rfv_r, na.rm = TRUE),
      rfv_h = sum(rfv_h, na.rm = TRUE)
    ) |>
    dplyr::ungroup()

  # Step 5: Merge the summed values back to the main dataset and remove duplicates
  result <- result |>
    dplyr::select(-dplyr::all_of(c("rfv_l", "rfv_r", "rfv_h"))) |>
    dplyr::left_join(rfv_sum, by = "chkey") |>
    dplyr::distinct()

  # Step 6: Return the final result
  return(result)
}


#--------------------------------------------------------------------------------------------------------
# Functions to run soil profile depth simulations

#' Query OSD Data and Convert Horizon Distinctness to Offset
#'
#' Retrieves Official Series Description (OSD) data for the soil series
#' (component names) present in `horizon_data`, extracts horizon distinctness,
#' and converts the distinctness codes into boundary-offset values via
#' `aqp::hzDistinctnessCodeToOffset()`.
#'
#' @param horizon_data Data frame containing SSURGO horizon data, with at
#'   least `compname` and `hzname` columns.
#'
#' @return A data frame containing the following columns:
#'         - `id`: Unique identifier for the soil profile (the component name).
#'         - `hzname`: The name of the soil horizon.
#'         - `distinctness`: The distinctness code for the horizon (e.g., A, C, G, D).
#'         - `genhz`: The generalized horizon designation.
#'         - `bound_sd`: The calculated offset value from the distinctness code.
#'
#' @keywords internal
query_osd_distinctness <- function(horizon_data) {
  # Subset SSURGO horizon data to keep only compname and hzname
  ssurgo_horizon_data <- horizon_data |>
    dplyr::select(compname, hzname)

  # Ensure 'genhz' column exists in ssurgo_horizon_data
  if (!"genhz" %in% colnames(ssurgo_horizon_data)) {
    ssurgo_horizon_data$genhz <- classify_genhz(ssurgo_horizon_data$hzname)
  }

  # Rename 'compname' to 'id'
  colnames(ssurgo_horizon_data)[colnames(ssurgo_horizon_data) == "compname"] <- "id"

  # Fetch OSD data for the list of soil component names, reusing cached per-series results
  # where available (see fetch_osd_horizons_cached()'s Known quirk note re: id casing) -
  # avoids redundant live soilDB::fetchOSD() calls when the same series appears in multiple
  # profiles/components within or across sessions (e.g. simulate_profile_depths_by_mukey()
  # calls this once per component, and components frequently share a series name).
  osd_horizon_data <- fetch_osd_horizons_cached(unique(horizon_data$compname))

  # Generate generalized horizon codes for the OSD data
  osd_horizon_data$genhz <- classify_genhz(osd_horizon_data$hzname)

  # Identify missing genhz values in ssurgo_horizon_data that are not in osd_horizon_data
  missing_genhz <- setdiff(ssurgo_horizon_data$genhz, osd_horizon_data$genhz)

  # If there are missing genhz values, append them with NA for distinctness
  if (length(missing_genhz) > 0) {
    missing_rows <- ssurgo_horizon_data[ssurgo_horizon_data$genhz %in% missing_genhz, c("id", "hzname", "genhz")]
    missing_rows$distinctness <- NA
    osd_horizon_data <- base::rbind(osd_horizon_data, missing_rows)
    message("Added missing genhz values: ", paste(missing_genhz, collapse = ", "))
  } else {
    message("No missing genhz values found.")
  }

  # Infill missing distinctness values using a custom function
  osd_horizon_data <- infill_missing_distinctness(osd_horizon_data)

  # Convert the distinctness codes to offset values
  osd_horizon_data$bound_sd <- aqp::hzDistinctnessCodeToOffset(osd_horizon_data$distinctness)

  # Rename the columns for clarity
  names(osd_horizon_data) <- c("id", "hzname", "distinctness", "genhz", "bound_sd")

  return(osd_horizon_data)
}

#' Attach OSD-Derived Boundary Distinctness to Horizon Data (per-genhz `bound_sd`)
#'
#' `query_osd_distinctness()`'s `bound_sd` (OSD boundary-distinctness offset: Abrupt/Clear/
#' Gradual/Diffuse converted to a numeric scale via `aqp::hzDistinctnessCodeToOffset()`) is
#' also used by `simulate_and_perturb_soil_profiles()` as `aqp::perturb()`'s `boundary.attr`
#' for horizon-depth perturbation. This function exposes the same genhz-averaged `bound_sd`
#' value on arbitrary horizon data (via the same
#' `dplyr::group_by(id, genhz) |> dplyr::summarise(bound_sd = mean(bound_sd, na.rm = TRUE))`
#' pattern), so it can also gate the vertical-correlation depth kernel.
#'
#' @param hz_data A data frame with at least `compname`, `hzname`, and `genhz` columns (the same
#'   shape `simulate_ssurgo_mapunit_draws()` builds right after `classify_genhz()`).
#'
#' @return `hz_data` with a `bound_sd` numeric column added/overwritten (one value per
#'   `compname`/`genhz` combination, `NA` for any combination the OSD lookup couldn't resolve).
#'   On OSD lookup failure (e.g. no network access - `query_osd_distinctness()` calls
#'   `soilDB::fetchOSD()`), degrades to an all-`NA` `bound_sd` column rather than erroring, so
#'   callers (and the kernel gating, which treats missing distinctness data as
#'   "no gating - fall back to the plain kernel") keep working without OSD access.
#'
#' @keywords internal
attach_osd_boundary_distinctness <- function(hz_data) {
  required_cols <- c("compname", "hzname", "genhz")
  if (!all(required_cols %in% names(hz_data))) {
    stop("hz_data must have compname, hzname, and genhz columns")
  }

  distinctness_data <- tryCatch({
    query_osd_distinctness(hz_data)
  }, error = function(e) {
    message("attach_osd_boundary_distinctness(): OSD distinctness lookup failed (",
            e$message, ") - proceeding with bound_sd = NA (no depth-kernel discontinuity gating).")
    NULL
  })

  if (is.null(distinctness_data)) {
    hz_data$bound_sd <- NA_real_
    return(hz_data)
  }

  # soilDB::fetchOSD() (called inside query_osd_distinctness() -> fetch_osd_horizons_cached())
  # returns its id column UPPER-CASED regardless of the requested case (documented quirk, see
  # fetch_osd_horizons_cached()'s own "Known quirk preserved" section) - the ORIGINAL caller of
  # query_osd_distinctness(), simulate_and_perturb_soil_profiles(), was documented as unaffected
  # by this because it joins by genhz only, never by id. This function DOES need to join by
  # compname (to get a per-cokey-relevant value, not just a global per-genhz average across every
  # series in the input), so it must match case-insensitively: title-case "Placentia" from
  # SSURGO and "PLACENTIA" from fetchOSD() must join.
  bound_lut <- distinctness_data |>
    dplyr::mutate(compname_upper = toupper(id)) |>
    dplyr::group_by(compname_upper, genhz) |>
    dplyr::summarise(bound_sd = mean(bound_sd, na.rm = TRUE), .groups = "drop")

  # left_join preserves hz_data's row count/order exactly (one bound_sd value broadcast across
  # every row sharing a compname/genhz combination) - a genhz/compname combination missing from
  # the OSD data (or dropped by NA-averaging when every OSD match for it was itself NA) yields
  # bound_sd = NA for those rows, the same graceful-missing-data contract as the rest of this
  # function.
  hz_data |>
    dplyr::mutate(compname_upper = toupper(compname)) |>
    dplyr::left_join(bound_lut, by = c("compname_upper", "genhz")) |>
    dplyr::select(-compname_upper)
}

#' Fetch OSD horizon distinctness data for a set of soil series names, with per-series caching
#'
#' `soilDB::fetchOSD()` is a live network call. Multiple profiles/components in the same
#' mukey/AOI frequently share the same series name (e.g. `simulate_profile_depths_by_mukey()`
#' calls `query_osd_distinctness()` once per component), so fetching unconditionally on every
#' call re-requests identical data redundantly. This wrapper caches each series' raw
#' `id`/`hzname`/`distinctness` rows on disk - the same `cache_get()`/`cache_set()` mechanism
#' `R/cache.R` already uses for AOI-keyed raster fetches, just with a series-name-keyed
#' cache key instead - fetching only the series not already cached.
#'
#' @section Known quirk preserved (not a bug introduced here): `soilDB::fetchOSD()` returns its
#'   `id` column upper-cased regardless of the requested case (e.g. requesting `"amador"`
#'   returns rows with `id == "AMADOR"`) - confirmed via a live call before writing this
#'   function. Splitting a combined multi-series fetch into per-series cache entries is done via
#'   case-insensitive matching on `id` to account for this. This casing difference was already
#'   present in `soilDB::fetchOSD()`'s output before this caching layer existed, and downstream
#'   code (`simulate_and_perturb_soil_profiles()`) was already unaffected by it, since it joins
#'   `query_osd_distinctness()`'s output by `genhz` only, never by `id`.
#'
#' @param compnames Character vector of soil series names, as would be passed directly to
#'   `soilDB::fetchOSD()`.
#' @return A data frame with columns `id`, `hzname`, `distinctness` - one row per horizon, across
#'   all requested series (any already-cached plus any newly fetched).
#' @keywords internal
fetch_osd_horizons_cached <- function(compnames) {
  compnames <- unique(compnames)
  keys <- paste0("osd_series_", vapply(toupper(compnames), digest::digest, character(1)))
  cached <- stats::setNames(lapply(keys, cache_get), compnames)

  missing <- compnames[vapply(cached, is.null, logical(1))]
  if (length(missing) > 0) {
    s <- soilDB::fetchOSD(missing)

    if (!all(c("distinctness", "top", "bottom", "hzname", "id") %in% names(s@horizons))) {
      stop("Necessary fields (distinctness, id, top, bottom, hzname) are not available in the OSD data.")
    }

    fetched <- s@horizons[, c("id", "hzname", "distinctness")]
    for (cn in missing) {
      rows <- fetched[toupper(fetched$id) == toupper(cn), , drop = FALSE]
      key <- paste0("osd_series_", digest::digest(toupper(cn)))
      cache_set(key, "osd_series", rows)
      cached[[cn]] <- rows
    }
  }

  do.call(rbind, cached)
}


#' Infill Missing Distinctness Values for Horizons
#'
#' Fills in missing `distinctness` values for common horizon names (O, A, B, C, R, Cr)
#' based on typical boundary characteristics observed in the field.
#'
#' ## Default Boundary Distinctness Assignments:
#' - **O Horizons (Organic Layers):** diffuse
#' - **A Horizons (Topsoil):** clear
#' - **B Horizons (Subsoil):** gradual
#' - **C Horizons (Parent Material):** gradual
#' - **Cr Horizons (Weathered Bedrock):** gradual
#' - **R Horizons (Bedrock):** abrupt
#'
#' @param horizon_data A data frame containing horizon data. It must include `hzname` and `distinctness` columns.
#'
#' @return A data frame with missing `distinctness` values infilled based on the horizon name or generalized horizon group.
#'
#' @keywords internal
infill_missing_distinctness <- function(horizon_data) {
  # List of default boundary distinctness values for common horizons
  default_distinctness <- list(
    "O"  = "diffuse",  # Organic horizon
    "A"  = "clear",    # Topsoil
    "E"  = "clear",    # Subsoil
    "B"  = "gradual",  # Subsoil
    "C"  = "gradual",  # Parent material
    "Cr" = "gradual",  # Weathered bedrock
    "R"  = "abrupt"    # Bedrock
  )

  # Assign a generalized horizon group (genhz) if it's not already present
  if (!"genhz" %in% colnames(horizon_data)) {
    horizon_data$genhz <- classify_genhz(horizon_data$hzname)
  }

  # Function to assign default distinctness based on horizon name or generalized horizon group
  fill_default_distinctness <- function(hzname, genhz, distinctness) {
    if (!is.na(distinctness)) {
      return(distinctness)  # If distinctness exists, return it
    }
    # Use default distinctness based on genhz first, then hzname if available
    if (!is.na(genhz) && genhz %in% names(default_distinctness)) {
      return(default_distinctness[[genhz]])
    } else if (!is.na(hzname) && hzname %in% names(default_distinctness)) {
      return(default_distinctness[[hzname]])
    } else {
      return(NA)  # Return NA if no match is found
    }
  }

  # Apply the default distinctness infilling using mapply
  horizon_data$distinctness <- mapply(
    fill_default_distinctness,
    horizon_data$hzname,
    horizon_data$genhz,
    horizon_data$distinctness
  )

  return(horizon_data)
}


#' Infill Missing Depth Variability
#'
#' Infills missing depth variability values in a horizon data frame by
#' replacing missing top and bottom depth bounds with the representative
#' depth value plus or minus 2, ensuring that the resulting values are not
#' negative.
#'
#' @param horizon_data A data frame containing horizon depth data. It must include the following columns:
#'   - `hzdept_r`: Representative top depth.
#'   - `hzdept_l`: Top low value (may be missing).
#'   - `hzdept_h`: Top high value (may be missing).
#'   - `hzdepb_r`: Representative bottom depth.
#'   - `hzdepb_l`: Bottom low value (may be missing).
#'   - `hzdepb_h`: Bottom high value (may be missing).
#'
#' @return A data frame with missing depth variability values filled in. Missing top low values are replaced with
#'   `hzdept_r - 2` and missing top high values with `hzdept_r + 2`. Similarly, missing bottom low values are replaced
#'   with `hzdepb_r - 2` and missing bottom high values with `hzdepb_r + 2`. All values are ensured to be non-negative.
#'
#' @keywords internal
infill_missing_depth_variability <- function(horizon_data) {
  hw <- RANGE_FALLBACK_HALFWIDTH

  # Missing low/high depth bounds -> representative -/+ hw cm, ensuring non-negative results.
  horizon_data$hzdept_l <- ifelse(
    is.na(horizon_data$hzdept_l),
    pmax(horizon_data$hzdept_r - hw, 0),
    horizon_data$hzdept_l
  )

  horizon_data$hzdept_h <- ifelse(
    is.na(horizon_data$hzdept_h),
    pmax(horizon_data$hzdept_r + hw, 0),
    horizon_data$hzdept_h
  )

  horizon_data$hzdepb_l <- ifelse(
    is.na(horizon_data$hzdepb_l),
    pmax(horizon_data$hzdepb_r - hw, 0),
    horizon_data$hzdepb_l
  )

  horizon_data$hzdepb_h <- ifelse(
    is.na(horizon_data$hzdepb_h),
    pmax(horizon_data$hzdepb_r + hw, 0),
    horizon_data$hzdepb_h
  )

  return(horizon_data)
}


#' Top-down Simulation of Soil Profile
#'
#' Simulates a soil profile by calculating the top and bottom depths for each
#' horizon in a top-down manner. The first horizon always starts at 0 cm, and
#' the top of each subsequent horizon is set equal to the bottom depth of the
#' previous horizon. The bottom depth of each horizon is simulated using a
#' triangular distribution via `tri_dist()`, defined by a lower bound
#' (`hzdepb_l`), a representative value (`hzdepb_r`), and an upper bound
#' (`hzdepb_h`).
#'
#' @param horizon_data A data frame containing horizon data. It should include the following columns:
#'   - `hzname`: Name of the horizon.
#'   - `hzdepb_l`: Lower bound for the bottom depth of the horizon.
#'   - `hzdepb_r`: Representative bottom depth of the horizon.
#'   - `hzdepb_h`: Upper bound for the bottom depth of the horizon.
#'
#' @return A data frame with the following columns:
#'   - `hzname`: The name of the horizon.
#'   - `top`: The simulated top depth of the horizon.
#'   - `bottom`: The simulated bottom depth of the horizon.
#'
#' @keywords internal
simulate_soil_profile_top_down <- function(horizon_data) {
  n <- nrow(horizon_data)

  # Initialize a result data frame with hzname, top, and bottom columns
  result <- data.frame(
    hzname = horizon_data$hzname,
    top = rep(NA, n),
    bottom = rep(NA, n)
  )

  # Set the top of the first horizon to 0
  result$top[1] <- 0

  for (i in 1:n) {
    if (i > 1) {
      # Set the top of the current horizon to be the bottom of the previous horizon
      result$top[i] <- result$bottom[i - 1]
    }

    # Determine bounds for bottom depth simulation
    bottom_low  <- max(horizon_data$hzdepb_l[i], result$top[i])
    bottom_high <- max(bottom_low, horizon_data$hzdepb_h[i])
    bottom_mode <- min(max(horizon_data$hzdepb_r[i], bottom_low), bottom_high)

    # Generate a single random sample using the tri_dist function
    bottom_depth <- tri_dist(n = 1, a = bottom_low, b = bottom_high, c = bottom_mode)
    bottom_depth <- round(bottom_depth)

    # Ensure the bottom depth is not less than the top depth
    if (bottom_depth < result$top[i]) {
      bottom_depth <- result$top[i]
    }
    result$bottom[i] <- bottom_depth
  }

  return(result)
}


#' Bottom-up Simulation of Soil Profile
#'
#' Simulates a soil profile from the bottom upward by calculating the top and
#' bottom depths for each horizon. The simulation starts from the bottom
#' horizon and proceeds upward, setting the bottom depth of each horizon to
#' the top depth of the horizon immediately below. For each horizon, the
#' bottom and top depths are simulated using a triangular distribution via
#' `tri_dist()`. The function ensures that there are no gaps between horizons
#' and that the top of the uppermost horizon is set to 0.
#'
#' @param horizon_data A data frame containing horizon data. The data frame must include the following columns:
#'   - `hzname`: The name of the horizon.
#'   - `hzdept_l`: The lower bound for the top depth of the horizon.
#'   - `hzdept_r`: The representative top depth of the horizon.
#'   - `hzdept_h`: The upper bound for the top depth of the horizon.
#'   - `hzdepb_l`: The lower bound for the bottom depth of the horizon.
#'   - `hzdepb_r`: The representative bottom depth of the horizon.
#'   - `hzdepb_h`: The upper bound for the bottom depth of the horizon.
#'
#' @return A data frame with the following columns:
#'   - `hzname`: The name of the horizon.
#'   - `top`: The simulated top depth of the horizon.
#'   - `bottom`: The simulated bottom depth of the horizon.
#'
#' @keywords internal
simulate_soil_profile_bottom_up <- function(horizon_data) {
  n <- nrow(horizon_data)

  # Initialize a result data frame with hzname, top, and bottom columns
  result <- data.frame(
    hzname = horizon_data$hzname,
    top = rep(NA, n),
    bottom = rep(NA, n)
  )

  # Process horizons from bottom to top
  for (i in n:1) {
    # Simulate bottom depth using tri_dist
    bottom_low  <- horizon_data$hzdepb_l[i]
    bottom_high <- horizon_data$hzdepb_h[i]
    bottom_mode <- horizon_data$hzdepb_r[i]

    # Generate a single random sample using tri_dist
    bottom_depth <- tri_dist(n = 1, a = bottom_low, b = bottom_high, c = bottom_mode)
    bottom_depth <- round(bottom_depth)

    if (i == n) {
      # For the bottom horizon, simulate top depth
      top_low  <- horizon_data$hzdept_l[i]
      top_high <- min(horizon_data$hzdept_h[i], bottom_depth)
      top_mode <- min(max(horizon_data$hzdept_r[i], top_low), top_high)

      # Generate a single random sample using tri_dist
      top_depth <- tri_dist(n = 1, a = top_low, b = top_high, c = top_mode)
      top_depth <- round(top_depth)

      # Ensure top depth is not greater than bottom depth
      if (top_depth > bottom_depth) {
        top_depth <- bottom_depth
      }
    } else {
      # For other horizons, set bottom depth to the top depth of the horizon below
      bottom_depth <- result$top[i + 1]

      # Simulate top depth using tri_dist
      top_low  <- horizon_data$hzdept_l[i]
      top_high <- min(horizon_data$hzdept_h[i], bottom_depth)
      top_mode <- min(max(horizon_data$hzdept_r[i], top_low), top_high)

      # Generate a single random sample using tri_dist
      top_depth <- tri_dist(n = 1, a = top_low, b = top_high, c = top_mode)
      top_depth <- round(top_depth)

      # Ensure top depth is not greater than bottom depth
      if (top_depth > bottom_depth) {
        top_depth <- bottom_depth
      }
    }

    # Assign the simulated depths to the result
    result$top[i]    <- top_depth
    result$bottom[i] <- bottom_depth
  }

  # Ensure the top depth of the uppermost horizon is 0
  result$top[1] <- 0

  # Adjust the bottom depth of the first horizon if necessary
  if (result$bottom[1] < result$top[1]) {
    result$bottom[1] <- result$top[1]
  }

  # Ensure no gaps between horizons
  for (i in 1:(n - 1)) {
    if (result$bottom[i] != result$top[i + 1]) {
      result$bottom[i] <- result$top[i + 1]
    }
  }

  return(result)
}


#' Simulate Soil Profile Thickness
#'
#' Simulates the thickness of soil horizons for a given soil profile using
#' both top-down and bottom-up simulation methods. It performs multiple
#' simulations and calculates the standard deviation of horizon thickness
#' (i.e., the difference between representative bottom and top depths)
#' across the simulations.
#'
#' @param horizon_data A dataframe containing horizon data. It must include the following columns:
#'        - `hzname`: The name of each horizon (e.g., A, B, C).
#'        - `hzdept_r`: The representative top depth of each horizon.
#'        - `hzdepb_r`: The representative bottom depth of each horizon.
#'        - `hzdept_l`: The low estimate for the top depth of each horizon.
#'        - `hzdept_h`: The high estimate for the top depth of each horizon.
#'        - `hzdepb_l`: The low estimate for the bottom depth of each horizon.
#'        - `hzdepb_h`: The high estimate for the bottom depth of each horizon.
#' @param n_simulations Integer, number of simulations to perform (default is 500).
#'
#' @return A dataframe containing the following columns:
#'         - `hzname`: The horizon name.
#'         - `top`: The representative top depth of the horizon (from `hzdept_r`).
#'         - `bottom`: The representative bottom depth of the horizon (from `hzdepb_r`).
#'         - `thickness_sd`: The standard deviation of the horizon thickness (bottom - top) across the simulations.
#'
#' @keywords internal
simulate_soil_profile_thickness <- function(horizon_data, n_simulations = 500) {
  # Step 1: Infill missing depth variability values
  horizon_data <- infill_missing_depth_variability(horizon_data)

  results_list <- list()

  # Loop over the number of simulations
  for (sim in 1:n_simulations) {
    # Randomly select the simulation direction: top-down or bottom-up
    if (stats::runif(1) < 0.5) {
      # Simulate profile using the top-down method
      result <- simulate_soil_profile_top_down(horizon_data)
      method <- "Top-Down"
    } else {
      # Simulate profile using the bottom-up method
      result <- simulate_soil_profile_bottom_up(horizon_data)
      method <- "Bottom-Up"
    }

    # Calculate thickness for each horizon as the difference between bottom and top depths
    result$Thickness <- result$bottom - result$top

    # Record the simulation method and simulation number
    result$Method <- method
    result$Simulation <- sim

    # Store the current simulation result
    results_list[[sim]] <- result
  }

  # Combine all simulation results into a single dataframe
  all_results <- do.call(rbind, results_list)

  # Reorder columns for clarity: Simulation, Method, hzname, top, bottom, Thickness
  all_results <- all_results[, c("Simulation", "Method", "hzname", "top", "bottom", "Thickness")]

  # Group results by horizon name and calculate the standard deviation of thickness for each horizon
  horizon_sd <- all_results |>
    dplyr::group_by(hzname) |>
    dplyr::summarise(thickness_sd = round(stats::sd(Thickness), 2)) |>
    dplyr::ungroup()

  # Merge the representative top and bottom depths from horizon_data
  summarized_results <- merge(
    horizon_data[, c("hzname", "hzdept_r", "hzdepb_r")],
    horizon_sd,
    by = "hzname"
  )

  # Rename columns for consistency: hzdept_r -> top, hzdepb_r -> bottom
  summarized_results <- summarized_results |>
    dplyr::rename(
      top = hzdept_r,
      bottom = hzdepb_r
    )

  # Return the summarized results
  return(summarized_results)
}


#' Simulate and Perturb Soil Profiles
#'
#' Takes a `SoilProfileCollection` object and perturbs its horizons to
#' simulate variability in soil profile thickness and boundary depths. The
#' simulation is performed in two steps: first, horizon thickness is
#' perturbed based on simulated thickness variability (using a top-down or
#' bottom-up simulation method), and then the boundary depths are perturbed
#' using distinctness-derived offsets. If the profile contains only one
#' horizon, the function simply replicates the profile.
#'
#' @section Known limitation:
#' This function requires `soil_profile`'s horizons to already carry a
#' `sim_comppct` column (it derives the number of simulations to run from
#' `unique(horizons(soil_profile)$sim_comppct)`). `simulate_component_composition()`
#' (`R/core-simulation.R`) produces this column, but at component
#' (`cokey`) grain, not horizon grain - callers must
#' `dplyr::left_join(horizon_data, simulate_component_composition(component_data), by = "cokey")`
#' before calling this function; it will still error with a missing-column
#' condition if that join hasn't been done first.
#'
#' @param soil_profile A SoilProfileCollection object containing soil profile data.
#'
#' @return A SoilProfileCollection object with perturbed horizon depths.
#'
#' @keywords internal
simulate_and_perturb_soil_profiles <- function(soil_profile) {
  # Step 1: Extract horizons from the soil profile and select relevant columns
  soil_profile@horizons <- aqp::horizons(soil_profile) |>
    dplyr::select(mukey, id, cokey, compname, hzname, sim_comppct,
                  hzdept_l, hzdept_r, hzdept_h, hzdepb_l, hzdepb_r, hzdepb_h, hzthk_l, sim_comppct)

  horizon_data <- aqp::horizons(soil_profile)
  horizon_data$genhz <- classify_genhz(horizon_data$hzname)
  n_simulations <- unique(horizon_data$sim_comppct)  # Set number of simulations based on sim_comppct

  # Check if there's only one horizon (e.g., R horizon only)
  if (nrow(horizon_data) == 1) {
    message("Only one horizon in the profile. Skipping perturbation for this profile.")

    # Return the original soil profile replicated n_simulations times (no perturbation needed)
    replicated_profiles <- list()
    for (i in 1:n_simulations) {
      # Create a copy of the soil profile and assign a unique profile ID
      new_profile <- soil_profile
      aqp::profile_id(new_profile) <- paste0(aqp::profile_id(new_profile), "_sim_", i)
      replicated_profiles[[i]] <- new_profile
    }
    simulated_profiles <- aqp::combine(replicated_profiles)
    return(simulated_profiles)
  }

  # Step 1 (continued): Simulate soil profile thickness
  simulated_thickness <- simulate_soil_profile_thickness(horizon_data, n_simulations)
  simulated_thickness$genhz <- classify_genhz(simulated_thickness$hzname)

  # Step 2: Query OSD distinctness and get bound_sd values
  distinctness_data <- query_osd_distinctness(horizon_data)
  bound_lut <- distinctness_data |>
    dplyr::group_by(id, genhz) |>
    dplyr::summarise(bound_sd = mean(bound_sd, na.rm = TRUE)) |>
    dplyr::ungroup()

  # Step 3: Merge the simulated thickness data and distinctness data
  combined_data <- merge(simulated_thickness, bound_lut, by = "genhz") |>
    dplyr::arrange(top)

  # Step 4: Add thickness standard deviation and bound_sd to the soil profile horizons
  horizon_data <- horizon_data |>
    dplyr::left_join(combined_data |> dplyr::select(hzname, thickness_sd, bound_sd) |> dplyr::distinct(), by = "hzname")
  soil_profile@horizons <- horizon_data |>
    dplyr::select(mukey, id, cokey, compname, hzname, sim_comppct, hzdept_r, hzdepb_r, thickness_sd, bound_sd)

  # Step 3: First perturbation (horizon thickness)
  perturbed_profiles_thickness <- aqp::perturb(
    soil_profile,
    n = n_simulations,
    thickness.attr = "thickness_sd"
  )

  # Step 4: Second perturbation (boundary distinctness)
  list_of_profiles <- list()
  for (i in seq_along(perturbed_profiles_thickness)) {
    single_profile <- perturbed_profiles_thickness[i]
    perturbed_profile <- aqp::perturb(
      single_profile,
      n = 1,
      boundary.attr = "bound_sd",
      id = single_profile@site$id
    )
    list_of_profiles[[i]] <- perturbed_profile
  }

  simulated_profiles <- aqp::combine(list_of_profiles)

  # Step 5: Identify out-of-range horizons and adjust depths
  adjust_out_of_range_profiles <- function(simulated_profiles, horizon_data) {
    sim_horizons <- aqp::horizons(simulated_profiles) |>
      dplyr::select(id, mukey, cokey, compname, hzname, top = hzdept_r, bottom = hzdepb_r)

    merged_data <- merge(
      sim_horizons,
      horizon_data[, c("hzname", "hzdept_l", "hzdepb_h")],
      by = "hzname",
      all.x = TRUE
    )

    merged_data <- merged_data |>
      dplyr::rowwise() |>
      dplyr::mutate(
        top = dplyr::if_else(is.na(hzdept_l), top, pmax(top, hzdept_l)),
        bottom = dplyr::if_else(is.na(hzdepb_h), bottom, pmin(bottom, hzdepb_h))
      ) |>
      dplyr::ungroup()

    simulated_profiles_hzdata <- aqp::horizons(simulated_profiles) |>
      dplyr::left_join(merged_data |> dplyr::select(hzname, id, top, bottom), by = c("id", "hzname"))
    simulated_profiles@horizons <- simulated_profiles_hzdata

    return(simulated_profiles)
  }

  simulated_profiles <- adjust_out_of_range_profiles(simulated_profiles, horizon_data)

  # Step 6: Return the adjusted simulated profiles
  return(simulated_profiles)
}


#' Run Soil Profile Depth Simulations for a SoilProfileCollection
#'
#' Simulates soil profile depths for each profile in a `SoilProfileCollection`
#' by applying `simulate_and_perturb_soil_profiles()` to each individual
#' profile, then combines the results into a single `SoilProfileCollection`.
#'
#' @param soil_collection A SoilProfileCollection object containing soil profile data.
#' @param seed An integer value for setting the random seed (default is 123) for reproducible simulations.
#'
#' @return A SoilProfileCollection object containing the simulated soil profiles.
#'
#' @examples
#' \dontrun{
#'   simulated_collection <- simulate_profile_depths_by_collection(soil_collection, seed = 123)
#' }
#'
#' @export
simulate_profile_depths_by_collection <- function(soil_collection, seed = 123) {

  # Check if input is a valid SoilProfileCollection
  if (!inherits(soil_collection, "SoilProfileCollection")) {
    stop("Input must be a SoilProfileCollection.")
  }

  # Set seed for reproducibility
  set.seed(seed)

  # Initialize an empty list to store all simulated profiles
  all_simulated_profiles <- list()

  # Loop over each profile in the SoilProfileCollection
  for (i in seq_along(soil_collection)) {
    # Extract the current soil profile
    soil_profile <- soil_collection[i, ]

    # Simulate and perturb the current soil profile
    simulated_profiles <- simulate_and_perturb_soil_profiles(soil_profile)

    # Append the simulated profiles to the list
    all_simulated_profiles[[i]] <- simulated_profiles
  }

  # Combine all simulated profiles into a single SoilProfileCollection
  combined_simulated_profiles <- aqp::combine(all_simulated_profiles)

  # Return the combined simulated profiles
  return(combined_simulated_profiles)
}


#' Run Soil Profile Depth Simulations in Parallel for a SoilProfileCollection
#'
#' Simulates soil profile depths for each profile in a `SoilProfileCollection`
#' in parallel using the `future`/`future.apply` packages, then combines the
#' results into a single `SoilProfileCollection`.
#'
#' @section Note on parallel workers:
#' Dispatched via `run_parallel_lapply()` (`R/parallel.R`), this package's shared
#' `future`/`future.apply` helper. `future::multisession` workers are fresh R processes;
#' `future`/`globals` auto-detect that `simulate_single_profile()` calls a `soilSIM`-namespaced
#' function and attach the package in each worker automatically - this only works when `soilSIM`
#' is actually installed and attached in the calling session, not merely `devtools::load_all()`'d.
#' The caller's own `future::plan()` (if any) is restored afterward regardless of success/failure.
#'
#' @param soil_collection A SoilProfileCollection object containing soil profile data.
#' @param seed An integer value to set the random seed for reproducibility (default is 123).
#' @param n_cores Integer, number of cores to use for parallel processing (default is 6).
#'
#' @return A SoilProfileCollection object containing the simulated and perturbed soil profiles.
#'
#' @examples
#' \dontrun{
#'   simulated_profiles <- simulate_profile_depths_by_collection_parallel(
#'     soil_collection, seed = 123, n_cores = 6
#'   )
#'   print(simulated_profiles)
#' }
#'
#' @export
simulate_profile_depths_by_collection_parallel <- function(soil_collection, seed = 123, n_cores = 6) {

  # Check if input is a valid SoilProfileCollection
  if (!inherits(soil_collection, "SoilProfileCollection")) {
    stop("Input must be a SoilProfileCollection.")
  }

  # Set seed for reproducibility
  set.seed(seed)

  # Define a function to simulate a single soil profile
  simulate_single_profile <- function(i) {
    soil_profile <- soil_collection[i, ]  # Select the ith profile
    simulate_and_perturb_soil_profiles(soil_profile)
  }

  # Run the simulation in parallel via run_parallel_lapply() (R/parallel.R) - future_seed =
  # TRUE gives future's parallel-safe RNG streams, needed for this function's own seed = ...
  # reproducibility contract (simulate_and_perturb_soil_profiles() perturbs randomly).
  # catch_errors = FALSE: this function's own contract is "graceful NULL on error, no sequential
  # fallback" (unlike the other run_parallel_lapply() call sites), so error handling stays in this
  # function's own tryCatch below rather than the helper's fallback-to-lapply default.
  result <- tryCatch({
    all_simulated_profiles <- run_parallel_lapply(
      seq_along(soil_collection), simulate_single_profile,
      n_cores = n_cores, future_seed = TRUE,
      op_name = "soil profile depth simulation",
      catch_errors = FALSE
    )

    aqp::combine(all_simulated_profiles)
  }, error = function(e) {
    handle_workflow_error(e, "soil profile depth simulation", "warn")
    # Optionally return NULL to indicate failure
    return(NULL)
  })

  return(result)
}


#' Simulate Profile Depths by Mukey
#'
#' Queries soil data for a given map unit key (mukey) using the SSURGO
#' database, sets up the corresponding soil profile data, and then simulates
#' and perturbs the soil profiles for that mukey via
#' `simulate_and_perturb_soil_profiles()`.
#'
#' @param mukey A string or numeric value representing the map unit key to query.
#' @param n_simulations Integer, the number of triangular draws per component used to derive
#'   `sim_comppct` via `simulate_component_composition()` (default is 100).
#' @param seed An integer to set the random seed for reproducibility (default is 123).
#'
#' @return A SoilProfileCollection object containing the simulated and perturbed soil profiles.
#'
#' @keywords internal
simulate_profile_depths_by_mukey <- function(mukey, n_simulations = 100, seed = 123) {
  # Step 1: Query the data based on the mukey
  mu_data <- fetch_ssurgo_aws_data(mukeys = mukey)

  if (is.null(mu_data)) {
    stop("No data found for the provided mukey.")
  }

  # Step 1.5: Derive sim_comppct (simulate_component_composition() produces it at cokey grain) and join
  # it onto every horizon row by cokey - simulate_and_perturb_soil_profiles() requires it at
  # horizon grain. Mirrors the already-verified pattern in R/adapter-ssurgo-simulate.R's
  # simulate_ssurgo_mapunit_draws().
  component_data <- simulate_component_composition(mu_data, n_simulations = n_simulations)
  mu_data <- dplyr::left_join(
    mu_data, component_data[, c("cokey", "sim_comppct")],
    by = "cokey"
  )

  # Step 2: Set up the soil profile data
  set.seed(seed)  # For reproducibility
  # `id` (not `compname`) is the profile-ID column convention the rest of
  # this file assumes (query_osd_distinctness() renames compname -> id
  # itself, and adjust_out_of_range_profiles()/evaluate_simulated_depths()
  # both hardcode a literal `id` column) - derive it before calling
  # aqp::depths<- so downstream horizon selects don't fail with
  # "object 'id' not found".
  mu_data$id <- mu_data$compname
  aqp::depths(mu_data) <- id ~ hzdept_r + hzdepb_r
  aqp::hzdesgnname(mu_data) <- "hzname"

  # Initialize an empty list to store all simulated profiles
  all_simulated_profiles <- list()

  # Step 3: Loop over each profile in mu_data
  for (i in seq_along(mu_data)) {
    # Select the current profile
    soil_profile <- mu_data[i, ]

    # Step 4: Simulate and perturb the current soil profile
    simulated_profiles <- simulate_and_perturb_soil_profiles(soil_profile)

    # Append the simulated profiles to the list
    all_simulated_profiles[[i]] <- simulated_profiles
  }

  # Step 5: Combine all simulated profiles into a single SoilProfileCollection
  combined_simulated_profiles <- aqp::combine(all_simulated_profiles)

  # Step 6: Return the combined simulated profiles
  return(combined_simulated_profiles)
}


#' Evaluate Simulated Profile Depths
#'
#' Evaluates whether the simulated soil profile depths fall outside the
#' specified range. Extracts the simulated top and bottom depths from a
#' `SoilProfileCollection`, merges these with the original horizon data, and
#' flags horizons where the simulated top is less than the lower bound or
#' the simulated bottom exceeds the upper bound.
#'
#' @param simulated_profiles A SoilProfileCollection object containing simulated soil profiles.
#' @param horizon_data A data frame with original horizon data, including columns:
#'        - `hzname`: Horizon name.
#'        - `hzdept_l`: Lower bound for the top depth.
#'        - `hzdepb_h`: Upper bound for the bottom depth.
#'
#' @return A data frame containing rows (horizons) where the simulated depths are out of range.
#'
#' @keywords internal
evaluate_simulated_depths <- function(simulated_profiles, horizon_data) {
  # Extract horizons from the simulated profiles
  sim_horizons <- aqp::horizons(simulated_profiles) |>
    dplyr::select(mukey, cokey, compname, hzname, top = hzdept_r, bottom = hzdepb_r)

  # Merge simulated horizons with original horizon data for comparison
  merged_data <- merge(
    sim_horizons,
    horizon_data[, c("hzname", "hzdept_l", "hzdepb_h")],
    by = "hzname",
    all.x = TRUE
  )

  # Check if simulated top and bottom depths fall outside the specified range
  merged_data$out_of_range <- with(merged_data, (top < hzdept_l) | (bottom > hzdepb_h))

  # Return rows where simulated depths are out of range
  out_of_range_data <- merged_data[merged_data$out_of_range == TRUE, ]

  return(out_of_range_data)
}
