#' @title Component-Composition and Correlated-Triangular Property Simulation
#'
#' @description Component-composition simulation (`sim_component_comp()`), correlated
#'   triangular-distribution sampling (`simulate_correlated_triangular()`), and per-cokey
#'   flexible-property simulation (`simulate_cokey_generalized()`), plus two small standalone
#'   horizon-data helpers (`remove_organic_layer()`, `slice_and_aggregate_soil_data()`).
#'
#'   `simulate_correlated_triangular()` draws uncorrelated normals with `stats::rnorm()` (no
#'   `MASS` dependency), and `simulate_cokey_generalized()`'s texture step uses
#'   `R/distributions.R`'s `ilr_forward()`/`ilr_inverse()` (no `compositions` dependency).
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
#' @export
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
#' @export
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
#' @section Grain mismatch with `R/depth-simulation.R`:
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
sim_component_comp <- function(data, n_simulations = 1000) {

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
#' @export
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
#' @export
calculate_mode <- function(x) {
  ux <- sort(unique(x))
  counts <- tabulate(match(x, ux))
  ux[which.max(counts)]
}

#' Van Bemmelen-factor conversion from SSURGO organic matter % to an estimated soil organic
#' carbon %
#'
#' SSURGO's `chorizon` table reports organic matter (`om_l/_r/_h`), not organic carbon - but
#' `simulate_cokey_generalized()`'s simulated `soc` column (and the KSSL reference correlation
#' matrix's `soc` row/column, see `R/kssl-reference-correlations.R`'s `.kssl_property_name_map`
#' docs) both need an actual SOC-scale value to be fused against SOLUS100's real `soc` variable
#' (organic carbon). `SOC = OM / 1.724` is the standard Van Bemmelen-factor approximation (the
#' same factor `remarginalized_awc()`'s `soc_to_om` argument uses in the opposite direction,
#' `R/raster-fusion-bridge.R`). Applying this conversion needs no change to the KSSL correlation
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
#'   `sim_component_comp()`).
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
#' @export
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
      # when correctly matched (see R/kssl-reference-correlations.R) - chol() inside
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
          ilr1_lrh = c(min(ilr1_vals), calculate_mode(ilr1_vals), max(ilr1_vals)),
          ilr2_lrh = c(min(ilr2_vals), calculate_mode(ilr2_vals), max(ilr2_vals))
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
