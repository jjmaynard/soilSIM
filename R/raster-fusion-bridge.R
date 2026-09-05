#' @title Per-Pixel Ensemble Re-Marginalization (raster-fusion <-> tabular MC bridge)
#'
#' @description Transforms the per-mukey joint Monte Carlo ensemble from
#'   `extract_mukey_joint_ensemble()` (`R/ssurgo-simulation.R`) onto the SOLUS grid so that each
#'   pixel's marginal distributions match `run_stage1_fusion()`'s SOLUS-fused posterior, while
#'   preserving the ensemble's empirical copula (cross-property and cross-depth-window rank
#'   structure). This is the "per-pixel" answer to wiring raster fusion into the modelling
#'   framework - see `planning-docs/RASTER_FUSION_PERPIXEL_ENSEMBLE_DESIGN.md`.
#'
#'   `zonal_distribution_from_posterior()` is the cheap mukey-collapsed comparison arm (feeds
#'   `observed_data_by_mukey` in `generate_monte_carlo_realizations()`); `remarginalize_ensemble_to_posterior()`
#'   is the per-pixel transform; `remarginalized_awc()` is the first derived-quantity consumer
#'   (available water capacity, via `saxton_rawls_raster()` since SOLUS100 has no water-retention
#'   variable). Nothing here touches either existing pipeline - purely additive.
#' @name raster_fusion_bridge
NULL

#' Probabilities behind a `run_stage1_fusion()` posterior's percentile layers
#' @param posterior A `run_stage1_fusion()`-style `list(percentiles = <named P.. rasters>, ...)`.
#' @return Numeric probabilities in (0, 1), ascending, parsed from the `"P<pct>"` layer names.
#' @keywords internal
posterior_percentile_probs <- function(posterior) {
  nm <- names(posterior$percentiles)
  if (is.null(nm) || !all(grepl("^P[0-9]+$", nm))) {
    stop("posterior_percentile_probs(): posterior$percentiles must be a list named 'P1'/'P50'/... .")
  }
  p <- as.numeric(sub("^P", "", nm)) / 100
  if (is.unsorted(p)) stop("posterior_percentile_probs(): percentile layers are not in ascending order.")
  p
}

#' Zonal (per-mukey) Reduction of a Fused Posterior Raster - Benchmark Baseline
#'
#' Collapses `run_stage1_fusion()`'s per-pixel posterior to one distribution per mukey by taking
#' `terra::zonal()` means of its percentile layers, in the `low`/`rep`/`high` shape
#' `fuse_observed_data_into_priors()`'s `observed_data_by_mukey` accepts. **Percentile-based only -
#' there is deliberately no point-estimate ("mean") path**, which would narrow the fed-back prior
#' and understate uncertainty (see the design doc's Problem A redraft).
#'
#' Its sole purpose is the A.6 benchmark: the mukey-collapsed comparison arm for
#' `remarginalize_ensemble_to_posterior()`, run through the *same* tabular Monte Carlo machinery
#' via `observed_data_by_mukey`. Kept `@keywords internal` until that benchmark says whether the
#' cheap path is adequate for production (then promote) or not (then delete).
#'
#' @param posterior A `run_stage1_fusion()`-style `list(percentiles = , ...)`.
#' @param mukey_raster A categorical mukey `terra::SpatRaster` (any grid - resampled to the
#'   posterior grid via nearest-neighbour internally).
#' @param probs Length-3 ascending probabilities for the `low`/`rep`/`high` triplet; each must be
#'   one of the posterior's own percentile-layer probabilities (no interpolation).
#' @return Named list keyed by mukey code, each element `c(low =, rep =, high =)`.
#' @keywords internal
zonal_distribution_from_posterior <- function(posterior, mukey_raster,
                                              probs = c(0.05, 0.5, 0.95)) {
  if (length(probs) != 3 || is.unsorted(probs)) {
    stop("zonal_distribution_from_posterior(): `probs` must be 3 ascending probabilities.")
  }
  post_probs <- posterior_percentile_probs(posterior)
  if (!all(probs %in% post_probs)) {
    stop(sprintf("zonal_distribution_from_posterior(): every `probs` value must be an exact posterior percentile (have %s).",
                 paste(post_probs, collapse = ", ")))
  }
  layer_for <- function(p) posterior$percentiles[[paste0("P", round(p * 100))]]
  ref <- layer_for(probs[1])

  mk <- terra::resample(mukey_raster, ref, method = "near")
  mk_labels <- terra::cats(mk)[[1]]
  label_col <- names(mk_labels)[2]

  reduce_one <- function(p) {
    z <- terra::zonal(layer_for(p), mk, fun = "mean", na.rm = TRUE)
    stats::setNames(z[[2]], as.character(z[[1]]))
  }
  lo <- reduce_one(probs[1]); rp <- reduce_one(probs[2]); hi <- reduce_one(probs[3])

  mukeys <- Reduce(intersect, list(names(lo), names(rp), names(hi)))
  # zonal() keys by the factor label when the raster is categorical; map back to the mukey code.
  code_of <- function(key) {
    if (!is.null(mk_labels) && key %in% as.character(mk_labels[[1]])) {
      as.character(mk_labels[[label_col]][match(key, as.character(mk_labels[[1]]))])
    } else {
      key
    }
  }
  stats::setNames(
    lapply(mukeys, function(k) c(low = unname(lo[k]), rep = unname(rp[k]), high = unname(hi[k]))),
    vapply(mukeys, code_of, character(1))
  )
}

#' Invert a per-pixel posterior CDF at a per-mukey probability, vectorized over the grid
#'
#' For a scalar-per-mukey probability raster `u_ras`, returns the raster whose every cell is that
#' cell's posterior value at its `u`, by piecewise-linear interpolation on the posterior's
#' percentile knots with `rule = 2` constant extrapolation - the inverse-direction analogue of
#' `sim_linear_cdf_batch()`.
#' @param pct Named list of percentile-value `SpatRaster`s (ascending).
#' @param post_probs Matching ascending probabilities.
#' @param u_ras A `SpatRaster` of probabilities in (0, 1). **May have many layers** - each of the
#'   single-layer `pct` rasters is recycled across them, so one call inverts a whole realization
#'   stack in a fixed `~4 * (length(post_probs) - 1)` terra ops regardless of `nlyr(u_ras)`.
#' @return A `SpatRaster` of interpolated posterior values, `nlyr(u_ras)` layers.
#' @keywords internal
invert_posterior_cdf_raster <- function(pct, post_probs, u_ras) {
  k <- length(post_probs)
  # Work on a placeholder where u is NA (arithmetic on NA in terra doesn't reliably propagate to
  # the final ifel), then restore NA at the end via terra::mask().
  u <- terra::ifel(is.na(u_ras), 0.5, u_ras)
  out <- pct[[1]] + u * 0                           # rule = 2 low clamp, broadcast to nlyr(u_ras)
  for (i in seq_len(k - 1L)) {
    t <- (u - post_probs[i]) / (post_probs[i + 1L] - post_probs[i])
    seg_val <- pct[[i]] * (1 - t) + pct[[i + 1L]] * t
    in_seg <- (u >= post_probs[i]) & (u < post_probs[i + 1L])
    out <- terra::ifel(in_seg, seg_val, out)
  }
  out <- terra::ifel(u >= post_probs[k], pct[[k]], out)
  terra::mask(out, u_ras)                           # NA wherever the input u was NA
}

#' Saxton-Rawls water retention (field capacity, wilting point), raster / stack-native
#'
#' The vectorized, `terra`-native counterpart of `calculate_saxton_rawls_single()` - identical
#' equations and clamps, but every operation is `terra` `Arith`/`Math`/`clamp`/`ifel` so it runs
#' on multi-layer realization stacks in one pass. Used by `remarginalized_awc()`'s default
#' (`method = "saxton_rawls"`) path, because SOLUS100 publishes **no** water-retention variable -
#' `wr_3b`/`wr_15b` can never be fused directly, only derived from the fusable
#' sand/silt/clay/`dbovendry`/`soc`/`fragvol`.
#'
#' Because `remarginalized_awc()` re-marginalizes `sand`/`silt`/`clay` to their own fused
#' posteriors independently, a realization's texture triple may no longer sum to 100. This
#' reproduces `calculate_saxton_rawls_single()`'s renormalization: where `|sand + silt + clay - 100|
#' > 5`, the three are rescaled to sum to 100 before the equations run (the equations themselves use
#' only the sand and clay fractions).
#'
#' @param sand,clay,silt,db,rfv,om Percent (or g/cm^3 for `db`) `SpatRaster`s. Any may be single
#'   or multi-layer; single-layer inputs are recycled.
#' @return `list(fc =, wp =)` - volumetric % `SpatRaster`s, `wp < fc` enforced.
#' @keywords internal
saxton_rawls_raster <- function(sand, clay, silt, db, rfv, om) {
  sc <- terra::clamp(sand, 0, 100)
  cc <- terra::clamp(clay, 0, 100)
  sic <- terra::clamp(silt, 0, 100)
  tot <- sc + cc + sic
  renorm <- abs(tot - 100) > 5                       # matches calculate_saxton_rawls_single()
  s  <- terra::ifel(renorm, sc / tot, sc / 100)
  cl <- terra::ifel(renorm, cc / tot, cc / 100)
  dbc <- terra::clamp(db, 0.6, 2.5)
  rfvf <- terra::clamp(rfv, 0, 95) / 100
  omf <- terra::clamp(om, 0.1, 50) / 100

  theta_s <- -0.251 * s + 0.195 * cl + 0.011 * omf + 0.006 * s * omf -
    0.027 * cl * omf + 0.452 * s * cl + 0.299

  A <- exp(log(33) + 1.54 * s + 0.95 * cl + 0.025 * omf - 0.351 * s * omf -
             0.023 * cl * omf - 0.427 * s * cl + 0.015 * s^2 * cl^2)
  fc_g <- theta_s * (A / (A + 0.31))^0.17

  B <- exp(log(1500) + 0.02 * cl^2 + 0.14 * s - 0.0002 * s^2 * cl -
             0.002 * cl^2 * s - 0.0002 * cl^2 * omf + 0.0003 * cl^2 * s * omf)
  wp_g <- theta_s * (B / (B + 0.31))^0.17

  fc_v <- terra::clamp(fc_g * dbc * (1 - rfvf) * 100, 3, 65)
  wp_v <- terra::clamp(wp_g * dbc * (1 - rfvf) * 100, 1, 45)
  wp_v <- terra::ifel(wp_v >= fc_v, fc_v * 0.6, wp_v)
  list(fc = fc_v, wp = wp_v)
}

#' Per-Pixel Ensemble Re-Marginalization to a Fused Posterior
#'
#' Rank-preserving marginal transform: for each property and depth window, maps every retained
#' source realization through `F_prior` (the mukey ensemble's empirical CDF) then `F_posterior^-1`
#' (the pixel's fused posterior, piecewise-linear on its percentile knots). Each realization keeps
#' its rank position in every property and window, so the ensemble's empirical copula
#' (cross-property, cross-depth rank correlations) carries over; only the marginals become the
#' pixel's posterior.
#'
#' @section Documented approximations:
#' The copula / correlation structure stays mukey-level (no sub-mukey spatially-varying
#' dependence); within-depth-window vertical shape comes from the source realization (window
#' aggregates only - no horizon push-down here); the copula is assumed invariant under
#' re-marginalization and is the SSURGO tabular one (KSSL + joint-copula vertical correlation) -
#' SOLUS carries no joint information, so nothing is re-estimated; the preserved correlation is
#' **rank (Spearman)**, not Pearson; composition (comppct) weighting stays mukey-level (already in
#' the source ensemble's differential replicate counts). Still strictly better than zonal
#' aggregation. See `docs/09_multi_source_raster_fusion_pipeline.md` -> "Statistical structure of
#' the per-pixel bridge" for how this shapes a derived-quantity uncertainty band.
#'
#' @param mukey_ensemble An `extract_mukey_joint_ensemble()` result.
#' @param posterior_by_property_window Named list `[[property]][[window]]` of `run_stage1_fusion()`
#'   posteriors (each `list(percentiles = , ...)`), on a common grid. Only properties present in
#'   BOTH this and the ensemble, and windows present in both, are produced.
#' @param n_out Max source realizations to carry per pixel (default 250). If a mukey has more
#'   ensemble replicates, the first `n_out` are used (a fixed subsample - same rows across every
#'   property/window, so the copula is preserved).
#' @param summarize `TRUE` (default) returns per-pixel percentile rasters (`probs`); `FALSE`
#'   returns the raw `n_kept`-layer realization stacks (for downstream use, e.g. `remarginalized_awc()`).
#' @param probs Percentile probabilities for the `summarize = TRUE` output.
#' @return `list(percentiles = [[property]][[window]] = <named list of P.. SpatRasters>)` when
#'   `summarize`, else `list(ensemble = [[property]][[window]] = <n_kept-layer SpatRaster>)`. Plus
#'   `properties`, `window_names`, `n_out`, `n_kept`.
#' @seealso `extract_mukey_joint_ensemble()`, `zonal_distribution_from_posterior()`, `remarginalized_awc()`
#' @export
remarginalize_ensemble_to_posterior <- function(mukey_ensemble, posterior_by_property_window,
                                                 n_out = 250, summarize = TRUE,
                                                 probs = c(0.05, 0.25, 0.5, 0.75, 0.95)) {
  ens <- mukey_ensemble$by_mukey
  props <- intersect(mukey_ensemble$properties, names(posterior_by_property_window))
  if (length(props) == 0) {
    stop("remarginalize_ensemble_to_posterior(): no property is present in both the ensemble and posterior_by_property_window.")
  }

  # ONE retained realization count across every mukey AND window, so realization r is the same
  # source draw everywhere - preserves the copula across properties and across depth windows.
  n_kept <- min(
    n_out,
    min(vapply(ens, function(e) min(vapply(e$windows, nrow, integer(1))), integer(1)))
  )
  if (n_kept < 1L) stop("remarginalize_ensemble_to_posterior(): the ensemble has no retained realizations.")

  mukey_codes <- suppressWarnings(as.numeric(names(ens)))
  ok_code <- !is.na(mukey_codes)

  out_pct <- list()
  out_ens <- list()

  for (p in props) {
    windows <- intersect(mukey_ensemble$window_names, names(posterior_by_property_window[[p]]))
    for (w in windows) {
      posterior <- posterior_by_property_window[[p]][[w]]
      post_probs <- posterior_percentile_probs(posterior)
      pct <- posterior$percentiles

      # levels(mk) <- NULL makes the cell values the numeric mukey codes themselves, matching the
      # terra::subst(from = <numeric codes>) broadcast pattern used in R/raster-fusion.R.
      mk <- terra::resample(mukey_ensemble$mukey_raster, pct[[1]], method = "near")
      levels(mk) <- NULL

      # [n_mukey x n_kept] rank-probabilities (empirical F_prior), one row per mukey.
      # na.last = "keep": a mukey with a missing value for this property/window (e.g. some cokey's
      # texture sim failed) keeps that realization NA rather than letting base rank()'s default
      # na.last = TRUE rank the NAs 1..n and fabricate a probability for them.
      u_mat <- do.call(rbind, lapply(names(ens), function(mkc) {
        v <- ens[[mkc]]$windows[[w]][seq_len(n_kept), p]
        nf <- sum(is.finite(v))
        if (nf == 0L) return(rep(NA_real_, n_kept))
        (rank(v, ties.method = "average", na.last = "keep") - 0.5) / nf
      }))

      # ONE terra::subst -> an n_kept-layer per-cell u stack; ONE CDF inversion over the 9 knots.
      u_stack <- terra::subst(mk, from = mukey_codes[ok_code],
                              to = u_mat[ok_code, , drop = FALSE], others = NA_real_)
      # terra::subst() with an NA in `to` LEAVES the cell at its raw mukey code rather than NA-ing
      # it; every real rank-probability is strictly in (0, 1), so clamp anything else to NA (mukey
      # codes are >> 1) - otherwise invert_posterior_cdf_raster() rule-2-clamps it to the P99 layer.
      u_stack <- terra::ifel(u_stack > 0 & u_stack < 1, u_stack, NA_real_)
      inv <- invert_posterior_cdf_raster(pct, post_probs, u_stack)
      names(inv) <- paste0("r", seq_len(n_kept))

      if (summarize) {
        qs <- terra::quantile(inv, probs = probs, na.rm = TRUE)
        out_pct[[p]][[w]] <- stats::setNames(
          lapply(seq_along(probs), function(i) qs[[i]]),
          paste0("P", round(probs * 100))
        )
      } else {
        out_ens[[p]][[w]] <- inv
      }
    }
  }

  res <- if (summarize) list(percentiles = out_pct) else list(ensemble = out_ens)
  c(res, list(properties = props, window_names = mukey_ensemble$window_names,
              n_out = n_out, n_kept = n_kept))
}

#' Per-Pixel Available Water Capacity from a Re-Marginalized Ensemble
#'
#' Plant-available water capacity (cm, summed over the ensemble's depth windows), per realization,
#' then summarized per pixel. Deliberately does NOT use `calculate_aws_df()` (ROSETTA + van
#' Genuchten): that needs a live network POST and its own Monte Carlo per call, so it cannot run
#' per (pixel, realization).
#'
#' @section Method:
#' - **`"saxton_rawls"`** (default): `SOLUS100 publishes no water-retention variable`, so `wr_3b`/
#'   `wr_15b` can never be fused directly. Instead the fusable
#'   `sand_total`/`silt_total`/`clay_total`/`db` (plus optional `soc`, `rfv`) are re-marginalized
#'   to their per-pixel fused posteriors, then `saxton_rawls_raster()` derives field capacity /
#'   wilting point per (pixel, realization). `AWC = sum_windows (fc - wp)/100 * thickness` (the
#'   Saxton-Rawls conversion already applies the rock-fragment correction internally).
#' - **`"direct"`**: uses `wr_3b`/`wr_15b` posteriors supplied in `posterior_by_property_window`
#'   directly (`AWC = sum_windows (wr_3b - wr_15b)/100 * thickness * (1 - rfv/100)`). Only usable
#'   if the caller has a water-retention likelihood from some non-SOLUS source.
#'
#' @section Interpreting the per-pixel AWC distribution:
#' `AWC` is computed **per source realization** and combined by `terra::quantile()` per pixel, so
#' `awc_cm$P95 - awc_cm$P5` is a genuine 90% credible interval. Because
#' `remarginalize_ensemble_to_posterior()` threads one realization index through every property and
#' every window (rank-preserving), each realization's inputs are cross-property *and* vertically
#' consistent - a "wet" profile is wet at all depths. What the band does and does not capture:
#' \itemize{
#'   \item **Marginal** uncertainty of each input property: the full SSURGO x SOLUS fused posterior
#'     spread per pixel, propagated through Saxton-Rawls. Captured.
#'   \item **Dependence** (cross-property, cross-depth): the SSURGO tabular copula
#'     (KSSL + joint-copula vertical correlation), *held fixed*. SOLUS carries no joint information,
#'     so nothing is re-estimated. The preserved correlation is **rank (Spearman)**, not Pearson.
#'   \item **Pedotransfer-function error**: NOT propagated. Saxton-Rawls is applied deterministically
#'     per realization; the +/-15\% band `calculate_saxton_rawls_single()` returns is unused.
#'   \item **Sub-mukey variation in dependence**: none. Every pixel in a mukey shares the source
#'     ensemble; only the fused marginals vary pixel-to-pixel.
#' }
#' See `docs/09_multi_source_raster_fusion_pipeline.md` -> "Statistical structure of the per-pixel
#' bridge" for the full treatment.
#'
#' @param mukey_ensemble An `extract_mukey_joint_ensemble()` result.
#' @param posterior_by_property_window As in `remarginalize_ensemble_to_posterior()`, keyed by the
#'   ensemble's property (sim-column) names. `"saxton_rawls"` needs
#'   `sand_total`/`silt_total`/`clay_total`/`db` (+ optional `soc`, `rfv`); `"direct"` needs
#'   `wr_3b`/`wr_15b` (+ optional `rfv`).
#' @param method `"saxton_rawls"` (default) or `"direct"` - see Method.
#' @param n_out,probs As in `remarginalize_ensemble_to_posterior()`.
#' @param rock_fragment If `TRUE` (default) apply the rock-fragment correction when an `rfv`
#'   posterior is available.
#' @param soc_to_om Multiplier converting the `soc` posterior (organic carbon %) to organic matter
#'   % for Saxton-Rawls (default 1.724, the Van Bemmelen factor). Only used by `"saxton_rawls"`.
#' @param tile_rows Optional. Process the AOI in horizontal row-strips, releasing each strip's
#'   realization stack before the next. Per-pixel quantiles are cell-independent, so a tiled run is
#'   numerically identical to `tile_rows = NULL`. `NULL` (default) / `>= nrow` disables tiling.
#' @param restriction_depth Optional single-layer `terra::SpatRaster` of restriction depth in cm -
#'   typically `fetch_solus_restriction_depth(aoi_vect)` (MULTI_PROPERTY_FUSION_PLAN.md task S2),
#'   already right-censoring-guarded (cells with no detected restriction are `Inf`, never
#'   `SOLUS_RESTRICTION_CENSOR_CM` itself). `NULL` (default) preserves today's behavior exactly -
#'   every window contributes its full nominal thickness, regardless of bedrock. When supplied,
#'   resampled once (bilinear - a continuous depth value, not a category) onto the working grid,
#'   and each window's contribution is scaled by its **effective** thickness
#'   `pmax(0, pmin(bottom, restriction_depth) - top)` instead of the nominal `bottom - top`: a
#'   window entirely above the restriction keeps its full thickness, one straddling it gets partial
#'   credit, one entirely below it contributes exactly **0** (not `NA` - `NA` would incorrectly
#'   blank the whole pixel's AWC, including valid shallower windows, rather than just this window's
#'   contribution). A pixel with `restriction_depth` missing/`NA` (e.g. no coverage at the AOI edge)
#'   is treated the same as `Inf` - no truncation - the same "no evidence -> don't truncate"
#'   convention this pipeline already uses for other SSURGO/SOLUS coverage gaps, rather than
#'   propagating `NA` into the AWC total.
#' @return `list(awc_cm = <named list of P.. SpatRasters>, n_kept =, windows_used =, n_tiles =)`.
#'   AWC clamped at 0.
#' @seealso `remarginalize_ensemble_to_posterior()`, `saxton_rawls_raster()`, `calculate_aws_df()`,
#'   `fetch_solus_restriction_depth()`
#' @export
remarginalized_awc <- function(mukey_ensemble, posterior_by_property_window,
                               method = c("saxton_rawls", "direct"),
                               n_out = 250, probs = c(0.05, 0.25, 0.5, 0.75, 0.95),
                               rock_fragment = TRUE, soc_to_om = 1.724, tile_rows = NULL,
                               restriction_depth = NULL) {
  method <- match.arg(method)
  req <- if (method == "saxton_rawls") c("sand_total", "silt_total", "clay_total", "db")
         else c("wr_3b", "wr_15b")
  opt <- if (method == "saxton_rawls") c("soc", "rfv") else c("rfv")

  if (!all(req %in% mukey_ensemble$properties)) {
    stop(sprintf("remarginalized_awc(method = '%s'): the ensemble must carry %s.",
                 method, paste(req, collapse = "/")))
  }
  if (!all(req %in% names(posterior_by_property_window))) {
    stop(sprintf("remarginalized_awc(method = '%s'): posterior_by_property_window needs %s posteriors.",
                 method, paste(req, collapse = "/")))
  }
  keys <- intersect(c(req, opt), names(posterior_by_property_window))
  if (!rock_fragment) keys <- setdiff(keys, "rfv")
  pbpw <- posterior_by_property_window[keys]

  block <- if (method == "saxton_rawls") {
    function(me, p, rd) .awc_block(me, p, n_out, probs, method = "saxton_rawls", soc_to_om = soc_to_om,
                                   restriction_depth = rd)
  } else {
    function(me, p, rd) .awc_block(me, p, n_out, probs, method = "direct", restriction_depth = rd)
  }

  w1 <- intersect(mukey_ensemble$window_names, names(pbpw[[req[1]]]))
  if (length(w1) == 0) stop("remarginalized_awc(): no window has a posterior for the required properties.")
  ref <- pbpw[[req[1]]][[w1[1]]]$percentiles[[1]]

  # S2: resample restriction_depth onto the working grid ONCE here (not per window inside
  # .awc_block(), not per tile below) - a continuous depth value, bilinear, matching how every
  # other continuous SOLUS/SSURGO raster in this pipeline is resampled onto a working grid. NA
  # cells (e.g. no coverage at the AOI edge) are treated as Inf - "no evidence, don't truncate" -
  # rather than propagating NA into the AWC total.
  restriction_depth_aligned <- if (is.null(restriction_depth)) {
    NULL
  } else {
    rd <- terra::resample(restriction_depth, ref, method = "bilinear")
    terra::ifel(is.na(rd), Inf, rd)
  }

  if (is.null(tile_rows) || tile_rows >= terra::nrow(ref)) {
    return(c(block(mukey_ensemble, pbpw, restriction_depth_aligned), list(n_tiles = 1L)))
  }

  nr <- terra::nrow(ref)
  ymax0 <- terra::ymax(ref); yres <- terra::yres(ref)
  x0 <- terra::xmin(ref); x1 <- terra::xmax(ref)
  starts <- seq(1L, nr, by = as.integer(tile_rows))

  blocks <- lapply(starts, function(s) {
    e <- min(s + as.integer(tile_rows) - 1L, nr)
    ext_b <- terra::ext(x0, x1, ymax0 - e * yres, ymax0 - (s - 1L) * yres)
    me_b <- mukey_ensemble
    me_b$mukey_raster <- terra::crop(mukey_ensemble$mukey_raster, ext_b, snap = "out")
    pbpw_b <- lapply(pbpw, function(prop) lapply(prop, function(post) {
      post$percentiles <- lapply(post$percentiles, terra::crop, y = ext_b, snap = "near")
      post
    }))
    rd_b <- if (is.null(restriction_depth_aligned)) NULL else terra::crop(restriction_depth_aligned, ext_b, snap = "near")
    block(me_b, pbpw_b, rd_b)
  })

  merged <- stats::setNames(lapply(seq_along(probs), function(j) {
    do.call(terra::merge, lapply(blocks, function(b) b$awc_cm[[j]]))
  }), paste0("P", round(probs * 100)))

  list(awc_cm = merged, n_kept = blocks[[1]]$n_kept,
       windows_used = blocks[[1]]$windows_used, n_tiles = length(starts))
}

# One-block (whole-grid) AWC core; called directly, or once per row-strip when tile_rows is set.
# restriction_depth (S2): NULL (default) = today's unchanged behavior, every window uses its full
# nominal thickness `th`. Otherwise a single-layer SpatRaster (already resampled onto the working
# grid, already NA->Inf'd, by the caller) - each window's thickness becomes the per-pixel effective
# thickness pmax(0, pmin(bottom, restriction_depth) - top), via two terra::clamp() calls (clamp
# with only `upper` set is exactly pmin(x, upper); with only `lower` set, exactly pmax(x, lower)).
.awc_block <- function(mukey_ensemble, pbpw, n_out, probs, method, soc_to_om = 1.724,
                       restriction_depth = NULL) {
  rm_res <- remarginalize_ensemble_to_posterior(mukey_ensemble, pbpw, n_out = n_out, summarize = FALSE)
  e <- rm_res$ensemble

  awc <- NULL
  used <- character(0)
  for (i in seq_along(mukey_ensemble$window_names)) {
    w <- mukey_ensemble$window_names[i]
    win <- mukey_ensemble$depth_windows[[i]]
    top <- win[1]; bottom <- win[2]
    th <- if (is.null(restriction_depth)) {
      bottom - top
    } else {
      terra::clamp(terra::clamp(restriction_depth, upper = bottom, values = TRUE) - top,
                   lower = 0, values = TRUE)
    }

    if (method == "saxton_rawls") {
      if (any(vapply(c("sand_total", "silt_total", "clay_total", "db"),
                     function(p) is.null(e[[p]]) || is.null(e[[p]][[w]]), logical(1)))) next
      tmpl <- e[["sand_total"]][[w]][[1]]
      om  <- if (!is.null(e[["soc"]]) && !is.null(e[["soc"]][[w]])) e[["soc"]][[w]] * soc_to_om
             else terra::setValues(tmpl, 2)
      rfv <- if (!is.null(e[["rfv"]]) && !is.null(e[["rfv"]][[w]])) e[["rfv"]][[w]]
             else terra::setValues(tmpl, 0)
      sr <- saxton_rawls_raster(e[["sand_total"]][[w]], e[["clay_total"]][[w]],
                                e[["silt_total"]][[w]], e[["db"]][[w]], rfv, om)
      layer <- (sr$fc - sr$wp) / 100 * th
    } else {
      if (is.null(e[["wr_3b"]][[w]]) || is.null(e[["wr_15b"]][[w]])) next
      layer <- (e[["wr_3b"]][[w]] - e[["wr_15b"]][[w]]) / 100 * th
      if (!is.null(e[["rfv"]]) && !is.null(e[["rfv"]][[w]])) layer <- layer * (1 - e[["rfv"]][[w]] / 100)
    }

    awc <- if (is.null(awc)) layer else awc + layer
    used <- c(used, w)
  }
  if (is.null(awc)) stop("remarginalized_awc(): no usable window.")

  awc <- terra::clamp(awc, lower = 0, values = TRUE)
  qs <- terra::quantile(awc, probs = probs, na.rm = TRUE)
  list(
    awc_cm = stats::setNames(lapply(seq_along(probs), function(j) qs[[j]]), paste0("P", round(probs * 100))),
    n_kept = rm_res$n_kept,
    windows_used = used
  )
}
