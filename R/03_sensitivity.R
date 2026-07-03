# =============================================================================
# R/03_sensitivity.R
#
# Sensitivity analysis for the threshold age a*:
#   (a) one-at-a-time (OAT) variation of each intensity scale parameter
#       (a12, a21, a1d, a2d) by +/-25% and +/-50% around baseline;
#   (b) variation of the disabled/healthy excess-mortality ratio
#       r = a2d / a1d over the range [1.5, 4.0], with a1d held fixed.
#
# Requires config.R, R/01_model.R, and R/02_run_baseline.R to already be
# sourced (the latter for `baseline_results`, used as the OAT reference).
# =============================================================================

#' Threshold ages for all three initial states under a given parameter set
#'
#' @param intensities A list of intensity functions from
#'   [construct_intensities()].
#' @return A named numeric vector of a* values, one per entry of
#'   `INITIAL_STATES`.
threshold_ages_by_state <- function(intensities) {
  transitions <- build_transition_matrices(intensities, ages)
  purrr::map_dbl(
    INITIAL_STATES,
    ~ compute_lifespan_measures(transitions, .x, intensities, ages)$a_star
  )
}

# ---- (a) One-at-a-time parameter variation -----------------------------------
oat_parameters <- c("a12", "a21", "a1d", "a2d")
oat_labels <- c(
  a12 = "Disability incidence scale (a[12])",
  a21 = "Recovery scale (a[21])",
  a1d = "Healthy mortality scale (a[1+])",
  a2d = "Disabled mortality scale (a[2+])"
)
oat_multipliers <- c(0.50, 0.75, 1.00, 1.25, 1.50)
oat_mult_labels <- c("-50%", "-25%", "Baseline", "+25%", "+50%")

oat_grid <- tidyr::expand_grid(
  param = oat_parameters,
  mult_idx = seq_along(oat_multipliers)
)

sensitivity_oat <- purrr::pmap_dfr(oat_grid, function(param, mult_idx) {
  params <- BASELINE_PARAMS
  params[[param]] <- BASELINE_PARAMS[[param]] * oat_multipliers[mult_idx]
  a_star <- threshold_ages_by_state(do.call(construct_intensities, params))

  tibble::tibble(
    parameter   = oat_labels[[param]],
    param_short = param,
    multiplier  = oat_multipliers[mult_idx],
    mult_label  = oat_mult_labels[mult_idx],
    is_baseline = oat_multipliers[mult_idx] == 1.00,
    a_healthy   = a_star["healthy"],
    a_disabled  = a_star["disabled"],
    a_mixed     = a_star["mixed"]
  )
}) |>
  tidyr::pivot_longer(c(a_healthy, a_disabled, a_mixed),
                       names_to = "init_state", values_to = "a_star") |>
  dplyr::mutate(
    init_state = factor(
      init_state,
      levels = c("a_healthy", "a_disabled", "a_mixed"),
      labels = c("Healthy (1,0)", "Disabled (0,1)", "Mixed (0.75,0.25)")
    ),
    mult_label = factor(mult_label, levels = oat_mult_labels)
  )

# ---- (b) Disabled/healthy mortality ratio -------------------------------------
mortality_ratios <- c(1.5, 2.0, 2.5, 3.0, 3.5, 4.0)

sensitivity_ratio <- purrr::map_dfr(mortality_ratios, function(r) {
  params <- BASELINE_PARAMS
  params$a2d <- BASELINE_PARAMS$a1d * r   # a1d held fixed at its calibrated value
  a_star <- threshold_ages_by_state(do.call(construct_intensities, params))

  tibble::tibble(
    ratio      = r,
    a_healthy  = a_star["healthy"],
    a_disabled = a_star["disabled"],
    a_mixed    = a_star["mixed"]
  )
}) |>
  tidyr::pivot_longer(c(a_healthy, a_disabled, a_mixed),
                       names_to = "init_state", values_to = "a_star") |>
  dplyr::mutate(
    init_state = factor(
      init_state,
      levels = c("a_healthy", "a_disabled", "a_mixed"),
      labels = c("Healthy (1,0)", "Disabled (0,1)", "Mixed (0.75,0.25)")
    )
  )

save_table(
  sensitivity_oat, "TableS1", "sensitivity_oat",
  caption = "One-at-a-time sensitivity of threshold age a* to +/-25%/+/-50% parameter variation.",
  supplementary = TRUE
)
save_table(
  sensitivity_ratio, "TableS2", "sensitivity_ratio",
  caption = "Sensitivity of threshold age a* to the disabled/healthy mortality ratio (1.5-4.0).",
  supplementary = TRUE
)

# ---- Summary tables printed at the end of the pipeline (see run_analysis.R) --
sensitivity_oat_summary <- sensitivity_oat |>
  dplyr::group_by(parameter, init_state) |>
  dplyr::summarise(
    min_a = round(min(a_star), 1),
    max_a = round(max(a_star), 1),
    range = round(max(a_star) - min(a_star), 1),
    .groups = "drop"
  )

sensitivity_ratio_summary <- sensitivity_ratio |>
  dplyr::group_by(init_state) |>
  dplyr::summarise(range = round(max(a_star) - min(a_star), 1), .groups = "drop")

# ---- (c) Joint sensitivity: mortality ratio x disability incidence -----------
# Complements the OAT panel (Fig 4A) by varying two parameters simultaneously,
# so interaction effects — invisible under one-at-a-time variation — show up
# as curvature/non-parallel contours rather than being assumed away.
joint_ratios      <- seq(1.5, 4.0, by = 0.5)
joint_multipliers <- seq(0.5, 1.5, by = 0.25)

sensitivity_joint <- tidyr::expand_grid(ratio = joint_ratios, mult = joint_multipliers) |>
  purrr::pmap_dfr(function(ratio, mult) {
    params <- BASELINE_PARAMS
    params$a2d <- BASELINE_PARAMS$a1d * ratio
    params$a12 <- BASELINE_PARAMS$a12 * mult
    a_star <- threshold_ages_by_state(do.call(construct_intensities, params))
    tibble::tibble(ratio = ratio, mult = mult, a_healthy = a_star["healthy"],
                    a_disabled = a_star["disabled"], a_mixed = a_star["mixed"])
  }) |>
  tidyr::pivot_longer(c(a_healthy, a_disabled, a_mixed),
                       names_to = "init_state", values_to = "a_star") |>
  dplyr::mutate(init_state = factor(
    init_state, levels = c("a_healthy", "a_disabled", "a_mixed"),
    labels = c("Healthy (1,0)", "Disabled (0,1)", "Mixed (0.75,0.25)")
  ))

save_table(
  sensitivity_joint, "TableS5", "sensitivity_joint",
  caption = "Joint sensitivity of a* to simultaneous variation in the mortality ratio and disability-incidence scale.",
  supplementary = TRUE
)
