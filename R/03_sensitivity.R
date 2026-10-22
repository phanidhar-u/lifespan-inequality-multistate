# =============================================================================
# R/03_sensitivity.R
#
# Sensitivity analysis for threshold age a*:
#   (a) One-at-a-time (OAT) variation of intensity level parameters
#       (a12, a21, a1d, a2d) by +/-25% and +/-50% around baseline;
#   (b) Variation of disabled/healthy excess-mortality ratio r = a2d / a1d in [1.5, 4.0];
#   (c) Joint sensitivity over mortality ratio x disability incidence.
# =============================================================================

threshold_ages_by_state <- function(intensities) {
  transitions <- build_transition_matrices(intensities, ages)
  purrr::map_dbl(
    INITIAL_STATES,
    ~ compute_lifespan_measures(transitions, .x, intensities, ages)$a_star
  )
}

# ---- (a) One-at-a-time (OAT) parameter variation -----------------------------
oat_parameters <- c("a12", "a21", "a1d", "a2d")
oat_labels <- c(
  a12 = "Disability incidence scale (a[12])",
  a21 = "Recovery scale (a[21])",
  a1d = "Healthy mortality scale (a[1*'\u2020'])",
  a2d = "Disabled mortality scale (a[2*'\u2020'])"
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

# ---- (b) Disabled/healthy mortality ratio r ----------------------------------
mortality_ratios <- c(1.5, 2.0, 2.5, 3.0, 3.5, 4.0)

sensitivity_ratio <- purrr::map_dfr(mortality_ratios, function(r) {
  params <- BASELINE_PARAMS
  params$a2d <- BASELINE_PARAMS$a1d * r   # a1d held fixed at calibrated level
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

# ---- (c) Joint sensitivity grid: ratio x incidence ---------------------------
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

# Save sensitivity tables
save_table(
  sensitivity_oat, "TableS1_sensitivity_oat",
  caption = "One-at-a-time sensitivity of threshold age a* to +/-25%/+/-50% parameter variation.",
  supplementary = TRUE
)

save_table(
  sensitivity_ratio, "TableS2_sensitivity_ratio",
  caption = "Sensitivity of threshold age a* to excess mortality ratio (1.5 - 4.0).",
  supplementary = TRUE
)

save_table(
  sensitivity_joint, "TableS5_sensitivity_joint",
  caption = "Joint sensitivity of a* to simultaneous variation in mortality ratio and disability-incidence scale.",
  supplementary = TRUE
)
