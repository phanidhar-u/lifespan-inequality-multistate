# =============================================================================
# R/07_robustness.R
#
# Numerical robustness and calibration checks:
#   (a) Grid-resolution convergence (step size dx from 1.0 down to 0.05 years)
#   (b) Starting-age sensitivity (a0 from 45 to 60 years)
#   (c) Calibration verification against CDC/NCHS 2019 life table targets
# =============================================================================

# ---- (a) Grid-resolution convergence -----------------------------------------
grid_steps <- c(1.0, 0.5, 0.25, 0.1, 0.05)

robustness_grid <- purrr::map_dfr(grid_steps, function(step) {
  ages_r <- seq(AGE_MIN, AGE_MAX, by = step)
  transitions_r <- build_transition_matrices(baseline_intensities, ages_r)
  a_star <- purrr::map_dbl(
    INITIAL_STATES,
    ~ compute_lifespan_measures(transitions_r, .x, baseline_intensities, ages_r)$a_star
  )
  tibble::tibble(
    dx         = step,
    a_healthy  = a_star["healthy"],
    a_disabled = a_star["disabled"],
    a_mixed    = a_star["mixed"]
  )
})

save_table(
  robustness_grid, "TableS3_grid_resolution_convergence",
  caption = "Convergence of threshold age a* as the age-grid step dx is refined.",
  supplementary = TRUE
)

# ---- (b) Starting-age sensitivity --------------------------------------------
start_ages <- c(45, 50, 55, 60)

robustness_start_age <- purrr::map_dfr(start_ages, function(a0) {
  ages_r <- seq(a0, AGE_MAX, by = AGE_STEP)
  transitions_r <- build_transition_matrices(baseline_intensities, ages_r)
  a_star <- purrr::map_dbl(
    INITIAL_STATES,
    ~ compute_lifespan_measures(transitions_r, .x, baseline_intensities, ages_r)$a_star
  )
  tibble::tibble(
    start_age  = a0,
    a_healthy  = a_star["healthy"],
    a_disabled = a_star["disabled"],
    a_mixed    = a_star["mixed"]
  )
})

save_table(
  robustness_start_age, "TableS4_starting_age_sensitivity",
  caption = "Sensitivity of threshold age a* to the model's origin (starting) age.",
  supplementary = TRUE
)

# ---- (c) Calibration check against US 2019 targets ---------------------------
calibration_log <- c(
  strrep("=", 65),
  "CALIBRATION VERIFICATION (CDC/NCHS 2019 Period Life Table Targets)",
  strrep("=", 65),
  sprintf("%-30s  %.2f yrs  (target: %.1f)",
          "e0ms, healthy-start",  baseline_results$healthy$e0,  CALIBRATION_TARGETS$e0_healthy_start),
  sprintf("%-30s  %.2f yrs",
          "e0ms, disabled-start", baseline_results$disabled$e0),
  sprintf("%-30s  %.2f yrs",
          "e0ms, mixed-start",    baseline_results$mixed$e0),
  sprintf("%-30s  %.5f  (target: %.5f)",
          "mu1d at age 50", baseline_intensities$mu1d(50), CALIBRATION_TARGETS$mu1d_at_50),
  sprintf("%-30s  %.5f  (target: %.5f)",
          "mu1d at age 70", baseline_intensities$mu1d(70), CALIBRATION_TARGETS$mu1d_at_70),
  sprintf("%-30s  %.5f  (target: %.5f)",
          "mu1d at age 80", baseline_intensities$mu1d(80), CALIBRATION_TARGETS$mu1d_at_80)
)

writeLines(calibration_log, file.path(DIR_TAB_SUPP, "calibration_check.txt"))
