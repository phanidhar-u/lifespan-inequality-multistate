# =============================================================================
# R/06_robustness.R
#
# Robustness checks beyond the OAT/ratio sensitivity in R/03_sensitivity.R:
#   (a) Numerical convergence of a* to the age-grid step size dx
#   (b) Sensitivity of a* to the model's starting (left-truncation) age
#
# Requires config.R, R/00_io.R, and R/01_model.R to already be sourced.
# =============================================================================

# ---- (a) Grid-resolution convergence ------------------------------------------
grid_steps <- c(1, 0.5, 0.25, 0.1, 0.05)   # AGE_STEP = 0.25 is the baseline

robustness_grid <- purrr::map_dfr(grid_steps, function(step) {
  ages_r <- seq(AGE_MIN, AGE_MAX, by = step)
  transitions_r <- build_transition_matrices(baseline_intensities, ages_r)
  a_star <- purrr::map_dbl(
    INITIAL_STATES,
    ~ compute_lifespan_measures(transitions_r, .x, baseline_intensities, ages_r)$a_star
  )
  tibble::tibble(dx = step, a_healthy = a_star["healthy"],
                  a_disabled = a_star["disabled"], a_mixed = a_star["mixed"])
})

save_table(
  robustness_grid, "TableS3", "grid_resolution_convergence",
  caption = "Convergence of threshold age a* as the age-grid step dx is refined from 1 to 0.05 years.",
  supplementary = TRUE
)

# ---- (b) Starting-age (left-truncation) sensitivity ---------------------------
# Intensities are parameterised relative to age 50 (see construct_intensities()
# in R/01_model.R), so no re-parameterisation is needed here — only the grid
# on which the model is evaluated changes.
start_ages <- c(45, 50, 55, 60)

robustness_start_age <- purrr::map_dfr(start_ages, function(a0) {
  ages_r <- seq(a0, AGE_MAX, by = AGE_STEP)
  transitions_r <- build_transition_matrices(baseline_intensities, ages_r)
  a_star <- purrr::map_dbl(
    INITIAL_STATES,
    ~ compute_lifespan_measures(transitions_r, .x, baseline_intensities, ages_r)$a_star
  )
  tibble::tibble(start_age = a0, a_healthy = a_star["healthy"],
                  a_disabled = a_star["disabled"], a_mixed = a_star["mixed"])
})

save_table(
  robustness_start_age, "TableS4", "starting_age_sensitivity",
  caption = "Sensitivity of threshold age a* to the model's starting (left-truncation) age.",
  supplementary = TRUE
)
