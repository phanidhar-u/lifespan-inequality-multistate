# =============================================================================
# R/02_run_baseline.R
#
# Runs the calibrated baseline model for the three initial health
# compositions (Healthy, Disabled, Mixed) and constructs Table 1.
# =============================================================================

baseline_intensities <- do.call(construct_intensities, BASELINE_PARAMS)
baseline_transitions <- build_transition_matrices(baseline_intensities, ages)

baseline_results <- purrr::map(
  INITIAL_STATES,
  ~ compute_lifespan_measures(baseline_transitions, .x, baseline_intensities, ages)
)

table1_df <- tibble::tibble(
  `Initial state`                  = c("Healthy (1,0)", "Disabled (0,1)", "Mixed (0.75,0.25)"),
  `e0ms (yrs)`                     = round(purrr::map_dbl(baseline_results, "e0"), 1),
  `edagger (yrs)`                  = round(purrr::map_dbl(baseline_results, "edagger"), 1),
  `Hms`                            = round(purrr::map_dbl(baseline_results, "H"), 3),
  `a* (yrs)`                       = round(purrr::map_dbl(baseline_results, "a_star"), 1),
  `Hbar(a*)`                       = round(purrr::map_dbl(baseline_results, "Hbar_at_astar"), 3),
  `1 + Hms + ln S(a*)`             = round(purrr::map_dbl(baseline_results, "rhs_at_astar"), 3)
)

save_table(
  table1_df, "Table1_baseline_results",
  caption = "Multistate demographic quantities and threshold ages by initial state distribution (Table 1 in manuscript)."
)
