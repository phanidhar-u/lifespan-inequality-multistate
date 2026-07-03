# =============================================================================
# R/05_calibration_check.R
#
# Verifies that the calibrated baseline model reproduces the target
# quantities reported in config.R (US CDC/NCHS 2019 period life table,
# Arias & Xu 2022). Writes a plain-text log for the reproducibility
# appendix and prints the same summary to the console.
#
# Requires config.R, R/01_model.R, and R/02_run_baseline.R to already be
# sourced.
# =============================================================================

calibration_log <- c(
  strrep("=", 65),
  "CALIBRATION VERIFICATION",
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

writeLines(calibration_log)
writeLines(calibration_log, file.path(DIR_TABLES, "calibration_check.txt"))
