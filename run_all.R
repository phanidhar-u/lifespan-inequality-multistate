# =============================================================================
# run_all.R
#
# Master reproduction pipeline for:
#   Upadhyaya, Aburto, Kumar, & Singh, "Threshold Age for Multistate Life Table
#   Entropy: A Perturbation Approach"
#
# Usage:
#   Rscript run_all.R
#   (or source("run_all.R") from the replication_code/ root directory)
#
# Generates all theoretical figures, numerical illustration figures,
# multi-indicator comparison figures, and tables.
# =============================================================================

cat("\n[1/7] Initializing configuration and global settings...\n")
source("config.R")
source("R/00_io.R")
source("R/01_multistate_core.R")

cat("[2/7] Running baseline calibrated model and Table 1...\n")
source("R/02_run_baseline.R")

cat("[3/7] Running sensitivity analysis (OAT, mortality ratio, joint grid)...\n")
source("R/03_sensitivity.R")
source("R/07_robustness.R")

cat("[4/7] Generating theoretical mechanics figures (Section 3)...\n")
source("R/04_fig_theory.R")

# Figure 1: Mechanics of localized killing perturbation (Lemmas 2 & 3)
fig1 <- generate_fig1_perturbation_mechanics(baseline_intensities, ages)
save_figure(
  fig1, "Fig1_perturbation_mechanics",
  alt_names = c("fig_perturbation_mechanics"),
  caption   = "Theoretical mechanics of age-local killing perturbation at age a = 75 (Lemmas 2 and 3).",
  width     = 9.5, height = 4.2
)

# Figure 2: Three-channel geometric decomposition of life disparity sensitivity (Prop. 4)
fig2 <- generate_fig2_prop4_structural_anatomy(baseline_intensities, ages)
save_figure(
  fig2, "Fig2_prop4_structural_anatomy",
  alt_names = c("fig_prop4_structural_anatomy"),
  caption   = "Three-channel geometric decomposition of multistate life disparity sensitivity at age a = 75 (Proposition 4).",
  width     = 6.8, height = 4.3
)

cat("[5/7] Generating numerical illustration figures (Section 4)...\n")
source("R/05_fig_numerical.R")

# Figure 3: Transition intensities
fig3 <- plot_fig3_intensities(baseline_intensities, ages)
save_figure(
  fig3, "Fig3_transition_intensities",
  alt_names = c("fig3_intensities", "fig1_transition_intensities"),
  caption   = "Gompertz-type transition intensity schedules on logarithmic scale.",
  width     = 8.5, height = 4.8
)

# Figure 4: State occupancy probabilities and aggregate survival S^ms(x)
fig4 <- plot_fig4_occupancy(baseline_transitions, ages)
save_figure(
  fig4, "Fig4_state_occupancy",
  alt_names = c("fig_state_occupancy_all_cases", "fig2_state_occupancy"),
  caption   = "State-occupancy probabilities and aggregate multistate survival across initial health compositions.",
  width     = 13.0, height = 4.2
)

# Figure 5: Verification of threshold condition across initial distributions
fig5 <- plot_fig5_threshold(baseline_results, ages)
save_figure(
  fig5, "Fig5_threshold_condition",
  alt_names = c("fig_mult_threshold_multistate", "fig3_threshold_condition"),
  caption   = "Verification of the threshold condition (eq. 30) for three initial state distributions.",
  width     = 7.5, height = 9.8
)

# Figure 6: Sensitivity analysis (OAT and excess mortality ratio)
fig6 <- plot_fig6_sensitivity(sensitivity_oat, sensitivity_ratio, baseline_results)
save_figure(
  fig6, "Fig6_sensitivity_analysis",
  alt_names = c("fig4_sensitivity_analysis"),
  caption   = "Sensitivity of threshold age a* to one-at-a-time scaling and excess mortality ratio.",
  width     = 10.0, height = 10.5
)

# Supplementary Figure S4: Joint sensitivity heatmap
figS4 <- plot_figS4_joint_sensitivity(sensitivity_joint)
save_figure(
  figS4, "FigS4_joint_sensitivity_heatmap",
  caption       = "Joint sensitivity heatmap of a* across mortality ratio and disability incidence.",
  width         = 11.0, height = 4.0, supplementary = TRUE
)

cat("[6/7] Generating multi-indicator threshold comparison (Appendix D)...\n")
source("R/06_fig_indicators.R")
ind_res <- generate_indicator_comparison(baseline_transitions, baseline_intensities, ages)

save_figure(
  ind_res$plot, "FigD1_threshold_comparison",
  alt_names = c("fig_threshold_comparison_compact"),
  caption   = "Distinct threshold ages for five lifespan inequality indicators across initial compositions (Appendix D).",
  width     = 8.0, height = 5.0
)

cat("[7/7] Writing reproduction manifest and summary output...\n")
write_manifest()
writeLines(capture.output(sessionInfo()), file.path(DIR_OUT, "session_info.txt"))

# Summary tables to console
cat("\n", strrep("=", 72), "\n")
cat("TABLE 1: MULTISTATE DEMOGRAPHIC QUANTITIES AND THRESHOLD AGES\n")
cat(strrep("=", 72), "\n")
print(as.data.frame(table1_df))

cat("\n", strrep("=", 72), "\n")
cat("TABLE D.1: THRESHOLD AGES ACROSS 5 LIFESPAN INEQUALITY INDICATORS\n")
cat(strrep("=", 72), "\n")
print(as.data.frame(ind_res$table))

cat("\nAll outputs successfully generated into 'output/figures/' and 'output/tables/'!\n")
cat("Manifest generated at 'output/manifest.md'.\n")
