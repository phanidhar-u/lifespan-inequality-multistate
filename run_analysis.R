# =============================================================================
# run_analysis.R
#
# Reproduces every figure and table in:
#   Upadhyaya & Kumar, "Threshold Age of Lifespan Inequality in a Multistate
#   (Healthy-Disabled) Framework"
#
# Usage:
#   Rscript run_analysis.R
#   (or source() interactively from the repository root)
#
# Outputs are written to output/figures/ and output/tables/.
# =============================================================================

source("config.R")
source("R/00_io.R")
source("R/01_model.R")
source("R/02_run_baseline.R")
source("R/03_sensitivity.R")
source("R/06_robustness.R")
source("R/04_figures.R")
source("R/05_calibration_check.R")

# ---- Build and save figures ---------------------------------------------------
# To add a new figure later (main text or supplementary), write a plotting
# function in R/04_figures.R, call it here, then save_figure() it — no other
# file needs to change. Supplementary figures just pass supplementary = TRUE
# and a "FigS#" label, which also routes them to output/figures/supplementary/.

fig1 <- plot_intensities(baseline_intensities, ages)
save_figure(fig1, "Fig1", "transition_intensities",
            caption = "Age-varying transition intensities (model parameters).",
            width = 9, height = 5)

fig2 <- plot_occupancy(baseline_transitions, ages)
save_figure(fig2, "Fig2", "state_occupancy",
            caption = "State occupancy and survival by initial health composition.",
            width = 14, height = 4.5)

fig3 <- plot_threshold(baseline_results, ages)
save_figure(fig3, "Fig3", "threshold_condition",
            caption = "Threshold-age condition (eq. 30) by initial health composition.",
            width = 7, height = 10)

fig4 <- plot_sensitivity(sensitivity_oat, sensitivity_ratio, baseline_results)
save_figure(fig4, "Fig4", "sensitivity_analysis",
            caption = "Sensitivity of threshold age a* to parameter variation.",
            width = 10, height = 11)

figS4 <- plot_joint_sensitivity(sensitivity_joint)
save_figure(figS4, "FigS4", "joint_sensitivity_heatmap",
            caption = "Joint sensitivity of a* to the mortality ratio and disability-incidence scale.",
            width = 11, height = 4, supplementary = TRUE)

# ---- Console summary ----------------------------------------------------------
cat("\n", strrep("=", 65), "\n")
cat("TABLE 1: MULTISTATE DEMOGRAPHIC QUANTITIES\n")
cat(strrep("=", 65), "\n")
print(as.data.frame(table1))

cat("\n", strrep("=", 65), "\n")
cat("SENSITIVITY SUMMARY: range of a* over +/-50% parameter variation\n")
cat(strrep("=", 65), "\n")
print(as.data.frame(sensitivity_oat_summary))

cat("\nMortality-ratio sensitivity (range across ratio 1.5 to 4.0):\n")
print(as.data.frame(sensitivity_ratio_summary))

cat("\nAll figures written to '", DIR_FIGURES, "/{main,supplementary}/'\n", sep = "")
cat("All tables written to '",  DIR_TABLES,  "/{main,supplementary}/'\n", sep = "")

# ---- Manifest + session info for the reproducibility appendix ------------------
write_manifest()
writeLines(capture.output(sessionInfo()), "output/session_info.txt")
cat("Manifest written to 'output/manifest.md'\n")
