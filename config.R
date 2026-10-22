# =============================================================================
# config.R
#
# Global configuration for the multistate threshold-age analysis.
# Upadhyaya, Aburto, Kumar, & Singh, "Threshold Age for Multistate Life Table
# Entropy: A Perturbation Approach"
#
# Isolating the age grid, calibrated baseline parameters, and paths here
# ensures full auditability without hardcoded magic numbers in analysis scripts.
# =============================================================================

suppressPackageStartupMessages({
  library(Matrix)
  library(expm)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(purrr)
  library(readr)
  library(ggplot2)
  library(patchwork)
  library(scales)
})

# ---- Age grid ---------------------------------------------------------------
AGE_MIN  <- 50
AGE_MAX  <- 110
AGE_STEP <- 0.25
ages     <- seq(AGE_MIN, AGE_MAX, by = AGE_STEP)
N_AGES   <- length(ages)

# ---- Calibrated baseline intensity parameters -------------------------------
# Disability transition parameters (a12, b12, a21, b21) follow the Health and
# Retirement Study (HRS) literature (Lievre et al. 2003; Crimmins et al. 2016).
#
# Mortality parameters (a1d, b1d, a2d, b2d) are calibrated so that healthy-start
# life expectancy at age 50 approximates the 2019 US period life table value
# (Arias & Xu 2022, CDC/NCHS). Excess mortality ratio a2d / a1d = 2.5 at baseline.
BASELINE_PARAMS <- list(
  a12 = 0.0100, b12 = 0.0450,   # healthy -> disabled (incidence)
  a21 = 0.0600, b21 = 0.0200,   # disabled -> healthy (recovery)
  a1d = 0.0030, b1d = 0.0800,   # healthy -> dead
  a2d = 0.0075, b2d = 0.0800    # disabled -> dead
)

# ---- Calibration reference targets (CDC 2019) -------------------------------
CALIBRATION_TARGETS <- list(
  e0_healthy_start = 32.3,   # years, US CDC 2019 period life table
  mu1d_at_50       = 0.00416,
  mu1d_at_70       = 0.01965,
  mu1d_at_80       = 0.04720
)

# ---- Initial health compositions examined throughout manuscript -------------
INITIAL_STATES <- list(
  healthy  = c(1.00, 0.00),
  disabled = c(0.00, 1.00),
  mixed    = c(0.75, 0.25)
)

INITIAL_STATE_LABELS <- c(
  "healthy"  = "Healthy (1,0)",
  "disabled" = "Disabled (0,1)",
  "mixed"    = "Mixed (0.75,0.25)"
)

# ---- Colorblind-safe palettes matching paper figures -------------------------
COLOR_PALETTE_3 <- c(
  "Healthy (1,0)"      = "#1b6ca8",
  "Disabled (0,1)"     = "#c0392b",
  "Mixed (0.75,0.25)"  = "#27ae60",
  "Mixed (0.75, 0.25)" = "#27ae60"
)

COLOR_PALETTE_INTENSITIES <- c(
  "#0072B2", "#009E73", "#D55E00", "#CC79A7"
)

# ---- Output directories ------------------------------------------------------
DIR_OUT         <- "output"
DIR_FIG_MAIN    <- file.path(DIR_OUT, "figures", "main")
DIR_FIG_SUPP    <- file.path(DIR_OUT, "figures", "supplementary")
DIR_TAB_MAIN    <- file.path(DIR_OUT, "tables", "main")
DIR_TAB_SUPP    <- file.path(DIR_OUT, "tables", "supplementary")

for (d in c(DIR_FIG_MAIN, DIR_FIG_SUPP, DIR_TAB_MAIN, DIR_TAB_SUPP)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}
