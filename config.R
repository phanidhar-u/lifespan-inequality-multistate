# =============================================================================
# config.R
#
# Global configuration for the multistate threshold-age analysis.
#
# Isolating the age grid and calibrated parameter set here keeps R/*.R free
# of magic numbers and makes the calibration auditable in one place.
#
# Upadhyaya & Kumar, "Threshold Age of Lifespan Inequality in a Multistate
# (Healthy-Disabled) Framework"
# =============================================================================

suppressPackageStartupMessages({
  library(expm)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(purrr)
  library(readr)
  library(ggplot2)
  library(patchwork)
})

# ---- Age grid ---------------------------------------------------------------
AGE_MIN <- 50
AGE_MAX <- 110
AGE_STEP <- 0.25
ages <- seq(AGE_MIN, AGE_MAX, by = AGE_STEP)

# ---- Calibrated baseline intensity parameters --------------------------------
# Disability-transition parameters (a12, b12, a21, b21) follow the Health and
# Retirement Study literature (Lievre, Brouard & Heathcote 2003; Crimmins,
# Zhang & Saito 2016) and are held fixed throughout.
#
# Mortality parameters (a1d, b1d, a2d, b2d) were recalibrated so that the
# healthy-start life expectancy at age 50, e_0^ms(50 | pi_0 = (1,0)),
# approximates the 2019 US period life table value of 32.3 years
# (Arias & Xu 2022, CDC/NCHS National Vital Statistics Reports). The
# disabled/healthy excess-mortality ratio a2d / a1d = 2.5 at baseline is
# consistent with Rogers, Rogers & Belanger (1990), Maddox et al. (1994),
# and Cai, Hayward & Saito (2006).
BASELINE_PARAMS <- list(
  a12 = 0.0100, b12 = 0.0450,   # healthy -> disabled (incidence)
  a21 = 0.0600, b21 = 0.0200,   # disabled -> healthy (recovery)
  a1d = 0.0030, b1d = 0.0800,   # healthy -> dead
  a2d = 0.0075, b2d = 0.0800    # disabled -> dead
)

# Calibration targets reported for verification in 05_calibration_check.R
CALIBRATION_TARGETS <- list(
  e0_healthy_start = 32.3,   # yrs, US CDC 2019
  mu1d_at_50 = 0.00416,
  mu1d_at_70 = 0.01965,
  mu1d_at_80 = 0.04720
)

# ---- Initial health compositions examined throughout the paper --------------
INITIAL_STATES <- list(
  healthy  = c(1.00, 0.00),
  disabled = c(0.00, 1.00),
  mixed    = c(0.75, 0.25)
)

# ---- Output paths -------------------------------------------------------------
DIR_FIGURES <- "output/figures"
DIR_TABLES  <- "output/tables"
dir.create(DIR_FIGURES, recursive = TRUE, showWarnings = FALSE)
dir.create(DIR_TABLES,  recursive = TRUE, showWarnings = FALSE)
