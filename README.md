# Replication Code: Threshold Age of Lifespan Inequality in a Multistate Framework

Reproducible code accompanying:

> **Phanidhar Upadhyaya & Kaushalendra Kumar**, *"Threshold Age of Lifespan Inequality in a Multistate Framework."*

This repository implements the continuous-age multistate (Healthy / Disabled / Dead) framework with age-varying Gompertz transition intensities, evaluated via product integration over a discretized age grid ($\Delta x = 0.25$ years). It reproduces every theoretical figure, numerical illustration figure, multi-indicator comparison, and table in the manuscript.

---

## Quick Start

### 1. Requirements

The code requires **R (>= 4.0.0)** and the following packages:

```r
install.packages(c(
  "Matrix", "expm", "dplyr", "tidyr", "tibble",
  "purrr", "readr", "ggplot2", "patchwork", "scales"
))
```

### 2. Execute Full Pipeline

Run the master script from the repository root:

```bash
Rscript run_all.R
```

Or within an interactive R session:

```r
source("run_all.R")
```

This single command generates all figures (both publication-ready vector **PDF** and **300 DPI PNG** preview) and tables (**CSV**) in `output/`, along with a complete `output/manifest.md`.

---

## Mapping: Manuscript Items to Output Files

| Manuscript Item | Description | Code File | Output Filename(s) |
| :--- | :--- | :--- | :--- |
| **Figure 1** | Local perturbation mechanics at $a = 75$ (Lemmas 2 & 3) | `R/04_fig_theory.R` | `output/figures/main/fig_perturbation_mechanics.{pdf,png}`<br>`output/figures/main/Fig1_perturbation_mechanics.{pdf,png}` |
| **Figure 2** | Three-channel geometric decomposition of $e^\dagger_{\mathrm{ms}}$ sensitivity (Proposition 4) | `R/04_fig_theory.R` | `output/figures/main/fig_prop4_structural_anatomy.{pdf,png}`<br>`output/figures/main/Fig2_prop4_structural_anatomy.{pdf,png}` |
| **Figure 3** | Transition intensity schedules ($\mu_{12}, \mu_{21}, \mu_{1\dagger}, \mu_{2\dagger}$) on log scale | `R/05_fig_numerical.R` | `output/figures/main/fig3_intensities.{pdf,png}`<br>`output/figures/main/Fig3_transition_intensities.{pdf,png}` |
| **Figure 4** | State-occupancy probabilities and aggregate survival $S^{\mathrm{ms}}(x)$ | `R/05_fig_numerical.R` | `output/figures/main/fig_state_occupancy_all_cases.{pdf,png}`<br>`output/figures/main/Fig4_state_occupancy.{pdf,png}` |
| **Figure 5** | Verification of threshold condition $\bar{H}^{\mathrm{ms}}(a) = 1 + H_{\mathrm{ms}} + \ln S^{\mathrm{ms}}(a)$ | `R/05_fig_numerical.R` | `output/figures/main/fig_mult_threshold_multistate.{pdf,png}`<br>`output/figures/main/Fig5_threshold_condition.{pdf,png}` |
| **Figure 6** | Sensitivity of threshold age $a^*$ to parameter variation (OAT and mortality ratio) | `R/05_fig_numerical.R` | `output/figures/main/fig4_sensitivity_analysis.{pdf,png}`<br>`output/figures/main/Fig6_sensitivity_analysis.{pdf,png}` |
| **Figure D.1** | Distinct threshold ages for 5 lifespan inequality indicators across initial compositions | `R/06_fig_indicators.R` | `output/figures/main/fig_threshold_comparison_compact.{pdf,png}`<br>`output/figures/main/FigD1_threshold_comparison.{pdf,png}` |
| **Figure S4** | Joint sensitivity heatmap ($a^*$ over mortality ratio $\times$ disability incidence) | `R/05_fig_numerical.R` | `output/figures/supplementary/FigS4_joint_sensitivity_heatmap.{pdf,png}` |
| **Table 1** | Baseline demographic quantities ($e_0, e^\dagger, H, a^*, \text{LHS}, \text{RHS}$) | `R/02_run_baseline.R` | `output/tables/main/Table1_baseline_results.csv` |
| **Table D.1** | Multistate threshold conditions & numerical roots for 5 inequality indicators | `R/06_fig_indicators.R` | `output/tables/main/TableD1_indicator_thresholds.csv` |
| **Table S1–S5** | Sensitivity analysis tables (OAT, mortality ratio, joint grid) | `R/03_sensitivity.R` | `output/tables/supplementary/TableS*.csv` |

*Note: For maximum convenience when compiling the LaTeX manuscript, all figure files are saved with their direct manuscript names (e.g. `fig3_intensities.pdf`, `fig_perturbation_mechanics.pdf`) as well as sequential compendium names.*

---

## Directory Structure

```text
.
├── config.R                     # Central parameters, age grid, palettes, output directories
├── run_all.R                    # Master reproduction script
├── DESCRIPTION                  # R compendium dependency manifest
├── LICENSE                      # MIT license
├── README.md                    # Documentation
├── R/
│   ├── 00_io.R                  # Figure/table export helpers (PDF + PNG) & manifest writer
│   ├── 01_multistate_core.R     # Gompertz hazards, product-integral engine, demographic measures
│   ├── 02_run_baseline.R        # Baseline model & Table 1 assembly
│   ├── 03_sensitivity.R         # OAT, mortality ratio, and joint sensitivity analyses
│   ├── 04_fig_theory.R          # Theoretical mechanics plots (Figures 1 & 2)
│   ├── 05_fig_numerical.R       # Numerical illustration plots (Figures 3, 4, 5, 6, S4)
│   └── 06_fig_indicators.R      # Multi-indicator analysis (Figure D.1 & Table D.1)
└── output/                      # Generated by run_all.R
    ├── figures/
    │   ├── main/                # Publication-ready figures (.pdf and .png)
    │   └── supplementary/       # Supplementary figures
    ├── tables/
    │   ├── main/                # Main-text tables (.csv)
    │   └── supplementary/       # Supplementary tables (.csv)
    └── manifest.md              # Auto-generated catalog of all artifacts
```

---

## Calibrated Model Summary

- **Transient states**: Healthy (1), Disabled (2); Absorbing state: Dead ($\dagger$).
- **Gompertz schedules** ($x \ge 50$):
  $$\mu_{12}(x) = 0.010\,e^{0.045(x-50)}, \quad \mu_{21}(x) = 0.060\,e^{-0.020(x-50)}$$
  $$\mu_{1\dagger}(x) = 0.003\,e^{0.080(x-50)}, \quad \mu_{2\dagger}(x) = 0.0075\,e^{0.080(x-50)}$$
- **Baseline calibration target**: Healthy-start life expectancy at age 50 $e_0^{\mathrm{ms}} \approx 32.3$ years approximates the 2019 US period life table (Arias & Xu 2022, CDC/NCHS).
- **Disability excess mortality ratio**: $\mu_{2\dagger}(x) / \mu_{1\dagger}(x) = 2.5$.
