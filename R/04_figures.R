# =============================================================================
# R/04_figures.R
#
# Figure-generating functions for the paper:
#   Fig 1  Transition intensities
#   Fig 2  State occupancy and survival, by initial health composition
#   Fig 3  Threshold-age condition (eq. 30), by initial health composition
#   Fig 4  Sensitivity analysis (OAT + mortality-ratio panels)
#
# Each function returns a ggplot/patchwork object; run_analysis.R handles
# saving. Keeping construction and saving separate makes the figures easy
# to reuse or tweak interactively (e.g. in an R Markdown supplement).
# =============================================================================

okabe_ito_4 <- c("#0072B2", "#009E73", "#D55E00", "#CC79A7")

#' Figure 1: age-varying transition intensities on the log scale
plot_intensities <- function(intensities, ages) {
  df <- tibble::tibble(
    age  = ages,
    mu1d = intensities$mu1d(ages),
    mu12 = intensities$mu12(ages),
    mu2d = intensities$mu2d(ages),
    mu21 = intensities$mu21(ages)
  ) |>
    tidyr::pivot_longer(-age, names_to = "rate", values_to = "value") |>
    dplyr::mutate(rate = factor(
      rate,
      levels = c("mu1d", "mu12", "mu2d", "mu21"),
      labels = c(
        "mu[1*'\u2020']~':'~Healthy %->% Dead",
        "mu[12]~':'~Healthy %->% Disabled",
        "mu[2*'\u2020']~':'~Disabled %->% Dead",
        "mu[21]~':'~Disabled %->% Healthy"
      )
    ))

  ggplot2::ggplot(df, ggplot2::aes(age, value, colour = rate, linetype = rate)) +
    ggplot2::geom_line(linewidth = 1.1) +
    ggplot2::scale_y_log10(
      breaks = c(0.003, 0.01, 0.03, 0.1, 0.3),
      labels = c("0.003", "0.010", "0.030", "0.100", "0.300")
    ) +
    ggplot2::scale_colour_manual(values = okabe_ito_4, labels = scales::label_parse()) +
    ggplot2::scale_linetype_manual(
      values = c("solid", "dashed", "solid", "dashed"),
      labels = scales::label_parse()
    ) +
    ggplot2::labs(
      title = "Transition intensities (model parameters)",
      x = "Age", y = "Intensity (log scale)", colour = NULL, linetype = NULL
    ) +
    ggplot2::theme_classic(base_size = 12) +
    ggplot2::theme(legend.position = "bottom", legend.text = ggplot2::element_text(size = 10))
}

#' A single state-occupancy panel for one initial composition
occupancy_panel <- function(pi0, title_str, transition_matrices, ages) {
  pi0_mat <- matrix(pi0, ncol = 1)
  state_probs <- do.call(rbind, lapply(transition_matrices, function(P) as.numeric(t(pi0_mat) %*% P)))

  df <- tibble::tibble(age = ages, p1 = state_probs[, 1], p2 = state_probs[, 2]) |>
    dplyr::mutate(S = p1 + p2)
  df_long <- tidyr::pivot_longer(df, c(p1, p2), names_to = "state", values_to = "prob")

  ggplot2::ggplot() +
    ggplot2::geom_line(data = df_long, ggplot2::aes(age, prob, colour = state), linewidth = 1.2) +
    ggplot2::geom_line(data = df, ggplot2::aes(age, S), colour = "black",
                        linetype = "dashed", linewidth = 1) +
    ggplot2::scale_colour_manual(
      values = c(p1 = "#0072B2", p2 = "#D55E00"),
      labels = c("Healthy", "Disabled")
    ) +
    ggplot2::coord_cartesian(ylim = c(0, 1)) +
    ggplot2::labs(title = title_str, x = "Age", y = "Probability", colour = NULL) +
    ggplot2::theme_classic(base_size = 12) +
    ggplot2::theme(plot.title = ggplot2::element_text(face = "bold", size = 11),
                   legend.position = "top")
}

#' Figure 2: state occupancy and survival across the three initial compositions
plot_occupancy <- function(transition_matrices, ages) {
  panels <- wrap_plots(
    occupancy_panel(INITIAL_STATES$healthy,  "Started healthy (1,0)",          transition_matrices, ages),
    occupancy_panel(INITIAL_STATES$disabled, "Started disabled (0,1)",         transition_matrices, ages),
    occupancy_panel(INITIAL_STATES$mixed,    "Started mixed (0.75,0.25)",      transition_matrices, ages),
    nrow = 1
  )

  panels + patchwork::plot_annotation(
    title = "State occupancy and survival under different initial health compositions",
    subtitle = "Solid lines: state-specific probabilities; dashed: overall survival S^ms(x)",
    caption = paste(
      "Source: Simulated multistate model with age-varying transition intensities.",
      "\nProbabilities derived from the product-integral of the transition intensity matrix."
    ),
    theme = ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = 13),
      plot.subtitle = ggplot2::element_text(size = 11),
      plot.caption = ggplot2::element_text(size = 9, colour = "grey40")
    )
  )
}

threshold_panel_colours <- list(
  healthy  = list(lhs = "#1b6ca8", rhs = "#e07b39"),
  disabled = list(lhs = "#c0392b", rhs = "#8e44ad"),
  mixed    = list(lhs = "#27ae60", rhs = "#d35400")
)

threshold_panel_titles <- list(
  healthy  = expression(bold(pi)[0] * " = (1, 0)   [Healthy start]"),
  disabled = expression(bold(pi)[0] * " = (0, 1)   [Disabled start]"),
  mixed    = expression(bold(pi)[0] * " = (0.75, 0.25)   [Mixed]")
)

#' A single threshold-condition panel for one initial composition
threshold_panel <- function(result, state_name, ages) {
  df <- tibble::tibble(age = ages, Hbar = result$Hbar, rhs = result$rhs)
  cols <- threshold_panel_colours[[state_name]]
  a_label <- paste0("italic(a)^'*' == ", round(result$a_star, 1))
  H_label <- paste0("italic(H)[ms] == ", round(result$H, 3))

  ggplot2::ggplot(df, ggplot2::aes(x = age)) +
    ggplot2::geom_line(ggplot2::aes(y = Hbar, colour = "LHS"), linewidth = 1.2, na.rm = TRUE) +
    ggplot2::geom_line(ggplot2::aes(y = rhs, colour = "RHS"),
                        linetype = "dashed", linewidth = 1.2, na.rm = TRUE) +
    ggplot2::geom_vline(xintercept = result$a_star, linetype = "dotted",
                         linewidth = 0.9, colour = "grey30") +
    ggplot2::annotate("text", x = result$a_star + 0.6,
                       y = max(result$Hbar, na.rm = TRUE) * 0.88,
                       label = a_label, parse = TRUE, hjust = -0.3, vjust = -0.3, size = 3.4) +
    ggplot2::annotate("text", x = 55, y = min(result$rhs, na.rm = TRUE) * 0.55,
                       label = H_label, parse = TRUE, hjust = 0, size = 3.1, colour = "grey45") +
    ggplot2::scale_colour_manual(
      values = c(LHS = cols$lhs, RHS = cols$rhs),
      labels = c(
        expression(bar(italic(H))^{ms}(italic(a))),
        expression(1 + italic(H)[ms] + ln ~ italic(S)^{ms}(italic(a)))
      )
    ) +
    ggplot2::coord_cartesian(xlim = c(50, 108)) +
    ggplot2::labs(title = threshold_panel_titles[[state_name]], x = "Age", y = "Value", colour = NULL) +
    ggplot2::theme_classic(base_size = 11) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = 10),
      legend.position = "bottom",
      legend.text = ggplot2::element_text(size = 9),
      legend.key.width = ggplot2::unit(1.2, "cm")
    )
}

#' Figure 3: threshold-age condition, stacked across the three compositions
plot_threshold <- function(results, ages) {
  panels <- threshold_panel(results$healthy, "healthy", ages) /
    threshold_panel(results$disabled, "disabled", ages) /
    threshold_panel(results$mixed, "mixed", ages)

  panels + patchwork::plot_annotation(
    title = expression(
      "Threshold condition: " * bar(italic(H))^{ms}(italic(a)^"*") *
        " = 1 + " * italic(H)[ms] * " + ln " * italic(S)^{ms}(italic(a)^"*")
    ),
    subtitle = "Solid = conditional entropy past age a.  Dashed = RHS of eq. (30).  Dotted = a*.",
    theme = ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = 12),
      plot.subtitle = ggplot2::element_text(size = 9, colour = "grey40")
    )
  )
}

state_colours_3 <- c(
  "Healthy (1,0)"      = "#1b6ca8",
  "Disabled (0,1)"     = "#c0392b",
  "Mixed (0.75,0.25)"  = "#27ae60"
)

#' Figure 4: sensitivity analysis (OAT panel + mortality-ratio panel)
plot_sensitivity <- function(sensitivity_oat, sensitivity_ratio, baseline_results) {
  reference_lines <- tibble::tibble(
    init_state = factor(
      c("Healthy (1,0)", "Disabled (0,1)", "Mixed (0.75,0.25)"),
      levels = c("Healthy (1,0)", "Disabled (0,1)", "Mixed (0.75,0.25)")
    ),
    ref = c(baseline_results$healthy$a_star, baseline_results$disabled$a_star, baseline_results$mixed$a_star)
  )

  panel_a <- ggplot2::ggplot(
    sensitivity_oat,
    ggplot2::aes(x = mult_label, y = a_star, colour = init_state, group = init_state)
  ) +
    ggplot2::geom_hline(data = reference_lines, ggplot2::aes(yintercept = ref, colour = init_state),
                         linetype = "dotted", linewidth = 0.5, alpha = 0.6) +
    ggplot2::geom_line(linewidth = 1.0) +
    ggplot2::geom_point(ggplot2::aes(size = is_baseline)) +
    ggplot2::scale_size_manual(values = c("FALSE" = 2.0, "TRUE" = 3.5), guide = "none") +
    ggplot2::scale_colour_manual(values = state_colours_3) +
    ggplot2::facet_wrap(~parameter, nrow = 2) +
    ggplot2::labs(
      title = "A.  One-at-a-time parameter variation",
      x = "Change from baseline",
      y = expression("Threshold age " * italic(a)^"*" * " (years)"),
      colour = "Initial state"
    ) +
    ggplot2::theme_classic(base_size = 11) +
    ggplot2::theme(
      strip.text = ggplot2::element_text(face = "bold", size = 9),
      legend.position = "bottom",
      plot.title = ggplot2::element_text(face = "bold", size = 11)
    )

  panel_b <- ggplot2::ggplot(
    sensitivity_ratio,
    ggplot2::aes(x = ratio, y = a_star, colour = init_state, group = init_state)
  ) +
    ggplot2::geom_vline(xintercept = 2.5, linetype = "dashed", colour = "grey55", linewidth = 0.7) +
    ggplot2::annotate("text", x = 2.5, y = max(sensitivity_ratio$a_star) - 0.3,
                       label = "Baseline\nratio = 2.5", size = 2.9, hjust = -0.08, colour = "grey40") +
    ggplot2::geom_line(linewidth = 1.1) +
    ggplot2::geom_point(size = 2.8) +
    ggplot2::scale_colour_manual(values = state_colours_3) +
    ggplot2::scale_x_continuous(breaks = c(1.5, 2.0, 2.5, 3.0, 3.5, 4.0)) +
    ggplot2::labs(
      title = expression("B.  Disabled/healthy mortality ratio (" * mu[{2 * "\u2020"}] * " / " * mu[{1 * "\u2020"}] * ")"),
      x = "Mortality ratio",
      y = expression("Threshold age " * italic(a)^"*" * " (years)"),
      colour = "Initial state"
    ) +
    ggplot2::theme_classic(base_size = 11) +
    ggplot2::theme(legend.position = "bottom", plot.title = ggplot2::element_text(face = "bold", size = 11))

  (panel_a / panel_b) + patchwork::plot_annotation(
    title = expression("Sensitivity of threshold age " * italic(a)^"*" * " to parameter variation"),
    caption = paste(
      "Panel A: each parameter varied by -50%, -25%, baseline, +25%, +50%; all others held fixed.",
      "\nLarge points mark the baseline value; dotted horizontal lines show baseline a*.",
      "\nPanel B: excess mortality ratio varied with mu_{1\u2020} held fixed at its calibrated value."
    ),
    theme = ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = 13),
      plot.caption = ggplot2::element_text(size = 8.5, colour = "grey40")
    )
  )
}

#' FigS4: joint sensitivity heatmap of a* over mortality ratio x disability incidence
plot_joint_sensitivity <- function(sensitivity_joint) {
  ggplot2::ggplot(sensitivity_joint, ggplot2::aes(x = ratio, y = mult, fill = a_star)) +
    ggplot2::geom_tile() +
    ggplot2::geom_contour(ggplot2::aes(z = a_star), colour = "white", linewidth = 0.3, alpha = 0.6) +
    ggplot2::facet_wrap(~init_state, nrow = 1) +
    ggplot2::scale_fill_viridis_c(name = expression(italic(a)^"*")) +
    ggplot2::labs(
      title = expression("Joint sensitivity of threshold age " * italic(a)^"*"),
      subtitle = "Simultaneous variation in the disabled/healthy mortality ratio and disability-incidence scale",
      x = expression("Mortality ratio  " * mu[{2 * "\u2020"}] / mu[{1 * "\u2020"}]),
      y = expression("Disability incidence multiplier (" * a[12] * ")")
    ) +
    ggplot2::theme_classic(base_size = 11) +
    ggplot2::theme(
      strip.text = ggplot2::element_text(face = "bold", size = 9),
      plot.title = ggplot2::element_text(face = "bold", size = 12),
      plot.subtitle = ggplot2::element_text(size = 9, colour = "grey40")
    )
}
