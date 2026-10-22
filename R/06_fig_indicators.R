# =============================================================================
# R/06_fig_indicators.R
#
# Generates the multi-indicator comparison analysis from Appendix D:
#   - Evaluates threshold conditions for 5 lifespan inequality indicators:
#       1. Multistate Life disparity (e†)
#       2. Multistate Life table entropy (H)
#       3. Multistate Lifespan variance (σ²)
#       4. Multistate Coefficient of variation (CV)
#       5. Multistate Drewnowski's index (D)
#   - Produces Table D.1 (TableD1_indicator_thresholds.csv)
#   - Produces Figure D.1 (fig_threshold_comparison_compact.pdf / .png)
# =============================================================================

compute_all_indicator_thresholds <- function(transition_matrices, pi0, intensities, ages_vec = ages) {
  n       <- length(ages_vec)
  dx      <- ages_vec[2] - ages_vec[1]
  one2    <- matrix(c(1, 1), ncol = 1)
  pi0_mat <- matrix(pi0, ncol = 1)
  mu_vec  <- function(x) matrix(c(intensities$mu1d(x), intensities$mu2d(x)), ncol = 1)

  # Multistate aggregate survival S^ms(x) and death density f(x)
  S <- pmax(vapply(transition_matrices, function(P) as.numeric(t(pi0_mat) %*% P %*% one2), numeric(1)), 1e-14)
  f <- pmax(vapply(seq_len(n), function(i) as.numeric(t(pi0_mat) %*% transition_matrices[[i]] %*% mu_vec(ages_vec[i])), numeric(1)), 0)

  # Life expectancy and remaining life expectancy
  e0     <- sum((S[-n] + S[-1]) / 2 * dx)
  tail_S <- rev(cumsum(rev(c((S[-n] + S[-1]) / 2 * dx, 0))))
  e0x    <- tail_S / S
  e0x[n] <- 0

  # Life disparity and entropy
  edagger <- sum((f[-n] * e0x[-n] + f[-1] * e0x[-1]) / 2 * dx)
  H_ms    <- edagger / e0

  # Cumulative attrition and conditional entropy
  Lambda   <- -log(S)
  Hbar_num <- rev(cumsum(rev(c((f[-n] * e0x[-n] + f[-1] * e0x[-1]) / 2 * dx, 0))))
  Hbar     <- ifelse(S < 1e-10 | e0x < 1e-10, NA_real_, (Hbar_num / S) / e0x)

  # Lifespan variable t = x - AGE_MIN
  t_ages <- ages_vec - AGE_MIN

  # Variance of remaining lifespan sigma2 = E[t^2] - e0^2
  Et2    <- sum((t_ages[-n]^2 * f[-n] + t_ages[-1]^2 * f[-1]) / 2 * dx)
  sigma2 <- Et2 - e0^2

  # Stationary population mean age above a: a_bar_ms(a) = int_a^inf t S(t) dt / int_a^inf S(t) dt
  tS_integrand <- t_ages * S
  tail_tS      <- rev(cumsum(rev(c((tS_integrand[-n] + tS_integrand[-1]) / 2 * dx, 0))))
  a_bar_ms     <- tail_tS / tail_S

  # Drewnowski index
  S2         <- S^2
  D_ms_total <- sum((S2[-n] + S2[-1]) / 2 * dx) / e0
  tail_S2    <- rev(cumsum(rev(c((S2[-n] + S2[-1]) / 2 * dx, 0))))
  D_ms_a     <- ifelse(S < 1e-10 | tail_S < 1e-10, NA_real_, (1 / S) * (tail_S2 / tail_S))

  # Root finder
  find_root <- function(lhs, rhs) {
    diff <- lhs - rhs
    idx  <- which(!is.na(diff) & !is.na(dplyr::lag(diff)) &
                    sign(diff) != sign(dplyr::lag(diff)))[1]
    if (is.na(idx) || idx == 1) idx <- which.min(abs(diff))
    ages_vec[idx]
  }

  list(
    e0           = e0,
    edagger      = edagger,
    H_ms         = H_ms,
    sigma2       = sigma2,
    D_ms         = D_ms_total,
    a_edagger    = find_root(Lambda + Hbar, rep(1, n)),
    a_entropy    = find_root(Lambda + Hbar, rep(1 + H_ms, n)),
    a_var        = find_root(a_bar_ms,      rep(e0, n)),
    a_cv         = find_root(a_bar_ms,      rep(e0 + sigma2 / e0, n)),
    a_drewnowski = find_root(2 * S * D_ms_a, rep(D_ms_total, n))
  )
}

generate_indicator_comparison <- function(baseline_transitions, baseline_intensities, ages_vec = ages) {
  cases <- list(
    "Healthy (1,0)"      = INITIAL_STATES$healthy,
    "Disabled (0,1)"     = INITIAL_STATES$disabled,
    "Mixed (0.75, 0.25)" = INITIAL_STATES$mixed
  )

  results_by_case <- lapply(cases, function(pi0) {
    compute_all_indicator_thresholds(baseline_transitions, pi0, baseline_intensities, ages_vec)
  })

  indicators_display <- c(
    "Life disparity\n(e\u2020)",
    "Life table\nentropy (H)",
    "Lifespan\nvariance (\u03c3\u00b2)",
    "Coeff. of\nvariation (CV)",
    "Drewnowski's\nindex (D)"
  )

  df_indicators <- dplyr::bind_rows(lapply(names(results_by_case), function(nm) {
    r <- results_by_case[[nm]]
    tibble::tibble(
      init_state = nm,
      indicator  = indicators_display,
      a_star     = c(r$a_edagger, r$a_entropy, r$a_var, r$a_cv, r$a_drewnowski)
    )
  })) |>
    dplyr::mutate(
      init_state = factor(init_state, levels = names(cases)),
      indicator  = factor(indicator,  levels = indicators_display)
    )

  # Table D.1 formatted for export
  table_d1 <- df_indicators |>
    dplyr::mutate(indicator = gsub("\n", " ", indicator)) |>
    tidyr::pivot_wider(names_from = init_state, values_from = a_star)

  save_table(
    table_d1, "TableD1_indicator_thresholds",
    caption = "Distinct threshold ages for five lifespan inequality indicators across initial compositions (Table D.1)."
  )

  # Plot Figure D.1 matching publication layout
  cb3    <- c("Healthy (1,0)" = "#1b6ca8", "Disabled (0,1)" = "#c0392b", "Mixed (0.75, 0.25)" = "#27ae60")
  shape3 <- c("Healthy (1,0)" = 16, "Disabled (0,1)" = 17, "Mixed (0.75, 0.25)" = 15)

  fig_d1 <- ggplot2::ggplot(df_indicators,
                            ggplot2::aes(x = indicator, y = a_star,
                                         colour = init_state, shape = init_state,
                                         group = init_state)) +
    ggplot2::geom_line(linewidth = 0.8, alpha = 0.55) +
    ggplot2::geom_point(size = 4.2) +
    ggplot2::scale_colour_manual(values = cb3) +
    ggplot2::scale_shape_manual(values = shape3) +
    ggplot2::scale_y_continuous(
      breaks = scales::pretty_breaks(n = 8),
      limits = c(
        floor(min(df_indicators$a_star)   - 1.5),
        ceiling(max(df_indicators$a_star) + 1.5)
      )
    ) +
    ggplot2::labs(
      x       = "Lifespan inequality indicator",
      y       = expression("Threshold age " * italic(a)^"*" * " (years)"),
      colour  = "Initial composition",
      shape   = "Initial composition"
    ) +
    ggplot2::theme_classic(base_size = 12) +
    ggplot2::theme(
      legend.position    = "bottom",
      legend.title       = ggplot2::element_text(size = 10),
      legend.text        = ggplot2::element_text(size = 10),
      axis.text.x        = ggplot2::element_text(size = 9.5, lineheight = 1.2),
      panel.grid.major.y = ggplot2::element_line(colour = "grey88", linewidth = 0.35),
      panel.grid.minor.y = ggplot2::element_line(colour = "grey93", linewidth = 0.20)
    )

  list(df = df_indicators, table = table_d1, plot = fig_d1)
}
