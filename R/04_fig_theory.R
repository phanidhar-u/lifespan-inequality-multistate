# =============================================================================
# R/04_fig_theory.R
#
# Generates the theoretical figures featured in Section 3 of the paper:
#   - Figure 1: Mechanics of localized killing perturbation at a = 75 (Lemmas 2 & 3)
#   - Figure 2: Three-channel geometric decomposition of disparity sensitivity (Prop. 4)
# =============================================================================

#' Generate Figure 1: Multiplicative survival shock and remaining LE depression
generate_fig1_perturbation_mechanics <- function(baseline_intensities, ages_vec = ages) {
  n       <- length(ages_vec)
  dx      <- ages_vec[2] - ages_vec[1]
  one2    <- matrix(c(1, 1), ncol = 1)
  pi0_mat <- matrix(c(1, 0), ncol = 1)   # Healthy start

  PL_base <- build_transition_matrices(baseline_intensities, ages_vec)
  S_base  <- pmax(vapply(PL_base, function(P) as.numeric(t(pi0_mat) %*% P %*% one2), numeric(1)), 1e-14)

  tail_S_base <- rev(cumsum(rev(c((S_base[-n] + S_base[-1]) / 2 * dx, 0))))
  e_bar_base  <- tail_S_base / S_base
  e_bar_base[n] <- 0

  a_pert  <- 75
  eps_vis <- 0.35
  idx_a   <- which.min(abs(ages_vec - a_pert))

  S_pert <- S_base
  S_pert[(idx_a + 1):n] <- S_base[(idx_a + 1):n] * exp(-eps_vis)

  e_bar_pert <- numeric(n)
  for (i in seq_len(n)) {
    if (ages_vec[i] >= a_pert) {
      e_bar_pert[i] <- e_bar_base[i]
    } else {
      int_x_a   <- sum((S_base[i:(idx_a - 1)] + S_base[(i + 1):idx_a]) / 2 * dx)
      int_a_inf <- tail_S_base[idx_a]
      e_bar_pert[i] <- (int_x_a + exp(-eps_vis) * int_a_inf) / S_base[i]
    }
  }

  df_S <- tibble::tibble(age = ages_vec, S_unpert = S_base, S_pert = S_pert)
  poly_lost_LE <- tibble::tibble(
    age = c(ages_vec[idx_a:n], rev(ages_vec[idx_a:n])),
    S   = c(S_base[idx_a:n], rev(S_pert[idx_a:n]))
  )
  s_jump_val <- S_base[idx_a] * exp(-eps_vis)

  # Panel (a): Survivorship Shock (Lemma 2)
  p_a <- ggplot2::ggplot() +
    ggplot2::geom_polygon(data = poly_lost_LE, ggplot2::aes(x = age, y = S), fill = "grey85", colour = NA) +
    ggplot2::geom_line(data = df_S |> dplyr::filter(age >= a_pert), ggplot2::aes(x = age, y = S_unpert),
                       colour = "black", linetype = "dashed", linewidth = 0.75) +
    ggplot2::geom_line(data = df_S |> dplyr::filter(age <= a_pert), ggplot2::aes(x = age, y = S_pert),
                       colour = "black", linewidth = 0.85) +
    ggplot2::geom_line(data = df_S |> dplyr::filter(age > a_pert), ggplot2::aes(x = age, y = S_pert),
                       colour = "black", linewidth = 0.85) +
    ggplot2::geom_vline(xintercept = a_pert, linetype = "dotted", colour = "grey40", linewidth = 0.5) +
    ggplot2::geom_segment(ggplot2::aes(x = a_pert, y = s_jump_val, xend = a_pert, yend = S_base[idx_a]),
                          colour = "black", linewidth = 0.5) +
    ggplot2::annotate("point", x = a_pert, y = S_base[idx_a], size = 2.0, colour = "black", shape = 16) +
    ggplot2::annotate("point", x = a_pert, y = s_jump_val, size = 2.0, colour = "black", shape = 1) +
    ggplot2::annotate("text", x = a_pert, y = 0.03, label = "italic(a) == 75", parse = TRUE,
                      hjust = 0.5, size = 3.3) +
    ggplot2::labs(tag = "(a)", x = "Age x (years)",
                  y = expression("Multistate survivorship " * italic(S)^{ms}(italic(x)))) +
    ggplot2::scale_x_continuous(breaks = seq(50, 110, by = 10), limits = c(50, 110)) +
    ggplot2::scale_y_continuous(breaks = seq(0, 1.0, by = 0.2), limits = c(0, 1.04), expand = c(0, 0)) +
    ggplot2::theme_classic(base_size = 11) +
    ggplot2::theme(
      plot.tag    = ggplot2::element_text(face = "bold", size = 11),
      axis.title  = ggplot2::element_text(size = 10.5),
      axis.text   = ggplot2::element_text(size = 9.5, colour = "black"),
      plot.margin = ggplot2::margin(t = 8, r = 16, b = 8, l = 8)
    )

  # Panel (b): Remaining Life Expectancy (Lemma 3)
  df_ebar <- tibble::tibble(age = ages_vec, ebar_unpert = e_bar_base, ebar_pert = e_bar_pert)
  poly_lost_ebar <- tibble::tibble(
    age  = c(ages_vec[1:idx_a], rev(ages_vec[1:idx_a])),
    ebar = c(e_bar_base[1:idx_a], rev(e_bar_pert[1:idx_a]))
  )
  ebar_jump_val <- (tail_S_base[idx_a] * exp(-eps_vis)) / S_base[idx_a]

  p_b <- ggplot2::ggplot() +
    ggplot2::geom_polygon(data = poly_lost_ebar, ggplot2::aes(x = age, y = ebar), fill = "grey85", colour = NA) +
    ggplot2::geom_line(data = df_ebar |> dplyr::filter(age <= a_pert), ggplot2::aes(x = age, y = ebar_unpert),
                       colour = "black", linetype = "dashed", linewidth = 0.75) +
    ggplot2::geom_line(data = df_ebar |> dplyr::filter(age < a_pert), ggplot2::aes(x = age, y = ebar_pert),
                       colour = "black", linewidth = 0.85) +
    ggplot2::geom_line(data = df_ebar |> dplyr::filter(age >= a_pert), ggplot2::aes(x = age, y = ebar_pert),
                       colour = "black", linewidth = 0.85) +
    ggplot2::geom_vline(xintercept = a_pert, linetype = "dotted", colour = "grey40", linewidth = 0.5) +
    ggplot2::geom_segment(ggplot2::aes(x = a_pert, y = ebar_jump_val, xend = a_pert, yend = e_bar_base[idx_a]),
                          colour = "black", linewidth = 0.5) +
    ggplot2::annotate("point", x = a_pert, y = e_bar_base[idx_a], size = 2.0, colour = "black", shape = 16) +
    ggplot2::annotate("point", x = a_pert, y = ebar_jump_val, size = 2.0, colour = "black", shape = 1) +
    ggplot2::annotate("text", x = a_pert, y = 1.2, label = "italic(a) == 75", parse = TRUE,
                      hjust = 0.5, size = 3.3) +
    ggplot2::labs(tag = "(b)", x = "Age x (years)",
                  y = expression("Remaining life expectancy " * bar(italic(e))^{ms}(italic(x)) * " (years)")) +
    ggplot2::scale_x_continuous(breaks = seq(50, 110, by = 10), limits = c(50, 110)) +
    ggplot2::scale_y_continuous(breaks = seq(0, 35, by = 5), limits = c(0, 36), expand = c(0, 0)) +
    ggplot2::theme_classic(base_size = 11) +
    ggplot2::theme(
      plot.tag    = ggplot2::element_text(face = "bold", size = 11),
      axis.title  = ggplot2::element_text(size = 10.5),
      axis.text   = ggplot2::element_text(size = 9.5, colour = "black"),
      plot.margin = ggplot2::margin(t = 8, r = 16, b = 8, l = 8)
    )

  fig1 <- p_a + p_b + patchwork::plot_layout(ncol = 2)
  fig1
}

#' Generate Figure 2: Three-channel decomposition of disparity sensitivity (Prop. 4)
generate_fig2_prop4_structural_anatomy <- function(baseline_intensities, ages_vec = ages) {
  n       <- length(ages_vec)
  dx      <- ages_vec[2] - ages_vec[1]
  one2    <- matrix(c(1, 1), ncol = 1)
  pi0_mat <- matrix(c(1, 0), ncol = 1)

  mu_vec <- function(x) matrix(c(baseline_intensities$mu1d(x), baseline_intensities$mu2d(x)), ncol = 1)

  PL_base <- build_transition_matrices(baseline_intensities, ages_vec)
  S_base  <- pmax(vapply(PL_base, function(P) as.numeric(t(pi0_mat) %*% P %*% one2), numeric(1)), 1e-14)
  f_base  <- pmax(vapply(seq_len(n), function(i) as.numeric(t(pi0_mat) %*% PL_base[[i]] %*% mu_vec(ages_vec[i])), numeric(1)), 0)

  tail_S_base <- rev(cumsum(rev(c((S_base[-n] + S_base[-1]) / 2 * dx, 0))))
  e_bar_base  <- tail_S_base / S_base
  e_bar_base[n] <- 0

  a_pert  <- 75
  eps_vis <- 0.35
  idx_a   <- which.min(abs(ages_vec - a_pert))

  S_pert <- S_base
  S_pert[(idx_a + 1):n] <- S_base[(idx_a + 1):n] * exp(-eps_vis)

  e_bar_pert <- numeric(n)
  for (i in seq_len(n)) {
    if (ages_vec[i] >= a_pert) {
      e_bar_pert[i] <- e_bar_base[i]
    } else {
      int_x_a   <- sum((S_base[i:(idx_a - 1)] + S_base[(i + 1):idx_a]) / 2 * dx)
      int_a_inf <- tail_S_base[idx_a]
      e_bar_pert[i] <- (int_x_a + exp(-eps_vis) * int_a_inf) / S_base[i]
    }
  }

  f_pert <- f_base
  f_pert[(idx_a + 1):n] <- f_base[(idx_a + 1):n] * exp(-eps_vis)

  disp_integrand_base <- f_base * e_bar_base
  disp_integrand_pert <- f_pert * e_bar_pert

  df_integrand <- tibble::tibble(
    age        = ages_vec,
    int_unpert = disp_integrand_base,
    int_pert   = disp_integrand_pert
  )

  poly_c1 <- tibble::tibble(
    age = c(ages_vec[1:idx_a], rev(ages_vec[1:idx_a])),
    y   = c(disp_integrand_base[1:idx_a], rev(disp_integrand_pert[1:idx_a]))
  )
  poly_c2 <- tibble::tibble(
    age = c(ages_vec[idx_a:n], rev(ages_vec[idx_a:n])),
    y   = c(disp_integrand_base[idx_a:n], rev(disp_integrand_pert[idx_a:n]))
  )

  spike_height <- 0.10

  p_prop4 <- ggplot2::ggplot() +
    ggplot2::geom_polygon(data = poly_c1, ggplot2::aes(x = age, y = y), fill = "grey85", colour = NA) +
    ggplot2::geom_polygon(data = poly_c2, ggplot2::aes(x = age, y = y), fill = "grey75", colour = NA) +
    ggplot2::geom_line(data = df_integrand, ggplot2::aes(x = age, y = int_unpert),
                       colour = "black", linetype = "dashed", linewidth = 0.75) +
    ggplot2::geom_line(data = df_integrand |> dplyr::filter(age <= a_pert), ggplot2::aes(x = age, y = int_pert),
                       colour = "black", linewidth = 0.85) +
    ggplot2::geom_line(data = df_integrand |> dplyr::filter(age >= a_pert), ggplot2::aes(x = age, y = int_pert),
                       colour = "black", linewidth = 0.85) +
    ggplot2::geom_vline(xintercept = a_pert, linetype = "dotted", colour = "grey40", linewidth = 0.5) +
    ggplot2::geom_segment(ggplot2::aes(x = a_pert, y = disp_integrand_pert[idx_a],
                                       xend = a_pert, yend = disp_integrand_pert[idx_a] + spike_height),
                          arrow = ggplot2::arrow(length = ggplot2::unit(0.16, "cm"), type = "closed"),
                          colour = "black", linewidth = 0.85) +
    ggplot2::annotate("point", x = a_pert, y = disp_integrand_pert[idx_a], size = 2.0, colour = "black", shape = 16) +
    ggplot2::annotate("text", x = a_pert, y = 0.012, label = "italic(a) == 75", parse = TRUE,
                      hjust = 0.5, size = 3.3) +
    ggplot2::labs(
      x = "Age x (years)",
      y = expression("Life disparity integrand " * italic(d)^{ms}(italic(x)) * bar(italic(e))^{ms}(italic(x)))
    ) +
    ggplot2::scale_x_continuous(breaks = seq(50, 110, by = 10), limits = c(50, 110)) +
    ggplot2::scale_y_continuous(breaks = seq(0, 0.40, by = 0.10), limits = c(0, 0.44), expand = c(0, 0)) +
    ggplot2::theme_classic(base_size = 11) +
    ggplot2::theme(
      axis.title  = ggplot2::element_text(size = 10.5),
      axis.text   = ggplot2::element_text(size = 9.5, colour = "black"),
      plot.margin = ggplot2::margin(t = 12, r = 16, b = 10, l = 10)
    )

  p_prop4
}
