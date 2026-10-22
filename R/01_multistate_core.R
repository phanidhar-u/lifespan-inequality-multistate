# =============================================================================
# R/01_multistate_core.R
#
# Core continuous-age multistate engine:
#   - Gompertz transition intensities (Healthy / Disabled / Dead)
#   - Product integral P(50, x) via matrix exponentials
#   - Multistate survival S^ms(x), death density f(x), life expectancy e0,
#     remaining life expectancy ebar(x), life disparity e_dagger, entropy H,
#     conditional entropy Hbar(a), and threshold age a*.
# =============================================================================

#' Construct age-varying Gompertz transition intensities
construct_intensities <- function(a12, b12, a21, b21, a1d, b1d, a2d, b2d) {
  list(
    mu12 = function(x) a12 * exp( b12 * (x - AGE_MIN)),
    mu21 = function(x) a21 * exp(-b21 * (x - AGE_MIN)),
    mu1d = function(x) a1d * exp( b1d * (x - AGE_MIN)),
    mu2d = function(x) a2d * exp( b2d * (x - AGE_MIN))
  )
}

#' Assemble the 2x2 infinitesimal generator at age x
transition_generator <- function(x, intensities) {
  matrix(
    c(-(intensities$mu12(x) + intensities$mu1d(x)),  intensities$mu12(x),
       intensities$mu21(x),                          -(intensities$mu21(x) + intensities$mu2d(x))),
    nrow = 2, byrow = TRUE
  )
}

#' Build product integral transition matrices P(50, x) across grid ages
build_transition_matrices <- function(intensities, ages_vec = ages) {
  n  <- length(ages_vec)
  dx <- ages_vec[2] - ages_vec[1]

  P_list <- vector("list", n)
  P_list[[1]] <- diag(2)

  for (i in seq(2, n)) {
    step <- expm::expm(transition_generator(ages_vec[i - 1], intensities) * dx)
    P_list[[i]] <- pmax(P_list[[i - 1]] %*% step, 0)
  }
  P_list
}

#' Find crossing root between LHS and RHS curves
find_crossing <- function(lhs, rhs, ages_vec = ages) {
  diff_vec <- lhs - rhs
  # Find sign change
  idx <- which(!is.na(diff_vec) & !is.na(dplyr::lag(diff_vec)) &
                 sign(diff_vec) != sign(dplyr::lag(diff_vec)))[1]

  if (is.na(idx) || idx <= 1) {
    idx <- which.min(abs(diff_vec))
    return(ages_vec[idx])
  }

  # Linear interpolation for sub-grid precision
  x1 <- ages_vec[idx - 1]
  x2 <- ages_vec[idx]
  y1 <- diff_vec[idx - 1]
  y2 <- diff_vec[idx]

  root_approx <- x1 - y1 * (x2 - x1) / (y2 - y1)
  root_approx
}

#' Compute demographic and threshold quantities for a given initial distribution
compute_lifespan_measures <- function(transition_matrices, pi0, intensities, ages_vec = ages) {
  n    <- length(ages_vec)
  dx   <- ages_vec[2] - ages_vec[1]
  one2 <- matrix(c(1, 1), ncol = 1)
  pi0_mat <- matrix(pi0, ncol = 1)

  mu_vec <- function(x) matrix(c(intensities$mu1d(x), intensities$mu2d(x)), ncol = 1)

  # Multistate aggregate survival S^ms(x) = pi0^T P(50, x) 1
  S <- vapply(
    transition_matrices,
    function(P) as.numeric(t(pi0_mat) %*% P %*% one2),
    numeric(1)
  )
  S <- pmax(S, 1e-14)

  # Death density f(x) = pi0^T P(50, x) mu_d(x)
  f <- vapply(
    seq_len(n),
    function(i) as.numeric(t(pi0_mat) %*% transition_matrices[[i]] %*% mu_vec(ages_vec[i])),
    numeric(1)
  )
  f <- pmax(f, 0)

  # Total life expectancy from origin age 50
  e0 <- sum((S[-n] + S[-1]) / 2 * dx)

  # Remaining life expectancy at each age: ebar^ms(x) = int_x^inf S(t)dt / S(x)
  tail_S <- rev(cumsum(rev(c((S[-n] + S[-1]) / 2 * dx, 0))))
  e0x <- tail_S / S
  e0x[n] <- 0

  # Multistate life disparity: e_dagger_ms = int_0^inf f(x) ebar^ms(x) dx
  edagger <- sum((f[-n] * e0x[-n] + f[-1] * e0x[-1]) / 2 * dx)

  # Multistate life table entropy: H_ms = e_dagger / e0
  H <- edagger / e0

  # Cumulative induced attrition: Lambda^ms(a) = -ln S^ms(a)
  Lambda <- -log(S)

  # Conditional entropy past age a: Hbar^ms(a) = e_dagger,>a / (S(a) * ebar(a))
  Hbar <- vapply(seq_len(n), function(i) {
    if (S[i] < 1e-10 || e0x[i] < 1e-10) return(NA_real_)
    (sum(f[i:n] * e0x[i:n]) * dx / S[i]) / e0x[i]
  }, numeric(1))

  # Threshold condition (eq. 30 in text): Hbar^ms(a) = 1 + H_ms + ln S^ms(a)
  rhs <- 1 + H + log(S)

  # Discrete grid crossing for reporting
  idx_discrete <- which.min(abs(Hbar - rhs))
  a_star_grid  <- ages_vec[idx_discrete]

  # Subgrid continuous root
  a_star_cont  <- find_crossing(Hbar, rhs, ages_vec)

  list(
    pi0            = pi0,
    S              = S,
    f              = f,
    e0             = e0,
    e0x            = e0x,
    edagger        = edagger,
    H              = H,
    Lambda         = Lambda,
    Hbar           = Hbar,
    rhs            = rhs,
    a_star         = a_star_grid,
    a_star_exact   = a_star_cont,
    Hbar_at_astar  = Hbar[idx_discrete],
    rhs_at_astar   = rhs[idx_discrete]
  )
}
