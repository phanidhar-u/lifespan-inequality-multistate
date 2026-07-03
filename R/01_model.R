# =============================================================================
# R/01_model.R
#
# Core two-state (Healthy / Disabled / Dead) multistate model:
#   - construct_intensities(): age-varying Gompertz-type transition intensities
#   - build_transition_matrices(): discretised product-integral P(50, x)
#   - compute_lifespan_measures(): survival, lifespan entropy, and the
#     threshold-age condition of eq. (30) in the paper
#
# State labels: 1 = Healthy, 2 = Disabled, "d" = Dead (absorbing).
# =============================================================================

#' Construct age-varying transition intensities
#'
#' Builds the four Gompertz-type intensity functions governing the
#' Healthy-Disabled-Dead process: mu12 (disability incidence), mu21
#' (recovery), mu1d and mu2d (cause-specific mortality). All intensities
#' are parameterised relative to age 50, the left endpoint of the age grid.
#'
#' @param a12,b12 Scale and shape of the healthy-to-disabled intensity.
#' @param a21,b21 Scale and shape of the disabled-to-healthy intensity.
#' @param a1d,b1d Scale and shape of healthy-state mortality.
#' @param a2d,b2d Scale and shape of disabled-state mortality.
#' @return A named list of four functions of age, `mu12`, `mu21`, `mu1d`,
#'   `mu2d`.
construct_intensities <- function(a12, b12, a21, b21, a1d, b1d, a2d, b2d) {
  list(
    mu12 = function(x) a12 * exp( b12 * (x - 50)),
    mu21 = function(x) a21 * exp(-b21 * (x - 50)),
    mu1d = function(x) a1d * exp( b1d * (x - 50)),
    mu2d = function(x) a2d * exp( b2d * (x - 50))
  )
}

#' Assemble the infinitesimal generator at a single age
#'
#' @param x Age at which the generator is evaluated.
#' @param intensities A list of intensity functions from
#'   [construct_intensities()].
#' @return A 2x2 generator matrix `A(x)` with rows summing to (approximately)
#'   zero once mortality outflow is included.
transition_generator <- function(x, intensities) {
  matrix(
    c(-(intensities$mu12(x) + intensities$mu1d(x)),  intensities$mu12(x),
       intensities$mu21(x),                          -(intensities$mu21(x) + intensities$mu2d(x))),
    nrow = 2, byrow = TRUE
  )
}

#' Build the discretised product-integral P(age_min, x) for every grid age
#'
#' Approximates the product-integral (Aalen-Johansen estimator) of the
#' generator matrix via piecewise-constant matrix exponentials on a fine
#' age grid, following standard multistate life table methodology
#' (Andersen et al. 1993; Schoen 1988).
#'
#' @param intensities A list of intensity functions from
#'   [construct_intensities()].
#' @param ages Numeric vector of grid ages, evenly spaced.
#' @return A list of length `length(ages)`, where element `i` is the 2x2
#'   transition matrix P(ages[1], ages[i]).
build_transition_matrices <- function(intensities, ages) {
  n  <- length(ages)
  dx <- ages[2] - ages[1]

  matrices <- vector("list", n)
  matrices[[1]] <- diag(2)

  for (i in seq(2, n)) {
    step <- expm::expm(transition_generator(ages[i - 1], intensities) * dx)
    matrices[[i]] <- pmax(matrices[[i - 1]] %*% step, 0)
  }
  matrices
}

#' Compute survival, entropy, and threshold-age quantities
#'
#' Given a discretised set of transition matrices and an initial health
#' composition `pi0`, this returns the full set of demographic quantities
#' used in the paper: multistate survival S^ms(x), the density of deaths
#' f(x), life expectancy e_0^ms, life-table entropy H_ms, the
#' age-conditional entropy \eqn{\bar H^{ms}(a)}, and the threshold age
#' a* solving eq. (30):
#'   \eqn{\bar H^{ms}(a^*) = 1 + H_{ms} + \ln S^{ms}(a^*)}.
#'
#' @param transition_matrices Output of [build_transition_matrices()].
#' @param pi0 Numeric length-2 vector, the initial (age-50) distribution
#'   across the Healthy and Disabled states. Must sum to 1.
#' @param intensities A list of intensity functions from
#'   [construct_intensities()], used to evaluate cause-specific mortality.
#' @param ages Numeric vector of grid ages matching `transition_matrices`.
#' @return A named list with elements `S`, `f`, `e0`, `e0x`, `edagger`, `H`,
#'   `Hbar`, `rhs`, `a_star`, `Hbar_at_astar`, `rhs_at_astar`.
compute_lifespan_measures <- function(transition_matrices, pi0, intensities, ages) {
  n  <- length(ages)
  dx <- ages[2] - ages[1]
  one2 <- matrix(c(1, 1), ncol = 1)
  pi0_mat <- matrix(pi0, ncol = 1)

  mortality_vector <- function(x) matrix(c(intensities$mu1d(x), intensities$mu2d(x)), ncol = 1)

  # Multistate survival S^ms(x) = pi0' P(50, x) 1
  S <- vapply(
    transition_matrices,
    function(P) as.numeric(t(pi0_mat) %*% P %*% one2),
    numeric(1)
  )
  S <- pmax(S, 1e-14)

  # Density of deaths f(x) = pi0' P(50, x) mu(x)
  f <- vapply(
    seq_len(n),
    function(i) as.numeric(t(pi0_mat) %*% transition_matrices[[i]] %*% mortality_vector(ages[i])),
    numeric(1)
  )
  f <- pmax(f, 0)

  # Life expectancy at age 50 (trapezoidal quadrature)
  e0 <- sum((S[-n] + S[-1]) / 2 * dx)

  # Remaining life expectancy at each age, e_0^ms(x)
  tail_area <- rev(cumsum(rev(c((S[-n] + S[-1]) / 2 * dx, 0))))
  e0x <- tail_area / S
  e0x[n] <- 0

  # Keyfitz-type lifespan entropy, H_ms = e-dagger / e0
  edagger <- sum((f[-n] * e0x[-n] + f[-1] * e0x[-1]) / 2 * dx)
  H <- edagger / e0

  # Age-conditional entropy bar{H}^{ms}(a), evaluated at every grid age
  Hbar <- vapply(seq_len(n), function(i) {
    if (S[i] < 1e-10 || e0x[i] < 1e-10) return(NA_real_)
    (sum(f[i:n] * e0x[i:n]) * dx / S[i]) / e0x[i]
  }, numeric(1))

  # Threshold condition (eq. 30): bar{H}^{ms}(a) = 1 + H_ms + ln S^{ms}(a)
  rhs <- 1 + H + log(S)
  idx <- which.min(abs(Hbar - rhs))
  a_star <- ages[idx]

  list(
    S = S, f = f, e0 = e0, e0x = e0x, edagger = edagger, H = H,
    Hbar = Hbar, rhs = rhs, a_star = a_star,
    Hbar_at_astar = Hbar[idx], rhs_at_astar = rhs[idx]
  )
}
