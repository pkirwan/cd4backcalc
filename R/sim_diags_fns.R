#' Generate infections according to specified pattern
#'
#' @param nquar Number of time quarters.
#' @param nlatst Number of latent states.
#' @param pattern Pattern type for infections ("varying", "varying_2", "constant", "increasing").
#' @param target_inf Target total number of infections to simulate (may be higher/lower with Poisson variation).
#' @param dist Logical indicating whether to add Poisson variation.
#' @return Matrix of new infections over time.
generate_infections <- function(
  nquar,
  nlatst,
  pattern,
  target_inf = 90000,
  dist = TRUE
) {
  # check pattern
  stopifnot(pattern %in% c("varying", "varying_1", "varying_2", "constant", "increasing"))

  ann_inf <- target_inf / nquar
  infs <- if (pattern == "varying") {
    c(
      seq(1.2, 2, length.out = round(nquar * 0.1)),
      seq(2, 1, length.out = round(nquar * 7 / 15)),
      seq(1, 0.5, length.out = round(nquar * 2 / 15)),
      seq(0.5, 0.4, length.out = round(nquar * 0.3))
    ) *
      ann_inf
  } else if (pattern == "varying_1") {
    c(
      seq(0.5, 0.6, length.out = round(nquar * 0.1)),
      seq(0.6, 0.6, length.out = round(nquar * 4 / 15)),
      seq(0.6, 0.5, length.out = round(nquar * 0.2)),
      seq(0.5, 0.4, length.out = round(nquar * 2 / 15)),
      seq(0.4, 0.45, length.out = round(nquar * 0.3))
    ) *
      ann_inf
  } else if (pattern == "varying_2") {
    c(
      seq(0.8, 1.6, length.out = round(nquar * 0.6)),
      seq(1.6, 0.9, length.out = round(nquar * 0.4))
    ) * ann_inf
  } else if (pattern == "constant") {
    rep(ann_inf, nquar)
  } else if (pattern == "increasing") {
    seq(0.5, 1.5, length.out = nquar) * ann_inf
  }
  infs <- stats::loess(infs ~ seq_along(infs), span = 0.4)$fitted
  infs <- if (dist) stats::rpois(nquar, infs) else infs
  h <- matrix(
    c(infs, rep(0, nquar * (nlatst))),
    nrow = nquar,
    ncol = nlatst + 1,
    byrow = FALSE
  )
  return(h)
}

#' Generate migration patterns
#'
#' @param nquar Number of time quarters.
#' @param pattern Pattern type for omega ("constant", "increasing", "step_change_1", "step_change_2", "none").
#' @param target_mig Target total number of migrations to simulate (may be higher/lower with Poisson variation).
#' @param dist Logical indicating whether to add Poisson variation.
#' @return Vector of omega values over time.
generate_omega <- function(nquar, pattern, target_mig = 12000, dist = TRUE) {
  # check pattern
  stopifnot(pattern %in% c("constant", "increasing", "step_change_1", "step_change_2", "none"))
  ann_mig <- target_mig / nquar
  omega <- if (pattern == "constant") {
    rep(ann_mig, nquar)
  } else if (pattern == "increasing") {
    seq(0.5, 1.5, length.out = nquar) * ann_mig
  } else if (pattern == "step_change_1") {
    c(
      seq(0.4, 0.4, length.out = round(nquar * 0.1)),
      seq(0.4, 1.2, length.out = round(nquar * 0.1)),
      seq(1.2, 0.4, length.out = round(nquar * 0.1)),
      seq(0.4, 1.3, length.out = round(nquar * 0.55)),
      seq(1.35, 2.5, length.out = round(nquar * 0.15))
    ) *
      ann_mig
  } else if (pattern == "step_change_2") {
    c(
      seq(0.7, 1.0, length.out = round(nquar * 0.65)),
      seq(1.0, 2.8, length.out = round(nquar * 0.20)),
      seq(2.8, 3.0, length.out = round(nquar * 0.15))
    ) * ann_mig
  } else if (pattern == "none") {
    rep(0, nquar)
  }
  omega <- stats::loess(omega ~ seq_along(omega), span = 0.4)$fitted
  omega <- if (dist) stats::rpois(nquar, omega) else omega
  return(omega)
}

#' Generate migration probability patterns
#'
#' @param nquar Number of time quarters.
#' @param nlatst Number of latent states.
#' @param pattern Pattern type for mu ("constant", "differing", "varying").
#' @param dist Logical indicating whether to add logit-scale noise.
#' @return Matrix of mu values over time.
generate_mu <- function(nquar, nlatst, pattern, dist = FALSE) {
  # check pattern
  stopifnot(pattern %in% c("constant", "differing", "varying"))
  if (pattern == "constant") {
    mu <- matrix(0.2, nrow = nquar, ncol = nlatst + 1)
  } else {
    mu <- if (pattern == "differing") {
      matrix(
        c(0.1, 0.2, 0.2, 0.1, 0.4),
        nrow = nquar,
        ncol = nlatst + 1,
        byrow = TRUE
      )
    } else if (pattern == "varying") {
      matrix(
        c(
          seq(0.15, 0.3, length.out = nquar),
          seq(0.2, 0.22, length.out = nquar),
          seq(0.2, 0.15, length.out = nquar),
          seq(0.10, 0.08, length.out = nquar),
          seq(0.2, 0.4, length.out = nquar)
        ),
        nrow = nquar,
        ncol = nlatst + 1
      )
    }
  }

  if (dist) {
    sigma <- 0.1
    # Add noise on logit scale with compartment nlatst as reference,
    # matching Stan's softmax parameterization
    log_mu <- log(mu)
    logits <- sweep(log_mu, 1, log_mu[, nlatst])
    noise_idx <- setdiff(seq_len(nlatst + 1), nlatst)
    logits[, noise_idx] <- logits[, noise_idx] +
      matrix(
        stats::rnorm(nquar * length(noise_idx), 0, sigma),
        nrow = nquar
      )
    exp_logits <- exp(logits)
    mu <- exp_logits / rowSums(exp_logits)
  } else {
    mu <- mu / rowSums(mu)
  }

  return(mu)
}

#' Generate proportion UK patterns
#'
#' @param nquar Number of time quarters.
#' @param pattern Pattern type for proportion UK ("constant", "increasing").
#' @param dist Logical indicating whether to add normal variation.
#' @return Vector of proportion UK values over time.
generate_ratio_abroad_uk <- function(nquar, pattern, dist = TRUE) {
  # check pattern
  stopifnot(pattern %in% c("constant", "increasing"))
  if (dist) {
    sigma <- 0.08
    if (pattern == "constant") {
      log_ratio_abroad_uk <- stats::rnorm(nquar, rep(log(0.2), nquar), sigma)
    } else if (pattern == "increasing") {
      log_ratio_abroad_uk <- stats::rnorm(
        nquar,
        seq(log(0.1), log(0.4), length.out = nquar),
        sigma
      )
    }
    # Transform back to proportion scale
    ratio_abroad_uk <- exp(log_ratio_abroad_uk)
  } else if (pattern == "constant") {
    ratio_abroad_uk <- rep(0.2, nquar)
  } else if (pattern == "increasing") {
    ratio_abroad_uk <- exp(seq(log(0.1), log(0.4), length.out = nquar))
  }
  return(ratio_abroad_uk)
}

#' Generate migration combining omega and mu
#'
#' @param nquar Number of time quarters.
#' @param nlatst Number of latent states.
#' @param omega Vector of omega values over time.
#' @param mu Matrix of mu values over time.
#' @param dist Logical indicating whether to use multinomial distribution.
#' @return Matrix of migration counts over time.
generate_migration <- function(nquar, nlatst, omega, mu, dist = TRUE) {
  m <- matrix(0, nrow = nquar, ncol = nlatst + 1)
  if (dist) {
    for (i in 1:nquar) {
      m[i, ] <- stats::rmultinom(1, omega[i], mu[i, ])
    }
  } else {
    for (i in 1:nquar) {
      m[i, ] <- omega[i] * mu[i, ]
    }
  }
  return(m)
}

#' Simulate AR(1) trajectory on the logit scale
#'
#' Produces a time series with the same correlation structure as the Stan AR(1)
#' random-walk prior used for diagnosis probabilities. AR(1) residuals are added
#' on top of a deterministic trend and then mean-centred so that the stochastic
#' trajectory does not drift systematically from the deterministic trend. This
#' ensures that \code{dist = TRUE} and \code{dist = FALSE} simulations produce
#' comparable epidemic trajectories (e.g. undiagnosed prevalence), while still
#' exhibiting realistic smooth quarter-to-quarter variation.
#'
#' @param nquar Number of time quarters.
#' @param trend_logit Numeric vector of length \code{nquar} giving the
#'   deterministic mean trajectory on the logit scale.
#' @param rho AR(1) autocorrelation coefficient.
#' @param vardelta Innovation variance (step variance per quarter).
#' @return Numeric vector of length \code{nquar} on the logit scale.
#' @keywords internal
simulate_ar1_logit <- function(nquar, trend_logit, rho = 0.86, vardelta = exp(-4.0)) {
  resid <- numeric(nquar)
  resid[1] <- sqrt(vardelta / (1 - rho^2)) * stats::rnorm(1)
  for (t in 2:nquar) {
    resid[t] <- rho * resid[t - 1] + sqrt(vardelta) * stats::rnorm(1)
  }
  # Centre residuals so the stochastic trajectory has the same mean as the
  # deterministic trend. Without this, a large initial draw (amplified by the
  # stationary-distribution scale factor 1/sqrt(1-rho^2) ~ 2 for rho=0.86)
  # persists throughout the simulation and causes systematic divergence in
  # derived quantities such as undiagnosed prevalence.
  resid <- resid - mean(resid)
  trend_logit + resid
}

#' Generate diagnosis probability patterns
#'
#' @param nquar Number of time quarters.
#' @param nage Number of age groups.
#' @param pattern Pattern type for diagnosis probabilities ("constant", "varying").
#' @param dist Logical indicating whether to add normal variation.
#' @param rita Logical indicating whether to generate an additional diagnosis probability for recently acquired HIV.
#' @param age_diag_effect Numeric scalar controlling the strength of age variation in diagnosis probabilities. A value of 0 means no age effect; positive values increase diagnosis probability with age. Applied as a multiplicative factor on the logit scale.
#' @return A 3D array of diagnosis probabilities with dimensions \code{[nquar, nage, nstates]}.
generate_diagnosis_probs <- function(
  nquar,
  nage = 52,
  pattern,
  dist = TRUE,
  rita = FALSE,
  age_diag_effect = 0
) {
  # check pattern
  stopifnot(pattern %in% c("constant", "varying"))
  if (dist) {
    # Simulate from an AR(1) on the logit scale, matching the Stan AR(1) prior.
    # rho and vardelta are set to the prior mode / tightened prior centre so that
    # the simulated trajectories are drawn from a distribution the model can fit.
    # Stationary SD = sqrt(vardelta / (1 - rho^2)) ~ 0.26 logit units.
    rho_sim <- 0.86       # prior mode of beta(7, 2)
    vardelta_sim <- exp(-4.0)  # tightened prior centre of mu_log_var

    if (pattern == "constant") {
      trend <- matrix(
        c(logit(0.055), logit(0.065), logit(0.13), logit(0.35)),
        nrow = nquar, ncol = 4, byrow = TRUE
      )
      if (rita) trend <- cbind(rep(logit(0.15), nquar), trend)
    } else if (pattern == "varying") {
      start_values <- c(logit(0.040), logit(0.050), logit(0.10), logit(0.28))
      mid_values   <- c(logit(0.055), logit(0.065), logit(0.13), logit(0.35))
      end_values   <- c(logit(0.065), logit(0.080), logit(0.16), logit(0.40))
      trend <- matrix(0, nrow = nquar, ncol = 4)
      for (k in 1:4) {
        trend[, k] <- c(
          seq(start_values[k], mid_values[k], length.out = nquar / 2),
          seq(mid_values[k],   end_values[k], length.out = nquar / 2)
        )
      }
      if (rita) {
        rita_trend <- c(
          seq(logit(0.12), logit(0.18), length.out = nquar / 2),
          seq(logit(0.18), logit(0.18), length.out = nquar / 2)
        )
        trend <- cbind(rita_trend, trend)
      }
    }

    nstates_sim <- ncol(trend)
    logit_d_probs <- matrix(0, nrow = nquar, ncol = nstates_sim)
    for (k in seq_len(nstates_sim)) {
      logit_d_probs[, k] <- simulate_ar1_logit(
        nquar, trend[, k], rho_sim, vardelta_sim
      )
    }
    d_probs <- invlogit(logit_d_probs)
  } else {
    if (pattern == "constant") {
      d_probs <- matrix(
        c(0.055, 0.065, 0.13, 0.35),
        nrow = nquar,
        ncol = 4,
        byrow = TRUE
      )
      if (rita) {
        d_probs <- cbind(
          rep(0.15, nquar),
          d_probs
        )
      }
    } else if (pattern == "varying") {
      # Amplitude kept within ~2x the AR(1) stationary SD (~0.28 logit units).
      start_values <- c(logit(0.040), logit(0.050), logit(0.10), logit(0.28))
      mid_values <- c(logit(0.055), logit(0.065), logit(0.13), logit(0.35))
      end_values <- c(logit(0.065), logit(0.080), logit(0.16), logit(0.40))
      d_probs <- matrix(
        invlogit(c(
          c(
            seq(start_values[1], mid_values[1], length.out = nquar / 2),
            seq(mid_values[1], end_values[1], length.out = nquar / 2)
          ),
          c(
            seq(start_values[2], mid_values[2], length.out = nquar / 2),
            seq(mid_values[2], end_values[2], length.out = nquar / 2)
          ),
          c(
            seq(start_values[3], mid_values[3], length.out = nquar / 2),
            seq(mid_values[3], end_values[3], length.out = nquar / 2)
          ),
          c(
            seq(start_values[4], mid_values[4], length.out = nquar / 2),
            seq(mid_values[4], end_values[4], length.out = nquar / 2)
          )
        )),
        nrow = nquar,
        ncol = 4
      )
      if (rita) {
        d_probs <- cbind(
          c(
            seq(0.12, 0.18, length.out = nquar / 2),
            seq(0.18, 0.18, length.out = nquar / 2)
          ),
          d_probs
        )
      }
    }
  }

  # Expand to 3D array [nquar, nage, nstates] with age-group-level effect.
  # Uses 4 age groups matching the Stan model:
  #   group 1: ages 1-10, group 2: ages 11-20, group 3: ages 21-30, group 4: ages 31+
  nstates <- ncol(d_probs)
  nagroup <- 4
  age_group_lookup <- pmin(nagroup, 1 + (seq_len(nage) - 1) %/% 10)

  # Centred group scale: e.g. for 4 groups -> c(-0.5, -0.167, 0.167, 0.5)
  group_scale <- (seq_len(nagroup) - 1) / (nagroup - 1) - 0.5
  # Logit-additive offset per age group; zero when age_diag_effect == 0
  group_offset <- age_diag_effect * group_scale # length nagroup

  d_probs_3d <- array(0, dim = c(nquar, nage, nstates))
  logit_base <- logit(d_probs) # [nquar, nstates]
  for (ag in seq_len(nagroup)) {
    idx <- which(age_group_lookup == ag)
    shifted <- invlogit(logit_base + group_offset[ag])
    for (a in idx) {
      d_probs_3d[, a, ] <- shifted
    }
  }

  return(d_probs_3d)
}

#' Generate age distributions
#'
#' Generates age distributions for a given number of age categories over the epidemic.
#'
#' @param nquar Number of time quarters.
#' @param nage Number of age groups.
#' @param pattern Pattern type for age distribution ("constant", "wedge", "varying_wedge").
#' @param peak_start Peak age during first time quarter.
#' @param peak_end Peak age during final time quarter.
#' @return Matrix of age distribution over time.
generate_age_dist <- function(
  nquar,
  nage = 52,
  pattern,
  peak_start = 10,
  peak_end = 30
) {
  # check pattern
  stopifnot(pattern %in% c("constant", "wedge", "varying_wedge"))

  if (pattern == "constant") {
    age_dist <- matrix(1, ncol = nquar, nrow = nage)
  } else if (pattern == "wedge") {
    age_dist <- matrix(
      c(
        seq(0.8, 1.4, length.out = ceiling(nage / 8 * 3)),
        seq(1.4, 0.8, length.out = floor(nage / 8 * 5))
      ),
      ncol = nquar,
      nrow = nage
    )
  } else if (pattern == "varying_wedge") {
    # currently doesn't pass the test of preserving totals to diagnoses
    peak_age <- round(seq(peak_start, peak_end, length.out = nquar))
    peak_age <- pmin(pmax(peak_age, 1L), nage - 1L)
    age_dist <- matrix(0, ncol = nquar, nrow = nage)
    for (i in 1:nquar) {
      age_dist[, i] <- c(
        seq(0.8, 1.4, length.out = peak_age[i]),
        seq(1.4, 0.8, length.out = nage - peak_age[i])
      )
    }
  }

  # When nage > 40, compensate for unequal age group widths so that wider groups
  # don't dominate narrower groups
  if (nage > 40) {
    group_sizes <- c(rep(10, 3), nage - 30)
    age_group_idx <- c(rep(1, 10), rep(2, 10), rep(3, 10), rep(4, nage - 30))
    compensation <- 1 / group_sizes[age_group_idx]
    age_dist <- age_dist * compensation
  }

  # normalize age distribution
  age_dist <- age_dist / colSums(age_dist)
  # transpose matrix
  return(t(age_dist))
}

#' Apply false recency rate misclassification to a RITA CD4 array
#'
#' Moves a fraction \code{frr} of long-standing (non-RITA-recent, states 2+)
#' diagnoses to the RITA-recent state (state 1), simulating false recency
#' classification errors in the RITA assay. Works on arrays of any rank whose
#' \emph{last} dimension indexes the RITA/CD4 states.
#'
#' @param cd4_arr Integer array whose last dimension indexes RITA/CD4 states
#' @param frr Scalar false recency rate in \eqn{[0, 1]}
#' @param dist Logical; if \code{TRUE} the number of misclassified cases per
#'   cell is drawn from a Binomial distribution, otherwise the deterministic
#'   expectation (rounded, preserving the total) is used.
#' @return An integer array with the same dimensions as \code{cd4_arr}.
#' @keywords internal
apply_frr <- function(cd4_arr, frr, dist = TRUE) {
  if (frr <= 0) return(cd4_arr)

  dims     <- dim(cd4_arr)
  n_states <- dims[length(dims)]
  n_cells  <- prod(dims[-length(dims)])

  # Reshape to (n_cells x n_states) for vectorised operations
  mat <- matrix(as.integer(cd4_arr), nrow = n_cells, ncol = n_states)

  # Total long-standing diagnoses per cell (states 2:n_states)
  n_long <- rowSums(mat[, 2:n_states, drop = FALSE])

  # Number misclassified as RITA-recent
  if (dist) {
    n_mis <- stats::rbinom(n_cells, n_long, frr)
  } else {
    n_mis <- round(n_long * frr)
  }

  # Subtract proportionally from long-standing states, add to recent (state 1)
  idx <- which(n_mis > 0 & n_long > 0)
  for (i in idx) {
    props              <- mat[i, 2:n_states] / n_long[i]
    reduction          <- round_preserve(n_mis[i] * props)
    mat[i, 2:n_states] <- mat[i, 2:n_states] - reduction
  }
  mat[, 1] <- mat[, 1] + n_mis

  array(mat, dim = dims)
}

#' Derive equilibrium initial state proportions
#'
#' Computes quasi-steady-state latent compartment proportions by solving the
#' linear no-aging system \code{(I - Q) lat = h}, where \code{Q} is the
#' single-age-group transition matrix built from the supplied diagnosis
#' probabilities and progression rates.  The result gives good starting values
#' for \code{init_props} or \code{init_props_abroad} when the undiagnosed stock
#' is assumed to be near equilibrium before the simulation window.
#'
#' The function mirrors the compartment logic in \code{\link{Q_mat}}: for the
#' native group (\code{abroad = FALSE}) all inflow enters state 1; for the
#' migrant stock (\code{abroad = TRUE}) inflow is distributed across latent
#' states proportionally to \code{mu_lat_t1}.
#'
#' @param d_probs_t1 Numeric vector of diagnosis probabilities at t = 1.
#'   Length 4 for CD4-only (\code{rita = FALSE}); length 5 for RITA
#'   (\code{rita = TRUE}): element 1 is the RITA (recent-infection) diagnosis
#'   probability, elements 2–5 are the CD4-stratified probabilities.
#' @param q Numeric vector of length 4: quarterly CD4 progression probabilities,
#'   matching the default \code{q = c(0.0976, 0.115, 0.116, 0.148)}.
#' @param mu_lat_t1 Numeric vector of relative inflow weights across latent
#'   states at t = 1.  Required when \code{abroad = TRUE}.  For RITA, the first
#'   three elements should be 0 (migrants enter long-standing states only).
#'   Will be normalised internally to sum to 1.
#' @param abroad Logical.  \code{FALSE} (default): solve for the native group
#'   with all inflow entering state 1.  \code{TRUE}: solve for the migrant stock
#'   with inflow distributed via \code{mu_lat_t1}.
#' @param scale Numeric or \code{NULL}.  When non-\code{NULL}, the output is
#'   scaled so that its first non-zero element equals \code{scale}.  Defaults to
#'   7 for the native group so the result is on the same scale as the built-in
#'   \code{init_props} default.  Defaults to \code{NULL} for the migrant group:
#'   the unscaled solve already gives values in the correct units for direct use
#'   as \code{init_props_abroad} (latent stock per unit of quarterly migration
#'   into latent states).
#' @param rita Logical.  Use the RITA compartment structure (7 latent states).
#' @return Numeric vector of length \code{nlatst} (4 for CD4-only, 7 for RITA).
#' @keywords internal
derive_init_props <- function(
  d_probs_t1 = c(0.055, 0.065, 0.13, 0.35),
  q          = c(0.0976, 0.115, 0.116, 0.148),
  mu_lat_t1  = NULL,
  abroad     = FALSE,
  scale      = if (!abroad) 7 else NULL,
  rita       = FALSE
) {
  nlatst <- if (rita) 7L else 4L

  if (rita) {
    if (length(d_probs_t1) != 5L) {
      stop("For rita = TRUE, d_probs_t1 must have length 5 (RITA + 4 CD4 states).")
    }
  } else {
    if (length(d_probs_t1) != 4L) {
      stop("For rita = FALSE, d_probs_t1 must have length 4.")
    }
  }
  if (length(q) != 4L) stop("q must have length 4.")
  if (abroad) {
    if (is.null(mu_lat_t1)) stop("mu_lat_t1 is required when abroad = TRUE.")
    if (length(mu_lat_t1) != nlatst) {
      stop(sprintf("mu_lat_t1 must have length %d (= nlatst).", nlatst))
    }
    if (sum(mu_lat_t1) == 0) stop("mu_lat_t1 must have at least one non-zero element.")
  }

  d <- d_probs_t1

  # Build Q matrix (mirrors Q_mat for a single representative age group) --------
  Q <- matrix(0, nlatst, nlatst)

  if (rita) {
    # Recent (RITA-positive) states 1–3: all diagnosed at the RITA rate d[1]
    Q[1L, 2L] <- (1 - d[1L]) * (1 - q[1L])   # state 1 → 2 (CD4 > 500, 3–6 mo)
    Q[1L, 3L] <- (1 - d[1L]) *      q[1L]    # state 1 → 3 (CD4 350–500, 3–6 mo)
    Q[2L, 4L] <- (1 - d[1L]) * (1 - q[1L])   # state 2 → 4 (CD4 > 500, long-standing)
    Q[2L, 5L] <- (1 - d[1L]) *      q[1L]    # state 2 → 5 (CD4 350–500, long-standing)
    Q[3L, 5L] <- (1 - d[1L]) * (1 - q[2L])   # state 3 → 5 (CD4 350–500, long-standing)
    Q[3L, 6L] <- (1 - d[1L]) *      q[2L]    # state 3 → 6 (CD4 200–350, long-standing)
    # Long-standing states 4–7: diagnosed at d[2:5]
    for (k in 4L:7L) {
      dk        <- d[k - 2L]
      qk        <- q[k - 3L]
      Q[k, k]   <- (1 - dk) * (1 - qk)
    }
    for (i in 4L:6L) {
      Q[i, i + 1L] <- (1 - d[i - 2L]) * q[i - 3L]
    }
  } else {
    for (k in 1L:4L) {
      Q[k, k] <- (1 - d[k]) * (1 - q[k])
    }
    for (k in 1L:3L) {
      Q[k, k + 1L] <- (1 - d[k]) * q[k]
    }
  }

  # Inflow vector h --------------------------------------------------------------
  h <- if (abroad) {
    mu_lat_t1 / sum(mu_lat_t1)
  } else {
    c(1, rep(0, nlatst - 1L))
  }

  # Solve (I - t(Q)) lat = h  (steady-state: lat = t(Q) %*% lat + h) ------------
  lat_ss <- solve(diag(nlatst) - t(Q), h)

  # Optional scaling -------------------------------------------------------------
  if (!is.null(scale)) {
    first_nz <- which(lat_ss > 0)[1L]
    lat_ss   <- lat_ss * (scale / lat_ss[first_nz])
  }

  lat_ss
}
