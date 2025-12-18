
#' Generate infections according to specified pattern
#'
#' @param nquar Number of time quarters.
#' @param nlatst Number of latent states.
#' @param pattern Pattern type for infections ("varying", "constant", "increasing").
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
  stopifnot(pattern %in% c("varying", "constant", "increasing"))

  ann_inf <- target_inf / nquar
  infs <- dplyr::case_when(
    pattern == "varying" ~ c(
      seq(1.2, 2, length.out = round(nquar * 0.1)),
      seq(2, 1.4, length.out = round(nquar * 4 / 15)),
      seq(1.4, 1, length.out = round(nquar * 0.2)),
      seq(1, 0.3, length.out = round(nquar * 2 / 15)),
      seq(0.3, 0.1, length.out = round(nquar * 0.3))
    ) *
      ann_inf,
    pattern == "constant" ~ rep(ann_inf, nquar),
    pattern == "increasing" ~ seq(0.5, 1.5, length.out = nquar) * ann_inf
  )
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
#' @param pattern Pattern type for omega ("constant", "increasing", "step_change", "none").
#' @param target_mig Target total number of migrations to simulate (may be higher/lower with Poisson variation).
#' @param dist Logical indicating whether to add Poisson variation.
#' @return Vector of omega values over time.
generate_omega <- function(nquar, pattern, target_mig = 12000, dist = TRUE) {
  # check pattern
  stopifnot(pattern %in% c("constant", "increasing", "step_change", "none"))
  ann_mig = target_mig / nquar
  omega <- dplyr::case_when(
    pattern == "constant" ~ rep(ann_mig, nquar),
    pattern == "increasing" ~ seq(0.5, 1.5, length.out = nquar) * ann_mig,
    pattern == "step_change" ~ c(
      seq(0.4, 0.4, length.out = round(nquar * 0.1)),
      seq(0.4, 1.2, length.out = round(nquar * 0.1)),
      seq(1.2, 0.4, length.out = round(nquar * 0.1)),
      seq(0.4, 1.3, length.out = round(nquar * 0.6)),
      seq(1.3, 4.5, length.out = round(nquar * 0.1))
    ) *
      ann_mig,
    pattern == "none" ~ rep(0, nquar)
  )
  omega <- if (dist) stats::rpois(nquar, omega) else omega
  return(omega)
}

#' Generate migration probability patterns
#'
#' @param nquar Number of time quarters.
#' @param nlatst Number of latent states.
#' @param pattern Pattern type for mu ("constant", "differing", "varying").
#' @return Matrix of mu values over time.
generate_mu <- function(nquar, nlatst, pattern) {
  # check pattern
  stopifnot(pattern %in% c("constant", "differing", "varying"))
  if (pattern == "constant") {
    mu <- matrix(0.2, nrow = nquar, ncol = nlatst + 1)
  } else {
    mu <- dplyr::case_when(
      pattern == "differing" ~ matrix(
        c(0.1, 0.2, 0.2, 0.1, 0.4),
        nrow = nquar,
        ncol = nlatst + 1,
        byrow = TRUE
      ),
      pattern == "varying" ~ matrix(
        c(
          seq(0.1, 0.3, length.out = nquar),
          seq(0.2, 0.22, length.out = nquar),
          seq(0.25, 0.15, length.out = nquar),
          seq(0.3, 0.05, length.out = nquar),
          seq(0.2, 0.4, length.out = nquar)
        ),
        nrow = nquar,
        ncol = nlatst + 1
      )
    )
  }
  # normalize mu
  mu <- mu / rowSums(mu)
  return(mu)
}

#' Generate proportion UK patterns
#'
#' @param nquar Number of time quarters.
#' @param pattern Pattern type for proportion UK ("constant", "increasing").
#' @param dist Logical indicating whether to add normal variation.
#' @return Vector of proportion UK values over time.
generate_prop_uk <- function(nquar, pattern, dist = TRUE) {
  # check pattern
  stopifnot(pattern %in% c("constant", "increasing"))
  if (dist) {
    sigma <- 0.08
    if (pattern == "constant") {
      log_prop_uk <- stats::rnorm(nquar, rep(log(0.2), nquar), sigma)
    } else if (pattern == "increasing") {
      log_prop_uk <- stats::rnorm(
        nquar,
        seq(log(0.1), log(0.4), length.out = nquar),
        sigma
      )
    }
    # Transform back to proportion scale
    prop_uk <- exp(log_prop_uk)
  } else if (pattern == "constant") {
    prop_uk <- rep(0.2, nquar)
  } else if (pattern == "increasing") {
    prop_uk <- exp(seq(log(0.1), log(0.4), length.out = nquar))
  }
  return(prop_uk)
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

#' Generate diagnosis probability patterns
#'
#' @param nquar Number of time quarters.
#' @param pattern Pattern type for diagnosis probabilities ("constant", "varying").
#' @param dist Logical indicating whether to add normal variation.
#' @param rita Logical indicating whether to generate an additional diagnosis probability for recently acquired HIV.
#' @return Matrix of diagnosis probabilities over time.
generate_diagnosis_probs <- function(
  nquar,
  pattern,
  dist = TRUE,
  rita = FALSE
) {
  # check pattern
  stopifnot(pattern %in% c("constant", "varying"))
  if (dist) {
    logit_d_probs <- matrix(0, nrow = nquar, ncol = 4)
    sigma <- 0.1

    if (pattern == "constant") {
      logit_d_probs[, 1] <- stats::rnorm(nquar, logit(0.05), sigma)
      logit_d_probs[, 2] <- stats::rnorm(nquar, logit(0.08), sigma)
      logit_d_probs[, 3] <- stats::rnorm(nquar, logit(0.11), sigma)
      logit_d_probs[, 4] <- stats::rnorm(nquar, logit(0.25), sigma)
      if (rita) {
        logit_d_probs <- cbind(
          stats::rnorm(nquar, logit(0.15), sigma),
          logit_d_probs
        )
      }
    } else if (pattern == "varying") {
      # Initialize first value
      start_values <- c(logit(0.02), logit(0.05), logit(0.10), logit(0.25))
      mid_values <- c(logit(0.03), logit(0.04), logit(0.07), logit(0.35))
      end_values <- c(logit(0.10), logit(0.08), logit(0.07), logit(0.3))

      logit_d_probs[, 1] <- stats::rnorm(
        nquar,
        c(
          seq(start_values[1], mid_values[1], length.out = nquar / 2),
          seq(mid_values[1], end_values[1], length.out = nquar / 2)
        ),
        sigma
      )
      logit_d_probs[, 2] <- stats::rnorm(
        nquar,
        c(
          seq(start_values[2], mid_values[2], length.out = nquar / 2),
          seq(mid_values[2], end_values[2], length.out = nquar / 2)
        ),
        sigma
      )
      logit_d_probs[, 3] <- stats::rnorm(
        nquar,
        c(
          seq(start_values[3], mid_values[3], length.out = nquar / 2),
          seq(mid_values[3], end_values[3], length.out = nquar / 2)
        ),
        sigma
      )
      logit_d_probs[, 4] <- stats::rnorm(
        nquar,
        c(
          seq(start_values[4], mid_values[4], length.out = nquar / 2),
          seq(mid_values[4], end_values[4], length.out = nquar / 2)
        ),
        sigma
      )
      if (rita) {
        logit_d_probs <- cbind(
          seq(logit(0.1), logit(0.2), length.out = nquar),
          logit_d_probs
        )
      }
    }
    d_probs <- invlogit(logit_d_probs)
  } else {
    if (pattern == "constant") {
      d_probs <- matrix(
        c(0.05, 0.08, 0.11, 0.25),
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
      # Initialize first value
      start_values <- c(logit(0.02), logit(0.05), logit(0.10), logit(0.25))
      mid_values <- c(logit(0.03), logit(0.04), logit(0.07), logit(0.35))
      end_values <- c(logit(0.10), logit(0.08), logit(0.07), logit(0.3))
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
          seq(0.1, 0.2, length.out = nquar),
          d_probs
        )
      }
    }
  }
  return(d_probs)
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
        seq(0.2, 1.8, length.out = nage / 2),
        seq(1.8, 0.2, length.out = nage / 2)
      ),
      ncol = nquar,
      nrow = nage
    )
  } else if (pattern == "varying_wedge") {
    peak_age = stats::rpois(
      nquar,
      (seq(peak_start, peak_end, length.out = nquar))
    )
    age_dist <- matrix(0, ncol = nquar, nrow = nage)
    for (i in 1:nquar) {
      age_dist[, i] <- c(
        seq(0.2, 1.8, length.out = peak_age[i]),
        seq(1.8, 0.2, length.out = nage - peak_age[i])
      )
    }
  }
  # normalize age distribution
  age_dist <- age_dist / colSums(age_dist)
  # transpose matrix
  return(t(age_dist))
}
