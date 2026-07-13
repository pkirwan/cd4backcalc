#' Calibrate initial prevalence from fitted draws
#'
#' Back-calculates improved initial prevalence values from quarter-1 observed
#' diagnoses and the fitted quarter-1 diagnosis probabilities. The helper works
#' with age-dependent and age-independent models, with or without migration, and
#' for both CD4-only and CD4-plus-RITA variants.
#'
#' The returned values are posterior summaries of baseline latent prevalence that
#' would reproduce the first-quarter diagnoses under the fitted diagnosis model.
#' For RITA models, the first three latent "recent" states are not separately
#' identified from quarter-1 data; by default they are split using the relative
#' weights in the currently supplied initial prevalence.
#'
#' @param object Either a list returned by [run_backcalc()] containing `fit` and
#'   `stan_data`, or a fitted model object with a `draws(format = "matrix")`
#'   method.
#' @param stan_data Optional Stan data list. Supply this when `object` is a fit
#'   rather than the full list returned by [run_backcalc()].
#' @param summary Posterior summary to return. One of `"median"` or `"mean"`.
#' @param use_hiv_total Logical; if `TRUE`, quarter-1 HIV totals are allocated
#'   across CD4/RITA categories using the observed category proportions. This is
#'   useful when CD4 completeness is below 100%.
#' @param smooth_age Logical; if `TRUE`, age-specific summaries are smoothed over
#'   age with loess after summarising across draws.
#' @param smooth_span Span passed to [stats::loess()] when `smooth_age = TRUE`.
#' @param recent_split How to split the unidentified recent latent pool in RITA
#'   models. `"current"` preserves the relative weights in the current
#'   `init_prev` input; `"equal"` splits it equally across states 1:3.
#' @param min_prob Lower bound used when inverting fitted diagnosis
#'   probabilities.
#' @param init_floor Small positive floor applied to calibrated initial
#'   prevalence values to avoid non-finite or tiny positive values that can
#'   cause Stan initialisation failures; exact structural zeros are preserved.
#' @param return_draws Logical; if `TRUE`, include draw-level calibrated values
#'   in the output.
#'
#' @return A list containing calibrated initial prevalence values, an updated
#'   `stan_data` list ready for refitting, and optional draw-level arrays.
#' @examples
#' # See ?run_backcalc for a full end-to-end workflow.
#' invisible(NULL)
#' @export
calibrate_init_prev <- function(
  object,
  stan_data = NULL,
  summary = c("median", "mean"),
  use_hiv_total = TRUE,
  smooth_age = TRUE,
  smooth_span = 0.5,
  recent_split = c("current", "equal"),
  min_prob = 1e-6,
  init_floor = 1e-8,
  return_draws = FALSE
) {
  summary <- match.arg(summary)
  recent_split <- match.arg(recent_split)

  resolved <- .resolve_init_prev_inputs(object, stan_data)
  fit <- resolved$fit
  stan_data <- resolved$stan_data

  if (is.null(stan_data$model_config)) {
    stop("stan_data$model_config is required.")
  }

  config <- stan_data$model_config
  fit_mat <- fit$draws(format = "matrix")
  n_draw <- nrow(fit_mat)
  n_latent <- if (isTRUE(config$rita)) 7L else 4L

  if (isTRUE(config$age)) {
    native_draws <- array(0, dim = c(n_draw, n_latent, stan_data$nage))
    migrant_draws <- if (isTRUE(config$migration)) {
      array(0, dim = c(n_draw, n_latent, stan_data$nage))
    } else {
      NULL
    }

    ad_idx <- .make_ad_t1_indices(fit_mat, config)
    mis_idx <- if (isTRUE(config$migration)) {
      .make_ad_misclassified_indices(fit_mat, stan_data$nage)
    } else {
      NULL
    }

    for (i in seq_len(n_draw)) {
      frr_i <- .extract_frr_draw(fit_mat, i, config)
      m_lat_i <- .extract_m_lat_t1(fit_mat, i, config)

      for (a in seq_len(stan_data$nage)) {
        ag <- .age_to_group(a)
        q4 <- stan_data$q_age[4, a]

        native_hiv <- if (isTRUE(config$migration)) {
          stan_data$HIV_age[1, a, 1]
        } else {
          stan_data$HIV_age[1, a]
        }
        native_aids <- if (isTRUE(config$migration)) {
          stan_data$AIDS_age[1, a, 1]
        } else {
          stan_data$AIDS_age[1, a]
        }
        native_cd4 <- if (isTRUE(config$migration)) {
          as.numeric(stan_data$CD4_age[1, a, 1, ])
        } else {
          as.numeric(stan_data$CD4_age[1, a, ])
        }
        native_current <- stan_data$init_prev_age[, a]
        native_d <- as.numeric(fit_mat[i, ad_idx$native[ag, ], drop = TRUE])
        native_diag <- .quarter1_diag_counts(
          hiv_total = native_hiv,
          cd4_counts = native_cd4,
          current_init = native_current,
          rita = config$rita,
          frr = frr_i,
          use_hiv_total = use_hiv_total
        )
        native_draws[i, , a] <- .solve_t1_latent(
          diag_counts = native_diag,
          aids_total = native_aids,
          d = native_d,
          q4 = q4,
          rita = config$rita,
          frr = frr_i,
          current_init = native_current,
          recent_split = recent_split,
          min_prob = min_prob
        )

        if (isTRUE(config$migration)) {
          migrant_hiv <- stan_data$HIV_age[1, a, 2]
          migrant_aids <- stan_data$AIDS_age[1, a, 2]
          migrant_cd4 <- as.numeric(stan_data$CD4_age[1, a, 2, ])
          migrant_current <- stan_data$init_prev_abroad_age[, a]
          migrant_d <- as.numeric(fit_mat[i, ad_idx$migrant[ag, ], drop = TRUE])
          migrant_diag <- .quarter1_diag_counts(
            hiv_total = migrant_hiv,
            cd4_counts = migrant_cd4,
            current_init = migrant_current,
            rita = config$rita,
            frr = frr_i,
            use_hiv_total = use_hiv_total
          )
          mis_total <- fit_mat[i, mis_idx[a]]
          migrant_diag <- pmax(
            migrant_diag - .misclassified_observed_counts(mis_total, m_lat_i, config$rita, frr_i),
            0
          )
          migrant_draws[i, , a] <- .solve_t1_latent(
            diag_counts = migrant_diag,
            aids_total = migrant_aids,
            d = migrant_d,
            q4 = q4,
            rita = config$rita,
            frr = frr_i,
            current_init = migrant_current,
            recent_split = recent_split,
            min_prob = min_prob
          )
        }
      }
    }

    init_prev_age <- .summarize_init_draws(native_draws, summary)
    if (isTRUE(smooth_age)) {
      init_prev_age <- .smooth_init_matrix(init_prev_age, smooth_span)
    }
    init_prev_age <- .sanitize_init_values(init_prev_age, floor = init_floor)

    out <- list(
      init_prev_age = init_prev_age,
      stan_data = utils::modifyList(stan_data, list(init_prev_age = init_prev_age))
    )

    if (isTRUE(config$migration)) {
      init_prev_abroad_age <- .summarize_init_draws(migrant_draws, summary)
      if (isTRUE(smooth_age)) {
        init_prev_abroad_age <- .smooth_init_matrix(init_prev_abroad_age, smooth_span)
      }
      init_prev_abroad_age <- .sanitize_init_values(init_prev_abroad_age, floor = init_floor)
      out$init_prev_abroad_age <- init_prev_abroad_age
      out$stan_data$init_prev_abroad_age <- init_prev_abroad_age
    }

    if (isTRUE(return_draws)) {
      out$draws <- list(init_prev_age = native_draws)
      if (isTRUE(config$migration)) {
        out$draws$init_prev_abroad_age <- migrant_draws
      }
    }

    return(out)
  }

  native_draws <- matrix(0, nrow = n_draw, ncol = n_latent)
  migrant_draws <- if (isTRUE(config$migration)) {
    matrix(0, nrow = n_draw, ncol = n_latent)
  } else {
    NULL
  }

  ai_idx <- .make_ai_t1_indices(fit_mat, config)

  for (i in seq_len(n_draw)) {
    frr_i <- .extract_frr_draw(fit_mat, i, config)
    m_lat_i <- .extract_m_lat_t1(fit_mat, i, config)

    native_current <- stan_data$init_prev
    native_diag <- .quarter1_diag_counts(
      hiv_total = if (isTRUE(config$migration)) stan_data$HIV[1, 1] else stan_data$HIV[1],
      cd4_counts = if (isTRUE(config$migration)) as.numeric(stan_data$CD4[1, 1, ]) else as.numeric(stan_data$CD4[1, ]),
      current_init = native_current,
      rita = config$rita,
      frr = frr_i,
      use_hiv_total = use_hiv_total
    )
    native_draws[i, ] <- .solve_t1_latent(
      diag_counts = native_diag,
      aids_total = if (isTRUE(config$migration)) stan_data$AIDS[1, 1] else stan_data$AIDS[1],
      d = as.numeric(fit_mat[i, ai_idx$native, drop = TRUE]),
      q4 = stan_data$q[4],
      rita = config$rita,
      frr = frr_i,
      current_init = native_current,
      recent_split = recent_split,
      min_prob = min_prob
    )

    if (isTRUE(config$migration)) {
      migrant_current <- stan_data$init_prev_abroad
      migrant_diag <- .quarter1_diag_counts(
        hiv_total = stan_data$HIV[1, 2],
        cd4_counts = as.numeric(stan_data$CD4[1, 2, ]),
        current_init = migrant_current,
        rita = config$rita,
        frr = frr_i,
        use_hiv_total = use_hiv_total
      )
      migrant_diag <- pmax(
        migrant_diag - .misclassified_observed_counts(fit_mat[i, ai_idx$mis], m_lat_i, config$rita, frr_i),
        0
      )
      migrant_draws[i, ] <- .solve_t1_latent(
        diag_counts = migrant_diag,
        aids_total = stan_data$AIDS[1, 2],
        d = as.numeric(fit_mat[i, ai_idx$migrant, drop = TRUE]),
        q4 = stan_data$q[4],
        rita = config$rita,
        frr = frr_i,
        current_init = migrant_current,
        recent_split = recent_split,
        min_prob = min_prob
      )
    }
  }

  init_prev <- as.numeric(.summarize_init_draws(native_draws, summary))
  init_prev <- as.numeric(.sanitize_init_values(init_prev, floor = init_floor))
  out <- list(
    init_prev = init_prev,
    stan_data = utils::modifyList(stan_data, list(init_prev = init_prev))
  )

  if (isTRUE(config$migration)) {
    init_prev_abroad <- as.numeric(.summarize_init_draws(migrant_draws, summary))
    init_prev_abroad <- as.numeric(.sanitize_init_values(init_prev_abroad, floor = init_floor))
    out$init_prev_abroad <- init_prev_abroad
    out$stan_data$init_prev_abroad <- init_prev_abroad
  }

  if (isTRUE(return_draws)) {
    out$draws <- list(init_prev = native_draws)
    if (isTRUE(config$migration)) {
      out$draws$init_prev_abroad <- migrant_draws
    }
  }

  out
}

.resolve_init_prev_inputs <- function(object, stan_data) {
  if (is.list(object) && !is.null(object$fit) && !is.null(object$stan_data)) {
    return(list(fit = object$fit, stan_data = object$stan_data))
  }

  if (is.null(stan_data)) {
    stop("Supply either a run_backcalc() result or both fit and stan_data.")
  }

  list(fit = object, stan_data = stan_data)
}

.make_ai_t1_indices <- function(fit_mat, config) {
  cn <- colnames(fit_mat)
  diag_states <- if (isTRUE(config$rita)) 5L else 4L
  native_prefix <- if (isTRUE(config$migration)) "d_nat" else "d"
  native <- vapply(
    seq_len(diag_states),
    function(k) match(sprintf("%s[%d,1]", native_prefix, k), cn),
    integer(1)
  )

  if (anyNA(native)) {
    stop("Could not find quarter-1 native diagnosis probabilities in the fit draws.")
  }

  out <- list(native = native)

  if (isTRUE(config$migration)) {
    migrant <- vapply(
      seq_len(diag_states),
      function(k) match(sprintf("d_mig[%d,1]", k), cn),
      integer(1)
    )
    if (anyNA(migrant)) {
      stop("Could not find quarter-1 migrant diagnosis probabilities in the fit draws.")
    }
    mis_idx <- match("misclassified_dx[1,2]", cn)
    if (is.na(mis_idx)) {
      stop("Could not find quarter-1 misclassified_dx[1,2] in the fit draws.")
    }
    out$migrant <- migrant
    out$mis <- mis_idx
  }

  out
}

.make_ad_t1_indices <- function(fit_mat, config) {
  cn <- colnames(fit_mat)
  diag_states <- if (isTRUE(config$rita)) 5L else 4L
  native_name <- if (isTRUE(config$migration)) "d_nat" else "d"
  native <- matrix(NA_integer_, nrow = 4, ncol = diag_states)

  for (ag in 1:4) {
    for (k in seq_len(diag_states)) {
      native[ag, k] <- match(sprintf("%s[%d,1,%d]", native_name, ag, k), cn)
    }
  }

  if (anyNA(native)) {
    stop("Could not find quarter-1 age-dependent diagnosis probabilities in the fit draws.")
  }

  out <- list(native = native)
  if (isTRUE(config$migration)) {
    migrant <- matrix(NA_integer_, nrow = 4, ncol = diag_states)
    for (ag in 1:4) {
      for (k in seq_len(diag_states)) {
        migrant[ag, k] <- match(sprintf("d_mig[%d,1,%d]", ag, k), cn)
      }
    }
    if (anyNA(migrant)) {
      stop("Could not find quarter-1 age-dependent migrant diagnosis probabilities in the fit draws.")
    }
    out$migrant <- migrant
  }

  out
}

.make_ad_misclassified_indices <- function(fit_mat, nage) {
  cn <- colnames(fit_mat)
  idx <- vapply(
    seq_len(nage),
    function(a) match(sprintf("misclassified_dx[1,%d]", a), cn),
    integer(1)
  )
  if (anyNA(idx)) {
    stop("Could not find quarter-1 age-specific misclassified_dx values in the fit draws.")
  }
  idx
}

.extract_frr_draw <- function(fit_mat, draw_id, config) {
  if (!isTRUE(config$rita)) {
    return(0)
  }

  idx <- match("frr", colnames(fit_mat))
  if (is.na(idx)) {
    stop("Could not find frr in the fit draws for a RITA model.")
  }

  as.numeric(fit_mat[draw_id, idx])
}

.extract_m_lat_t1 <- function(fit_mat, draw_id, config) {
  if (!isTRUE(config$migration)) {
    return(NULL)
  }

  cn <- colnames(fit_mat)
  idx <- vapply(
    1:4,
    function(k) match(sprintf("m[%d,1]", k), cn),
    integer(1)
  )
  if (anyNA(idx)) {
    stop("Could not find m[1:4,1] in the fit draws.")
  }

  vals <- as.numeric(fit_mat[draw_id, idx])
  total <- sum(vals)
  if (total <= 0) {
    rep(0.25, 4)
  } else {
    vals / total
  }
}

.quarter1_diag_counts <- function(
  hiv_total,
  cd4_counts,
  current_init,
  rita,
  frr,
  use_hiv_total
) {
  cd4_counts <- as.numeric(cd4_counts)
  hiv_total <- as.numeric(hiv_total)
  total_cd4 <- sum(cd4_counts)

  if (hiv_total <= 0 && total_cd4 <= 0) {
    return(rep(0, length(cd4_counts)))
  }

  if (total_cd4 > 0) {
    weights <- cd4_counts / total_cd4
  } else {
    weights <- .default_observed_weights(current_init, rita, frr)
  }

  if (!isTRUE(use_hiv_total) && total_cd4 > 0) {
    return(cd4_counts)
  }

  hiv_total * weights
}

.default_observed_weights <- function(current_init, rita, frr) {
  current_init <- pmax(as.numeric(current_init), 0)

  if (!isTRUE(rita)) {
    raw <- current_init[1:4]
  } else {
    raw <- c(sum(current_init[1:3]), current_init[4:7])
    raw <- c(
      raw[1] + frr * sum(raw[2:5]),
      (1 - frr) * raw[2:5]
    )
  }

  total <- sum(raw)
  if (total <= 0) {
    rep(1 / length(raw), length(raw))
  } else {
    raw / total
  }
}

.misclassified_observed_counts <- function(mis_total, m_lat, rita, frr) {
  mis_total <- as.numeric(mis_total)
  if (mis_total <= 0) {
    return(rep(0, if (isTRUE(rita)) 5 else 4))
  }

  if (!isTRUE(rita)) {
    return(mis_total * m_lat)
  }

  c(mis_total * frr, mis_total * (1 - frr) * m_lat)
}

.solve_t1_latent <- function(
  diag_counts,
  aids_total,
  d,
  q4,
  rita,
  frr,
  current_init,
  recent_split,
  min_prob
) {
  diag_counts <- pmax(as.numeric(diag_counts), 0)
  aids_total <- max(as.numeric(aids_total), 0)
  d <- pmin(pmax(as.numeric(d), min_prob), 1 - min_prob)
  q4 <- pmin(pmax(as.numeric(q4), 0), 1)

  if (!isTRUE(rita)) {
    out <- numeric(4)
    out[1:3] <- diag_counts[1:3] / d[1:3]
    out[4] <- (diag_counts[4] + aids_total) / (d[4] + (1 - d[4]) * q4)
    return(out)
  }

  out <- numeric(7)
  nonrecent_scale <- max(1 - frr, min_prob)
  nonrecent_true <- diag_counts[2:5] / nonrecent_scale
  hiv_true <- sum(diag_counts)
  recent_true <- max(hiv_true - sum(nonrecent_true), 0)
  recent_total <- recent_true / d[1]
  recent_weights <- .recent_split_weights(current_init[1:3], recent_split)

  out[1:3] <- recent_total * recent_weights
  out[4] <- nonrecent_true[1] / d[2]
  out[5] <- nonrecent_true[2] / d[3]
  out[6] <- nonrecent_true[3] / d[4]
  out[7] <- (nonrecent_true[4] + aids_total) / (d[5] + (1 - d[5]) * q4)
  out
}

.recent_split_weights <- function(current_recent, recent_split) {
  if (identical(recent_split, "equal")) {
    return(rep(1 / 3, 3))
  }

  current_recent <- pmax(as.numeric(current_recent), 0)
  total <- sum(current_recent)
  if (total <= 0) {
    rep(1 / 3, 3)
  } else {
    current_recent / total
  }
}

.summarize_init_draws <- function(draws, summary) {
  fun <- if (identical(summary, "mean")) mean else stats::median
  if (length(dim(draws)) == 2L) {
    apply(draws, 2, fun)
  } else {
    apply(draws, c(2, 3), fun)
  }
}

.sanitize_init_values <- function(x, floor) {
  floor <- as.numeric(floor)
  if (!is.finite(floor) || floor <= 0) {
    stop("init_floor must be a finite positive number.")
  }

  vals <- as.numeric(x)
  zero_mask <- vals == 0
  vals[!is.finite(vals)] <- floor
  vals[vals < 0] <- floor
  vals[vals > 0 & vals < floor] <- floor
  vals[zero_mask] <- 0
  dim(vals) <- dim(x)
  dimnames(vals) <- dimnames(x)
  vals
}

.smooth_init_matrix <- function(init_mat, span) {
  init_mat <- as.matrix(init_mat)
  if (ncol(init_mat) < 4) {
    return(init_mat)
  }

  age_idx <- seq_len(ncol(init_mat))
  out <- init_mat
  for (k in seq_len(nrow(init_mat))) {
    y <- pmax(as.numeric(init_mat[k, ]), 0)
    if (sum(y > 0) < 3) {
      next
    }
    fit <- try(
      stats::loess(log1p(y) ~ age_idx, span = span, degree = 1),
      silent = TRUE
    )
    if (inherits(fit, "try-error")) {
      next
    }
    pred <- suppressWarnings(stats::predict(fit, newdata = data.frame(age_idx = age_idx)))
    pred <- ifelse(is.finite(pred), expm1(pred), y)
    pred <- pmax(pred, 0)
    if (sum(pred) > 0 && sum(y) > 0) {
      pred <- pred * sum(y) / sum(pred)
    }
    out[k, ] <- pred
  }
  out
}

.age_to_group <- function(age_idx) {
  min(4L, 1L + (age_idx - 1L) %/% 10L)
}