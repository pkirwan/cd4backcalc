#' Post-processing model fit
#'
#' This function decides which post-processing routine needs to be run
#' and passes the relevant data and arguments.
#'
#' @param fit Fitted stan model object and stan_data list.
#' @param stan_data Data list used in stan model fitting.
#' @param start_yr Starting year for posterior predictive quantities.
#' @param nworkers Number of parallel workers to use.
#' @param chain_results Logical indicating whether to return results by chain.
#' @return List of posterior predictive quantities.
#' @examples
#' \donttest{
#' sim_diags <- simulate_diagnoses(sim_seed = 123)
#' stan_data <- cd4backcalc:::model_type_prep(sim_diags, model_seed = 123)
#' fit <- cd4backcalc:::fit_stan_model(
#'    stan_data,
#'    iter_warmup = 10,
#'    iter_sampling = 10,
#'    chains = 1,
#'    refresh = 0,
#'    messages = FALSE
#' )
#' hiv_list <- cd4backcalc:::postproc(fit, stan_data, start_yr = 0, nworkers = 1)
#' }
#' @keywords internal
postproc <- function(
  fit,
  stan_data,
  start_yr = 0,
  nworkers = 1,
  chain_results = FALSE
) {
  # Run the appropriate post-processing routine
  if (stan_data$model_config$age) {
    postproc_ad(fit, stan_data, start_yr, nworkers, chain_results)
  } else {
    postproc_ai(fit, stan_data, start_yr, nworkers, chain_results)
  }
}

#' Post-process age-independent model fit
#'
#' @param fit Fitted stan model object.
#' @param stan_data Data list used in stan model fitting.
#' @param start_yr Starting year for posterior predictive quantities.
#' @param nworkers Number of parallel workers to use.
#' @param chain_results Logical indicating whether to return results by chain.
#' @return List of posterior predictive quantities.
#' @keywords internal
postproc_ai <- function(
  fit,
  stan_data,
  start_yr = 1995,
  nworkers = 1,
  chain_results = FALSE
) {
  # Use multisession on Windows (multicore is not supported)
  if (.Platform$OS.type == "windows") {
    future::plan(future::multisession, workers = nworkers)
  } else {
    future::plan(future::multicore, workers = nworkers)
  }
  on.exit(future::plan(future::sequential), add = TRUE)

  # increase memory limit
  old_future_max_size <- getOption("future.globals.maxSize")
  options(future.globals.maxSize = 3000 * 1024^2)
  on.exit(options(future.globals.maxSize = old_future_max_size), add = TRUE)

  # generate matrix object
  fit_mat <- fit$draws(format = "matrix")

  # Precompute column indices once
  idx <- make_fit_indices(fit_mat, stan_data$model_config, stan_data$nquar)

  if (chain_results) {
    iter_per_chain <- nrow(fit_mat) / fit$num_chains()
    draw_splits <- rep(1:fit$num_chains(), each = iter_per_chain)

    hiv_list <- list()
    for (ch in 1:fit$num_chains()) {
      draw_ids <- which(draw_splits == ch)

      inf_dist <- future.apply::future_sapply(
        draw_ids,
        function(ii) infect_dist_by_year_dx(stan_data, fit_mat, ii, idx),
        future.seed = TRUE,
        future.globals = c("stan_data", "fit_mat", "idx")
      )
      prev_dist <- future.apply::future_sapply(
        draw_ids,
        function(ii) prev_by_year_infection(stan_data, fit_mat, ii, idx),
        future.seed = TRUE,
        future.globals = c("stan_data", "fit_mat", "idx")
      )

      ppq <- build_ppq_list(
        stan_data = stan_data,
        prev_dist = prev_dist,
        inf_dist = inf_dist,
        fit_slice = fit_mat[draw_ids, , drop = FALSE],
        start_yr = start_yr
      )

      # Store results by chain
      hiv_list[[paste0("incidence_", ch)]] <- ppq$incidence
      hiv_list[[paste0("incidence_annual_", ch)]] <- ppq$incidence_annual
      hiv_list[[paste0("undiag_prev_", ch)]] <- ppq$undiag_prev
      hiv_list[[paste0("diag_prob_", ch)]] <- ppq$diag_prob
      hiv_list[[paste0("diagnoses_", ch)]] <- ppq$diagnoses

      if (stan_data$model_config$migration) {
        hiv_list[[paste0("all_migration_", ch)]] <- ppq$all_migration
        hiv_list[[paste0("undiag_migration_", ch)]] <- ppq$undiag_migration
        hiv_list[[paste0("migration_prob_", ch)]] <- ppq$migration_prob
        hiv_list[[paste0("transfers_", ch)]] <- ppq$transfers
        hiv_list[[paste0("transfers_detected_", ch)]] <- ppq$transfers_detected
        hiv_list[[paste0("ratio_abroad_uk_", ch)]] <- ppq$ratio_abroad_uk
        hiv_list[[paste0("detect_prob_", ch)]] <- ppq$detect_prob
      }
    }
  } else {
    draw_ids <- seq_len(nrow(fit_mat))

    inf_dist <- future.apply::future_sapply(
      draw_ids,
      function(ii) infect_dist_by_year_dx(stan_data, fit_mat, ii, idx),
      future.seed = TRUE,
      future.globals = c("stan_data", "fit_mat", "idx")
    )
    prev_dist <- future.apply::future_sapply(
      draw_ids,
      function(ii) prev_by_year_infection(stan_data, fit_mat, ii, idx),
      future.seed = TRUE,
      future.globals = c("stan_data", "fit_mat", "idx")
    )

    hiv_list <- build_ppq_list(
      stan_data = stan_data,
      prev_dist = prev_dist,
      inf_dist = inf_dist,
      fit_slice = fit_mat,
      start_yr = start_yr
    )
  }

  diagnostics <- fit$diagnostic_summary()
  hiv_list$conv.chk <- diagnostics$ebfmi > 0.3
  hiv_list$n.divg <- sum(diagnostics$num_divergent)
  hiv_list$time <- fit$time()$total

  return(hiv_list)
}

#' Make fit indices
#'
#' Creates a list of indices for extracting parameters from the fit matrix.
#'
#' @param fit_mat Matrix of posterior samples from stan model.
#' @param config Configuration list used in stan model fitting.
#' @param nquar Number of quarters in the model.
#' @return List of indices for various parameters in the fit matrix.
make_fit_indices <- function(fit_mat, config, nquar) {
  cn <- colnames(fit_mat)

  # Migration models expose d_nat[k,t] and d_mig[k,t] instead of d[k,t].
  # Use d_nat as the representative UK-born diagnosis probabilities.
  dp <- if (isTRUE(config$migration)) "d_nat" else "d"

  idx <- list(
    d1 = grep(sprintf("^%s\\[1,", dp), cn),
    d2 = grep(sprintf("^%s\\[2,", dp), cn),
    d3 = grep(sprintf("^%s\\[3,", dp), cn),
    d4 = grep(sprintf("^%s\\[4,", dp), cn),
    d5 = grep(sprintf("^%s\\[5,", dp), cn),
    d1_mig = grep("^d_mig\\[1,", cn),
    d2_mig = grep("^d_mig\\[2,", cn),
    d3_mig = grep("^d_mig\\[3,", cn),
    d4_mig = grep("^d_mig\\[4,", cn),
    d5_mig = grep("^d_mig\\[5,", cn),
    lambda = grep("^lambda\\[", cn),
    gamma = grep("^gamma\\[", cn),
    p = grep("^p\\[", cn),
    m = grep("^m\\[", cn),
    mis_dx = grep("^misclassified_dx\\[[0-9]+,2\\]$", cn),
    frr = match("frr", cn),
    under_rep = match("under_rep", cn)
  )

  if (
    any(c(length(idx$d1), length(idx$d2), length(idx$d3), length(idx$d4)) == 0L)
  ) {
    stop("Could not find d[1..4,] in fit matrix column names.")
  }
  if (config$rita && length(idx$d5) == 0L) {
    stop(
      "config$rita=TRUE but could not find d[5,] in fit matrix column names."
    )
  }
  if (length(idx$lambda) == 0L) {
    stop("Could not find lambda[...] in fit matrix column names.")
  }

  if (config$migration) {
    if (length(idx$gamma) == 0L) {
      stop("config$migration=TRUE but gamma[...] not found.")
    }
    if (length(idx$p) == 0L) {
      stop("config$migration=TRUE but p[...] not found.")
    }
    if (length(idx$m) == 0L) stop("config$migration=TRUE but m[...] not found.")
    if (length(idx$mis_dx) != nquar) {
      stop("config$migration=TRUE but misclassified_dx[*,2] not found or wrong length.")
    }
    if (config$rita && is.na(idx$frr)) {
      stop("config$rita=TRUE but frr not found in fit matrix column names.")
    }
  }

  idx
}

split_ai_misclassified_observed <- function(fit_mat, i, idx, temp_m, config, nquar) {
  nstage <- config$ndiagst - 1
  out <- list(
    total = rep(0, nquar),
    cd4 = matrix(0, nrow = nstage, ncol = nquar)
  )

  if (!isTRUE(config$migration)) {
    return(out)
  }

  out$total <- as.numeric(fit_mat[i, idx$mis_dx])

  m_lat_raw <- temp_m[, 1:4, drop = FALSE]
  m_lat <- m_lat_raw
  lat_total <- rowSums(m_lat_raw)
  nz <- lat_total > 0
  if (any(nz)) {
    m_lat[nz, ] <- sweep(m_lat_raw[nz, , drop = FALSE], 1, lat_total[nz], "/")
  }
  if (any(!nz)) {
    m_lat[!nz, ] <- matrix(0.25, nrow = sum(!nz), ncol = 4)
  }

  if (isTRUE(config$rita)) {
    frr <- as.numeric(fit_mat[i, idx$frr])
    out$cd4[1, ] <- out$total * frr
    out$cd4[2:5, ] <- t(sweep(m_lat, 1, out$total * (1 - frr), "*"))
  } else {
    out$cd4 <- t(sweep(m_lat, 1, out$total, "*"))
  }

  out
}

#' Build posterior predictive quantities list
#'
#' Constructs a list of posterior predictive quantities from the provided distributions.
#'
#' @param stan_data Data list used in stan model fitting.
#' @param prev_dist List of prevalence distributions.
#' @param inf_dist List of infection distributions.
#' @param fit_slice Matrix of posterior samples from stan model.
#' @param start_yr Starting year for posterior predictive quantities.
#' @return List of posterior predictive quantities.
build_ppq_list <- function(
  stan_data,
  prev_dist,
  inf_dist,
  fit_slice,
  start_yr
) {
  out <- list(
    incidence = incidence.df(
      stan_data,
      prev_dist,
      start.yr = start_yr,
      annual = FALSE
    ),
    incidence_annual = incidence.df(
      stan_data,
      prev_dist,
      start.yr = start_yr,
      annual = TRUE
    ),
    undiag_prev = undiag.prev.df(
      stan_data,
      prev_dist,
      inf_dist,
      start.yr = start_yr
    ),
    diag_prob = diag.prob.df(
      stan_data,
      fit_slice,
      start.yr = start_yr
    ),
    diagnoses = diagnoses.df(
      stan_data,
      prev_dist,
      start.yr = start_yr
    )
  )
  if (stan_data$model_config$migration) {
    out <- c(
      out,
      list(
        diag_prob_mig = diag.prob.df(
          stan_data,
          fit_slice,
          start.yr = start_yr,
          group = "mig"
        ),
        all_migration = total.migration.df(
          stan_data,
          prev_dist,
          start.yr = start_yr,
          annual = FALSE
        ),
        undiag_migration = undx.migration.df(
          stan_data,
          prev_dist,
          start.yr = start_yr,
          annual = FALSE
        ),
        migration_prob = migration.prob.df(
          stan_data,
          fit_slice,
          start.yr = start_yr
        ),
        transfers = transfers.df(
          stan_data,
          prev_dist,
          start.yr = start_yr
        ),
        transfers_detected = transfers_detected.df(
          stan_data,
          prev_dist,
          fit_slice,
          start.yr = start_yr
        ),
        ratio_abroad_uk = ratio_abroad_uk.df(
          stan_data,
          fit_slice,
          start.yr = start_yr
        ),
        detect_prob = detect.prob.df(
          stan_data,
          fit_slice,
          start.yr = start_yr
        )
      )
    )
  }

  out
}

# Internal helpers for age-independent posterior summaries.
.posterior_summary_vec <- function(x) {
  x <- as.numeric(x)
  qs <- stats::quantile(
    x,
    probs = c(0.025, 0.5, 0.975),
    na.rm = TRUE,
    names = FALSE,
    type = 7
  )
  out <- c(qs[1], qs[2], mean(x, na.rm = TRUE), qs[3])
  names(out) <- c("2.5%", "50%", "mean", "97.5%")
  out
}

#' Infection distribution by year of diagnosis
#'
#' Calculates the distribution of infections by year of diagnosis.
#'
#' @param stan_data Data list used in stan model fitting.
#' @param fit_mat Matrix of posterior samples from stan model.
#' @param i Index of the posterior sample to process.
#' @param idx Named list of column indices from [make_fit_indices()].
#' @return List containing HIV and AIDS diagnosis infection distributions.
infect_dist_by_year_dx <- function(stan_data, fit_mat, i, idx) {
  if (i %% 10 == 0) {
    cat("Iteration Number = ", i, "\n")
  }
  config <- stan_data$model_config

  output <- list()

  # Extract dimensions and initialize variables
  n_CD4 <- nrow(stan_data$CD4)
  nquar <- stan_data$nquar
  if (config$migration) {
    N_sample_uk <- rowSums(stan_data$CD4[, 1, ])
    N_sample_abroad <- rowSums(stan_data$CD4[, 2, ])
  } else {
    N_sample <- rowSums(stan_data$CD4)
  }

  # Get initial prevalence
  init_prev <- if (!("init_prev" %in% names(stan_data))) {
    rep(0, config$nlatst)
  } else {
    stan_data$init_prev
  }

  init_prev_abroad <- if (!("init_prev_abroad" %in% names(stan_data))) {
    rep(0, config$nlatst)
  } else {
    stan_data$init_prev_abroad
  }

  # Get diagnosis probabilities
  d1 <- as.numeric(fit_mat[i, idx$d1])
  d2 <- as.numeric(fit_mat[i, idx$d2])
  d3 <- as.numeric(fit_mat[i, idx$d3])
  d4 <- as.numeric(fit_mat[i, idx$d4])

  # Get infections and migrations
  temp_h <- exp(as.numeric(fit_mat[i, idx$lambda]))

  if (config$migration) {
    temp_omega <- exp(as.numeric(fit_mat[i, idx$gamma]))
    temp_p <- as.numeric(fit_mat[i, idx$p])
    temp_m <- matrix(as.numeric(fit_mat[i, idx$m]), nrow = nquar, byrow = TRUE)
    temp_om <- t(temp_m[, 1:4] * temp_omega)
    misclassified_obs <- split_ai_misclassified_observed(
      fit_mat,
      i,
      idx,
      temp_m,
      config,
      nquar
    )
    d1_mig <- as.numeric(fit_mat[i, idx$d1_mig])
    d2_mig <- as.numeric(fit_mat[i, idx$d2_mig])
    d3_mig <- as.numeric(fit_mat[i, idx$d3_mig])
    d4_mig <- as.numeric(fit_mat[i, idx$d4_mig])
  }

  # Calculate diagnoses
  if (config$rita) {
    d5 <- as.numeric(fit_mat[i, idx$d5])
    diagnoses <- forward_calculate_rita(
      init_prev,
      temp_h,
      d1 = d1,
      d2 = d2,
      d3 = d3,
      d4 = d4,
      d5 = d5,
      q = stan_data$q
    )
  } else {
    diagnoses <- forward_calculate(
      init_prev,
      temp_h,
      d1 = d1,
      d2 = d2,
      d3 = d3,
      d4 = d4,
      q = stan_data$q
    )
  }

  if (config$migration) {
    if (config$rita) {
      d5 <- as.numeric(fit_mat[i, idx$d5])
      d5_mig <- as.numeric(fit_mat[i, idx$d5_mig])
      diagnoses_b <- forward_calculate_rita(
        init_prev_abroad,
        temp_h * temp_p,
        d1 = d1_mig,
        d2 = d2_mig,
        d3 = d3_mig,
        d4 = d4_mig,
        d5 = d5_mig,
        q = stan_data$q
      )
      diagnoses_m <- forward_calculate_m_rita(
        rep(0, length(init_prev_abroad)),
        temp_om,
        d1 = d1_mig,
        d2 = d2_mig,
        d3 = d3_mig,
        d4 = d4_mig,
        d5 = d5_mig,
        q = stan_data$q
      )
    } else {
      diagnoses_b <- forward_calculate(
        init_prev_abroad,
        temp_h * temp_p,
        d1 = d1_mig,
        d2 = d2_mig,
        d3 = d3_mig,
        d4 = d4_mig,
        q = stan_data$q
      )
      diagnoses_m <- forward_calculate_m(
        rep(0, length(init_prev_abroad)),
        temp_om,
        d1 = d1_mig,
        d2 = d2_mig,
        d3 = d3_mig,
        d4 = d4_mig,
        q = stan_data$q
      )
    }

    # Combine diagnoses
    diagnoses_combined <- list(
      HIV.diagnoses = diagnoses$HIV.diagnoses +
        diagnoses_m$HIV.diagnoses +
        diagnoses_b$HIV.diagnoses,
      AIDS.diagnoses = diagnoses$AIDS.diagnoses +
        diagnoses_m$AIDS.diagnoses +
        diagnoses_b$AIDS.diagnoses,
      prevalence = diagnoses$prevalence +
        diagnoses_m$prevalence +
        diagnoses_b$prevalence,
      expected.p = diagnoses$expected.p +
        diagnoses_m$expected.p +
        diagnoses_b$expected.p
    )

    diagnoses_combined_uk <- list(
      HIV.diagnoses = diagnoses$HIV.diagnoses + diagnoses_b$HIV.diagnoses,
      AIDS.diagnoses = diagnoses$AIDS.diagnoses + diagnoses_b$AIDS.diagnoses,
      prevalence = diagnoses$prevalence + diagnoses_b$prevalence,
      expected.p = diagnoses$expected.p + diagnoses_b$expected.p
    )

    diagnoses_combined_m <- list(
      HIV.diagnoses = diagnoses_m$HIV.diagnoses,
      AIDS.diagnoses = diagnoses_m$AIDS.diagnoses,
      prevalence = diagnoses_m$prevalence,
      expected.p = diagnoses_m$expected.p
    )
  }

  # Calculate distributions
  if (config$migration) {
    temp_HIV_dx_uk <- temp_diags_fun(diagnoses_combined_uk)
    d_to_i_HIV_uk <- diag_to_infec_wrap(temp_HIV_dx_uk)
    d_to_i_AIDS_uk <- normalize_mat(t(diagnoses_combined_uk$AIDS.diagnoses))

    temp_HIV_dx_abroad <- temp_diags_fun(diagnoses_combined_m)
    d_to_i_HIV_abroad <- diag_to_infec_wrap(temp_HIV_dx_abroad)
    d_to_i_AIDS_abroad <- normalize_mat(t(diagnoses_combined_m$AIDS.diagnoses))
  } else {
    temp_HIV_dx <- temp_diags_fun(diagnoses)
    d_to_i_HIV <- diag_to_infec_wrap(temp_HIV_dx)
    d_to_i_AIDS <- normalize_mat(t(diagnoses$AIDS.diagnoses))
  }

  # Calculate unreported AIDS
  if (!is.na(idx$under_rep)) {
    under.report <- 1 - as.numeric(fit_mat[i, idx$under_rep])
    under.report <- c(rep(0, 88), rep(under.report, nquar - 88))
  } else {
    under.report <- rep(0, nquar)
  }

  if (config$migration) {
    unreported.AIDS <- stats::rpois(
      nquar,
      under.report * diagnoses_combined_uk$AIDS.diagnoses[1, ]
    )
    for (j in 2:nquar) {
      unreported.AIDS <- rbind(
        unreported.AIDS,
        stats::rpois(
          nquar,
          under.report * diagnoses_combined_uk$AIDS.diagnoses[j, ]
        )
      )
    }
  } else {
    unreported.AIDS <- stats::rpois(
      nquar,
      under.report * diagnoses$AIDS.diagnoses[1, ]
    )
    for (j in 2:nquar) {
      unreported.AIDS <- rbind(
        unreported.AIDS,
        stats::rpois(nquar, under.report * diagnoses$AIDS.diagnoses[j, ])
      )
    }
  }

  # Number of HIV diagnoses with no CD4
  if (config$migration) {
    Dx.no.CD4.uk <- stan_data$HIV[, 1] - c(rep(0, nquar - n_CD4), N_sample_uk)
    Dx.no.CD4.abroad <- stan_data$HIV[, 2] -
      c(rep(0, nquar - n_CD4), N_sample_abroad)
    Dx.no.CD4.b <- Dx.no.CD4.abroad *
      colSums(diagnoses_b$HIV.diagnoses) /
      (colSums(diagnoses_b$HIV.diagnoses) + colSums(diagnoses_m$HIV.diagnoses))
    Dx.no.CD4.m <- Dx.no.CD4.abroad *
      colSums(diagnoses_m$HIV.diagnoses) /
      (colSums(diagnoses_b$HIV.diagnoses) + colSums(diagnoses_m$HIV.diagnoses))

    temp_uk <- sapply(2:nquar, function(x) {
      if (Dx.no.CD4.uk[x] > 0 && sum(d_to_i_HIV_uk[,, x]) > 0) {
        stats::rmultinom(
          1,
          Dx.no.CD4.uk[x] + Dx.no.CD4.b[x],
          d_to_i_HIV_uk[,, x]
        )
      } else {
        matrix(0, nrow = nquar * (config$ndiagst - 1), ncol = 1)
      }
    })

    temp_abroad <- sapply(2:nquar, function(x) {
      if (Dx.no.CD4.m[x] > 0 && sum(d_to_i_HIV_abroad[,, x]) > 0) {
        stats::rmultinom(1, Dx.no.CD4.m[x], d_to_i_HIV_abroad[,, x])
      } else {
        matrix(0, nrow = nquar * (config$ndiagst - 1), ncol = 1)
      }
    })

    temp <- array(
      c(matrix(0, config$ndiagst - 1, nquar), temp_uk + temp_abroad),
      dim = c(config$ndiagst - 1, nquar, nquar)
    )
  } else {
    Dx.no.CD4 <- stan_data$HIV - c(rep(0, nquar - n_CD4), N_sample)

    temp <- sapply(2:nquar, function(x) {
      stats::rmultinom(1, Dx.no.CD4[x], d_to_i_HIV[,, x])
    })
    temp <- array(
      c(matrix(0, config$ndiagst - 1, nquar), temp),
      dim = c(config$ndiagst - 1, nquar, nquar)
    )
  }
  output$HIV.Dx.Infect.Dist <- temp

  # Simulating AIDS diagnoses
  tempvec <- if (config$migration) {
    sapply(5:nquar, function(x) {
      uk_aids <- if (stan_data$AIDS[x, 1] > 0 && sum(d_to_i_AIDS_uk[x, ]) > 0) {
        stats::rmultinom(1, stan_data$AIDS[x, 1], d_to_i_AIDS_uk[x, ])
      } else {
        matrix(0, nrow = nquar, ncol = 1)
      }
      abroad_aids <- if (stan_data$AIDS[x, 2] > 0 && sum(d_to_i_AIDS_abroad[x, ]) > 0) {
        stats::rmultinom(1, stan_data$AIDS[x, 2], d_to_i_AIDS_abroad[x, ])
      } else {
        matrix(0, nrow = nquar, ncol = 1)
      }
      uk_aids + abroad_aids
    })
  } else {
    sapply(5:nquar, function(x) {
      stats::rmultinom(1, stan_data$AIDS[x], d_to_i_AIDS[x, ])
    })
  }

  output$AIDS.Dx.Infect.Dist <- cbind(0, 0, 0, 0, tempvec) + unreported.AIDS

  if (config$rita) {
    prev_for_dist <- if (config$migration) diagnoses_combined$prevalence else diagnoses$prevalence
    collapsed_s1 <- colSums(prev_for_dist[, 1:3, ])
    dist_s1 <- t(collapsed_s1) / colSums(collapsed_s1)

    for (state in 1) {
      d_to_i_temp_1 <- diag_to_infec_wrap1(prev_for_dist, 1) *
        dist_s1[, 1]
      d_to_i_temp_2 <- diag_to_infec_wrap1(prev_for_dist, 2) *
        dist_s1[, 2]
      d_to_i_temp_3 <- diag_to_infec_wrap1(prev_for_dist, 3) *
        dist_s1[, 3]
      d_to_i_temp <- d_to_i_temp_1 + d_to_i_temp_2 + d_to_i_temp_3

      HIV.diags <- if (config$migration) {
        c(rep(0, nquar - n_CD4), rowSums(stan_data$CD4[,, 1]))
      } else {
        c(rep(0, nquar - n_CD4), stan_data$CD4[, 1])
      }
      if (config$migration) {
        HIV.diags <- pmax(0, HIV.diags - misclassified_obs$cd4[1, ])
      }

      tempvec1 <- sapply((state + 1):nquar, function(x) {
        stats::rmultinom(1, HIV.diags[x], d_to_i_temp[x - state, ])
      })
      output$HIV.Dx.Infect.Dist[1, , ] <- output$HIV.Dx.Infect.Dist[1, , ] +
        cbind(matrix(0, nquar, nquar - ncol(tempvec1)), tempvec1)
    }

    for (state in 4:config$nlatst) {
      d_to_i_temp <- diag_to_infec_wrap1(prev_for_dist, state)

      HIV.diags <- if (config$migration) {
        c(rep(0, nquar - n_CD4), rowSums(stan_data$CD4[,, state - 2]))
      } else {
        c(rep(0, nquar - n_CD4), stan_data$CD4[, state - 2])
      }
      if (config$migration) {
        HIV.diags <- pmax(0, HIV.diags - misclassified_obs$cd4[state - 2, ])
      }

      tempvec1 <- sapply((state + 1):nquar, function(x) {
        stats::rmultinom(1, HIV.diags[x], d_to_i_temp[x - state, ])
      })
      output$HIV.Dx.Infect.Dist[state - 2, , ] <- output$HIV.Dx.Infect.Dist[
        state - 2,
        ,
      ] +
        cbind(matrix(0, nquar, nquar - ncol(tempvec1)), tempvec1)
    }
  } else {
    for (state in 1:(config$ndiagst - 1)) {
      d_to_i_temp <- if (config$migration) {
        diag_to_infec_wrap1(diagnoses_combined$prevalence, state)
      } else {
        diag_to_infec_wrap1(diagnoses$prevalence, state)
      }

      HIV.diags <- if (config$migration) {
        c(rep(0, nquar - n_CD4), rowSums(stan_data$CD4[,, state]))
      } else {
        c(rep(0, nquar - n_CD4), stan_data$CD4[, state])
      }
      if (config$migration) {
        HIV.diags <- pmax(0, HIV.diags - misclassified_obs$cd4[state, ])
      }

      tempvec1 <- sapply((state + 1):nquar, function(x) {
        stats::rmultinom(1, HIV.diags[x], d_to_i_temp[x - state, ])
      })
      output$HIV.Dx.Infect.Dist[state, , ] <- output$HIV.Dx.Infect.Dist[
        state,
        ,
      ] +
        cbind(matrix(0, nquar, nquar - ncol(tempvec1)), tempvec1)
    }
  }

  return(output)
}

#' Prevalence distribution by year of infection
#'
#' Calculates the prevalence distribution by year of infection,
#' as estimated in the final quarter of the model.
#'
#' @param stan_data Data list used in stan model fitting.
#' @param fit_mat Matrix of posterior samples from stan model.
#' @param i Index of the posterior sample to process.
#' @param idx Named list of column indices from [make_fit_indices()].
#' @return List containing posterior predictive quantities.
prev_by_year_infection <- function(stan_data, fit_mat, i, idx) {
  if (i %% 10 == 0) {
    cat("Iteration Number = ", i, "\n")
  }
  config <- stan_data$model_config

  output <- list()

  n_CD4 <- nrow(stan_data$CD4)
  nquar <- stan_data$nquar
  if (config$migration) {
    N_sample_uk <- rowSums(stan_data$CD4[, 1, ])
    N_sample_abroad <- rowSums(stan_data$CD4[, 2, ])
  } else {
    N_sample <- rowSums(stan_data$CD4)
  }

  # Get initial prevalence
  init_prev <- if (!("init_prev" %in% names(stan_data))) {
    rep(0, config$nlatst) ### Allows modification for initial prevalence
  } else {
    stan_data$init_prev
  }

  init_prev_abroad <- if (!("init_prev_abroad" %in% names(stan_data))) {
    rep(0, config$nlatst)
  } else {
    stan_data$init_prev_abroad
  }

  # Get diagnosis probabilities
  d1 <- as.numeric(fit_mat[i, idx$d1])
  d2 <- as.numeric(fit_mat[i, idx$d2])
  d3 <- as.numeric(fit_mat[i, idx$d3])
  d4 <- as.numeric(fit_mat[i, idx$d4])

  # Get infections and migrations
  temp_h <- exp(as.numeric(fit_mat[i, idx$lambda]))

  if (config$migration) {
    temp_omega <- exp(as.numeric(fit_mat[i, idx$gamma]))
    temp_p <- as.numeric(fit_mat[i, idx$p])
    temp_m <- matrix(as.numeric(fit_mat[i, idx$m]), nrow = nquar, byrow = TRUE)
    temp_om <- t(temp_m[, 1:4] * temp_omega)
    misclassified_obs <- split_ai_misclassified_observed(
      fit_mat,
      i,
      idx,
      temp_m,
      config,
      nquar
    )
    d1_mig <- as.numeric(fit_mat[i, idx$d1_mig])
    d2_mig <- as.numeric(fit_mat[i, idx$d2_mig])
    d3_mig <- as.numeric(fit_mat[i, idx$d3_mig])
    d4_mig <- as.numeric(fit_mat[i, idx$d4_mig])
  }

  # Calculate diagnoses
  if (config$rita) {
    d5 <- as.numeric(fit_mat[i, idx$d5])
    diagnoses <- forward_calculate_rita(
      init_prev,
      temp_h,
      d1 = d1,
      d2 = d2,
      d3 = d3,
      d4 = d4,
      d5 = d5,
      q = stan_data$q
    )
  } else {
    diagnoses <- forward_calculate(
      init_prev,
      temp_h,
      d1 = d1,
      d2 = d2,
      d3 = d3,
      d4 = d4,
      q = stan_data$q
    )
  }

  if (config$migration) {
    if (config$rita) {
      d5 <- as.numeric(fit_mat[i, idx$d5])
      d5_mig <- as.numeric(fit_mat[i, idx$d5_mig])
      diagnoses_b <- forward_calculate_rita(
        init_prev_abroad,
        temp_h * temp_p,
        d1 = d1_mig,
        d2 = d2_mig,
        d3 = d3_mig,
        d4 = d4_mig,
        d5 = d5_mig,
        q = stan_data$q
      )
      diagnoses_m <- forward_calculate_m_rita(
        rep(0, 7),
        temp_om,
        d1 = d1_mig,
        d2 = d2_mig,
        d3 = d3_mig,
        d4 = d4_mig,
        d5 = d5_mig,
        q = stan_data$q
      )
    } else {
      diagnoses_b <- forward_calculate(
        init_prev_abroad,
        temp_h * temp_p,
        d1 = d1_mig,
        d2 = d2_mig,
        d3 = d3_mig,
        d4 = d4_mig,
        q = stan_data$q
      )
      diagnoses_m <- forward_calculate_m(
        rep(0, 4),
        temp_om,
        d1 = d1_mig,
        d2 = d2_mig,
        d3 = d3_mig,
        d4 = d4_mig,
        q = stan_data$q
      )
    }

    # Combine diagnoses
    diagnoses_combined <- list(
      HIV.diagnoses = diagnoses$HIV.diagnoses +
        diagnoses_m$HIV.diagnoses +
        diagnoses_b$HIV.diagnoses,
      AIDS.diagnoses = diagnoses$AIDS.diagnoses +
        diagnoses_m$AIDS.diagnoses +
        diagnoses_b$AIDS.diagnoses,
      prevalence = diagnoses$prevalence +
        diagnoses_m$prevalence +
        diagnoses_b$prevalence,
      expected.p = diagnoses$expected.p +
        diagnoses_m$expected.p +
        diagnoses_b$expected.p
    )

    diagnoses_combined_uk <- list(
      HIV.diagnoses = diagnoses$HIV.diagnoses + diagnoses_b$HIV.diagnoses,
      AIDS.diagnoses = diagnoses$AIDS.diagnoses + diagnoses_b$AIDS.diagnoses,
      prevalence = diagnoses$prevalence + diagnoses_b$prevalence,
      expected.p = diagnoses$expected.p + diagnoses_b$expected.p
    )

    diagnoses_combined_m <- list(
      HIV.diagnoses = diagnoses_m$HIV.diagnoses,
      AIDS.diagnoses = diagnoses_m$AIDS.diagnoses,
      prevalence = diagnoses_m$prevalence,
      expected.p = diagnoses_m$expected.p
    )
  }

  if (!config$migration) {
    diagnoses_combined <- diagnoses
  }

  output$AIDS.diagnoses <- matrix(
    stats::rpois(
      length(diagnoses_combined$AIDS.diagnoses),
      diagnoses_combined$AIDS.diagnoses
    ),
    nquar
  )

  output$HIV.diagnoses <- matrix(
    stats::rpois(
      length(diagnoses_combined$HIV.diagnoses),
      diagnoses_combined$HIV.diagnoses
    ),
    nquar
  )

  if (config$migration) {
    # Keep transfers in observed HIV diagnoses, but do not back-project them
    # through latent infection histories.
    output$misclassified.HIV <- stats::rpois(nquar, misclassified_obs$total)
  }

  if (config$migration) {
    output$total.migration <- stats::rpois(length(temp_omega), temp_omega)
    output$undx.migration <- stats::rpois(length(temp_omega), colSums(temp_om))
    output$transfers <- stats::rpois(
      length(temp_omega),
      temp_omega * temp_m[, 5]
    )
  }
  # Calculate distributions
  if (config$migration) {
    temp_HIV_dx_uk <- temp_diags_fun(diagnoses_combined_uk)
    d_to_i_HIV_uk <- diag_to_infec_wrap(temp_HIV_dx_uk)
    d_to_i_AIDS_uk <- normalize_mat(t(diagnoses_combined_uk$AIDS.diagnoses))

    temp_HIV_dx_abroad <- temp_diags_fun(diagnoses_combined_m)
    d_to_i_HIV_abroad <- diag_to_infec_wrap(temp_HIV_dx_abroad)
    d_to_i_AIDS_abroad <- normalize_mat(t(diagnoses_combined_m$AIDS.diagnoses))
  } else {
    temp_HIV_dx <- temp_diags_fun(diagnoses)
    d_to_i_HIV <- diag_to_infec_wrap(temp_HIV_dx)
    d_to_i_AIDS <- normalize_mat(t(diagnoses$AIDS.diagnoses))
  }

  # Calculate unreported AIDS
  if (!is.na(idx$under_rep)) {
    under.report <- 1 - as.numeric(fit_mat[i, idx$under_rep])
    under.report <- c(rep(0, 88), rep(under.report, nquar - 88))
  } else {
    under.report <- rep(0, nquar)
  }
  unreported.AIDS <- stats::rpois(
    nquar,
    under.report * diagnoses_combined$AIDS.diagnoses[1, ]
  )
  for (j in 2:nquar) {
    unreported.AIDS <- rbind(
      unreported.AIDS,
      stats::rpois(nquar, under.report * diagnoses_combined$AIDS.diagnoses[j, ])
    )
  }

  output$unreported.AIDS <- rowSums(unreported.AIDS)

  # In migration models the observed HIV/CD4 diagnoses include misclassified
  # arrivals who were already diagnosed abroad. Reconstructing undiagnosed
  # prevalence from diagnosis timings back-projects those cases into pre-arrival
  # history and inflates the early quarters. Preserve the direct latent
  # prevalence trajectory so undiagnosed prevalence can be read off without
  # that back-imputation.
  output$undiagnosed_total <- vapply(seq_len(nquar), function(tt) {
    sum(diagnoses_combined$prevalence[tt, 1:config$nlatst, ])
  }, numeric(1))

  # it appears that the undiagnosed component will only hold data for the most recent year
  output$undiagnosed <- matrix(
    stats::rpois(
      n = prod(dim(diagnoses_combined$prevalence[
        nquar,
        1:config$nlatst,
      ])),
      lambda = diagnoses_combined$prevalence[nquar, 1:config$nlatst, ]
    ),
    nrow = config$nlatst
  )

  # Number of HIV diagnoses with no CD4
  if (config$migration) {
    Dx.no.CD4.uk <- stan_data$HIV[, 1] - c(rep(0, nquar - n_CD4), N_sample_uk)
    Dx.no.CD4.abroad <- stan_data$HIV[, 2] -
      c(rep(0, nquar - n_CD4), N_sample_abroad)
    Dx.no.CD4.b <- Dx.no.CD4.abroad *
      colSums(diagnoses_b$HIV.diagnoses) /
      (colSums(diagnoses_b$HIV.diagnoses) + colSums(diagnoses_m$HIV.diagnoses))
    Dx.no.CD4.m <- Dx.no.CD4.abroad *
      colSums(diagnoses_m$HIV.diagnoses) /
      (colSums(diagnoses_b$HIV.diagnoses) + colSums(diagnoses_m$HIV.diagnoses))

    temp_uk <- sapply(2:nquar, function(x) {
      if (Dx.no.CD4.uk[x] > 0 && sum(d_to_i_HIV_uk[,, x]) > 0) {
        stats::rmultinom(
          1,
          Dx.no.CD4.uk[x] + Dx.no.CD4.b[x],
          d_to_i_HIV_uk[,, x]
        )
      } else {
        matrix(0, nrow = nquar * (config$ndiagst - 1), ncol = 1)
      }
    })

    temp_abroad <- sapply(2:nquar, function(x) {
      if (Dx.no.CD4.m[x] > 0 && sum(d_to_i_HIV_abroad[,, x]) > 0) {
        stats::rmultinom(1, Dx.no.CD4.m[x], d_to_i_HIV_abroad[,, x])
      } else {
        matrix(0, nrow = nquar * (config$ndiagst - 1), ncol = 1)
      }
    })

    temp <- array(
      c(matrix(0, config$ndiagst - 1, nquar), temp_uk + temp_abroad),
      dim = c(config$ndiagst - 1, nquar, nquar)
    )
  } else {
    Dx.no.CD4 <- stan_data$HIV - c(rep(0, nquar - n_CD4), N_sample)

    temp <- sapply(2:nquar, function(x) {
      stats::rmultinom(1, Dx.no.CD4[x], d_to_i_HIV[,, x])
    })
    temp <- array(
      c(matrix(0, config$ndiagst - 1, nquar), temp),
      dim = c(config$ndiagst - 1, nquar, nquar)
    )
  }

  output$HIV.Dx.no.CD4 <- apply(temp, 1:2, sum)
  output$HIV.Dx.CD4 <- matrix(0, config$ndiagst - 1, nquar) ## Initialization needed later

  # Simulating AIDS diagnoses
  tempvec <- if (config$migration) {
    sapply(5:nquar, function(x) {
      uk_aids <- if (stan_data$AIDS[x, 1] > 0 && sum(d_to_i_AIDS_uk[x, ]) > 0) {
        stats::rmultinom(1, stan_data$AIDS[x, 1], d_to_i_AIDS_uk[x, ])
      } else {
        matrix(0, nrow = nquar, ncol = 1)
      }
      abroad_aids <- if (stan_data$AIDS[x, 2] > 0 && sum(d_to_i_AIDS_abroad[x, ]) > 0) {
        stats::rmultinom(1, stan_data$AIDS[x, 2], d_to_i_AIDS_abroad[x, ])
      } else {
        matrix(0, nrow = nquar, ncol = 1)
      }
      uk_aids + abroad_aids
    })
  } else {
    sapply(5:nquar, function(x) {
      stats::rmultinom(1, stan_data$AIDS[x], d_to_i_AIDS[x, ])
    })
  }

  output$reported.AIDS <- rowSums(tempvec)

  if (config$rita) {
    collapsed_s1 <- colSums(diagnoses_combined$prevalence[, 1:3, ])
    dist_s1 <- t(collapsed_s1) / colSums(collapsed_s1)

    for (state in 1) {
      d_to_i_temp_1 <- diag_to_infec_wrap1(diagnoses_combined$prevalence, 1) *
        dist_s1[, 1]
      d_to_i_temp_2 <- diag_to_infec_wrap1(diagnoses_combined$prevalence, 2) *
        dist_s1[, 2]
      d_to_i_temp_3 <- diag_to_infec_wrap1(diagnoses_combined$prevalence, 3) *
        dist_s1[, 3]
      d_to_i_temp <- d_to_i_temp_1 + d_to_i_temp_2 + d_to_i_temp_3

      HIV.diags <- if (config$migration) {
        c(rep(0, nquar - n_CD4), rowSums(stan_data$CD4[,, 1]))
      } else {
        c(rep(0, nquar - n_CD4), stan_data$CD4[, 1])
      }
      if (config$migration) {
        HIV.diags <- pmax(0, HIV.diags - misclassified_obs$cd4[1, ])
      }

      tempvec1 <- sapply((state + 1):nquar, function(x) {
        stats::rmultinom(1, HIV.diags[x], d_to_i_temp[x - state, ])
      })
      output$HIV.Dx.CD4[1, ] <- output$HIV.Dx.CD4[1, ] + rowSums(tempvec1)
    }

    for (state in 4:config$nlatst) {
      d_to_i_temp <- diag_to_infec_wrap1(diagnoses_combined$prevalence, state)

      HIV.diags <- if (config$migration) {
        c(rep(0, nquar - n_CD4), rowSums(stan_data$CD4[,, state - 2]))
      } else {
        c(rep(0, nquar - n_CD4), stan_data$CD4[, state - 2])
      }
      if (config$migration) {
        HIV.diags <- pmax(0, HIV.diags - misclassified_obs$cd4[state - 2, ])
      }

      tempvec1 <- sapply((state + 1):nquar, function(x) {
        stats::rmultinom(1, HIV.diags[x], d_to_i_temp[x - state, ])
      })
      output$HIV.Dx.CD4[state - 2, ] <- rowSums(tempvec1)
    }
  } else {
    for (state in 1:(config$ndiagst - 1)) {
      d_to_i_temp <- diag_to_infec_wrap1(diagnoses_combined$prevalence, state)
      HIV.diags <- if (config$migration) {
        c(
          rep(0, nquar - n_CD4),
          rowSums(stan_data$CD4[,, state])
        )
      } else {
        c(
          rep(0, nquar - n_CD4),
          stan_data$CD4[, state]
        )
      }
      if (config$migration) {
        HIV.diags <- pmax(0, HIV.diags - misclassified_obs$cd4[state, ])
      }

      tempvec1 <- sapply((state + 1):nquar, function(x) {
        stats::rmultinom(1, HIV.diags[x], d_to_i_temp[x - state, ])
      })
      output$HIV.Dx.CD4[state, ] <- rowSums(tempvec1)
    }
  }
  return(output)
}

#' Incidence dataframe
#'
#' Calculates incidence dataframe from prevalence distribution.
#'
#' @param stan_data Data list used in stan model fitting.
#' @param prev_dist Prevalence distribution by year of infection.
#' @param start.yr Starting year for the dataframe.
#' @param annual Logical indicating if annual incidence is required.
#' @return Dataframe of incidence estimates.
incidence.df <- function(
  stan_data,
  prev_dist,
  start.yr = 1978,
  annual = FALSE
) {
  nyr <- stan_data$nquar / 4

  ### Total infections
  ## - Those undx at end of surveillance period (in any states)
  ## - Those Dx with HIV and no CD4 (prev_dist["HIV.Dx.no.CD4",]) - actually the CD4 is "imputed"
  ## - Those Dx with HIV and CD4 (prev_dist["HIV.Dx.CD4",])
  ## - Those Dx with AIDS and reported (prev_dist["unreported.AIDS",])
  ## - Those Dx with AIDS and not reported (prev_dist["reported.AIDS",])

  infections.total <- sapply(prev_dist["undiagnosed", ], function(x) {
    colSums(x)
  }) +
    sapply(prev_dist["HIV.Dx.no.CD4", ], function(x) apply(x, 2, sum)) +
    sapply(prev_dist["HIV.Dx.CD4", ], function(x) apply(x, 2, sum)) +
    simplify2array(prev_dist["unreported.AIDS", ]) +
    simplify2array(prev_dist["reported.AIDS", ])

  ## - Minus those infections acquired abroad and not diagnosed on arrival
  infections.total <- if ("undx.migration" %in% rownames(prev_dist)) {
    infections.total - simplify2array(prev_dist["undx.migration", ])
  } else {
    infections.total
  }

  # Adjustments to remove initial prevalence from first quarter of data
  # usually results in zero infections
  infections.total <- as.matrix(infections.total)
  init_prev_total <- sum(stan_data$init_prev)
  if ("init_prev_abroad" %in% names(stan_data)) {
    init_prev_total <- init_prev_total + sum(stan_data$init_prev_abroad)
  }
  infections.total[1, ] <- pmax(0, infections.total[1, ] - init_prev_total)

  if (annual == FALSE) {
    q.infections.total <- apply(
      infections.total,
      1,
      .posterior_summary_vec
    )
    df <- as.data.frame(t(q.infections.total))
    df$year <- seq(start.yr, start.yr + nyr - 0.25, by = 0.25)
    df <- df |> tidyr::pivot_longer(-"year", names_to = "quartile")
  } else {
    infections.years <- apply(infections.total, 2, function(x) {
      tapply(x, rep(1:(stan_data$nquar / 4), each = 4), sum)
    })
    q.infections.years <- apply(
      infections.years,
      1,
      .posterior_summary_vec
    )
    df <- as.data.frame(t(q.infections.years))
    df$year <- seq(start.yr, start.yr + nyr - 0.25, by = 1)
    df <- df |> tidyr::pivot_longer(-"year", names_to = "quartile")
  }
  return(df)
}

#' Migration dataframe
#'
#' Calculates migration dataframe from prevalence distribution.
#'
#' @param stan_data Data list used in stan model fitting.
#' @param prev_dist Prevalence distribution by year of infection.
#' @param start.yr Starting year for the dataframe.
#' @param annual Logical indicating if annual migration is required.
#' @return Dataframe of migration estimates.
total.migration.df <- function(
  stan_data,
  prev_dist,
  start.yr = 1978,
  annual = FALSE
) {
  nyr <- stan_data$nquar / 4

  migration.total <- simplify2array(prev_dist["total.migration", ])

  if (annual == FALSE) {
    q.migration.total <- apply(
      migration.total,
      1,
      .posterior_summary_vec
    )
    df <- as.data.frame(t(q.migration.total))
    df$year <- seq(start.yr, start.yr + nyr - 0.25, by = 0.25)
    df <- df |> tidyr::pivot_longer(-"year", names_to = "quartile")
  } else {
    migration.years <- apply(migration.total, 2, function(x) {
      tapply(x, rep(1:(stan_data$nquar / 4), each = 4), sum)
    })
    q.migration.years <- apply(
      migration.years,
      1,
      .posterior_summary_vec
    )
    df <- as.data.frame(t(q.migration.years))
    df$year <- seq(start.yr, start.yr + nyr - 0.25, by = 1)
    df <- df |> tidyr::pivot_longer(-"year", names_to = "quartile")
  }
  return(df)
}

#' Transfers dataframe
#'
#' Calculates transfers dataframe from prevalence distribution.
#'
#' @param stan_data Data list used in stan model fitting.
#' @param prev_dist Prevalence distribution by year of infection.
#' @param start.yr Starting year for the dataframe.
#' @param annual Logical indicating if annual transfers is required.
#' @return Dataframe of transfer estimates.
transfers.df <- function(
  stan_data,
  prev_dist,
  start.yr = 1978,
  annual = FALSE
) {
  nyr <- stan_data$nquar / 4

  transfers <- simplify2array(prev_dist["transfers", ])

  if (annual == FALSE) {
    q.transfers <- apply(
      transfers,
      1,
      .posterior_summary_vec
    )
    df <- as.data.frame(t(q.transfers))
    df$year <- seq(start.yr, start.yr + nyr - 0.25, by = 0.25)
    df <- df |> tidyr::pivot_longer(-"year", names_to = "quartile")
  } else {
    transfers.years <- apply(transfers, 2, function(x) {
      tapply(x, rep(1:(stan_data$nquar / 4), each = 4), sum)
    })
    q.transfers.years <- apply(
      transfers.years,
      1,
      .posterior_summary_vec
    )
    df <- as.data.frame(t(q.transfers.years))
    df$year <- seq(start.yr, start.yr + nyr - 0.25, by = 1)
    df <- df |> tidyr::pivot_longer(-"year", names_to = "quartile")
  }
  return(df)
}

#' Detected transfers dataframe
#'
#' Calculates detected transfers dataframe from prevalence distribution and
#' detection probabilities.
#'
#' @param stan_data Data list used in stan model fitting.
#' @param prev_dist Prevalence distribution by year of infection.
#' @param fit_mat Matrix of posterior samples from stan model.
#' @param start.yr Starting year for the dataframe.
#' @param annual Logical indicating if annual detected transfers is required.
#' @return Dataframe of detected transfer estimates.
transfers_detected.df <- function(
  stan_data,
  prev_dist,
  fit_mat,
  start.yr = 1978,
  annual = FALSE
) {
  nyr <- stan_data$nquar / 4
  transfers <- as.matrix(simplify2array(prev_dist["transfers", ]))
  if (nrow(transfers) != stan_data$nquar && ncol(transfers) == stan_data$nquar) {
    transfers <- t(transfers)
  }

  detect.ind <- grep("^detect_prob\\[", colnames(fit_mat))
  detect_draws <- if (length(detect.ind) == 0) {
    matrix(
      rep(stan_data$detect_prob, ncol(transfers)),
      nrow = stan_data$nquar,
      ncol = ncol(transfers)
    )
  } else {
    t(as.matrix(fit_mat[, detect.ind, drop = FALSE]))
  }

  detected_transfers <- transfers * detect_draws

  if (annual == FALSE) {
    q.detected_transfers <- apply(
      detected_transfers,
      1,
      .posterior_summary_vec
    )
    df <- as.data.frame(t(q.detected_transfers))
    df$year <- seq(start.yr, start.yr + nyr - 0.25, by = 0.25)
    df <- df |> tidyr::pivot_longer(-"year", names_to = "quartile")
  } else {
    detected_transfers.years <- apply(detected_transfers, 2, function(x) {
      tapply(x, rep(1:(stan_data$nquar / 4), each = 4), sum)
    })
    q.detected_transfers.years <- apply(
      detected_transfers.years,
      1,
      .posterior_summary_vec
    )
    df <- as.data.frame(t(q.detected_transfers.years))
    df$year <- seq(start.yr, start.yr + nyr - 0.25, by = 1)
    df <- df |> tidyr::pivot_longer(-"year", names_to = "quartile")
  }

  df
}

#' Undiagnosed migration dataframe
#'
#' Calculates undiagnosed migration inflow (omega * (1 - m5)) from prevalence distribution.
#'
#' @param stan_data Data list used in stan model fitting.
#' @param prev_dist Prevalence distribution by year of infection.
#' @param start.yr Starting year for the dataframe.
#' @param annual Logical indicating if annual undiagnosed migration is required.
#' @return Dataframe of undiagnosed migration estimates.
undx.migration.df <- function(
  stan_data,
  prev_dist,
  start.yr = 1978,
  annual = FALSE
) {
  nyr <- stan_data$nquar / 4

  undx_migration <- simplify2array(prev_dist["undx.migration", ])

  if (annual == FALSE) {
    q.undx <- apply(
      undx_migration,
      1,
      .posterior_summary_vec
    )
    df <- as.data.frame(t(q.undx))
    df$year <- seq(start.yr, start.yr + nyr - 0.25, by = 0.25)
    df <- df |> tidyr::pivot_longer(-"year", names_to = "quartile")
  } else {
    undx_years <- apply(undx_migration, 2, function(x) {
      tapply(x, rep(1:(stan_data$nquar / 4), each = 4), sum)
    })
    q.undx_years <- apply(
      undx_years,
      1,
      .posterior_summary_vec
    )
    df <- as.data.frame(t(q.undx_years))
    df$year <- seq(start.yr, start.yr + nyr - 0.25, by = 1)
    df <- df |> tidyr::pivot_longer(-"year", names_to = "quartile")
  }

  df
}

#' Undiagnosed prevalence dataframe
#'
#' Calculates undiagnosed prevalence dataframe from prevalence and infection distributions.
#'
#' @param stan_data Data list used in stan model fitting.
#' @param prev_dist Prevalence distribution by year of infection.
#' @param inf_dist Infection distribution by diagnosis.
#' @param start.yr Starting year for the dataframe.
#' @param annual Logical indicating if annual undiagnosed prevalence is required.
#' @return Dataframe of undiagnosed prevalence estimates.
undiag.prev.df <- function(
  stan_data,
  prev_dist,
  inf_dist,
  start.yr = 1978,
  annual = FALSE
) {
  nyr <- stan_data$nquar / 4
  use_direct_prev <- "undiagnosed_total" %in% rownames(prev_dist)

  if (use_direct_prev) {
    undiagnosed.prevalence <- as.matrix(
      simplify2array(prev_dist["undiagnosed_total", ])
    )
    if (nrow(undiagnosed.prevalence) != stan_data$nquar) {
      undiagnosed.prevalence <- t(undiagnosed.prevalence)
    }
  } else {
    ## Undiagnosed prevalence is equal to the number of people
    ## that are undiagnosed by the end of the surveillance period
    undiagnosed.prevalence <- sapply(prev_dist["undiagnosed", ], function(x) {
      cumsum(colSums(x))
    })

    ### Plus all those that were diagnosed (AIDS or HIV)
    ### between the infection time and dx time
    ### these are counted by the prev_fun(1)
    undiagnosed.prevalence.AIDS <- sapply(
      inf_dist["AIDS.Dx.Infect.Dist", ],
      function(x) t(prev_fun(x))
    )
    undiagnosed.prevalence.HIV <- sapply(
      inf_dist["HIV.Dx.Infect.Dist", ],
      function(x) t(prev_fun1(x))
    )

    undiagnosed.prevalence <- undiagnosed.prevalence +
      undiagnosed.prevalence.AIDS +
      undiagnosed.prevalence.HIV
  }

  q.undiagnosed.prevalence <- apply(
    undiagnosed.prevalence,
    1,
    .posterior_summary_vec
  )[, 1:(stan_data$nquar)]
  df <- as.data.frame(t(q.undiagnosed.prevalence))
  df$year <- seq(start.yr, start.yr + nyr - 0.25, by = 0.25)
  df <- df |> tidyr::pivot_longer(-"year", names_to = "quartile")

  if (annual == TRUE) {
    df <- df |>
      dplyr::filter(.data$year %% 1 == 0.75) |>
      dplyr::mutate(year = .data$year - 0.75)
  }

  return(df)
}

#' Diagnosis probabilities dataframe
#'
#' Calculates diagnosis probabilities dataframe from posterior samples.
#'
#' @param stan_data Data list used in stan model fitting.
#' @param fit_mat Matrix of posterior samples from stan model.
#' @param start.yr Starting year for the dataframe.
#' @param group Character indicating which diagnosis probabilities to extract
#'   (`"nat"` for UK-born, `"mig"` for born abroad).
#' @return Dataframe of diagnosis probabilities estimates.
diag.prob.df <- function(stan_data, fit_mat, start.yr = 1978, group = c("nat", "mig")) {
  group <- match.arg(group)
  nyr <- stan_data$nquar / 4

  # Migration models expose d_nat[k,t] / d_mig[k,t]; non-migration models use d[k,t].
  # group = "nat" selects d_nat (or d for non-migration); group = "mig" selects d_mig.
  dp <- if (any(grepl("^d_nat\\[", colnames(fit_mat)))) {
    if (group == "mig") "d_mig" else "d_nat"
  } else {
    "d"
  }

  ### Quantiles of diagnosis
  d1.q <- apply(
    fit_mat[, grep(sprintf("^%s\\[1,", dp), colnames(fit_mat))],
    2,
    .posterior_summary_vec
  )
  d2.q <- apply(
    fit_mat[, grep(sprintf("^%s\\[2,", dp), colnames(fit_mat))],
    2,
    .posterior_summary_vec
  )
  d3.q <- apply(
    fit_mat[, grep(sprintf("^%s\\[3,", dp), colnames(fit_mat))],
    2,
    .posterior_summary_vec
  )
  d4.q <- apply(
    fit_mat[, grep(sprintf("^%s\\[4,", dp), colnames(fit_mat))],
    2,
    .posterior_summary_vec
  )

  year <- seq(start.yr, start.yr + nyr - 0.25, by = 0.25)

  if (any(grepl(sprintf("^%s\\[5,", dp), colnames(fit_mat)))) {
    d5.q <- apply(
        fit_mat[, grep(sprintf("^%s\\[5,", dp), colnames(fit_mat))],
        2,
        .posterior_summary_vec
      )

    d1.q <- t(d1.q) |>
      as.data.frame() |>
      cbind(year, strata = "Recently acquired") |>
      tidyr::pivot_longer(-c("year", "strata"), names_to = "quartile")
    d2.q <- t(d2.q) |>
      as.data.frame() |>
      cbind(year, strata = "CD4>=500") |>
      tidyr::pivot_longer(-c("year", "strata"), names_to = "quartile")
    d3.q <- t(d3.q) |>
      as.data.frame() |>
      cbind(year, strata = "CD4 350-499") |>
      tidyr::pivot_longer(-c("year", "strata"), names_to = "quartile")
    d4.q <- t(d4.q) |>
      as.data.frame() |>
      cbind(year, strata = "CD4 200-349") |>
      tidyr::pivot_longer(-c("year", "strata"), names_to = "quartile")
    d5.q <- t(d5.q) |>
      as.data.frame() |>
      cbind(year, strata = "CD4<200") |>
      tidyr::pivot_longer(-c("year", "strata"), names_to = "quartile")
    df <- rbind(d1.q, d2.q, d3.q, d4.q, d5.q)
  } else {
    d1.q <- t(d1.q) |>
      as.data.frame() |>
      cbind(year, strata = "CD4>=500") |>
      tidyr::pivot_longer(-c("year", "strata"), names_to = "quartile")
    d2.q <- t(d2.q) |>
      as.data.frame() |>
      cbind(year, strata = "CD4 350-499") |>
      tidyr::pivot_longer(-c("year", "strata"), names_to = "quartile")
    d3.q <- t(d3.q) |>
      as.data.frame() |>
      cbind(year, strata = "CD4 200-349") |>
      tidyr::pivot_longer(-c("year", "strata"), names_to = "quartile")
    d4.q <- t(d4.q) |>
      as.data.frame() |>
      cbind(year, strata = "CD4<200") |>
      tidyr::pivot_longer(-c("year", "strata"), names_to = "quartile")
    df <- rbind(d1.q, d2.q, d3.q, d4.q)
  }
  return(df)
}

#' Diagnoses dataframe
#'
#' Calculates diagnoses dataframe from prevalence distribution.
#'
#' @param stan_data Data list used in stan model fitting.
#' @param prev_dist Prevalence distribution by year of infection.
#' @param start.yr Starting year for the dataframe.
#' @param annual Logical indicating if annual diagnoses is required.
#' @return Dataframe of diagnoses estimates.
diagnoses.df <- function(
  stan_data,
  prev_dist,
  start.yr = 1978,
  annual = FALSE
) {
  nyr <- (stan_data$nquar) / 4

  diagnosed.AIDS <- sapply(prev_dist["AIDS.diagnoses", ], function(x) {
    colSums(x)
  })

  diagnosed.HIV <- sapply(prev_dist["HIV.diagnoses", ], function(x) {
    colSums(x)
  })
  if ("misclassified.HIV" %in% rownames(prev_dist)) {
    diagnosed.HIV <- diagnosed.HIV + matrix(
      simplify2array(prev_dist["misclassified.HIV", ]),
      nrow = stan_data$nquar
    )
  }

  q.diagnosed <- apply(
    diagnosed.HIV + diagnosed.AIDS,
    1,
    .posterior_summary_vec
  )[, 1:(stan_data$nquar)]
  df <- as.data.frame(t(q.diagnosed))
  df$year <- seq(start.yr, start.yr + nyr - 0.25, by = 0.25)
  df <- df |> tidyr::pivot_longer(-"year", names_to = "quartile")

  if (annual == TRUE) {
    df <- df |>
      dplyr::filter(.data$year %% 1 == 0.75) |>
      dplyr::mutate(year = .data$year - 0.75)
  }

  return(df)
}

#' Migration probabilities dataframe
#'
#' Calculates migration probabilities dataframe from posterior samples.
#'
#' @param stan_data Data list used in stan model fitting.
#' @param fit_mat Matrix of posterior samples from stan model.
#' @param start.yr Starting year for the dataframe.
#' @return Dataframe of migration probabilities estimates.
migration.prob.df <- function(stan_data, fit_mat, start.yr = 1978) {
  nyr <- stan_data$nquar / 4

  m1.ind <- grep("^m\\[1", colnames(fit_mat))
  m2.ind <- grep("^m\\[2", colnames(fit_mat))
  m3.ind <- grep("^m\\[3", colnames(fit_mat))
  m4.ind <- grep("^m\\[4", colnames(fit_mat))
  m5.ind <- grep("^m\\[5", colnames(fit_mat))

  ### Quantiles of diagnosis
  m1.q <- apply(
    fit_mat[, m1.ind],
    2,
    .posterior_summary_vec
  )
  m2.q <- apply(
    fit_mat[, m2.ind],
    2,
    .posterior_summary_vec
  )
  m3.q <- apply(
    fit_mat[, m3.ind],
    2,
    .posterior_summary_vec
  )
  m4.q <- apply(
    fit_mat[, m4.ind],
    2,
    .posterior_summary_vec
  )
  m5.q <- apply(
    fit_mat[, m5.ind],
    2,
    .posterior_summary_vec
  )
  year <- seq(start.yr, start.yr + nyr - 0.25, by = 0.25)

  m1.q <- t(m1.q) |>
    as.data.frame() |>
    cbind(year, strata = "CD4>=500") |>
    tidyr::pivot_longer(-c("year", "strata"), names_to = "quartile")
  m2.q <- t(m2.q) |>
    as.data.frame() |>
    cbind(year, strata = "CD4 350-499") |>
    tidyr::pivot_longer(-c("year", "strata"), names_to = "quartile")
  m3.q <- t(m3.q) |>
    as.data.frame() |>
    cbind(year, strata = "CD4 200-349") |>
    tidyr::pivot_longer(-c("year", "strata"), names_to = "quartile")
  m4.q <- t(m4.q) |>
    as.data.frame() |>
    cbind(year, strata = "CD4<200") |>
    tidyr::pivot_longer(-c("year", "strata"), names_to = "quartile")
  m5.q <- t(m5.q) |>
    as.data.frame() |>
    cbind(year, strata = "Transfers of care") |>
    tidyr::pivot_longer(-c("year", "strata"), names_to = "quartile")
  df <- rbind(m1.q, m2.q, m3.q, m4.q, m5.q)
  return(df)
}

#' UK-acquired incidence-ratio dataframe
#'
#' Calculates the born-abroad:UK-born ratio of UK-acquired incidence from
#' posterior samples.
#'
#' @param stan_data Data list used in stan model fitting.
#' @param fit_mat Matrix of posterior samples from stan model.
#' @param start.yr Starting year for the dataframe.
#' @return Data frame of born-abroad:UK-born UK-acquired incidence-ratio
#'   estimates.
ratio_abroad_uk.df <- function(stan_data, fit_mat, start.yr = 1978) {
  nyr <- stan_data$nquar / 4

  puk.ind <- grep("^p\\[", colnames(fit_mat))

  ### Quantiles of diagnosis
  if (length(puk.ind) == 0) {
    puk.q <- matrix(rep(stan_data$p, each = 4), nrow = 4)
    rownames(puk.q) <- c("2.5%", "50%", "mean", "97.5%")
  } else {
    puk.q <- apply(
      fit_mat[, puk.ind],
      2,
      .posterior_summary_vec
    )
  }

  year <- seq(start.yr, start.yr + nyr - 0.25, by = 0.25)

  df <- t(puk.q) |>
    as.data.frame() |>
    cbind(year) |>
    tidyr::pivot_longer(-"year", names_to = "quartile")

  return(df)
}

#' Detection probabilities dataframe
#'
#' Calculates detection probabilities dataframe from posterior samples.
#'
#' @param stan_data Data list used in stan model fitting.
#' @param fit_mat Matrix of posterior samples from stan model.
#' @param start.yr Starting year for the dataframe.
#' @return Dataframe of detection probabilities estimates.
detect.prob.df <- function(stan_data, fit_mat, start.yr = 1978) {
  nyr <- stan_data$nquar / 4

  detect.ind <- grep("^detect_prob\\[", colnames(fit_mat))

  ### Quantiles of diagnosis
  if (length(detect.ind) == 0) {
    detect.q <- matrix(rep(stan_data$detect_prob, each = 4), nrow = 4)
    rownames(detect.q) <- c("2.5%", "50%", "mean", "97.5%")
  } else {
    detect.q <- apply(
      fit_mat[, detect.ind],
      2,
      .posterior_summary_vec
    )
  }

  year <- seq(start.yr, start.yr + nyr - 0.25, by = 0.25)

  df <- t(detect.q) |>
    as.data.frame() |>
    cbind(year) |>
    tidyr::pivot_longer(-"year", names_to = "quartile")

  return(df)
}
