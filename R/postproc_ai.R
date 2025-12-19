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
#' sim_diags <- simulate_diagnoses(sim_seed = 123)
#' stan_data <- model_type_prep(sim_diags, model_seed = 123)
#' fit <- fit_stan_model(
#'    stan_data,
#'    iter = 10,
#'    chains = 1,
#'    refresh = 0,
#'    messages = FALSE
#' )
#' hiv_list <- postproc(fit, stan_data, start_yr = 0, nworkers = 1)
#' @export
postproc <- function(
  fit,
  stan_data,
  start_yr = 0,
  nworkers = 1,
  chain_results = FALSE
) {
  # Run the appropriate post-processing routine
  if (stan_data$config$age) {
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
#' @export
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
  options(future.globals.maxSize = 3000 * 1024^2)

  # generate matrix object
  fit_mat <- as.matrix(fit)

  # make fit indices
  idx <- make_fit_indices(fit_mat, stan_data$config, stan_data$nquar)

  if (chain_results) {
    nchains <- as.integer(fit@sim$chains)
    draw_splits <- split_draws(nrow(fit_mat), nchains)

    hiv_list <- list()
    for (ch in seq_along(draw_splits)) {
      draw_ids <- draw_splits[[ch]]

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

      if (stan_data$config$migration) {
        hiv_list[[paste0("migration_", ch)]] <- ppq$migration
        hiv_list[[paste0("migration_prob_", ch)]] <- ppq$migration_prob
        hiv_list[[paste0("transfers_", ch)]] <- ppq$transfers
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

  idx <- list(
    d1 = grep("^d\\[1,", cn),
    d2 = grep("^d\\[2,", cn),
    d3 = grep("^d\\[3,", cn),
    d4 = grep("^d\\[4,", cn),
    d5 = grep("^d\\[5,", cn),
    lambda = grep("^lambda\\[", cn),
    gamma = grep("^gamma\\[", cn),
    p = grep("^p\\[", cn),
    m = grep("^m\\[", cn),
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
  }

  idx
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
  if (stan_data$config$migration) {
    out <- c(
      out,
      list(
        migration = migration.df(
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

#' Split draws by chain
#'
#' Splits the total number of draws into groups corresponding to each chain.
#'
#' @param n_draws Total number of draws.
#' @param nchains Number of chains.
#' @return List of indices for each chain.
split_draws <- function(n_draws, nchains) {
  nchains <- max(1L, as.integer(nchains))
  if (nchains > 1) {
    split(
      seq_len(n_draws),
      cut(seq_len(n_draws), breaks = nchains, labels = FALSE)
    )
  } else {
    list(seq_len(n_draws))
  }
}

#' Infection distribution by year of diagnosis
#'
#' Calculates the distribution of infections by year of diagnosis.
#'
#' @param stan_data Data list used in stan model fitting.
#' @param fit_mat Matrix of posterior samples from stan model.
#' @param i Index of the posterior sample to process.
#' @param idx Additional index parameter for processing.
#' @return List containing HIV and AIDS diagnosis infection distributions.
infect_dist_by_year_dx <- function(stan_data, fit_mat, i, idx) {
  if (i %% 10 == 0) {
    cat("Iteration Number = ", i, "\n")
  }
  config <- stan_data$config

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
  d1 <- fit_mat[i, idx$d1]
  d2 <- fit_mat[i, idx$d2]
  d3 <- fit_mat[i, idx$d3]
  d4 <- fit_mat[i, idx$d4]

  # Get infections and migrations
  temp_h <- exp(fit_mat[i, idx$lambda])
  if (config$migration) {
    temp_omega <- exp(fit_mat[i, idx$gamma])
    temp_p <- fit_mat[i, idx$p]
    temp_m <- matrix(
      fit_mat[i, idx$m],
      nrow = nquar,
      byrow = TRUE
    )
    temp_om <- t(temp_m[, 1:4] * temp_omega)
  }

  # Calculate diagnoses
  if (config$rita) {
    d5 <- fit_mat[i, idx$d5]
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
      d5 <- fit_mat[i, idx$d5]
      diagnoses_b <- forward_calculate_rita(
        init_prev_abroad / 2,
        temp_h * temp_p,
        d1 = d1,
        d2 = d2,
        d3 = d3,
        d4 = d4,
        d5 = d5,
        q = stan_data$q
      )
      diagnoses_m <- forward_calculate_m_rita(
        init_prev_abroad / 2,
        temp_om,
        d1 = d1,
        d2 = d2,
        d3 = d3,
        d4 = d4,
        d5 = d5,
        q = stan_data$q
      )
    } else {
      diagnoses_b <- forward_calculate(
        init_prev_abroad / 2,
        temp_h * temp_p,
        d1 = d1,
        d2 = d2,
        d3 = d3,
        d4 = d4,
        q = stan_data$q
      )
      diagnoses_m <- forward_calculate_m(
        init_prev_abroad / 2,
        temp_om,
        d1 = d1,
        d2 = d2,
        d3 = d3,
        d4 = d4,
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
    d_to_i_AIDS_uk <- normalize_mat(t(diagnoses_combined_uk$AIDS))

    temp_HIV_dx_abroad <- temp_diags_fun(diagnoses_combined_m)
    d_to_i_HIV_abroad <- diag_to_infec_wrap(temp_HIV_dx_abroad)
    d_to_i_AIDS_abroad <- normalize_mat(t(diagnoses_combined_m$AIDS))
  } else {
    temp_HIV_dx <- temp_diags_fun(diagnoses)
    d_to_i_HIV <- diag_to_infec_wrap(temp_HIV_dx)
    d_to_i_AIDS <- normalize_mat(t(diagnoses$AIDS))
  }

  # Calculate unreported AIDS
  if (!is.na(idx$under_rep)) {
    under.report <- 1 - fit_mat[i, idx$under_rep]
    under.report <- c(rep(0, 88), rep(under.report, nquar - 88))
  } else {
    under.report <- rep(0, nquar)
  }

  if (config$migration) {
    unreported.AIDS <- stats::rpois(
      nquar,
      under.report * diagnoses_combined_uk$AIDS[1, ]
    )
    for (j in 2:nquar) {
      unreported.AIDS <- rbind(
        unreported.AIDS,
        stats::rpois(nquar, under.report * diagnoses_combined_uk$AIDS[j, ])
      )
    }
  } else {
    unreported.AIDS <- stats::rpois(nquar, under.report * diagnoses$AIDS[1, ])
    for (j in 2:nquar) {
      unreported.AIDS <- rbind(
        unreported.AIDS,
        stats::rpois(nquar, under.report * diagnoses$AIDS[j, ])
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
      if (sum(stan_data$AIDS[x, ]) > 0 && sum(d_to_i_AIDS_uk[x, ]) > 0) {
        stats::rmultinom(1, sum(stan_data$AIDS[x, ]), d_to_i_AIDS_uk[x, ])
      } else {
        matrix(0, nrow = nquar, ncol = 1)
      }
    })
  } else {
    sapply(5:nquar, function(x) {
      stats::rmultinom(1, stan_data$AIDS[x], d_to_i_AIDS[x, ])
    })
  }

  output$AIDS.Dx.Infect.Dist <- cbind(0, 0, 0, 0, tempvec) + unreported.AIDS

  if (config$rita) {
    collapsed_s1 <- colSums(diagnoses$prevalence[, 1:3, ])
    dist_s1 <- t(collapsed_s1) / colSums(collapsed_s1)

    for (state in 1) {
      #prev_list <- array(sapply(1:nquar, function(x) colSums(diagnoses$prevalence[x, 1:3, ])), dim = c(nquar, 1, nquar))
      d_to_i_temp_1 <- diag_to_infec_wrap1(diagnoses$prevalence, 1) *
        dist_s1[, 1]
      d_to_i_temp_2 <- diag_to_infec_wrap1(diagnoses$prevalence, 2) *
        dist_s1[, 2]
      d_to_i_temp_3 <- diag_to_infec_wrap1(diagnoses$prevalence, 3) *
        dist_s1[, 3]
      d_to_i_temp <- d_to_i_temp_1 + d_to_i_temp_2 + d_to_i_temp_3

      HIV.diags <- if (config$migration) {
        c(rep(0, nquar - n_CD4), rowSums(stan_data$CD4[,, 1]))
      } else {
        c(rep(0, nquar - n_CD4), stan_data$CD4[, 1])
      }

      tempvec1 <- sapply((state + 1):nquar, function(x) {
        stats::rmultinom(1, HIV.diags[x], d_to_i_temp[x - state, ])
      })
      output$HIV.Dx.Infect.Dist[1, , ] <- output$HIV.Dx.Infect.Dist[1, , ] +
        cbind(matrix(0, nquar, nquar - ncol(tempvec1)), tempvec1)
    }

    for (state in 4:config$nlatst) {
      d_to_i_temp <- diag_to_infec_wrap1(diagnoses$prevalence, state)

      HIV.diags <- if (config$migration) {
        c(rep(0, nquar - n_CD4), rowSums(stan_data$CD4[,, state - 2]))
      } else {
        c(rep(0, nquar - n_CD4), stan_data$CD4[, state - 2])
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
#' @param idx Additional index parameter for processing.
#' @return List containing posterior predictive quantities.
prev_by_year_infection <- function(stan_data, fit_mat, i, idx) {
  if (i %% 10 == 0) {
    cat("Iteration Number = ", i, "\n")
  }
  config <- stan_data$config

  output <- list()

  n_CD4 <- nrow(stan_data$CD4)
  nquar <- stan_data$nquar
  N_sample <- rowSums(stan_data$CD4)

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
  d1 <- fit_mat[i, idx$d1]
  d2 <- fit_mat[i, idx$d2]
  d3 <- fit_mat[i, idx$d3]
  d4 <- fit_mat[i, idx$d4]

  # Get infections and migrations
  temp_h <- exp(fit_mat[i, idx$lambda])

  if (config$migration) {
    temp_omega <- exp(fit_mat[i, idx$gamma])
    temp_p <- fit_mat[i, idx$p]
    temp_m <- matrix(
      fit_mat[i, idx$m],
      nrow = nquar,
      byrow = TRUE
    )
    temp_om <- t(temp_m[, 1:4] * temp_omega)
  }

  # Calculate diagnoses
  if (config$rita) {
    d5 <- fit_mat[i, idx$d5]
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
      d5 <- fit_mat[i, idx$d5]
      diagnoses_b <- forward_calculate_rita(
        init_prev_abroad,
        temp_h * temp_p,
        d1 = d1,
        d2 = d2,
        d3 = d3,
        d4 = d4,
        d5 = d5,
        q = stan_data$q
      )
      diagnoses_m <- forward_calculate_m_rita(
        rep(0, 7),
        temp_om,
        d1 = d1,
        d2 = d2,
        d3 = d3,
        d4 = d4,
        d5 = d5,
        q = stan_data$q
      )
    } else {
      diagnoses_b <- forward_calculate(
        init_prev_abroad,
        temp_h * temp_p,
        d1 = d1,
        d2 = d2,
        d3 = d3,
        d4 = d4,
        q = stan_data$q
      )
      diagnoses_m <- forward_calculate_m(
        rep(0, 4),
        temp_om,
        d1 = d1,
        d2 = d2,
        d3 = d3,
        d4 = d4,
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
    output$total.migration <- stats::rpois(length(temp_omega), temp_omega)
    output$undx.migration <- stats::rpois(length(temp_omega), colSums(temp_om))
    output$transfers <- stats::rpois(
      length(temp_omega),
      temp_omega * temp_m[, 5]
    )
  }
  # output$expected.p <- matrix(rpois(length(diagnoses$expected.p),diagnoses$expected.p),4,nquar,nquar)
  # Calculate distributions
  if (config$migration) {
    temp_HIV_dx_uk <- temp_diags_fun(diagnoses_combined_uk)
    d_to_i_HIV_uk <- diag_to_infec_wrap(temp_HIV_dx_uk)
    d_to_i_AIDS_uk <- normalize_mat(t(diagnoses_combined_uk$AIDS))

    temp_HIV_dx_abroad <- temp_diags_fun(diagnoses_combined_m)
    d_to_i_HIV_abroad <- diag_to_infec_wrap(temp_HIV_dx_abroad)
    d_to_i_AIDS_abroad <- normalize_mat(t(diagnoses_combined_m$AIDS))
  } else {
    temp_HIV_dx <- temp_diags_fun(diagnoses)
    d_to_i_HIV <- diag_to_infec_wrap(temp_HIV_dx)
    d_to_i_AIDS <- normalize_mat(t(diagnoses$AIDS))
  }

  # Calculate unreported AIDS
  if (!is.na(idx$under_rep)) {
    under.report <- 1 - fit_mat[i, idx$under_rep]
    under.report <- c(rep(0, 88), rep(under.report, nquar - 88))
  } else {
    under.report <- rep(0, nquar)
  }
  unreported.AIDS <- stats::rpois(
    nquar,
    under.report * diagnoses_combined$AIDS[1, ]
  )
  for (j in 2:nquar) {
    unreported.AIDS <- rbind(
      unreported.AIDS,
      stats::rpois(nquar, under.report * diagnoses_combined$AIDS[j, ])
    )
  }

  output$unreported.AIDS <- rowSums(unreported.AIDS)

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
  Dx.no.CD4 <- if (config$migration) {
    rowSums(stan_data$HIV) - c(rep(0, nquar - n_CD4), N_sample)
  } else {
    stan_data$HIV - c(rep(0, nquar - n_CD4), N_sample)
  }

  temp <- sapply(2:nquar, function(x) {
    if (Dx.no.CD4[x] > 0 && sum(d_to_i_HIV[,, x]) > 0) {
      stats::rmultinom(1, Dx.no.CD4[x], d_to_i_HIV[,, x])
    } else {
      matrix(0, nrow = nquar * (config$ndiagst - 1), ncol = 1)
    }
  })
  temp <- array(
    c(matrix(0, config$ndiagst - 1, nquar), temp),
    dim = c(config$ndiagst - 1, nquar, nquar)
  )

  output$HIV.Dx.no.CD4 <- apply(temp, 1:2, sum)
  output$HIV.Dx.CD4 <- matrix(0, config$ndiagst - 1, nquar) ## Initialization needed later

  # Simulating AIDS diagnoses
  tempvec <- if (config$migration) {
    sapply(5:nquar, function(x) {
      if (sum(stan_data$AIDS[x, ]) > 0 && sum(d_to_i_AIDS_uk[x, ]) > 0) {
        stats::rmultinom(1, sum(stan_data$AIDS[x, ]), d_to_i_AIDS_uk[x, ])
      } else {
        matrix(0, nrow = nquar, ncol = 1)
      }
    })
  } else {
    sapply(5:nquar, function(x) {
      stats::rmultinom(1, stan_data$AIDS[x], d_to_i_AIDS[x, ])
    })
  }

  output$reported.AIDS <- rowSums(tempvec)

  if (config$rita) {
    collapsed_s1 <- colSums(diagnoses$prevalence[, 1:3, ])
    dist_s1 <- t(collapsed_s1) / colSums(collapsed_s1)

    for (state in 1) {
      #prev_list <- array(sapply(1:nquar, function(x) colSums(diagnoses$prevalence[x, 1:3, ])), dim = c(nquar, 1, nquar))
      d_to_i_temp_1 <- diag_to_infec_wrap1(diagnoses$prevalence, 1) *
        dist_s1[, 1]
      d_to_i_temp_2 <- diag_to_infec_wrap1(diagnoses$prevalence, 2) *
        dist_s1[, 2]
      d_to_i_temp_3 <- diag_to_infec_wrap1(diagnoses$prevalence, 3) *
        dist_s1[, 3]
      d_to_i_temp <- d_to_i_temp_1 + d_to_i_temp_2 + d_to_i_temp_3

      HIV.diags <- if (config$migration) {
        c(rep(0, nquar - n_CD4), rowSums(stan_data$CD4[,, 1]))
      } else {
        c(rep(0, nquar - n_CD4), stan_data$CD4[, 1])
      }

      tempvec1 <- sapply((state + 1):nquar, function(x) {
        stats::rmultinom(1, HIV.diags[x], d_to_i_temp[x - state, ])
      })
      output$HIV.Dx.CD4[1, ] <- output$HIV.Dx.CD4[1, ] + rowSums(tempvec1)
    }

    for (state in 4:config$nlatst) {
      d_to_i_temp <- diag_to_infec_wrap1(diagnoses$prevalence, state)

      HIV.diags <- if (config$migration) {
        c(rep(0, nquar - n_CD4), rowSums(stan_data$CD4[,, state - 2]))
      } else {
        c(rep(0, nquar - n_CD4), stan_data$CD4[, state - 2])
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

  # remove initial prevalence from first quarter
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
      stats::quantile,
      probs = c(0.025, 0.5, 0.975)
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
      stats::quantile,
      probs = c(0.025, 0.5, 0.975)
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
migration.df <- function(
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
      stats::quantile,
      probs = c(0.025, 0.5, 0.975)
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
      stats::quantile,
      probs = c(0.025, 0.5, 0.975)
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
      stats::quantile,
      probs = c(0.025, 0.5, 0.975)
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
      stats::quantile,
      probs = c(0.025, 0.5, 0.975)
    )
    df <- as.data.frame(t(q.transfers.years))
    df$year <- seq(start.yr, start.yr + nyr - 0.25, by = 1)
    df <- df |> tidyr::pivot_longer(-"year", names_to = "quartile")
  }
  return(df)
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
  # start.yr <- 1995
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

  q.undiagnosed.prevalence <- apply(
    undiagnosed.prevalence,
    1,
    stats::quantile,
    probs = c(0.025, 0.5, 0.975)
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
#' @return Dataframe of diagnosis probabilities estimates.
diag.prob.df <- function(stan_data, fit_mat, start.yr = 1978) {
  nyr <- stan_data$nquar / 4

  ### Quantiles of diagnosis
  d1.q <- apply(
    fit_mat[, grep("^d\\[1,", colnames(fit_mat))],
    2,
    stats::quantile,
    probs = c(0.025, 0.5, 0.975)
  )
  d2.q <- apply(
    fit_mat[, grep("^d\\[2,", colnames(fit_mat))],
    2,
    stats::quantile,
    probs = c(0.025, 0.5, 0.975)
  )
  d3.q <- apply(
    fit_mat[, grep("^d\\[3,", colnames(fit_mat))],
    2,
    stats::quantile,
    probs = c(0.025, 0.5, 0.975)
  )
  d4.q <- apply(
    fit_mat[, grep("^d\\[4,", colnames(fit_mat))],
    2,
    stats::quantile,
    probs = c(0.025, 0.5, 0.975)
  )

  year <- seq(start.yr, start.yr + nyr - 0.25, by = 0.25)

  if (any(grepl("^d\\[5,", colnames(fit_mat)))) {
    d5.ind <-
      d5.q <- apply(
        fit_mat[, grep("^d\\[5,", colnames(fit_mat))],
        2,
        stats::quantile,
        probs = c(0.025, 0.5, 0.975)
      )

    d1.q <- t(d1.q) |>
      as.data.frame() |>
      cbind(year, strata = "Recently acquired") |>
      tidyr::pivot_longer(-c("year", "strata"), names_to = "quartile")
    d2.q <- t(d2.q) |>
      as.data.frame() |>
      cbind(year, strata = "CD4>500") |>
      tidyr::pivot_longer(-c("year", "strata"), names_to = "quartile")
    d3.q <- t(d3.q) |>
      as.data.frame() |>
      cbind(year, strata = "CD4 350-500") |>
      tidyr::pivot_longer(-c("year", "strata"), names_to = "quartile")
    d4.q <- t(d4.q) |>
      as.data.frame() |>
      cbind(year, strata = "CD4 200-350") |>
      tidyr::pivot_longer(-c("year", "strata"), names_to = "quartile")
    d5.q <- t(d5.q) |>
      as.data.frame() |>
      cbind(year, strata = "CD4<200") |>
      tidyr::pivot_longer(-c("year", "strata"), names_to = "quartile")
    df <- rbind(d1.q, d2.q, d3.q, d4.q, d5.q)
  } else {
    d1.q <- t(d1.q) |>
      as.data.frame() |>
      cbind(year, strata = "CD4>500") |>
      tidyr::pivot_longer(-c("year", "strata"), names_to = "quartile")
    d2.q <- t(d2.q) |>
      as.data.frame() |>
      cbind(year, strata = "CD4 350-500") |>
      tidyr::pivot_longer(-c("year", "strata"), names_to = "quartile")
    d3.q <- t(d3.q) |>
      as.data.frame() |>
      cbind(year, strata = "CD4 200-350") |>
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

  q.diagnosed <- apply(
    diagnosed.HIV + diagnosed.AIDS,
    1,
    stats::quantile,
    probs = c(0.025, 0.5, 0.975)
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
    stats::quantile,
    probs = c(0.025, 0.5, 0.975)
  )
  m2.q <- apply(
    fit_mat[, m2.ind],
    2,
    stats::quantile,
    probs = c(0.025, 0.5, 0.975)
  )
  m3.q <- apply(
    fit_mat[, m3.ind],
    2,
    stats::quantile,
    probs = c(0.025, 0.5, 0.975)
  )
  m4.q <- apply(
    fit_mat[, m4.ind],
    2,
    stats::quantile,
    probs = c(0.025, 0.5, 0.975)
  )
  m5.q <- apply(
    fit_mat[, m5.ind],
    2,
    stats::quantile,
    probs = c(0.025, 0.5, 0.975)
  )
  year <- seq(start.yr, start.yr + nyr - 0.25, by = 0.25)

  m1.q <- t(m1.q) |>
    as.data.frame() |>
    cbind(year, strata = "CD4>500") |>
    tidyr::pivot_longer(-c("year", "strata"), names_to = "quartile")
  m2.q <- t(m2.q) |>
    as.data.frame() |>
    cbind(year, strata = "CD4 350-500") |>
    tidyr::pivot_longer(-c("year", "strata"), names_to = "quartile")
  m3.q <- t(m3.q) |>
    as.data.frame() |>
    cbind(year, strata = "CD4 200-350") |>
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

#' Prop UK dataframe
#'
#' Calculates proportion UK dataframe from posterior samples.
#'
#' @param stan_data Data list used in stan model fitting.
#' @param fit_mat Matrix of posterior samples from stan model.
#' @param start.yr Starting year for the dataframe.
#' @return Dataframe of proportion UK estimates.
ratio_abroad_uk.df <- function(stan_data, fit_mat, start.yr = 1978) {
  nyr <- stan_data$nquar / 4

  puk.ind <- grep("^p\\[", colnames(fit_mat))

  ### Quantiles of diagnosis
  if (length(puk.ind) == 0) {
    puk.q <- stan_data$p
  } else {
    puk.q <- apply(
      fit_mat[, puk.ind],
      2,
      stats::quantile,
      probs = c(0.025, 0.5, 0.975)
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
    detect.q <- stan_data$detect_prob
  } else {
    detect.q <- apply(
      fit_mat[, detect.ind],
      2,
      stats::quantile,
      probs = c(0.025, 0.5, 0.975)
    )
  }

  year <- seq(start.yr, start.yr + nyr - 0.25, by = 0.25)

  df <- t(detect.q) |>
    as.data.frame() |>
    cbind(year) |>
    tidyr::pivot_longer(-"year", names_to = "quartile")

  return(df)
}
