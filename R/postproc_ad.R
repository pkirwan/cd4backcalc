#' Post-process the age-dependent model fit
#'
#' Post-processes the fitted age-dependent back-calculation model to extract
#' various summaries and diagnostics.
#'
#' @param fit Fitted stan model object.
#' @param stan_data Data list used in stan model fitting.
#' @param start_yr Starting year for posterior predictive quantities.
#' @param nworkers Number of workers for parallel processing.
#' @param chain_results Logical indicating if chain-specific results are required.
#' @return List of post-processed results including summaries and diagnostics.
postproc_ad <- function(
  fit,
  stan_data,
  start_yr = 0,
  nworkers = 1,
  chain_results = FALSE
) {
  ### MODEL SETUP ###
  out <- list()
  ind_list <- list()
  N.Poiss <- 300

  # Useful constants
  config <- stan_data$config
  nquar <- stan_data$nquar
  if (stan_data$nquar %% 4 != 0) {
    stop("Number of quarters must be divisible by 4")
  }
  nage <- stan_data$nage

  # Define prevalent states
  prev_states <- c("all", 1:4)
  if (config$rita) {
    prev_states <- c(prev_states, 5:7)
  }

  # Define diagnosed states
  diag_states <- list(
    `CD4>=500` = 1,
    `CD4 350-499` = 2,
    `CD4 200-349` = 3,
    `CD4<200` = 4
  )
  if (config$rita) {
    diag_states <- list(
      `Recently acquired` = 1,
      `CD4>=500` = 2,
      `CD4 350-499` = 3,
      `CD4 200-349` = 4,
      `CD4<200` = 5
    )
  }

  # Define age bands
  age_bands <- list(
    all = c(15, 66),
    `1524` = c(15, 24),
    `2534` = c(25, 34),
    `3544` = c(35, 44),
    `45` = c(45, 66)
  )

  # Define quartiles
  q_levels <- c("2.5%", "50%", "mean", "97.5%")

  ### Fit summaries ###
  fit.mat <- as.matrix(fit)
  fit.summ <- rstan::summary(fit)[[1]]

  # Flags for the model considered
  ageprobs.flag <- ifelse(
    length(grep("^d\\[", colnames(fit.mat))) == nquar * 4,
    FALSE,
    TRUE
  )
  ageprobs.spl.flag <- ifelse(
    length(grep("^d\\[", colnames(fit.mat))) == nquar * 4 * nage,
    TRUE,
    FALSE
  )
  urep.flag <- ifelse(
    ("under_rep" %in% colnames(fit.mat)) == TRUE,
    TRUE,
    FALSE
  )

  ## Lack of convergence and sampling diagnostics
  out$conv.chk <- length(fit.summ[, 10][fit.summ[, 10] > 1.05]) == 0
  sampler.params <- rstan::get_sampler_params(fit, inc_warmup = FALSE)
  out$n.divg <- sum(sapply(sampler.params, function(x) x[, 5]))
  out$time <- max(apply(rstan::get_elapsed_time(fit), 1, sum))

  ### INCIDENCE ESTIMATES ###

  # Splines - beta: spline pars, gamma = X*beta
  beta.ind <- grep("beta_inf", colnames(fit.mat))
  tau.ind <- grep("tau", colnames(fit.mat))
  x <- stan_data$X_inf

  beta_draws <- fit.mat[, beta.ind, drop = FALSE]
  if (!is.matrix(x)) {
    x <- as.matrix(x)
  }
  infs <- beta_draws %*% t(x)

  if (length(tau.ind) > 1) {
    tau <- apply(fit.mat[, tau.ind], 2, get_quantiles)
  } else {
    tau <- get_quantiles(fit.mat[, tau.ind])
  }

  # Posterior means of parameter of interest
  out$tau <- tau
  alpha.ind <- grep("alpha", colnames(fit.mat))
  sigma2.ind <- grep("vardelta", colnames(fit.mat))

  if (length(alpha.ind) != 0) {
    out$alpha <- fit.mat[, alpha.ind]
  }
  if (length(sigma2.ind) != 0) {
    out$sigma2 <- fit.mat[, sigma2.ind]
  }

  # Posterior mean for all time and ages
  out$infs.all <- apply(infs, 2, get_quantiles)

  # Extract incidence estimates
  out$incidence <- infs.fct.qt(
    infs,
    age.st = age_bands$all[1],
    age.end = age_bands$all[2],
    nquar,
    nage = nage,
    start_yr = start_yr
  )

  out$annual_incidence <- infs.fct(
    infs,
    age.st = age_bands$all[1],
    age.end = age_bands$all[2],
    nquar,
    type = "time",
    nage = nage,
    start_yr = start_yr
  )

  # Extract incidence estimates by age band
  out$incidence_age <- purrr::imap(age_bands, \(age, bn) {
    infs.fct.qt(
      infs,
      age[1],
      age[2],
      nquar,
      nage = nage,
      start_yr = start_yr
    ) |>
      dplyr::mutate(age_band = factor(bn, levels = names(age_bands)))
  }) |>
    purrr::list_rbind()

  out$annual_incidence_age <- purrr::imap(age_bands, \(age, bn) {
    infs.fct(
      infs,
      age[1],
      age[2],
      nquar,
      type = "time",
      nage = nage,
      start_yr = start_yr
    ) |>
      dplyr::mutate(age_band = factor(bn, levels = names(age_bands)))
  }) |>
    purrr::list_rbind()

  ### DIAGNOSIS PROBABILITIES ###

  # Indices for diagnosis parameters by age band
  ind_list$d1.ind.1524 <- grep("^d\\[1,1", colnames(fit.mat))
  ind_list$d2.ind.1524 <- grep("^d\\[1,2", colnames(fit.mat))
  ind_list$d3.ind.1524 <- grep("^d\\[1,3", colnames(fit.mat))
  ind_list$d4.ind.1524 <- grep("^d\\[1,4", colnames(fit.mat))
  ind_list$d1.ind.2534 <- grep("^d\\[2,1", colnames(fit.mat))
  ind_list$d2.ind.2534 <- grep("^d\\[2,2", colnames(fit.mat))
  ind_list$d3.ind.2534 <- grep("^d\\[2,3", colnames(fit.mat))
  ind_list$d4.ind.2534 <- grep("^d\\[2,4", colnames(fit.mat))
  ind_list$d1.ind.3544 <- grep("^d\\[3,1", colnames(fit.mat))
  ind_list$d2.ind.3544 <- grep("^d\\[3,2", colnames(fit.mat))
  ind_list$d3.ind.3544 <- grep("^d\\[3,3", colnames(fit.mat))
  ind_list$d4.ind.3544 <- grep("^d\\[3,4", colnames(fit.mat))
  ind_list$d1.ind.45 <- grep("^d\\[4,1", colnames(fit.mat))
  ind_list$d2.ind.45 <- grep("^d\\[4,2", colnames(fit.mat))
  ind_list$d3.ind.45 <- grep("^d\\[4,3", colnames(fit.mat))
  ind_list$d4.ind.45 <- grep("^d\\[4,4", colnames(fit.mat))
  if (config$rita) {
    ind_list$d5.ind.1524 <- grep("^d\\[1,5", colnames(fit.mat))
    ind_list$d5.ind.2534 <- grep("^d\\[2,5", colnames(fit.mat))
    ind_list$d5.ind.3544 <- grep("^d\\[3,5", colnames(fit.mat))
    ind_list$d5.ind.45 <- grep("^d\\[4,5", colnames(fit.mat))
  }

  # Extract diagnosis probabilities by age band
  diag_grid <- tidyr::crossing(
    age_band = names(age_bands)[2:5],
    state = 1:length(diag_states)
  )

  out$diag_prob_age <- purrr::pmap(
    list(age_band = diag_grid$age_band, state = diag_grid$state),
    function(age_band, state) {
      idx_name <- sprintf("d%i.ind.%s", state, age_band)
      idx <- ind_list[[idx_name]]

      q <- apply(fit.mat[, idx, drop = FALSE], 2, get_quantiles)

      tibble::as_tibble(t(q)) |>
        dplyr::mutate(
          age_band = factor(age_band, levels = names(age_bands)),
          strata = factor(
            state,
            levels = 1:length(diag_states),
            labels = names(diag_states)
          ),
          year = seq(
            from = start_yr,
            length.out = nquar,
            by = 0.25
          )
        )
    }
  ) |>
    purrr::list_rbind() |>
    tidyr::pivot_longer(
      cols = dplyr::all_of(q_levels),
      names_to = "quartile",
      values_to = "value"
    )

  # Delta parameters (age progression rates)
  delta1.ind <- grep("delta_1", colnames(fit.mat))
  delta2.ind <- grep("delta_2", colnames(fit.mat))
  delta3.ind <- grep("delta_3", colnames(fit.mat))
  delta4.ind <- grep("delta_4", colnames(fit.mat))

  out$delta1.q <- apply(fit.mat[, delta1.ind], 2, get_quantiles)
  out$delta2.q <- apply(fit.mat[, delta2.ind], 2, get_quantiles)
  out$delta3.q <- apply(fit.mat[, delta3.ind], 2, get_quantiles)
  out$delta4.q <- apply(fit.mat[, delta4.ind], 2, get_quantiles)

  if (config$rita) {
    delta5.ind <- grep("delta_5", colnames(fit.mat))
    out$delta5.q <- apply(fit.mat[, delta5.ind], 2, get_quantiles)
  }

  ### GENERATE PREVALENCE AND DIAGNOSIS ARRAYS ###

  # Choose the correct Cpp functions
  if (config$rita) {
    g.fct <- gof_iter_fct_rita
    p.fct <- prev_iter_fct_rita
  } else {
    g.fct <- gof_iter_fct
    p.fct <- prev_iter_fct
  }

  # Convert q and init_prev to matrices for easier handling in Cpp
  q_mat <- matrix(stan_data$q_age, nrow = 4, ncol = nage)
  prev_mat <- matrix(
    stan_data$init_prev_age,
    nrow = nrow(stan_data$init_prev_age),
    ncol = nage
  )

  # Expected diagnoses
  exp.diags <- lapply(1:nrow(fit.mat), function(i) {
    g.fct(
      exp(infs[i, ]),
      d.mat.iter(
        fit.mat,
        i,
        age.dx.flag = ageprobs.flag,
        age.dx.spl.flag = ageprobs.spl.flag,
        rita.flag = config$rita,
        nt = nquar,
        nage = nage,
        ind_list = ind_list
      ),
      q_mat,
      prev_mat,
      nquar,
      nage
    )
  })

  # Expected prevalence
  exp.prev <- lapply(1:nrow(fit.mat), function(i) {
    p.fct(
      exp(infs[i, ]),
      d.mat.iter(
        fit.mat,
        i,
        age.dx.flag = ageprobs.flag,
        age.dx.spl.flag = ageprobs.spl.flag,
        rita.flag = config$rita,
        nt = nquar,
        nage = nage,
        ind_list = ind_list
      ),
      q_mat,
      prev_mat,
      nquar,
      nage
    )
  })

  # List to array (easier to handle)
  exp.diags.arr <- simplify2array(exp.diags)
  exp.prev.arr <- simplify2array(exp.prev)

  ### EXPECTED DIAGNOSES ###

  # Expected number of CD4
  nCD4 <- apply(stan_data$CD4_age, c(1, 2), sum)

  # denom dims: (nt, nage, nsim)
  cd4_denom <- apply(
    exp.diags.arr[,, seq_len(config$ndiagst - 1), , drop = FALSE],
    c(1, 2, 4),
    sum
  )

  # stage proportions p(stage | diagnosed), dims: (nt, nage, stage, nsim)
  cd4_num <- exp.diags.arr[,, seq_len(config$ndiagst - 1), , drop = FALSE]
  cd4_p_all <- sweep(cd4_num, MARGIN = c(1, 2, 4), STATS = cd4_denom, FUN = "/")
  cd4_p_all[!is.finite(cd4_p_all)] <- 0

  # expected CD4 counts by stage
  CD4_all <- sweep(cd4_p_all, MARGIN = c(1, 2), STATS = nCD4, FUN = "*")

  # extract stages into the original variables expected later
  CD4.st1 <- CD4_all[,, 1, ]
  CD4.st2 <- CD4_all[,, 2, ]
  CD4.st3 <- CD4_all[,, 3, ]
  CD4.st4 <- CD4_all[,, 4, ]
  if (isTRUE(config$rita)) {
    CD4.st5 <- CD4_all[,, 5, ]
  } else {
    CD4.st5 <- NULL
  }

  ### Expected HIV diagnoses
  if (config$rita) {
    HIV <- exp.diags.arr[,, 1, ] +
      exp.diags.arr[,, 2, ] +
      exp.diags.arr[,, 3, ] +
      exp.diags.arr[,, 4, ] +
      exp.diags.arr[,, 5, ]
  } else {
    HIV <- exp.diags.arr[,, 1, ] +
      exp.diags.arr[,, 2, ] +
      exp.diags.arr[,, 3, ] +
      exp.diags.arr[,, 4, ]
  }

  ### Expected AIDS diagnoses
  ### Depending on whether under-reporting present
  ### (With real-data always)
  aids_slice <- ifelse(config$rita, 6, 5)
  if (urep.flag) {
    urep.vec <- fit.mat[, "under_rep"]
    AIDS <- simplify2array(lapply(1:nrow(fit.mat), function(i) {
      urep.fct(
        urep.vec[i],
        exp.diags.arr[,, aids_slice, i],
        urep.flag = urep.flag,
        nt = nquar
      )
    }))
  } else {
    AIDS <- exp.diags.arr[,, aids_slice, ]
  }

  ### GOODNESS OF FIT ###
  gof_types <- list(
    HIV = "HIV",
    AIDS = "AIDS",
    `CD4.1` = "CD4.st1",
    `CD4.2` = "CD4.st2",
    `CD4.3` = "CD4.st3",
    `CD4.4` = "CD4.st4"
  )
  if (config$rita) {
    gof_types[["CD4.5"]] <- "CD4.st5"
  }

  grid <- tidyr::crossing(
    base_name = names(gof_types),
    bn = names(age_bands)
  )

  gof_results <- purrr::pmap(
    list(base_name = grid$base_name, bn = grid$bn),
    function(base_name, bn) {
      type_arg <- gof_types[[base_name]]
      age <- age_bands[[bn]]

      gof.fct(
        type_arg,
        age[1],
        age[2],
        HIV,
        AIDS,
        CD4.st1,
        CD4.st2,
        CD4.st3,
        CD4.st4,
        CD4.st5,
        nquar
      )
    }
  )

  names(gof_results) <- paste0("gof.", grid$base_name, grid$bn)
  out <- c(out, gof_results)

  ### UNDIAGNOSED PREVALENCE ###

  # undiag_prev overall
  out$undiag_prev <- prev.fct(
    exp.prev.arr,
    "all",
    age.st = age_bands$all[1],
    age.end = age_bands$all[2],
    config$rita,
    nquar,
    start_yr
  )

  # undiag_prev by state
  out$undiag_prev_state <- purrr::map(prev_states, \(state) {
    prev.fct(
      exp.prev.arr,
      state,
      age.st = age_bands$all[1],
      age.end = age_bands$all[2],
      config$rita,
      nquar,
      start_yr
    ) |>
      dplyr::mutate(
        strata = factor(
          as.character(state),
          levels = as.character(prev_states)
        )
      )
  }) |>
    purrr::list_rbind()

  # undiag_prev by age
  out$undiag_prev_age <- purrr::imap(age_bands, \(age, bn) {
    prev.fct(
      exp.prev.arr,
      "all",
      age[1],
      age[2],
      config$rita,
      nquar,
      start_yr
    ) |>
      dplyr::mutate(age_band = factor(bn, levels = names(age_bands)))
  }) |>
    purrr::list_rbind()

  # undiag_prev by age and state
  # build all combinations and map once
  prev_grid <- tidyr::expand_grid(
    state = prev_states,
    bn = names(age_bands)
  )

  out$undiag_prev_age_state <- purrr::pmap(
    list(state = prev_grid$state, bn = prev_grid$bn),
    \(state, bn) {
      age <- age_bands[[bn]]
      prev.fct(
        exp.prev.arr,
        state,
        age[1],
        age[2],
        config$rita,
        nquar,
        start_yr
      ) |>
        dplyr::mutate(
          age_band = factor(bn, levels = names(age_bands)),
          strata = factor(
            as.character(state),
            levels = as.character(prev_states)
          )
        )
    }
  ) |>
    purrr::list_rbind()

  ### POSTERIOR PREDICTIVE DISTRIBUTIONS ###

  HIV.m <- apply(HIV, c(1, 2), mean)

  HIV.p <- apply(HIV.m, c(1, 2), function(x) stats::rpois(N.Poiss, lambda = x))

  AIDS.m <- apply(AIDS, c(1, 2), mean)

  AIDS.p <- apply(AIDS.m, c(1, 2), function(x) {
    stats::rpois(N.Poiss, lambda = x)
  })

  # Probability of diagnosis by stage

  # denom dims: (nt, nage, nsim)
  denom <- apply(
    exp.diags.arr[,, seq_len(config$ndiagst - 1), , drop = FALSE],
    c(1, 2, 4),
    sum
  )

  # p_all dims: (nt, nage, stage, nsim)
  p_num <- exp.diags.arr[,, seq_len(config$ndiagst - 1), , drop = FALSE]
  p_all <- sweep(p_num, MARGIN = c(1, 2, 4), STATS = denom, FUN = "/")
  p_all[!is.finite(p_all)] <- 0

  # extract each stage
  p1 <- p_all[,, 1, ]
  p2 <- p_all[,, 2, ]
  p3 <- p_all[,, 3, ]
  p4 <- p_all[,, 4, ]
  if (config$rita) {
    p5 <- p_all[,, 5, ]
  }

  p1.m <- apply(p1, c(1, 2), mean)
  p2.m <- apply(p2, c(1, 2), mean)
  p3.m <- apply(p3, c(1, 2), mean)
  p4.m <- apply(p4, c(1, 2), mean)
  if (config$rita) {
    p5.m <- apply(p5, c(1, 2), mean)
  }

  if (config$rita) {
    CD4.Poiss <- array(0, dim = c(nquar, nage, 5, N.Poiss))
  } else {
    CD4.Poiss <- array(0, dim = c(nquar, nage, 4, N.Poiss))
  }

  if (config$rita) {
    ### if quarterly there cannot be CD4 for a=1, at t=5, t=9,...
    ### as anyway nCD4[,1] = 0 then ignore and put 0
    for (t in 1:nquar) {
      for (a in 2:nage) {
        CD4.Poiss[t, a, , ] <- stats::rmultinom(
          N.Poiss,
          size = nCD4[t, a],
          prob = c(p1.m[t, a], p2.m[t, a], p3.m[t, a], p4.m[t, a], p5.m[t, a])
        )
      }
    }
  } else {
    for (t in 1:nquar) {
      for (a in 2:nage) {
        CD4.Poiss[t, a, , ] <- stats::rmultinom(
          N.Poiss,
          size = nCD4[t, a],
          prob = c(p1.m[t, a], p2.m[t, a], p3.m[t, a], p4.m[t, a])
        )
      }
    }
  }

  # CD4 total over stages: dims (nt, nage, N.Poiss)
  CD4.total <- apply(CD4.Poiss, c(1, 2, 4), sum)

  # Create tidy data frames for posterior predictive summaries
  out$HIV_ann_age <- tidy_bind(
    HIV.p,
    age_bands,
    q_levels,
    nquar,
    start_yr,
    qt.flag = FALSE
  )
  out$AIDS_ann_age <- tidy_bind(
    AIDS.p,
    age_bands,
    q_levels,
    nquar,
    start_yr,
    qt.flag = FALSE
  )
  out$CD4_ann_age <- tidy_bind(
    CD4.total,
    age_bands,
    q_levels,
    nquar,
    start_yr,
    qt.flag = FALSE,
    CD4 = TRUE
  )

  out$HIV_age <- tidy_bind(HIV.p, age_bands, q_levels, nquar, start_yr)
  out$AIDS_age <- tidy_bind(AIDS.p, age_bands, q_levels, nquar, start_yr)
  out$CD4_age <- tidy_bind(
    CD4.total,
    age_bands,
    q_levels,
    nquar,
    start_yr,
    CD4 = TRUE
  )

  # Combine diagnoses
  out$diagnoses_ann_age <- dplyr::left_join(
    out$HIV_ann_age,
    out$AIDS_ann_age,
    by = c("age_band", "year", "quartile")
  ) |>
    dplyr::mutate(value = .data$value.x + .data$value.y) |>
    dplyr::select(-c("value.x", "value.y"))

  out$diagnoses_age <- dplyr::left_join(
    out$HIV_age,
    out$AIDS_age,
    by = c("age_band", "year", "quartile")
  ) |>
    dplyr::mutate(value = .data$value.x + .data$value.y) |>
    dplyr::select(-c("value.x", "value.y"))

  out$diagnoses <- out$diagnoses_age |> dplyr::filter(.data$age_band == "all")

  # ### Time to diagnosis stuff
  # if(!ageprobs.flag & !quar.flag & !ageprobs.spl.flag){
  #   stop("Time to diagnosis not coded up yet")
  # }else if(ageprobs.flag & !quar.flag & !ageprobs.spl.flag){
  #   stop("Time to diagnosis not coded up yet")
  # }else if(!ageprobs.flag & quar.flag & !ageprobs.spl.flag){
  #   Rcpp::sourceCpp('HIV backcalculation/Age dependency/MSMdata16/TimeToDxFct_aidiags.cpp')
  # }else if(ageprobs.flag & quar.flag & !ageprobs.spl.flag){
  #   Rcpp::sourceCpp('HIV backcalculation/Age dependency/MSMdata16/TimeToDxFct_acdiags.cpp')
  # }else if(ageprobs.flag & quar.flag & ageprobs.spl.flag){
  #   Rcpp::sourceCpp('HIV backcalculation/Age dependency/MSMdata16/TimeToDxFct_adiags.cpp')
  # }else{
  #   stop("Flags somehow uncorrectly specified")
  # }
  # d.ind <- grep("d\\[",colnames(fit.mat))
  # d.mat <- fit.mat[,d.ind]
  #
  # ### the call depends on model... maybe include nac in fct
  # t.to.d <- t_to_d(nt, nage, 4, n_ind_sim=20, d.mat, stan_data$q)
  # snaps <- snap(nt, nage, 4, n_ind_sim=20, t_max=150,  d.mat, stan_data$q)
  #
  # t.to.d.dist <- apply(t.to.d[[1]], c(1,2), quantile, probs=c(0.1,0.25,0.5,0.75,0.9))
  #
  # t.to.d.dist[t.to.d.dist==(nt+1)] <- NA
  # ### Scaling quarters to years
  # out$t.to.d.dist  <- t.to.d.dist/4
  #
  # snaps.dist <- apply(snaps[[1]], c(1,2), quantile, probs=c(0.1,0.25,0.5,0.75,0.9))
  #
  # ### Scaling quarters to years
  # out$snaps.dist  <- snaps.dist/4
  #

  return(out)
}

#' Helper function to bind tidy data frames
#'
#' @param list List of HIV posterior predictive matrices.
#' @param age_bands Age bands list.
#' @param q_levels Quartile levels.
#' @param nquar Number of quarters.
#' @param start_yr Starting year.
#' @param qt.flag Logical indicating if quarterly output is required.
#' @param CD4 Logical indicating if CD4 data is being processed.
#' @return Combined tidy data frame.
tidy_bind <- function(
  list,
  age_bands,
  q_levels,
  nquar,
  start_yr,
  qt.flag = TRUE,
  CD4 = FALSE
) {
  out_template <- tibble::tibble(
    age_band = factor(character(0), levels = names(age_bands)),
    year = numeric(0),
    quartile = factor(character(0), levels = q_levels),
    value = numeric(0)
  )

  # Create grouping variable for annual sums
  grp_year <- rep(seq_len(nquar / 4), each = 4)

  # Define years vector
  if (qt.flag) {
    yrs <- start_yr + (seq_len(nquar) - 1) / 4
  } else {
    yrs <- start_yr + seq_len(nquar / 4) - 1
  }

  # Loop over age bands
  for (bn in names(age_bands)) {
    age <- age_bands[[bn]]
    age_range <- (age[1]:age[2]) - 14

    if (qt.flag) {
      if (CD4) {
        tmp <- apply(list[, age_range, , drop = FALSE], 3, function(m) {
          rowSums(m)
        })
      } else {
        tmp <- apply(list[,, age_range, drop = FALSE], 1, function(m) {
          rowSums(m)
        })
      }
    } else {
      if (CD4) {
        tmp <- apply(list[, age_range, , drop = FALSE], 3, function(m) {
          as.numeric(rowsum(rowSums(m), grp_year, reorder = FALSE))
        })
      } else {
        tmp <- apply(list[,, age_range, drop = FALSE], 1, function(m) {
          as.numeric(rowsum(rowSums(m), grp_year, reorder = FALSE))
        })
      }
    }

    q <- apply(tmp, 1, get_quantiles)

    df1 <- tibble::as_tibble(t(q)) |>
      dplyr::mutate(
        age_band = factor(bn),
        year = yrs
      ) |>
      tidyr::pivot_longer(
        cols = -c("age_band", "year"),
        names_to = "quartile",
        values_to = "value"
      ) |>
      dplyr::mutate(
        quartile = factor(.data$quartile)
      )
    out_template <- dplyr::bind_rows(out_template, df1)
  }
  return(out_template)
}

#' Get initial prevalence matrix to provide to age-dependent model
#'
#' @param config Configuration list.
#' @param exp.prev.arr Expected prevalence array.
#' @param start_yr Starting year.
#' @return Initial prevalence matrix.
get_init_prev <- function(config, exp.prev.arr, start_yr) {
  ### Initial prevalence matrix
  # Prev.mat for specific years
  if (config$rita) {
    prev.mat.2011 <- rbind(
      apply(exp.prev.arr[1, , 1, ], 1, get_quantiles)[2, ],
      apply(exp.prev.arr[1, , 2, ], 1, get_quantiles)[2, ],
      apply(exp.prev.arr[1, , 3, ], 1, get_quantiles)[2, ],
      apply(exp.prev.arr[1, , 4, ], 1, get_quantiles)[2, ],
      apply(exp.prev.arr[1, , 5, ], 1, get_quantiles)[2, ],
      apply(exp.prev.arr[1, , 6, ], 1, get_quantiles)[2, ],
      apply(exp.prev.arr[1, , 7, ], 1, get_quantiles)[2, ]
    )
  } else if (start_yr == 1995) {
    prev.mat.1995 <- rbind(
      apply(exp.prev.arr[1, , 1, ], 1, get_quantiles)[2, ],
      apply(exp.prev.arr[1, , 2, ], 1, get_quantiles)[2, ],
      apply(exp.prev.arr[1, , 3, ], 1, get_quantiles)[2, ],
      apply(exp.prev.arr[1, , 4, ], 1, get_quantiles)[2, ]
    )
    prev.mat.2011 <- rbind(
      apply(exp.prev.arr[65, , 1, ], 1, get_quantiles)[2, ],
      apply(exp.prev.arr[65, , 2, ], 1, get_quantiles)[2, ],
      apply(exp.prev.arr[65, , 3, ], 1, get_quantiles)[2, ],
      apply(exp.prev.arr[65, , 4, ], 1, get_quantiles)[2, ]
    )
  } else if (start_yr == 2011) {
    prev.mat.2011 <- rbind(
      apply(exp.prev.arr[1, , 1, ], 1, get_quantiles)[2, ],
      apply(exp.prev.arr[1, , 2, ], 1, get_quantiles)[2, ],
      apply(exp.prev.arr[1, , 3, ], 1, get_quantiles)[2, ],
      apply(exp.prev.arr[1, , 4, ], 1, get_quantiles)[2, ]
    )
  }
}
