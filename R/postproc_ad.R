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
  out <- list()

  ## Useful constants
  config <- stan_data$config

  ## nt is the the number of period considered
  ## i.e. nquar if model quarterly
  ## or nyr if model yearly
  ## is.null allows distinction between quarter and years
  nt <- stan_data$nquar
  nage <- stan_data$nage

  ### To store indices of various parameters
  ind_list <- list()

  ### Specifications
  # out$data <- 1 # tentatively set data.rm <- 0
  out$nt <- nt
  # out$data.id <- data.id

  ### Fit summaries
  fit.mat <- as.matrix(fit)
  fit.summ <- rstan::summary(fit)[[1]]

  ### Flags for the model considered
  ## Is this quarterly? Yearly?
  ## AgeDx Probs? or not
  quar.flag <- ifelse(is.null(stan_data$nquar) == FALSE, TRUE, FALSE)
  ageprobs.flag <- ifelse(
    length(grep("^d\\[", colnames(fit.mat))) == nt * 4,
    FALSE,
    TRUE
  )
  ageprobs.spl.flag <- ifelse(
    length(grep("^d\\[", colnames(fit.mat))) == nt * 4 * nage,
    TRUE,
    FALSE
  )
  delta.flag <- ifelse(
    length(grep("delta1", colnames(fit.mat))) > 0,
    TRUE,
    FALSE
  )
  urep.flag <- ifelse(("under_rep" %in% colnames(fit.mat)) == TRUE, TRUE, FALSE)

  out$quar.flag <- quar.flag
  out$ageprobs.flag <- ageprobs.flag
  out$ageprobs.spl.flag <- ageprobs.spl.flag
  out$urep.flag <- urep.flag

  ### General information
  ## Lack of convergence
  out$conv.chk <- length(fit.summ[, 10][fit.summ[, 10] > 1.05]) == 0

  ## Divergent stuff
  sampler.params <- rstan::get_sampler_params(fit, inc_warmup = FALSE)
  out$n.divg <- sum(sapply(sampler.params, function(x) x[, 5]))

  ## Time taken
  out$time <- max(apply(rstan::get_elapsed_time(fit), 1, sum))

  ## Infections (yearly)
  ## Splines - beta: spline pars, gamma = X*beta
  ## infs matrix with log-infs for every simulation
  ## has same structure for GP and splines

  beta.ind <- grep("beta", colnames(fit.mat))
  lambda.ind <- grep("lambda", colnames(fit.mat))
  x <- stan_data$X

  beta_draws <- fit.mat[, beta.ind, drop = FALSE]
  if (!is.matrix(x)) {
    x <- as.matrix(x)
  }
  infs <- beta_draws %*% t(x)

  if (length(lambda.ind) > 1) {
    lambda <- apply(fit.mat[, lambda.ind], 2, get_quantiles)
  } else {
    ## if length==1 not a matrix thus does not work
    lambda <- get_quantiles(fit.mat[, lambda.ind])
  }

  ### Posterior mean of smoothing parameter of interest
  out$lambda <- lambda
  ### Posterior mean of alpha and sigma coefficients of interest
  alpha.ind <- grep("alpha", colnames(fit.mat))
  sigma2.ind <- grep("vardelta", colnames(fit.mat))

  if (length(alpha.ind) != 0) {
    out$alpha <- fit.mat[, alpha.ind]
  }
  if (length(sigma2.ind) != 0) {
    out$sigma2 <- fit.mat[, sigma2.ind]
  }

  ## Posterior mean for all time and ages
  out$infs.all <- apply(infs, 2, get_quantiles)

  ## Define age bands
  age_bands <- list(
    all = c(15, 66),
    `1524` = c(15, 24),
    `2534` = c(25, 34),
    `3544` = c(35, 44),
    `45` = c(45, 66)
  )

  band_suffix <- function(bn) {
    if (identical(bn, "all")) "" else paste0(".", bn)
  }

  for (bn in names(age_bands)) {
    ages <- age_bands[[bn]]
    suf <- band_suffix(bn)

    out[[paste0("infs", suf)]] <- infs.fct(
      infs,
      ages[1],
      ages[2],
      nt,
      type = "time",
      qt.flag = quar.flag,
      nage = nage
    )

    out[[paste0("infs.qt", suf)]] <- infs.fct.qt(
      infs,
      ages[1],
      ages[2],
      nt,
      qt.flag = quar.flag,
      nage = nage
    )

    out[[paste0("peak", suf)]] <- peak.fct(
      infs,
      ages[1],
      ages[2],
      nt,
      last.yr = 10,
      qt.flag = quar.flag,
      nage = nage
    )

    for (last_yr in c(2, 3)) {
      out[[paste0("trend.", last_yr, suf)]] <- trend.fct(
        infs,
        ages[1],
        ages[2],
        nt,
        last.yr = last_yr,
        qt.flag = quar.flag,
        nage = nage
      )
    }
  }

  ### Diagnoses probabilities
  if (!ageprobs.flag) {
    ind_list$d1.ind <- grep("^d\\[1,", colnames(fit.mat))
    ind_list$d2.ind <- grep("^d\\[2,", colnames(fit.mat))
    ind_list$d3.ind <- grep("^d\\[3,", colnames(fit.mat))
    ind_list$d4.ind <- grep("^d\\[4,", colnames(fit.mat))

    out$d1.q <- apply(fit.mat[, ind_list$d1.ind], 2, get_quantiles)
    out$d2.q <- apply(fit.mat[, ind_list$d2.ind], 2, get_quantiles)
    out$d3.q <- apply(fit.mat[, ind_list$d3.ind], 2, get_quantiles)
    out$d4.q <- apply(fit.mat[, ind_list$d4.ind], 2, get_quantiles)
  } else {
    if (!ageprobs.spl.flag) {
      ind_list$d1.ind.y <- grep("^d\\[1,1", colnames(fit.mat))
      ind_list$d2.ind.y <- grep("^d\\[1,2", colnames(fit.mat))
      ind_list$d3.ind.y <- grep("^d\\[1,3", colnames(fit.mat))
      ind_list$d4.ind.y <- grep("^d\\[1,4", colnames(fit.mat))
      ind_list$d1.ind.my <- grep("^d\\[2,1", colnames(fit.mat))
      ind_list$d2.ind.my <- grep("^d\\[2,2", colnames(fit.mat))
      ind_list$d3.ind.my <- grep("^d\\[2,3", colnames(fit.mat))
      ind_list$d4.ind.my <- grep("^d\\[2,4", colnames(fit.mat))
      ind_list$d1.ind.mo <- grep("^d\\[3,1", colnames(fit.mat))
      ind_list$d2.ind.mo <- grep("^d\\[3,2", colnames(fit.mat))
      ind_list$d3.ind.mo <- grep("^d\\[3,3", colnames(fit.mat))
      ind_list$d4.ind.mo <- grep("^d\\[3,4", colnames(fit.mat))
      ind_list$d1.ind.o <- grep("^d\\[4,1", colnames(fit.mat))
      ind_list$d2.ind.o <- grep("^d\\[4,2", colnames(fit.mat))
      ind_list$d3.ind.o <- grep("^d\\[4,3", colnames(fit.mat))
      ind_list$d4.ind.o <- grep("^d\\[4,4", colnames(fit.mat))

      out$d1.q.y <- apply(fit.mat[, ind_list$d1.ind.y], 2, get_quantiles)
      out$d2.q.y <- apply(fit.mat[, ind_list$d2.ind.y], 2, get_quantiles)
      out$d3.q.y <- apply(fit.mat[, ind_list$d3.ind.y], 2, get_quantiles)
      out$d4.q.y <- apply(fit.mat[, ind_list$d4.ind.y], 2, get_quantiles)
      out$d1.q.my <- apply(fit.mat[, ind_list$d1.ind.my], 2, get_quantiles)
      out$d2.q.my <- apply(fit.mat[, ind_list$d2.ind.my], 2, get_quantiles)
      out$d3.q.my <- apply(fit.mat[, ind_list$d3.ind.my], 2, get_quantiles)
      out$d4.q.my <- apply(fit.mat[, ind_list$d4.ind.my], 2, get_quantiles)
      out$d1.q.mo <- apply(fit.mat[, ind_list$d1.ind.mo], 2, get_quantiles)
      out$d2.q.mo <- apply(fit.mat[, ind_list$d2.ind.mo], 2, get_quantiles)
      out$d3.q.mo <- apply(fit.mat[, ind_list$d3.ind.mo], 2, get_quantiles)
      out$d4.q.mo <- apply(fit.mat[, ind_list$d4.ind.mo], 2, get_quantiles)
      out$d1.q.o <- apply(fit.mat[, ind_list$d1.ind.o], 2, get_quantiles)
      out$d2.q.o <- apply(fit.mat[, ind_list$d2.ind.o], 2, get_quantiles)
      out$d3.q.o <- apply(fit.mat[, ind_list$d3.ind.o], 2, get_quantiles)
      out$d4.q.o <- apply(fit.mat[, ind_list$d4.ind.o], 2, get_quantiles)

      if (config$rita) {
        ind_list$d5.ind.y <- grep("^d\\[1,5", colnames(fit.mat))
        ind_list$d5.ind.my <- grep("^d\\[2,5", colnames(fit.mat))
        ind_list$d5.ind.mo <- grep("^d\\[3,5", colnames(fit.mat))
        ind_list$d5.ind.o <- grep("^d\\[4,5", colnames(fit.mat))

        out$d5.q.y <- apply(fit.mat[, ind_list$d5.ind.y], 2, get_quantiles)
        out$d5.q.my <- apply(fit.mat[, ind_list$d5.ind.my], 2, get_quantiles)
        out$d5.q.mo <- apply(fit.mat[, ind_list$d5.ind.mo], 2, get_quantiles)
        out$d5.q.o <- apply(fit.mat[, ind_list$d5.ind.o], 2, get_quantiles)
      }
    } else {
      ind_list$d1.ind <- grep("^d\\[1,", colnames(fit.mat))
      ind_list$d2.ind <- grep("^d\\[2,", colnames(fit.mat))
      ind_list$d3.ind <- grep("^d\\[3,", colnames(fit.mat))
      ind_list$d4.ind <- grep("^d\\[4,", colnames(fit.mat))
      d1.q <- apply(fit.mat[, ind_list$d1.ind], 2, get_quantiles)
      d2.q <- apply(fit.mat[, ind_list$d2.ind], 2, get_quantiles)
      d3.q <- apply(fit.mat[, ind_list$d3.ind], 2, get_quantiles)
      d4.q <- apply(fit.mat[, ind_list$d4.ind], 2, get_quantiles)
      out$d1.q <- array(d1.q, dim = c(4, nt, nage)) ## array (c(2.5,50,mean,97.5),t,a)
      out$d2.q <- array(d2.q, dim = c(4, nt, nage)) ## array (c(2.5,50,mean,97.5),t,a)
      out$d3.q <- array(d3.q, dim = c(4, nt, nage)) ## array (c(2.5,50,mean,97.5),t,a)
      out$d4.q <- array(d4.q, dim = c(4, nt, nage)) ## array (c(2.5,50,mean,97.5),t,a)
    }
  }

  if (delta.flag) {
    delta1.ind <- grep("delta1", colnames(fit.mat))
    delta2.ind <- grep("delta2", colnames(fit.mat))
    delta3.ind <- grep("delta3", colnames(fit.mat))
    delta4.ind <- grep("delta4", colnames(fit.mat))

    out$delta1.q <- apply(fit.mat[, delta1.ind], 2, get_quantiles)
    out$delta2.q <- apply(fit.mat[, delta2.ind], 2, get_quantiles)
    out$delta3.q <- apply(fit.mat[, delta3.ind], 2, get_quantiles)
    out$delta4.q <- apply(fit.mat[, delta4.ind], 2, get_quantiles)

    if (config$rita) {
      delta5.ind <- grep("delta5", colnames(fit.mat))
      out$delta5.q <- apply(fit.mat[, delta5.ind], 2, get_quantiles)
    }
  } else {
    # deltas are the same for every age group so only derive once
    out$delta1.q <- apply(
      log((fit.mat[, ind_list$d1.ind.y]) / (1 - fit.mat[, ind_list$d1.ind.y])) -
        fit.mat[, alpha.ind][, 1],
      2,
      get_quantiles
    )
    out$delta2.q <- apply(
      log((fit.mat[, ind_list$d2.ind.y]) / (1 - fit.mat[, ind_list$d2.ind.y])) -
        fit.mat[, alpha.ind][, 2],
      2,
      get_quantiles
    )
    out$delta3.q <- apply(
      log((fit.mat[, ind_list$d3.ind.y]) / (1 - fit.mat[, ind_list$d3.ind.y])) -
        fit.mat[, alpha.ind][, 3],
      2,
      get_quantiles
    )
    out$delta4.q <- apply(
      log((fit.mat[, ind_list$d4.ind.y]) / (1 - fit.mat[, ind_list$d4.ind.y])) -
        fit.mat[, alpha.ind][, 4],
      2,
      get_quantiles
    )
  }

  ## Need to set environment of d.mat.iter
  ## to be global environment (otherwise too many arguments to take care of ...)
  #  environment(d.mat.iter) <- environment()

  ## GOODNESS OF FIT
  # Choose the correct Cpp functions
  if (config$rita) {
    g.fct <- gof_iter_fct_rita
    p.fct <- prev_iter_fct_rita
  } else {
    g.fct <- gof_iter_fct
    p.fct <- prev_iter_fct
  }

  # convert q and init_prev to matrices for easier handling in Cpp
  q_mat <- matrix(stan_data$q_age, nrow = 4, ncol = nage)
  prev_mat <- matrix(
    stan_data$init_prev_age,
    nrow = nrow(stan_data$init_prev_age),
    ncol = nage
  )

  exp.diags <- lapply(1:nrow(fit.mat), function(i) {
    g.fct(
      exp(infs[i, ]),
      d.mat.iter(
        fit.mat,
        i,
        age.dx.flag = ageprobs.flag,
        age.dx.spl.flag = ageprobs.spl.flag,
        rita.flag = config$rita,
        nt = nt,
        nage = nage,
        ind_list = ind_list
      ),
      q_mat,
      prev_mat,
      nt,
      nage
    )
  })
  exp.prev <- lapply(1:nrow(fit.mat), function(i) {
    p.fct(
      exp(infs[i, ]),
      d.mat.iter(
        fit.mat,
        i,
        age.dx.flag = ageprobs.flag,
        age.dx.spl.flag = ageprobs.spl.flag,
        rita.flag = config$rita,
        nt = nt,
        nage = nage,
        ind_list = ind_list
      ),
      q_mat,
      prev_mat,
      nt,
      nage
    )
  })

  ### List to array (easier to handle)
  exp.diags.arr <- simplify2array(exp.diags)
  exp.prev.arr <- simplify2array(exp.prev)

  ### Prev.mat

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

  ### Expected number of CD4
  nCD4 <- apply(stan_data$CD4_age, c(1, 2), sum)
  if (config$rita) {
    CD4.st1 <- simplify2array(lapply(1:nrow(fit.mat), function(i) {
      nCD4 *
        (exp.diags.arr[,, 1, i] /
          (exp.diags.arr[,, 1, i] +
            exp.diags.arr[,, 2, i] +
            exp.diags.arr[,, 3, i] +
            exp.diags.arr[,, 4, i] +
            exp.diags.arr[,, 5, i]))
    }))
    CD4.st2 <- simplify2array(lapply(1:nrow(fit.mat), function(i) {
      nCD4 *
        (exp.diags.arr[,, 2, i] /
          (exp.diags.arr[,, 1, i] +
            exp.diags.arr[,, 2, i] +
            exp.diags.arr[,, 3, i] +
            exp.diags.arr[,, 4, i] +
            exp.diags.arr[,, 5, i]))
    }))
    CD4.st3 <- simplify2array(lapply(1:nrow(fit.mat), function(i) {
      nCD4 *
        (exp.diags.arr[,, 3, i] /
          (exp.diags.arr[,, 1, i] +
            exp.diags.arr[,, 2, i] +
            exp.diags.arr[,, 3, i] +
            exp.diags.arr[,, 4, i] +
            exp.diags.arr[,, 5, i]))
    }))
    CD4.st4 <- simplify2array(lapply(1:nrow(fit.mat), function(i) {
      nCD4 *
        (exp.diags.arr[,, 4, i] /
          (exp.diags.arr[,, 1, i] +
            exp.diags.arr[,, 2, i] +
            exp.diags.arr[,, 3, i] +
            exp.diags.arr[,, 4, i] +
            exp.diags.arr[,, 5, i]))
    }))
    CD4.st5 <- simplify2array(lapply(1:nrow(fit.mat), function(i) {
      nCD4 *
        (exp.diags.arr[,, 5, i] /
          (exp.diags.arr[,, 1, i] +
            exp.diags.arr[,, 2, i] +
            exp.diags.arr[,, 3, i] +
            exp.diags.arr[,, 4, i] +
            exp.diags.arr[,, 5, i]))
    }))
  } else {
    CD4.st1 <- simplify2array(lapply(1:nrow(fit.mat), function(i) {
      nCD4 *
        (exp.diags.arr[,, 1, i] /
          (exp.diags.arr[,, 1, i] +
            exp.diags.arr[,, 2, i] +
            exp.diags.arr[,, 3, i] +
            exp.diags.arr[,, 4, i]))
    }))
    CD4.st2 <- simplify2array(lapply(1:nrow(fit.mat), function(i) {
      nCD4 *
        (exp.diags.arr[,, 2, i] /
          (exp.diags.arr[,, 1, i] +
            exp.diags.arr[,, 2, i] +
            exp.diags.arr[,, 3, i] +
            exp.diags.arr[,, 4, i]))
    }))
    CD4.st3 <- simplify2array(lapply(1:nrow(fit.mat), function(i) {
      nCD4 *
        (exp.diags.arr[,, 3, i] /
          (exp.diags.arr[,, 1, i] +
            exp.diags.arr[,, 2, i] +
            exp.diags.arr[,, 3, i] +
            exp.diags.arr[,, 4, i]))
    }))
    CD4.st4 <- simplify2array(lapply(1:nrow(fit.mat), function(i) {
      nCD4 *
        (exp.diags.arr[,, 4, i] /
          (exp.diags.arr[,, 1, i] +
            exp.diags.arr[,, 2, i] +
            exp.diags.arr[,, 3, i] +
            exp.diags.arr[,, 4, i]))
    }))
  }

  ### Remember that because of quarterly
  ### dynamics no CD4 dx can occur
  ### at times 5,9,13,... for age 1
  ### Thus replace all na by 0
  if (quar.flag) {
    CD4.st1[which(is.na(CD4.st1) == TRUE)] <- 0
    CD4.st2[which(is.na(CD4.st2) == TRUE)] <- 0
    CD4.st3[which(is.na(CD4.st3) == TRUE)] <- 0
    CD4.st4[which(is.na(CD4.st4) == TRUE)] <- 0
    if (config$rita) {
      CD4.st5[which(is.na(CD4.st5) == TRUE)] <- 0
    }
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
  if (urep.flag) {
    urep.vec <- fit.mat[, "under_rep"]
    AIDS <- simplify2array(lapply(1:nrow(fit.mat), function(i) {
      urep.fct(
        urep.vec[i],
        exp.diags.arr[,, 5, i],
        urep.flag = urep.flag,
        qt.flag = quar.flag,
        nt = nt
      )
    }))
  } else if (config$rita) {
    AIDS <- exp.diags.arr[,, 6, ]
  } else {
    AIDS <- exp.diags.arr[,, 5, ]
  }

  ### Storing summary statistics
  if (!config$rita) {
    CD4.st5 <- NULL
  }

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

  for (base_name in names(gof_types)) {
    type_arg <- gof_types[[base_name]]

    for (bn in names(age_bands)) {
      ages <- age_bands[[bn]]
      suf <- band_suffix(bn)

      out[[paste0(base_name, suf)]] <- gof.fct(
        type_arg,
        ages[1],
        ages[2],
        HIV,
        AIDS,
        CD4.st1,
        CD4.st2,
        CD4.st3,
        CD4.st4,
        CD4.st5,
        nt,
        quar.flag
      )
    }
  }

  ### Prevalence
  prev_states <- list(
    prev = "all",
    `prev.st1` = 1,
    `prev.st2` = 2,
    `prev.st3` = 3,
    `prev.st4` = 4
  )
  if (config$rita) {
    prev_states <- c(
      prev_states,
      list(`prev.st5` = 5, `prev.st6` = 6, `prev.st7` = 7)
    )
  }

  for (base_name in names(prev_states)) {
    state <- prev_states[[base_name]]

    for (bn in names(age_bands)) {
      ages <- age_bands[[bn]]
      suf <- band_suffix(bn)

      out[[paste0(base_name, suf)]] <- prev.fct(
        exp.prev.arr,
        state,
        ages[1],
        ages[2],
        config$rita
      )
    }
  }

  ### Getting data posterior predictive distributions
  N.Poiss <- 300

  HIV.m <- apply(HIV, c(1, 2), mean)

  HIV.p <- apply(HIV.m, c(1, 2), function(x) stats::rpois(N.Poiss, lambda = x))

  AIDS.m <- apply(AIDS, c(1, 2), mean)

  AIDS.p <- apply(AIDS.m, c(1, 2), function(x) {
    stats::rpois(N.Poiss, lambda = x)
  })

  if (config$rita) {
    p1 <- simplify2array(lapply(1:nrow(fit.mat), function(i) {
      (exp.diags.arr[,, 1, i] /
        (exp.diags.arr[,, 1, i] +
          exp.diags.arr[,, 2, i] +
          exp.diags.arr[,, 3, i] +
          exp.diags.arr[,, 4, i] +
          exp.diags.arr[,, 5, i]))
    }))
    p2 <- simplify2array(lapply(1:nrow(fit.mat), function(i) {
      (exp.diags.arr[,, 2, i] /
        (exp.diags.arr[,, 1, i] +
          exp.diags.arr[,, 2, i] +
          exp.diags.arr[,, 3, i] +
          exp.diags.arr[,, 4, i] +
          exp.diags.arr[,, 5, i]))
    }))
    p3 <- simplify2array(lapply(1:nrow(fit.mat), function(i) {
      (exp.diags.arr[,, 3, i] /
        (exp.diags.arr[,, 1, i] +
          exp.diags.arr[,, 2, i] +
          exp.diags.arr[,, 3, i] +
          exp.diags.arr[,, 4, i] +
          exp.diags.arr[,, 5, i]))
    }))
    p4 <- simplify2array(lapply(1:nrow(fit.mat), function(i) {
      (exp.diags.arr[,, 4, i] /
        (exp.diags.arr[,, 1, i] +
          exp.diags.arr[,, 2, i] +
          exp.diags.arr[,, 3, i] +
          exp.diags.arr[,, 4, i] +
          exp.diags.arr[,, 5, i]))
    }))
    p5 <- simplify2array(lapply(1:nrow(fit.mat), function(i) {
      (exp.diags.arr[,, 5, i] /
        (exp.diags.arr[,, 1, i] +
          exp.diags.arr[,, 2, i] +
          exp.diags.arr[,, 3, i] +
          exp.diags.arr[,, 4, i] +
          exp.diags.arr[,, 5, i]))
    }))
  } else {
    p1 <- simplify2array(lapply(1:nrow(fit.mat), function(i) {
      (exp.diags.arr[,, 1, i] /
        (exp.diags.arr[,, 1, i] +
          exp.diags.arr[,, 2, i] +
          exp.diags.arr[,, 3, i] +
          exp.diags.arr[,, 4, i]))
    }))
    p2 <- simplify2array(lapply(1:nrow(fit.mat), function(i) {
      (exp.diags.arr[,, 2, i] /
        (exp.diags.arr[,, 1, i] +
          exp.diags.arr[,, 2, i] +
          exp.diags.arr[,, 3, i] +
          exp.diags.arr[,, 4, i]))
    }))
    p3 <- simplify2array(lapply(1:nrow(fit.mat), function(i) {
      (exp.diags.arr[,, 3, i] /
        (exp.diags.arr[,, 1, i] +
          exp.diags.arr[,, 2, i] +
          exp.diags.arr[,, 3, i] +
          exp.diags.arr[,, 4, i]))
    }))
    p4 <- simplify2array(lapply(1:nrow(fit.mat), function(i) {
      (exp.diags.arr[,, 4, i] /
        (exp.diags.arr[,, 1, i] +
          exp.diags.arr[,, 2, i] +
          exp.diags.arr[,, 3, i] +
          exp.diags.arr[,, 4, i]))
    }))
  }

  if (quar.flag) {
    p1[which(is.na(p1) == TRUE)] <- 0
    p2[which(is.na(p2) == TRUE)] <- 0
    p3[which(is.na(p3) == TRUE)] <- 0
    p4[which(is.na(p4) == TRUE)] <- 0
    if (config$rita) {
      p5[which(is.na(p5) == TRUE)] <- 0
    }
  }

  p1.m <- apply(p1, c(1, 2), mean)
  p2.m <- apply(p2, c(1, 2), mean)
  p3.m <- apply(p3, c(1, 2), mean)
  p4.m <- apply(p4, c(1, 2), mean)
  if (config$rita) {
    p5.m <- apply(p5, c(1, 2), mean)
  }

  if (config$rita) {
    CD4.Poiss <- array(0, dim = c(nt, nage, 5, N.Poiss))
  } else {
    CD4.Poiss <- array(0, dim = c(nt, nage, 4, N.Poiss))
  }

  if (quar.flag && config$rita) {
    ### if quarterly there cannot be CD4 for a=1, at t=5, t=9,...
    ### as anyway nCD4[,1] = 0 then ignore and put 0
    for (t in 1:nt) {
      for (a in 2:nage) {
        CD4.Poiss[t, a, , ] <- stats::rmultinom(
          N.Poiss,
          size = nCD4[t, a],
          prob = c(p1.m[t, a], p2.m[t, a], p3.m[t, a], p4.m[t, a], p5.m[t, a])
        )
      }
    }
  } else if (quar.flag) {
    for (t in 1:nt) {
      for (a in 2:nage) {
        CD4.Poiss[t, a, , ] <- stats::rmultinom(
          N.Poiss,
          size = nCD4[t, a],
          prob = c(p1.m[t, a], p2.m[t, a], p3.m[t, a], p4.m[t, a])
        )
      }
    }
  } else {
    for (t in 1:nt) {
      for (a in 1:nage) {
        CD4.Poiss[t, a, , ] <- stats::rmultinom(
          N.Poiss,
          size = nCD4[t, a],
          prob = c(p1.m[t, a], p2.m[t, a], p3.m[t, a], p4.m[t, a])
        )
      }
    }
  }

  age_idx <- lapply(age_bands, function(ages) {
    .age_to_idx(ages[1], ages[2], nage = nage)
  })

  # If quarterly model: "non-Quar" outputs are yearly aggregates; else identity
  time_agg_year <- make_time_agg(quar.flag, nt)
  time_agg_quar <- function(ts) ts

  # --- HIV / AIDS Poisson summaries (yearly-or-quarterly aggregated) ---
  for (bn in names(age_idx)) {
    suf <- band_suffix(bn)
    out[[paste0("HIV", suf, ".Poiss")]] <- summarise_poiss(
      HIV.p,
      age_idx[[bn]],
      time_agg_year
    )
    out[[paste0("AIDS", suf, ".Poiss")]] <- summarise_poiss(
      AIDS.p,
      age_idx[[bn]],
      time_agg_year
    )
  }

  # --- CD4 multinomial summaries (yearly-or-quarterly aggregated) ---
  nstages <- if (config$rita) 5 else 4
  npoiss <- dim(CD4.Poiss)[4]

  for (s in seq_len(nstages)) {
    # CD4.Poiss dims: (nt, nage, stage, N.Poiss)
    stage_arr <- array(CD4.Poiss[,, s, ], dim = c(nt, nage, npoiss))

    for (bn in names(age_idx)) {
      suf <- band_suffix(bn)
      out[[paste0("CD4.", s, suf, ".Mult")]] <- summarise_mult(
        stage_arr,
        age_idx[[bn]],
        time_agg_year
      )
    }
  }

  # --- Quarterly (non-aggregated) outputs (only when quar.flag == TRUE) ---
  if (quar.flag) {
    for (bn in names(age_idx)) {
      suf <- band_suffix(bn)
      out[[paste0("HIV", suf, ".Quar.Poiss")]] <- summarise_poiss(
        HIV.p,
        age_idx[[bn]],
        time_agg_quar
      )
      out[[paste0("AIDS", suf, ".Quar.Poiss")]] <- summarise_poiss(
        AIDS.p,
        age_idx[[bn]],
        time_agg_quar
      )
    }

    for (s in seq_len(nstages)) {
      stage_arr <- array(CD4.Poiss[,, s, ], dim = c(nt, nage, npoiss))

      for (bn in names(age_idx)) {
        suf <- band_suffix(bn)
        out[[paste0("CD4.", s, suf, ".Quar.Mult")]] <- summarise_mult(
          stage_arr,
          age_idx[[bn]],
          time_agg_quar
        )
      }
    }
  }

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

make_time_agg <- function(quar_flag, nt_local) {
  if (isTRUE(quar_flag)) {
    if (nt_local %% 4 != 0) {
      stop("nt must be divisible by 4 when quar.flag=TRUE")
    }
    grp <- rep(seq_len(nt_local / 4), each = 4)
    function(ts) as.numeric(rowsum(ts, grp, reorder = FALSE))
  } else {
    function(ts) ts
  }
}

summarise_poiss <- function(p_arr, age_cols, time_agg) {
  # p_arr dims from apply(rpois): (N.Poiss, nt, nage)
  # for each Poiss draw: matrix(nt, nage) -> rowSums over age_cols -> time_agg
  tmp <- apply(p_arr[,, age_cols, drop = FALSE], 1, function(m) {
    time_agg(rowSums(m))
  })
  apply(tmp, 1, get_quantiles)
}

summarise_mult <- function(stage_arr, age_cols, time_agg) {
  # stage_arr dims: (nt, nage, N.Poiss)
  tmp <- apply(stage_arr[, age_cols, , drop = FALSE], 3, function(m) {
    time_agg(rowSums(m))
  })
  apply(tmp, 1, get_quantiles)
}
