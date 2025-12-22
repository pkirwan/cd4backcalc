#' Run CD4 back-calculation model
#'
#' This function runs the CD4 back-calculation model using Stan and runs the post-processing.
#'
#' @param stan_data A list of data formatted for Stan.
#' @param migration Logical indicating whether to include migration in the model.
#' @param rita Logical indicating whether to include RITA data in the model.
#' @param age Logical indicating whether to include age-specific data in the model.
#' @param inf_model Infection model type (1 = spline, 2 = random walk, 3 = GP).
#' @param diag_model Diagnosis model type (1 = spline, 2 = random walk).
#' @param model_seed Seed for model fitting (will be set randomly if NULL).
#' @param iter Number of iterations for Stan sampling.
#' @param chains Number of chains for Stan sampling.
#' @param refresh Number of iterations between progress messages.
#' @param messages Logical indicating whether to show messages during sampling.
#' @param start_yr Starting year for posterior predictive quantities.
#' @param nworkers Number of parallel workers to use.
#' @param chain_results Logical indicating whether to return results by chain.
#' @param return_fit Logical indicating whether to return the fit object
#' @return List of posterior predictive quantities and (optionlly) the stanfit object.
#' @examples
#' sim_diags <- simulate_diagnoses(sim_seed = 123)
#' hiv_list <- run_backcalc(sim_diags, iter = 5, chains = 1, messages = FALSE, nworkers = 1)
#' @export
run_backcalc <- function(
  stan_data,
  migration = FALSE,
  rita = FALSE,
  age = FALSE,
  inf_model = 1,
  diag_model = 1,
  model_seed = NULL,
  iter = 2000,
  chains = 4,
  refresh = 100,
  messages = TRUE,
  start_yr = 0,
  nworkers = 1,
  chain_results = FALSE,
  return_fit = TRUE
) {
  # Check if sim_diags object has been provided
  if ("stan_data" %in% names(stan_data)) {
    sim_data <- stan_data$sim_data
    stan_data <- stan_data$stan_data
  } else {
    sim_data <- NULL
  }

  # Prepare spline basis and augment stan_data
  stan_data <- model_type_prep(
    stan_data = stan_data,
    migration = migration,
    rita = rita,
    age = age,
    inf_model = inf_model,
    diag_model = diag_model,
    model_seed = model_seed
  )

  # fit model
  fit <- fit_stan_model(
    stan_data = stan_data,
    model_seed = model_seed,
    iter = iter,
    chains = chains,
    refresh = refresh,
    messages = messages
  )

  gc() # clean up memory before post-processing

  # post-process
  hiv_list <- postproc(fit, stan_data, start_yr, nworkers, chain_results)

  if (return_fit) {
    hiv_list$fit <- fit
  }
  hiv_list$stan_data <- stan_data
  hiv_list$sim_data <- sim_data

  return(hiv_list)
}

#' Run CD4 back-calculation model
#'
#' This function runs the CD4 back-calculation model using Stan.
#'
#' @param stan_data A list of data formatted for Stan.
#' @param model_seed Seed for model fitting (will be set to stan_data$config$model_seed if NULL).
#' @param iter Number of iterations for Stan sampling.
#' @param chains Number of chains for Stan sampling.
#' @param refresh Number of iterations between progress messages.
#' @param messages Logical indicating whether to show messages during sampling.
#' @return A stanfit object containing the results of the model fitting.
#' @examples
#' sim_diags <- simulate_diagnoses(sim_seed = 123)
#' stan_data <- model_type_prep(sim_diags)
#' fit <- fit_stan_model(stan_data, iter = 10, chains = 1, refresh = 0, messages = FALSE)
#' @export
fit_stan_model <- function(
  stan_data,
  model_seed = NULL,
  iter = 2000,
  chains = 4,
  refresh = 100,
  messages = TRUE
) {
  # rstan options
  options(mc.cores = parallelly::availableCores())
  rstan::rstan_options(auto_write = TRUE)

  # model seed
  if (is.null(model_seed)) {
    model_seed <- stan_data$config$model_seed
  }

  # Check stan_data validity
  confirm_stan_data(stan_data)

  # Load appropriate model
  model <- load_stan_model(stan_data$config)

  # Make pars_save
  pars_save <- make_pars_save(stan_data$config)

  # Sample from model
  fit <- rstan::sampling(
    model,
    data = stan_data,
    iter = iter,
    chains = chains,
    seed = model_seed,
    pars = pars_save,
    control = list(
      adapt_delta = 0.95,
      max_treedepth = 12
    ),
    save_warmup = FALSE,
    refresh = refresh,
    verbose = messages,
    show_messages = messages
  )

  return(fit)
}

#' Define configuration for fitting back-calculation model
#'
#' This function defines the configuration settings for fitting the CD4 back-calculation model.
#'
#' @param migration Logical indicating whether to include migration in the model.
#' @param rita Logical indicating whether to include RITA data in the model.
#' @param age Logical indicating whether to include age-specific data in the model.
#' @param inf_model Infection model type (1 = spline, 2 = random walk, 3 = GP).
#' @param diag_model Diagnosis model type (1 = spline, 2 = random walk).
#' @param model_seed Seed for model fitting.
#' @return A list containing configuration settings for model fitting.
stan_model_config <- function(
  migration = FALSE,
  rita = FALSE,
  age = FALSE,
  inf_model = 1,
  diag_model = 1,
  model_seed = NULL
) {
  if ((inf_model != 1 || diag_model != 1) && (migration || rita || age)) {
    stop(
      "The migration, rita, or age-specific models only support spline model configurations. Please use inf_model = 1 and diag_model = 1."
    )
  }
  if (is.null(model_seed)) {
    model_seed <- sample.int(.Machine$integer.max, 1)
  }
  config <- list(
    inf_model = inf_model,
    diag_model = diag_model,
    migration = migration,
    rita = rita,
    age = age,
    model_seed = model_seed
  )
  return(config)
}

#' Confirm stan_data validity
#'
#' This function checks that the stan_data list contains the expected elements
#' based on the provided configuration.
#'
#' @param stan_data A list of data formatted for Stan.
#' @return NULL. Throws an error if checks fail.
confirm_stan_data <- function(stan_data) {
  if (stan_data$config$age && stan_data$config$migration) {
    stopifnot(
      "Incorrect dimensions in stan_data, HIV_age should be (nquar, nage, 2)" = all(
        dim(stan_data$HIV_age) == c(stan_data$nquar, stan_data$nage, 2)
      )
    )
    stopifnot(
      "Incorrect dimensions in stan_data, AIDS_age should be (nquar, nage, 2)" = all(
        dim(stan_data$AIDS_age) == c(stan_data$nquar, stan_data$nage, 2)
      )
    )
    if (stan_data$config$rita) {
      stopifnot(
        "Incorrect dimensions in stan_data, CD4_age should be (nquar, nage, 2, 5)" = all(
          dim(stan_data$CD4_age) == c(stan_data$nquar, stan_data$nage, 2, 5)
        )
      )
    } else {
      stopifnot(
        "Incorrect dimensions in stan_data, CD4_age should be (nquar, nage, 2, 4)" = all(
          dim(stan_data$CD4_age) == c(stan_data$nquar, stan_data$nage, 2, 4)
        )
      )
    }
    stopifnot(
      "Incorrect dimensions in stan_data, ABROAD_age should be (nquar, nage)" = all(
        dim(stan_data$ABROAD_age) == c(stan_data$nquar, stan_data$nage)
      )
    )
  } else if (stan_data$config$age) {
    stopifnot(
      "Incorrect dimensions in stan_data, HIV_age should be (nquar, nage)" = all(
        dim(stan_data$HIV_age) == c(stan_data$nquar, stan_data$nage)
      )
    )
    if (stan_data$config$rita) {
      stopifnot(
        "Incorrect dimensions in stan_data, CD4_age should be (nquar, nage, 5)" = all(
          dim(stan_data$CD4_age) == c(stan_data$nquar, stan_data$nage, 5)
        )
      )
    }
  } else if (stan_data$config$migration) {
    stopifnot(
      "Incorrect dimensions in stan_data, HIV should be (nquar, 2)" = all(
        dim(stan_data$HIV) == c(stan_data$nquar, 2)
      )
    )
    stopifnot(
      "Incorrect dimensions in stan_data, AIDS should be (nquar, 2)" = all(
        dim(stan_data$AIDS) == c(stan_data$nquar, 2)
      )
    )
    stopifnot(
      "Incorrect length in stan_data, ABROAD should be (nquar)" = length(
        stan_data$ABROAD
      ) ==
        stan_data$nquar
    )
    if (stan_data$config$rita) {
      stopifnot(
        "Incorrect dimensions in stan_data, CD4 should be (nquar, 2, 5)" = all(
          dim(stan_data$CD4) == c(stan_data$nquar, 2, 5)
        )
      )
    } else {
      stopifnot(
        "Incorrect dimensions in stan_data, CD4 should be (nquar, 2, 4)" = all(
          dim(stan_data$CD4) == c(stan_data$nquar, 2, 4)
        )
      )
    }
  } else if (stan_data$config$rita) {
    stopifnot(
      "Incorrect dimensions in stan_data, CD4 should be (nquar, 5)" = all(
        dim(stan_data$CD4) == c(stan_data$nquar, 5)
      )
    )
  } else {
    stopifnot(
      "Incorrect dimensions in stan_data, CD4 should be (nquar, 4)" = all(
        dim(stan_data$CD4) == c(stan_data$nquar, 4)
      )
    )
  }
  return(NULL)
}

#' Choose Stan model based on configuration
#'
#' This function chooses the appropriate Stan model based on the provided configuration.
#'
#' @param config Configuration object containing model settings.
#' @return The Stan model object.
load_stan_model <- function(config) {
  if (config$age) {
    if (config$migration && config$rita) {
      model <- stanmodels$AD_ptens_migration_rita
    } else if (config$migration) {
      model <- stanmodels$AD_ptens_migration
    } else if (config$rita) {
      model <- stanmodels$AD_ptens_rita
    } else {
      model <- stanmodels$AD_ptens
    }
  } else if (config$migration && config$rita) {
    model <- stanmodels$AI_spl_migration_rita
  } else if (config$migration) {
    model <- stanmodels$AI_spl_migration
  } else if (config$rita) {
    model <- stanmodels$AI_spl_rita
  } else {
    model <- stanmodels$AI_spl_rw_gp
  }
  return(model)
}

#' Prepare spline basis for back-calculation model
#'
#' This function prepares the spline basis and augments the stan_data list
#' with the design matrices for the back-calculation model.
#'
#' @param stan_data A list of data formatted for Stan.
#' @param migration Logical indicating whether to include migration in the model.
#' @param rita Logical indicating whether to include RITA data in the model.
#' @param age Logical indicating whether to include age-specific data in the model.
#' @param inf_model Infection model type (1 = spline, 2 = random walk, 3 = GP).
#' @param diag_model Diagnosis model type (1 = spline, 2 = random walk).
#' @param model_seed Seed for model fitting.
#' @param basis Basis type for the spline (default is "ts" for thin plate regression spline).
#' @param knots Number of knots for the spline (default is 10).
#' @param penalty Order of the penalty for the spline (default is 2).
#' @return The augmented stan_data list with spline design matrices and prior information.
#' @export
model_type_prep <- function(
  stan_data,
  migration = FALSE,
  rita = FALSE,
  age = FALSE,
  inf_model = 1,
  diag_model = 1,
  model_seed = NULL,
  basis = "ts",
  knots = 10,
  penalty = 2
) {
  # Check if sim_diags object has been provided
  if ("stan_data" %in% names(stan_data)) {
    sim_diags <- stan_data$sim_diags
    stan_data <- stan_data$stan_data
  }

  # Add model_seed if NULL
  if (is.null(model_seed)) {
    model_seed <- sample.int(.Machine$integer.max, 1)
  }

  # Create stan_data$config if it doesn't exist
  if (is.null(stan_data$config)) {
    stan_data$config <- stan_model_config(
      migration,
      rita,
      age,
      inf_model,
      diag_model,
      model_seed
    )
  } else {
    stan_data$config$inf_model <- inf_model
    stan_data$config$diag_model <- diag_model
    stan_data$config$model_seed <- model_seed
  }

  # Generate fake data for jagam to create the spline basis
  jags_data <- list(
    qts = 1:stan_data$nquar,
    D = stats::runif(stan_data$nquar)
  )

  # Creates the spline objects using a bug file with some stan code
  # diagonalize = TRUE makes a reparameterization so that the priors are iid normal
  create_spline <- function(k) {
    mgcv::jagam(
      D ~ s(qts, bs = basis, k = k, m = penalty),
      family = stats::gaussian,
      data = jags_data,
      file = tempfile(),
      diagonalize = TRUE
    )
  }

  # Diagnosis model spline
  jagam_out_diag <- create_spline(knots - 4)
  stan_data$X_diag <- jagam_out_diag$jags.data$X
  stan_data$npars_diag <- ncol(stan_data$X_diag)

  if (stan_data$config$age) {
    # Age-specific infection model spline
    jags_data_age <- list(
      age = rep(1:stan_data$nage, each = stan_data$nquar),
      qts = rep(1:stan_data$nquar, times = stan_data$nage),
      D = stats::runif(stan_data$nquar * stan_data$nage)
    )

    jagam_out_age <- mgcv::jagam(
      D ~ te(
        qts,
        age,
        bs = c("ps", "ps"),
        k = c(10, 8),
        m = list(c(2, 1), c(2, 1))
      ),
      family = stats::gaussian,
      data = jags_data_age,
      file = tempfile(),
      diagonalize = TRUE
    )

    # Augment stan_data with design matrix and prior info
    stan_data$X_inf <- jagam_out_age$jags.data$X
    stan_data$npars_inf <- ncol(stan_data$X_inf)
    stan_data$S1 <- jagam_out_age$jags.data$S1
  } else {
    # Infection model spline (age-independent)
    jagam_out <- create_spline(knots)
    stan_data$X_inf <- jagam_out$jags.data$X
    stan_data$npars_inf <- ncol(stan_data$X_inf)

    # For GP model
    stan_data$sigma_sq <- 0.001 ## fixed nugget for gp
    range01 <- function(x) {
      (x - min(x)) / (max(x) - min(x))
    } ## scales inputs of GP to 0-1.
    stan_data$x <- range01(1:stan_data$nquar)
  }

  if (stan_data$config$migration) {
    # Migration model splines
    jagam_out_mig <- create_spline(knots)
    jagam_out_cob <- create_spline(knots - 4)
    jagam_out_mp <- create_spline(knots - 2)

    # Augment stan_data with design matrix and prior info
    stan_data$X_mig <- jagam_out_mig$jags.data$X
    stan_data$npars_mig <- ncol(stan_data$X_mig)

    stan_data$X_cob <- jagam_out_cob$jags.data$X
    stan_data$npars_cob <- ncol(stan_data$X_cob)

    stan_data$X_mp <- jagam_out_mp$jags.data$X
    stan_data$npars_mp <- ncol(stan_data$X_mp)
  }

  # Add model types to stan_data
  stan_data$inf_model <- ifelse(
    is.null(stan_data$config$inf_model),
    1,
    stan_data$config$inf_model
  )
  stan_data$diag_model <- ifelse(
    is.null(stan_data$config$diag_model),
    1,
    stan_data$config$diag_model
  )

  return(stan_data)
}

#' Make parameters to save based on configuration
#'
#' This function generates a vector of parameter names to save from the Stan model
#' based on the provided configuration.
#'
#' @param config Configuration object containing model settings.
#' @return A character vector of parameter names to save.
make_pars_save <- function(config) {
  # basic parameters to save
  pars_save <- c(
    "lambda",
    "d",
    "log_lik",
    "exp_HIV_dx",
    "exp_AIDS_dx",
    "exp_p"
  )

  if (config$age) {
    # Add age-specific parameters
    pars_save <- c(
      pars_save,
      "tau",
      "beta_inf",
      "delta_1",
      "delta_2",
      "delta_3",
      "delta_4"
    )
    if (!config$migration) {
      pars_save <- c(
        pars_save,
        "vardelta",
        "alpha"
      )
    }
  } else {
    pars_save <- c(
      pars_save,
      "prev_mat"
    )
  }

  if (config$migration) {
    # Add migration-specific parameters
    pars_save <- c(
      pars_save,
      "gamma",
      "omega",
      "p",
      "m",
      "detect_prob",
      "log_lik_abroad",
      "misclassified_dx"
    )
  }

  return(pars_save)
}
