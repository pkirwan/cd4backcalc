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
#' @param iter_warmup Number of warmup iterations for Stan sampling.
#' @param iter_sampling Number of sampling iterations for Stan sampling.
#' @param adapt_delta Adapt delta parameter for Stan sampling (default is 0.95).
#' @param max_treedepth Max tree depth for Stan sampling (default is 12).
#' @param chains Number of chains for Stan sampling.
#' @param refresh Number of iterations between progress messages.
#' @param messages Logical indicating whether to show messages during sampling.
#' @param start_yr Starting year for posterior predictive quantities.
#' @param nworkers Number of parallel workers to use.
#' @param chain_results Logical indicating whether to return results by chain.
#' @param return_fit Logical indicating whether to return the fit object
#' @param checkpoint_dir Directory for saving/loading checkpoints.
#'   If \code{NULL} (default), no checkpointing is performed.
#' @param iter_per_chunk Number of sampling iterations per chunk when
#'   checkpointing.  If \code{NULL}, all sampling iterations are run in a
#'   single chunk.
#' @return List of posterior predictive quantities and (optionally) the stanfit object.
#' @examples
#' sim_diags <- simulate_diagnoses(sim_seed = 123)
#' hiv_list <- run_backcalc(
#'  sim_diags, iter_warmup = 10, iter_sampling = 10, chains = 1, messages = FALSE, nworkers = 1
#' )
#' @export
run_backcalc <- function(
  stan_data,
  migration = NULL,
  rita = NULL,
  age = NULL,
  inf_model = 1,
  diag_model = 1,
  model_seed = NULL,
  iter_warmup = 1000,
  iter_sampling = 1000,
  chains = 4,
  adapt_delta = 0.95,
  max_treedepth = 12,
  refresh = 100,
  messages = TRUE,
  start_yr = 0,
  nworkers = 1,
  chain_results = FALSE,
  return_fit = TRUE,
  checkpoint_dir = NULL,
  iter_per_chunk = NULL
) {
  # Check if sim_diags object has been provided
  if ("stan_data" %in% names(stan_data)) {
    sim_data <- stan_data$sim_data
    stan_data <- stan_data$stan_data
  } else {
    sim_data <- NULL
  }

  # Prepare spline basis and augment stan_data with model_config
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
    iter_warmup = iter_warmup,
    iter_sampling = iter_sampling,
    adapt_delta = adapt_delta,
    max_treedepth = max_treedepth,
    chains = chains,
    refresh = refresh,
    messages = messages,
    checkpoint_dir = checkpoint_dir,
    iter_per_chunk = iter_per_chunk
  )

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
#' When \code{checkpoint_dir} is provided, fitting is split into a warmup
#' phase followed by one or more sampling chunks of size
#' \code{iter_per_chunk}.
#' After each chunk a checkpoint is written so that a timed-out SLURM job
#' can be resubmitted and will resume from the last completed chunk.
#'
#' @param stan_data A list of data formatted for Stan.
#' @param model_seed Seed for model fitting (will be set to stan_data$model_config$model_seed if NULL).
#' @param iter_warmup Number of warmup iterations for Stan sampling.
#' @param iter_sampling Number of sampling iterations for Stan sampling.
#' @param chains Number of chains for Stan sampling.
#' @param adapt_delta Adapt delta parameter for Stan sampling (default is 0.95).
#' @param max_treedepth Max tree depth for Stan sampling (default is 12).
#' @param refresh Number of iterations between progress messages.
#' @param messages Logical indicating whether to show messages during sampling.
#' @param checkpoint_dir Directory for saving/loading checkpoints.
#'   If \code{NULL} (default), no checkpointing is performed.
#' @param iter_per_chunk Number of sampling iterations per chunk when
#'   checkpointing.  If \code{NULL}, all sampling iterations are run in a
#'   single chunk.
#' @return A CmdStanFit object (without checkpointing) or a
#'   \code{cd4backcalc_combined_fit} object (with checkpointing).
#' @examples
#' \dontrun{
#' sim_diags <- simulate_diagnoses(sim_seed = 123)
#' stan_data <- model_type_prep(sim_diags)
#' fit <- fit_stan_model(
#'  stan_data, iter_warmup = 10, iter_sampling = 10, chains = 1, refresh = 0, messages = FALSE
#' )
#' }
#' @keywords internal
fit_stan_model <- function(
  stan_data,
  model_seed = NULL,
  iter_warmup = 1000,
  iter_sampling = 1000,
  chains = 4,
  adapt_delta = 0.95,
  max_treedepth = 12,
  refresh = 100,
  messages = TRUE,
  checkpoint_dir = NULL,
  iter_per_chunk = NULL
) {
  options(mc.cores = parallelly::availableCores())

  # model seed
  if (is.null(model_seed)) {
    model_seed <- stan_data$model_config$model_seed
  }

  # Check stan_data validity
  confirm_stan_data(stan_data)

  # Load appropriate model
  model <- load_stan_model(stan_data$model_config)

  model_data <- stan_data[
    !names(stan_data) %in% c("sim_config", "model_config")
  ]

  # ---- Non-checkpoint path (original behaviour) ----
  if (is.null(checkpoint_dir)) {
    fit <- model$sample(
      data = model_data,
      iter_warmup = iter_warmup,
      iter_sampling = iter_sampling,
      chains = chains,
      seed = model_seed,
      adapt_delta = adapt_delta,
      max_treedepth = max_treedepth,
      refresh = refresh,
      show_exceptions = messages,
      show_messages = messages
    )
    return(fit)
  }

  # ---- Checkpoint path ----
  if (is.null(iter_per_chunk)) {
    iter_per_chunk <- iter_sampling
  }

  cp <- load_checkpoint(checkpoint_dir)
  if (!is.null(cp)) cp <- pad_checkpoint_chains(cp, chains)
  draws_list <- if (!is.null(cp)) cp$draws_list else list()
  diagnostics_list <- if (!is.null(cp)) cp$diagnostics_list else list()
  total_time <- if (!is.null(cp)) cp$total_time else 0

  # Phase 1: Warmup in chunks (only if no checkpoint or still in warmup)
  if (is.null(cp) || cp$phase == "warmup") {
    warmup_done <- if (!is.null(cp)) cp$warmup_done else 0L

    while (warmup_done < iter_warmup) {
      this_chunk <- min(iter_per_chunk, iter_warmup - warmup_done)
      message(
        "Checkpoint: warmup chunk ",
        warmup_done + 1, "-", warmup_done + this_chunk,
        " of ", iter_warmup, "..."
      )

      # First warmup chunk: fresh start; subsequent chunks: seed from
      # previous adaptation state but keep adaptation enabled so the
      # sampler continues tuning.
      # For resumed chunks we skip the init_buffer (step size already
      # tuned) and, for non-final chunks, also skip the term_buffer so
      # that all iterations go to mass-matrix adaptation windows.
      is_last_warmup <- (warmup_done + this_chunk >= iter_warmup)

      if (warmup_done == 0) {
        fit <- model$sample(
          data = model_data,
          iter_warmup = this_chunk,
          iter_sampling = 0,
          chains = chains,
          seed = model_seed,
          adapt_delta = adapt_delta,
          max_treedepth = max_treedepth,
          refresh = refresh,
          show_exceptions = messages,
          show_messages = messages,
          save_warmup = TRUE,
          sig_figs = 18
        )
      } else {
        fit <- model$sample(
          data = model_data,
          iter_warmup = this_chunk,
          iter_sampling = 0,
          chains = chains,
          seed = model_seed,
          adapt_delta = adapt_delta,
          max_treedepth = max_treedepth,
          init_buffer = 0,
          term_buffer = if (is_last_warmup) NULL else 0,
          inv_metric = cp$inv_metric,
          step_size = cp$step_size,
          init = cp$init,
          refresh = refresh,
          show_exceptions = messages,
          show_messages = messages,
          save_warmup = TRUE,
          sig_figs = 18
        )
      }

      warmup_done <- warmup_done + this_chunk
      total_time <- total_time + fit$time()$total

      if (warmup_done < iter_warmup) {
        # More warmup to go — stay in warmup phase
        save_checkpoint(
          checkpoint_dir, fit, "warmup", warmup_done, 0L,
          draws_list, diagnostics_list, total_time
        )
      } else {
        # Warmup complete — transition to sampling phase
        save_checkpoint(
          checkpoint_dir, fit, "sampling", warmup_done, 0L,
          draws_list, diagnostics_list, total_time
        )
      }
      cp <- load_checkpoint(checkpoint_dir)
      cp <- pad_checkpoint_chains(cp, chains)
      message("Checkpoint: ", warmup_done, "/", iter_warmup, " warmup iterations complete.")
    }
  }

  # Phase 2: Sampling in chunks
  sampling_done <- cp$sampling_done

  while (sampling_done < iter_sampling) {
    this_chunk <- min(iter_per_chunk, iter_sampling - sampling_done)
    message(
      "Checkpoint: sampling chunk ",
      sampling_done + 1, "-", sampling_done + this_chunk,
      " of ", iter_sampling, "..."
    )

    fit <- model$sample(
      data = model_data,
      iter_warmup = 0,
      iter_sampling = this_chunk,
      chains = chains,
      seed = model_seed,
      adapt_engaged = FALSE,
      inv_metric = cp$inv_metric,
      step_size = cp$step_size,
      init = cp$init,
      refresh = refresh,
      show_exceptions = messages,
      show_messages = messages,
      sig_figs = 18
    )

    sampling_done <- sampling_done + this_chunk
    draws_list <- append(draws_list, list(fit$draws()))
    diagnostics_list <- append(
      diagnostics_list, list(fit$diagnostic_summary())
    )
    total_time <- total_time + fit$time()$total

    save_checkpoint(
      checkpoint_dir, fit, "sampling", warmup_done = iter_warmup,
      sampling_done, draws_list, diagnostics_list, total_time
    )
    cp <- load_checkpoint(checkpoint_dir)
    cp <- pad_checkpoint_chains(cp, chains)
    message("Checkpoint: ", sampling_done, "/", iter_sampling, " sampling iterations complete.")
  }

  # Combine draws from all chunks into a single fit-like object
  if (length(draws_list) == 1) {
    return(fit)
  }
  create_combined_fit(fit, draws_list, diagnostics_list, chains, total_time)
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
    model_seed = model_seed,
    ndiagst = ifelse(rita, 6, 5),
    nlatst = ifelse(rita, 7, 4)
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
  if (stan_data$model_config$age && stan_data$model_config$migration) {
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
    if (stan_data$model_config$rita) {
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
  } else if (stan_data$model_config$age) {
    stopifnot(
      "Incorrect dimensions in stan_data, HIV_age should be (nquar, nage)" = all(
        dim(stan_data$HIV_age) == c(stan_data$nquar, stan_data$nage)
      )
    )
    if (stan_data$model_config$rita) {
      stopifnot(
        "Incorrect dimensions in stan_data, CD4_age should be (nquar, nage, 5)" = all(
          dim(stan_data$CD4_age) == c(stan_data$nquar, stan_data$nage, 5)
        )
      )
    }
  } else if (stan_data$model_config$migration) {
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
    if (stan_data$model_config$rita) {
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
  } else if (stan_data$model_config$rita) {
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

#' Determine the CmdStan executable cache directory
#'
#' This helper picks a writable location for compiled CmdStan executables.
#' It prefers a job-local temporary directory to avoid failures from shared
#' home-directory caches on batch systems, and allows an explicit override via
#' the `CD4BACKCALC_CMDSTAN_DIR` environment variable.
#'
#' @return The absolute path to a writable directory for CmdStan executables.
#' @noRd
stan_model_cache_dir <- function() {
  override_dir <- Sys.getenv("CD4BACKCALC_CMDSTAN_DIR", unset = "")
  if (nzchar(override_dir)) {
    dir.create(override_dir, recursive = TRUE, showWarnings = FALSE)
    return(normalizePath(override_dir, winslash = "/", mustWork = FALSE))
  }

  temp_cache_dir <- file.path(tempdir(), "cd4backcalc-cmdstan")
  dir.create(temp_cache_dir, recursive = TRUE, showWarnings = FALSE)

  if (dir.exists(temp_cache_dir) && file.access(temp_cache_dir, 2) == 0) {
    return(normalizePath(temp_cache_dir, winslash = "/", mustWork = FALSE))
  }

  user_cache_dir <- file.path(
    tools::R_user_dir("cd4backcalc", which = "cache"),
    "cmdstan"
  )
  dir.create(user_cache_dir, recursive = TRUE, showWarnings = FALSE)
  normalizePath(user_cache_dir, winslash = "/", mustWork = FALSE)
}

#' Choose Stan model based on configuration
#'
#' This function chooses the appropriate Stan model based on the provided
#' configuration. The model is compiled at runtime using cmdstanr, with the
#' compiled binary cached so subsequent calls are fast.
#'
#' @param config Configuration object containing model settings.
#' @return The Stan model object (a CmdStanModel).
load_stan_model <- function(config) {
  name <- if (config$age) {
    if (config$migration && config$rita) {
      "AD_ptens_migration_rita_softmax"
    } else if (config$migration) {
      "AD_ptens_migration_softmax"
    } else if (config$rita) {
      "AD_ptens_rita"
    } else {
      "AD_ptens"
    }
  } else if (config$migration && config$rita) {
    "AI_spl_migration_rita"
  } else if (config$migration) {
    "AI_spl_migration"
  } else if (config$rita) {
    "AI_spl_rita"
  } else if (
    config$diag_model == 2 || config$inf_model == 2 || config$inf_model == 3
  ) {
    "AI_spl_rw_gp"
  } else {
    "AI_spl"
  }

  # Find the .stan file shipped with the package
  stan_file <- system.file(
    "stan",
    paste0(name, ".stan"),
    package = "cd4backcalc"
  )
  if (!nzchar(stan_file) || !file.exists(stan_file)) {
    stop(
      "Stan model file '",
      name,
      ".stan' not found in the installed package. ",
      "The package may need to be reinstalled.",
      call. = FALSE
    )
  }

  # Prefer a job-local cache to avoid failures from shared home-directory caches
  # on batch systems, while still allowing an explicit override when needed.
  cmdstanr::cmdstan_model(
    stan_file = stan_file,
    dir = stan_model_cache_dir()
  )
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
#' @keywords internal
model_type_prep <- function(
  stan_data,
  migration = NULL,
  rita = NULL,
  age = NULL,
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

  # Set migration, rita, age based on stan_data if NULL
  # If NULL and sim_config values are available, use those; otherwise default to FALSE
  if (is.null(migration)) {
    migration <- ifelse(
      !is.null(stan_data$sim_config$migration),
      stan_data$sim_config$migration,
      FALSE
    )
  }
  if (is.null(rita)) {
    rita <- ifelse(
      !is.null(stan_data$sim_config$rita),
      stan_data$sim_config$rita,
      FALSE
    )
  }
  if (is.null(age)) {
    age <- ifelse(
      !is.null(stan_data$sim_config$age),
      stan_data$sim_config$age,
      FALSE
    )
  }

  # Create stan_data$model_config
  stan_data$model_config <- stan_model_config(
    migration,
    rita,
    age,
    inf_model,
    diag_model,
    model_seed
  )

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
  jagam_out_diag <- create_spline(knots - 2)
  stan_data$X_diag <- jagam_out_diag$jags.data$X
  stan_data$npars_diag <- ncol(stan_data$X_diag)

  if (stan_data$model_config$age) {
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

  if (stan_data$model_config$migration) {
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
    is.null(stan_data$model_config$inf_model),
    1,
    stan_data$model_config$inf_model
  )
  stan_data$diag_model <- ifelse(
    is.null(stan_data$model_config$diag_model),
    1,
    stan_data$model_config$diag_model
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
#' @keywords internal
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
