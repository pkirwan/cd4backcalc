#' Run CD4 back-calculation model
#'
#' This function runs the CD4 back-calculation model using Stan and runs the post-processing.
#'
#' @param stan_data A list of data formatted for Stan.
#' @param migration Logical indicating whether to include migration in the model.
#' @param reduced_migration Logical indicating whether to use the reduced
#'   age-dependent migration model.
#' @param rita Logical indicating whether to include RITA data in the model.
#' @param age Logical indicating whether to include age-specific data in the model.
#' @param inf_model Infection model type (1 = spline, 2 = random walk, 3 = GP).
#' @param diag_model Diagnosis model type (1 = spline, 2 = AR(1)).
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
#' @param init Initial values passed to Stan. When checkpointing, these are
#'   only used when starting a new run; resumed checkpointed fits use the
#'   saved checkpoint state.
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
  iter_per_chunk = NULL,
  init = NULL,
  reduced_migration = FALSE
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
    reduced_migration = reduced_migration,
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
    iter_per_chunk = iter_per_chunk,
    init = init
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
#' @param init Initial values passed to Stan. When checkpointing, these are
#'   only used when starting a new run; resumed checkpointed fits use the
#'   saved checkpoint state.
#' @return A CmdStanFit object (without checkpointing) or a
#'   \code{cd4backcalc_combined_fit} object (with checkpointing).
#' @examples
#' \donttest{
#' sim_diags <- simulate_diagnoses(sim_seed = 123)
#' stan_data <- cd4backcalc:::model_type_prep(sim_diags)
#' fit <- cd4backcalc:::fit_stan_model(
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
  iter_per_chunk = NULL,
  init = NULL
) {
  old_mc_cores <- getOption("mc.cores")
  options(mc.cores = parallelly::availableCores())
  on.exit(options(mc.cores = old_mc_cores), add = TRUE)

  # Helper: safely extract total time from a CmdStanFit even when some
  # chains failed (fit$time() can error with mismatched row counts).
  safe_fit_time <- function(fit) {
    tryCatch(fit$time()$total, error = function(e) NA_real_)
  }

  # model seed
  if (is.null(model_seed)) {
    model_seed <- stan_data$model_config$model_seed
  }

  # Check stan_data validity
  confirm_stan_data(stan_data)

  # Load appropriate model
  model <- load_stan_model(stan_data$model_config)

  model_data <- stan_data[
    !names(stan_data) %in% c("sim_config", "model_config", "ABROAD_age", "gamma_obs_grouped")
  ]

  # ---- Non-checkpoint path (original behaviour) ----
  if (is.null(checkpoint_dir)) {
    fit <- model$sample(
      data = model_data,
      iter_warmup = iter_warmup,
      iter_sampling = iter_sampling,
      chains = chains,
      seed = model_seed,
      init = init,
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

  message("Checkpointing enabled. Checkpoint directory: ", checkpoint_dir)
  cp <- load_checkpoint(checkpoint_dir)
  if (!is.null(cp)) cp <- pad_checkpoint_chains(cp, chains)
  if (!is.null(init) && !is.null(cp)) {
    warning(
      "Ignoring user-supplied init because an existing checkpoint was found; resuming from the saved checkpoint state.",
      call. = FALSE
    )
  }
  cmdstan_output_dir <- file.path(checkpoint_dir, "cmdstan-output")
  dir.create(cmdstan_output_dir, recursive = TRUE, showWarnings = FALSE)
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
          init = init,
          output_dir = cmdstan_output_dir,
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
          output_dir = cmdstan_output_dir,
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
      total_time <- total_time + safe_fit_time(fit)

      # Check for failed chains
      n_ok <- length(fit$metadata()$id)
      if (n_ok < chains) {
        warning(
          n_ok, " of ", chains, " chains completed warmup chunk. ",
          "Use fit$output() to inspect failures."
        )
      }

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
      output_dir = cmdstan_output_dir,
      adapt_engaged = FALSE,
      inv_metric = cp$inv_metric,
      step_size = cp$step_size,
      init = cp$init,
      refresh = refresh,
      show_exceptions = messages,
      show_messages = messages,
      sig_figs = 18
    )

    # Check for failed chains — error early rather than silently
    # continuing with mismatched draws.
    n_ok <- length(fit$metadata()$id)
    if (n_ok < chains) {
      stop(
        n_ok, " of ", chains, " chains completed the sampling chunk. ",
        "Use fit$output() to inspect failures. ",
        "The previous checkpoint is still on disk; resubmit to retry.",
        call. = FALSE
      )
    }

    sampling_done <- sampling_done + this_chunk

    # Always accumulate draws so the checkpoint contains all draws to date.
    # This ensures that a resumed run can combine draws from all chunks,
    # even when a prior run covered only a single chunk.
    draws_list <- append(draws_list, list(fit$draws()))
    diagnostics_list <- append(
      diagnostics_list, list(fit$diagnostic_summary())
    )
    total_time <- total_time + safe_fit_time(fit)

    save_checkpoint(
      checkpoint_dir, fit, "sampling", warmup_done = iter_warmup,
      sampling_done, draws_list, diagnostics_list, total_time
    )
    cp <- load_checkpoint(checkpoint_dir)
    cp <- pad_checkpoint_chains(cp, chains)
    message("Checkpoint: ", sampling_done, "/", iter_sampling, " sampling iterations complete.")
  }

  # Single chunk (no prior draws in checkpoint, only one chunk run):
  # return the raw CmdStanFit so callers get the native object.
  # Multiple chunks: combine into a single fit-like object.
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
#' @param diag_model Diagnosis model type (1 = spline, 2 = AR(1)).
#' @param model_seed Seed for model fitting.
#' @param reduced_migration Logical indicating whether to use the reduced age-dependent migration model.
#' @return A list containing configuration settings for model fitting.
stan_model_config <- function(
  migration = FALSE,
  rita = FALSE,
  age = FALSE,
  inf_model = 1,
  diag_model = 1,
  model_seed = NULL,
  reduced_migration = FALSE
) {
  if (reduced_migration && (!migration || !age || rita)) {
    stop("reduced_migration requires migration = TRUE, age = TRUE, and rita = FALSE.")
  }
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
    reduced_migration = reduced_migration,
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
    stopifnot(
      "Incorrect dimensions in stan_data, ABROAD_group should be (nquar, 4)" = all(
        dim(stan_data$ABROAD_group) == c(stan_data$nquar, 4L)
      )
    )
    stopifnot(
      "ABROAD_group row sums must equal ABROAD" = all(
        rowSums(stan_data$ABROAD_group) == stan_data$ABROAD
      )
    )
    stopifnot(
      "Incorrect length in stan_data, gamma_obs should be (nquar)" =
        length(stan_data$gamma_obs) == stan_data$nquar
    )
    stopifnot(
      "Incorrect dimensions in stan_data, gamma_obs_group should be (nquar, 4)" = all(
        dim(stan_data$gamma_obs_group) == c(stan_data$nquar, 4L)
      )
    )
    gamma_group_available <- rowSums(stan_data$gamma_obs_group >= 0)
    stopifnot(
      "Each row of gamma_obs_group must be either all -1 or have all 4 age-group values" =
        all(gamma_group_available %in% c(0L, 4L))
    )
    stopifnot(
      "n_detect_steps must match the length of detect_step_quarter" =
        identical(
          as.integer(stan_data$n_detect_steps),
          length(stan_data$detect_step_quarter)
        )
    )
    stopifnot(
      "gamma_obs_drift_sd must be a single non-negative number" =
        length(stan_data$gamma_obs_drift_sd) == 1L &&
          is.finite(stan_data$gamma_obs_drift_sd) &&
          stan_data$gamma_obs_drift_sd >= 0
    )
    stopifnot(
      "gamma_obs_slope_sd must be a single non-negative number" =
        length(stan_data$gamma_obs_slope_sd) == 1L &&
          is.finite(stan_data$gamma_obs_slope_sd) &&
          stan_data$gamma_obs_slope_sd >= 0
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
    stopifnot(
      "detect_logit_obs and detect_obs_index must have the same length" =
        length(stan_data$detect_logit_obs) == length(stan_data$detect_obs_index)
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
      "AD_ptens_migration_rita"
    } else if (config$migration && isTRUE(config$reduced_migration)) {
      "AD_ptens_migration_reduced"
    } else if (config$migration) {
      "AD_ptens_migration"
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
  cache_dir <- stan_model_cache_dir()
  message("CmdStan model cache dir: ", cache_dir)
  message("R tempdir(): ", tempdir())
  cmdstanr::cmdstan_model(
    stan_file = stan_file,
    dir = cache_dir
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
#' @param diag_model Diagnosis model type (1 = spline, 2 = AR(1)).
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
  penalty = 2,
  reduced_migration = FALSE
) {
  # Check if sim_diags object has been provided
  if ("stan_data" %in% names(stan_data)) {
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

  if (migration) {
    if (is.null(stan_data$detect_obs_sd)) {
      stan_data$detect_obs_sd <- 0
    }
    if (is.null(stan_data$detect_logit_obs)) {
      stan_data$detect_logit_obs <- numeric(0)
    }
    if (is.null(stan_data$detect_obs_index)) {
      if (length(stan_data$detect_logit_obs) == 0) {
        stan_data$detect_obs_index <- integer(0)
      } else if (length(stan_data$detect_logit_obs) == stan_data$nquar) {
        stan_data$detect_obs_index <- seq_len(stan_data$nquar)
      } else {
        stop(
          "detect_obs_index must be supplied when detect_logit_obs is sparse."
        )
      }
    }
    stan_data$detect_obs_index <- as.integer(stan_data$detect_obs_index)

    if (length(stan_data$detect_logit_obs) != length(stan_data$detect_obs_index)) {
      stop("detect_logit_obs and detect_obs_index must have the same length.")
    }
    if (length(stan_data$detect_obs_index) > 0 &&
        any(stan_data$detect_obs_index < 1 | stan_data$detect_obs_index > stan_data$nquar)) {
      stop("detect_obs_index entries must lie between 1 and nquar.")
    }

    # n_detect_obs is derived, not stored, so always computed fresh here.
    stan_data$n_detect_obs <- length(stan_data$detect_obs_index)

    # Quarters at which ascertainment of a previous diagnosis abroad changes
    # level (e.g. a surveillance reporting change). Supplied per surveillance
    # system; none by default.
    if (is.null(stan_data$detect_step_quarter)) {
      stan_data$detect_step_quarter <- integer(0)
    }
    stan_data$detect_step_quarter <- sort(unique(as.integer(
      stan_data$detect_step_quarter
    )))
    if (
      length(stan_data$detect_step_quarter) > 0 &&
        any(
          stan_data$detect_step_quarter < 2 |
            stan_data$detect_step_quarter > stan_data$nquar
        )
    ) {
      stop("detect_step_quarter entries must lie between 2 and nquar.")
    }
    stan_data$n_detect_steps <- length(stan_data$detect_step_quarter)

    # For age + migration models, pre-aggregate the single-year ABROAD_age into
    # the 4 diagnosis age groups used by the Stan multinomial likelihood.  This
    # is the sufficient statistic the model actually fits; the raw ABROAD_age is
    # retained in stan_data for R-side diagnostics but is excluded from the Stan
    # input (see fit_stan_model).  If ABROAD_group is already present (e.g.
    # supplied directly by a real-data pipeline) it is left unchanged.
    if (isTRUE(age) || isTRUE(stan_data$sim_config$age)) {
      if (is.null(stan_data$ABROAD_group)) {
        if (is.null(stan_data$ABROAD_age)) {
          stop("stan_data must contain ABROAD_age (or ABROAD_group) for age + migration models.")
        }
        nage_prep <- stan_data$nage
        agl <- pmin(4L, 1L + (seq_len(nage_prep) - 1L) %/% 10L)
        abroad_mat <- as.matrix(stan_data$ABROAD_age)
        grp <- matrix(0L, nrow = stan_data$nquar, ncol = 4L)
        for (ag in 1:4) {
          idx <- which(agl == ag)
          grp[, ag] <- as.integer(rowSums(abroad_mat[, idx, drop = FALSE]))
        }
        stan_data$ABROAD_group <- grp
      }
    }

    # gamma_obs: EXTERNAL auxiliary observations for log migration volume.
    # Support a dense scalar rolling-annual log series (gamma_obs) plus an
    # optional dense (nquar x 4) grouped matrix (gamma_obs_group) on the same
    # age groups as ABROAD_group. If gamma_obs_grouped is supplied, it is used
    # only to convert inactive grouped quarters into -1 sentinel rows.
    # Use only with evidence independent of the HIV/AIDS/CD4/ABROAD counts
    # (e.g. ONS migration flows). Do NOT pass the posterior of another model fit
    # to these same counts: that double-uses the data and yields an overconfident
    # posterior.
    if (is.null(stan_data$gamma_obs_sd)) {
      stan_data$gamma_obs_sd <- 0
    }
    if (is.null(stan_data$gamma_obs) && !is.null(stan_data$gamma_obs_mean)) {
      stan_data$gamma_obs <- stan_data$gamma_obs_mean
    }
    if (is.null(stan_data$gamma_obs_group) && !is.null(stan_data$gamma_obs_mean_group)) {
      stan_data$gamma_obs_group <- stan_data$gamma_obs_mean_group
    }
    if (is.null(stan_data$gamma_obs)) {
      stan_data$gamma_obs <- rep(0, stan_data$nquar)
    }
    if (length(stan_data$gamma_obs) == 0 && stan_data$gamma_obs_sd == 0) {
      stan_data$gamma_obs <- rep(0, stan_data$nquar)
    }
    if (length(stan_data$gamma_obs) != stan_data$nquar) {
      stop("gamma_obs must have length nquar.")
    }
    if (is.null(stan_data$gamma_obs_group)) {
      stan_data$gamma_obs_group <- matrix(-1, nrow = stan_data$nquar, ncol = 4L)
    }
    if (length(stan_data$gamma_obs_group) == 0 && stan_data$gamma_obs_sd == 0) {
      stan_data$gamma_obs_group <- matrix(-1, nrow = stan_data$nquar, ncol = 4L)
    }
    stan_data$gamma_obs_group <- as.matrix(stan_data$gamma_obs_group)
    if (!all(dim(stan_data$gamma_obs_group) == c(stan_data$nquar, 4L))) {
      stop("gamma_obs_group must have dimensions (nquar, 4).")
    }
    if (!is.null(stan_data$gamma_obs_grouped)) {
      if (is.logical(stan_data$gamma_obs_grouped)) {
        if (length(stan_data$gamma_obs_grouped) == 1L) {
          stan_data$gamma_obs_grouped <- rep(as.integer(stan_data$gamma_obs_grouped), stan_data$nquar)
        } else {
          stan_data$gamma_obs_grouped <- as.integer(stan_data$gamma_obs_grouped)
        }
      }
      if (length(stan_data$gamma_obs_grouped) != stan_data$nquar) {
        stop("gamma_obs_grouped must have length nquar.")
      }
      if (!all(stan_data$gamma_obs_grouped %in% c(0L, 1L))) {
        stop("gamma_obs_grouped must contain only 0/1 values.")
      }
      if (any(stan_data$gamma_obs_grouped == 1L)) {
        grp_rows <- stan_data$gamma_obs_group[stan_data$gamma_obs_grouped == 1L, , drop = FALSE]
        if (!all(rowSums(grp_rows >= 0) == 4L)) {
          stop("Rows with gamma_obs_grouped == 1 must provide all 4 gamma_obs_group values.")
        }
      }
      stan_data$gamma_obs_group[stan_data$gamma_obs_grouped == 0L, ] <- -1
    }
    gamma_group_available <- rowSums(stan_data$gamma_obs_group >= 0)
    if (!all(gamma_group_available %in% c(0L, 4L))) {
      stop("Each row of gamma_obs_group must be either all -1 or have all 4 age-group values.")
    }
    if (stan_data$gamma_obs_sd > 0 && any(stan_data$gamma_obs < 0)) {
      stop("gamma_obs must be provided for every quarter when gamma_obs_sd > 0.")
    }
    stan_data$gamma_obs_grouped <- NULL

    # Per-quarter random-walk sd on the log ratio of HIV-positive migration to
    # the external all-cause index. Set comparable to gamma_obs_sd so the index
    # still constrains short-run shape; the long-run change in that ratio is
    # carried by the model's own gamma_obs_slope.
    if (is.null(stan_data$gamma_obs_drift_sd)) {
      stan_data$gamma_obs_drift_sd <- 0.02
    }
    if (
      length(stan_data$gamma_obs_drift_sd) != 1L ||
        !is.finite(stan_data$gamma_obs_drift_sd) ||
        stan_data$gamma_obs_drift_sd < 0
    ) {
      stop("gamma_obs_drift_sd must be a single non-negative number.")
    }

    # Prior sd on the linear trend of the offset between gamma_obs and
    # log(omega); default 1 preserves the old std_normal() prior. Tighten
    # (e.g. towards 0) when the external index and HIV+ migration share are
    # expected to move together, so the offset cannot silently absorb a
    # trend that belongs to omega.
    if (is.null(stan_data$gamma_obs_slope_sd)) {
      stan_data$gamma_obs_slope_sd <- 1
    }
    if (
      length(stan_data$gamma_obs_slope_sd) != 1L ||
        !is.finite(stan_data$gamma_obs_slope_sd) ||
        stan_data$gamma_obs_slope_sd < 0
    ) {
      stop("gamma_obs_slope_sd must be a single non-negative number.")
    }
  }

  # Create stan_data$model_config
  stan_data$model_config <- stan_model_config(
    migration = migration,
    rita = rita,
    age = age,
    inf_model = inf_model,
    diag_model = diag_model,
    model_seed = model_seed,
    reduced_migration = reduced_migration
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
    if (!reduced_migration) {
      jagam_out_mp <- create_spline(knots - 4)
    }

    # Augment stan_data with design matrix and prior info
    stan_data$X_mig <- jagam_out_mig$jags.data$X
    stan_data$npars_mig <- ncol(stan_data$X_mig)

    stan_data$X_cob <- jagam_out_cob$jags.data$X
    stan_data$npars_cob <- ncol(stan_data$X_cob)

    if (!reduced_migration) {
      stan_data$X_mp <- if (stan_data$model_config$age) {
        jagam_out_mp$jags.data$X[, -1, drop = FALSE]
      } else {
        jagam_out_mp$jags.data$X
      }
      stan_data$npars_mp <- ncol(stan_data$X_mp)
    } else {
      stan_data$X_mp <- NULL
      stan_data$npars_mp <- NULL
    }
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
