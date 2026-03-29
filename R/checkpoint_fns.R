#' Create initial values for Stan from a fitted model
#'
#' This function extracts the last iteration's parameter values from a fitted Stan model.
#' It returns a list of initial values for each chain, which can be used to initialize a new Stan model.
#'
#' @param fit A fitted Stan model object (e.g., from rstan or cmdstanr).
#' @return A list of initial values for each chain.
#' @export
create_stan_inits <- function(fit) {
  # Try to include warmup draws (needed when iter_sampling = 0).

  # Fall back to sampling-only draws when warmup was not saved.
  draws <- tryCatch(
    as.array(fit$draws(inc_warmup = TRUE)),
    error = function(e) as.array(fit$draws())
  )
  n_iter <- dim(draws)[1]
  n_chains <- dim(draws)[2]
  all_vars <- dimnames(draws)[[3]]

  # Exclude internal Stan variables
  param_vars <- all_vars[!grepl("^lp__|__$", all_vars)]

  # Parse base names and index strings once
  base_names <- sub("\\[.*$", "", param_vars)
  idx_strs <- ifelse(
    grepl("\\[", param_vars),
    sub("^[^\\[]*\\[(.*)\\]$", "\\1", param_vars),
    NA_character_
  )

  lapply(seq_len(n_chains), function(chain_id) {
    vals <- as.numeric(draws[n_iter, chain_id, param_vars])
    unique_params <- unique(base_names)

    result <- lapply(unique_params, function(pname) {
      mask <- base_names == pname
      pvals <- vals[mask]
      pidx <- idx_strs[mask]

      # Scalar parameter
      if (length(pvals) == 1 && is.na(pidx[1])) {
        return(pvals[1])
      }

      indices <- strsplit(pidx, ",")
      ndim <- length(indices[[1]])

      if (ndim == 1) {
        # Vector
        idx <- as.integer(unlist(indices))
        vec <- numeric(max(idx))
        vec[idx] <- pvals
        return(vec)
      }

      # Multi-dimensional array
      idx_mat <- do.call(rbind, lapply(indices, as.integer))
      dims <- apply(idx_mat, 2, max)
      arr <- array(NA_real_, dim = dims)
      for (i in seq_len(nrow(idx_mat))) {
        arr[matrix(idx_mat[i, ], nrow = 1)] <- pvals[i]
      }
      arr
    })
    names(result) <- unique_params
    result
  })
}

#' Save a checkpoint during Stan model fitting
#'
#' Saves the current state of a Stan model fit so that sampling can be
#' resumed in a subsequent run. The checkpoint includes the adaptation
#' state (inverse metric and step size), the last iteration's parameter
#' values (for chain initialisation), and any accumulated draws.
#'
#' @param checkpoint_dir Directory in which to save the checkpoint.
#' @param fit A CmdStanFit object from the most recent sampling call.
#' @param phase Current phase: \code{"warmup"} or \code{"sampling"}.
#' @param warmup_done Total warmup iterations completed so far.
#' @param sampling_done Total sampling iterations completed so far.
#' @param draws_list List of \code{posterior::draws_array} objects from
#'   completed sampling chunks.
#' @param diagnostics_list List of diagnostic summaries from each chunk.
#' @param total_time Cumulative wall-clock time in seconds across all chunks.
#' @return The checkpoint list (invisibly).
#' @noRd
save_checkpoint <- function(checkpoint_dir, fit, phase, warmup_done,
                            sampling_done, draws_list = list(),
                            diagnostics_list = list(), total_time = 0) {
  dir.create(checkpoint_dir, recursive = TRUE, showWarnings = FALSE)
  checkpoint <- list(
    phase = phase,
    warmup_done = warmup_done,
    sampling_done = sampling_done,
    inv_metric = fit$inv_metric(matrix = FALSE),
    step_size = fit$metadata()$step_size_adaptation,
    init = create_stan_inits(fit),
    draws_list = draws_list,
    diagnostics_list = diagnostics_list,
    total_time = total_time
  )
  checkpoint_file <- file.path(checkpoint_dir, "checkpoint.RData")
  # Write to temp file then rename for atomicity
  tmp <- tempfile(tmpdir = checkpoint_dir, fileext = ".RData")
  save(checkpoint, file = tmp)
  file.rename(tmp, checkpoint_file)
  invisible(checkpoint)
}

#' Load a checkpoint from disk
#'
#' @param checkpoint_dir Directory containing a \code{checkpoint.RData} file.
#' @return The checkpoint list, or \code{NULL} if no checkpoint exists.
#' @noRd
load_checkpoint <- function(checkpoint_dir) {
  checkpoint_file <- file.path(checkpoint_dir, "checkpoint.RData")
  if (!file.exists(checkpoint_file)) return(NULL)
  env <- new.env(parent = emptyenv())
  load(checkpoint_file, envir = env)
  env$checkpoint
}

#' Create a combined fit-like object from multiple sampling chunks
#'
#' After chunked sampling, this function binds draws from all chunks and
#' returns an object with the same interface as a CmdStanFit (the subset
#' used by \code{\link{postproc}} and \code{\link{plot_diagnostics}}).
#'
#' @param last_fit The CmdStanFit from the final sampling chunk.
#' @param draws_list List of \code{posterior::draws_array} objects.
#' @param diagnostics_list List of diagnostic summaries from each chunk.
#' @param n_chains Number of chains.
#' @param total_time Cumulative wall-clock time in seconds.
#' @return An object of class \code{cd4backcalc_combined_fit}.
#' @noRd
create_combined_fit <- function(last_fit, draws_list, diagnostics_list,
                                n_chains, total_time) {
  all_draws <- posterior::bind_draws(draws_list, along = "iteration")

  # Aggregate diagnostics: sum counts, keep last EBFMI as approximation
  agg_diag <- diagnostics_list[[1]]
  if (length(diagnostics_list) > 1) {
    for (i in 2:length(diagnostics_list)) {
      d <- diagnostics_list[[i]]
      agg_diag$num_divergent <- agg_diag$num_divergent + d$num_divergent
      agg_diag$num_max_treedepth <-
        agg_diag$num_max_treedepth + d$num_max_treedepth
      agg_diag$ebfmi <- d$ebfmi
    }
  }

  structure(
    list(
      draws = function(format = "array", ...) {
        if (format == "matrix") {
          posterior::as_draws_matrix(all_draws)
        } else if (format == "df" || format == "data.frame") {
          posterior::as_draws_df(all_draws)
        } else {
          all_draws
        }
      },
      num_chains = function() n_chains,
      diagnostic_summary = function(...) agg_diag,
      time = function() list(total = total_time),
      cmdstan_diagnose = function() last_fit$cmdstan_diagnose(),
      inv_metric = function(...) last_fit$inv_metric(...),
      metadata = function() last_fit$metadata()
    ),
    class = "cd4backcalc_combined_fit"
  )
}
