test_that("save_checkpoint / load_checkpoint round-trips correctly", {
  skip_if_no_cmdstan()
  sim_diags <- simulate_diagnoses(sim_seed = 42)
  stan_data <- model_type_prep(sim_diags, model_seed = 42)
  checkpoint_dir <- file.path(tempdir(), "cp_roundtrip")
  on.exit(unlink(checkpoint_dir, recursive = TRUE), add = TRUE)

  # Run a tiny warmup-only fit to get a real fit object
  suppressWarnings({
    fit <- fit_stan_model(
      stan_data,
      iter_warmup = 5,
      iter_sampling = 5,
      chains = 1,
      refresh = 0,
      messages = FALSE
    )
  })

  save_checkpoint(
    checkpoint_dir, fit, "sampling",
    warmup_done = 5, sampling_done = 3,
    draws_list = list(fit$draws()),
    diagnostics_list = list(fit$diagnostic_summary()),
    total_time = 1.5
  )

  cp <- load_checkpoint(checkpoint_dir)
  expect_type(cp, "list")
  expect_equal(cp$phase, "sampling")
  expect_equal(cp$warmup_done, 5)
  expect_equal(cp$sampling_done, 3)
  expect_type(cp$inv_metric, "list")
  expect_type(cp$step_size, "double")
  expect_type(cp$init, "list")
  expect_length(cp$draws_list, 1)
  expect_length(cp$diagnostics_list, 1)
  expect_equal(cp$total_time, 1.5)
})

test_that("load_checkpoint returns NULL when no checkpoint exists", {
  dir <- file.path(tempdir(), "cp_nonexistent")
  expect_null(load_checkpoint(dir))
})

# ---- create_stan_inits ----

test_that("create_stan_inits returns one init list per chain", {
  skip_if_no_cmdstan()
  sim_diags <- simulate_diagnoses(sim_seed = 42)
  stan_data <- model_type_prep(sim_diags, model_seed = 42)
  suppressWarnings({
    fit <- fit_stan_model(
      stan_data,
      iter_warmup = 5,
      iter_sampling = 5,
      chains = 2,
      refresh = 0,
      messages = FALSE
    )
  })

  inits <- create_stan_inits(fit)
  expect_type(inits, "list")
  expect_length(inits, 2)
  # Each element should be a named list of parameter values
  expect_type(inits[[1]], "list")
  expect_true(length(names(inits[[1]])) > 0)
  # No internal Stan variables should be present
  expect_false(any(grepl("lp__", names(inits[[1]]))))
})

# ---- fit_stan_model with checkpointing ----

test_that("fit_stan_model with checkpoint_dir produces valid fit", {
  skip_if_no_cmdstan()
  sim_diags <- simulate_diagnoses(sim_seed = 42)
  stan_data <- model_type_prep(sim_diags, model_seed = 42)
  checkpoint_dir <- file.path(tempdir(), "cp_fit_test")
  on.exit(unlink(checkpoint_dir, recursive = TRUE), add = TRUE)

  suppressWarnings({
    fit <- fit_stan_model(
      stan_data,
      iter_warmup = 5,
      iter_sampling = 10,
      chains = 1,
      refresh = 0,
      messages = FALSE,
      checkpoint_dir = checkpoint_dir,
      iter_per_chunk = 5
    )
  })

  # Should return a combined fit since iter_per_chunk < iter_sampling
  expect_s3_class(fit, "cd4backcalc_combined_fit")

  # draws() should work in both formats
  draws_mat <- fit$draws(format = "matrix")
  expect_true(is.matrix(draws_mat))
  expect_true(nrow(draws_mat) > 0)

  draws_arr <- fit$draws()
  expect_s3_class(draws_arr, "draws_array")

  # Other interface methods
  expect_equal(fit$num_chains(), 1)
  expect_type(fit$diagnostic_summary(), "list")
  expect_type(fit$time(), "list")
  expect_true(fit$time()$total > 0)
})

test_that("checkpointed fit works with postproc", {
  skip_if_no_cmdstan()
  sim_diags <- simulate_diagnoses(sim_seed = 42)
  stan_data <- model_type_prep(sim_diags, model_seed = 42)
  checkpoint_dir <- file.path(tempdir(), "cp_postproc_test")
  on.exit(unlink(checkpoint_dir, recursive = TRUE), add = TRUE)

  suppressWarnings({
    fit <- fit_stan_model(
      stan_data,
      iter_warmup = 5,
      iter_sampling = 10,
      chains = 1,
      refresh = 0,
      messages = FALSE,
      checkpoint_dir = checkpoint_dir,
      iter_per_chunk = 5
    )
  })

  hiv_list <- postproc(fit, stan_data, start_yr = 0, nworkers = 1)
  expect_true("diag_prob" %in% names(hiv_list))
})

# ---- Resumption ----

test_that("fit_stan_model resumes from a checkpoint", {
  skip_if_no_cmdstan()
  sim_diags <- simulate_diagnoses(sim_seed = 42)
  stan_data <- model_type_prep(sim_diags, model_seed = 42)
  checkpoint_dir <- file.path(tempdir(), "cp_resume_test")
  on.exit(unlink(checkpoint_dir, recursive = TRUE), add = TRUE)

  # Run first chunk only (5 of 10 sampling iters)
  suppressWarnings({
    fit1 <- fit_stan_model(
      stan_data,
      iter_warmup = 5,
      iter_sampling = 5,
      chains = 1,
      refresh = 0,
      messages = FALSE,
      checkpoint_dir = checkpoint_dir,
      iter_per_chunk = 5
    )
  })

  # Checkpoint should exist
  cp <- load_checkpoint(checkpoint_dir)
  expect_equal(cp$phase, "sampling")
  expect_equal(cp$sampling_done, 5)

  # Now resume asking for 10 total — should pick up where we left off
  # and run another 5
  suppressWarnings({
    fit2 <- fit_stan_model(
      stan_data,
      iter_warmup = 5,
      iter_sampling = 10,
      chains = 1,
      refresh = 0,
      messages = FALSE,
      checkpoint_dir = checkpoint_dir,
      iter_per_chunk = 5
    )
  })

  # Should have combined draws from both chunks
  expect_s3_class(fit2, "cd4backcalc_combined_fit")
  draws <- fit2$draws(format = "matrix")
  expect_equal(nrow(draws), 10)
})

# ---- run_backcalc with checkpointing (end-to-end) ----

test_that("run_backcalc with checkpointing produces full results", {
  skip_if_no_cmdstan()
  sim_diags <- simulate_diagnoses(sim_seed = 42)
  stan_data <- model_type_prep(sim_diags, model_seed = 42)
  checkpoint_dir <- file.path(tempdir(), "cp_e2e_test")
  on.exit(unlink(checkpoint_dir, recursive = TRUE), add = TRUE)

  suppressWarnings({
    hiv_list <- run_backcalc(
      stan_data,
      model_seed = 42,
      iter_warmup = 10,
      iter_sampling = 10,
      chains = 1,
      refresh = 0,
      messages = FALSE,
      nworkers = 1,
      checkpoint_dir = checkpoint_dir,
      iter_per_chunk = 5
    )
  })

  expect_true("diag_prob" %in% names(hiv_list))
  expect_s3_class(hiv_list$fit, "cd4backcalc_combined_fit")

  # Cleanup happened to leave checkpoint behind (caller responsibility)
  expect_true(file.exists(file.path(checkpoint_dir, "checkpoint.RData")))
})

# ---- Single chunk returns CmdStanMCMC directly ----

test_that("checkpoint with single chunk returns CmdStanMCMC", {
  skip_if_no_cmdstan()
  sim_diags <- simulate_diagnoses(sim_seed = 42)
  stan_data <- model_type_prep(sim_diags, model_seed = 42)
  checkpoint_dir <- file.path(tempdir(), "cp_single_chunk")
  on.exit(unlink(checkpoint_dir, recursive = TRUE), add = TRUE)

  suppressWarnings({
    fit <- fit_stan_model(
      stan_data,
      iter_warmup = 5,
      iter_sampling = 5,
      chains = 1,
      refresh = 0,
      messages = FALSE,
      checkpoint_dir = checkpoint_dir,
      iter_per_chunk = 10
    )
  })

  # Single chunk — should return the raw CmdStanMCMC, not a combined fit
  expect_r6_class(fit, "CmdStanMCMC")
})
