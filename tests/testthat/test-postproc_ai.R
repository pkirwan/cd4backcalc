test_that("cd4-only post-processing works", {
  # Simulate diagnoses
  sim_config <- define_sim_config(nquar = 60, seed = 123)
  sim_diags <- simulate_diagnoses(sim_config)
  # Prepare spline basis and augment stan_data
  stan_data <- model_type_prep(sim_diags$stan_data, sim_config)
  # Fit back-calculation model
  suppressWarnings({
    fit <- fit_backcalc(
      stan_data,
      sim_config,
      iter = 20,
      chains = 1,
      refresh = 0,
      messages = FALSE
    )
  })
  # Post-process fitted model to obtain posterior predictive
  hiv_list <- postproc(
    fit,
    stan_data,
    sim_config,
    start_yr = 0,
    nworkers = 1
  )
  expect_true("diag_prob" %in% names(hiv_list))
})

test_that("migration post-processing works", {
  # Simulate diagnoses
  sim_config <- define_sim_config(nquar = 60, migration = TRUE, seed = 123)
  sim_diags <- simulate_diagnoses(sim_config)
  # Prepare spline basis and augment stan_data
  stan_data <- model_type_prep(sim_diags$stan_data, sim_config)
  # Fit back-calculation model
  suppressWarnings({
    fit <- fit_backcalc(
      stan_data,
      sim_config,
      iter = 20,
      chains = 1,
      refresh = 0,
      messages = FALSE
    )
  })
  # Post-process fitted model to obtain posterior predictive
  hiv_list <- postproc(
    fit,
    stan_data,
    sim_config,
    start_yr = 0,
    nworkers = 1
  )
  expect_true("migration_prob" %in% names(hiv_list))
})

test_that("rita post-processing works", {
  # Simulate diagnoses
  sim_config <- define_sim_config(nquar = 60, rita = TRUE, seed = 123)
  sim_diags <- simulate_diagnoses(sim_config)
  # Prepare spline basis and augment stan_data
  stan_data <- model_type_prep(sim_diags$stan_data, sim_config)
  # Fit back-calculation model
  suppressWarnings({
    fit <- fit_backcalc(
      stan_data,
      sim_config,
      iter = 20,
      chains = 1,
      refresh = 0,
      messages = FALSE
    )
  })
  # Post-process fitted model to obtain posterior predictive
  hiv_list <- postproc(
    fit,
    stan_data,
    sim_config,
    start_yr = 0,
    nworkers = 1
  )
  expect_true("diag_prob" %in% names(hiv_list))
})

test_that("age-specific post-processing works", {
  sim_config <- define_sim_config(seed = 123, nquar = 60, age = TRUE)
  sim_diags <- simulate_diagnoses(sim_config)
  # Prepare spline basis and augment stan_data
  stan_data <- model_type_prep(sim_diags$stan_data, sim_config)
  # Fit back-calculation model
  suppressWarnings({
    fit <- fit_backcalc(
      stan_data,
      sim_config,
      iter = 10,
      chains = 1,
      refresh = 0,
      messages = FALSE
    )
  })
  # Post-process fitted model to obtain posterior predictive
  hiv_list <- postproc(
    fit,
    stan_data,
    sim_config,
    start_yr = 0,
    nworkers = 1
  )
  expect_true("infs.qt" %in% names(hiv_list))
})

test_that("age-specific rita back-calculation runs without error", {
  sim_config <- define_sim_config(seed = 123, nquar = 60, age = TRUE, rita = TRUE)
  sim_diags <- simulate_diagnoses(sim_config)
  # Prepare spline basis and augment stan_data
  stan_data <- model_type_prep(sim_diags$stan_data, sim_config)
  # Fit back-calculation model
  suppressWarnings({
    fit <- fit_backcalc(
      stan_data,
      sim_config,
      iter = 10,
      chains = 1,
      refresh = 1,
      messages = TRUE
    )
  })
  # Post-process fitted model to obtain posterior predictive
  hiv_list <- postproc(
    fit,
    stan_data,
    sim_config,
    start_yr = 0,
    nworkers = 1
  )
  expect_true("infs.qt" %in% names(hiv_list))
})