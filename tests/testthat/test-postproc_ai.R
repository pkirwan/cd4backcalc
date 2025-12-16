test_that("cd4-only post-processing works", {
  # Simulate diagnoses
  stan_data <- simulate_diagnoses(nquar = 60, seed = 123)
  # Fit back-calculation model
  suppressWarnings({
    fit <- fit_stan_model(
      stan_data,
      iter = 20,
      chains = 1,
      refresh = 0,
      messages = FALSE
    )
  })
  expect_s4_class(fit, "stanfit")

  # Post-process fitted model to obtain posterior predictive
  hiv_list <- postproc(
    fit,
    stan_data,
    start_yr = 0,
    nworkers = 1
  )
  expect_true("diag_prob" %in% names(hiv_list))
})

test_that("migration post-processing works", {
  # Simulate diagnoses
  stan_data <- simulate_diagnoses(nquar = 60, migration = TRUE, seed = 123)
  # Fit back-calculation model
  suppressWarnings({
    fit <- fit_stan_model(
      stan_data,
      iter = 20,
      chains = 1,
      refresh = 0,
      messages = FALSE
    )
  })
  expect_s4_class(fit, "stanfit")

  # Post-process fitted model to obtain posterior predictive
  hiv_list <- postproc(
    fit,
    stan_data,
    start_yr = 0,
    nworkers = 1
  )
  expect_true("migration_prob" %in% names(hiv_list))
})

test_that("rita post-processing works", {
  # Simulate diagnoses
  stan_data <- simulate_diagnoses(nquar = 60, rita = TRUE, seed = 123)
  # Fit back-calculation model
  suppressWarnings({
    fit <- fit_stan_model(
      stan_data,
      iter = 20,
      chains = 1,
      refresh = 0,
      messages = FALSE
    )
  })
  expect_s4_class(fit, "stanfit")

  # Post-process fitted model to obtain posterior predictive
  hiv_list <- postproc(
    fit,
    stan_data,
    start_yr = 0,
    nworkers = 1
  )
  expect_true("diag_prob" %in% names(hiv_list))
})

test_that("age-specific post-processing works", {
  stan_data <- simulate_diagnoses(seed = 123, nquar = 60, age = TRUE)
  # Fit back-calculation model
  suppressWarnings({
    fit <- fit_stan_model(
      stan_data,
      iter = 2,
      chains = 1,
      refresh = 0,
      messages = FALSE
    )
  })
  expect_s4_class(fit, "stanfit")

  # Post-process fitted model to obtain posterior predictive
  hiv_list <- postproc(
    fit,
    stan_data,
    start_yr = 0,
    nworkers = 1
  )
  expect_true("infs.qt" %in% names(hiv_list))
})

test_that("age-specific rita back-calculation runs without error", {
  stan_data <- simulate_diagnoses(
    seed = 123,
    nquar = 60,
    age = TRUE,
    rita = TRUE
  )
  # Fit back-calculation model
  suppressWarnings({
    fit <- fit_stan_model(
      stan_data,
      iter = 10,
      chains = 1,
      refresh = 0,
      messages = FALSE
    )
  })
  expect_s4_class(fit, "stanfit")

  # Post-process fitted model to obtain posterior predictive
  hiv_list <- postproc(
    fit,
    stan_data,
    start_yr = 0,
    nworkers = 1
  )
  expect_true("infs.qt" %in% names(hiv_list))
})
