test_that("cd4-only back-calculation produces a stanfit object", {
  # Simulate diagnoses
  sim_config <- define_sim_config(seed = 123, nquar = 30)
  sim_diags <- simulate_diagnoses(sim_config)
  # Basic CD4-only model can be fit with splines, random walk, or GP
  stan_config <- stan_model_config(
    inf_model = 1, # 1 spline, 2 random walk, 3 gp
    diag_model = 1, # 1 spline, 2 random walk
    sim_config = sim_config
  )
  # Prepare spline basis and augment stan_data
  stan_data <- model_type_prep(sim_diags$stan_data, stan_config)
  # Fit back-calculation model
  suppressWarnings({
    fit <- fit_backcalc(
      stan_data,
      stan_config,
      iter = 10,
      chains = 1,
      refresh = 0,
      messages = FALSE
    )
  })
  expect_s4_class(fit, "stanfit")
})

test_that("migration back-calculation produces a stanfit object", {
  sim_config <- define_sim_config(seed = 123, migration = TRUE)
  sim_diags <- simulate_diagnoses(sim_config)
  # Prepare spline basis and augment stan_data
  stan_data <- model_type_prep(sim_diags$stan_data, sim_config)
  # 1 chain and 10 iterations for testing
  # Suppress Stan warnings due to low number of iterations
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
  expect_s4_class(fit, "stanfit")
})

test_that("rita back-calculation produces a stanfit object", {
  sim_config <- define_sim_config(seed = 123, rita = TRUE)
  sim_diags <- simulate_diagnoses(sim_config)
  # Prepare spline basis and augment stan_data
  stan_data <- model_type_prep(sim_diags$stan_data, sim_config)
  # 1 chain and 10 iterations for testing
  # Suppress Stan warnings due to low number of iterations
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
  expect_s4_class(fit, "stanfit")
})

test_that("age-specific back-calculation runs without error", {
  sim_config <- define_sim_config(seed = 123, age = TRUE)
  sim_diags <- simulate_diagnoses(sim_config)
  # Prepare spline basis and augment stan_data
  stan_data <- model_type_prep(sim_diags$stan_data, sim_config)
  # 1 chain and 10 iterations for testing
  # Suppress Stan warnings due to low number of iterations
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
  expect_s4_class(fit, "stanfit")
})

test_that("age-specific rita back-calculation runs without error", {
  sim_config <- define_sim_config(seed = 123, age = TRUE, rita = TRUE)
  sim_diags <- simulate_diagnoses(sim_config)
  # Prepare spline basis and augment stan_data
  stan_data <- model_type_prep(sim_diags$stan_data, sim_config)
  # 1 chain and 10 iterations for testing
  # Suppress Stan warnings due to low number of iterations
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
  expect_s4_class(fit, "stanfit")
})
