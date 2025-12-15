test_that("mse_list is generated correctly", {
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
  # Generate plots
  model_list <- extract_simulated_quantities(sim_diags, sim_config)
  expect_true(model_list$incidence[1,1] == 1478)
  mse_list <- make_mse_list(hiv_list, sim_diags, sim_config)
  expect_true(mse_list$incidence_bias[1,8] == -1546)
})

test_that("plots are generated correctly", {
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
  # Generate plots
  mse_list <- make_mse_list(hiv_list, sim_diags, sim_config)
  p1 <- mse_plot(mse_list$incidence_mse)
  expect_s3_class(p1, "ggplot")
  p2 <- bias_plot(mse_list$incidence_bias)
  expect_s3_class(p2, "ggplot")
})
