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
  # plot diagnostics
  p0 <- plot_diagnostics(fit)
  expect_s3_class(p0, "ggplot")
  # Generate plots of key quantities
  p1 <- plot_estimates(hiv_list, sim_diags, sim_config, "incidence")
  p2 <- plot_estimates(hiv_list, sim_diags, sim_config, "undiag_prev")
  p3 <- plot_estimates(hiv_list, sim_diags, sim_config, "diag_prob")
  p4 <- plot_estimates(hiv_list, sim_diags, sim_config, "diagnoses")
  expect_s3_class(p1, "ggplot")
  expect_s3_class(p2, "ggplot")
  expect_s3_class(p3, "ggplot")
  expect_s3_class(p4, "ggplot")
})

test_that("plots are generated correctly for migration model", {
  # Simulate diagnoses
  sim_config <- define_sim_config(nquar = 60, seed = 123, migration = TRUE)
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
  # plot diagnostics
  p0 <- plot_diagnostics(fit)
  expect_s3_class(p0, "ggplot")
  # Generate plots of key quantities
  p1 <- plot_estimates(hiv_list, sim_diags, sim_config, "migration")
  p2 <- plot_estimates(hiv_list, sim_diags, sim_config, "migration_prob")
  p3 <- plot_estimates(hiv_list, sim_diags, sim_config, "transfers")
  p4 <- plot_estimates(hiv_list, sim_diags, sim_config, "prop_uk")
  p5 <- plot_estimates(hiv_list, sim_diags, sim_config, "detect_prob")
  expect_s3_class(p1, "ggplot")
  expect_s3_class(p2, "ggplot")
  expect_s3_class(p3, "ggplot")
  expect_s3_class(p4, "ggplot")
  expect_s3_class(p5, "ggplot")
})
