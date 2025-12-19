test_that("plots are generated correctly", {
  # Simulate diagnoses
  stan_data <- simulate_diagnoses(nquar = 60, sim_seed = 123)
  # Fit back-calculation model
  suppressWarnings({
    hiv_list <- run_backcalc(
      stan_data,
      model_seed = 123,
      iter = 20,
      chains = 1,
      refresh = 0,
      messages = FALSE,
      nworkers = 1
    )
  })
  # plot diagnostics
  p0 <- plot_diagnostics(hiv_list)
  expect_s3_class(p0, "ggplot")
  # Generate plots of key quantities
  p1 <- plot_estimates(hiv_list, "incidence")
  p2 <- plot_estimates(hiv_list, "undiag_prev")
  p3 <- plot_estimates(hiv_list, "diag_prob")
  p4 <- plot_estimates(hiv_list, "diagnoses")
  expect_s3_class(p1, "ggplot")
  expect_s3_class(p2, "ggplot")
  expect_s3_class(p3, "ggplot")
  expect_s3_class(p4, "ggplot")
})

test_that("plots are generated correctly for migration model", {
  # Simulate diagnoses
  stan_data <- simulate_diagnoses(nquar = 60, sim_seed = 123, migration = TRUE)
  # Fit back-calculation model
  suppressWarnings({
    hiv_list <- run_backcalc(
      stan_data,
      model_seed = 123,
      iter = 20,
      chains = 1,
      refresh = 0,
      messages = FALSE,
      nworkers = 1
    )
  })
  # plot diagnostics
  p0 <- plot_diagnostics(hiv_list)
  expect_s3_class(p0, "ggplot")
  # Generate plots of key quantities
  p1 <- plot_estimates(hiv_list, "migration")
  p2 <- plot_estimates(hiv_list, "migration_prob")
  p3 <- plot_estimates(hiv_list, "transfers")
  p4 <- plot_estimates(hiv_list, "prop_uk")
  p5 <- plot_estimates(hiv_list, "detect_prob")
  expect_s3_class(p1, "ggplot")
  expect_s3_class(p2, "ggplot")
  expect_s3_class(p3, "ggplot")
  expect_s3_class(p4, "ggplot")
  expect_s3_class(p5, "ggplot")
})

test_that("plots work for lists", {
  # Simulate diagnoses
  stan_data <- simulate_diagnoses(nquar = 60, sim_seed = 123, migration = TRUE)
  # Fit back-calculation model
  hiv_list <- list()
  for (i in 1:2) {
    suppressWarnings({
      hiv_list[[i]] <- run_backcalc(
        stan_data,
        model_seed = 123 + i,
        iter = 20,
        chains = 1,
        refresh = 0,
        messages = FALSE,
        nworkers = 1
      )
    })
  }
  # plot diagnostics
  p0 <- plot_diagnostics(hiv_list)
  expect_s3_class(p0, "ggplot")
  # Generate plots of key quantities
  p1 <- plot_estimates(hiv_list, "incidence")
  p2 <- plot_estimates(hiv_list, "undiag_prev")
  p3 <- plot_estimates(hiv_list, "diag_prob")
  p4 <- plot_estimates(hiv_list, "diagnoses")
  expect_s3_class(p1, "ggplot")
  expect_s3_class(p2, "ggplot")
  expect_s3_class(p3, "ggplot")
  expect_s3_class(p4, "ggplot")
    # Generate plots of key quantities
  p5 <- plot_estimates(hiv_list, "migration")
  p6 <- plot_estimates(hiv_list, "migration_prob")
  p7 <- plot_estimates(hiv_list, "transfers")
  p8 <- plot_estimates(hiv_list, "prop_uk")
  p9 <- plot_estimates(hiv_list, "detect_prob")
  expect_s3_class(p5, "ggplot")
  expect_s3_class(p6, "ggplot")
  expect_s3_class(p7, "ggplot")
  expect_s3_class(p8, "ggplot")
  expect_s3_class(p9, "ggplot")
})