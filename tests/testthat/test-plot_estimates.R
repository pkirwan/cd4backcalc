test_that("simulated plots are generated correctly", {
  # Simulate diagnoses
  stan_data <- simulate_diagnoses(nquar = 40, sim_seed = 123)
  # Generate plots of key quantities
  p1 <- plot_simulations(stan_data, "incidence")
  p2 <- plot_simulations(stan_data, "undiag_prev")
  p3 <- plot_simulations(stan_data, "diag_prob")
  p4 <- plot_simulations(stan_data, "diagnoses")
  expect_s3_class(p1, "ggplot")
  expect_s3_class(p2, "ggplot")
  expect_s3_class(p3, "ggplot")
  expect_s3_class(p4, "ggplot")
})

test_that("simulated plots are generated correctly for migration_prob", {
  # Simulate diagnoses
  stan_data <- simulate_diagnoses(nquar = 40, sim_seed = 123, migration = TRUE, m_pattern = "varying")
  # Generate plots of key quantities
  p1 <- plot_simulations(stan_data, "migration_prob")
  expect_s3_class(p1, "ggplot")
  stan_data <- simulate_diagnoses(nquar = 40, sim_seed = 123, migration = TRUE, rita = TRUE, m_pattern = "varying")
  # Generate plots of key quantities
  p2 <- plot_simulations(stan_data, "migration_prob")
  expect_s3_class(p2, "ggplot")
})

test_that("simulated age plots are generated correctly", {
  # Simulate diagnoses
  stan_data <- simulate_diagnoses(nquar = 40, sim_seed = 123, age = TRUE)
  # Generate plots of key quantities
  p1 <- plot_simulations(stan_data, "incidence_age")
  p2 <- plot_simulations(stan_data, "undiag_prev_age")
  expect_s3_class(p1, "ggplot")
  expect_s3_class(p2, "ggplot")
})

test_that("plots are generated correctly", {
  # Simulate diagnoses
  stan_data <- simulate_diagnoses(nquar = 40, sim_seed = 123)
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

test_that("chained plots are generated correctly", {
  # Simulate diagnoses
  stan_data <- simulate_diagnoses(nquar = 40, sim_seed = 123)
  # Fit back-calculation model
  suppressWarnings({
    hiv_list <- run_backcalc(
      stan_data,
      model_seed = 123,
      iter = 20,
      chains = 2,
      refresh = 0,
      messages = FALSE,
      nworkers = 1,
      chain_results = TRUE
    )
  })
  # Generate plots of key quantities
  p1 <- plot_estimates(hiv_list, "incidence", chain_results = TRUE)
  expect_s3_class(p1, "ggplot")
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
  p1 <- plot_estimates(hiv_list, "undiag_migration")
  p2 <- plot_estimates(hiv_list, "migration_prob")
  p3 <- plot_estimates(hiv_list, "transfers")
  p4 <- plot_estimates(hiv_list, "ratio_abroad_uk")
  p5 <- plot_estimates(hiv_list, "detect_prob")
  expect_s3_class(p1, "ggplot")
  expect_s3_class(p2, "ggplot")
  expect_s3_class(p3, "ggplot")
  expect_s3_class(p4, "ggplot")
  expect_s3_class(p5, "ggplot")
})

test_that("plots work for migration model lists", {
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
  p5 <- plot_estimates(hiv_list, "undiag_migration")
  p6 <- plot_estimates(hiv_list, "migration_prob")
  p7 <- plot_estimates(hiv_list, "transfers")
  p8 <- plot_estimates(hiv_list, "ratio_abroad_uk")
  p9 <- plot_estimates(hiv_list, "detect_prob")
  expect_s3_class(p5, "ggplot")
  expect_s3_class(p6, "ggplot")
  expect_s3_class(p7, "ggplot")
  expect_s3_class(p8, "ggplot")
  expect_s3_class(p9, "ggplot")
})

test_that("plots work for rita model lists", {
  # Simulate diagnoses
  stan_data <- simulate_diagnoses(nquar = 40, sim_seed = 123, rita = TRUE)
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
})


test_that("plots work for migration and rita model lists", {
  # Simulate diagnoses
  stan_data <- simulate_diagnoses(nquar = 40, sim_seed = 123, rita = TRUE, migration = TRUE, m_pattern = "varying")
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
  p1 <- bias_plot(hiv_list, "incidence")
  p2 <- plot_estimates(hiv_list, "migration_prob")
  expect_s3_class(p1, "ggplot")
  expect_s3_class(p2, "ggplot")
})

test_that("plots work for age-dependent model", {
  # Simulate diagnoses
  stan_data <- simulate_diagnoses(nquar = 40, sim_seed = 123, age = TRUE)
  # Fit back-calculation model
  suppressWarnings({
    hiv_list <- run_backcalc(
      stan_data,
      model_seed = 123,
      iter = 10,
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
  p3 <- plot_estimates(hiv_list, "diagnoses")
  p4 <- plot_estimates(hiv_list, "incidence_age")
  p5 <- plot_estimates(hiv_list, "undiag_prev_age")
  p6 <- plot_estimates(hiv_list, "diagnoses_age")
  p7 <- plot_estimates(hiv_list, "diag_prob_age")
  expect_s3_class(p1, "ggplot")
  expect_s3_class(p2, "ggplot")
  expect_s3_class(p3, "ggplot")
  expect_s3_class(p4, "ggplot")
  expect_s3_class(p5, "ggplot")
  expect_s3_class(p6, "ggplot")
  expect_s3_class(p7, "ggplot")
})

test_that("plots work for age-dependent rita model", {
  # Simulate diagnoses
  stan_data <- simulate_diagnoses(nquar = 40, sim_seed = 123, age = TRUE, rita = TRUE)
  # Fit back-calculation model
  suppressWarnings({
    hiv_list <- run_backcalc(
      stan_data,
      model_seed = 123,
      iter = 10,
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
  p3 <- plot_estimates(hiv_list, "diagnoses")
  p4 <- plot_estimates(hiv_list, "incidence_age")
  p5 <- plot_estimates(hiv_list, "undiag_prev_age")
  p6 <- plot_estimates(hiv_list, "diagnoses_age")
  p7 <- plot_estimates(hiv_list, "diag_prob_age")
  expect_s3_class(p1, "ggplot")
  expect_s3_class(p2, "ggplot")
  expect_s3_class(p3, "ggplot")
  expect_s3_class(p4, "ggplot")
  expect_s3_class(p5, "ggplot")
  expect_s3_class(p6, "ggplot")
  expect_s3_class(p7, "ggplot")
})


test_that("start_year and overlay of real values works", {
  # Simulate diagnoses
  stan_data <- simulate_diagnoses(nquar = 40, sim_seed = 123, age = TRUE, rita = TRUE, migration = TRUE)
  # Fit back-calculation model
  suppressWarnings({
    hiv_list <- run_backcalc(
      stan_data,
      model_seed = 123,
      iter = 10,
      chains = 1,
      refresh = 0,
      messages = FALSE,
      nworkers = 1
    )
  })
  # Remove sim_data element from list to mimic real data scenario
  hiv_list <- hiv_list[names(hiv_list) != "sim_data"]
  # Generate plots of key quantities
  p1 <- plot_estimates(hiv_list, "diagnoses", start_year = 2000)
  p2 <- plot_estimates(hiv_list, "diagnoses_age", start_year = 2000)
  p3 <- plot_estimates(hiv_list, "transfers", start_year = 2000)
  expect_s3_class(p1, "ggplot")
  expect_s3_class(p2, "ggplot")
  expect_s3_class(p3, "ggplot")
})