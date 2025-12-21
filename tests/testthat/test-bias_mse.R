test_that("extract_simulated_quantities correctly extracts values ai and ad", {
  # Simulate diagnoses
  stan_data <- simulate_diagnoses(nquar = 60, sim_seed = 123)
  sim_estimates <- extract_simulated_quantities(stan_data)
  expect_equal(
    sim_estimates$incidence$value[1],
    stan_data$sim_data$infs[1, 1]
  )

  stan_data <- simulate_diagnoses(nquar = 60, sim_seed = 123, migration = TRUE)
  sim_estimates <- extract_simulated_quantities(stan_data)
  expect_equal(
    sim_estimates$migration$value[1],
    sum(stan_data$sim_data$m[1, 1:4])
  )

  stan_data <- simulate_diagnoses(nquar = 60, sim_seed = 123, age = TRUE)
  sim_estimates <- extract_simulated_quantities(stan_data)
  expect_equal(
    sim_estimates$diagnoses$value[1],
    stan_data$stan_data$HIV[1] + stan_data$stan_data$AIDS[1]
  )
  expect_equal(
    sim_estimates$diagnoses_age$value[1],
    sum(
      stan_data$stan_data$HIV_age[1, 1:10] +
        stan_data$stan_data$AIDS_age[1, 1:10]
    )
  )
  expect_equal(
    sim_estimates$undiag_prev_age$value[1],
    sum(stan_data$sim_data$epi_uk$lat_age[1:10, 1, ]) +
      sum(stan_data$sim_data$epi_uk_abroad$lat_age[1:10, 1, ])
  )
})

test_that("mse_list and plots are generated correctly", {
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
  # Test extraction and plotting functions
  mse_list <- make_mse_list(hiv_list)
  expect_true(mse_list$incidence_bias[1, 8] == -1546)
  p1 <- mse_plot(hiv_list, quantity = "incidence")
  expect_s3_class(p1, "ggplot")
  p2 <- bias_plot(hiv_list, quantity = "incidence")
  expect_s3_class(p2, "ggplot")
})

test_that("mse_list and plots are generated correctly for migration model", {
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
  # Test extraction and plotting functions
  p1 <- mse_plot(hiv_list, quantity = "migration")
  expect_s3_class(p1, "ggplot")
  p2 <- bias_plot(hiv_list, quantity = "migration")
  expect_s3_class(p2, "ggplot")
})

test_that("mse_list and plots work for lists", {
  # Simulate diagnoses
  stan_data <- simulate_diagnoses(nquar = 60, sim_seed = 123, rita = TRUE)
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
  mse_list <- make_mse_list(hiv_list)
  expect_true(mse_list$incidence_bias[1, 8] == -997)
  # Test extraction and plotting functions
  p1 <- mse_plot(hiv_list, quantity = "incidence")
  expect_s3_class(p1, "ggplot")
  p2 <- bias_plot(hiv_list, quantity = "incidence")
  expect_s3_class(p2, "ggplot")
})
