test_that("mse_list and plots are generated correctly", {
  # Simulate diagnoses
  stan_data <- simulate_diagnoses(nquar = 60, sim_seed = 123)
  # Fit back-calculation model
  suppressWarnings({
    hiv_list <- run_backcalc(
      stan_data,
      model_seed = 123,
      iter_warmup = 10,
      iter_sampling = 10,
      chains = 1,
      refresh = 0,
      messages = FALSE,
      nworkers = 1
    )
  })
  # Test extraction and plotting functions
  mse_list <- make_mse_list(hiv_list)
  # expect_true(mse_list$incidence_bias[1, 8] == -1146)
  p1 <- mse_plot(hiv_list, quantity = "incidence")
  expect_s3_class(p1, "ggplot")
  p2 <- bias_plot(hiv_list, quantity = "incidence")
  expect_s3_class(p2, "ggplot")
})

test_that("mse_list and plots are generated correctly for dual fits", {
  # Simulate diagnoses
  stan_data <- simulate_diagnoses(nquar = 60, sim_seed = 123)
  # Fit back-calculation model
  suppressWarnings({
    hiv_list <- run_backcalc(
      stan_data,
      model_seed = 123,
      iter_warmup = 10,
      iter_sampling = 10,
      chains = 1,
      refresh = 0,
      messages = FALSE,
      nworkers = 1
    )
  })
  # Fit back-calculation model
  suppressWarnings({
    hiv_list_2 <- run_backcalc(
      stan_data,
      model_seed = 1234,
      iter_warmup = 10,
      iter_sampling = 10,
      chains = 1,
      refresh = 0,
      messages = FALSE,
      nworkers = 1
    )
  })
  # Test extraction and plotting functions
  mse_list <- make_mse_list(hiv_list)
  # expect_true(mse_list$incidence_bias[1, 8] == -1146)
  p1 <- mse_plot(hiv_list, quantity = "incidence", hiv_list_2 = hiv_list_2, model_1_name = "Fit 1", model_2_name = "Fit 2")
  expect_s3_class(p1, "ggplot")
  p2 <- bias_plot(hiv_list, quantity = "incidence", hiv_list_2 = hiv_list_2, model_1_name = "Fit 1", model_2_name = "Fit 2")
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
      iter_warmup = 10,
      iter_sampling = 10,
      chains = 1,
      refresh = 0,
      messages = FALSE,
      nworkers = 1
    )
  })
  # Test extraction and plotting functions
  p1 <- mse_plot(hiv_list, quantity = "undiag_migration")
  expect_s3_class(p1, "ggplot")
  p2 <- bias_plot(hiv_list, quantity = "undiag_migration")
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
        iter_warmup = 10,
      iter_sampling = 10,
        chains = 1,
        refresh = 0,
        messages = FALSE,
        nworkers = 1
      )
    })
  }
  mse_list <- make_mse_list(hiv_list)
  # don't check specific value on Mac and Linux CI
  # expect_true(mse_list$incidence_bias[1, 8] == -1108)
  
  # Test extraction and plotting functions
  p1 <- mse_plot(hiv_list, quantity = "incidence")
  expect_s3_class(p1, "ggplot")
  p2 <- bias_plot(hiv_list, quantity = "incidence")
  expect_s3_class(p2, "ggplot")
})
