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
  p1 <- mse_plot(list(hiv_list), quantity = "incidence")
  expect_s3_class(p1, "ggplot")
  p2 <- bias_plot(list(hiv_list), quantity = "incidence")
  expect_s3_class(p2, "ggplot")
})

test_that("make_mse_list uses posterior means and fit-specific truths", {
  make_mock_fit <- function(truth, estimate_mean, estimate_median) {
    zero_matrix <- matrix(0, nrow = 3, ncol = 1)
    zero_migration <- matrix(0, nrow = 3, ncol = 2)
    list(
      incidence = tibble::tibble(
        year = c(0, 0.25, 0.5),
        `50%` = estimate_median,
        mean = estimate_mean
      ),
      stan_data = list(
        nquar = 3,
        nage = 1,
        HIV = c(0, 0, 0),
        AIDS = c(0, 0, 0),
        model_config = list(migration = FALSE),
        sim_config = list(
          migration = FALSE,
          rita = FALSE,
          nlatst = 2
        )
      ),
      sim_data = list(
        epi_uk = list(lat = zero_matrix, diags = zero_matrix),
        epi_uk_abroad = list(lat = zero_matrix, diags = zero_matrix),
        epi_abroad = list(diags = zero_matrix),
        infs = cbind(truth),
        m = zero_migration,
        omega = c(0, 0, 0),
        d_probs = array(0, dim = c(3, 1, 4)),
        misclassified = c(0, 0, 0)
      )
    )
  }

  fit_1 <- make_mock_fit(
    truth = c(0, 8, 18),
    estimate_mean = c(0, 10, 20),
    estimate_median = c(0, 100, 200)
  )
  fit_2 <- make_mock_fit(
    truth = c(0, 11, 21),
    estimate_mean = c(0, 13, 23),
    estimate_median = c(0, 130, 230)
  )

  mse_list <- make_mse_list(list(fit_1, fit_2))

  expect_equal(mse_list$incidence_bias$year, c(0.25, 0.5))
  expect_equal(mse_list$incidence_bias$bias, c(2, 2))
  expect_equal(mse_list$incidence_mse$mse, c(4, 4))
  expect_equal(mse_list$incidence_mse$n, c(2L, 2L))

  filtered_mse <- make_mse_list(list(fit_1, fit_2), years = 0.5)
  expect_equal(filtered_mse$incidence_bias$year, 0.5)
  expect_equal(filtered_mse$incidence_bias$bias, 2)
  expect_equal(filtered_mse$incidence_mse$mse, c(4, 4))
  expect_equal(filtered_mse$incidence_mse$n, c(1L, 1L))
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
  p1 <- mse_plot(list("Fit 1" = hiv_list, "Fit 2" = hiv_list_2), quantity = "incidence")
  expect_s3_class(p1, "ggplot")
  p2 <- bias_plot(list("Fit 1" = hiv_list, "Fit 2" = hiv_list_2), quantity = "incidence")
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
  p1 <- mse_plot(list(hiv_list), quantity = "undiag_migration")
  expect_s3_class(p1, "ggplot")
  p2 <- bias_plot(list(hiv_list), quantity = "undiag_migration")
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

  # Test extraction and plotting functions
  p1 <- mse_plot(list(hiv_list), quantity = "incidence")
  expect_s3_class(p1, "ggplot")
  p2 <- bias_plot(list(hiv_list), quantity = "incidence")
  expect_s3_class(p2, "ggplot")
})
