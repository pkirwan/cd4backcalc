test_that("cd4-only post-processing works", {
  # Simulate diagnoses
  sim_diags <- simulate_diagnoses(
    nquar = 60,
    sim_seed = 123,
    cd4_completeness = 0.8
  )
  stan_data <- model_type_prep(sim_diags, model_seed = 123)
  # Fit back-calculation model
  suppressWarnings({
    fit <- fit_stan_model(
      stan_data,
      iter_warmup = 10,
      iter_sampling = 10,
      chains = 2,
      refresh = 0,
      messages = FALSE
    )
  })
  expect_r6_class(fit, "CmdStanMCMC")

  # Post-process fitted model to obtain posterior predictive
  hiv_list <- postproc(
    fit,
    stan_data,
    start_yr = 0,
    nworkers = 1
  )
  expect_true("diag_prob" %in% names(hiv_list))

  # Chained post-processing
  hiv_list <- postproc(
    fit,
    stan_data,
    start_yr = 0,
    nworkers = 1,
    chain_results = TRUE
  )
  expect_true("diag_prob_1" %in% names(hiv_list))
})

test_that("migration post-processing works", {
  # Simulate diagnoses
  sim_diags <- simulate_diagnoses(
    nquar = 60,
    migration = TRUE,
    sim_seed = 123,
    cd4_completeness = 0.8
  )
  stan_data <- model_type_prep(sim_diags, model_seed = 123)
  # Fit back-calculation model
  suppressWarnings({
    fit <- fit_stan_model(
      stan_data,
      iter_warmup = 10,
      iter_sampling = 10,
      chains = 1,
      refresh = 0,
      messages = FALSE
    )
  })
  expect_r6_class(fit, "CmdStanMCMC")

  # Post-process fitted model to obtain posterior predictive
  hiv_list <- postproc(
    fit,
    stan_data,
    start_yr = 0,
    nworkers = 1
  )
  expect_true("migration_prob" %in% names(hiv_list))
  expect_true("migrant_hiv_prevalence" %in% names(hiv_list))
})

test_that("undiag_prev uses direct latent prevalence", {
  stan_data <- list(
    nquar = 4,
    model_config = list(migration = FALSE)
  )

  prev_dist <- sapply(
    list(
      list(
        undiagnosed_total = c(10, 9, 8, 7),
        undiagnosed = matrix(1000, nrow = 1)
      ),
      list(
        undiagnosed_total = c(14, 13, 12, 11),
        undiagnosed = matrix(2000, nrow = 1)
      )
    ),
    identity
  )

  inf_dist <- sapply(
    list(
      list(
        AIDS.Dx.Infect.Dist = matrix(1000, nrow = 4, ncol = 4),
        HIV.Dx.Infect.Dist = array(1000, dim = c(1, 4, 4))
      ),
      list(
        AIDS.Dx.Infect.Dist = matrix(2000, nrow = 4, ncol = 4),
        HIV.Dx.Infect.Dist = array(2000, dim = c(1, 4, 4))
      )
    ),
    identity
  )

  undiag_prev <- undiag.prev.df(
    stan_data,
    prev_dist,
    inf_dist,
    start.yr = 0
  )

  median_values <- undiag_prev |>
    dplyr::filter(.data$quartile == "50%") |>
    dplyr::pull(.data$value)

  expect_equal(median_values, c(12, 11, 10, 9))
})

test_that("RITA forward calculation transitions initial prevalence at quarter 1", {
  init_prev <- c(10, 20, 30, 40, 50, 60, 70)
  h <- c(100, 100)
  d1 <- c(0.1, 0.1)
  d2 <- c(0.2, 0.2)
  d3 <- c(0.3, 0.3)
  d4 <- c(0.4, 0.4)
  d5 <- c(0.5, 0.5)
  q <- c(0.1, 0.2, 0.3, 0.4)

  out <- forward_calculate_rita(
    init_prev,
    h,
    d1,
    d2,
    d3,
    d4,
    d5,
    q
  )

  expected <-
    sum(init_prev[1:3]) * (1 - d1[1]) +
    init_prev[4] * (1 - d2[1]) +
    init_prev[5] * (1 - d3[1]) +
    init_prev[6] * (1 - d4[1]) +
    init_prev[7] * (1 - d5[1]) * (1 - q[4]) +
    h[1]

  expect_equal(sum(out$prevalence[1, 1:7, 1]), expected)
})

test_that("rita post-processing works", {
  # Simulate diagnoses
  sim_diags <- simulate_diagnoses(
    nquar = 60,
    rita = TRUE,
    sim_seed = 123,
    cd4_completeness = 0.8
  )
  stan_data <- model_type_prep(sim_diags, model_seed = 123)
  # Fit back-calculation model
  suppressWarnings({
    fit <- fit_stan_model(
      stan_data,
      iter_warmup = 10,
      iter_sampling = 10,
      chains = 1,
      refresh = 0,
      messages = FALSE
    )
  })
  expect_r6_class(fit, "CmdStanMCMC")

  # Post-process fitted model to obtain posterior predictive
  hiv_list <- postproc(
    fit,
    stan_data,
    start_yr = 0,
    nworkers = 1
  )
  expect_true("diag_prob" %in% names(hiv_list))
})

test_that("migration and rita post-processing works", {
  # Simulate diagnoses
  sim_diags <- simulate_diagnoses(
    nquar = 60,
    migration = TRUE,
    rita = TRUE,
    sim_seed = 123,
    cd4_completeness = 0.8
  )
  stan_data <- model_type_prep(sim_diags, model_seed = 123)
  # Fit back-calculation model
  suppressWarnings({
    fit <- fit_stan_model(
      stan_data,
      iter_warmup = 10,
      iter_sampling = 10,
      chains = 1,
      refresh = 0,
      messages = FALSE
    )
  })
  expect_r6_class(fit, "CmdStanMCMC")

  # Post-process fitted model to obtain posterior predictive
  hiv_list <- postproc(
    fit,
    stan_data,
    start_yr = 0,
    nworkers = 1
  )
  expect_true("migration_prob" %in% names(hiv_list))
})

test_that("age-specific post-processing works", {
  sim_diags <- simulate_diagnoses(
    sim_seed = 123,
    nquar = 40,
    age = TRUE,
    cd4_completeness = 0.8
  )
  stan_data <- model_type_prep(sim_diags, model_seed = 123)
  # Fit back-calculation model
  suppressWarnings({
    fit <- fit_stan_model(
      stan_data,
      iter_warmup = 5,
      iter_sampling = 5,
      chains = 1,
      refresh = 0,
      messages = FALSE
    )
  })
  expect_r6_class(fit, "CmdStanMCMC")

  # Post-process fitted model to obtain posterior predictive
  hiv_list <- postproc(
    fit,
    stan_data,
    start_yr = 0,
    nworkers = 1
  )
  expect_true("incidence" %in% names(hiv_list))
  expect_true("undiag_prev" %in% names(hiv_list))
  expect_true("incidence_age" %in% names(hiv_list))
  expect_true("diag_prob_age" %in% names(hiv_list))
  expect_true("diagnoses_age" %in% names(hiv_list))
})

test_that("age-specific rita post-processing works", {
  sim_diags <- simulate_diagnoses(
    sim_seed = 123,
    nquar = 40,
    age = TRUE,
    rita = TRUE,
    cd4_completeness = 0.8
  )
  stan_data <- model_type_prep(sim_diags, model_seed = 123)
  # Fit back-calculation model
  suppressWarnings({
    fit <- fit_stan_model(
      stan_data,
      iter_warmup = 5,
      iter_sampling = 5,
      chains = 1,
      refresh = 0,
      messages = FALSE
    )
  })
  expect_r6_class(fit, "CmdStanMCMC")

  # Post-process fitted model to obtain posterior predictive
  hiv_list <- postproc(
    fit,
    stan_data,
    start_yr = 0,
    nworkers = 1
  )
  expect_true("incidence" %in% names(hiv_list))
  expect_true("undiag_prev" %in% names(hiv_list))
})

test_that("age-specific migration model fit works", {
  sim_diags <- simulate_diagnoses(
    sim_seed = 123,
    nquar = 40,
    age = TRUE,
    migration = TRUE,
    cd4_completeness = 0.8
  )
  stan_data <- model_type_prep(sim_diags, model_seed = 123)
  # Fit back-calculation model
  suppressWarnings({
    fit <- fit_stan_model(
      stan_data,
      iter_warmup = 5,
      iter_sampling = 5,
      chains = 1,
      refresh = 0,
      messages = FALSE
    )
  })
  expect_r6_class(fit, "CmdStanMCMC")

  # Post-process fitted model to obtain posterior predictive
  hiv_list <- postproc(
    fit,
    stan_data,
    start_yr = 0,
    nworkers = 1
  )
  expect_true("incidence" %in% names(hiv_list))
  expect_true("undiag_prev" %in% names(hiv_list))
  expect_true("migration_prob" %in% names(hiv_list))
})

test_that("age-specific rita migration model fit works", {
  sim_diags <- simulate_diagnoses(
    sim_seed = 123,
    nquar = 40,
    age = TRUE,
    migration = TRUE,
    rita = TRUE,
    cd4_completeness = 0.8
  )
  stan_data <- model_type_prep(sim_diags, model_seed = 123)
  # Fit back-calculation model
  suppressWarnings({
    fit <- fit_stan_model(
      stan_data,
      iter_warmup = 5,
      iter_sampling = 5,
      chains = 1,
      refresh = 0,
      messages = FALSE
    )
  })
  expect_r6_class(fit, "CmdStanMCMC")

  # Post-process fitted model to obtain posterior predictive
  hiv_list <- postproc(
    fit,
    stan_data,
    start_yr = 0,
    nworkers = 1
  )
  expect_true("incidence" %in% names(hiv_list))
  expect_true("undiag_prev" %in% names(hiv_list))
  expect_true("migration_prob" %in% names(hiv_list))
})
