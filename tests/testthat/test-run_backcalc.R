test_that("cd4-only back-calculation produces a stanfit object", {
  # Simulate diagnoses
  stan_data <- simulate_diagnoses(sim_seed = 123)
  # Test run_backcalc
  suppressWarnings({
    hiv_list <- run_backcalc(
      stan_data,
      model_seed = 123,
      inf_model = 1, # 1 spline, 2 random walk, 3 gp
      diag_model = 1, # 1 spline, 2 random walk
      iter = 10,
      chains = 1,
      refresh = 0,
      messages = FALSE,
      nworkers = 1
    )
  })
  expect_s4_class(hiv_list$fit, "stanfit")
  expect_true("diag_prob" %in% names(hiv_list))
})

test_that("migration back-calculation produces a stanfit object", {
  stan_data <- simulate_diagnoses(sim_seed = 123, migration = TRUE)
  # Test run_backcalc
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
  expect_s4_class(hiv_list$fit, "stanfit")
  expect_true("diag_prob" %in% names(hiv_list))
})

test_that("rita back-calculation produces a stanfit object", {
  stan_data <- simulate_diagnoses(sim_seed = 123, rita = TRUE)
  # Test run_backcalc
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
  expect_s4_class(hiv_list$fit, "stanfit")
  expect_true("diag_prob" %in% names(hiv_list))
})

test_that("age-specific back-calculation runs without error", {
  stan_data <- simulate_diagnoses(sim_seed = 123, nquar = 40, age = TRUE)
  # Test run_backcalc
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
  expect_s4_class(hiv_list$fit, "stanfit")
  expect_true("d1.q.y" %in% names(hiv_list))
})

test_that("age-specific rita back-calculation runs without error", {
  stan_data <- simulate_diagnoses(
    sim_seed = 123,
    nquar = 40,
    age = TRUE,
    rita = TRUE
  )
  # Test run_backcalc
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
  expect_s4_class(hiv_list$fit, "stanfit")
  expect_true("d1.q.y" %in% names(hiv_list))
})
