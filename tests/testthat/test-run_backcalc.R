# for local testing use Sys.setenv(CD4BACKCALC_CMDSTAN_DIR = cmdstanr::cmdstan_default_install_path()) to avoid re-compiling the model each time

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
      iter_warmup = 5,
      iter_sampling = 5,
      chains = 1,
      refresh = 0,
      messages = FALSE,
      nworkers = 1
    )
  })
  expect_r6_class(hiv_list$fit, "CmdStanMCMC")
  expect_true("diag_prob" %in% names(hiv_list))
})

test_that("run_backcalc accepts user-supplied inits", {
  stan_data <- simulate_diagnoses(sim_seed = 123)

  suppressWarnings({
    initial_fit <- run_backcalc(
      stan_data,
      model_seed = 123,
      iter_warmup = 5,
      iter_sampling = 5,
      chains = 1,
      refresh = 0,
      messages = FALSE,
      nworkers = 1
    )
  })

  inits <- create_stan_inits(initial_fit$fit)

  suppressWarnings({
    hiv_list <- run_backcalc(
      stan_data,
      model_seed = 123,
      iter_warmup = 5,
      iter_sampling = 5,
      chains = 1,
      refresh = 0,
      messages = FALSE,
      nworkers = 1,
      init = inits
    )
  })

  expect_r6_class(hiv_list$fit, "CmdStanMCMC")
  expect_true("diag_prob" %in% names(hiv_list))
})

test_that("back-calculation can cope with missing data", {
  # Simulate diagnoses
  stan_data <- simulate_diagnoses(sim_seed = 123, cd4_completeness = 0.5)
  # Test run_backcalc
  suppressWarnings({
    hiv_list <- run_backcalc(
      stan_data,
      model_seed = 123,
      inf_model = 1, # 1 spline, 2 random walk, 3 gp
      diag_model = 1, # 1 spline, 2 random walk
      iter_warmup = 5,
      iter_sampling = 5,
      chains = 1,
      refresh = 0,
      messages = FALSE,
      nworkers = 1
    )
  })
  expect_r6_class(hiv_list$fit, "CmdStanMCMC")
  expect_true("diag_prob" %in% names(hiv_list))

  # Simulate diagnoses
  stan_data <- simulate_diagnoses(sim_seed = 123, rita = TRUE, rita_completeness = 0.5)
  # Test run_backcalc
  suppressWarnings({
    hiv_list <- run_backcalc(
      stan_data,
      model_seed = 123,
      inf_model = 1, # 1 spline, 2 random walk, 3 gp
      diag_model = 1, # 1 spline, 2 random walk
      iter_warmup = 5,
      iter_sampling = 5,
      chains = 1,
      refresh = 0,
      messages = FALSE,
      nworkers = 1
    )
  })
  expect_r6_class(hiv_list$fit, "CmdStanMCMC")
  expect_true("diag_prob" %in% names(hiv_list))
})

test_that("migration back-calculation produces a stanfit object", {
  stan_data <- simulate_diagnoses(sim_seed = 123, migration = TRUE)
  # Test run_backcalc
  suppressWarnings({
    hiv_list <- run_backcalc(
      stan_data,
      model_seed = 123,
      iter_warmup = 5,
      iter_sampling = 5,
      chains = 1,
      refresh = 0,
      messages = FALSE,
      nworkers = 1
    )
  })
  expect_r6_class(hiv_list$fit, "CmdStanMCMC")
  expect_true("diag_prob" %in% names(hiv_list))
})

test_that("rita back-calculation produces a stanfit object", {
  stan_data <- simulate_diagnoses(sim_seed = 123, rita = TRUE)
  # Test run_backcalc
  suppressWarnings({
    hiv_list <- run_backcalc(
      stan_data,
      model_seed = 123,
      iter_warmup = 5,
      iter_sampling = 5,
      chains = 1,
      refresh = 0,
      messages = FALSE,
      nworkers = 1
    )
  })
  expect_r6_class(hiv_list$fit, "CmdStanMCMC")
  expect_true("diag_prob" %in% names(hiv_list))
})

test_that("rita data with cd4-only back-calculation produces a stanfit object", {
  stan_data <- simulate_diagnoses(sim_seed = 123, rita = TRUE)
  stan_data <- sim_remodel(stan_data, rita = FALSE)
  # Test run_backcalc
  suppressWarnings({
    hiv_list <- run_backcalc(
      stan_data,
      rita = FALSE,
      model_seed = 123,
      iter_warmup = 5,
      iter_sampling = 5,
      chains = 1,
      refresh = 0,
      messages = FALSE,
      nworkers = 1
    )
  })
  expect_r6_class(hiv_list$fit, "CmdStanMCMC")
  expect_equal(nrow(hiv_list$diag_prob), 120 * 4 * 4) # 120 quarters, 4 summaries (including mean), 4 cd4 states
})


test_that("age-specific back-calculation runs without error", {
  stan_data <- simulate_diagnoses(sim_seed = 123, nquar = 40, age = TRUE)
  # Test run_backcalc
  suppressWarnings({
    hiv_list <- run_backcalc(
      stan_data,
      model_seed = 123,
      iter_warmup = 5,
      iter_sampling = 5,
      chains = 1,
      refresh = 0,
      messages = FALSE,
      nworkers = 1
    )
  })
  expect_r6_class(hiv_list$fit, "CmdStanMCMC")
  expect_true("diag_prob_age" %in% names(hiv_list))
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
      iter_warmup = 5,
      iter_sampling = 5,
      chains = 1,
      refresh = 0,
      messages = FALSE,
      nworkers = 1
    )
  })
  expect_r6_class(hiv_list$fit, "CmdStanMCMC")
  expect_true("diag_prob_age" %in% names(hiv_list))
})

test_that("age-specific migration back-calculation runs without error", {
  stan_data <- simulate_diagnoses(
    sim_seed = 123,
    nquar = 40,
    age = TRUE,
    migration = TRUE,
    sim_type = "combo_3",
    eta_sd = 0.3,
    detect_obs_sd = 0.3,
    detect_obs_times = seq(1, 40, by = 20)
  )
  # Test run_backcalc
  suppressWarnings({
    hiv_list <- run_backcalc(
      stan_data,
      model_seed = 123,
      iter_warmup = 5,
      iter_sampling = 5,
      chains = 1,
      refresh = 0,
      messages = FALSE,
      nworkers = 1
    )
  })
  expect_r6_class(hiv_list$fit, "CmdStanMCMC")
  expect_true("diag_prob_age" %in% names(hiv_list))
})

test_that("age-specific migration and rita back-calculation runs without error", {
  stan_data <- simulate_diagnoses(
    sim_seed = 123,
    nquar = 40,
    age = TRUE,
    migration = TRUE,
    rita = TRUE,
    sim_type = "combo_3"
  )
  # Test run_backcalc
  suppressWarnings({
    hiv_list <- run_backcalc(
      stan_data,
      model_seed = 123,
      iter_warmup = 5,
      iter_sampling = 5,
      chains = 1,
      refresh = 0,
      messages = FALSE,
      nworkers = 1
    )
  })
  expect_r6_class(hiv_list$fit, "CmdStanMCMC")
  expect_true("diag_prob_age" %in% names(hiv_list))
})
