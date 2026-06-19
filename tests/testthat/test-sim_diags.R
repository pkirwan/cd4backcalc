test_that("migration simulation produces expected HIV data", {
  sim_diags <- simulate_diagnoses(sim_seed = 123, migration = TRUE)
  expect_equal(
    sim_diags$sim_data$omega[1],
    sum(sim_diags$sim_data$omega_age[1, ])
  )
})

test_that("age-specific migration simulation produces expected HIV data", {
  sim_diags <- simulate_diagnoses(sim_seed = 123, age = TRUE, migration = TRUE)
  expect_equal(
    sim_diags$stan_data$ABROAD,
    rowSums(sim_diags$stan_data$ABROAD_age)
  )
})

test_that("sparse detection observations follow requested spacing", {
  sim_diags <- simulate_diagnoses(
    sim_seed = 123,
    migration = TRUE,
    detect_obs_sd = 0.3,
    detect_obs_times = seq(1, 120, by = 20)
  )
  expect_equal(
    sim_diags$stan_data$detect_obs_index,
    c(1L, 21L, 41L, 61L, 81L, 101L)
  )
  expect_length(sim_diags$stan_data$detect_logit_obs, 6)
  expect_null(sim_diags$stan_data$n_detect_obs)
})

test_that("age-specific migration and rita simulation produces expected HIV data", {
  sim_diags <- simulate_diagnoses(
    sim_seed = 123,
    age = TRUE,
    migration = TRUE,
    rita = TRUE
  )
  expect_equal(
    sim_diags$stan_data$ABROAD,
    rowSums(sim_diags$stan_data$ABROAD_age)
  )
})

test_that("age-specific migration simulation produces expected AIDS data", {
  sim_diags <- simulate_diagnoses(
    sim_seed = 123,
    age = TRUE,
    migration = TRUE,
    sim_type = "combo_3"
  )
  expect_equal(
    sim_diags$stan_data$AIDS[, 1],
    rowSums(sim_diags$stan_data$AIDS_age[,, 1])
  )
  expect_equal(
    sim_diags$stan_data$AIDS[, 2],
    rowSums(sim_diags$stan_data$AIDS_age[,, 2])
  )
})

test_that("age-specific migration simulation with gamma produces expected external migration data", {
  sim_diags <- simulate_diagnoses(
    sim_seed = 123,
    age = TRUE,
    migration = TRUE,
    sim_type = "combo_3",
    gamma_obs_sd = 0.1
  )
  expect_equal(
    sim_diags$stan_data$gamma_obs_sd,
    0.1
  )
  expect_equal(
    length(sim_diags$stan_data$gamma_obs_mean),
    sim_diags$stan_data$nquar
  )
})
