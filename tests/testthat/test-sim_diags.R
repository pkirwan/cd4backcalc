test_that("rita completeness by state simulation produces expected HIV data", {
  sim_diags <- simulate_diagnoses(
    sim_seed = 123,
    age = TRUE,
    rita = TRUE,
    rita_completeness_by_state = c(0.8, 0.8, 0.7, 0.6, 0.5)
  )
  expect_true(all(
    rowSums(sim_diags$stan_data$CD4) < 0.8 * sim_diags$stan_data$HIV
  ))
})

test_that("migration simulation produces expected HIV data", {
  sim_diags <- simulate_diagnoses(sim_seed = 123, migration = TRUE)
  expect_equal(
    sim_diags$sim_data$omega[1],
    sum(sim_diags$sim_data$omega_age[1, ])
  )
})

test_that("migration simulation preserves scientific input contracts", {
  sim_diags <- simulate_diagnoses(
    nquar = 20,
    sim_seed = 123,
    migration = TRUE,
    dist = FALSE,
    eta_sd = 0
  )

  expect_equal(sim_diags$stan_data$eta_sd, 0)
  expect_equal(
    exp(sim_diags$stan_data$eta_mean),
    sim_diags$sim_data$ratio_abroad_uk
  )
  expect_equal(rowSums(sim_diags$sim_data$mu), rep(1, 20))
  expect_equal(rowSums(sim_diags$sim_data$m), sim_diags$sim_data$omega)
})

test_that("migration probability patterns follow the Stan decomposition", {
  expected_cd4 <- list(
    constant = rep(0.25, 4),
    differing = c(1, 2, 2, 1) / 6,
    varying = c(2, 2, 2, 1) / 7
  )
  expected_transfer <- list(
    constant = rep(0.2, 20),
    differing = rep(0.4, 20),
    varying = seq(0.2, 0.4, length.out = 20)
  )

  for (pattern in names(expected_cd4)) {
    mu <- generate_mu(20, 4, pattern, dist = FALSE)
    conditional_cd4 <-
      mu[, 1:4, drop = FALSE] / rowSums(mu[, 1:4, drop = FALSE])
    expect_equal(
      unname(conditional_cd4),
      matrix(expected_cd4[[pattern]], 20, 4, byrow = TRUE)
    )
    expect_equal(mu[, 5], expected_transfer[[pattern]])

    set.seed(123)
    stochastic_mu <- generate_mu(20, 4, pattern, dist = TRUE)
    stochastic_conditional_cd4 <-
      stochastic_mu[, 1:4, drop = FALSE] /
        rowSums(stochastic_mu[, 1:4, drop = FALSE])
    expect_equal(
      unname(stochastic_conditional_cd4),
      matrix(expected_cd4[[pattern]], 20, 4, byrow = TRUE)
    )
    expect_equal(rowSums(stochastic_mu), rep(1, 20))
    expect_false(isTRUE(all.equal(
      stochastic_mu[, 5],
      expected_transfer[[pattern]]
    )))
  }
})

test_that("transfer detection preserves totals and age allocation", {
  sim_diags <- simulate_diagnoses(
    nquar = 40,
    nage = 52,
    sim_seed = 123,
    age = TRUE,
    migration = TRUE,
    m_pattern = "varying"
  )

  expect_equal(
    sim_diags$stan_data$ABROAD + sim_diags$sim_data$misclassified,
    sim_diags$sim_data$m[, 5]
  )
  expect_equal(
    sim_diags$stan_data$ABROAD,
    rowSums(sim_diags$stan_data$ABROAD_age)
  )
})

test_that("recent-arrival CD4 data are annual age-by-state counts", {
  sim_diags <- simulate_diagnoses(
    nquar = 40,
    nage = 52,
    sim_seed = 123,
    age = TRUE,
    migration = TRUE,
    m_pattern = "varying",
    recent_arrival_completeness = 0.2
  )

  expect_equal(sim_diags$stan_data$nyear, 10L)
  expect_equal(dim(sim_diags$stan_data$recent_arrival_cd4), c(10L, 4L, 4L))
  expect_true(sum(sim_diags$stan_data$recent_arrival_cd4) > 0)
  expect_true(all(sim_diags$stan_data$recent_arrival_cd4 >= 0))
})

test_that("migration model preparation defaults missing recent-arrival data", {
  sim_diags <- simulate_diagnoses(
    nquar = 20,
    sim_seed = 123,
    age = TRUE,
    migration = TRUE,
    m_pattern = "varying"
  )
  sim_diags$stan_data$nyear <- NULL
  sim_diags$stan_data$recent_arrival_cd4 <- NULL

  stan_data <- model_type_prep(sim_diags)

  expect_equal(stan_data$nyear, 5L)
  expect_equal(stan_data$recent_arrival_cd4, array(0L, c(5L, 4L, 4L)))
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
    gamma_obs_sd = 0.1,
    gamma_obs_scale_trend = 0.1,
    gamma_obs_drift_sd = 0.05,
    abroad_age_trend = 0.1,
    detect_step_quarter = 30,
    detect_step_size = 0.1
  )
  expect_equal(
    sim_diags$stan_data$gamma_obs_sd,
    0.1
  )
  expect_equal(
    length(sim_diags$stan_data$gamma_obs),
    sim_diags$stan_data$nquar
  )
  expect_equal(sim_diags$stan_data$gamma_obs_drift_sd, 0.05)
  expect_gt(
    stats::sd(diff(log(sim_diags$sim_data$migrant_hiv_prevalence))),
    0
  )
  t_std <- (seq_len(sim_diags$stan_data$nquar) -
    0.5 * (sim_diags$stan_data$nquar + 1)) / sim_diags$stan_data$nquar
  expected_log_prevalence <-
    -log(sim_diags$stan_data$sim_config$gamma_obs_scale) - 0.1 * t_std
  expect_equal(
    mean(
      log(sim_diags$sim_data$migrant_hiv_prevalence) -
        expected_log_prevalence
    ),
    0,
    tolerance = 1e-12
  )
})

test_that("gamma_obs_grouped vector controls grouped external migration by quarter", {
  gamma_grouped_schedule <- c(rep(FALSE, 20), rep(TRUE, 20), rep(FALSE, 20))
  sim_diags <- simulate_diagnoses(
    sim_seed = 123,
    age = TRUE,
    migration = TRUE,
    nquar = 60,
    sim_type = "combo_3",
    gamma_obs_grouped = gamma_grouped_schedule,
    gamma_obs_sd = 0.1
  )

  expect_equal(length(sim_diags$stan_data$gamma_obs), sim_diags$stan_data$nquar)
  expect_equal(
    dim(sim_diags$stan_data$gamma_obs_group),
    c(sim_diags$stan_data$nquar, 4L)
  )
  expect_false("gamma_obs_grouped" %in% names(sim_diags$stan_data))

  expect_true(all(
    sim_diags$stan_data$gamma_obs_group[!gamma_grouped_schedule, ] == -1
  ))
  expect_true(all(
    sim_diags$stan_data$gamma_obs_group[gamma_grouped_schedule, ] > -1
  ))
})
