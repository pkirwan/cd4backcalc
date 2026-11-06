test_that("rita completeness by state simulation produces expected HIV data", {
  state_completeness <- c(0.8, 0.8, 0.7, 0.6, 0.5)
  complete_diags <- simulate_diagnoses(
    sim_seed = 123,
    age = TRUE,
    rita = TRUE,
    dist = FALSE
  )
  sim_diags <- simulate_diagnoses(
    sim_seed = 123,
    age = TRUE,
    rita = TRUE,
    dist = FALSE,
    rita_completeness_by_state = state_completeness
  )
  expected_cd4 <- round(sweep(
    complete_diags$stan_data$CD4,
    2,
    state_completeness,
    "*"
  ))

  expect_equal(sim_diags$stan_data$CD4, expected_cd4)
  expect_equal(sim_diags$stan_data$HIV, complete_diags$stan_data$HIV)
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
    mu <- generate_mu(20, 4, pattern)
    conditional_cd4 <-
      mu[, 1:4, drop = FALSE] / rowSums(mu[, 1:4, drop = FALSE])
    expect_equal(
      unname(conditional_cd4),
      matrix(expected_cd4[[pattern]], 20, 4, byrow = TRUE)
    )
    expect_equal(mu[, 5], expected_transfer[[pattern]])
    expect_equal(rowSums(mu), rep(1, 20))
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

test_that("recent-arrival CD4 data follow migration model dimensions", {
  for (age in c(FALSE, TRUE)) {
    for (rita in c(FALSE, TRUE)) {
      sim_diags <- simulate_diagnoses(
        nquar = 40,
        nage = 52,
        sim_seed = 123,
        age = age,
        rita = rita,
        migration = TRUE,
        m_pattern = "varying",
        recent_arrival_completeness = 0.2
      )

      expected_dim <- if (age) c(10L, 4L, 4L) else c(10L, 4L)
      expect_equal(sim_diags$stan_data$nyear, 10L)
      expect_equal(dim(sim_diags$stan_data$recent_arrival_cd4), expected_dim)
      expect_equal(
        dim(sim_diags$stan_data$recent_arrival_cd4_total),
        expected_dim
      )
      expect_true(sum(sim_diags$stan_data$recent_arrival_cd4) > 0)
      expect_true(all(sim_diags$stan_data$recent_arrival_cd4 >= 0))
      expect_true(all(
        sim_diags$stan_data$recent_arrival_cd4 <=
          sim_diags$stan_data$recent_arrival_cd4_total
      ))
    }
  }
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
  sim_diags$stan_data$recent_arrival_cd4_total <- NULL

  stan_data <- model_type_prep(sim_diags)

  expect_equal(stan_data$nyear, 5L)
  expect_equal(stan_data$recent_arrival_cd4, array(0L, c(5L, 4L, 4L)))
  expect_equal(
    stan_data$recent_arrival_cd4_total,
    array(0L, c(5L, 4L, 4L))
  )

  sim_diags$stan_data$sim_config$age <- FALSE
  sim_diags$stan_data$recent_arrival_cd4 <- NULL
  sim_diags$stan_data$recent_arrival_cd4_total <- NULL
  stan_data <- model_type_prep(sim_diags)
  expect_equal(stan_data$recent_arrival_cd4, matrix(0L, 5L, 4L))
  expect_equal(stan_data$recent_arrival_cd4_total, matrix(0L, 5L, 4L))
})

test_that("recent-arrival CD4 data require their conditional denominator", {
  sim_diags <- simulate_diagnoses(
    nquar = 40,
    sim_seed = 123,
    migration = TRUE,
    recent_arrival_completeness = 0.5
  )
  sim_diags$stan_data$recent_arrival_cd4_total <- NULL

  expect_error(
    model_type_prep(sim_diags),
    "recent_arrival_cd4_total is required"
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
