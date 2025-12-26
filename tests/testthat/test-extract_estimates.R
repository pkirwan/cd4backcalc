test_that("extract_simulated_quantity correctly extracts values for ai", {
  # Simulate diagnoses
  stan_data <- simulate_diagnoses(nquar = 60, sim_seed = 123)
  sim_estimates <- extract_simulated_quantity(stan_data, quantity = "incidence")
  expect_equal(
    sim_estimates$value[1],
    stan_data$sim_data$infs[1, 1]
  )
  sim_estimates <- extract_simulated_quantity(
    stan_data,
    quantity = "undiag_prev"
  )
  expect_equal(
    sim_estimates$value[1],
    sum(
      stan_data$sim_data$epi_uk$lat[1, ] +
        stan_data$sim_data$epi_uk_abroad$lat[1, ]
    )
  )
  sim_estimates <- extract_simulated_quantity(stan_data, quantity = "diagnoses")
  expect_equal(
    sim_estimates$value[1],
    stan_data$stan_data$HIV[1] + stan_data$stan_data$AIDS[1]
  )
  sim_estimates <- extract_simulated_quantity(stan_data, quantity = "diag_prob")
  expect_equal(
    sim_estimates$value[1],
    stan_data$sim_data$d_probs[1, 1]
  )
  sim_estimates <- extract_simulated_quantity(
    stan_data,
    quantity = "cd4_strata"
  )
  expect_equal(
    sim_estimates$value[1],
    stan_data$stan_data$CD4[1, 1]
  )

  stan_data <- simulate_diagnoses(nquar = 60, sim_seed = 123, migration = TRUE)
  sim_estimates <- extract_simulated_quantity(
    stan_data,
    quantity = "undiag_migration"
  )
  expect_equal(
    sim_estimates$value[1],
    sum(stan_data$sim_data$m[1, 1:4]) -
      (stan_data$stan_data$ABROAD[1] *
        ((1 / stan_data$sim_data$detect_prob[1]) - 1))
  )
  sim_estimates <- extract_simulated_quantity(
    stan_data,
    quantity = "all_migration"
  )
  expect_equal(
    sim_estimates$value[1],
    sum(stan_data$sim_data$omega[1])
  )
  sim_estimates <- extract_simulated_quantity(
    stan_data,
    quantity = "migration_prob"
  )
  expect_equal(
    sim_estimates$value[1],
    stan_data$sim_data$mu[1, 1]
  )
  sim_estimates <- extract_simulated_quantity(
    stan_data,
    quantity = "ratio_abroad_uk"
  )
  expect_equal(
    sim_estimates$value[1],
    stan_data$sim_data$ratio_abroad_uk[1]
  )
  sim_estimates <- extract_simulated_quantity(stan_data, quantity = "transfers")
  expect_equal(
    sim_estimates$value[1],
    stan_data$stan_data$ABROAD[1] * (1 / stan_data$sim_data$detect_prob[1])
  )
  sim_estimates <- extract_simulated_quantity(
    stan_data,
    quantity = "detect_prob"
  )
  expect_equal(
    sim_estimates$value[1],
    stan_data$sim_data$detect_prob[1]
  )

  stan_data <- simulate_diagnoses(
    nquar = 60,
    sim_seed = 123,
    rita = TRUE,
    migration = TRUE
  )
  sim_estimates <- extract_simulated_quantity(
    stan_data,
    quantity = "migration_prob"
  )
})

test_that("extract_simulated_quantity correctly extracts values for ad", {
  stan_data <- simulate_diagnoses(nquar = 60, sim_seed = 123, age = TRUE)
  sim_estimates <- extract_simulated_quantity(
    stan_data,
    quantity = "incidence_age"
  )
  expect_equal(
    sim_estimates$value[1],
    sum(stan_data$sim_data$infs_age[1:10, 1, 1])
  )
  sim_estimates <- extract_simulated_quantity(
    stan_data,
    quantity = "undiag_prev_age"
  )
  expect_equal(
    sim_estimates$value[1],
    sum(stan_data$sim_data$epi_uk$lat_age[1:10, 1, ]) +
      sum(stan_data$sim_data$epi_uk_abroad$lat_age[1:10, 1, ])
  )
  sim_estimates <- extract_simulated_quantity(
    stan_data,
    quantity = "diagnoses_age"
  )
  expect_equal(
    sim_estimates$value[1],
    sum(
      stan_data$stan_data$HIV_age[1, 1:10] +
        stan_data$stan_data$AIDS_age[1, 1:10]
    )
  )
  sim_estimates <- extract_simulated_quantity(
    stan_data,
    quantity = "diag_prob_age"
  )
  expect_equal(
    sim_estimates$value[1],
    stan_data$sim_data$d_probs[1, 1]
  )
  # age and migration
  stan_data <- simulate_diagnoses(
    nquar = 60,
    sim_seed = 123,
    age = TRUE,
    migration = TRUE
  )
  sim_estimates <- extract_simulated_quantity(
    stan_data,
    quantity = "undiag_migration_age"
  )
  expect_equal(
    sim_estimates$value[1],
    sum(stan_data$sim_data$m_age[1:10, 1, 1:4]) -
      (sum(stan_data$stan_data$ABROAD_age[1, 1:10]) *
        ((1 / stan_data$sim_data$detect_prob[1]) - 1))
  )
  sim_estimates <- extract_simulated_quantity(
    stan_data,
    quantity = "all_migration_age"
  )
  expect_equal(
    sim_estimates$value[1],
    sum(stan_data$sim_data$omega_age[1:10, 1])
  )
})
