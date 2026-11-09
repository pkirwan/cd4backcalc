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
    stan_data$sim_data$d_probs[1, 1, 1]
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
    sum(stan_data$sim_data$m[1, 1:4])
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
    stan_data$sim_data$m[1, 1] / sum(stan_data$sim_data$m[1, ])
  )
  expected_estimates <- extract_simulated_quantity(
    stan_data,
    quantity = "migration_prob",
    deterministic = TRUE
  )
  expect_equal(
    expected_estimates$value[expected_estimates$strata == "Transfers of care"],
    stan_data$sim_data$mu[, 5]
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
    quantity = "transfers_detected"
  )
  expect_equal(
    sim_estimates$value[1],
    stan_data$stan_data$ABROAD[1]
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
  # Check the age-stratified diag_prob output has the expected structure
  expect_true(nrow(sim_estimates) > 0)
  expect_true("age_band" %in% names(sim_estimates))
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
    sum(stan_data$sim_data$m_age[1:10, 1, 1:4])
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

test_that("observed Stan quantities are extracted for plotting", {
  sim_diags <- simulate_diagnoses(
    nquar = 40,
    sim_seed = 123,
    age = TRUE,
    migration = TRUE,
    gamma_obs_sd = 0.1,
    gamma_obs_grouped = TRUE
  )
  quantities <- c(
    "aids_diagnoses",
    "aids_diagnoses_age",
    "transfers_detected_age",
    "gamma_obs",
    "gamma_obs_group",
    "cd4_diagnoses",
    "cd4_diagnoses_age"
  )

  extracted <- lapply(
    quantities,
    function(quantity) extract_simulated_quantity(sim_diags, quantity)
  )

  expect_true(all(vapply(extracted, nrow, integer(1)) > 0))
  expect_true(all(vapply(
    extracted,
    function(x) all(c("value", "year", "group") %in% names(x)),
    logical(1)
  )))
  expect_equal(
    extracted[[which(quantities == "aids_diagnoses")]]$value[1],
    sum(sim_diags$stan_data$AIDS[1, ])
  )
  expect_equal(
    extracted[[which(quantities == "aids_diagnoses_age")]]$value[1],
    sum(sim_diags$stan_data$AIDS_age[1, 1:10, ])
  )
  expect_equal(
    extracted[[which(quantities == "transfers_detected_age")]]$value[1],
    sim_diags$stan_data$ABROAD_group[1, 1]
  )
  expect_equal(
    extracted[[which(quantities == "gamma_obs")]]$value[1],
    exp(sim_diags$stan_data$gamma_obs[1])
  )
  expect_equal(
    extracted[[which(quantities == "gamma_obs_group")]]$value[1],
    exp(sim_diags$stan_data$gamma_obs_group[1, 1])
  )
  expect_false(any(extracted[[which(quantities == "gamma_obs_group")]]$value < 0))
  expect_equal(
    extracted[[which(quantities == "cd4_diagnoses")]]$value[1],
    sum(sim_diags$stan_data$CD4[1, , 1])
  )
  expect_equal(
    extracted[[which(quantities == "cd4_diagnoses_age")]]$value[1],
    sum(sim_diags$stan_data$CD4_age[1, 1:10, , 1])
  )
})
