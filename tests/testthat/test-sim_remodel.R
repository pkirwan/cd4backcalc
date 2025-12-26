test_that("basic re-modelling process works", {
  # Simulate diagnoses
  sim_diags_rita <- simulate_diagnoses(nquar = 60, sim_seed = 123, rita = TRUE)
  sim_diags_remod <- sim_remodel(sim_diags_rita, rita = FALSE)

  # confirm remodeled stan data has no rita states
  expect_equal(ncol(sim_diags_remod$stan_data$CD4), 4)
  expect_equal(length(sim_diags_remod$stan_data$init_prev), 4)
  # check total CD4 diagnoses preserved
  expect_equal(
    sum(sim_diags_rita$stan_data$CD4),
    sum(sim_diags_remod$stan_data$CD4)
  )
  expect_equal(
    sum(sim_diags_rita$stan_data$CD4[3, ]),
    sum(sim_diags_remod$stan_data$CD4[3, ])
  )
  expect_equal(
    sum(sim_diags_rita$stan_data$CD4),
    sum(sim_diags_remod$stan_data$HIV)
  )
  expect_equal(
    rowSums(sim_diags_remod$stan_data$CD4),
    sim_diags_remod$stan_data$HIV
  )
  # check total age-specific CD4 diagnoses preserved
  expect_equal(
    sum(sim_diags_rita$stan_data$CD4_age),
    sum(sim_diags_remod$stan_data$CD4_age)
  )
  expect_equal(
    sum(sim_diags_rita$stan_data$CD4_age[10, , ]),
    sum(sim_diags_remod$stan_data$CD4_age[10, , ])
  )
  expect_equal(
    sum(sim_diags_rita$stan_data$CD4_age[1, 1, ]),
    sum(sim_diags_remod$stan_data$CD4_age[1, 1, ])
  )
  expect_equal(
    rowSums(colSums(sim_diags_remod$stan_data$CD4_age)),
    colSums(sim_diags_remod$stan_data$HIV_age)
  )
})
