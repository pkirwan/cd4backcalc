test_that("basic simulation produces expected HIV data", {
  sim_diags <- simulate_diagnoses(sim_seed = 123)
  expect_equal(
    sim_diags$stan_data$HIV[1],
    606
  )
})

test_that("migration simulation produces expected HIV data", {
  sim_diags <- simulate_diagnoses(sim_seed = 123, migration = TRUE)
  expect_equal(
    sim_diags$stan_data$HIV[1, 1],
    599
  )
  expect_equal(
    sim_diags$sim_data$omega[1],
    sum(sim_diags$sim_data$omega_age[1, ])
  )
})

test_that("rita simulation produces expected HIV data", {
  sim_diags <- simulate_diagnoses(sim_seed = 123, rita = TRUE)
  expect_equal(
    sim_diags$stan_data$HIV[1],
    614
  )
})

test_that("age-specific simulation produces expected HIV data", {
  sim_diags <- simulate_diagnoses(sim_seed = 123, age = TRUE)
  expect_equal(
    sim_diags$stan_data$HIV[1],
    606
  )
  expect_equal(
    sim_diags$stan_data$HIV_age[1],
    17
  )
})

test_that("age-specific rita simulation produces expected HIV data", {
  sim_diags <- simulate_diagnoses(sim_seed = 123, age = TRUE, rita = TRUE)
  expect_equal(
    sim_diags$stan_data$HIV[1],
    614
  )
  expect_equal(
    sim_diags$stan_data$HIV_age[1],
    16
  )
})

test_that("age-specific migration simulation produces expected HIV data", {
  sim_diags <- simulate_diagnoses(sim_seed = 123, age = TRUE, migration = TRUE)
  expect_equal(
    sim_diags$stan_data$HIV[1, 1],
    599
  )
  expect_equal(
    sim_diags$stan_data$HIV_age[1, 1, 1],
    14
  )
  expect_equal(
    sim_diags$stan_data$ABROAD,
    rowSums(sim_diags$stan_data$ABROAD_age)
  )
  expect_equal(
    sum(sim_diags$stan_data$ABROAD),
    1880
  )
})

test_that("age-specific migration and rita simulation produces expected HIV data", {
  sim_diags <- simulate_diagnoses(
    sim_seed = 123,
    age = TRUE,
    migration = TRUE,
    rita = TRUE
  )
  expect_equal(
    sim_diags$stan_data$HIV[1, 1],
    584
  )
  expect_equal(
    sim_diags$stan_data$HIV_age[1, 1, 1],
    14
  )
  expect_equal(
    sim_diags$stan_data$ABROAD,
    rowSums(sim_diags$stan_data$ABROAD_age)
  )
  expect_equal(
    sum(sim_diags$stan_data$ABROAD),
    1880
  )
})

test_that("different sim_type produces expected results", {
  sim_diags <- simulate_diagnoses(
    nquar = 40,
    sim_seed = 123,
    sim_type = "combo_3"
  )
  expect_equal(
    sim_diags$stan_data$HIV[1],
    1683
  )
})

test_that("cd4_completeness and rita_completeness parameters affect simulated data as expected", {
  sim_diags <- simulate_diagnoses(
    nquar = 40,
    sim_seed = 123,
    sim_type = "combo_3",
    cd4_completeness = 1
  )
  expect_equal(
    sim_diags$stan_data$CD4[1],
    383
  )

  sim_diags <- simulate_diagnoses(
    nquar = 40,
    sim_seed = 123,
    sim_type = "combo_3",
    cd4_completeness = 0.5
  )
  expect_equal(
    sim_diags$stan_data$CD4[1],
    166
  )
  expect_equal(
    sim_diags$stan_data$CD4_age[1],
    3
  )

  sim_diags_rita <- simulate_diagnoses(
    nquar = 40,
    sim_seed = 123,
    sim_type = "combo_3",
    rita = TRUE,
    cd4_completeness = 1,
    rita_completeness = 1
  )
  expect_equal(
    sim_diags_rita$stan_data$CD4[1],
    501
  )
  expect_equal(
    sim_diags_rita$stan_data$CD4_age[1],
    8
  )

  sim_diags_rita_cd4 <- simulate_diagnoses(
    sim_seed = 123,
    nquar = 40,
    sim_type = "combo_3",
    rita = TRUE,
    cd4_completeness = 1,
    rita_completeness = 0.5
  )
  expect_equal(
    sim_diags_rita_cd4$stan_data$CD4[1],
    264
  )
  expect_equal(
    sim_diags_rita_cd4$stan_data$CD4_age[1],
    4
  )

  sim_diags_rita_cd4 <- simulate_diagnoses(
    sim_seed = 123,
    nquar = 40,
    sim_type = "combo_3",
    rita = TRUE,
    cd4_completeness = 0.5,
    rita_completeness = 1
  )
  expect_equal(
    sim_diags_rita_cd4$stan_data$CD4[1],
    264
  )
  expect_equal(
    sim_diags_rita_cd4$stan_data$CD4_age[1],
    4
  )

  sim_diags_rita_cd4 <- simulate_diagnoses(
    sim_seed = 123,
    nquar = 40,
    sim_type = "combo_3",
    rita = TRUE,
    cd4_completeness = 0.5,
    rita_completeness = 0.5
  )
  expect_equal(
    sim_diags_rita_cd4$stan_data$CD4[1],
    141
  )
  expect_equal(
    sim_diags_rita_cd4$stan_data$CD4_age[1],
    2
  )
})
