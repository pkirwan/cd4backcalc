test_that("basic simulation produces expected HIV data", {
  sim_diags <- simulate_diagnoses(sim_seed = 123)
  expect_equal(
    sim_diags$stan_data$HIV[1],
    722
  )
})

test_that("migration simulation produces expected HIV data", {
  sim_diags <- simulate_diagnoses(sim_seed = 123, migration = TRUE)
  expect_equal(
    sim_diags$stan_data$HIV[1, 1],
    718
  )
})

test_that("rita simulation produces expected HIV data", {
  sim_diags <- simulate_diagnoses(sim_seed = 123, rita = TRUE)
  expect_equal(
    sim_diags$stan_data$HIV[1],
    1082
  )
})

test_that("age-specific simulation produces expected HIV data", {
  sim_diags <- simulate_diagnoses(sim_seed = 123, age = TRUE)
  expect_equal(
    sim_diags$stan_data$HIV[1],
    722
  )
  expect_equal(
    sim_diags$stan_data$HIV_age[1],
    14
  )
})

test_that("age-specific rita simulation produces expected HIV data", {
  sim_diags <- simulate_diagnoses(sim_seed = 123, age = TRUE, rita = TRUE)
  expect_equal(
    sim_diags$stan_data$HIV[1],
    1082
  )
  expect_equal(
    sim_diags$stan_data$HIV_age[1],
    21
  )
})

test_that("age-specific migration simulation produces expected HIV data", {
  sim_diags <- simulate_diagnoses(sim_seed = 123, age = TRUE, migration = TRUE)
  expect_equal(
    sim_diags$stan_data$HIV[1, 1],
    718
  )
  expect_equal(
    sim_diags$stan_data$HIV_age[1, 1, 1],
    14
  )
})

test_that("age-specific migration and rita simulation produces expected HIV data", {
  sim_diags <- simulate_diagnoses(sim_seed = 123, age = TRUE, migration = TRUE, rita = TRUE)
  expect_equal(
    sim_diags$stan_data$HIV[1, 1],
    901
  )
  expect_equal(
    sim_diags$stan_data$HIV_age[1, 1, 1],
    18
  )
})
