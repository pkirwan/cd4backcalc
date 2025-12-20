test_that("basic simulation produces expected HIV data", {
  sim_diags <- simulate_diagnoses(sim_seed = 123)
  expect_equal(
    sim_diags$stan_data$HIV[1],
    604
  )
})

test_that("migration simulation produces expected HIV data", {
  sim_diags <- simulate_diagnoses(sim_seed = 123, migration = TRUE)
  expect_equal(
    sim_diags$stan_data$HIV[1, 1],
    601
  )
})

test_that("rita simulation produces expected HIV data", {
  sim_diags <- simulate_diagnoses(sim_seed = 123, rita = TRUE)
  expect_equal(
    sim_diags$stan_data$HIV[1],
    906
  )
})

test_that("age-specific simulation produces expected HIV data", {
  sim_diags <- simulate_diagnoses(sim_seed = 123, age = TRUE)
  expect_equal(
    sim_diags$stan_data$HIV[1],
    604
  )
  expect_equal(
    sim_diags$stan_data$HIV_age[1],
    11
  )
})

test_that("age-specific rita simulation produces expected HIV data", {
  sim_diags <- simulate_diagnoses(sim_seed = 123, age = TRUE, rita = TRUE)
  expect_equal(
    sim_diags$stan_data$HIV[1],
    906
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
    601
  )
  expect_equal(
    sim_diags$stan_data$HIV_age[1, 1, 1],
    11
  )
})

test_that("age-specific migration and rita simulation produces expected HIV data", {
  sim_diags <- simulate_diagnoses(sim_seed = 123, age = TRUE, migration = TRUE, rita = TRUE)
  expect_equal(
    sim_diags$stan_data$HIV[1, 1],
    755
  )
  expect_equal(
    sim_diags$stan_data$HIV_age[1, 1, 1],
    14
  )
})
