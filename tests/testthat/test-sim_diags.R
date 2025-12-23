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
    614
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
    614
  )
  expect_equal(
    sim_diags$stan_data$HIV_age[1],
    10
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
  expect_equal(
    sim_diags$stan_data$ABROAD,
    rowSums(sim_diags$stan_data$ABROAD_age)
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
    555
  )
  expect_equal(
    sim_diags$stan_data$HIV_age[1, 1, 1],
    10
  )
  expect_equal(
    sim_diags$stan_data$ABROAD,
    rowSums(sim_diags$stan_data$ABROAD_age)
  )
})

test_that("different sim_type produces expected results", {
  sim_diags <- simulate_diagnoses(nquar = 40, sim_seed = 123, sim_type = "combo_3")
  expect_equal(
    sim_diags$stan_data$HIV[1],
    1686
  )
})
