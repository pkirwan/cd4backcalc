test_that("basic simulation produces expected HIV data", {
  sim_config <- define_sim_config(seed = 123)
  basic_simulation <- simulate_diagnoses(sim_config)
  expect_equal(
    basic_simulation$stan_data$HIV[1],
    722
  )
})

test_that("migration simulation produces expected HIV data", {
  sim_config <- define_sim_config(seed = 123, migration = TRUE)
  migration_simulation <- simulate_diagnoses(sim_config)
  expect_equal(
    migration_simulation$stan_data$HIV[1, 1],
    718
  )
})

test_that("rita simulation produces expected HIV data", {
  sim_config <- define_sim_config(seed = 123, rita = TRUE)
  rita_simulation <- simulate_diagnoses(sim_config)
  expect_equal(
    rita_simulation$stan_data$HIV[1],
    1082
  )
})

test_that("age-specific simulation produces expected HIV data", {
  sim_config <- define_sim_config(seed = 123, age = TRUE)
  age_simulation <- simulate_diagnoses(sim_config)
  expect_equal(
    age_simulation$stan_data$HIV[1],
    722
  )
  expect_equal(
    age_simulation$stan_data$HIV_age[1],
    13
  )
})

test_that("age-specific rita simulation produces expected HIV data", {
  sim_config <- define_sim_config(seed = 123, age = TRUE, rita = TRUE)
  age_rita_simulation <- simulate_diagnoses(sim_config)
  expect_equal(
    age_rita_simulation$stan_data$HIV[1],
    1082
  )
  expect_equal(
    age_rita_simulation$stan_data$HIV_age[1],
    21
  )
})

test_that("age-specific migration simulation produces expected HIV data", {
  sim_config <- define_sim_config(seed = 123, age = TRUE, migration = TRUE)
  age_migration_simulation <- simulate_diagnoses(sim_config)
  expect_equal(
    age_migration_simulation$stan_data$HIV[1, 1],
    718
  )
  expect_equal(
    age_migration_simulation$stan_data$HIV_age[1, 1, 1],
    14
  )
})

test_that("age-specific migration and rita simulation produces expected HIV data", {
  sim_config <- define_sim_config(seed = 123, age = TRUE, migration = TRUE, rita = TRUE)
  age_migration_rita_simulation <- simulate_diagnoses(sim_config)
  expect_equal(
    age_migration_rita_simulation$stan_data$HIV[1, 1],
    901
  )
  expect_equal(
    age_migration_rita_simulation$stan_data$HIV_age[1, 1, 1],
    18
  )
})
