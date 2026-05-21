make_fake_fit <- function(draw_matrix) {
  list(draws = function(format = "matrix") draw_matrix)
}

test_that("calibrate_init_prev recovers age-independent CD4-only baseline", {
  draw_matrix <- matrix(
    c(0.2, 0.3, 0.4, 0.5),
    nrow = 1,
    dimnames = list(
      NULL,
      c("d[1,1]", "d[2,1]", "d[3,1]", "d[4,1]")
    )
  )

  stan_data <- list(
    HIV = c(17, 0, 0, 0),
    AIDS = c(2, 0, 0, 0),
    CD4 = rbind(c(2, 3, 4, 8), c(0, 0, 0, 0), c(0, 0, 0, 0), c(0, 0, 0, 0)),
    q = c(0.1, 0.1, 0.1, 0.25),
    init_prev = c(1, 1, 1, 1),
    nquar = 4,
    model_config = list(age = FALSE, migration = FALSE, rita = FALSE)
  )

  out <- calibrate_init_prev(make_fake_fit(draw_matrix), stan_data, smooth_age = FALSE)

  expect_equal(out$init_prev, c(10, 10, 10, 16), tolerance = 1e-10)
  expect_equal(out$stan_data$init_prev, out$init_prev, tolerance = 1e-10)
})

test_that("calibrate_init_prev recovers age-independent migration RITA baselines", {
  native_init <- c(4, 4, 4, 10, 8, 10, 20)
  migrant_init <- c(1, 2, 3, 5, 4, 3, 10)
  frr <- 0.1

  native_obs <- c(7.5, 1.8, 1.8, 2.7, 7.2)
  migrant_obs <- c(4.95, 1.575, 1.98, 2.295, 6.3)

  draw_matrix <- matrix(
    c(
      0.5, 0.2, 0.25, 0.3, 0.4,
      0.6, 0.25, 0.3, 0.35, 0.5,
      5,
      0.1, 0.2, 0.3, 0.4,
      frr
    ),
    nrow = 1,
    dimnames = list(
      NULL,
      c(
        "d_nat[1,1]", "d_nat[2,1]", "d_nat[3,1]", "d_nat[4,1]", "d_nat[5,1]",
        "d_mig[1,1]", "d_mig[2,1]", "d_mig[3,1]", "d_mig[4,1]", "d_mig[5,1]",
        "misclassified_dx[1,2]",
        "m[1,1]", "m[2,1]", "m[3,1]", "m[4,1]",
        "frr"
      )
    )
  )

  cd4 <- array(0, dim = c(4, 2, 5))
  cd4[1, 1, ] <- native_obs
  cd4[1, 2, ] <- migrant_obs

  stan_data <- list(
    HIV = rbind(c(sum(native_obs), sum(migrant_obs)), c(0, 0), c(0, 0), c(0, 0)),
    AIDS = rbind(c(3, 1.25), c(0, 0), c(0, 0), c(0, 0)),
    CD4 = cd4,
    q = c(0.1, 0.2, 0.3, 0.25),
    init_prev = native_init,
    init_prev_abroad = migrant_init,
    nquar = 4,
    model_config = list(age = FALSE, migration = TRUE, rita = TRUE)
  )

  out <- calibrate_init_prev(make_fake_fit(draw_matrix), stan_data, smooth_age = FALSE)

  expect_equal(out$init_prev, native_init, tolerance = 1e-10)
  expect_equal(out$init_prev_abroad, migrant_init, tolerance = 1e-10)
})

test_that("calibrate_init_prev recovers age-dependent migration CD4-only baselines", {
  nage <- 12
  native_init <- matrix(0, nrow = 4, ncol = nage)
  migrant_init <- matrix(0, nrow = 4, ncol = nage)
  native_init[, 1] <- c(10, 10, 10, 16)
  native_init[, 11] <- c(8, 9, 10, 14)
  migrant_init[, 1] <- c(6, 7, 8, 12)
  migrant_init[, 11] <- c(5, 6, 7, 11)

  native_obs_age1 <- c(2, 3, 4, 8)
  native_obs_age11 <- c(2, 3.15, 4.5, 7.7)
  migrant_base_age1 <- c(1.32, 2.24, 3.36, 6.24)
  migrant_base_age11 <- c(1.35, 2.22, 3.29, 6.27)
  migrant_obs_age1 <- migrant_base_age1 + c(0.1, 0.2, 0.3, 0.4)
  migrant_obs_age11 <- migrant_base_age11 + c(0.2, 0.4, 0.6, 0.8)

  draw_values <- c(
    0.2, 0.3, 0.4, 0.5,
    0.25, 0.35, 0.45, 0.55,
    0.2, 0.3, 0.4, 0.5,
    0.2, 0.3, 0.4, 0.5,
    0.22, 0.32, 0.42, 0.52,
    0.27, 0.37, 0.47, 0.57,
    0.22, 0.32, 0.42, 0.52,
    0.22, 0.32, 0.42, 0.52,
    c(1, rep(0, 9), 2, 0),
    0.1, 0.2, 0.3, 0.4
  )

  draw_names <- c(
    unlist(lapply(1:4, function(ag) sprintf("d_nat[%d,1,%d]", ag, 1:4))),
    unlist(lapply(1:4, function(ag) sprintf("d_mig[%d,1,%d]", ag, 1:4))),
    sprintf("misclassified_dx[1,%d]", seq_len(nage)),
    sprintf("m[%d,1]", 1:4)
  )

  draw_matrix <- matrix(draw_values, nrow = 1, dimnames = list(NULL, draw_names))

  hiv_age <- array(0, dim = c(4, nage, 2))
  aids_age <- array(0, dim = c(4, nage, 2))
  cd4_age <- array(0, dim = c(4, nage, 2, 4))

  hiv_age[1, 1, 1] <- sum(native_obs_age1)
  hiv_age[1, 11, 1] <- sum(native_obs_age11)
  hiv_age[1, 1, 2] <- sum(migrant_obs_age1)
  hiv_age[1, 11, 2] <- sum(migrant_obs_age11)

  aids_age[1, 1, 1] <- 2
  aids_age[1, 11, 1] <- 1.575
  aids_age[1, 1, 2] <- 1.44
  aids_age[1, 11, 2] <- 1.1825

  cd4_age[1, 1, 1, ] <- native_obs_age1
  cd4_age[1, 11, 1, ] <- native_obs_age11
  cd4_age[1, 1, 2, ] <- migrant_obs_age1
  cd4_age[1, 11, 2, ] <- migrant_obs_age11

  q_age <- matrix(rep(c(0.1, 0.1, 0.1, 0.25), nage), nrow = 4)

  stan_data <- list(
    HIV_age = hiv_age,
    AIDS_age = aids_age,
    CD4_age = cd4_age,
    q_age = q_age,
    init_prev_age = native_init,
    init_prev_abroad_age = migrant_init,
    nage = nage,
    nquar = 4,
    model_config = list(age = TRUE, migration = TRUE, rita = FALSE)
  )

  out <- calibrate_init_prev(make_fake_fit(draw_matrix), stan_data, smooth_age = FALSE)

  expect_equal(out$init_prev_age, native_init, tolerance = 1e-10)
  expect_equal(out$init_prev_abroad_age, migrant_init, tolerance = 1e-10)
})