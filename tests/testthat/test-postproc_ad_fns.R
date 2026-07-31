test_that("total UK-acquired incidence applies the ratio by quarter and age", {
  incidence_uk <- matrix(
    c(
      10, 20, 30, 40, 50, 60,
      5, 10, 15, 20, 25, 30
    ),
    nrow = 2,
    byrow = TRUE
  )
  ratio_abroad_uk <- matrix(
    c(
      0, 1, 3,
      1, 0, 0.5
    ),
    nrow = 2,
    byrow = TRUE
  )

  result <- .total_uk_acquired_log_incidence(
    log(incidence_uk),
    ratio_abroad_uk,
    nquar = 3,
    nage = 2
  )

  quarter_index <- rep(1:3, times = 2)
  expected <- incidence_uk * (1 + ratio_abroad_uk[, quarter_index])
  expect_equal(exp(result), expected)
})

test_that("total UK-acquired incidence rejects invalid ratios", {
  log_incidence_uk <- matrix(log(1:6), nrow = 1)

  expect_error(
    .total_uk_acquired_log_incidence(
      log_incidence_uk,
      matrix(c(0, -1, 1), nrow = 1),
      nquar = 3,
      nage = 2
    ),
    "finite non-negative"
  )
})