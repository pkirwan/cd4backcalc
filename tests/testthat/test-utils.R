test_that("logit works", {
  expect_equal(logit(0.5), 0)
})

test_that("invlogit works", {
  expect_equal(invlogit(0), 0.5)
})

test_that("round_preserve works", {
  x <- c(1.2, 2.5, 3.3)
  rounded <- round_preserve(x)
  expect_equal(sum(rounded), round(sum(x)))
})

test_that("round_matrix_preserve works", {
  mat <- matrix(c(1.2, 2.5, 3.3, 4.7), nrow = 2)
  rounded_mat <- round_matrix_preserve(mat)
  expect_equal(sum(rounded_mat), round(sum(mat)))
})

test_that("round_array_preserve works", {
  array <- array(c(1.2, 2.5, 3.3, 4.7, 5.6, 6.1), dim = c(2, 3))
  rounded_array <- round_array_preserve(array)
  expect_equal(sum(unlist(rounded_array)), round(sum(unlist(array))))
})

test_that("round_fields works", {
  epi_obj <- list(
    lat = matrix(c(1.2, 2.5, 3.3, 4.7), nrow = 2),
    diags = matrix(c(2.1, 3.4, 1.6, 5.8), nrow = 2)
  )
  rounded_epi <- round_fields(epi_obj, fields = c("lat", "diags"))
  expect_equal(sum(rounded_epi$lat), round(sum(epi_obj$lat)))
  expect_equal(sum(rounded_epi$diags), round(sum(epi_obj$diags)))
})
