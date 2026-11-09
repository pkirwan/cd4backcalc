skip_if_no_cmdstan <- function() {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("cmdstanr")
  testthat::skip_if(
    !tryCatch(
      {
        cmdstanr::cmdstan_path()
        TRUE
      },
      error = function(e) FALSE
    ),
    "CmdStan is not installed"
  )
}
