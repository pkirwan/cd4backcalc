test_that("stan_model_cache_dir prefers explicit override", {
  original_override <- Sys.getenv("CD4BACKCALC_CMDSTAN_DIR", unset = NA)
  original_tmpdir <- Sys.getenv("TMPDIR", unset = NA)

  override_dir <- file.path(tempdir(), "cmdstan-override")

  on.exit({
    if (is.na(original_override)) {
      Sys.unsetenv("CD4BACKCALC_CMDSTAN_DIR")
    } else {
      Sys.setenv(CD4BACKCALC_CMDSTAN_DIR = original_override)
    }

    if (is.na(original_tmpdir)) {
      Sys.unsetenv("TMPDIR")
    } else {
      Sys.setenv(TMPDIR = original_tmpdir)
    }
  }, add = TRUE)

  Sys.setenv(CD4BACKCALC_CMDSTAN_DIR = override_dir)
  Sys.setenv(TMPDIR = file.path(tempdir(), "ignored-tmpdir"))

  expect_equal(cd4backcalc:::stan_model_cache_dir(), normalizePath(override_dir, winslash = "/", mustWork = FALSE))
  expect_true(dir.exists(override_dir))
})

test_that("stan_model_cache_dir falls back to tempdir cache", {
  original_override <- Sys.getenv("CD4BACKCALC_CMDSTAN_DIR", unset = NA)

  on.exit({
    if (is.na(original_override)) {
      Sys.unsetenv("CD4BACKCALC_CMDSTAN_DIR")
    } else {
      Sys.setenv(CD4BACKCALC_CMDSTAN_DIR = original_override)
    }
  }, add = TRUE)

  Sys.unsetenv("CD4BACKCALC_CMDSTAN_DIR")

  cache_dir <- cd4backcalc:::stan_model_cache_dir()

  expect_equal(
    cache_dir,
    normalizePath(
      file.path(tempdir(), "cd4backcalc-cmdstan"),
      winslash = "/",
      mustWork = FALSE
    )
  )
  expect_true(dir.exists(cache_dir))
})
