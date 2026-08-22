test_that("migration models disable melding when eta_sd is zero", {
  model_names <- c(
    "AD_ptens_migration.stan",
    "AD_ptens_migration_reduced.stan",
    "AD_ptens_migration_rita.stan",
    "AI_spl_migration.stan",
    "AI_spl_migration_rita.stan"
  )

  purrr::walk(model_names, function(model_name) {
    model_path <- system.file("stan", model_name, package = "cd4backcalc")
    expect_true(nzchar(model_path), info = model_name)
    model_code <- paste(readLines(model_path, warn = FALSE), collapse = "\n")
    expect_match(model_code, "if \\(eta_sd > 0\\)", info = model_name)
  })
})