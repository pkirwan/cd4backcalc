libs <- file.path(R_PACKAGE_DIR, "libs", R_ARCH)
dir.create(libs, recursive = TRUE, showWarnings = FALSE)
for (file in c("symbols.rds", Sys.glob(paste0("*", SHLIB_EXT)))) {
  if (file.exists(file)) {
    file.copy(file, file.path(libs, file))
  }
}
inst_stan <- file.path("..", "inst", "stan")
if (dir.exists(inst_stan)) {
  warning(
    "Stan models in inst/stan/ are deprecated in {instantiate} ",
    ">= 0.0.4.9001 (2024-01-03). Please put them in src/stan/ instead."
  )
  if (file.exists("stan")) {
    warning("src/stan/ already exists. Not copying models from inst/stan/.")
  } else {
    message("Copying inst/stan/ to src/stan/.")
    fs::dir_copy(path = inst_stan, new_path = "stan")
  }
}
bin <- file.path(R_PACKAGE_DIR, "bin")
if (!file.exists(bin)) {
  dir.create(bin, recursive = TRUE, showWarnings = FALSE)
}
bin_stan <- file.path(bin, "stan")

# Compile Stan models in a temp directory to avoid permission issues
# that occur when compiling directly in the package install directory.
# Copying files into the install directory can strip execute bits on Linux.
tmp_stan <- file.path(tempdir(), "stan_compile")
if (dir.exists(tmp_stan)) unlink(tmp_stan, recursive = TRUE)
fs::dir_copy(path = "stan", new_path = tmp_stan)
callr::r(
  func = function(tmp_stan) {
    instantiate::stan_package_compile(
      models = instantiate::stan_package_model_files(path = tmp_stan)
    )
  },
  args = list(tmp_stan = tmp_stan),
  show = TRUE,
  stderr = "2>&1"
)

# Copy compiled output to the final package directory
fs::dir_copy(path = tmp_stan, new_path = bin_stan)

# Restore execute permissions on compiled Stan binaries (file copy can
# strip the execute bit on Linux/macOS)
if (.Platform$OS.type == "unix") {
  stan_bins <- list.files(bin_stan, full.names = TRUE)
  stan_bins <- stan_bins[!grepl("\\.", basename(stan_bins))]
  for (b in stan_bins) {
    Sys.chmod(b, mode = "0755")
  }
}
