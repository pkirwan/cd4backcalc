# Pre-render figures for vignette articles.
#
# Run from repository root:
#   Rscript tools/make-figures.R
#
# Optional environment variables:
#   CD4BC_ITER_WARMUP, CD4BC_ITER_SAMPLING, CD4BC_CHAINS,
#   CD4BC_AGE_ITER_WARMUP, CD4BC_AGE_ITER_SAMPLING
#
# The script writes static image files into vignettes/articles/figures/
# so article .Rmd files can include images without running model-fitting code.

if (!requireNamespace("cd4backcalc", quietly = TRUE)) {
  stop("Package 'cd4backcalc' must be installed before running this script.")
}
if (!requireNamespace("ggplot2", quietly = TRUE)) {
  stop("Package 'ggplot2' must be installed before running this script.")
}
if (!requireNamespace("patchwork", quietly = TRUE)) {
  stop("Package 'patchwork' must be installed before running this script.")
}

library(cd4backcalc)
library(ggplot2)
library(patchwork)

out_dir <- file.path("vignettes", "articles", "figures")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
load_saved_models <- TRUE

readme_out_dir <- file.path("man", "figures")
dir.create(readme_out_dir, recursive = TRUE, showWarnings = FALSE)

iter_warmup <- as.integer(Sys.getenv("CD4BC_ITER_WARMUP", "1000"))
iter_sampling <- as.integer(Sys.getenv("CD4BC_ITER_SAMPLING", "1000"))
chains <- as.integer(Sys.getenv("CD4BC_CHAINS", "4"))

age_iter_warmup <- as.integer(Sys.getenv("CD4BC_AGE_ITER_WARMUP", "200"))
age_iter_sampling <- as.integer(Sys.getenv("CD4BC_AGE_ITER_SAMPLING", "200"))

fit_args <- list(
  chains = chains,
  iter_warmup = iter_warmup,
  iter_sampling = iter_sampling,
  model_seed = 123
)

fit_args_age <- list(
  chains = chains,
  iter_warmup = age_iter_warmup,
  iter_sampling = age_iter_sampling,
  model_seed = 123
)

save_plot <- function(plot_obj, name, width = 8, height = 4.8, dpi = 320) {
  png_path <- file.path(out_dir, paste0(name, ".png"))
  pdf_path <- file.path(out_dir, paste0(name, ".pdf"))

  ggsave(
    filename = png_path,
    plot = plot_obj,
    width = width,
    height = height,
    dpi = dpi
  )
  ggsave(filename = pdf_path, plot = plot_obj, width = width, height = height)

  message("Saved: ", png_path)
  message("Saved: ", pdf_path)
}

save_readme_plot <- function(
  plot_obj,
  name,
  width = 10,
  height = 6,
  dpi = 320
) {
  png_path <- file.path(readme_out_dir, paste0(name, ".png"))

  ggsave(
    filename = png_path,
    plot = plot_obj,
    width = width,
    height = height,
    dpi = dpi
  )
  message("Saved: ", png_path)
}

save_estimate_plot <- function(fit_obj, quantity, name, ...) {
  p <- plot_estimates(fit_obj, quantity = quantity, ...)
  save_plot(p, name)
}

save_minimal_plot <- function(fit_obj, quantity, name, ...) {
  p <- plot_estimates(fit_obj, quantity = quantity, ...) +
    theme_gray(14)
  save_plot(p, name)
}

save_sim_plot <- function(sim_obj, quantity, name, ...) {
  p <- plot_simulations(sim_obj, quantity = quantity, ...)
  save_plot(p, name)
}

fit_model <- function(stan_data, age = FALSE, ...) {
  args <- if (age) fit_args_age else fit_args
  extra <- list(...)
  args[names(extra)] <- extra
  args <- c(list(stan_data = stan_data), args)
  do.call(run_backcalc, args)
}

healthy_fit <- function(object, quantity = NULL, max_truth_ratio = Inf) {
  diagnostics_ok <- !is.null(object$n.divg) &&
    object$n.divg == 0 &&
    !is.null(object$conv.chk) &&
    all(object$conv.chk)

  if (!diagnostics_ok || is.null(quantity)) {
    return(diagnostics_ok)
  }

  estimates <- get_estimates(object, quantity = quantity)
  truth <- extract_simulated_quantity(object, quantity = quantity)
  all(is.finite(estimates$mean)) &&
    max(estimates$mean) <= max_truth_ratio * max(truth$value)
}

load_or_build_cached <- function(
  cache_file,
  object_name,
  build_fn,
  label,
  validate_fn = function(object) TRUE
) {
  cache_path <- file.path(out_dir, cache_file)

  if (load_saved_models && file.exists(cache_path)) {
    env <- new.env(parent = emptyenv())
    load(cache_path, envir = env)

    if (!exists(object_name, envir = env, inherits = FALSE)) {
      stop("Cache file does not contain expected object: ", object_name)
    }

    object <- get(object_name, envir = env, inherits = FALSE)
    if (validate_fn(object)) {
      message("Loaded saved ", label, ".")
      return(object)
    }

    message("Saved ", label, " failed validation. Rebuilding...")
  }

  message("No saved ", label, " found. Building...")
  obj <- build_fn()

  if (!validate_fn(obj)) {
    stop("Newly built ", label, " failed validation.")
  }

  assign(object_name, obj)
  save(list = object_name, file = cache_path)

  obj
}

message("Preparing simulated datasets...")
set.seed(140926)

message("Generating figures from README.Rmd...")
sim_readme <- simulate_diagnoses(h_pattern = "increasing")
fit_readme <- load_or_build_cached(
  cache_file = "README-fit-sim.RData",
  object_name = "fit_readme",
  build_fn = function() run_backcalc(sim_readme),
  label = "model fit for README simulation"
)

save_readme_plot(
  plot_simulations(sim_readme, quantity = "incidence"),
  "README-simulate-data-1"
)
save_readme_plot(plot_diagnostics(fit_readme), "README-simulate-data-2")
save_readme_plot(
  bias_plot(list(fit_readme), quantity = "incidence"),
  "README-simulate-data-3"
)
save_readme_plot(
  plot_estimates(fit_readme, quantity = "incidence"),
  "README-simulate-data-4"
)

sim_cd4 <- simulate_diagnoses(sim_type = "combo_3", sim_seed = 123)
sim_age <- simulate_diagnoses(sim_type = "combo_3", age = TRUE, sim_seed = 124)
sim_rita <- simulate_diagnoses(
  sim_type = "combo_3",
  rita = TRUE,
  sim_seed = 125
)
sim_mig <- simulate_diagnoses(
  sim_type = "combo_3",
  migration = TRUE,
  sim_seed = 126
)
sim_rita_mig <- simulate_diagnoses(
  sim_type = "combo_4",
  rita = TRUE,
  migration = TRUE,
  sim_seed = 127
)

message("Generating figures from vignettes/cd4backcalc.Rmd...")
save_sim_plot(sim_cd4, "incidence", "getting-started-sim-incidence")

fit_cd4 <- load_or_build_cached(
  cache_file = "getting-started-fit-cd4.RData",
  object_name = "fit_cd4",
  build_fn = function() fit_model(sim_cd4),
  label = "model fit for CD4-only model"
)
save_plot(
  plot_diagnostics(fit_cd4, ntraces = 4),
  "getting-started-diagnostics",
  width = 10,
  height = 8
)
save_plot(
  bias_plot(list(fit_cd4), quantity = "incidence"),
  "getting-started-bias-incidence"
)
p1 <- plot_estimates(fit_cd4, quantity = "incidence")
p2 <- plot_estimates(fit_cd4, quantity = "undiag_prev")
p3 <- plot_estimates(fit_cd4, quantity = "diag_prob")
p4 <- plot_estimates(fit_cd4, quantity = "diagnoses")
save_plot(
  (p1 + p2) /
    p3 /
    p4 +
    plot_layout(guides = "collect") &
    theme(legend.position = "bottom"),
  "getting-started-estimates-patchwork",
  width = 10,
  height = 13
)
save_minimal_plot(fit_cd4, "incidence", "getting-started-est-incidence-minimal")

message("Generating figures from vignettes/articles/model-types.Rmd...")
p1 <- plot_estimates(fit_cd4, quantity = "incidence")
p2 <- plot_estimates(fit_cd4, quantity = "undiag_prev")
p3 <- plot_estimates(fit_cd4, quantity = "diag_prob")
save_plot(
  (p1 + p2) /
    p3 +
    plot_layout(guides = "collect") &
    theme(legend.position = "bottom"),
  "model-types-cd4-patchwork",
  width = 10,
  height = 8
)

fit_age <- load_or_build_cached(
  cache_file = "model-types-fit-age.RData",
  object_name = "fit_age",
  build_fn = function() fit_model(sim_age, age = TRUE, return_fit = FALSE),
  label = "model fit for age-dependent model"
)
p1 <- plot_estimates(fit_age, quantity = "incidence_age")
p2 <- plot_estimates(fit_age, quantity = "undiag_prev_age")
p3 <- plot_estimates(fit_age, quantity = "diagnoses_age")
p4 <- plot_estimates(fit_age, quantity = "incidence")
save_plot(
  p1 /
    p2 /
    p3 /
    p4 +
    plot_layout(guides = "collect") &
    theme(legend.position = "bottom"),
  "model-types-age-patchwork",
  width = 10,
  height = 14
)

fit_rita <- load_or_build_cached(
  cache_file = "model-types-fit-rita.RData",
  object_name = "fit_rita",
  build_fn = function() fit_model(sim_rita, rita = TRUE, return_fit = FALSE),
  label = "model fit for RITA model"
)
p1 <- plot_estimates(fit_rita, quantity = "incidence")
p2 <- plot_estimates(fit_rita, quantity = "undiag_prev")
p3 <- plot_estimates(fit_rita, quantity = "diag_prob")
save_plot(
  (p1 + p2) /
    p3 +
    plot_layout(guides = "collect") &
    theme(legend.position = "bottom"),
  "model-types-rita-patchwork",
  width = 10,
  height = 8
)

fit_mig <- load_or_build_cached(
  cache_file = "model-types-fit-migration.RData",
  object_name = "fit_mig",
  build_fn = function() {
    fit_model(sim_mig, migration = TRUE, return_fit = FALSE, init = 0)
  },
  label = "model fit for migration model",
  validate_fn = healthy_fit
)
p1 <- plot_estimates(fit_mig, quantity = "undiag_migration")
p2 <- plot_estimates(fit_mig, quantity = "all_migration")
p3 <- plot_estimates(fit_mig, quantity = "ratio_abroad_uk")
p4 <- plot_estimates(fit_mig, quantity = "detect_prob")
p5 <- plot_estimates(fit_mig, quantity = "migration_prob")
p6 <- plot_estimates(fit_mig, quantity = "diag_prob_mig")
save_plot(
  (p1 + p2) /
    (p3 + p4) /
    p5 /
    p6 +
    plot_layout(guides = "collect") &
    theme(legend.position = "bottom"),
  "model-types-migration-patchwork",
  width = 10,
  height = 12
)

fit_rita_mig <- load_or_build_cached(
  cache_file = "model-types-fit-rita-migration.RData",
  object_name = "fit_rita_mig",
  build_fn = function() {
    fit_model(sim_rita_mig, rita = TRUE, migration = TRUE, return_fit = FALSE)
  },
  label = "model fit for RITA + migration model"
)
p1 <- plot_estimates(fit_rita_mig, quantity = "incidence")
p2 <- plot_estimates(fit_rita_mig, quantity = "undiag_migration")
p3 <- plot_estimates(fit_rita_mig, quantity = "diag_prob")
p4 <- plot_estimates(fit_rita_mig, quantity = "migration_prob")
save_plot(
  (p1 + p2) /
    p3 /
    p4 +
    plot_layout(guides = "collect") &
    theme(legend.position = "bottom"),
  "model-types-rita-migration-patchwork",
  width = 10,
  height = 8
)

fit_cd4_from_rita <- load_or_build_cached(
  cache_file = "model-types-fit-remodel-cd4.RData",
  object_name = "fit_cd4_from_rita",
  build_fn = function() {
    sim_cd4_from_rita <- sim_remodel(sim_rita, rita = FALSE)
    fit_model(sim_cd4_from_rita, rita = FALSE, return_fit = FALSE)
  },
  label = "model fit for CD4-only model fitted to RITA-simulated data"
)

message("Generating figures from vignettes/articles/analysing-results.Rmd...")
p1 <- plot_simulations(sim_mig, quantity = "incidence")
p2 <- plot_simulations(sim_mig, quantity = "undiag_prev")
p3 <- plot_simulations(sim_mig, quantity = "all_migration")
p4 <- plot_simulations(sim_mig, quantity = "undiag_migration")
p5 <- plot_simulations(sim_mig, quantity = "diag_prob")
save_plot(
  (p1 + p2) /
    (p3 + p4) /
    p5 +
    plot_layout(guides = "collect") &
    theme(legend.position = "bottom"),
  "analysing-results-sim-patchwork",
  width = 10,
  height = 12
)

p1 <- plot_estimates(fit_mig, quantity = "incidence")
p2 <- plot_estimates(fit_mig, quantity = "undiag_prev")
p3 <- plot_estimates(fit_mig, quantity = "diag_prob")
save_plot(
  (p1 + p2) /
    p3 +
    plot_layout(guides = "collect") &
    theme(legend.position = "bottom"),
  "analysing-results-estimates-patchwork",
  width = 10,
  height = 8
)
save_estimate_plot(
  fit_mig,
  "incidence",
  "analysing-results-est-incidence-generative",
  sim_comparison = "generative"
)

p1 <- bias_plot(list(fit_mig), quantity = "incidence")
p2 <- mse_plot(list(fit_mig), quantity = "incidence")

save_plot(
  p1 + p2,
  "analysing-results-bias-mse-incidence",
  width = 10,
  height = 6
)

inc_df <- get_estimates(fit_mig, quantity = "incidence")
head(inc_df)

sim_inc <- extract_simulated_quantity(sim_mig, quantity = "incidence")
head(sim_inc)

real_mig <- sim_mig
real_mig$stan_data$model_config <- real_mig$stan_data$sim_config
obs <- extract_real_data(real_mig, quantity = "diagnoses")
head(obs)

compare_two_est <- plot_estimates(
  fit_rita,
  hiv_list_2 = fit_cd4_from_rita,
  quantity = "incidence",
  model_1_name = "CD4 + RITA",
  model_2_name = "CD4 only"
)
compare_two_bias <- bias_plot(
  list("CD4 + RITA" = fit_rita, "CD4 only" = fit_cd4_from_rita),
  quantity = "incidence",
  colours = c("CD4 + RITA" = "mediumorchid", "CD4 only" = "coral"),
  offset = 15
)

save_plot(
  compare_two_est /
    compare_two_bias,
  "analysing-results-compare-two-patchwork",
  width = 10,
  height = 10
)

models <- load_or_build_cached(
  cache_file = "analysing-results-compare-smoothing-models.RData",
  object_name = "models",
  build_fn = function() {
    list(
      "Spline incidence" = fit_model(
        sim_cd4,
        inf_model = 1,
        diag_model = 1,
        model_seed = 123,
        return_fit = FALSE
      ),
      "Random-walk incidence" = fit_model(
        sim_cd4,
        inf_model = 2,
        diag_model = 1,
        model_seed = 123,
        return_fit = FALSE
      ),
      "Gaussian-process incidence" = fit_model(
        sim_cd4,
        inf_model = 3,
        diag_model = 1,
        model_seed = 123,
        return_fit = FALSE
      )
    )
  },
  label = "model fits for incidence-smoothing comparison",
  validate_fn = function(fits) {
    all(vapply(
      fits,
      healthy_fit,
      logical(1),
      quantity = "incidence",
      max_truth_ratio = 10
    ))
  }
)
compare_list_est <- plot_estimates(models, quantity = "incidence")
compare_list_bias <- bias_plot(models, quantity = "incidence", offset = 6)

save_plot(
  compare_list_est /
    compare_list_bias,
  "analysing-results-compare-smoothing-patchwork",
  width = 10,
  height = 10
)

fit_chain <- load_or_build_cached(
  cache_file = "analysing-results-chain-results.RData",
  object_name = "fit_chain",
  build_fn = function() {
    fit_model(
      sim_cd4,
      chain_results = TRUE,
      return_fit = FALSE,
      iter_warmup = 50,
      iter_sampling = 50,
      chains = 3
    )
  },
  label = "model fit for chain results"
)
save_estimate_plot(
  fit_chain,
  "incidence",
  "analysing-results-chain-results-incidence",
  chain_results = TRUE,
  colours = c("#0072B2", "#D55E00", "#009E73")
)

message("Completed figure generation for vignette examples.")
