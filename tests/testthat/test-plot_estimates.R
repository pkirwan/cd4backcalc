test_that("migration plot titles describe their estimands", {
  expect_identical(
    get_plot_title("ratio_abroad_uk"),
    "born-abroad:UK-born ratio of UK-acquired incidence"
  )
  expect_identical(
    get_plot_title("diag_prob_mig"),
    "HIV diagnosis probabilities (born abroad), by CD4 strata"
  )
})

test_that("plot labels are inherited from model list names", {
  estimates <- tibble::tibble(
    year = 0:2,
    `2.5%` = 1:3,
    `50%` = 2:4,
    mean = 2:4,
    `97.5%` = 3:5
  )
  models <- list(
    "Spline incidence" = list(incidence = estimates),
    "Gaussian-process incidence" = list(incidence = estimates)
  )

  plot <- plot_estimates(models, quantity = "incidence")

  expect_identical(unique(plot$data$group), names(models))
})

test_that("simulated plots are generated correctly", {
  # Simulate diagnoses
  stan_data <- simulate_diagnoses(
    nquar = 40,
    sim_seed = 123,
    age = TRUE,
    migration = TRUE,
    age_diag_effect = 0.5,
    sim_type = "combo_3",
    gamma_obs_sd = 0.1,
    gamma_obs_scale_trend = 0.1,
    abroad_age_trend = 0.5,
    detect_step_quarter = 30,
    detect_step_size = 0.1
  )
  # Generate plots of key quantities
  p1 <- plot_simulations(stan_data, "incidence")
  p2 <- plot_simulations(stan_data, "undiag_prev")
  p3 <- plot_simulations(stan_data, "diag_prob_age")
  p4 <- plot_simulations(stan_data, "diagnoses")
  p5 <- plot_simulations(stan_data, "undiag_migration_age")
  p6 <- plot_simulations(stan_data, "migration_prob")
  p7 <- plot_simulations(stan_data, "transfers")
  p8 <- plot_simulations(stan_data, "transfers_detected")
  p9 <- plot_simulations(stan_data, "ratio_abroad_uk")
  p10 <- plot_simulations(stan_data, "detect_prob")
  p11 <- plot_simulations(stan_data, "migrant_hiv_prevalence")
  expect_s3_class(p1, "ggplot")
  expect_s3_class(p2, "ggplot")
  expect_s3_class(p3, "ggplot")
  expect_s3_class(p4, "ggplot")
  expect_s3_class(p5, "ggplot")
  expect_s3_class(p6, "ggplot")
  expect_s3_class(p7, "ggplot")
  expect_s3_class(p8, "ggplot")
  expect_s3_class(p9, "ggplot")
  expect_s3_class(p10, "ggplot")
  expect_s3_class(p11, "ggplot")
})

test_that("observed Stan quantities can be plotted", {
  sim_diags <- simulate_diagnoses(
    nquar = 40,
    sim_seed = 123,
    age = TRUE,
    migration = TRUE,
    sim_type = "combo_3",
    gamma_obs_sd = 0.1,
    gamma_obs_grouped = TRUE
  )

  quantities <- c(
    "aids_diagnoses",
    "aids_diagnoses_age",
    "transfers_detected_age",
    "gamma_obs",
    "gamma_obs_group",
    "cd4_diagnoses",
    "cd4_diagnoses_age",
    "recent_arrival_cd4"
  )
  plots <- lapply(quantities, function(quantity) {
    plot_simulations(sim_diags, quantity)
  })

  expect_true(all(vapply(plots, inherits, logical(1), what = "ggplot")))

  ai_sim_diags <- simulate_diagnoses(
    nquar = 40,
    sim_seed = 123,
    migration = TRUE
  )
  recent_arrival <- extract_simulated_quantity(
    ai_sim_diags,
    "recent_arrival_cd4"
  )
  expect_equal(nrow(recent_arrival), 40L)
  expect_equal(
    recent_arrival$value,
    as.vector(ai_sim_diags$stan_data$recent_arrival_cd4)
  )
  expect_s3_class(plot_simulations(ai_sim_diags, "recent_arrival_cd4"), "ggplot")
})

test_that("observed quantities can be plotted without sim_config", {
  sim_diags <- simulate_diagnoses(
    nquar = 40,
    sim_seed = 123,
    age = TRUE,
    migration = TRUE,
    gamma_obs_sd = 0.1,
    gamma_obs_grouped = TRUE
  )
  real_data <- sim_diags
  real_data$stan_data$sim_config <- NULL
  real_data$stan_data$model_config <- list(
    migration = TRUE,
    rita = FALSE,
    age = TRUE
  )
  real_data$sim_data <- NULL

  quantities <- c(
    "diagnoses",
    "diagnoses_age",
    "aids_diagnoses",
    "aids_diagnoses_age",
    "transfers_detected",
    "transfers_detected_age",
    "gamma_obs",
    "gamma_obs_group",
    "cd4_diagnoses",
    "cd4_diagnoses_age",
    "recent_arrival_cd4"
  )
  plots <- lapply(quantities, function(quantity) {
    plot_simulations(real_data, quantity)
  })

  expect_true(all(vapply(plots, inherits, logical(1), what = "ggplot")))
  expect_identical(
    plots[[1]]$labels$title,
    "Observed HIV diagnoses"
  )
})

test_that("simulated plots are generated correctly", {
  # Simulate diagnoses
  stan_data <- simulate_diagnoses(
    nquar = 40,
    sim_seed = 123,
    migration = TRUE,
    sim_type = "combo_4",
    eta_sd = 0.3,
    detect_obs_sd = 0.3,
    detect_obs_times = seq(1, 40, by = 20)
  )
  # Generate plots of key quantities
  p1 <- plot_simulations(stan_data, "incidence")
  p2 <- plot_simulations(stan_data, "undiag_prev")
  p3 <- plot_simulations(stan_data, "diag_prob")
  p4 <- plot_simulations(stan_data, "diagnoses")
  p5 <- plot_simulations(stan_data, "undiag_migration")
  p6 <- plot_simulations(stan_data, "migration_prob")
  p7 <- plot_simulations(stan_data, "transfers")
  p8 <- plot_simulations(stan_data, "transfers_detected")
  p9 <- plot_simulations(stan_data, "ratio_abroad_uk")
  p10 <- plot_simulations(stan_data, "detect_prob") +
    ggplot2::geom_point(
      data = data.frame(
        year = (stan_data$stan_data$detect_obs_index - 1) / 4,
        value = invlogit(stan_data$stan_data$detect_logit_obs),
        group = "detect_prob_obs"
      ),
      color = "red",
      size = 2
    )
  p11 <- plot_simulations(stan_data, "diag_prob_mig")
  expect_s3_class(p1, "ggplot")
  expect_s3_class(p2, "ggplot")
  expect_s3_class(p3, "ggplot")
  expect_s3_class(p4, "ggplot")
  expect_s3_class(p5, "ggplot")
  expect_s3_class(p6, "ggplot")
  expect_s3_class(p7, "ggplot")
  expect_s3_class(p8, "ggplot")
  expect_s3_class(p9, "ggplot")
  expect_s3_class(p10, "ggplot")
  expect_s3_class(p11, "ggplot")
})

test_that("simulated plots are generated correctly for migration_prob", {
  # Simulate diagnoses
  stan_data <- simulate_diagnoses(
    nquar = 40,
    sim_seed = 123,
    migration = TRUE,
    m_pattern = "varying"
  )
  # Generate plots of key quantities
  p1 <- plot_simulations(stan_data, "migration_prob")
  expect_s3_class(p1, "ggplot")
  stan_data <- simulate_diagnoses(
    nquar = 40,
    sim_seed = 123,
    migration = TRUE,
    rita = TRUE,
    m_pattern = "varying"
  )
  # Generate plots of key quantities
  p2 <- plot_simulations(stan_data, "migration_prob")
  expect_s3_class(p2, "ggplot")
})

test_that("simulated age plots are generated correctly", {
  # Simulate diagnoses
  stan_data <- simulate_diagnoses(
    nquar = 40,
    sim_seed = 123,
    age = TRUE,
    migration = TRUE,
    sim_type = "combo_3"
  )
  # Generate plots of key quantities
  p1 <- plot_simulations(stan_data, "incidence_age")
  p2 <- plot_simulations(stan_data, "undiag_prev_age")
  p3 <- plot_simulations(stan_data, "diag_prob_age")
  p4 <- plot_simulations(stan_data, "diag_prob_mig_age")
  expect_s3_class(p1, "ggplot")
  expect_s3_class(p2, "ggplot")
  expect_s3_class(p3, "ggplot")
  expect_s3_class(p4, "ggplot")
})

test_that("plots are generated correctly", {
  # Simulate diagnoses
  stan_data <- simulate_diagnoses(nquar = 40, sim_seed = 123)
  # Fit back-calculation model
  suppressWarnings({
    hiv_list <- run_backcalc(
      stan_data,
      model_seed = 123,
      iter_warmup = 10,
      iter_sampling = 10,
      chains = 1,
      refresh = 0,
      messages = FALSE,
      nworkers = 1
    )
  })
  # plot diagnostics
  p0 <- plot_diagnostics(hiv_list)
  expect_s3_class(p0, "ggplot")
  # Generate plots of key quantities
  p1 <- plot_estimates(hiv_list, "incidence")
  p2 <- plot_estimates(hiv_list, "undiag_prev")
  p3 <- plot_estimates(hiv_list, "diag_prob")
  p4 <- plot_estimates(hiv_list, "diagnoses")
  expect_s3_class(p1, "ggplot")
  expect_s3_class(p2, "ggplot")
  expect_s3_class(p3, "ggplot")
  expect_s3_class(p4, "ggplot")
})

test_that("dual plots are generated correctly", {
  # Simulate diagnoses
  stan_data <- simulate_diagnoses(nquar = 40, sim_seed = 123)
  # Fit back-calculation model
  suppressWarnings({
    hiv_list <- run_backcalc(
      stan_data,
      model_seed = 123,
      iter_warmup = 10,
      iter_sampling = 10,
      chains = 1,
      refresh = 0,
      messages = FALSE,
      nworkers = 1
    )
  })
  suppressWarnings({
    hiv_list_2 <- run_backcalc(
      stan_data,
      model_seed = 1234,
      iter_warmup = 10,
      iter_sampling = 10,
      chains = 1,
      refresh = 0,
      messages = FALSE,
      nworkers = 1
    )
  })
  # plot diagnostics
  p0 <- plot_diagnostics(hiv_list)
  expect_s3_class(p0, "ggplot")
  # Generate plots of key quantities
  p1 <- plot_estimates(
    hiv_list,
    "incidence",
    hiv_list_2 = hiv_list_2,
    model_1_name = "Model A",
    model_2_name = "Model B"
  )
  p2 <- plot_estimates(
    hiv_list,
    "undiag_prev",
    hiv_list_2 = hiv_list_2,
    model_1_name = "Model A",
    model_2_name = "Model B"
  )
  p3 <- plot_estimates(
    hiv_list,
    "diag_prob",
    hiv_list_2 = hiv_list_2,
    model_1_name = "Model A",
    model_2_name = "Model B"
  )
  p4 <- plot_estimates(
    hiv_list,
    "diagnoses",
    hiv_list_2 = hiv_list_2,
    model_1_name = "Model A",
    model_2_name = "Model B"
  )
  expect_s3_class(p1, "ggplot")
  expect_s3_class(p2, "ggplot")
  expect_s3_class(p3, "ggplot")
  expect_s3_class(p4, "ggplot")
})

test_that("chained plots are generated correctly", {
  # Simulate diagnoses
  stan_data <- simulate_diagnoses(nquar = 40, sim_seed = 123)
  # Fit back-calculation model
  suppressWarnings({
    hiv_list <- run_backcalc(
      stan_data,
      model_seed = 123,
      iter_warmup = 10,
      iter_sampling = 10,
      chains = 2,
      refresh = 0,
      messages = FALSE,
      nworkers = 1,
      chain_results = TRUE
    )
  })
  # Generate plots of key quantities
  p1 <- plot_estimates(hiv_list, "incidence", chain_results = TRUE)
  expect_s3_class(p1, "ggplot")
})

test_that("chain plot colours can be customised", {
  groups <- c("Model 1, Chain 1", "Model 1, Chain 2")

  expect_identical(
    cd4backcalc:::make_model_colors(
      groups,
      chain_results = TRUE,
      colours = c("navy", "firebrick")
    ),
    stats::setNames(c("navy", "firebrick"), groups)
  )
  expect_identical(
    cd4backcalc:::make_model_colors(
      groups,
      chain_results = TRUE,
      colours = c(
        "Model 1, Chain 2" = "firebrick",
        "Model 1, Chain 1" = "navy"
      )
    ),
    stats::setNames(c("navy", "firebrick"), groups)
  )
})

test_that("plots are generated correctly for migration model", {
  # Simulate diagnoses
  stan_data <- simulate_diagnoses(nquar = 60, sim_seed = 123, migration = TRUE)
  # Fit back-calculation model
  suppressWarnings({
    hiv_list <- run_backcalc(
      stan_data,
      model_seed = 123,
      iter_warmup = 10,
      iter_sampling = 10,
      chains = 1,
      refresh = 0,
      messages = FALSE,
      nworkers = 1
    )
  })
  # plot diagnostics
  p0 <- plot_diagnostics(hiv_list)
  expect_s3_class(p0, "ggplot")
  # Generate plots of key quantities
  p1 <- plot_estimates(hiv_list, "undiag_migration")
  p2 <- plot_estimates(hiv_list, "migration_prob")
  p3 <- plot_estimates(hiv_list, "transfers")
  p4 <- plot_estimates(hiv_list, "transfers_detected")
  p5 <- plot_estimates(hiv_list, "ratio_abroad_uk")
  p6 <- plot_estimates(hiv_list, "detect_prob")
  expect_s3_class(p1, "ggplot")
  expect_s3_class(p2, "ggplot")
  expect_s3_class(p3, "ggplot")
  expect_s3_class(p4, "ggplot")
  expect_s3_class(p5, "ggplot")
  expect_s3_class(p6, "ggplot")
})

test_that("plots work for migration model lists", {
  # Simulate diagnoses
  stan_data <- simulate_diagnoses(nquar = 60, sim_seed = 123, migration = TRUE)
  # Fit back-calculation model
  hiv_list <- list()
  for (i in 1:2) {
    suppressWarnings({
      hiv_list[[i]] <- run_backcalc(
        stan_data,
        model_seed = 123 + i,
        iter_warmup = 10,
        iter_sampling = 10,
        chains = 1,
        refresh = 0,
        messages = FALSE,
        nworkers = 1
      )
    })
  }
  # plot diagnostics
  p0 <- plot_diagnostics(hiv_list)
  expect_s3_class(p0, "ggplot")
  # Generate plots of key quantities
  p1 <- plot_estimates(hiv_list, "incidence")
  p2 <- plot_estimates(hiv_list, "undiag_prev")
  p3 <- plot_estimates(hiv_list, "diag_prob")
  p4 <- plot_estimates(hiv_list, "diagnoses")
  expect_s3_class(p1, "ggplot")
  expect_s3_class(p2, "ggplot")
  expect_s3_class(p3, "ggplot")
  expect_s3_class(p4, "ggplot")
  # Generate plots of key quantities
  p5 <- plot_estimates(hiv_list, "undiag_migration")
  p6 <- plot_estimates(hiv_list, "migration_prob")
  p7 <- plot_estimates(hiv_list, "transfers")
  p8 <- plot_estimates(hiv_list, "transfers_detected")
  p9 <- plot_estimates(hiv_list, "ratio_abroad_uk")
  p10 <- plot_estimates(hiv_list, "detect_prob")
  expect_s3_class(p5, "ggplot")
  expect_s3_class(p6, "ggplot")
  expect_s3_class(p7, "ggplot")
  expect_s3_class(p8, "ggplot")
  expect_s3_class(p9, "ggplot")
  expect_s3_class(p10, "ggplot")
})

test_that("plots work for rita model lists", {
  # Simulate diagnoses
  stan_data <- simulate_diagnoses(nquar = 40, sim_seed = 123, rita = TRUE)
  # Fit back-calculation model
  hiv_list <- list()
  for (i in 1:2) {
    suppressWarnings({
      hiv_list[[i]] <- run_backcalc(
        stan_data,
        model_seed = 123 + i,
        iter_warmup = 10,
        iter_sampling = 10,
        chains = 1,
        refresh = 0,
        messages = FALSE,
        nworkers = 1
      )
    })
  }
  # plot diagnostics
  p0 <- plot_diagnostics(hiv_list)
  expect_s3_class(p0, "ggplot")
  # Generate plots of key quantities
  p1 <- plot_estimates(hiv_list, "incidence")
  p2 <- plot_estimates(hiv_list, "undiag_prev")
  p3 <- plot_estimates(hiv_list, "diag_prob")
  p4 <- plot_estimates(hiv_list, "diagnoses")
  expect_s3_class(p1, "ggplot")
  expect_s3_class(p2, "ggplot")
  expect_s3_class(p3, "ggplot")
  expect_s3_class(p4, "ggplot")
})


test_that("plots work for migration and rita model bias plots", {
  # Simulate diagnoses
  stan_data <- simulate_diagnoses(
    nquar = 40,
    sim_seed = 123,
    rita = TRUE,
    migration = TRUE,
    m_pattern = "varying"
  )
  # Fit back-calculation model
  hiv_list <- list()
  for (i in 1:2) {
    suppressWarnings({
      hiv_list[[i]] <- run_backcalc(
        stan_data,
        model_seed = 123 + i,
        iter_warmup = 10,
        iter_sampling = 10,
        chains = 1,
        refresh = 0,
        messages = FALSE,
        nworkers = 1
      )
    })
  }
  # plot diagnostics
  p0 <- plot_diagnostics(hiv_list)
  expect_s3_class(p0, "ggplot")
  # Generate plots of key quantities
  p1 <- bias_plot(list(hiv_list), "incidence")
  p2 <- plot_estimates(hiv_list, "migration_prob")
  expect_s3_class(p1, "ggplot")
  expect_s3_class(p2, "ggplot")
})

test_that("plots work for age-dependent model", {
  # Simulate diagnoses
  stan_data <- simulate_diagnoses(nquar = 40, sim_seed = 123, age = TRUE)
  # Fit back-calculation model
  suppressWarnings({
    hiv_list <- run_backcalc(
      stan_data,
      model_seed = 123,
      iter_warmup = 5,
      iter_sampling = 5,
      chains = 1,
      refresh = 0,
      messages = FALSE,
      nworkers = 1
    )
  })
  # plot diagnostics
  p0 <- plot_diagnostics(hiv_list)
  expect_s3_class(p0, "ggplot")
  # Generate plots of key quantities
  p1 <- plot_estimates(hiv_list, "incidence")
  p2 <- plot_estimates(hiv_list, "undiag_prev")
  p3 <- plot_estimates(hiv_list, "diagnoses")
  p4 <- plot_estimates(hiv_list, "incidence_age")
  p5 <- plot_estimates(hiv_list, "undiag_prev_age")
  p6 <- plot_estimates(hiv_list, "diagnoses_age")
  p7 <- plot_estimates(hiv_list, "diag_prob_age")
  expect_s3_class(p1, "ggplot")
  expect_s3_class(p2, "ggplot")
  expect_s3_class(p3, "ggplot")
  expect_s3_class(p4, "ggplot")
  expect_s3_class(p5, "ggplot")
  expect_s3_class(p6, "ggplot")
  expect_s3_class(p7, "ggplot")
})

test_that("plots work for age-dependent rita model", {
  # Simulate diagnoses
  stan_data <- simulate_diagnoses(
    nquar = 40,
    sim_seed = 123,
    age = TRUE,
    rita = TRUE
  )
  # Fit back-calculation model
  suppressWarnings({
    hiv_list <- run_backcalc(
      stan_data,
      model_seed = 123,
      iter_warmup = 5,
      iter_sampling = 5,
      chains = 1,
      refresh = 0,
      messages = FALSE,
      nworkers = 1
    )
  })
  # plot diagnostics
  p0 <- plot_diagnostics(hiv_list)
  expect_s3_class(p0, "ggplot")
  # Generate plots of key quantities
  p1 <- plot_estimates(hiv_list, "incidence")
  p2 <- plot_estimates(hiv_list, "undiag_prev")
  p3 <- plot_estimates(hiv_list, "diagnoses")
  p4 <- plot_estimates(hiv_list, "incidence_age")
  p5 <- plot_estimates(hiv_list, "undiag_prev_age")
  p6 <- plot_estimates(hiv_list, "diagnoses_age")
  p7 <- plot_estimates(hiv_list, "diag_prob_age")
  expect_s3_class(p1, "ggplot")
  expect_s3_class(p2, "ggplot")
  expect_s3_class(p3, "ggplot")
  expect_s3_class(p4, "ggplot")
  expect_s3_class(p5, "ggplot")
  expect_s3_class(p6, "ggplot")
  expect_s3_class(p7, "ggplot")
})


test_that("start_year and overlay of real values works", {
  # Simulate diagnoses
  stan_data <- simulate_diagnoses(
    nquar = 40,
    sim_seed = 123,
    age = TRUE,
    rita = TRUE,
    migration = TRUE
  )
  # Fit back-calculation model
  suppressWarnings({
    hiv_list <- run_backcalc(
      stan_data,
      model_seed = 123,
      iter_warmup = 5,
      iter_sampling = 5,
      chains = 1,
      refresh = 0,
      messages = FALSE,
      nworkers = 1
    )
  })
  # Remove sim_data element from list to mimic real data scenario
  hiv_list <- hiv_list[names(hiv_list) != "sim_data"]
  # Generate plots of key quantities
  p1 <- plot_estimates(hiv_list, "diagnoses", start_year = 2000)
  p2 <- plot_estimates(hiv_list, "diagnoses_age", start_year = 2000)
  p3 <- plot_estimates(hiv_list, "transfers", start_year = 2000)
  p4 <- plot_estimates(hiv_list, "transfers_detected", start_year = 2000)
  p5 <- plot_estimates(hiv_list, "migrant_hiv_prevalence", start_year = 2000)
  expect_s3_class(p1, "ggplot")
  expect_s3_class(p2, "ggplot")
  expect_s3_class(p3, "ggplot")
  expect_s3_class(p4, "ggplot")
  expect_s3_class(p5, "ggplot")
})
