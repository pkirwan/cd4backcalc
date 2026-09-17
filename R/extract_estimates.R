#' Get simulation data for plotting.
#'
#' @param sim_list A list of simulation results.
#' @param var The variable to extract.
#' @param var2 The sub-variable to extract.
#' @param plot_years Number of years in the plot.
#' @return A data frame of simulation results for plotting.
#' @keywords internal
get_sims <- function(sim_list, var, var2, plot_years = 29.75) {
  sim_list <- purrr::compact(sim_list)
  if (length(sim_list) == 0) return(tibble::tibble())
  # make sim_list a list of lists if not already
  if (!("stan_data" %in% names(sim_list[[1]]))) {
    sim_list <- list(sim_list)
  }
  # extract estimates
  estimates <- purrr::map(seq_along(sim_list), function(i) {
    x <- sim_list[[i]]
    x[[var]][[var2]] |>
      tibble::as_tibble() |>
      dplyr::mutate(
        year = seq(0, plot_years, by = 0.25),
        group = paste0("Simulation ", i)
      )
  }) |>
    purrr::list_rbind()
  return(estimates)
}

#' Get estimates from HIV model results for plotting.
#'
#' @param hiv_list A list of HIV model results.
#' @param quantity The variable to extract.
#' @return A data frame of model estimates for plotting.
#' @examples
#' # See ?run_backcalc for a full end-to-end workflow.
#' invisible(NULL)
#' @export
get_estimates <- function(hiv_list, quantity) {
  hiv_list <- purrr::compact(hiv_list)
  if (length(hiv_list) == 0) return(tibble::tibble())
  # make hiv_list a list of lists if not already
  if (!("incidence" %in% names(hiv_list[[1]])) && !("incidence_1" %in% names(hiv_list[[1]]))) {
    hiv_list <- list(hiv_list)
  }
  # extract estimates
  estimates <- purrr::map(seq_along(hiv_list), function(i) {
    x <- hiv_list[[i]]
    list_name <- names(hiv_list)[i]
    group_label <- if (!is.null(x$group)) {
      x$group
    } else if (!is.null(list_name) && nzchar(list_name)) {
      list_name
    } else {
      paste0("Model ", i)
    }

    estimate_tbl <- x[[quantity]] |>
      tibble::as_tibble()

    if (all(c("quartile", "value") %in% names(estimate_tbl))) {
      estimate_tbl <- estimate_tbl |>
        tidyr::pivot_wider(names_from = "quartile", values_from = "value")
    } else {
      quantile_cols <- intersect(c("2.5%", "50%", "mean", "97.5%"), names(estimate_tbl))
      if (length(quantile_cols) == 0) {
        stop(
          sprintf(
            "Estimate '%s' must contain either quartile/value columns or summary columns.",
            quantity
          )
        )
      }
    }

    estimate_tbl |>
      dplyr::mutate(group = group_label, iteration = i)
  }) |>
    purrr::list_rbind()

  # Remap age_band labels to match simulated quantity labels
  if ("age_band" %in% names(estimates)) {
    age_band_map <- c(
      "1524" = "15-24",
      "2534" = "25-34",
      "3544" = "35-44",
      "45" = "45+"
    )
    estimates <- estimates |>
      dplyr::mutate(
        age_band = dplyr::if_else(
          as.character(.data$age_band) %in% names(age_band_map),
          age_band_map[as.character(.data$age_band)],
          as.character(.data$age_band)
        )
      )
  }

  return(estimates)
}

make_grouped_age_data <- function(values, group = "Generative values") {
  age_group_labels <- c("15-24", "25-34", "35-44", "45+")
  nquar <- nrow(values)

  tibble::tibble(
    value = as.vector(values),
    year = rep(seq(0, nquar / 4 - 0.25, by = 0.25), 4),
    age_band = factor(
      rep(age_group_labels, each = nquar),
      levels = age_group_labels
    ),
    group = group
  )
}

make_age_band_data <- function(values, group = "Generative values") {
  age_group <- pmin(4L, 1L + (seq_len(ncol(values)) - 1L) %/% 10L)
  grouped_values <- vapply(
    1:4,
    function(i) rowSums(values[, age_group == i, drop = FALSE]),
    numeric(nrow(values))
  )

  make_grouped_age_data(grouped_values, group)
}

cd4_state_labels <- function(rita) {
  labels <- c("CD4>=500", "CD4 350-499", "CD4 200-349", "CD4<200")
  if (rita) c("Recently acquired", labels) else labels
}

make_cd4_data <- function(values, migration, rita, group) {
  if (migration) {
    values <- apply(values, c(1, length(dim(values))), sum)
  }
  state_labels <- cd4_state_labels(rita)
  colnames(values) <- state_labels

  tibble::as_tibble(values) |>
    dplyr::mutate(
      year = seq(0, nrow(values) / 4 - 0.25, by = 0.25),
      group = group
    ) |>
    tidyr::pivot_longer(
      cols = tidyselect::all_of(state_labels),
      names_to = "strata",
      values_to = "value"
    )
}

make_cd4_age_data <- function(values, migration, rita, group) {
  if (migration) {
    values <- apply(values, c(1, 2, 4), sum)
  }
  state_labels <- cd4_state_labels(rita)

  purrr::map(seq_along(state_labels), function(i) {
    make_age_band_data(values[,, i], group = group) |>
      dplyr::mutate(
        strata = factor(state_labels[i], levels = state_labels)
      )
  }) |>
    purrr::list_rbind()
}

make_recent_arrival_cd4_data <- function(stan_data, age, group) {
  state_labels <- cd4_state_labels(FALSE)
  age_bands <- if (age) c("15-24", "25-34", "35-44", "45+") else "all"
  nyear <- if (is.null(stan_data$nyear)) {
    stan_data$nquar %/% 4
  } else {
    stan_data$nyear
  }

  tidyr::expand_grid(
    strata = factor(state_labels, levels = state_labels),
    age_band = factor(age_bands, levels = age_bands),
    year = seq_len(nyear) - 1
  ) |>
    dplyr::mutate(
      value = as.vector(stan_data$recent_arrival_cd4),
      group = group
    )
}

#' Extract simulated quantities based on the provided configuration.
#'
#' @param sim_diags The simulated diagnosis data containing stan_data and sim_data.
#' @param quantity The quantity to extract.
#' @param deterministic Logical indicating whether to re-run the simulation deterministically (dist = FALSE) to obtain the expected values (lambda) rather than the stochastic realisation. Useful for bias/MSE calculations across multiple seeds (default is FALSE).
#' @return A list of model estimates including undiagnosed prevalence,
#'   incidence, migration, diagnosis probabilities, diagnoses, migration
#'   probabilities, the born-abroad:UK-born ratio of UK-acquired incidence,
#'   transfers of care, CD4 strata, and the HIV-positive share of the
#'   external migration index (migrant_hiv_prevalence).
#' @examples
#' # See ?run_backcalc for a full end-to-end workflow.
#' invisible(NULL)
#' @export
extract_simulated_quantity <- function(sim_diags, quantity, deterministic = FALSE) {
  if (deterministic) {
    sim_config <- sim_diags$stan_data$sim_config
    sim_diags <- simulate_diagnoses(
      h_pattern = sim_config$h_pattern,
      o_pattern = sim_config$o_pattern,
      d_pattern = sim_config$d_pattern,
      p_pattern = sim_config$p_pattern,
      m_pattern = sim_config$m_pattern,
      r_pattern = sim_config$r_pattern,
      a_pattern = sim_config$a_pattern,
      sim_seed = sim_config$sim_seed,
      dist = FALSE,
      migration = sim_config$migration,
      rita = sim_config$rita,
      age = sim_config$age,
      nquar = sim_config$nquar,
      nage = sim_config$nage,
      q = sim_config$q,
      target_inf = sim_config$target_inf,
      target_mig = sim_config$target_mig,
      abroad_prop = sim_config$abroad_prop,
      cd4_completeness = 1,
      rita_completeness = 1,
      age_diag_effect = sim_config$age_diag_effect,
      mig_diag_offset = sim_config$mig_diag_offset,
      frr = sim_config$frr,
      gamma_obs_grouped = sim_config$gamma_obs_grouped,
      gamma_obs_sd = sim_config$gamma_obs_sd,
      gamma_obs_scale = sim_config$gamma_obs_scale,
      abroad_age_trend = sim_config$abroad_age_trend,
      detect_step_quarter = sim_config$detect_step_quarter,
      detect_step_size = sim_config$detect_step_size,
      gamma_obs_scale_trend = sim_config$gamma_obs_scale_trend,
      gamma_obs_drift_sd = sim_config$gamma_obs_drift_sd,
      recent_arrival_completeness = sim_config$recent_arrival_completeness
    )
  }
  stan_data <- sim_diags$stan_data
  sim_data <- sim_diags$sim_data

  plot_years <- stan_data$nquar / 4 - 0.25
  age_group_labels <- c("15-24", "25-34", "35-44", "45+")

  # Create data frames for each basic quantity
  undiag_prev <- if (stan_data$sim_config$migration) {
    tibble::tibble(
      value = rowSums(sim_data$epi_uk$lat) +
        rowSums(sim_data$epi_uk_abroad$lat) +
        rowSums(sim_data$epi_abroad$lat),
      year = seq(0, plot_years, by = 0.25),
      group = "Generative values"
    )
  } else {
    tibble::tibble(
      value = rowSums(sim_data$epi_uk$lat) +
        rowSums(sim_data$epi_uk_abroad$lat),
      year = seq(0, plot_years, by = 0.25),
      group = "Generative values"
    )
  }
  incidence <- tibble::tibble(
    value = sim_data$infs[, 1],
    year = seq(0, plot_years, by = 0.25),
    group = "Generative values"
  )
  undiag_migration <- tibble::tibble(
    # m[, 1:nlatst] now contains only true undiagnosed migrants (misclassified
    # are handled separately via detect_prob split in sim_diags)
    value = rowSums(sim_data$m[, 1:stan_data$sim_config$nlatst]),
    year = seq(0, plot_years, by = 0.25),
    group = "Generative values"
  )
  all_migration <- tibble::tibble(
    value = sim_data$omega,
    year = seq(0, plot_years, by = 0.25),
    group = "Generative values"
  )
  total_diags <- sim_data$epi_uk$diags +
    sim_data$epi_uk_abroad$diags +
    sim_data$epi_abroad$diags
  total_diags_age <- sim_data$epi_uk$diags_age +
    sim_data$epi_uk_abroad$diags_age +
    sim_data$epi_abroad$diags_age
  diagnoses <- tibble::tibble(
    value = rowSums(total_diags) + sim_data$misclassified,
    year = seq(0, plot_years, by = 0.25),
    group = "Generative values"
  )
  diag_prob <- if (stan_data$sim_config$rita) {
    # Use middle age as representative for non-age-stratified summary
    mid_age <- ceiling(stan_data$nage / 2)
    tibble::tibble(
      value = as.vector(sim_data$d_probs[, mid_age, ]),
      year = rep(seq(0, plot_years, by = 0.25), 5),
      strata = rep(
        c(
          "Recently acquired",
          "CD4>=500",
          "CD4 350-499",
          "CD4 200-349",
          "CD4<200"
        ),
        each = stan_data$nquar
      ),
      group = "Generative values"
    ) |>
      dplyr::mutate(
        strata = factor(
          .data$strata,
          levels = c(
            "Recently acquired",
            "CD4>=500",
            "CD4 350-499",
            "CD4 200-349",
            "CD4<200"
          )
        )
      )
  } else {
    mid_age <- ceiling(stan_data$nage / 2)
    tibble::tibble(
      value = as.vector(sim_data$d_probs[, mid_age, ]),
      year = rep(seq(0, plot_years, by = 0.25), 4),
      strata = rep(
        c("CD4>=500", "CD4 350-499", "CD4 200-349", "CD4<200"),
        each = stan_data$nquar
      ),
      group = "Generative values"
    ) |>
      dplyr::mutate(
        strata = factor(
          .data$strata,
          levels = c("CD4>=500", "CD4 350-499", "CD4 200-349", "CD4<200")
        )
      )
  }

  if (quantity == "undiag_prev") {
    out <- undiag_prev
  } else if (quantity == "incidence") {
    out <- incidence
  } else if (quantity == "undiag_migration") {
    out <- undiag_migration
  } else if (quantity == "all_migration") {
    out <- all_migration
  } else if (quantity == "diagnoses") {
    out <- diagnoses
  } else if (quantity == "aids_diagnoses") {
    out <- tibble::tibble(
      value = total_diags[, stan_data$sim_config$ndiagst],
      year = seq(0, plot_years, by = 0.25),
      group = "Generative values"
    )
  } else if (quantity == "diag_prob") {
    out <- diag_prob
  } else if (quantity == "diag_prob_mig") {
    # Born-abroad diagnosis probabilities mirror diag_prob but use d_probs_mig
    mid_age <- ceiling(stan_data$nage / 2)
    if (stan_data$sim_config$rita) {
      out <- tibble::tibble(
        value = as.vector(sim_data$d_probs_mig[, mid_age, ]),
        year = rep(seq(0, plot_years, by = 0.25), 5),
        strata = rep(
          c("Recently acquired", "CD4>=500", "CD4 350-499", "CD4 200-349", "CD4<200"),
          each = stan_data$nquar
        ),
        group = "Generative values"
      ) |>
        dplyr::mutate(
          strata = factor(
            .data$strata,
            levels = c("Recently acquired", "CD4>=500", "CD4 350-499", "CD4 200-349", "CD4<200")
          )
        )
    } else {
      out <- tibble::tibble(
        value = as.vector(sim_data$d_probs_mig[, mid_age, ]),
        year = rep(seq(0, plot_years, by = 0.25), 4),
        strata = rep(
          c("CD4>=500", "CD4 350-499", "CD4 200-349", "CD4<200"),
          each = stan_data$nquar
        ),
        group = "Generative values"
      ) |>
        dplyr::mutate(
          strata = factor(
            .data$strata,
            levels = c("CD4>=500", "CD4 350-499", "CD4 200-349", "CD4<200")
          )
        )
    }
  } else if (quantity == "undiag_prev_age") {
    out <- tidyr::expand_grid(
      year = seq(0, plot_years, by = 0.25),
      age = 15:(stan_data$nage + 14)
    ) |>
      dplyr::mutate(
        value = as.vector(apply(
          sim_data$epi_uk$lat_age +
            sim_data$epi_uk_abroad$lat_age +
            if (stan_data$sim_config$migration) {
              sim_data$epi_abroad$lat_age
            } else {
              0
            },
          c(1, 2),
          sum
        )),
        group = "Generative values",
        age_band = cut(
          .data$age,
          breaks = c(15, 25, 35, 45, 100),
          labels = age_group_labels,
          right = FALSE
        )
      ) |>
      dplyr::group_by(.data$year, .data$age_band, .data$group) |>
      dplyr::summarise(value = sum(.data$value), .groups = "drop") |>
      dplyr::bind_rows(
        undiag_prev |>
          dplyr::mutate(
            age_band = "all"
          )
      )
  } else if (quantity == "incidence_age") {
    out <- tidyr::expand_grid(
      year = seq(0, plot_years, by = 0.25),
      age = 15:(stan_data$nage + 14)
    ) |>
      dplyr::mutate(
        value = as.vector(sim_data$infs_age[,, 1]),
        group = "Generative values",
        age_band = cut(
          .data$age,
          breaks = c(15, 25, 35, 45, 100),
          labels = age_group_labels,
          right = FALSE
        )
      ) |>
      dplyr::group_by(.data$year, .data$age_band, .data$group) |>
      dplyr::summarise(value = sum(.data$value), .groups = "drop") |>
      dplyr::bind_rows(
        incidence |>
          dplyr::mutate(
            age_band = "all"
          )
      )
  } else if (quantity == "undiag_migration_age") {
    out <- tidyr::expand_grid(
      year = seq(0, plot_years, by = 0.25),
      age = 15:(stan_data$nage + 14)
    ) |>
      dplyr::mutate(
        value = as.vector(
          apply(
            sim_data$m_age[,,
              1:stan_data$sim_config$nlatst
            ],
            c(1, 2),
            sum
          )
        ),
        group = "Generative values",
        age_band = cut(
          .data$age,
          breaks = c(15, 25, 35, 45, 100),
          labels = age_group_labels,
          right = FALSE
        )
      ) |>
      dplyr::group_by(.data$year, .data$age_band, .data$group) |>
      dplyr::summarise(value = sum(.data$value), .groups = "drop") |>
      dplyr::bind_rows(
        undiag_migration |>
          dplyr::mutate(
            age_band = "all"
          )
      )
  } else if (quantity == "all_migration_age") {
    out <- tidyr::expand_grid(
      year = seq(0, plot_years, by = 0.25),
      age = 15:(stan_data$nage + 14)
    ) |>
      dplyr::mutate(
        value = as.vector(sim_data$omega_age),
        group = "Generative values",
        age_band = cut(
          .data$age,
          breaks = c(15, 25, 35, 45, 100),
          labels = age_group_labels,
          right = FALSE
        )
      ) |>
      dplyr::group_by(.data$year, .data$age_band, .data$group) |>
      dplyr::summarise(value = sum(.data$value), .groups = "drop") |>
      dplyr::bind_rows(
        undiag_migration |>
          dplyr::mutate(
            age_band = "all"
          )
      )
  } else if (quantity %in% c("diag_prob_age", "diag_prob_mig_age")) {
    d_array <- if (quantity == "diag_prob_mig_age") sim_data$d_probs_mig else sim_data$d_probs
    state_labels <- if (stan_data$sim_config$rita) {
      c("Recently acquired", "CD4>=500", "CD4 350-499", "CD4 200-349", "CD4<200")
    } else {
      c("CD4>=500", "CD4 350-499", "CD4 200-349", "CD4<200")
    }
    # d_probs is constant within each age group, so pick the first age per group
    rep_ages <- c(1, 11, 21, 31) # first age index in each group
    out <- tidyr::expand_grid(
      year = seq(0, plot_years, by = 0.25),
      age_band = factor(age_group_labels, levels = age_group_labels),
      strata = factor(state_labels, levels = state_labels)
    ) |>
      dplyr::mutate(
        value = as.vector(aperm(d_array[, rep_ages, , drop = FALSE], c(3, 2, 1))),
        group = "Generative values"
      )
  } else if (quantity == "diagnoses_age") {
    out <- make_age_band_data(
      t(apply(total_diags_age, c(1, 2), sum) + sim_data$misclassified_age)
    ) |>
      dplyr::bind_rows(
        diagnoses |>
          dplyr::mutate(
            age_band = "all"
          )
      )
  } else if (quantity == "aids_diagnoses_age") {
    out <- make_age_band_data(
      t(total_diags_age[,, stan_data$sim_config$ndiagst])
    )
  } else if (quantity == "transfers_detected_age") {
    out <- make_grouped_age_data(stan_data$ABROAD_group)
  } else if (quantity == "gamma_obs") {
    out <- tibble::tibble(
      value = exp(stan_data$gamma_obs),
      year = seq(0, plot_years, by = 0.25),
      group = "Generative values"
    )
  } else if (quantity == "gamma_obs_group") {
    gamma_obs_group <- stan_data$gamma_obs_group
    gamma_obs_group[gamma_obs_group >= 0] <- exp(gamma_obs_group[gamma_obs_group >= 0])
    out <- make_grouped_age_data(gamma_obs_group) |>
      dplyr::filter(.data$value >= 0)
  } else if (quantity == "cd4_diagnoses_age") {
    out <- make_cd4_age_data(
      stan_data$CD4_age,
      migration = stan_data$sim_config$migration,
      rita = stan_data$sim_config$rita,
      group = "Generative values"
    )
  } else if (quantity == "recent_arrival_cd4") {
    out <- make_recent_arrival_cd4_data(
      stan_data,
      age = stan_data$sim_config$age,
      group = "Generative values"
    )
  } else if (quantity == "migration_prob") {
    # no recently acquired category in migration prob
    out <- tibble::tibble(
      value = if (stan_data$sim_config$rita) {
        as.vector(sim_data$mu[, 4:ncol(sim_data$mu)])
      } else {
        as.vector(sim_data$mu)
      },
      year = rep(seq(0, plot_years, by = 0.25), 5),
      strata = rep(
        c(
          "CD4>=500",
          "CD4 350-499",
          "CD4 200-349",
          "CD4<200",
          "Transfers of care"
        ),
        each = stan_data$nquar
      ),
      group = "Generative values"
    ) |>
      dplyr::mutate(
        value = dplyr::if_else(is.na(.data$value), 0, .data$value),
        strata = factor(
          .data$strata,
          levels = c(
            "CD4>=500",
            "CD4 350-499",
            "CD4 200-349",
            "CD4<200",
            "Transfers of care"
          )
        )
      )
  } else if (quantity == "ratio_abroad_uk") {
    out <- tibble::tibble(
      value = sim_data$ratio_abroad_uk,
      year = seq(0, plot_years, by = 0.25),
      group = "Generative values"
    )
  } else if (quantity == "transfers") {
    out <- tibble::tibble(
      value = stan_data$ABROAD * (1 / sim_data$detect_prob),
      year = seq(0, plot_years, by = 0.25),
      group = "Generative values"
    )
  } else if (quantity == "transfers_detected") {
    out <- tibble::tibble(
      value = stan_data$ABROAD,
      year = seq(0, plot_years, by = 0.25),
      group = "Generative values"
    )
  } else if (quantity %in% c("cd4_strata", "cd4_diagnoses")) {
    out <- make_cd4_data(
      stan_data$CD4,
      migration = stan_data$sim_config$migration,
      rita = stan_data$sim_config$rita,
      group = "Generative values"
    )
  } else if (quantity == "detect_prob") {
    out <- tibble::tibble(
      value = sim_data$detect_prob,
      year = seq(0, plot_years, by = 0.25),
      group = "Generative values"
    )
  } else if (quantity == "migrant_hiv_prevalence") {
    out <- tibble::tibble(
      value = sim_data$migrant_hiv_prevalence,
      year = seq(0, plot_years, by = 0.25),
      group = "Generative values"
    )
  }
  return(out)
}

#' Extract real data for plotting
#'
#' @param hiv_list A list of HIV model results containing `stan_data`.
#' @param quantity A string specifying the quantity to extract.
#' @return A data frame of observed values for plotting, or `NULL` when the
#'   requested quantity is not present in the real-data inputs.
#' @examples
#' # See ?run_backcalc for a full end-to-end workflow.
#' invisible(NULL)
#' @export
extract_real_data <- function(hiv_list, quantity) {
  stan_data <- hiv_list$stan_data
  plot_years <- stan_data$nquar / 4 - 0.25
  model_config <- stan_data$model_config
  if (is.null(model_config)) {
    stop("stan_data$model_config is required for real-data plotting.")
  }
  migration <- isTRUE(model_config$migration)
  rita <- isTRUE(model_config$rita)
  age <- isTRUE(model_config$age)

  diagnoses <- if (migration) {
    tibble::tibble(
      value = rowSums(stan_data$HIV) + rowSums(stan_data$AIDS),
      year = seq(0, plot_years, by = 0.25),
      group = "Observed values"
    )
  } else {
    tibble::tibble(
      value = stan_data$HIV + stan_data$AIDS,
      year = seq(0, plot_years, by = 0.25),
      group = "Observed values"
    )
  }
  if (quantity == "diagnoses") {
    return(diagnoses)
  } else if (quantity == "aids_diagnoses") {
    return(tibble::tibble(
      value = if (migration) rowSums(stan_data$AIDS) else stan_data$AIDS,
      year = seq(0, plot_years, by = 0.25),
      group = "Observed values"
    ))
  } else if (quantity == "diagnoses_age") {
    if (is.null(stan_data$HIV_age) || is.null(stan_data$AIDS_age)) {
      return(NULL)
    }
    out <- make_age_band_data(
      apply(stan_data$HIV_age, c(1, 2), sum) +
        apply(stan_data$AIDS_age, c(1, 2), sum),
      group = "Observed values"
    ) |>
      dplyr::bind_rows(
        diagnoses |>
          dplyr::mutate(
            age_band = "all"
          )
      )
  } else if (quantity == "aids_diagnoses_age") {
    if (is.null(stan_data$AIDS_age)) {
      return(NULL)
    }
    out <- make_age_band_data(
      apply(stan_data$AIDS_age, c(1, 2), sum),
      group = "Observed values"
    )
  } else if (quantity == "transfers_detected") {
    if (is.null(stan_data$ABROAD)) {
      return(NULL)
    }
    # ABROAD counts only transfers whose previous diagnosis was identified, so
    # they are comparable with transfers_detected, not with total transfers.
    out <- tibble::tibble(
      value = stan_data$ABROAD,
      year = seq(0, plot_years, by = 0.25),
      group = "Observed values"
    )
  } else if (quantity == "transfers_detected_age") {
    if (!age || is.null(stan_data$ABROAD_group)) {
      return(NULL)
    }
    out <- make_grouped_age_data(
      stan_data$ABROAD_group,
      group = "Observed values"
    )
  } else if (quantity == "gamma_obs") {
    if (!migration || is.null(stan_data$gamma_obs)) {
      return(NULL)
    }
    out <- tibble::tibble(
      value = exp(stan_data$gamma_obs),
      year = seq(0, plot_years, by = 0.25),
      group = "Observed values"
    )
  } else if (quantity == "gamma_obs_group") {
    if (!migration || is.null(stan_data$gamma_obs_group)) {
      return(NULL)
    }
    gamma_obs_group <- stan_data$gamma_obs_group
    gamma_obs_group[gamma_obs_group >= 0] <-
      exp(gamma_obs_group[gamma_obs_group >= 0])
    out <- make_grouped_age_data(
      gamma_obs_group,
      group = "Observed values"
    ) |>
      dplyr::filter(.data$value >= 0)
  } else if (quantity == "cd4_diagnoses") {
    if (is.null(stan_data$CD4)) {
      return(NULL)
    }
    out <- make_cd4_data(
      stan_data$CD4,
      migration = migration,
      rita = rita,
      group = "Observed values"
    )
  } else if (quantity == "cd4_diagnoses_age") {
    if (is.null(stan_data$CD4_age)) {
      return(NULL)
    }
    out <- make_cd4_age_data(
      stan_data$CD4_age,
      migration = migration,
      rita = rita,
      group = "Observed values"
    )
  } else if (quantity == "recent_arrival_cd4") {
    if (!migration || is.null(stan_data$recent_arrival_cd4)) {
      return(NULL)
    }
    out <- make_recent_arrival_cd4_data(
      stan_data,
      age = age,
      group = "Observed values"
    )
  } else if (quantity == "detect_prob") {
    if (
      is.null(stan_data$detect_obs_index) ||
        is.null(stan_data$detect_logit_obs) ||
        length(stan_data$detect_obs_index) == 0
    ) {
      return(NULL)
    }
    out <- tibble::tibble(
      value = invlogit(stan_data$detect_logit_obs),
      year = (stan_data$detect_obs_index - 1) / 4,
      group = "Observed values"
    )
  }
  return(out)
}
