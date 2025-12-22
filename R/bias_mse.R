#' Get simulation data for plotting.
#'
#' @param sim_list A list of simulation results.
#' @param var The variable to extract.
#' @param var2 The sub-variable to extract.
#' @param plot_years Number of years in the plot.
#' @return A data frame of simulation results for plotting.
#' @export
get_sims <- function(sim_list, var, var2, plot_years = 29.75) {
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
#' @param var The variable to extract.
#' @return A data frame of model estimates for plotting.
#' @export
get_estimates <- function(hiv_list, var) {
  # make hiv_list a list of lists if not already
  if (!("incidence" %in% names(hiv_list[[1]]))) {
    hiv_list <- list(hiv_list)
  }
  # extract estimates
  estimates <- purrr::map(seq_along(hiv_list), function(i) {
    x <- hiv_list[[i]]
    x[[var]] |>
      tibble::as_tibble() |>
      dplyr::mutate(group = paste0("Model ", i)) |>
      tidyr::pivot_wider(names_from = "quartile", values_from = "value")
  }) |>
    purrr::list_rbind()
  return(estimates)
}

#' Extract simulated quantities based on the provided configuration.
#'
#' @param sim_diags The simulated diagnosis data containing stan_data and sim_data.
#' @return A list of model estimates including undiagnosed prevalence, incidence, migration, diagnosis probabilities, diagnoses, migration probabilities, proportion of UK-acquired infections that are migrants, transfers of care, and CD4 strata.
#' @export
extract_simulated_quantities <- function(sim_diags) {
  stan_data <- sim_diags$stan_data
  sim_data <- sim_diags$sim_data

  plot_years <- stan_data$nquar / 4 - 0.25
  age_group_labels <- c("1524", "2534", "3544", "45")

  undiag_prev <- if (stan_data$config$migration) {
    tibble::tibble(
      value = rowSums(sim_data$epi_uk$lat) +
        rowSums(sim_data$epi_uk_abroad$lat) +
        rowSums(sim_data$epi_abroad$lat),
      year = seq(0, plot_years, by = 0.25),
      group = "Simulated values"
    )
  } else {
    tibble::tibble(
      value = rowSums(sim_data$epi_uk$lat) +
        rowSums(sim_data$epi_uk_abroad$lat),
      year = seq(0, plot_years, by = 0.25),
      group = "Simulated values"
    )
  }

  undiag_prev_age <- tidyr::expand_grid(
    year = seq(0, plot_years, by = 0.25),
    age = 15:(stan_data$nage + 14)
  ) |>
    dplyr::mutate(
      value = as.vector(apply(
        sim_data$epi_uk$lat_age +
          sim_data$epi_uk_abroad$lat_age +
          if (stan_data$config$migration) {
            sim_data$epi_abroad$lat_age
          } else {
            0
          },
        c(1, 2),
        sum
      )),
      group = "Simulated values",
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

  incidence <- tibble::tibble(
    value = sim_data$infs[, 1],
    year = seq(0, plot_years, by = 0.25),
    group = "Simulated values"
  )

  incidence_age <- tidyr::expand_grid(
    year = seq(0, plot_years, by = 0.25),
    age = 15:(stan_data$nage + 14)
  ) |>
    dplyr::mutate(
      value = as.vector(sim_data$infs_age[,, 1]),
      group = "Simulated values",
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

  migration <- tibble::tibble(
    value = rowSums(sim_data$m[, 1:stan_data$config$nlatst]),
    year = seq(0, plot_years, by = 0.25),
    group = "Simulated values"
  )

  migration_age <- tidyr::expand_grid(
    year = seq(0, plot_years, by = 0.25),
    age = 15:(stan_data$nage + 14)
  ) |>
    dplyr::mutate(
      value = as.vector(apply(
        sim_data$m_age[,,
          1:stan_data$config$nlatst
        ],
        c(1, 2),
        sum
      )),
      group = "Simulated values",
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
      migration |>
        dplyr::mutate(
          age_band = "all"
        )
    )

  diag_prob <- if (stan_data$config$rita) {
    tibble::tibble(
      value = as.vector(sim_data$d_probs),
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
      group = "Simulated values"
    ) |>
      dplyr::mutate(
        strata = factor(
          .data$strata,
          levels = c(
            "CD4<200",
            "CD4 200-349",
            "CD4 350-499",
            "CD4>=500",
            "Recently acquired"
          )
        )
      )
  } else {
    tibble::tibble(
      value = as.vector(sim_data$d_probs),
      year = rep(seq(0, plot_years, by = 0.25), 4),
      strata = rep(
        c("CD4>=500", "CD4 350-499", "CD4 200-349", "CD4<200"),
        each = stan_data$nquar
      ),
      group = "Simulated values"
    ) |>
      dplyr::mutate(
        strata = factor(
          .data$strata,
          levels = c("CD4<200", "CD4 200-349", "CD4 350-499", "CD4>=500")
        )
      )
  }

  diag_prob_age <- tidyr::crossing(
    diag_prob,
    age_band = age_group_labels
  )

  diagnoses <- if (stan_data$config$migration) {
    tibble::tibble(
      value = rowSums(stan_data$HIV) + rowSums(stan_data$AIDS),
      year = seq(0, plot_years, by = 0.25),
      group = "Simulated values"
    )
  } else {
    tibble::tibble(
      value = stan_data$HIV + stan_data$AIDS,
      year = seq(0, plot_years, by = 0.25),
      group = "Simulated values"
    )
  }

  diagnoses_age <- tidyr::expand_grid(
    # diagnoses are age x year
    age = 15:(stan_data$nage + 14),
    year = seq(0, plot_years, by = 0.25)
  ) |>
    dplyr::mutate(
      value = as.vector(
        apply(stan_data$HIV_age, c(1, 2), sum) +
          apply(stan_data$AIDS_age, c(1, 2), sum)
      ),
      group = "Simulated values",
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
      diagnoses |>
        dplyr::mutate(
          age_band = "all"
        )
    )

  migration_prob <- if (stan_data$config$rita) {
    tibble::tibble(
      value = as.vector(sweep(
        sim_data$m,
        1,
        ifelse(sim_data$omega == 0, NA_real_, sim_data$omega),
        "/"
      )),
      year = rep(seq(0, plot_years, by = 0.25), 8),
      strata = rep(
        c(
          "Recently acquired 0-3, CD4>500",
          "Recently acquired 3-6, CD4 >=500",
          "Recently acquired 0-3, CD4 350-499",
          "CD4>=500",
          "CD4 350-499",
          "CD4 200-349",
          "CD4<200",
          "Transfers of care"
        ),
        each = stan_data$nquar
      ),
      group = "Simulated values"
    ) |>
      dplyr::mutate(
        value = dplyr::if_else(is.na(.data$value), 0, .data$value),
        strata = factor(
          .data$strata,
          levels = c(
            "CD4<200",
            "CD4 200-349",
            "CD4 350-499",
            "CD4>=500",
            "Recently acquired 0-3, CD4 350-499",
            "Recently acquired 3-6, CD4 >=500",
            "Recently acquired 0-3, CD4>=500",
            "Transfers of care"
          )
        )
      )
  } else {
    tibble::tibble(
      value = as.vector(sweep(
        sim_data$m,
        1,
        ifelse(sim_data$omega == 0, NA_real_, sim_data$omega),
        "/"
      )),
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
      group = "Simulated values"
    ) |>
      dplyr::mutate(
        value = dplyr::if_else(is.na(.data$value), 0, .data$value),
        strata = factor(
          .data$strata,
          levels = c(
            "CD4<200",
            "CD4 200-349",
            "CD4 350-499",
            "CD4>=500",
            "Transfers of care"
          )
        )
      )
  }

  ratio_abroad_uk <- tibble::tibble(
    value = sim_data$ratio_abroad_uk,
    year = seq(0, plot_years, by = 0.25),
    group = "Simulated values"
  )

  transfers <- tibble::tibble(
    value = stan_data$ABROAD,
    year = seq(0, plot_years, by = 0.25),
    group = "Simulated values"
  )

  cd4_strata <- if (stan_data$config$rita) {
    if (stan_data$config$migration) {
      tibble::tibble(
        `Recently acquired` = rowSums(stan_data$CD4[,, 1]),
        `CD4>=500` = rowSums(stan_data$CD4[,, 2]),
        `CD4 350-499` = rowSums(stan_data$CD4[,, 3]),
        `CD4 200-349` = rowSums(stan_data$CD4[,, 4]),
        `CD4<200` = rowSums(stan_data$CD4[,, 5]),
        year = seq(0, plot_years, by = 0.25),
        group = "Simulated values"
      )
    } else {
      tibble::tibble(
        `Recently acquired` = stan_data$CD4[, 1],
        `CD4>=500` = stan_data$CD4[, 2],
        `CD4 350-499` = stan_data$CD4[, 3],
        `CD4 200-349` = stan_data$CD4[, 4],
        `CD4<200` = stan_data$CD4[, 5],
        year = seq(0, plot_years, by = 0.25),
        group = "Simulated values"
      )
    }
  } else if (stan_data$config$migration) {
    tibble::tibble(
      `CD4>=500` = rowSums(stan_data$CD4[,, 1]),
      `CD4 350-499` = rowSums(stan_data$CD4[,, 2]),
      `CD4 200-349` = rowSums(stan_data$CD4[,, 3]),
      `CD4<200` = rowSums(stan_data$CD4[,, 4]),
      year = seq(0, plot_years, by = 0.25),
      group = "Simulated values"
    )
  } else {
    tibble::tibble(
      `CD4>=500` = stan_data$CD4[, 1],
      `CD4 350-499` = stan_data$CD4[, 2],
      `CD4 200-349` = stan_data$CD4[, 3],
      `CD4<200` = stan_data$CD4[, 4],
      year = seq(0, plot_years, by = 0.25),
      group = "Simulated values"
    )
  }

  cd4_strata <- cd4_strata |>
    tidyr::pivot_longer(
      cols = tidyselect::starts_with(c("CD4", "Recently")),
      names_to = "strata",
      values_to = "value"
    )

  detect_prob <- tibble::tibble(
    value = sim_data$detect_prob,
    year = seq(0, plot_years, by = 0.25),
    group = "Simulated values"
  )

  return(
    list(
      "stan_data" = stan_data,
      "sim_data" = sim_data,
      "undiag_prev" = undiag_prev,
      "undiag_prev_age" = undiag_prev_age,
      "incidence" = incidence,
      "incidence_age" = incidence_age,
      "migration" = migration,
      "diag_prob" = diag_prob,
      "diag_prob_age" = diag_prob_age,
      "diagnoses" = diagnoses,
      "diagnoses_age" = diagnoses_age,
      "migration_prob" = migration_prob,
      "ratio_abroad_uk" = ratio_abroad_uk,
      "transfers" = transfers,
      "cd4_strata" = cd4_strata,
      "detect_prob" = detect_prob
    )
  )
}

#' Get bias between HIV model estimates and simulation results.
#'
#' @param hiv_list A list of HIV model results.
#' @param model_estimates The model estimates for comparison.
#' @param type The variable type to calculate bias for.
#' @return A data frame of bias values.
get_bias <- function(
  hiv_list,
  model_estimates,
  type
) {
  bias <- dplyr::left_join(
    get_estimates(hiv_list, type),
    model_estimates[[type]],
    by = c("year")
  ) |>
    dplyr::filter(.data$year > 0) |>
    dplyr::mutate(
      bias = .data$`50%` - .data$value
    ) |>
    dplyr::rename(group = "group.x")
}

#' Get bias between HIV model estimates and simulation results for strata.
#'
#' @param hiv_list A list of HIV model results.
#' @param model_estimates The model estimates for comparison.
#' @param type The variable type to calculate bias for.
#' @return A data frame of bias values for strata.
get_bias_strata <- function(
  hiv_list,
  model_estimates,
  type
) {
  bias <- dplyr::left_join(
    get_estimates(hiv_list, type),
    model_estimates[[type]],
    by = c("year", "strata")
  ) |>
    dplyr::filter(.data$year > 0) |>
    dplyr::mutate(
      bias = .data$`50%` - .data$value
    ) |>
    dplyr::rename(group = "group.x")
}

#' Get mean squared error (MSE) list for HIV model estimates.
#'
#' @param bias_df A data frame of bias values.
#' @return A data frame of MSE values.
get_mse <- function(bias_df) {
  bias_df |>
    dplyr::group_by(.data$group) |>
    dplyr::summarise(mse = mean(.data$bias^2), .groups = "drop")
}

#' Make a list of mean squared error (MSE) values for HIV model estimates.
#'
#' @param hiv_list A list of HIV model results.
#' @return A list of bias and MSE values for incidence, migration, proportion of UK-acquired infections that are migrants, and migration probabilities.
#' @export
make_mse_list <- function(
  hiv_list
) {
  # make hiv_list a list of lists if not already
  if (!("incidence" %in% names(hiv_list[[1]]))) {
    hiv_list <- list(hiv_list)
  }

  simulated_quantities <- extract_simulated_quantities(hiv_list[[1]])

  out <- list()

  out$incidence_bias <- get_bias(
    hiv_list,
    simulated_quantities,
    "incidence"
  )
  out$incidence_mse <- get_mse(out$incidence_bias)

  if (hiv_list[[1]]$stan_data$config$migration) {
    out$migration_bias <- get_bias(
      hiv_list,
      simulated_quantities,
      "migration"
    )
    out$ratio_abroad_uk_bias <- get_bias(
      hiv_list,
      simulated_quantities,
      "ratio_abroad_uk"
    )
    out$mp_bias <- get_bias_strata(
      hiv_list,
      simulated_quantities,
      "migration_prob"
    )

    out$migration_mse <- get_mse(out$migration_bias)
    out$ratio_abroad_uk_mse <- get_mse(out$ratio_abroad_uk_bias)
    out$mp_mse <- get_mse(out$mp_bias)
  }
  return(out)
}

#' Box plots of mean squared error (MSE).
#'
#' @param hiv_list A list returned from the `postproc` function containing model estimates.
#' @param quantity The quantity to plot MSE for (default is "incidence", other options are "migration", "ratio_abroad_uk", and "mp").
#' @param offset The offset for text labels.
#' @param lab A logical indicating whether to add labels.
#' @return A ggplot object of the MSE box plot.
#' @export
mse_plot <- function(
  hiv_list,
  quantity = "incidence",
  offset = 0.1,
  lab = TRUE
) {
  data <- make_mse_list(
    hiv_list
  )[[paste0(quantity, "_mse")]]

  # create box plot
  p1 <- data |>
    ggplot2::ggplot() +
    ggplot2::aes(x = .data$group, y = .data$mse) +
    ggplot2::geom_boxplot(
      coef = 1.5,
      alpha = 0.2,
      width = 0.4,
      outliers = FALSE
    ) +
    ggplot2::labs(
      x = "",
      y = "Squared difference"
    ) +
    ggplot2::scale_y_continuous(
      breaks = scales::pretty_breaks(),
      labels = scales::label_comma()
    ) +
    ggplot2::theme_minimal(15) +
    ggplot2::theme(
      legend.position = "none",
      panel.grid.minor = ggplot2::element_blank()
    )

  if (lab) {
    p1 <- p1 +
      ggplot2::geom_text(
        data = data |>
          dplyr::group_by(.data$group) |>
          dplyr::summarise(
            bias = mean(.data$mse),
            pos = stats::quantile(.data$mse, 0.75) + offset,
            # Format bias based on its magnitude
            formatted_bias = dplyr::case_when(
              abs(.data$bias) >= 1000 ~ format(
                round(.data$bias, 0),
                big.mark = ","
              ),
              abs(.data$bias) >= 1 ~ format(round(.data$bias, 1), nsmall = 1),
              TRUE ~ format(
                round(.data$bias, 5),
                nsmall = 5,
                scientific = FALSE
              )
            )
          ),
        ggplot2::aes(
          y = .data$pos,
          label = paste0("PMSE: ", .data$formatted_bias),
        ),
        nudge_x = 0.3,
        size = 4
      )
  }
  return(p1)
}

#' Box plots of bias.
#'
#' @param hiv_list A list returned from the `postproc` function containing model estimates.
#' @param quantity The quantity to plot bias for (default is "incidence", other options are "migration", "ratio_abroad_uk", and "mp").
#' @param offset The offset for text labels.
#' @param lab A logical indicating whether to add labels.
#' @return A ggplot object of the bias box plot.
#' @export
bias_plot <- function(
  hiv_list,
  quantity = "incidence",
  offset = 0.1,
  lab = TRUE
) {
  data <- make_mse_list(
    hiv_list
  )[[paste0(quantity, "_bias")]]

  p1 <- data |>
    ggplot2::ggplot() +
    ggplot2::aes(x = .data$group, y = .data$bias) +
    ggplot2::geom_hline(yintercept = 0, colour = "grey") +
    ggplot2::geom_boxplot(coef = Inf, alpha = 0.2, width = 0.4) +
    ggplot2::labs(
      x = "",
      y = "Difference"
    ) +
    ggplot2::scale_y_continuous(
      breaks = scales::pretty_breaks(),
      labels = scales::label_comma()
    ) +
    ggplot2::theme_minimal(15) +
    ggplot2::theme(
      legend.position = "none",
      panel.grid.minor = ggplot2::element_blank()
    )

  if (lab) {
    p1 <- p1 +
      ggplot2::geom_text(
        data = data |>
          dplyr::group_by(.data$group) |>
          dplyr::summarise(
            meanbias = mean(.data$bias),
            pos = stats::quantile(.data$bias, 0.75) + offset,
            # Format bias based on its magnitude
            formatted_bias = dplyr::case_when(
              abs(.data$meanbias) >= 1000 ~ format(
                round(.data$meanbias, 0),
                big.mark = ","
              ),
              abs(.data$meanbias) >= 1 ~ format(
                round(.data$meanbias, 1),
                nsmall = 1
              ),
              TRUE ~ format(
                round(.data$meanbias, 5),
                nsmall = 5,
                scientific = FALSE
              )
            )
          ),
        ggplot2::aes(
          y = .data$pos,
          label = paste0("Bias: ", .data$formatted_bias),
        ),
        nudge_x = 0.3,
        size = 4
      )
  }
  return(p1)
}
