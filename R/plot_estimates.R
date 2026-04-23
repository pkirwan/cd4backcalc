#' Plot simulated quantities
#'
#' Generate a ggplot2 object plotting simulated quantities over time.
#'
#' @param sim_diags A list containing the simulated diagnosis data objects sim_data and stan_data.
#' @param quantity A string specifying the quantity to plot. Options are "incidence", "undiag_prev", "undiag_migration", "all_migration", "diag_prob", "diagnoses", "migration_prob", "transfers", "ratio_abroad_uk", and "detect_prob".
#' @return A ggplot2 object.
#' @export
plot_simulations <- function(
  sim_diags,
  quantity
) {
  plot_title <- get_plot_title(quantity)
  # replace "Estimated" with "Simulated"
  plot_title <- paste0("Simulated ", plot_title)
  sim_estimates <- extract_simulated_quantity(sim_diags, quantity)
  sim_estimates_det <- extract_simulated_quantity(
    sim_diags,
    quantity,
    deterministic = TRUE
  )
  sim_estimates <- sim_estimates |>
    dplyr::mutate(group = "Sampled values") |>
    dplyr::bind_rows(
      sim_estimates_det |>
        dplyr::mutate(group = "Generative values")
    )

  # exclude "all" age_band from age-facetted plots
  if ("age_band" %in% names(sim_estimates) && grepl("_age$", quantity)) {
    sim_estimates <- sim_estimates |>
      dplyr::filter(.data$age_band != "all")
  }

  model_colors <- stats::setNames(
    c("mediumturquoise", "mediumorchid"),
    c("Sampled values", "Generative values")
  )

  # plot simulated values as points and expected values as lines
  p1 <- sim_estimates |>
    ggplot2::ggplot() +
    ggplot2::aes(
      .data$year,
      .data$value,
      color = .data$group,
      group = .data$group
    ) +
    ggplot2::geom_point(
      data = sim_estimates |>
        dplyr::filter(.data$group == "Sampled values"),
      size = 1.1,
      shape = 3
    ) +
    ggplot2::geom_line(
      data = sim_estimates |>
        dplyr::filter(.data$group == "Generative values")
    ) +
    ggplot2::labs(
      title = plot_title,
      x = "Year",
      # remove "by CD4 strata" and "by age" from y-axis label for better readability
      y = sub(
        "^(.)",
        "\\U\\1",
        gsub(
          ", by CD4 strata| by age|, by CD4 strata and age|Simulated | for transfers of care",
          "",
          get_plot_title(quantity)
        ),
        perl = TRUE
      ),
      color = ""
    ) +
    ggplot2::scale_color_manual(
      values = model_colors,
      name = ""
    ) +
    ggplot2::theme_classic() +
    ggplot2::scale_x_continuous(
      expand = ggplot2::expansion(mult = c(0.01, 0.01))
    ) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_line(color = "grey80"),
      legend.position = "bottom"
    )

  p1 <- do_facetting(p1, quantity)

  return(p1)
}

#' Plot estimated model quantities
#'
#' Generate a ggplot2 object plotting estimated quantities over time, optionally compared to the simulated values.
#'
#' @param hiv_list A list returned from the `postproc` function containing model estimates.
#' @param quantity A string specifying the quantity to plot. Options are "incidence", "undiag_prev", "undiag_migration", "all_migration", "diag_prob", "diagnoses", "migration_prob", "transfers", "ratio_abroad_uk", and "detect_prob".
#' @param start_year An integer specifying the starting year for the x-axis (default is 0).
#' @param chain_results A logical indicating whether to plot chain-by-chain estimates (default is FALSE).
#' @param hiv_list_2 An optional second list of model estimates for comparison (default is NULL).
#' @param model_1_name A string specifying the name of the first model for the legend (default is "Model 1").
#' @param model_2_name A string specifying the name of the second model for the legend (default is "Model 2").
#' @param credible_intervals A logical indicating whether to plot 95% credible intervals (default is TRUE).
#' @return A ggplot2 object.
#' @export
plot_estimates <- function(
  hiv_list,
  quantity,
  start_year = 0,
  chain_results = FALSE,
  hiv_list_2 = NULL,
  model_1_name = "Model 1",
  model_2_name = "Model 2",
  credible_intervals = TRUE
) {
  hiv_list <- purrr::compact(hiv_list)
  if (length(hiv_list) == 0) {
    stop("All entries in hiv_list are NULL.")
  }
  # make hiv_list a list of lists if not already
  if (
    !("incidence" %in% names(hiv_list[[1]])) &&
      !("incidence_1" %in% names(hiv_list[[1]]))
  ) {
    hiv_list <- list(hiv_list)
  }

  if (!is.null(hiv_list_2)) {
    hiv_list_2 <- purrr::compact(hiv_list_2)
    if (length(hiv_list_2) == 0) {
      hiv_list_2 <- NULL
    }
  }
  if (!is.null(hiv_list_2)) {
    if (
      !("incidence" %in% names(hiv_list_2[[1]])) &&
        !("incidence_1" %in% names(hiv_list_2[[1]]))
    ) {
      hiv_list_2 <- list(hiv_list_2)
    }
    # add group variable to hiv_list_2 estimates
    hiv_list_2 <- lapply(hiv_list_2, function(x) {
      x$group <- model_2_name
      return(x)
    })
    # add group variable to hiv_list estimates
    hiv_list <- lapply(hiv_list, function(x) {
      x$group <- model_1_name
      return(x)
    })
    # combine hiv_list and hiv_list_2 into one list
    hiv_list <- c(hiv_list, hiv_list_2)
  }

  plot_title <- paste0("Estimated ", get_plot_title(quantity))

  # chain-by-chain plotting if requested
  if (chain_results) {
    # derive the number of chains from hiv_list
    n_chains <- sum(grepl("incidence_[0-9]", names(hiv_list[[1]])))
    quantiles <- tibble::tibble()
    for (i in 1:n_chains) {
      quantiles <- quantiles |>
        dplyr::bind_rows(
          get_estimates(hiv_list, paste0(quantity, "_", i)) |>
            dplyr::mutate(chain = paste0("Chain ", i))
        )
    }
    quantiles <- quantiles |>
      dplyr::mutate(group = paste0(.data$group, ", ", .data$chain))
  } else {
    quantiles <- get_estimates(hiv_list, quantity)
  }

  # filter out zero time point for incidence and diagnosis plots
  # by construction these are always zero at time zero
  if (
    quantity %in% c("incidence", "incidence_age", "diagnoses", "diagnoses_age")
  ) {
    quantiles <- quantiles |>
      dplyr::filter(.data$year > 0)
  }

  if (
    quantity %in%
      c("diag_prob", "diag_prob_age", "diag_prob_mig", "diag_prob_mig_age")
  ) {
    quantiles <- quantiles |>
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
  } else if (
    quantity %in%
      c("migration_prob", "migration_prob_age")
  ) {
    quantiles <- quantiles |>
      dplyr::mutate(
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
  }

  # exclude "all" age_band from age-facetted plots
  if ("age_band" %in% names(quantiles) && grepl("_age$", quantity)) {
    quantiles <- quantiles |>
      dplyr::filter(.data$age_band != "all")
  }

  p1 <- quantiles |>
    ggplot2::ggplot() +
    ggplot2::aes(
      x = .data$year + start_year,
      y = .data$`50%`,
      color = .data$group,
      fill = .data$group,
      group = interaction(.data$group, .data$iteration)
    )

  if (credible_intervals) {
    p1 <- p1 +
      ggplot2::geom_line(
        ggplot2::aes(y = .data$`2.5%`),
        linewidth = 0.5,
        alpha = 0.5
      ) +
      ggplot2::geom_line(
        ggplot2::aes(y = .data$`97.5%`),
        linewidth = 0.5,
        alpha = 0.5
      ) +
      ggplot2::geom_ribbon(
        ggplot2::aes(ymin = .data$`2.5%`, ymax = .data$`97.5%`),
        alpha = 0.2
      )
  }

  p1 <- p1 +
    ggplot2::geom_line(linewidth = 1, alpha = 0.7) +
    ggplot2::labs(
      title = plot_title,
      x = "Year",
      # remove "by CD4 strata" and "by age" from y-axis label for better readability
      y = sub(
        "^(.)",
        "\\U\\1",
        gsub(
          ", by CD4 strata| by age|, by CD4 strata and age|Estimated | for transfers of care",
          "",
          get_plot_title(quantity)
        ),
        perl = TRUE
      ),
      color = "",
      fill = ""
    ) +
    ggplot2::theme_classic() +
    ggplot2::scale_x_continuous(
      expand = ggplot2::expansion(mult = c(0.01, 0.01))
    ) +
    ggplot2::theme(
      legend.position = "bottom",
      panel.grid.major = ggplot2::element_line(color = "grey80")
    )

  # define model colors
  if (!is.null(hiv_list_2)) {
    model_colors <- stats::setNames(
      c("mediumorchid", "coral", "mediumturquoise"),
      c(model_1_name, model_2_name, "Generative values")
    )
  } else {
    model_colors <- stats::setNames(
      c(rep("mediumorchid", length(unique(quantiles$group))), "mediumturquoise"),
      c(unique(quantiles$group), "Generative values")
    )
  }

  p1 <- p1 +
    ggplot2::scale_color_manual(
      values = model_colors,
      name = ""
    ) +
    ggplot2::scale_fill_manual(
      values = model_colors,
      name = ""
    )

  # add simulated values if data were simulated
  # use deterministic (expected) values so comparison is against lambda, not rpois(lambda)
  if ("sim_data" %in% names(hiv_list[[1]])) {
    sim_estimates <- extract_simulated_quantity(
      hiv_list[[1]],
      quantity,
      deterministic = TRUE
    )
    if ("age_band" %in% names(sim_estimates) && grepl("_age$", quantity)) {
      sim_estimates <- sim_estimates |>
        dplyr::filter(.data$age_band != "all")
    }
    p1 <- p1 +
      ggplot2::geom_point(
        data = sim_estimates,
        ggplot2::aes(y = .data$value, group = .data$group),
        size = 1.1,
        shape = 3
      )
  } else if (quantity %in% c("diagnoses", "diagnoses_age", "transfers")) {
    # add observed diagnoses if real data
    observed_data <- extract_real_data(hiv_list[[1]], quantity)
    if ("age_band" %in% names(observed_data) && grepl("_age$", quantity)) {
      observed_data <- observed_data |>
        dplyr::filter(.data$age_band != "all")
    }
    p1 <- p1 +
      ggplot2::geom_point(
        data = observed_data,
        ggplot2::aes(y = .data$value, group = .data$group),
        size = 1.1,
        shape = 3
      )
  }

  p1 <- do_facetting(p1, quantity)

  return(p1)
}

#' Plot diagnostics for fitted model
#'
#' Generate diagnostic plots for a fitted back-calculation model.
#'
#' @param hiv_list A list returned from the `run_backcalc` function containing model estimates and fit objects.
#' @param ntraces Number of parameter traces to plot for each diagnostic.
#' @return A ggplot2 object containing diagnostic plots.
#' @export
plot_diagnostics <- function(hiv_list, ntraces = 4) {
  hiv_list <- purrr::compact(hiv_list)
  if (length(hiv_list) == 0) {
    stop("All entries in hiv_list are NULL.")
  }
  # check if hiv_list is a list of results or a single result
  if ("fit" %in% names(hiv_list[[1]])) {
    fit <- hiv_list[[1]]$fit
    message(
      "Traceplots for the first iteration in the list, for other iterations please provide the appropriate fit object."
    )
  } else if ("fit" %in% names(hiv_list)) {
    fit <- hiv_list$fit
  } else {
    stop("No fitted model found in hiv_list.")
  }

  draws <- fit$draws(format = "matrix")

  # log incidence
  p1 <- bayesplot::mcmc_trace(
    draws,
    pars = colnames(draws)[sort(sample(
      grep("^lambda\\[", colnames(draws)),
      ntraces
    ))],
    facet_args = list(nrow = 1)
  ) +
    # fewer x-axis ticks for better readability
    ggplot2::scale_x_continuous(
      breaks = scales::pretty_breaks(n = 4)
    )
  # log posterior density
  p2 <- bayesplot::mcmc_trace(draws, pars = "lp__")

  if (any(grepl("gamma", colnames(draws)))) {
    # detection probabilities
    p3 <- bayesplot::mcmc_trace(
      draws,
      pars = colnames(draws)[sort(sample(
        grep("^detect_prob\\[", colnames(draws)),
        ntraces
      ))],
      facet_args = list(nrow = 1)
    )
    # proportion of UK-acquired who are migrants
    p4 <- bayesplot::mcmc_trace(
      draws,
      pars = colnames(draws)[sort(sample(
        grep("^p\\[", colnames(draws)),
        ntraces
      ))],
      facet_args = list(nrow = 1)
    )
    # log migration incidence
    p5 <- bayesplot::mcmc_trace(
      draws,
      pars = colnames(draws)[sort(sample(
        grep("^gamma\\[", colnames(draws)),
        ntraces
      ))],
      facet_args = list(nrow = 1)
    )
    # migration probabilities
    p6 <- bayesplot::mcmc_trace(
      draws,
      pars = colnames(draws)[sort(sample(
        grep("^m\\[", colnames(draws)),
        ntraces
      ))],
      facet_args = list(nrow = 1)
    )

    p0 <- patchwork::wrap_plots(
      p1,
      p2,
      p3,
      p4,
      p5,
      p6,
      ncol = 1,
      heights = c(2, 2, 2, 2, 2, 1)
    ) &
      ggplot2::theme(text = ggplot2::element_text(family = "sans"))
  } else {
    p0 <- patchwork::wrap_plots(
      p1,
      p2,
      ncol = 1,
      heights = c(2, 1)
    ) &
      ggplot2::theme(text = ggplot2::element_text(family = "sans"))
  }

  return(p0)
}

#' Helper function to make plot titles
#'
#' Generate appropriate plot titles based on quantity.
#'
#' @param quantity A string specifying the quantity.
#' @return A string containing the plot title.
get_plot_title <- function(quantity) {
  plot_title <- if (quantity == "incidence") {
    "HIV incidence"
  } else if (quantity == "undiag_prev") {
    "undiagnosed HIV prevalence"
  } else if (quantity == "undiag_migration") {
    "migration with undiagnosed HIV"
  } else if (quantity == "all_migration") {
    "overall migration with HIV (including diagnosed and undiagnosed)"
  } else if (quantity == "diag_prob") {
    "HIV diagnosis probabilities, by CD4 strata"
  } else if (quantity == "diagnoses") {
    "HIV diagnoses"
  } else if (quantity == "migration_prob") {
    "migration probabilities, by CD4 strata"
  } else if (quantity == "transfers") {
    "transfers of care"
  } else if (quantity == "ratio_abroad_uk") {
    "born-abroad:native-born ratio"
  } else if (quantity == "detect_prob") {
    "detection probability for transfers of care"
  } else if (quantity == "incidence_age") {
    "HIV incidence by age"
  } else if (quantity == "undiag_prev_age") {
    "undiagnosed HIV prevalence by age"
  } else if (quantity == "undiag_migration_age") {
    "migration with undiagnosed HIV by age"
  } else if (quantity == "diagnoses_age") {
    "HIV diagnoses by age"
  } else if (quantity == "diag_prob_age") {
    "HIV diagnosis probabilities, by CD4 strata and age"
  } else if (quantity == "diag_prob_mig") {
    "HIV diagnosis probabilities (migrants), by CD4 strata"
  } else if (quantity == "diag_prob_mig_age") {
    "HIV diagnosis probabilities (migrants), by CD4 strata and age"
  } else if (quantity == "incidence_trend") {
    "HIV incidence trend by age"
  } else if (quantity == "peak_incidence") {
    "year of peak HIV incidence by age"
  } else {
    stop("Invalid quantity specified.")
  }
  return(plot_title)
}

#' Helper function to do facetting for plots
#'
#' @param p1 A ggplot2 object.
#' @param quantity A string specifying the quantity.
#' @return A ggplot2 object.
do_facetting <- function(p1, quantity) {
  # custom y-axis for probabilities
  if (
    quantity %in%
      c(
        "diag_prob", "diag_prob_mig", "diag_prob_age", "diag_prob_mig_age",
        "ratio_abroad_uk", "detect_prob", "migration_prob", "migration_prob_age"
      )
  ) {
    p1 <- p1 +
      ggplot2::scale_y_continuous(
        labels = scales::percent,
        limits = c(0, NA),
        expand = ggplot2::expansion(mult = c(0, 0.05))
      )
  } else {
    p1 <- p1 +
      ggplot2::scale_y_continuous(
        limits = c(0, NA),
        expand = ggplot2::expansion(mult = c(0, 0.05)),
        labels = scales::comma
      )
  }

  # facet by strata for diagnosis and migration probabilities
  if (quantity %in% c("diag_prob", "diag_prob_mig", "migration_prob")) {
    p1 <- p1 +
      ggplot2::facet_grid(cols = ggplot2::vars(.data$strata), scales = "free_y")
  }

  # facet by age if age-dependent model
  if (
    quantity %in%
      c(
        "incidence_age",
        "undiag_prev_age",
        "undiag_migration_age",
        "diagnoses_age"
      )
  ) {
    p1 <- p1 +
      ggplot2::facet_grid(
        cols = ggplot2::vars(.data$age_band),
        scales = "free_y"
      )
  }

  # facet by strata and age for age-dependent diagnosis and migration probabilities
  if (quantity %in% c("diag_prob_age", "diag_prob_mig_age", "migration_prob_age")) {
    p1 <- p1 +
      ggplot2::facet_grid(
        rows = ggplot2::vars(.data$strata),
        cols = ggplot2::vars(.data$age_band),
        scales = "free_y"
      )
  }

  return(p1)
}
