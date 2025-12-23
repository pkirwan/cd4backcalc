#' Plot simulated quantities
#'
#' Generate a ggplot2 object plotting simulated quantities over time.
#'
#' @param sim_diags A list containing the simulated diagnosis data objects sim_data and stan_data.
#' @param quantity A string specifying the quantity to plot. Options are "incidence", "undiag_prev", "migration", "diag_prob", "diagnoses", "migration_prob", "transfers", "ratio_abroad_uk", and "detect_prob".
#' @return A ggplot2 object.
#' @export
plot_simulations <- function(
  sim_diags,
  quantity
) {
  plot_title <- get_plot_title(quantity)
  # replace "Estimated" with "Simulated"
  plot_title <- sub("Estimated", "Simulated", plot_title)
  sim_estimates <- extract_simulated_quantity(sim_diags, quantity)

  p1 <- sim_estimates |>
    ggplot2::ggplot() +
    ggplot2::aes(.data$year, .data$value, group = .data$group) +
    ggplot2::geom_point(color = "purple") +
    ggplot2::labs(
      title = plot_title,
      x = "Year",
      y = plot_title
    ) +
    ggplot2::theme_bw()

  p1 <- do_facetting(p1, quantity)

  return(p1)
}

#' Plot estimated model quantities
#'
#' Generate a ggplot2 object plotting estimated quantities over time, optionally compared to the simulated values.
#'
#' @param hiv_list A list returned from the `postproc` function containing model estimates.
#' @param quantity A string specifying the quantity to plot. Options are "incidence", "undiag_prev", "migration", "diag_prob", "diagnoses", "migration_prob", "transfers", "ratio_abroad_uk", and "detect_prob".
#' @param chain_results A logical indicating whether to plot chain-by-chain estimates (default is FALSE).
#' @return A ggplot2 object.
#' @export
plot_estimates <- function(
  hiv_list,
  quantity,
  chain_results = FALSE
) {
  # make hiv_list a list of lists if not already
  if (!("incidence" %in% names(hiv_list[[1]]))) {
    hiv_list <- list(hiv_list)
  }

  plot_title <- get_plot_title(quantity)

  # chain-by-chain plotting if requested
  if (chain_results) {
    quantiles <- dplyr::bind_rows(
      get_estimates(hiv_list, paste0(quantity, "_1")) |>
        dplyr::mutate(chain = "Chain 1"),
      get_estimates(hiv_list, paste0(quantity, "_2")) |>
        dplyr::mutate(chain = "Chain 2"),
      get_estimates(hiv_list, paste0(quantity, "_3")) |>
        dplyr::mutate(chain = "Chain 3"),
      get_estimates(hiv_list, paste0(quantity, "_4")) |>
        dplyr::mutate(chain = "Chain 4")
    ) |>
      dplyr::mutate(group = paste0(.data$group, ", ", .data$chain))
  } else {
    quantiles <- get_estimates(hiv_list, quantity)
  }

  p1 <- quantiles |>
    ggplot2::ggplot() +
    ggplot2::aes(.data$year, .data$`50%`, group = .data$group) +
    ggplot2::scale_linetype_manual(values = c(2, 1, 2)) +
    ggplot2::geom_line(
      ggplot2::aes(y = .data$`2.5%`),
      linetype = "dashed",
      color = "gray"
    ) +
    ggplot2::geom_line(
      ggplot2::aes(y = .data$`97.5%`),
      linetype = "dashed",
      color = "gray"
    ) +
    ggplot2::geom_line(color = "purple") +
    ggplot2::labs(
      title = plot_title,
      x = "Year",
      y = plot_title
    ) +
    ggplot2::theme_bw()

  # add simulated values if data were simulated
  # simulated values are the same for all iterations, so just use first
  if ("sim_data" %in% names(hiv_list[[1]])) {
    sim_estimates <- extract_simulated_quantity(hiv_list[[1]], quantity)
    p1 <- p1 +
      ggplot2::geom_point(
        data = sim_estimates,
        ggplot2::aes(y = .data$value),
        size = 0.4
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
plot_diagnostics <- function(hiv_list, ntraces = 6) {
  # check if hiv_list is a list of results or a single result
  if ("fit" %in% names(hiv_list[[1]])) {
    fit <- hiv_list[[1]]$fit
    message(
      "Traceplots will be shown for the first iteration of the list, for other iterations please provide the appropriate fit object."
    )
  } else if ("fit" %in% names(hiv_list)) {
    fit <- hiv_list$fit
  } else {
    stop("No fitted model found in hiv_list.")
  }

  # log incidence
  p1 <- rstan::traceplot(
    fit,
    pars = names(fit)[sort(sample(grep("^lambda\\[", names(fit)), ntraces))],
    nrow = 1
  )
  # log posterior density
  p2 <- rstan::traceplot(fit, pars = "lp__")

  if (any(grepl("gamma", names(fit)))) {
    # detection probabilities
    p3 <- rstan::traceplot(
      fit,
      pars = names(fit)[sort(sample(
        grep("^detect_prob\\[", names(fit)),
        ntraces
      ))],
      nrow = 1
    )
    # proportion of UK-acquired who are migrants
    p4 <- rstan::traceplot(
      fit,
      pars = names(fit)[sort(sample(grep("^p\\[", names(fit)), ntraces))],
      nrow = 1
    )
    # log migration incidence
    p5 <- rstan::traceplot(
      fit,
      pars = names(fit)[sort(sample(grep("^gamma\\[", names(fit)), ntraces))],
      nrow = 1
    )
    # migration probabilities
    p6 <- rstan::traceplot(
      fit,
      pars = names(fit)[sort(sample(grep("^m\\[", names(fit)), ntraces))],
      nrow = 1
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
    )
  } else {
    p0 <- patchwork::wrap_plots(
      p1,
      p2,
      ncol = 1,
      heights = c(2, 1)
    )
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
    "Estimated HIV incidence"
  } else if (quantity == "undiag_prev") {
    "Estimated undiagnosed HIV prevalence"
  } else if (quantity == "migration") {
    "Estimated migration with undiagnosed HIV"
  } else if (quantity == "diag_prob") {
    "Estimated HIV diagnosis probabilities, by CD4 strata"
  } else if (quantity == "diagnoses") {
    "Estimated HIV diagnoses"
  } else if (quantity == "migration_prob") {
    "Estimated migration probabilities, by CD4 strata"
  } else if (quantity == "transfers") {
    "Estimated transfers of care"
  } else if (quantity == "ratio_abroad_uk") {
    "Estimated born-abroad:native-born ratio"
  } else if (quantity == "detect_prob") {
    "Estimated detection probability for transfers of care"
  } else if (quantity == "incidence_age") {
    "Estimated HIV incidence by age"
  } else if (quantity == "undiag_prev_age") {
    "Estimated undiagnosed HIV prevalence by age"
  } else if (quantity == "migration_age") {
    "Estimated migration with undiagnosed HIV by age"
  } else if (quantity == "diagnoses_age") {
    "Estimated HIV diagnoses by age"
  } else if (quantity == "diag_prob_age") {
    "Estimated HIV diagnosis probabilities by CD4 strata and age"
  } else if (quantity == "incidence_trend") {
    "Estimated HIV incidence trend by age"
  } else if (quantity == "peak_incidence") {
    "Estimated year of peak HIV incidence by age"
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
      c("diag_prob", "ratio_abroad_uk", "detect_prob", "migration_prob")
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
        expand = ggplot2::expansion(mult = c(0, 0.05))
      )
  }

  # facet by strata for diagnosis and migration probabilities
  if (quantity %in% c("diag_prob", "migration_prob")) {
    p1 <- p1 +
      ggplot2::facet_wrap(~ .data$strata, scales = "free_y")
  }

  # facet by age if age-dependent model
  if (
    quantity %in%
      c(
        "incidence_age",
        "undiag_prev_age",
        "migration_age",
        "diagnoses_age"
      )
  ) {
    p1 <- p1 +
      ggplot2::facet_wrap(~ .data$age_band, scales = "free_y")
  }

  # facet by strata and age for age-dependent diagnosis and migration probabilities
  if (quantity %in% c("diag_prob_age", "migration_prob_age")) {
    p1 <- p1 +
      ggplot2::facet_wrap(
        ~ .data$strata + .data$age_band,
        scales = "free_y",
        ncol = 4
      )
  }

  return(p1)
}
