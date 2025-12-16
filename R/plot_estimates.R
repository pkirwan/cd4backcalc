#' Plot estimated model quantities
#'
#' Generate a ggplot2 object plotting estimated quantities over time, optionally compared to the simulated values.
#'
#' @param hiv_list A list returned from the `postproc` function containing model estimates.
#' @param quantity A string specifying the quantity to plot. Options are "incidence", "undiag_prev", "migration", "diag_prob", "diagnoses", "migration_prob", "transfers", "prop_uk", and "detect_prob".
#' @param chain_results A logical indicating whether to plot chain-by-chain estimates (default is FALSE).
#' @return A ggplot2 object.
#' @export
plot_estimates <- function(
  hiv_list,
  quantity,
  chain_results = FALSE
) {
  plot_title <- get_plot_title(quantity)

  # chain-by-chain plotting if requested
  if (chain_results) {
    hiv_list <- dplyr::bind_rows(
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
    ggplot2::aes(.data$year, .data$value, group = .data$group) +
    ggplot2::scale_linetype_manual(values = c(2, 1, 2)) +
    ggplot2::geom_line(ggplot2::aes(y = .data$`50%`), color = "purple") +
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
    ggplot2::labs(
      title = plot_title,
      x = "Year",
      y = plot_title
    ) +
    ggplot2::theme_bw()

  # add simulated values if data were simulated
  if ("sim_data" %in% names(hiv_list)) {
    sim_estimates <- extract_simulated_quantities(hiv_list)
    p1 <- p1 +
      ggplot2::geom_point(size = 0.4, data = sim_estimates[[quantity]])
  }

  # custom y-axis for probabilities
  if (
    quantity %in% c("diag_prob", "prop_uk", "detect_prob", "migration_prob")
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

  return(p1)
}

#' Plot diagnostics for fitted model
#'
#' Generate diagnostic plots for a fitted back-calculation model.
#'
#' @param fit A fitted model object returned from the `fit_stan_model` or `run_backcalc` function.
#' @param ntraces Number of parameter traces to plot for each diagnostic.
#' @return A ggplot2 object containing diagnostic plots.
#' @export
plot_diagnostics <- function(fit, ntraces = 6) {

  # check if hiv_list has been provided
  if ("fit" %in% names(fit)) {
    fit <- fit$fit
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
  } else if (quantity == "prop_uk") {
    "Estimated born-abroad:native-born ratio"
  } else if (quantity == "detect_prob") {
    "Estimated detection probability for transfers of care"
  } else {
    stop("Invalid quantity specified.")
  }
  return(plot_title)
}
