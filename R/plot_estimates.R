#' Plot simulated or observed quantities
#'
#' Generate a ggplot2 object plotting simulated quantities over time. When
#' `sim_data` is absent, quantities available directly from real-data Stan
#' inputs are plotted instead.
#'
#' @param sim_diags A list containing simulated diagnosis data objects `sim_data`
#'   and `stan_data`, or a real-data result containing `stan_data` and no
#'   `sim_data`.
#' @param quantity A string specifying the quantity to plot. Options include
#'   "incidence", "undiag_prev", "undiag_migration", "all_migration",
#'   "diag_prob", "diagnoses", "migration_prob", "transfers",
#'   "transfers_detected", "ratio_abroad_uk", "detect_prob",
#'   "migrant_hiv_prevalence", "aids_diagnoses", "aids_diagnoses_age",
#'   "cd4_diagnoses", "cd4_diagnoses_age", "transfers_detected_age",
#'   "gamma_obs", "gamma_obs_group", and "recent_arrival_cd4".
#' @return A ggplot2 object.
#' @examples
#' # See ?run_backcalc for a full end-to-end workflow.
#' invisible(NULL)
#' @export
plot_simulations <- function(
  sim_diags,
  quantity
) {
  is_real_data <- is.null(sim_diags$sim_data)
  plot_title <- paste0(
    if (is_real_data) "Observed " else "Simulated ",
    get_plot_title(quantity)
  )

  if (is_real_data) {
    sim_estimates <- extract_real_data(sim_diags, quantity)
    if (is.null(sim_estimates) || nrow(sim_estimates) == 0) {
      stop(
        sprintf(
          "Quantity '%s' is not available as observed data.",
          quantity
        ),
        call. = FALSE
      )
    }
  } else {
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
  }

  if (
    "age_band" %in% names(sim_estimates) &&
      (grepl("_age$", quantity) || quantity == "gamma_obs_group")
  ) {
    sim_estimates <- sim_estimates |>
      dplyr::filter(.data$age_band != "all")
  }

  model_colors <- if (is_real_data) {
    stats::setNames("mediumturquoise", "Observed values")
  } else {
    stats::setNames(
      c("mediumturquoise", "mediumorchid"),
      c("Sampled values", "Generative values")
    )
  }

  # plot simulated values as points and expected values as lines
  p1 <- sim_estimates |>
    ggplot2::ggplot() +
    ggplot2::aes(
      .data$year,
      .data$value,
      color = .data$group,
      group = .data$group
    ) +
    ggplot2::labs(
      title = plot_title,
      x = "Year",
      # remove "by CD4 strata" and "by age" from y-axis label for better readability
      y = sub(
        "^(.)",
        "\\U\\1",
        gsub(
          ", by CD4 strata| by age group|, by CD4 strata and age group|Simulated |Observed | for transfers of care",
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

  if (is_real_data) {
    p1 <- p1 +
      ggplot2::geom_point(
        size = 1.1,
        shape = 3
      )
  } else {
    p1 <- p1 +
      ggplot2::geom_point(
        data = sim_estimates |>
          dplyr::filter(.data$group == "Sampled values"),
        size = 1.1,
        shape = 3
      ) +
      ggplot2::geom_line(
        data = sim_estimates |>
          dplyr::filter(.data$group == "Generative values")
      )
  }

  p1 <- do_facetting(p1, quantity)

  return(p1)
}

#' Plot estimated model quantities
#'
#' Generate a ggplot2 object plotting estimated quantities over time, optionally compared to the simulated values.
#'
#' @param hiv_list A list returned from the `postproc` function containing model
#'   estimates, or a named list of model results. Names are used as plot labels.
#' @param quantity A string specifying the quantity to plot. Options are "incidence", "undiag_prev", "undiag_migration", "all_migration", "diag_prob", "diagnoses", "migration_prob", "transfers", "transfers_detected", "ratio_abroad_uk", "detect_prob", and "migrant_hiv_prevalence".
#' @param start_year An integer specifying the starting year for the x-axis (default is 0).
#' @param chain_results A logical indicating whether to plot chain-by-chain estimates (default is FALSE).
#' @param hiv_list_2 An optional second list of model estimates for comparison (default is NULL).
#' @param model_1_name A string specifying the name of the first model for the legend (default is "Model 1").
#' @param model_2_name A string specifying the name of the second model for the legend (default is "Model 2").
#' @param credible_intervals A logical indicating whether to plot 95% credible intervals (default is TRUE).
#' @param sim_comparison A string specifying which simulated values to compare against when `hiv_list` contains simulated data. Options are `"sampled"` (the stochastic realisation for the particular seed, the default) or `"generative"` (deterministic expected values).
#' @param colours Optional character vector of model or chain colours. Unnamed
#'   values are matched in legend order; named values must match legend labels.
#' @return A ggplot2 object.
#' @examples
#' # See ?run_backcalc for a full end-to-end workflow.
#' invisible(NULL)
#' @export
plot_estimates <- function(
  hiv_list,
  quantity,
  start_year = 0,
  chain_results = FALSE,
  hiv_list_2 = NULL,
  model_1_name = "Model 1",
  model_2_name = "Model 2",
  credible_intervals = TRUE,
  sim_comparison = c("sampled", "generative"),
  colours = NULL
) {
  model_1_name_supplied <- !missing(model_1_name)
  sim_comparison <- match.arg(sim_comparison)
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
  }
  if (model_1_name_supplied || !is.null(hiv_list_2)) {
    # add group variable to hiv_list estimates
    hiv_list <- lapply(hiv_list, function(x) {
      x$group <- model_1_name
      return(x)
    })
  }
  if (!is.null(hiv_list_2)) {
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

  model_groups <- unique(quantiles$group)
  model_linetypes <- make_model_linetypes(model_groups, chain_results)
  model_colors <- make_model_colors(
    model_groups,
    chain_results = chain_results,
    has_second_model = !is.null(hiv_list_2),
    model_1_name = model_1_name,
    model_2_name = model_2_name,
    colours = colours
  )
  legend_breaks <- model_groups

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
        ggplot2::aes(y = .data$`2.5%`, linetype = .data$group),
        linewidth = 0.5,
        alpha = 0.5
      ) +
      ggplot2::geom_line(
        ggplot2::aes(y = .data$`97.5%`, linetype = .data$group),
        linewidth = 0.5,
        alpha = 0.5
      ) +
      ggplot2::geom_ribbon(
        ggplot2::aes(ymin = .data$`2.5%`, ymax = .data$`97.5%`),
        alpha = 0.2
      )
  }

  p1 <- p1 +
    ggplot2::geom_line(
      ggplot2::aes(linetype = .data$group),
      linewidth = 1,
      alpha = 0.7
    ) +
    ggplot2::labs(
      title = plot_title,
      x = "Year",
      # remove "by CD4 strata" and "by age" from y-axis label for better readability
      y = sub(
        "^(.)",
        "\\U\\1",
        gsub(
          ", by CD4 strata| by age group|, by CD4 strata and age group|Estimated | for transfers of care",
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

  sim_ref_label <- if (sim_comparison == "sampled") "Sampled values" else "Generative values"

  # add simulated values if data were simulated
  if ("sim_data" %in% names(hiv_list[[1]])) {
    use_deterministic <- sim_comparison == "generative"
    model_colors <- c(model_colors, stats::setNames("mediumturquoise", sim_ref_label))
    legend_breaks <- c(legend_breaks, sim_ref_label)
    sim_estimates <- extract_simulated_quantity(
      hiv_list[[1]],
      quantity,
      deterministic = use_deterministic
    )
    if (!use_deterministic) {
      sim_estimates <- sim_estimates |>
        dplyr::mutate(group = "Sampled values")
    }
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
  } else if (quantity %in% c(
    "diagnoses",
    "diagnoses_age",
    "transfers_detected"
  )) {
    # add observed diagnoses if real data
    observed_data <- extract_real_data(hiv_list[[1]], quantity)
    if ("age_band" %in% names(observed_data) && grepl("_age$", quantity)) {
      observed_data <- observed_data |>
        dplyr::filter(.data$age_band != "all")
    }

    model_colors <- c(model_colors, "Observed values" = "mediumturquoise")
    legend_breaks <- c(legend_breaks, "Observed values")

    p1 <- p1 +
      ggplot2::geom_point(
        data = observed_data,
        ggplot2::aes(y = .data$value, group = .data$group),
        size = 1.1,
        shape = 3
      )
  }

  legend_linetypes <- stats::setNames(rep("blank", length(legend_breaks)), legend_breaks)
  legend_linetypes[model_groups] <- unname(model_linetypes[model_groups])
  legend_shapes <- stats::setNames(rep(NA, length(legend_breaks)), legend_breaks)
  legend_shapes[setdiff(legend_breaks, model_groups)] <- 3

  p1 <- p1 +
    ggplot2::scale_color_manual(
      values = model_colors,
      breaks = legend_breaks,
      name = ""
    ) +
    ggplot2::scale_fill_manual(
      values = model_colors,
      breaks = legend_breaks,
      name = ""
    ) +
    ggplot2::scale_linetype_manual(
      values = model_linetypes,
      breaks = model_groups,
      name = ""
    ) +
    ggplot2::guides(
      fill = "none",
      linetype = "none",
      color = ggplot2::guide_legend(
        override.aes = list(
          linetype = unname(legend_linetypes[legend_breaks]),
          shape = unname(legend_shapes[legend_breaks]),
          fill = NA
        )
      )
    )

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
#' @examples
#' # See ?run_backcalc for a full end-to-end workflow.
#' invisible(NULL)
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

  if (any(grepl("detect_prob", colnames(draws)))) {
    # detection probabilities
    p3 <- bayesplot::mcmc_trace(
      draws,
      pars = colnames(draws)[sort(sample(
        grep("^detect_prob\\[", colnames(draws)),
        ntraces
      ))],
      facet_args = list(nrow = 1)
    )
    # Born-abroad:UK-born ratio of UK-acquired incidence
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
      heights = c(1, 1, 1, 1, 1, 1)
    ) &
      ggplot2::theme(text = ggplot2::element_text(family = "sans"))
  } else {
    p0 <- patchwork::wrap_plots(
      p1,
      p2,
      ncol = 1,
      heights = c(1, 1)
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
    "overall migration with HIV"
  } else if (quantity == "diag_prob") {
    "HIV diagnosis probabilities, by CD4 strata"
  } else if (quantity == "diagnoses") {
    "HIV diagnoses"
  } else if (quantity == "aids_diagnoses") {
    "AIDS diagnoses"
  } else if (quantity == "cd4_diagnoses") {
    "CD4 diagnoses, by CD4 strata"
  } else if (quantity == "migration_prob") {
    "migration probabilities, by CD4 strata"
  } else if (quantity == "transfers") {
    "transfers of care"
  } else if (quantity == "transfers_detected") {
    "detected transfers of care"
  } else if (quantity == "ratio_abroad_uk") {
    "UK-acquired born-abroad:UK-born ratio"
  } else if (quantity == "detect_prob") {
    "detection probability for transfers of care"
  } else if (quantity == "migrant_hiv_prevalence") {
    "HIV prevalence among migrants"
  } else if (quantity == "incidence_age") {
    "HIV incidence by age group"
  } else if (quantity == "undiag_prev_age") {
    "undiagnosed HIV prevalence by age group"
  } else if (quantity == "undiag_migration_age") {
    "migration with undiagnosed HIV by age group"
  } else if (quantity == "diagnoses_age") {
    "HIV diagnoses by age group"
  } else if (quantity == "aids_diagnoses_age") {
    "AIDS diagnoses by age group"
  } else if (quantity == "cd4_diagnoses_age") {
    "CD4 diagnoses, by CD4 strata and age group"
  } else if (quantity == "transfers_detected_age") {
    "detected transfers of care by age group"
  } else if (quantity == "gamma_obs") {
    "external migration index"
  } else if (quantity == "gamma_obs_group") {
    "external migration index by age group"
  } else if (quantity == "recent_arrival_cd4") {
    "recent-arrival CD4 diagnoses, by CD4 strata and age group"
  } else if (quantity == "diag_prob_age") {
    "HIV diagnosis probabilities, by CD4 strata and age group"
  } else if (quantity == "diag_prob_mig") {
    "HIV diagnosis probabilities (born abroad), by CD4 strata"
  } else if (quantity == "diag_prob_mig_age") {
    "HIV diagnosis probabilities (born abroad), by CD4 strata and age group"
  } else if (quantity == "incidence_trend") {
    "HIV incidence trend by age group"
  } else if (quantity == "peak_incidence") {
    "year of peak HIV incidence by age group"
  } else {
    stop("Invalid quantity specified.")
  }
  return(plot_title)
}

# Internal helper to assign line types to plotted model groups.
make_model_linetypes <- function(groups, chain_results = FALSE) {
  if (length(groups) == 0) {
    return(stats::setNames(character(), character()))
  }

  base_linetypes <- c("solid", "dashed", "dotted", "dotdash", "longdash", "twodash")
  model_groups <- if (chain_results) {
    sub(", Chain [0-9]+$", "", groups)
  } else {
    groups
  }

  model_linetypes <- stats::setNames(
    rep(base_linetypes, length.out = length(unique(model_groups))),
    unique(model_groups)
  )

  stats::setNames(unname(model_linetypes[model_groups]), groups)
}

# Internal helper to assign colours to plotted model groups.
make_model_colors <- function(
  groups,
  chain_results = FALSE,
  has_second_model = FALSE,
  model_1_name = "Model 1",
  model_2_name = "Model 2",
  colours = NULL
) {
  if (length(groups) == 0) {
    return(stats::setNames(character(), character()))
  }

  if (!is.null(colours)) {
    if (!is.character(colours) || length(colours) != length(groups)) {
      stop("`colours` must contain one character value per model or chain.")
    }
    if (!is.null(names(colours)) && !setequal(names(colours), groups)) {
      stop("Named `colours` must match the model or chain legend labels.")
    }
    if (!is.null(names(colours))) {
      return(stats::setNames(unname(colours[groups]), groups))
    }
    return(stats::setNames(colours, groups))
  }

  if (!has_second_model) {
    return(stats::setNames(rep("mediumorchid", length(groups)), groups))
  }

  model_groups <- if (chain_results) {
    sub(", Chain [0-9]+$", "", groups)
  } else {
    groups
  }

  model_colors <- stats::setNames(
    c("mediumorchid", "coral"),
    c(model_1_name, model_2_name)
  )

  stats::setNames(unname(model_colors[model_groups]), groups)
}

# Internal helper to format CD4 facet labels.
label_cd4_strata <- function(labels) {
  ifelse(
    labels == "CD4>=500",
    "CD4>=500",
    encodeString(labels, quote = "'")
  )
  ifelse(
    labels == "CD4<200",
    "CD4<200",
    encodeString(labels, quote = "'")
  )
}

cd4_strata_labeller <- ggplot2::as_labeller(
  label_cd4_strata,
  default = ggplot2::label_parsed
)

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
        "ratio_abroad_uk", "detect_prob", "migration_prob", "migration_prob_age",
        "migrant_hiv_prevalence"
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
  if (quantity %in% c(
    "diag_prob", "diag_prob_mig", "migration_prob", "cd4_diagnoses"
  )) {
    p1 <- p1 +
      ggplot2::facet_grid(
        cols = ggplot2::vars(.data$strata),
        scales = "free_y",
        labeller = ggplot2::labeller(strata = cd4_strata_labeller)
      )
  }

  # facet by age if age-dependent model
  if (
    quantity %in%
      c(
        "incidence_age",
        "undiag_prev_age",
        "undiag_migration_age",
        "diagnoses_age",
        "aids_diagnoses_age",
        "transfers_detected_age",
        "gamma_obs_group"
      )
  ) {
    p1 <- p1 +
      ggplot2::facet_grid(
        cols = ggplot2::vars(.data$age_band),
        scales = "free_y"
      )
  }

  # facet by strata and age for age-dependent diagnosis and migration probabilities
  if (quantity %in% c(
    "diag_prob_age",
    "diag_prob_mig_age",
    "migration_prob_age",
    "cd4_diagnoses_age",
    "recent_arrival_cd4"
  )) {
    p1 <- p1 +
      ggplot2::facet_grid(
        rows = ggplot2::vars(.data$strata),
        cols = ggplot2::vars(.data$age_band),
        scales = "free_y",
        labeller = ggplot2::labeller(strata = cd4_strata_labeller)
      )
  }

  return(p1)
}
