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
  dplyr::left_join(
    get_estimates(hiv_list, type),
    model_estimates,
    by = c("year")
  ) |>
    # bias should not be calculated for year 0
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
  dplyr::left_join(
    get_estimates(hiv_list, type),
    model_estimates,
    by = c("year", "strata")
  ) |>
    # bias should not be calculated for year 0
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
#' @param sim_comparison A string specifying which simulated values to compare against. Options are `"sampled"` (the stochastic realisation for the particular seed, the default) or `"generative"` (deterministic expected values).
#' @return A list of bias and MSE values for incidence, migration, proportion of UK-acquired infections that are migrants, and migration probabilities.
#' @keywords internal
make_mse_list <- function(
  hiv_list,
  sim_comparison = c("sampled", "generative")
) {
  sim_comparison <- match.arg(sim_comparison)
  use_deterministic <- sim_comparison == "generative"
  hiv_list <- purrr::compact(hiv_list)
  if (length(hiv_list) == 0) {
    stop("All entries in hiv_list are NULL.")
  }
  # make hiv_list a list of lists if not already
  if (!("incidence" %in% names(hiv_list[[1]]))) {
    hiv_list <- list(hiv_list)
  }

  out <- list()

  out$incidence_bias <- get_bias(
    hiv_list,
    extract_simulated_quantity(
      hiv_list[[1]],
      quantity = "incidence",
      deterministic = use_deterministic
    ),
    "incidence"
  )
  out$incidence_mse <- get_mse(out$incidence_bias)

  if (hiv_list[[1]]$stan_data$model_config$migration) {
    out$undiag_migration_bias <- get_bias(
      hiv_list,
      extract_simulated_quantity(
        hiv_list[[1]],
        quantity = "undiag_migration",
        deterministic = use_deterministic
      ),
      "undiag_migration"
    )
    out$all_migration_bias <- get_bias(
      hiv_list,
      extract_simulated_quantity(
        hiv_list[[1]],
        quantity = "all_migration",
        deterministic = use_deterministic
      ),
      "all_migration"
    )
    out$ratio_abroad_uk_bias <- get_bias(
      hiv_list,
      extract_simulated_quantity(
        hiv_list[[1]],
        quantity = "ratio_abroad_uk",
        deterministic = use_deterministic
      ),
      "ratio_abroad_uk"
    )
    out$mp_bias <- get_bias_strata(
      hiv_list,
      extract_simulated_quantity(
        hiv_list[[1]],
        quantity = "migration_prob",
        deterministic = use_deterministic
      ),
      "migration_prob"
    )

    out$undiag_migration_mse <- get_mse(out$undiag_migration_bias)
    out$all_migration_mse <- get_mse(out$all_migration_bias)
    out$ratio_abroad_uk_mse <- get_mse(out$ratio_abroad_uk_bias)
    out$mp_mse <- get_mse(out$mp_bias)
  }
  return(out)
}

#' Box plots of mean squared error (MSE).
#'
#' @param hiv_lists A named list of model result lists, each returned from the
#'   `postproc` function. Names are used as model labels. If unnamed, models are
#'   labelled "Model 1", "Model 2", etc.
#' @param quantity The quantity to plot MSE for (default is `"incidence"`).
#'   Migration-model options are `"undiag_migration"`, `"ratio_abroad_uk"`, and
#'   `"mp"`.
#' @param offset The offset for text labels.
#' @param lab A logical indicating whether to add labels.
#' @param colours An optional named character vector of colours for each model.
#'   Names must match the names of `hiv_lists`. If `NULL` (default), ggplot's
#'   default discrete colour palette is used.
#' @param sim_comparison A string specifying which simulated values to compare
#'   against. Options are `"sampled"` (the stochastic realisation for the
#'   particular seed, the default) or `"generative"` (deterministic expected
#'   values).
#' @return A ggplot object of the MSE box plot.
#' @export
mse_plot <- function(
  hiv_lists,
  quantity = "incidence",
  offset = 0.1,
  lab = TRUE,
  colours = NULL,
  sim_comparison = c("sampled", "generative")
) {
  sim_comparison <- match.arg(sim_comparison)

  if (is.null(names(hiv_lists))) {
    names(hiv_lists) <- paste0("Model ", seq_along(hiv_lists))
  }

  data <- purrr::imap(hiv_lists, function(hl, name) {
    make_mse_list(hl, sim_comparison = sim_comparison)[[paste0(quantity, "_mse")]] |>
      dplyr::mutate(group = name)
  }) |>
    dplyr::bind_rows()

  # create box plot
  p1 <- data |>
    ggplot2::ggplot() +
    ggplot2::aes(x = .data$group, y = .data$mse, color = .data$group) +
    ggplot2::geom_boxplot(
      coef = 1.5,
      alpha = 0.2,
      width = 0.4,
      outliers = FALSE
    ) +
    ggplot2::labs(
      title = paste(
        "Mean Squared Error of estimated",
        get_plot_title(quantity)
      ),
      x = "",
      y = "Squared difference"
    ) +
    ggplot2::scale_y_continuous(
      breaks = scales::pretty_breaks(),
      labels = scales::label_comma()
    ) +
    ggplot2::theme_classic() +
    ggplot2::theme(
      legend.position = "none"
    )

  if (!is.null(colours)) {
    p1 <- p1 + ggplot2::scale_color_manual(values = colours)
  }

  if (lab) {
    p1 <- p1 +
      ggplot2::geom_text(
        data = data |>
          dplyr::group_by(.data$group) |>
          dplyr::summarise(
            bias = mean(.data$mse),
            pos = stats::quantile(.data$mse, 0.75) + offset,
            # Format bias based on its magnitude
            formatted_bias = if (abs(.data$bias) >= 1000) {
              format(
                round(.data$bias, 0),
                big.mark = ","
              )
            } else if (abs(.data$bias) >= 1) {
              format(round(.data$bias, 1), nsmall = 1)
            } else {
              format(round(.data$bias, 5), nsmall = 5, scientific = FALSE)
            }
          ),
        ggplot2::aes(
          y = .data$pos,
          label = paste0("PMSE: ", .data$formatted_bias),
        ),
        nudge_x = 0.4,
        size = 4
      )
  }
  return(p1)
}

#' Box plots of bias.
#'
#' @param hiv_lists A named list of model result lists, each returned from the
#'   `postproc` function. Names are used as model labels. If unnamed, models are
#'   labelled "Model 1", "Model 2", etc.
#' @param quantity The quantity to plot bias for (default is `"incidence"`).
#'   Migration-model options are `"undiag_migration"`, `"ratio_abroad_uk"`, and
#'   `"mp"`.
#' @param offset The offset for text labels.
#' @param lab A logical indicating whether to add labels.
#' @param colours An optional named character vector of colours for each model.
#'   Names must match the names of `hiv_lists`. If `NULL` (default), ggplot's
#'   default discrete colour palette is used.
#' @param sim_comparison A string specifying which simulated values to compare
#'   against. Options are `"sampled"` (the stochastic realisation for the
#'   particular seed, the default) or `"generative"` (deterministic expected
#'   values).
#' @return A ggplot object of the bias box plot.
#' @export
bias_plot <- function(
  hiv_lists,
  quantity = "incidence",
  offset = 0.1,
  lab = TRUE,
  colours = NULL,
  sim_comparison = c("sampled", "generative")
) {
  sim_comparison <- match.arg(sim_comparison)

  if (is.null(names(hiv_lists))) {
    names(hiv_lists) <- paste0("Model ", seq_along(hiv_lists))
  }

  data <- purrr::imap(hiv_lists, function(hl, name) {
    make_mse_list(hl, sim_comparison = sim_comparison)[[paste0(quantity, "_bias")]] |>
      dplyr::mutate(group = name)
  }) |>
    dplyr::bind_rows()

  p1 <- data |>
    ggplot2::ggplot() +
    ggplot2::aes(x = .data$group, y = .data$bias, color = .data$group) +
    ggplot2::geom_hline(yintercept = 0, colour = "grey") +
    ggplot2::geom_boxplot(coef = Inf, alpha = 0.2, width = 0.4) +
    ggplot2::labs(
      title = paste("Bias in estimated", get_plot_title(quantity)),
      x = "",
      y = "Difference"
    ) +
    ggplot2::scale_y_continuous(
      breaks = scales::pretty_breaks(),
      labels = scales::label_comma()
    ) +
    ggplot2::theme_classic() +
    ggplot2::theme(
      legend.position = "none"
    )

  if (!is.null(colours)) {
    p1 <- p1 + ggplot2::scale_color_manual(values = colours)
  }

  if (lab) {
    p1 <- p1 +
      ggplot2::geom_text(
        data = data |>
          dplyr::group_by(.data$group) |>
          dplyr::summarise(
            meanbias = mean(.data$bias),
            pos = stats::quantile(.data$bias, 0.75) + offset,
            # Format bias based on its magnitude
            formatted_bias = if (abs(.data$meanbias) >= 1000) {
              format(
                round(.data$meanbias, 0),
                big.mark = ","
              )
            } else if (abs(.data$meanbias) >= 1) {
              format(round(.data$meanbias, 1), nsmall = 1)
            } else {
              format(round(.data$meanbias, 5), nsmall = 5, scientific = FALSE)
            }
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
