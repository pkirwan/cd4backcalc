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
  bias <- dplyr::left_join(
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
#' @return A list of bias and MSE values for incidence, migration, proportion of UK-acquired infections that are migrants, and migration probabilities.
#' @keywords internal
make_mse_list <- function(
  hiv_list
) {
  # make hiv_list a list of lists if not already
  if (!("incidence" %in% names(hiv_list[[1]]))) {
    hiv_list <- list(hiv_list)
  }

  out <- list()

  out$incidence_bias <- get_bias(
    hiv_list,
    extract_simulated_quantity(hiv_list[[1]], quantity = "incidence"),
    "incidence"
  )
  out$incidence_mse <- get_mse(out$incidence_bias)

  if (hiv_list[[1]]$stan_data$model_config$migration) {
    out$undiag_migration_bias <- get_bias(
      hiv_list,
      extract_simulated_quantity(hiv_list[[1]], quantity = "undiag_migration"),
      "undiag_migration"
    )
    out$all_migration_bias <- get_bias(
      hiv_list,
      extract_simulated_quantity(hiv_list[[1]], quantity = "all_migration"),
      "all_migration"
    )
    out$ratio_abroad_uk_bias <- get_bias(
      hiv_list,
      extract_simulated_quantity(hiv_list[[1]], quantity = "ratio_abroad_uk"),
      "ratio_abroad_uk"
    )
    out$mp_bias <- get_bias_strata(
      hiv_list,
      extract_simulated_quantity(hiv_list[[1]], quantity = "migration_prob"),
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
#' @param hiv_list A list returned from the `postproc` function containing model estimates.
#' @param hiv_list_2 A second list returned from the `postproc` function for comparison (default is NULL).
#' @param quantity The quantity to plot MSE for (default is "incidence", other options are "undiag_migration", "ratio_abroad_uk", and "mp").
#' @param offset The offset for text labels.
#' @param lab A logical indicating whether to add labels.
#' @param model_1_name The name of the first model for labeling (default is "Model 1").
#' @param model_2_name The name of the second model for labeling (default is "Model 2").
#' @return A ggplot object of the MSE box plot.
#' @export
mse_plot <- function(
  hiv_list,
  hiv_list_2 = NULL,
  quantity = "incidence",
  offset = 0.1,
  lab = TRUE,
  model_1_name = "Model 1",
  model_2_name = "Model 2"
) {
  data <- make_mse_list(
    hiv_list
  )[[paste0(quantity, "_mse")]]

  if (!is.null(hiv_list_2)) {
    data_2 <- make_mse_list(
      hiv_list_2
    )[[paste0(quantity, "_mse")]]

    data <- dplyr::bind_rows(
      data |> dplyr::mutate(group = model_1_name),
      data_2 |> dplyr::mutate(group = model_2_name)
    )
  }

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

  if (!is.null(hiv_list_2)) {
    p1 <- p1 +
      ggplot2::scale_color_manual(
        values = c("mediumorchid3", "coral1"),
        name = "Model",
        labels = c(model_1_name, model_2_name)
      )
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
        nudge_x = 0.3,
        size = 4
      )
  }
  return(p1)
}

#' Box plots of bias.
#'
#' @param hiv_list A list returned from the `postproc` function containing model estimates.
#' @param quantity The quantity to plot bias for (default is "incidence", other options are "undiag_migration", "ratio_abroad_uk", and "mp").
#' @param offset The offset for text labels.
#' @param lab A logical indicating whether to add labels.
#' @param hiv_list_2 A second list returned from the `postproc` function for comparison (default is NULL).
#' @param model_1_name The name of the first model for labeling (default is "Model 1").
#' @param model_2_name The name of the second model for labeling (default is "Model 2").
#' @return A ggplot object of the bias box plot.
#' @export
bias_plot <- function(
  hiv_list,
  quantity = "incidence",
  offset = 0.1,
  lab = TRUE,
  hiv_list_2 = NULL,
  model_1_name = "Model 1",
  model_2_name = "Model 2"
) {
  data <- make_mse_list(
    hiv_list
  )[[paste0(quantity, "_bias")]]

  if (!is.null(hiv_list_2)) {
    data_2 <- make_mse_list(
      hiv_list_2
    )[[paste0(quantity, "_bias")]]

    data <- dplyr::bind_rows(
      data |> dplyr::mutate(group = model_1_name),
      data_2 |> dplyr::mutate(group = model_2_name)
    )
  }

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

  if (!is.null(hiv_list_2)) {
    p1 <- p1 +
      ggplot2::scale_color_manual(
        values = c("mediumorchid3", "coral1"),
        name = "Model",
        labels = c(model_1_name, model_2_name)
      )
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
