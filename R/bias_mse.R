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
#' @export
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

  if (hiv_list[[1]]$stan_data$config$migration) {
    out$migration_bias <- get_bias(
      hiv_list,
      extract_simulated_quantity(hiv_list[[1]], quantity = "migration"),
      "migration"
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
      title = paste0(
        "Mean Squared Error of estimated ",
        get_plot_title(quantity)
      ),
      x = "",
      y = "Squared difference"
    ) +
    ggplot2::scale_y_continuous(
      breaks = scales::pretty_breaks(),
      labels = scales::label_comma()
    ) +
    ggplot2::theme_bw() +
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
      title = paste("Bias in estimated ", get_plot_title(quantity)),
      x = "",
      y = "Difference"
    ) +
    ggplot2::scale_y_continuous(
      breaks = scales::pretty_breaks(),
      labels = scales::label_comma()
    ) +
    ggplot2::theme_bw() +
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
