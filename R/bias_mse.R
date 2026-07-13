#' Get bias between HIV model estimates and simulation results.
#'
#' @param hiv_list A list of HIV model results.
#' @param type The variable type to calculate bias for.
#' @param deterministic Logical indicating whether to compare against
#'   deterministic generative values.
#' @return A data frame of bias values.
get_bias <- function(
  hiv_list,
  type,
  deterministic = FALSE,
  years = NULL
) {
  purrr::imap(hiv_list, function(x, simulation) {
    estimates <- get_estimates(x, type)
    if (!("mean" %in% names(estimates))) {
      stop(sprintf(
        "Estimate '%s' must contain a 'mean' column for bias/MSE calculations.",
        type
      ))
    }

    join_cols <- intersect(c("year", "strata", "age_band"), names(estimates))
    joined <- estimates |>
      dplyr::select(dplyr::all_of(c(join_cols, "mean"))) |>
      dplyr::left_join(
        extract_simulated_quantity(
          x,
          quantity = type,
          deterministic = deterministic
        ) |>
          dplyr::select(dplyr::all_of(c(join_cols, "value"))),
        by = join_cols
      )

    if (any(is.na(joined$value))) {
      stop(sprintf(
        "Failed to match simulated values to '%s' estimates when calculating bias/MSE.",
        type
      ))
    }

    joined |>
      dplyr::filter(
        .data$year > 0,
        if (is.null(years)) TRUE else .data$year %in% years
      ) |>
      dplyr::mutate(
        simulation = as.character(simulation),
        bias = .data$mean - .data$value
      ) |>
      dplyr::select(dplyr::all_of(join_cols), "simulation", "bias")
  }) |>
    purrr::list_rbind()
}

#' Get empirical bias between HIV model estimates and simulation results.
#'
#' @param bias_df A data frame of raw bias values.
#' @return A data frame of empirical bias values averaged across simulations.
#' @keywords internal
summarise_bias <- function(bias_df) {
  group_cols <- intersect(c("year", "strata", "age_band"), names(bias_df))
  bias_df |>
    dplyr::group_by(dplyr::across(dplyr::all_of(group_cols))) |>
    dplyr::summarise(bias = mean(.data$bias), .groups = "drop")
}

#' Get mean squared error (MSE) list for HIV model estimates.
#'
#' @param bias_df A data frame of raw bias values.
#' @return A data frame of MSE values.
get_mse <- function(bias_df) {
  bias_df |>
    dplyr::group_by(.data$simulation) |>
    dplyr::summarise(
      mse = mean(.data$bias^2),
      n = dplyr::n(),
      .groups = "drop"
    )
}

format_bias_summary <- function(value) {
  if (abs(value) >= 1000) {
    format(round(value, 0), big.mark = ",")
  } else if (abs(value) >= 1) {
    format(round(value, 1), nsmall = 1)
  } else {
    format(round(value, 5), nsmall = 5, scientific = FALSE)
  }
}

add_bias_outputs <- function(
  out,
  out_name,
  hiv_list,
  quantity,
  deterministic,
  strip_total = FALSE,
  years = NULL
) {
  raw_bias <- get_bias(
    hiv_list,
    type = quantity,
    deterministic = deterministic,
    years = years
  )

  # Remove the aggregate "all" row when age-band-specific summaries are wanted,
  # to avoid double-counting the total in bias/MSE calculations.
  if (strip_total && "age_band" %in% names(raw_bias)) {
    raw_bias <- dplyr::filter(raw_bias, .data$age_band != "all")
  }

  out[[paste0(out_name, "_bias")]] <- summarise_bias(raw_bias)
  out[[paste0(out_name, "_mse")]] <- get_mse(raw_bias)
  out
}

.model_flag_enabled <- function(hiv_result, flag_name) {
  isTRUE(hiv_result$stan_data$model_config[[flag_name]]) ||
    isTRUE(hiv_result$stan_data$sim_config[[flag_name]])
}

#' Make a list of mean squared error (MSE) values for HIV model estimates.
#'
#' @param hiv_list A list of HIV model results.
#' @param sim_comparison A string specifying which simulated values to compare against. Options are `"sampled"` (the stochastic realisation for the particular seed, the default) or `"generative"` (deterministic expected values).
#' @param years Optional numeric/integer vector specifying which years to
#'   include in bias/MSE calculations. If `NULL` (default), all available years
#'   are used.
#' @return A list of bias and MSE values for incidence, migration, proportion of UK-acquired infections that are migrants, and migration probabilities.
#' @keywords internal
make_mse_list <- function(
  hiv_list,
  sim_comparison = c("sampled", "generative"),
  years = NULL
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

  out <- add_bias_outputs(
    out,
    out_name = "incidence",
    hiv_list = hiv_list,
    quantity = "incidence",
    deterministic = use_deterministic,
    years = years
  )

  if (.model_flag_enabled(hiv_list[[1]], "age")) {
    out <- add_bias_outputs(
      out,
      out_name = "incidence_age",
      hiv_list = hiv_list,
      quantity = "incidence_age",
      deterministic = use_deterministic,
      strip_total = TRUE,
      years = years
    )
  }

  if (.model_flag_enabled(hiv_list[[1]], "migration")) {
    out <- add_bias_outputs(
      out,
      out_name = "undiag_migration",
      hiv_list = hiv_list,
      quantity = "undiag_migration",
      deterministic = use_deterministic,
      years = years
    )
    out <- add_bias_outputs(
      out,
      out_name = "all_migration",
      hiv_list = hiv_list,
      quantity = "all_migration",
      deterministic = use_deterministic,
      years = years
    )
    out <- add_bias_outputs(
      out,
      out_name = "ratio_abroad_uk",
      hiv_list = hiv_list,
      quantity = "ratio_abroad_uk",
      deterministic = use_deterministic,
      years = years
    )
    out <- add_bias_outputs(
      out,
      out_name = "mp",
      hiv_list = hiv_list,
      quantity = "migration_prob",
      deterministic = use_deterministic,
      years = years
    )
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
#' @param years Optional numeric/integer vector specifying which years to
#'   include in bias/MSE calculations. If `NULL` (default), all available years
#'   are used.
#' @return A ggplot object of the MSE box plot.
#' @examples
#' # See ?run_backcalc for a full end-to-end workflow.
#' invisible(NULL)
#' @export
mse_plot <- function(
  hiv_lists,
  quantity = "incidence",
  offset = 0.1,
  lab = TRUE,
  colours = NULL,
  sim_comparison = c("sampled", "generative"),
  years = NULL
) {
  sim_comparison <- match.arg(sim_comparison)

  if (is.null(names(hiv_lists))) {
    names(hiv_lists) <- paste0("Model ", seq_along(hiv_lists))
  }

  data <- purrr::imap(hiv_lists, function(hl, name) {
    make_mse_list(hl, sim_comparison = sim_comparison, years = years)[[paste0(quantity, "_mse")]] |>
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
            pmse = stats::weighted.mean(.data$mse, .data$n),
            pos = stats::quantile(.data$mse, 0.75) + offset,
            formatted_pmse = format_bias_summary(.data$pmse)
          ),
        ggplot2::aes(
          y = .data$pos,
          label = paste0("PMSE: ", .data$formatted_pmse),
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
#' @param years Optional numeric/integer vector specifying which years to
#'   include in bias/MSE calculations. If `NULL` (default), all available years
#'   are used.
#' @return A ggplot object of the bias box plot.
#' @examples
#' # See ?run_backcalc for a full end-to-end workflow.
#' invisible(NULL)
#' @export
bias_plot <- function(
  hiv_lists,
  quantity = "incidence",
  offset = 0.1,
  lab = TRUE,
  colours = NULL,
  sim_comparison = c("sampled", "generative"),
  years = NULL
) {
  sim_comparison <- match.arg(sim_comparison)

  if (is.null(names(hiv_lists))) {
    names(hiv_lists) <- paste0("Model ", seq_along(hiv_lists))
  }

  data <- purrr::imap(hiv_lists, function(hl, name) {
    make_mse_list(hl, sim_comparison = sim_comparison, years = years)[[paste0(quantity, "_bias")]] |>
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
            formatted_bias = format_bias_summary(.data$meanbias)
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
