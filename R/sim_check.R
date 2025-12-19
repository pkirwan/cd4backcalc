
#' Validate simulation
#'
#' Validates the consistency of the simulated epidemic data.
#'
#' @param sim_data List of simulated epidemic data.
#' @return NULL. Throws an error if validation checks fail.
validate_sim_data <- function(sim_data) {

  total_uk_acq <- sim_data$infs[, 1]
  
  uk_born_uk_acq <- total_uk_acq / (1 + sim_data$ratio_abroad_uk)
  abroad_born_uk_acq <- total_uk_acq * sim_data$ratio_abroad_uk / (1 + sim_data$ratio_abroad_uk)

  stopifnot(
    "UK epidemic totals not preserved" = all(
      abs(
        rowSums(sim_data$epi_uk$lat) +
          cumsum(rowSums(sim_data$epi_uk$diags)) -
          sum(sim_data$lat_init_uk) -
          cumsum(uk_born_uk_acq)
      ) <
        15
    )
  )

  stopifnot(
    "UK abroad epidemic totals not preserved" = all(
      abs(
        rowSums(sim_data$epi_uk_abroad$lat) +
          cumsum(rowSums(sim_data$epi_uk_abroad$diags)) -
          sum(sim_data$lat_init_uk_abroad) -
          cumsum(abroad_born_uk_acq)
      ) <
        15
    )
  )

  stopifnot(
    "Abroad epidemic totals not preserved" = all(
      abs(
        rowSums(sim_data$epi_abroad$lat) +
          cumsum(rowSums(sim_data$epi_abroad$diags)) -
          sum(sim_data$lat_init_abroad) -
          cumsum(rowSums(sim_data$m[, 1:(sim_data$config$nlatst)]))
      ) <
        15
    )
  )

  return(NULL)
}

#' Validate Stan data
#'
#' Validates the consistency of the Stan data with the simulated epidemic data.
#'
#' @param stan_data List of data prepared for Stan model fitting.
#' @param sim_data List of simulated epidemic data.
#' @param config Configuration list defining simulation parameters.
#' @return NULL. Throws an error if validation checks fail.
validate_stan_data <- function(stan_data, sim_data, config) {
  if (config$migration) {
    stopifnot(
      "HIV diags different to sum of CD4 states" = sum(rowSums(
        stan_data$HIV
      )) ==
        sum(rowSums(stan_data$CD4))
    )
  } else {
    stopifnot(
      "HIV diags different to sum of CD4 states" = sum(stan_data$HIV) ==
        sum(rowSums(stan_data$CD4))
    )
  }

  stopifnot(
    "Infections not preserved to diagnoses" = abs(
      sum(sim_data$infs[,1]) +
        sum(stan_data$init_prev) +
        sum(stan_data$init_prev_abroad) +
        # abroad infections migrating to UK
        sum(sim_data$m[, 1:config$nlatst]) -
        # total diagnoses plus latent infections
        sum(stan_data$HIV) -
        sum(stan_data$AIDS) -
        sum(sim_data$epi_uk$lat[config$nquar, ]) -
        sum(sim_data$epi_uk_abroad$lat[config$nquar, ]) -
        sum(sim_data$epi_abroad$lat[config$nquar, ])
    ) <
      15
  )
  return(NULL)
}
