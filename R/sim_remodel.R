#' Remodel the simulated diagnosis data for fitting a different model type.
#'
#' This function takes the output from \code{simulate_diagnoses()} and
#' remodels the data to be suitable for fitting a different model type
#' (e.g., from RITA to CD4-only model).
#'
#' @param sim_diags A list output from \code{simulate_diagnoses()}.
#' @param rita Logical; if \code{TRUE}, remodels the data for fitting an RITA
#' model; if \code{FALSE}, remodels the data for fitting a CD4-only model.
#' @return A list containing the remodeled stan data and configuration.
#' @export
#' @examples
#' #' # Simulate diagnoses data using RITA model
#' sim_diags <- simulate_diagnoses(nquar = 60, rita = TRUE, sim_seed = 123)
#' # Remodel data for fitting CD4-only model
#' remodeled_data <- sim_remodel(sim_diags, rita = FALSE)
sim_remodel <- function(sim_diags, rita = TRUE) {
  # Extract simulated stan data and configuration
  stan_data <- sim_diags$stan_data
  sim_data <- sim_diags$sim_data

  # Sum across RITA states if remodeling to CD4-only model
  if (!rita) {
    # Remodel CD4 data

    # recently acquired have CD4 strata with proportions 0.59, 0.27, 0.11, 0.03
    cd4_dx <- matrix(0, nrow = stan_data$nquar, ncol = 4)
    cd4_dx[, 1] <- stan_data$CD4[, 2] + stan_data$CD4[, 1] * 0.59
    cd4_dx[, 2] <- stan_data$CD4[, 3] + stan_data$CD4[, 1] * 0.27
    cd4_dx[, 3] <- stan_data$CD4[, 4] + stan_data$CD4[, 1] * 0.11
    cd4_dx[, 4] <- stan_data$CD4[, 5] + stan_data$CD4[, 1] * 0.03

    cd4_dx <- round_matrix_preserve(cd4_dx, rowwise = TRUE)

    cd4_dx_age <- array(0, dim = c(stan_data$nquar, stan_data$nage, 4))
    cd4_dx_age[,, 1] <- stan_data$CD4_age[,, 2] +
      stan_data$CD4_age[,, 1] * 0.59
    cd4_dx_age[,, 2] <- stan_data$CD4_age[,, 3] +
      stan_data$CD4_age[,, 1] * 0.27
    cd4_dx_age[,, 3] <- stan_data$CD4_age[,, 4] +
      stan_data$CD4_age[,, 1] * 0.11
    cd4_dx_age[,, 4] <- stan_data$CD4_age[,, 5] +
      stan_data$CD4_age[,, 1] * 0.03

    # round_matrix_preserve applied for each quarter
    for (t in 1:stan_data$nquar) {
      cd4_dx_age[t, , ] <- round_matrix_preserve(
        cd4_dx_age[t, , ],
        rowwise = TRUE
      )
    }

    # Remodel initial prevalence data
    init_prev <- c(
      stan_data$init_prev[1] + stan_data$init_prev[2] + stan_data$init_prev[4],
      stan_data$init_prev[3] + stan_data$init_prev[5],
      stan_data$init_prev[6],
      stan_data$init_prev[7]
    )

    init_prev_abroad <- stan_data$init_prev_abroad[4:7]

    init_prev_age <- matrix(0, nrow = 4, ncol = stan_data$nage)
    init_prev_age[1, ] <-
      stan_data$init_prev_age[1, ] +
      stan_data$init_prev_age[2, ] +
      stan_data$init_prev_age[4, ]
    init_prev_age[2, ] <-
      stan_data$init_prev_age[3, ] +
      stan_data$init_prev_age[5, ]
    init_prev_age[3, ] <- stan_data$init_prev_age[6, ]
    init_prev_age[4, ] <- stan_data$init_prev_age[7, ]

    init_prev_abroad_age <- stan_data$init_prev_abroad_age[4:7, ]

    # Update simulated stan data
    stan_data$CD4 <- cd4_dx
    stan_data$CD4_age <- cd4_dx_age
    stan_data$init_prev <- init_prev
    stan_data$init_prev_age <- init_prev_age
    stan_data$init_prev_abroad <- init_prev_abroad
    stan_data$init_prev_abroad_age <- init_prev_abroad_age
  }

  return(list(
    "stan_data" = stan_data,
    "sim_data" = sim_data
  ))
}
