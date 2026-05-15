#' Define simulation configuration
#'
#' Defines the patterns for various parameters in the epidemic simulation.
#' The number of latent and diagnosis states will be inferred from the simulation specification.
#'
#' @param h_pattern Pattern type for infections ("constant", "increasing", "varying", "varying_2").
#' @param o_pattern Pattern type for omega ("constant", "increasing", "step_change_1", "step_change_2", "none").
#' @param d_pattern Pattern type for diagnosis probabilities ("constant", "varying").
#' @param p_pattern Pattern type for proportion UK ("constant", "increasing").
#' @param m_pattern Pattern type for mu ("constant", "differing", "varying").
#' @param r_pattern Pattern type for RITA data ("constant", "increasing", "none").
#' @param a_pattern Pattern type for age-based data ("constant", "wedge", "varying_wedge").
#' @param sim_type If provided, a preset type of simulation to define patterns for.
#' @param sim_seed Random seed for simulation. If NULL, a random seed is generated.
#' @param dist Logical indicating whether to add stochastic variation.
#' @param migration Logical indicating whether to include migration effects.
#' @param rita Logical indicating whether to simulate RITA (Recent Infection Testing Algorithm) data.
#' @param age Logical indicating whether to simulate age-stratified data.
#' @param nquar Number of time quarters.
#' @param nage Number of age groups.
#' @param q Vector of progression probabilities.
#' @param target_inf Target total number of HIV acquisitions.
#' @param target_mig Target total number of migrations.
#' @param init_props Initial state proportions.
#' @param init_props_abroad Initial state proportions for infections acquired abroad.
#' @param abroad_prop Proportion of infections acquired abroad at the outset.
#' @param cd4_completeness CD4 completeness.
#' @param rita_completeness RITA completeness.
#' @param age_diag_effect Numeric scalar controlling the strength of age variation in diagnosis probabilities. Applied as a logit-additive offset scaled linearly across ages.
#' @param mig_diag_offset Numeric vector of per-state logit-scale offsets applied to diagnosis probabilities for the migrant group (y=2). Length must equal the number of diagnosis states (4 for non-RITA, 5 for RITA). NULL (default) applies no offset.
#' @param frr False recency rate for RITA misclassification.
#' @param detect_obs_every Positive integer giving the spacing, in quarters, between external detection-probability observations used for simulation diagnostics.
#' @return Configuration list.
#' @keywords internal
define_sim_config <- function(
  h_pattern = "constant",
  o_pattern = "constant",
  d_pattern = "constant",
  p_pattern = "constant",
  m_pattern = "constant",
  r_pattern = "constant",
  a_pattern = "constant",
  sim_type = NULL,
  sim_seed = NULL,
  dist = TRUE,
  migration = FALSE,
  rita = FALSE,
  age = FALSE,
  nquar = 120, # number of quarters
  nage = 52, # number of age groups
  q = c(0.0976, 0.115, 0.116, 0.148), # progression probabilities
  target_inf = 90000, # target total number of HIV acquisitions
  target_mig = 12000, # target total number of migrations
  init_props = c(7, 3.5, 1.75, 0.5), # initial state proportions
  init_props_abroad = c(1.75, 2.25, 2.25, 1.25), # initial state proportions abroad
  abroad_prop = 0.2, # proportion acquired abroad at outset
  cd4_completeness = 1, # CD4 completeness
  rita_completeness = 1, # RITA completeness
  age_diag_effect = 0, # age effect on diagnosis probabilities
  mig_diag_offset = NULL, # per-state logit offsets for migrant diagnosis probs
  frr = 0.019, # false recency rate for RITA misclassification
  detect_obs_every = 1 # spacing between detect_prob observations
) {
  if (rita) {
    nlatst <- 7
    ndiagst <- 6
    init_props <- c(1.2, 1.0, 0.1, 4.8, 3.4, 1.75, 0.5)
    init_props_abroad <- c(0, 0, 0, 1.75, 2.25, 2.25, 1.25)
  } else {
    nlatst <- 4
    ndiagst <- 5
  }

  q_age <- matrix(rep(q, nage), nrow = nage, byrow = TRUE)

  if (is.null(sim_seed)) {
    sim_seed <- sample.int(.Machine$integer.max, 1)
  }

  if (
    length(detect_obs_every) != 1 ||
      is.na(detect_obs_every) ||
      detect_obs_every < 1 ||
      detect_obs_every != as.integer(detect_obs_every)
  ) {
    stop("detect_obs_every must be a positive integer.")
  }

  # Standard configuration parameters
  config <- list(
    h_pattern = h_pattern,
    o_pattern = o_pattern,
    d_pattern = d_pattern,
    p_pattern = p_pattern,
    m_pattern = m_pattern,
    r_pattern = r_pattern,
    a_pattern = a_pattern,
    sim_seed = sim_seed,
    dist = dist,
    migration = migration,
    rita = rita,
    age = age,
    nquar = nquar,
    nage = nage,
    nlatst = nlatst,
    ndiagst = ndiagst,
    q = q,
    q_age = q_age,
    target_inf = target_inf,
    target_mig = target_mig,
    init_props = init_props,
    init_props_abroad = init_props_abroad,
    abroad_prop = abroad_prop,
    cd4_completeness = cd4_completeness,
    rita_completeness = rita_completeness,
    age_diag_effect = age_diag_effect,
    mig_diag_offset = mig_diag_offset,
    frr = frr,
    detect_obs_every = as.integer(detect_obs_every)
  )

  if (!is.null(sim_type)) {
    config <- define_simulation_type(config, sim_type)
  }

  if (migration == FALSE) {
    config$o_pattern <- "none"
  }

  if (rita == FALSE) {
    config$r_pattern <- "none"
  }

  return(config)
}

#' Define simulation type
#'
#' This is a helper function to define various pre-specified simulation
#' types for epidemic simulations. Each configuration specifies patterns
#' for infections, migration, diagnosis probabilities, and other parameters.
#'
#' List of available simulation types:
#' - "constant": All parameters constant over time.
#' - "h_increasing": Infections increase over time.
#' - "h_varying": Infections vary over time.
#' - "o_increasing": Migration (omega) increases over time.
#' - "o_step_change": Migration (omega) has step changes over time.
#' - "d_varying": Diagnosis probabilities vary over time.
#' - "p_increasing": Proportion UK increases over time.
#' - "m_differing": Migration probabilities (mu) differ by diagnosed state.
#' - "m_varying": Migration probabilities (mu) vary over time.
#' - "combo_1": Combination of varying infections, step change migration, and varying diagnosis probabilities.
#' - "combo_2": Combination of varying infections, step change migration, varying diagnosis probabilities, and increasing proportion UK.
#' - "combo_3": Combination of early-peaking infections, increasing migration, varying diagnosis probabilities, increasing proportion UK, varying migration probabilities, and age effects.
#' - "combo_4": Late-peaking infections (varying_2) with a late upward step change in migration (step_change_2),
#'   representing a surge in arrivals from high-prevalence regions. Shares varying diagnosis probabilities,
#'   increasing proportion UK, varying migration probabilities, and age effects with combo_3.
#'
#' @param config Configuration list to modify.
#' @param sim_type Type of simulation to define patterns for.
#' @return Configuration list for the specified simulation type.
#' @keywords internal
define_simulation_type <- function(config, sim_type = "constant") {
  # check sim_type is valid
  valid_sim_types <- c(
    "constant",
    "h_increasing",
    "h_varying",
    "o_increasing",
    "o_step_change",
    "d_varying",
    "p_increasing",
    "m_differing",
    "m_varying",
    "combo_1",
    "combo_2",
    "combo_3",
    "combo_4"
  )
  if (!sim_type %in% valid_sim_types) {
    stop(
      "Invalid sim_type. Choose from: ",
      paste(valid_sim_types, collapse = ", ")
    )
  }

  if (sim_type == "h_increasing") {
    config$h_pattern <- "increasing"
  } else if (sim_type == "h_varying") {
    config$h_pattern <- "varying"
  } else if (sim_type == "o_increasing") {
    config$o_pattern <- "increasing"
  } else if (sim_type == "o_step_change") {
    config$o_pattern <- "step_change_1"
  } else if (sim_type == "d_varying") {
    config$d_pattern <- "varying"
  } else if (sim_type == "p_increasing") {
    config$p_pattern <- "increasing"
  } else if (sim_type == "m_differing") {
    config$m_pattern <- "differing"
  } else if (sim_type == "m_varying") {
    config$m_pattern <- "varying"
  } else if (sim_type == "combo_1") {
    config$h_pattern <- "varying"
    config$o_pattern <- "step_change"
    config$d_pattern <- "varying"
  } else if (sim_type == "combo_2") {
    config$h_pattern <- "varying"
    config$o_pattern <- "step_change"
    config$d_pattern <- "varying"
    config$p_pattern <- "increasing"
  } else if (sim_type == "combo_3") {
    config$h_pattern <- "varying"
    config$o_pattern <- "increasing"
    config$d_pattern <- "varying"
    config$p_pattern <- "increasing"
    config$m_pattern <- "varying"
    config$a_pattern <- "wedge"
    if (config$age) {
      config$age_diag_effect <- 0.5
    }
    if (config$migration) {
      if (config$rita) {
        config$mig_diag_offset <- c(-0.5, -0.5, -0.25, -0.2, -0.1)
      } else {
        config$mig_diag_offset <- c(-0.5, -0.25, -0.2, -0.1)
      }
    }
  } else if (sim_type == "combo_4") {
    config$h_pattern <- "varying_2"
    config$o_pattern <- "step_change_2"
    config$d_pattern <- "varying"
    config$p_pattern <- "increasing"
    config$m_pattern <- "varying"
    config$a_pattern <- "wedge"
    if (config$age) {
      config$age_diag_effect <- 0.5
    }
    if (config$migration) {
      if (config$rita) {
        config$mig_diag_offset <- c(-0.5, -0.5, -0.25, -0.2, -0.1)
      } else {
        config$mig_diag_offset <- c(-0.5, -0.25, -0.2, -0.1)
      }
    }
  }

  return(config)
}

#' Generate Q Matrix
#'
#' Generates the Q matrix for the epidemic model.
#' CD4-only latent states definition:
#' Latent state 1: CD4 > 500
#' Latent state 2: 500 < CD4 < 350
#' Latent state 3: 200 < CD4 < 350
#' Latent state 4: CD4 < 200
#'
#' RITA latent states definition:
#' Latent state 1: CD4 > 500, HIV acquired in the last three months
#' Latent state 2: CD4 > 500, HIV acquired between 3 and 6 months ago
#' Latent state 3: 500 < CD4 < 350, HIV acquired between 3 and 6 months ago
#' Latent state 4: CD4 > 500, HIV acquired more than 6 months ago
#' Latent state 5: 500 < CD4 < 350, HIV acquired more than 6 months ago
#' Latent state 6: 350 < CD4 < 200, HIV acquired more than 6 months ago
#' Latent state 7: CD4 < 200, HIV acquired more than 6 months ago
#'
#' @param t Time index.
#' @param d_probs Diagnosis probabilities, a 3D array with dimensions \code{[nquar, nage, nstates]}.
#' @param q_age Age-specific progression probabilities.
#' @param nlatst Number of latent states.
#' @param nage Number of age groups.
#' @param rita Logical indicating whether RITA-specific logic should be applied.
#' @return Q matrix.
#' @examples
#' \dontrun{
#' q_age <- matrix(rep(c(0.2, 0.3, 0.4, 0.5), 52), nrow = 52, byrow = TRUE)
#' d_probs <- array(0.1, dim = c(10, 52, 4))
#' Q_mat(t = 10, d_probs, q_age, nlatst = 4, nage = 52, rita = FALSE)
#' d_probs <- array(0.1, dim = c(10, 52, 5))
#' Q_mat(t = 10, d_probs, q_age, nlatst = 7, nage = 52, rita = TRUE)
#' }
#' @keywords internal
Q_mat <- function(t, d_probs, q_age, nlatst, nage = 52, rita = FALSE) {
  out <- array(0, dim = c(nage, nlatst, nlatst))

  for (a in 1:nage) {
    if (rita) {
      # RITA-specific logic
      diag(out[a, , ])[4:7] <- (1 - d_probs[t, a, 2:5]) * (1 - q_age[a, ])
      out[a, 1, 2] <- (1 - d_probs[t, a, 1]) * (1 - q_age[a, 1])
      out[a, 1, 3] <- (1 - d_probs[t, a, 1]) * q_age[a, 1]
      out[a, 2, 4] <- (1 - d_probs[t, a, 1]) * (1 - q_age[a, 1])
      out[a, 2, 5] <- (1 - d_probs[t, a, 1]) * q_age[a, 1]

      out[a, 3, 5] <- (1 - d_probs[t, a, 1]) * (1 - q_age[a, 2])
      out[a, 3, 6] <- (1 - d_probs[t, a, 1]) * q_age[a, 2] # allow progression to lower CD4 state
      for (i in 4:6) {
        out[a, i, i + 1] <- (1 - d_probs[t, a, i - 2]) * q_age[a, i - 3]
      }
    } else {
      # Standard logic
      diag(out[a, , ])[1:4] <- (1 - d_probs[t, a, ]) * (1 - q_age[a, ])
      for (i in 1:3) {
        out[a, i, i + 1] <- (1 - d_probs[t, a, i]) * q_age[a, i]
      }
    }
  }

  return(out)
}

#' Generate D Matrix
#'
#' Generates the D matrix for the epidemic model.
#' CD4-only diagnosed states definition:
#' Diagnosed state 1: CD4 > 500 at diagnosis
#' Diagnosed state 2: 500 < CD4 < 350 at diagnosis
#' Diagnosed state 3: 350 < CD4 < 200 at diagnosis
#' Diagnosed state 4: CD4 < 200 at diagnosis
#'
#' RITA diagnosed states definition:
#' Diagnosed state 1: RITA positive diagnosis
#' Diagnosed state 2: CD4 > 500 at diagnosis
#' Diagnosed state 3: 500 < CD4 < 350 at diagnosis
#' Diagnosed state 4: 350 < CD4 < 200 at diagnosis
#' Diagnosed state 5: CD4 < 200 at diagnosis
#'
#' @param t Time index.
#' @param d_probs Diagnosis probabilities, a 3D array with dimensions \code{[nquar, nage, nstates]}.
#' @param q_age Age-specific progression probabilities.
#' @param nlatst Number of latent states.
#' @param ndiagst Number of diagnosed states.
#' @param nage Number of age groups.
#' @param rita Logical indicating whether RITA-specific logic should be applied.
#' @return D matrix.
#' @examples
#' \dontrun{
#' q_age <- matrix(rep(c(0.2, 0.3, 0.4, 0.5), 52), nrow = 52, byrow = TRUE)
#' d_probs <- array(0.1, dim = c(10, 52, 4))
#' D_mat(t = 10, d_probs = d_probs, q_age, nlatst = 4, ndiagst = 5, nage = 52, rita = FALSE)
#' d_probs <- array(0.1, dim = c(10, 52, 5))
#' D_mat(t = 10, d_probs = d_probs, q_age, nlatst = 7, ndiagst = 6, nage = 52, rita = TRUE)
#' }
#' @keywords internal
D_mat <- function(t, d_probs, q_age, nlatst, ndiagst, nage = 52, rita = FALSE) {
  out <- array(0, dim = c(nage, nlatst, ndiagst))

  for (a in 1:nage) {
    if (rita) {
      # RITA-specific logic
      out[a, 1:3, 1] <- d_probs[t, a, 1]
      out[a, 4:7, 2:5] <- diag(d_probs[t, a, 2:5])
      out[a, 7, 6] <- (1 - d_probs[t, a, 5]) * q_age[a, 4]
    } else {
      # Standard logic
      diag(out[a, , ])[1:4] <- d_probs[t, a, ]
      out[a, 4, 5] <- (1 - d_probs[t, a, 4]) * q_age[a, 4]
    }
  }

  return(out)
}

#' Simulate epidemic over time
#'
#' Simulates the epidemic over a specified number of time quarters using the provided Q and D matrices, initial latent state, and new infections.
#' @param Q_l List of Q matrices for each time quarter.
#' @param D_l List of D matrices for each time quarter.
#' @param h Matrix of new infections over time.
#' @param lat_init Initial latent state vector.
#' @param config Configuration list containing simulation parameters.
#' @return A list containing matrices of latent and diagnosed states over time, as well as detailed quarterly records.
#' @examples
#' \dontrun{
#' config <- define_sim_config(nquar = 10)
#' q_age <- matrix(rep(c(0.0976, 0.115, 0.116, 0.148), 52), nrow = 52, byrow = TRUE)
#' d_probs <- matrix(rep(0.1, 40), nrow = 10)
#' Q_l <- lapply(1:10, function(x) {
#'    Q_mat(x, d_probs, q_age, 4)
#' })
#' D_l <- lapply(1:10, function(x) {
#'    D_mat(x, d_probs, q_age, 4, 5)
#' })
#' lat_init <- array(0, dim = c(52, 4))
#' lat_init[, 1] <- c(40, 30, 20, 10)
#' h <- array(matrix(rep(10, 20), nrow = 10), dim = c(52, 10, 4))
#' epi_mat(Q_l, D_l, h, lat_init, config)
#' }
#' @keywords internal
epi_mat <- function(Q_l, D_l, h, lat_init, config) {
  ### Epidemics loop

  # Matrices for age-stratified data
  lat_age <- array(0, dim = c(config$nage, config$nquar, config$nlatst))
  diags_age <- array(0, dim = c(config$nage, config$nquar, config$ndiagst))

  for (a in 1:config$nage) {
    lat_age[a, 1, ] <- t(Q_l[[1]][a, , ]) %*%
      lat_init[a, ] +
      h[a, 1, ]
    diags_age[a, 1, ] <- t(D_l[[1]][a, , ]) %*% lat_init[a, ]
  }

  # Simulate epidemic with aging
  for (t in 2:config$nquar) {
    # Apply epidemic dynamics for each age group
    for (a in 1:config$nage) {
      # Epidemic progression
      lat_age[a, t, ] <- t(Q_l[[t]][a, , ]) %*% lat_age[a, t - 1, ] + h[a, t, ]
      diags_age[a, t, ] <- t(D_l[[t]][a, , ]) %*% lat_age[a, t - 1, ]
    }

    # Age people every 4 quarters (1 year)
    if ((t - 1) %% 4 == 0 && t < config$nquar) {
      # Shift latent states up one age group, retain oldest group
      lat_age[config$nage - 1, t, ] <- lat_age[config$nage, t, ] +
        lat_age[config$nage - 1, t, ]
      lat_age[2:config$nage, t, ] <- lat_age[1:(config$nage - 1), t, ]
      lat_age[1, t, ] <- 0 # No one ages into youngest group
    }
  }

  # lat <- matrix(0, nrow = nquar, ncol = nlatst)
  # diags <- matrix(0, nrow = nquar, ncol = ndiagst)

  # # First iteration using matrix multiplication
  # lat[1, ] <- t(Q_l[[1]]) %*% lat_init + h[1, ]
  # diags[1, ] <- t(D_l[[1]]) %*% lat_init

  # # Main loop with vectorized operations
  # for (t in 2:nquar) {
  #   lat[t, ] <- (t(Q_l[[t]]) %*% lat[t - 1, ]) + h[t, ]
  #   diags[t, ] <- t(D_l[[t]]) %*% lat[t - 1, ]
  # }

  # # method to record when infections enter the latent state
  # lat_quar <- list()
  # diags_quar <- list()

  # lat_quar[[1]] <- matrix(0, nrow = nquar, ncol = nlatst)
  # diags_quar[[1]] <- matrix(0, nrow = nquar, ncol = ndiagst)
  # lat_quar[[1]][1, ] <- (t(Q_l[[1]]) %*% lat_init) + h[1, ]
  # diags_quar[[1]][1, ] <- t(D_l[[1]]) %*% lat_init

  # for (u in 1:nquar) {
  #   if (u > 1) {
  #     lat_quar[[u]] <- matrix(0, nrow = nquar, ncol = nlatst)
  #     diags_quar[[u]] <- matrix(0, nrow = nquar, ncol = ndiagst)
  #     lat_quar[[u]][u, ] <- h[u, ]
  #   }

  #   if (u < nquar) {
  #     for (t in (u + 1):nquar) {
  #       lat_quar[[u]][t, ] <- t(Q_l[[t]]) %*% lat_quar[[u]][t - 1, ]
  #       diags_quar[[u]][t, ] <- t(D_l[[t]]) %*% lat_quar[[u]][t - 1, ]
  #     }
  #   }
  # }

  out <- list(
    lat = apply(lat_age, c(2, 3), sum), # Sum over ages
    diags = apply(diags_age, c(2, 3), sum),
    lat_age = lat_age,
    diags_age = diags_age
  )
  return(out)
}

#' Simulate HIV Epidemic Data
#'
#' Simulates HIV epidemic data based on the provided configuration.
#'
#' @param h_pattern Pattern type for infections ("constant", "increasing", "varying", "varying_2").
#' @param o_pattern Pattern type for omega ("constant", "increasing", "step_change", "step_change_late", "none").
#' @param d_pattern Pattern type for diagnosis probabilities ("constant", "varying").
#' @param p_pattern Pattern type for proportion UK ("constant", "increasing").
#' @param m_pattern Pattern type for mu ("constant", "differing", "varying").
#' @param r_pattern Pattern type for RITA data ("constant", "increasing", "none").
#' @param a_pattern Pattern type for age-based data ("constant", "wedge", "varying_wedge").
#' @param sim_type If provided, a preset type of simulation to define patterns for.
#' @param sim_seed Random seed for reproducibility. If NULL, a random seed is generated.
#' @param dist Logical indicating whether to add stochastic variation.
#' @param migration Logical indicating whether to include migration effects.
#' @param rita Logical indicating whether to simulate RITA (Recent Infection Testing Algorithm) data.
#' @param age Logical indicating whether to simulate age-stratified data.
#' @param nquar Number of time quarters.
#' @param nage Number of age groups.
#' @param q Vector of progression probabilities.
#' @param target_inf Target total number of HIV acquisitions.
#' @param target_mig Target total number of migrations.
#' @param init_props Initial state proportions.
#' @param init_props_abroad Initial state proportions for infections acquired abroad.
#' @param abroad_prop Proportion of infections acquired abroad at the outset.
#' @param cd4_completeness CD4 completeness.
#' @param rita_completeness RITA completeness.
#' @param age_diag_effect Numeric scalar controlling the strength of age variation in diagnosis probabilities. Applied as a logit-additive offset scaled linearly across ages.
#' @param mig_diag_offset Numeric vector of per-state logit-scale offsets applied to diagnosis probabilities for the migrant group (y=2). Length must equal the number of diagnosis states (4 for non-RITA, 5 for RITA). NULL (default) applies no offset.
#' @param frr False recency rate for RITA misclassification.
#' @param eta_sd Standard deviation of the approximate melding prior for the log ratio of abroad-acquired to UK-acquired infections. Passed directly to Stan as data. Smaller values anchor the model more tightly to the external \code{ratio_abroad_uk} evidence. Default is 0.10, which matches the simulation noise on the log scale.
#' @param detect_obs_sd Standard deviation of detection probability.
#' @param detect_obs_every Positive integer giving the spacing, in quarters, between external detection-probability observations. Use values larger than 1 to weaken the auxiliary likelihood by observing only occasional linkage studies.
#' @return List containing simulated epidemic data.
#' @examples
#' simulate_diagnoses()
#' @export
simulate_diagnoses <- function(
  h_pattern = "constant",
  o_pattern = "constant",
  d_pattern = "constant",
  p_pattern = "constant",
  m_pattern = "constant",
  r_pattern = "constant",
  a_pattern = "constant",
  sim_type = NULL,
  sim_seed = NULL,
  dist = TRUE,
  migration = FALSE,
  rita = FALSE,
  age = FALSE,
  nquar = 120, # number of quarters
  nage = 52, # number of age groups
  q = c(0.0976, 0.115, 0.116, 0.148), # progression probabilities
  target_inf = 90000, # target total number of HIV acquisitions
  target_mig = 12000, # target total number of migrations
  init_props = c(7, 3.5, 1.75, 0.5), # initial state proportions
  init_props_abroad = c(1.75, 2.25, 2.25, 1.25), # initial state proportions abroad
  abroad_prop = 0.2, # proportion acquired abroad at outset
  cd4_completeness = 1, # CD4 completeness
  rita_completeness = 1, # RITA completeness
  age_diag_effect = 0, # age effect on diagnosis probabilities
  mig_diag_offset = NULL, # per-state logit offsets for migrant diagnosis probs
  frr = 0.019, # false recency rate for RITA misclassification
  eta_sd = 0.10, # sd of the melding prior for log(ratio_abroad_uk)
  detect_obs_sd = 0,
  detect_obs_every = 1
) {
  sim_config <- define_sim_config(
    h_pattern,
    o_pattern,
    d_pattern,
    p_pattern,
    m_pattern,
    r_pattern,
    a_pattern,
    sim_type,
    sim_seed,
    dist,
    migration,
    rita,
    age,
    nquar,
    nage,
    q,
    target_inf,
    target_mig,
    init_props,
    init_props_abroad,
    abroad_prop,
    cd4_completeness,
    rita_completeness,
    age_diag_effect,
    mig_diag_offset,
    frr,
    detect_obs_every
  )
  set.seed(sim_config$sim_seed)

  # Generate model components
  infs <- generate_infections(
    sim_config$nquar,
    sim_config$nlatst,
    sim_config$h_pattern,
    sim_config$target_inf,
    sim_config$dist
  )
  omega <- generate_omega(
    sim_config$nquar,
    sim_config$o_pattern,
    sim_config$target_mig,
    sim_config$dist
  )

  # ratio of born-abroad UK-acquired to UK-born UK-acquired; can exceed 1
  ratio_abroad_uk <- generate_ratio_abroad_uk(
    sim_config$nquar,
    sim_config$p_pattern,
    sim_config$dist
  )

  # mu is a matrix of dimension [nquar, nlatst + 1] giving migration probabilities
  # mu are the true migration probabilities
  mu <- generate_mu(
    sim_config$nquar,
    4,
    sim_config$m_pattern,
    sim_config$dist
  )

  # detection probabilities for those migrating already diagnosed
  # In Stan: detect_prob[t] = inv_logit(alpha_detect + 0.5 * X_mp * beta_detect)
  if (sim_config$dist) {
    detect_prob <- invlogit(stats::rnorm(
      sim_config$nquar,
      logit(seq(0.7, 0.8, length.out = sim_config$nquar)),
      0.1
    ))
  } else {
    detect_prob <- seq(0.7, 0.8, length.out = sim_config$nquar)
  }

  if (sim_config$rita) {
    # pad first three columns with zeros for RITA states
    # because migrants assumed to enter non-RITA states only
    mu <- cbind(
      matrix(0, nrow = sim_config$nquar, ncol = 3),
      mu
    )
  }

  # m is a matrix of dimension [nquar, nlatst + 1] with the last column being those
  # who migrate already diagnosed (total, before detection split).
  # Use true mu (not mu_detect) to match the Stan model, which handles
  # misclassification separately via detect_prob rather than redistributing
  # undetected abroad diagnoses into CD4 states.
  m <- generate_migration(
    sim_config$nquar,
    sim_config$nlatst,
    omega,
    mu,
    sim_config$dist
  )
  d_probs <- generate_diagnosis_probs(
    sim_config$nquar,
    sim_config$nage,
    sim_config$d_pattern,
    sim_config$dist,
    sim_config$rita,
    sim_config$age_diag_effect
  )

  # Build migrant diagnosis probabilities by applying per-state logit offsets.
  # Mirrors alpha_diag_mig in the Stan migration models: migrants (y=2) may
  # have a systematically different baseline diagnosis rate than natives (y=1)
  # while sharing the same time trend.
  mig_offset <- sim_config$mig_diag_offset
  nstates_diag <- dim(d_probs)[3]
  if (is.null(mig_offset) || all(mig_offset == 0)) {
    d_probs_mig <- d_probs
  } else {
    if (length(mig_offset) != nstates_diag) {
      stop(paste0(
        "mig_diag_offset must have length ",
        nstates_diag,
        " (one per diagnosis state); got length ",
        length(mig_offset),
        "."
      ))
    }
    d_probs_mig <- d_probs
    for (k in seq_len(nstates_diag)) {
      d_probs_mig[,, k] <- invlogit(logit(d_probs[,, k]) + mig_offset[k])
    }
  }

  # Generate initial prevalence
  lat_init_uk <- round(
    sum(infs[1, ]) / (1 + sim_config$abroad_prop) * sim_config$init_props
  )
  lat_init_uk_abroad <- round(
    sum(infs[1, ]) *
      (sim_config$abroad_prop / (1 + sim_config$abroad_prop)) *
      sim_config$init_props
  )
  lat_init_abroad <- round(
    sum(m[1, 1:(sim_config$nlatst)]) * sim_config$init_props_abroad
  )

  # Generate Q and D matrices for native group (y=1)
  Q_l <- lapply(1:sim_config$nquar, function(x) {
    Q_mat(
      x,
      d_probs,
      sim_config$q_age,
      sim_config$nlatst,
      sim_config$nage,
      sim_config$rita
    )
  })
  D_l <- lapply(1:sim_config$nquar, function(x) {
    D_mat(
      x,
      d_probs,
      sim_config$q_age,
      sim_config$nlatst,
      sim_config$ndiagst,
      sim_config$nage,
      sim_config$rita
    )
  })

  # Generate Q and D matrices for migrant group (y=2)
  Q_l_mig <- lapply(1:sim_config$nquar, function(x) {
    Q_mat(
      x,
      d_probs_mig,
      sim_config$q_age,
      sim_config$nlatst,
      sim_config$nage,
      sim_config$rita
    )
  })
  D_l_mig <- lapply(1:sim_config$nquar, function(x) {
    D_mat(
      x,
      d_probs_mig,
      sim_config$q_age,
      sim_config$nlatst,
      sim_config$ndiagst,
      sim_config$nage,
      sim_config$rita
    )
  })

  # Distribute across ages using UK age distribution
  # Dimensions: [age, quarter, state]
  agedist_inf <- generate_age_dist(
    sim_config$nquar,
    sim_config$nage,
    sim_config$a_pattern,
    10,
    30
  )

  # Split total uk-acq across into uk-born and born abroad
  total_uk_acq <- infs[, 1]
  infs_age <- array(
    0,
    dim = c(sim_config$nage, sim_config$nquar, sim_config$nlatst)
  )
  infs_age[,, 1] <- t(total_uk_acq * agedist_inf)

  born_uk_acq_uk <- total_uk_acq / (1 + ratio_abroad_uk)
  born_abroad_acq_uk <- total_uk_acq * ratio_abroad_uk / (1 + ratio_abroad_uk)

  # Distribute UK incidence across ages
  infs_uk_age <- array(
    0,
    dim = c(sim_config$nage, sim_config$nquar, sim_config$nlatst)
  )
  infs_uk_age[,, 1] <- t(born_uk_acq_uk * agedist_inf)

  infs_uk_abroad_age <- array(
    0,
    dim = c(sim_config$nage, sim_config$nquar, sim_config$nlatst)
  )
  infs_uk_abroad_age[,, 1] <- t(born_abroad_acq_uk * agedist_inf)

  # Distribute UK-born initial prevalence across ages
  lat_init_uk_age <- round_preserve(outer(agedist_inf[1, ], lat_init_uk))
  lat_init_uk_abroad_age <- round_preserve(outer(
    agedist_inf[1, ],
    lat_init_uk_abroad
  ))

  # Distribute abroad initial prevalence across ages
  agedist_abroad <- generate_age_dist(
    sim_config$nquar,
    sim_config$nage,
    sim_config$a_pattern,
    20,
    40
  )
  m_age <- array(
    0,
    dim = c(sim_config$nage, sim_config$nquar, (sim_config$nlatst + 1))
  )
  for (i in 1:(sim_config$nlatst + 1)) {
    m_age[,, i] <- t(m[, i] * agedist_abroad)
  }
  lat_init_abroad_age <- round_preserve(outer(
    agedist_abroad[1, ],
    lat_init_abroad
  ))

  # Useful to store omega by age for diagnostics
  omega_age <- omega * agedist_abroad

  # Split compartment 5 (total diagnosed migrants) into detected (ABROAD) and
  # misclassified, matching the Stan model's misclassification decomposition:
  #   ABROAD[t]          = omega[t] * m[5,t] * detect_prob[t]
  #   misclassified[t,a] = omega[t] * m[5,t] * (1 - detect_prob[t]) * age_prop[t,a]
  total_diag_age <- round_matrix_preserve(
    m_age[,, sim_config$nlatst + 1],
    colwise = TRUE
  ) # [nage, nquar], integer, column sums preserved

  if (sim_config$dist) {
    ABROAD_age_mat <- matrix(
      stats::rbinom(
        sim_config$nage * sim_config$nquar,
        as.integer(total_diag_age),
        rep(detect_prob, each = sim_config$nage)
      ),
      nrow = sim_config$nage,
      ncol = sim_config$nquar
    )
  } else {
    ABROAD_age_mat <- round_matrix_preserve(
      sweep(total_diag_age, 2, detect_prob, "*"),
      colwise = TRUE
    )
  }
  misclassified_age <- total_diag_age - ABROAD_age_mat # [nage, nquar]
  misclassified <- colSums(misclassified_age)

  # Run epidemic simulations
  epi_uk <- epi_mat(
    Q_l,
    D_l,
    infs_uk_age,
    lat_init_uk_age,
    sim_config
  )
  epi_uk_abroad <- epi_mat(
    Q_l_mig,
    D_l_mig,
    infs_uk_abroad_age,
    lat_init_uk_abroad_age,
    sim_config
  )
  epi_abroad <- epi_mat(
    Q_l_mig,
    D_l_mig,
    m_age[,, 1:sim_config$nlatst],
    lat_init_abroad_age,
    sim_config
  )

  sim_data <- list(
    "epi_uk" = epi_uk,
    "epi_uk_abroad" = epi_uk_abroad,
    "epi_abroad" = epi_abroad,
    "omega" = omega,
    "omega_age" = omega_age,
    "infs" = infs,
    "infs_age" = infs_age,
    "m" = m,
    "m_age" = m_age,
    "mu" = mu,
    "d_probs" = d_probs,
    "d_probs_mig" = d_probs_mig,
    "ratio_abroad_uk" = ratio_abroad_uk,
    "lat_init_uk" = lat_init_uk,
    "lat_init_uk_abroad" = lat_init_uk_abroad,
    "lat_init_abroad" = lat_init_abroad,
    "detect_prob" = detect_prob,
    "misclassified" = misclassified,
    "misclassified_age" = misclassified_age
  )

  # Perform validation checks before rounding
  validate_sim_data(sim_data)

  # Apply rounding which preserves sums
  sim_data$epi_uk <- round_fields(sim_data$epi_uk)
  sim_data$epi_uk_abroad <- round_fields(sim_data$epi_uk_abroad)
  sim_data$epi_abroad <- round_fields(sim_data$epi_abroad)

  # Prepare Stan data
  epi_uk <- sim_data$epi_uk
  epi_uk_abroad <- sim_data$epi_uk_abroad
  epi_abroad <- sim_data$epi_abroad
  born_abroad <- epi_abroad$diags + epi_uk_abroad$diags
  born_abroad_age <- epi_abroad$diags_age + epi_uk_abroad$diags_age

  # cd4_dx is an nquar x 2 x (ndiagst - 1) array
  cd4_dx <- array(0, dim = c(sim_config$nquar, 2, sim_config$ndiagst - 1))
  for (i in 1:(sim_config$ndiagst - 1)) {
    cd4_dx[, 1, i] <- epi_uk$diags[, i]
    cd4_dx[, 2, i] <- born_abroad[, i]
  }
  # cd4_dx_age is an [nage, nquar, 2, ndiagst - 1] array
  cd4_dx_age <- array(
    0,
    dim = c(sim_config$nage, sim_config$nquar, 2, sim_config$ndiagst - 1)
  )
  for (i in 1:(sim_config$ndiagst - 1)) {
    cd4_dx_age[,, 1, i] <- epi_uk$diags_age[,, i]
    cd4_dx_age[,, 2, i] <- born_abroad_age[,, i]
  }

  # Distribute misclassified abroad diagnoses across CD4 states for inclusion
  # in CD4 data. In reality misclassified cases have CD4 counts recorded,
  # we just don't know they were already diagnosed abroad. Use mu[1:nlatst]
  # proportions as proxy for their CD4 distribution at arrival.
  if (sim_config$rita) {
    mu_cd4_cols <- 4:sim_config$nlatst
    dx_cd4_cols <- 2:(sim_config$ndiagst - 1)
  } else {
    mu_cd4_cols <- 1:sim_config$nlatst
    dx_cd4_cols <- 1:(sim_config$ndiagst - 1)
  }

  n_cd4_cats <- length(dx_cd4_cols)
  misclassified_cd4_age <- array(
    0,
    dim = c(sim_config$nage, sim_config$nquar, sim_config$ndiagst - 1)
  )

  for (t in 1:sim_config$nquar) {
    mu_raw <- mu[t, mu_cd4_cols]
    mu_total <- sum(mu_raw)
    mu_cd4 <- if (mu_total > 0) {
      mu_raw / mu_total
    } else {
      rep(1 / n_cd4_cats, n_cd4_cats)
    }

    for (a in 1:sim_config$nage) {
      n <- misclassified_age[a, t]
      if (n > 0) {
        if (sim_config$dist) {
          misclassified_cd4_age[a, t, dx_cd4_cols] <-
            as.vector(stats::rmultinom(1, n, mu_cd4))
        } else {
          misclassified_cd4_age[a, t, dx_cd4_cols] <-
            round_preserve(n * mu_cd4)
        }
      }
    }
  }

  # Add misclassified CD4 counts to migrant CD4 data
  cd4_dx_age[,, 2, ] <- cd4_dx_age[,, 2, ] + misclassified_cd4_age
  for (i in 1:(sim_config$ndiagst - 1)) {
    cd4_dx[, 2, i] <- cd4_dx[, 2, i] + colSums(misclassified_cd4_age[,, i])
  }

  # aids_dx is an nquar x 2 matrix
  aids_dx <- matrix(
    c(epi_uk$diags[, sim_config$ndiagst], born_abroad[, sim_config$ndiagst]),
    nrow = sim_config$nquar
  )
  # aids_dx_age is an [nage, nquar, 2] array
  aids_dx_age <- array(0, dim = c(sim_config$nage, sim_config$nquar, 2))
  aids_dx_age[,, 1] <- epi_uk$diags_age[,, sim_config$ndiagst]
  aids_dx_age[,, 2] <- born_abroad_age[,, sim_config$ndiagst]

  # hix_dx is an nquar x 2 matrix
  hiv_dx <- matrix(
    c(epi_uk$diags[, 1], born_abroad[, 1]),
    nrow = sim_config$nquar
  ) +
    matrix(c(epi_uk$diags[, 2], born_abroad[, 2]), nrow = sim_config$nquar) +
    matrix(c(epi_uk$diags[, 3], born_abroad[, 3]), nrow = sim_config$nquar) +
    matrix(c(epi_uk$diags[, 4], born_abroad[, 4]), nrow = sim_config$nquar)
  # hiv_dx_age is an [nage, nquar, 2] array
  hiv_dx_age <- array(0, dim = c(sim_config$nage, sim_config$nquar, 2))
  hiv_dx_age[,, 1] <-
    epi_uk$diags_age[,, 1] +
    epi_uk$diags_age[,, 2] +
    epi_uk$diags_age[,, 3] +
    epi_uk$diags_age[,, 4]
  hiv_dx_age[,, 2] <-
    born_abroad_age[,, 1] +
    born_abroad_age[,, 2] +
    born_abroad_age[,, 3] +
    born_abroad_age[,, 4]

  if (sim_config$rita) {
    hiv_dx <- hiv_dx +
      matrix(c(epi_uk$diags[, 5], born_abroad[, 5]), nrow = sim_config$nquar)
    hiv_dx_age[,, 1] <- hiv_dx_age[,, 1] + epi_uk$diags_age[,, 5]
    hiv_dx_age[,, 2] <- hiv_dx_age[,, 2] + born_abroad_age[,, 5]
  }

  # Add misclassified abroad diagnoses directly to migrant HIV totals.
  # In the Stan model these bypass CD4 progression entirely and appear as

  # immediate HIV diagnoses (not through diagnosis probabilities).
  hiv_dx[, 2] <- hiv_dx[, 2] + misclassified
  hiv_dx_age[,, 2] <- hiv_dx_age[,, 2] + misclassified_age

  if (sim_config$migration) {
    init_prev <- lat_init_uk
    init_prev_age <- lat_init_uk_age
    init_prev_abroad <- lat_init_uk_abroad + lat_init_abroad
    init_prev_abroad_age <- lat_init_uk_abroad_age + lat_init_abroad_age
    hiv_dx_age <- aperm(hiv_dx_age, c(2, 1, 3))
    aids_dx_age <- aperm(aids_dx_age, c(2, 1, 3))
    cd4_dx_age <- aperm(cd4_dx_age, c(2, 1, 3, 4))
  } else {
    # sum over migration status
    hiv_dx <- rowSums(hiv_dx)
    aids_dx <- rowSums(aids_dx)
    cd4_dx <- apply(
      cd4_dx[,, 1:(sim_config$ndiagst - 1), drop = FALSE],
      c(1, 3),
      sum
    )
    hiv_dx_age <- apply(hiv_dx_age[,, 1:2, drop = FALSE], c(1, 2), sum)
    hiv_dx_age <- t(hiv_dx_age)
    aids_dx_age <- apply(aids_dx_age[,, 1:2, drop = FALSE], c(1, 2), sum)
    aids_dx_age <- t(aids_dx_age)
    cd4_dx_age <- apply(
      cd4_dx_age[,, 1:2, 1:(sim_config$ndiagst - 1), drop = FALSE],
      c(1, 2, 4),
      sum
    )
    cd4_dx_age <- aperm(cd4_dx_age, c(2, 1, 3))
    init_prev <- lat_init_uk + lat_init_uk_abroad
    init_prev_age <- lat_init_uk_age + lat_init_uk_abroad_age
    init_prev_abroad <- lat_init_abroad
    init_prev_abroad_age <- lat_init_abroad_age
  }

  # Apply false recency rate misclassification to RITA-classified data
  if (sim_config$rita && sim_config$frr > 0) {
    cd4_dx <- apply_frr(cd4_dx, sim_config$frr, sim_config$dist)
    cd4_dx_age <- apply_frr(cd4_dx_age, sim_config$frr, sim_config$dist)
  }

  # Generate noisy logit-scale detect_prob observations (data-linkage auxiliary data).
  # detect_obs_sd = 0 disables the auxiliary likelihood in Stan.
  # Larger detect_obs_every values mimic occasional linkage studies rather than
  # an observation every quarter.
  detect_obs_index <- if (detect_obs_sd > 0) {
    seq.int(1, sim_config$nquar, by = sim_config$detect_obs_every)
  } else {
    integer(0)
  }
  detect_logit_obs <- if (length(detect_obs_index) > 0) {
    stats::rnorm(
      length(detect_obs_index),
      logit(detect_prob[detect_obs_index]),
      detect_obs_sd
    )
  } else {
    numeric(0)
  }

  # ABROAD data: detected portion of diagnosed migrants
  ABROAD <- colSums(ABROAD_age_mat)
  ABROAD_age <- t(ABROAD_age_mat)

  # Make stan data list
  stan_data <- list(
    "HIV" = hiv_dx,
    "HIV_age" = hiv_dx_age,
    "AIDS" = aids_dx,
    "AIDS_age" = aids_dx_age,
    "CD4" = cd4_dx,
    "CD4_age" = cd4_dx_age,
    "ABROAD" = ABROAD,
    "ABROAD_age" = ABROAD_age,
    "nquar" = sim_config$nquar,
    "nage" = sim_config$nage,
    "q" = sim_config$q,
    "q_age" = t(sim_config$q_age),
    "init_prev" = init_prev,
    "init_prev_age" = t(init_prev_age),
    "init_prev_abroad" = init_prev_abroad,
    "init_prev_abroad_age" = t(init_prev_abroad_age),
    # log(p) is the appropriate transformation for the eta parameter in the migration model
    "eta_mean" = log(sim_data$ratio_abroad_uk),
    "eta_sd" = eta_sd,
    "detect_obs_index" = detect_obs_index,
    "detect_logit_obs" = detect_logit_obs,
    "detect_obs_sd" = detect_obs_sd,
    "sim_config" = sim_config
  )

  # Perform validation checks
  validate_stan_data(stan_data, sim_data, sim_config)

  # Adjust CD4 counts if completeness less than 100% using a normal approximation to the binomial distribution
  if (sim_config$cd4_completeness < 1) {
    stan_data$CD4 <- round(
      stan_data$CD4 *
        pmin(
          1,
          stats::rnorm(length(stan_data$CD4), sim_config$cd4_completeness, 0.05)
        )
    )
    stan_data$CD4_age <- round(
      stan_data$CD4_age *
        pmin(
          1,
          stats::rnorm(
            length(stan_data$CD4_age),
            sim_config$cd4_completeness,
            0.05
          )
        )
    )
  }
  # Adjust RITA counts if completeness less than 100% using a normal approximation to the binomial distribution
  if (sim_config$rita && sim_config$rita_completeness < 1) {
    stan_data$CD4 <- round(
      stan_data$CD4 *
        pmin(
          1,
          stats::rnorm(
            length(stan_data$CD4),
            sim_config$rita_completeness,
            0.05
          )
        )
    )
    stan_data$CD4_age <- round(
      stan_data$CD4_age *
        pmin(
          1,
          stats::rnorm(
            length(stan_data$CD4_age),
            sim_config$rita_completeness,
            0.05
          )
        )
    )
  }

  return(list(
    "stan_data" = stan_data,
    "sim_data" = sim_data
  ))
}
