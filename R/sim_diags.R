#' Define simulation configuration
#'
#' Defines the patterns for various parameters in the epidemic simulation.
#' The number of latent and diagnosis states will be inferred from the simulation specification.
#'
#' @param h_pattern Pattern type for infections ("constant", "increasing", "varying", "varying_2").
#' @param o_pattern Pattern type for omega ("constant", "increasing", "step_change_1", "step_change_2", "none").
#' @param d_pattern Pattern type for diagnosis probabilities ("constant", "varying").
#' @param p_pattern Pattern for the born-abroad:UK-born ratio of UK-acquired
#'   incidence ("constant", "increasing").
#' @param m_pattern Pattern type for mu ("constant", "differing", or
#'   "varying").
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
#' @param abroad_prop Proportion of infections acquired abroad at the outset.
#' @param cd4_completeness CD4 completeness.
#' @param rita_completeness RITA completeness.
#' @param rita_completeness_by_state Optional numeric vector of length 5
#'   giving state-specific completeness for RITA simulations. Elements
#'   correspond to RITA positive, CD4 500 or above, CD4 350--499,
#'   CD4 200--349, and CD4 below 200.
#'   Use lower values in lower CD4 states to simulate informative missingness.
#' @param age_diag_effect Numeric scalar controlling the strength of age variation in diagnosis probabilities. Applied as a logit-additive offset scaled linearly across ages.
#' @param mig_diag_offset Numeric vector of per-state logit-scale offsets applied
#'   to diagnosis probabilities for people born abroad (`y = 2`). Length must
#'   equal the number of diagnosis states (4 for non-RITA, 5 for RITA). `NULL`
#'   (default) applies no offset.
#' @param frr False recency rate for RITA misclassification.
#' @param detect_obs_times Integer vector of quarter indices (1-based) at which external detection-probability observations are generated.
#' @param gamma_obs_grouped Logical scalar or logical vector of length
#'   \code{nquar} controlling which quarters receive age-stratified external
#'   migration observations on the same 4 age groups as \code{ABROAD_group}.
#'   The scalar \code{gamma_obs} series is always simulated; quarters where
#'   \code{gamma_obs_grouped} is \code{FALSE} have \code{-1} rows in
#'   \code{gamma_obs_group}.
#' @param abroad_age_trend Total logit-scale change (start to end) in the
#'   15-24/25-34 age-group shares of born-abroad arrivals relative to the 45+
#'   reference group; 0 = time-constant age composition.
#' @param detect_step_quarter Quarter index (1-based) of an ascertainment
#'   step in detect_prob (e.g. a surveillance reporting change). \code{NULL}
#'   (default) means no step.
#' @param detect_step_size Logit-scale jump in detect_prob applied from
#'   \code{detect_step_quarter} onward.
#' @param gamma_obs_scale_trend Total log-scale change (start to end) in the
#'   HIV-positive share of the external migration index (gamma_obs_scale);
#'   0 = constant share.
#' @param gamma_obs_drift_sd Quarterly random-walk innovation SD on the log
#'   ratio of the external migration index to HIV-positive migration. This
#'   controls stochastic variation in migrant HIV prevalence and is passed to
#'   the matching Stan model. Set to 0 to disable.
#' @param recent_arrival_completeness Proportion of eligible diagnoses with
#'   avidity testing done and a recorded year of arrival. Combined
#'   multiplicatively with \code{cd4_completeness} (a recorded CD4 value is
#'   also required) to give the overall completeness of the recent-arrival
#'   CD4 auxiliary data. Set to 0 to disable.
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
  nquar = 40, # number of quarters
  nage = 52, # number of age groups
  q = c(0.0976, 0.115, 0.116, 0.148), # progression probabilities
  target_inf = 30000, # target total number of HIV acquisitions
  target_mig = 6000, # target total number of migrations
  abroad_prop = 0.2, # proportion acquired abroad at outset
  cd4_completeness = 1, # CD4 completeness
  rita_completeness = 1, # RITA completeness
  rita_completeness_by_state = NULL, # optional state-specific RITA completeness
  age_diag_effect = 0, # age effect on diagnosis probabilities
  mig_diag_offset = NULL, # per-state logit offsets for migrant diagnosis probs
  frr = 0.019, # false recency rate for RITA misclassification
  detect_obs_times = NULL, # quarter indices for detect_prob observations
  gamma_obs_grouped = FALSE, # whether external migration data are age-stratified
  gamma_obs_sd = 0, # sd of the external migration-index trend likelihood (0 disables)
  gamma_obs_scale = 500, # ratio of total external migration to HIV+ migration volume
  abroad_age_trend = 0, # logit-scale drift in born-abroad arrival age composition
  detect_step_quarter = NULL, # quarter of an ascertainment step in detect_prob
  detect_step_size = 0, # logit-scale size of that step
  gamma_obs_scale_trend = 0, # log-scale drift in gamma_obs_scale
  gamma_obs_drift_sd = 0.02, # quarterly log-scale random-walk innovation sd
  recent_arrival_completeness = 0
) {
  if (rita) {
    nlatst <- 7
    ndiagst <- 6
  } else {
    nlatst <- 4
    ndiagst <- 5
  }

  q_age <- matrix(rep(q, nage), nrow = nage, byrow = TRUE)

  if (is.null(sim_seed)) {
    sim_seed <- sample.int(.Machine$integer.max, 1)
  }

  if (!is.null(detect_obs_times)) {
    if (
      length(detect_obs_times) == 0 ||
        anyNA(detect_obs_times) ||
        any(detect_obs_times < 1) ||
        any(detect_obs_times != as.integer(detect_obs_times)) ||
        any(detect_obs_times > nquar)
    ) {
      stop(
        "detect_obs_times must be a non-empty integer vector with values in 1:nquar."
      )
    }
  }

  if (!is.null(rita_completeness_by_state)) {
    if (
      !is.numeric(rita_completeness_by_state) ||
        anyNA(rita_completeness_by_state) ||
        length(rita_completeness_by_state) != 5 ||
        any(rita_completeness_by_state < 0) ||
        any(rita_completeness_by_state > 1)
    ) {
      stop(
        "rita_completeness_by_state must be a numeric vector of length 5 with values in [0, 1]."
      )
    }
  }

  if (
    !is.logical(gamma_obs_grouped) ||
      anyNA(gamma_obs_grouped) ||
      !(length(gamma_obs_grouped) %in% c(1L, nquar))
  ) {
    stop(
      "gamma_obs_grouped must be a logical scalar or a logical vector of length nquar."
    )
  }

  if (!is.null(detect_step_quarter)) {
    if (
      length(detect_step_quarter) != 1 ||
        is.na(detect_step_quarter) ||
        detect_step_quarter != as.integer(detect_step_quarter) ||
        detect_step_quarter < 2 ||
        detect_step_quarter > nquar
    ) {
      stop("detect_step_quarter must be NULL or a single integer in 2:nquar.")
    }
  }

  if (
    length(recent_arrival_completeness) != 1L ||
      !is.finite(recent_arrival_completeness) ||
      recent_arrival_completeness < 0 ||
      recent_arrival_completeness > 1
  ) {
    stop("recent_arrival_completeness must be a number between 0 and 1.")
  }
  if (
    length(gamma_obs_drift_sd) != 1L ||
      !is.finite(gamma_obs_drift_sd) ||
      gamma_obs_drift_sd < 0
  ) {
    stop("gamma_obs_drift_sd must be a single non-negative number.")
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
    abroad_prop = abroad_prop,
    cd4_completeness = cd4_completeness,
    rita_completeness = rita_completeness,
    rita_completeness_by_state = rita_completeness_by_state,
    age_diag_effect = age_diag_effect,
    mig_diag_offset = mig_diag_offset,
    frr = frr,
    detect_obs_times = if (is.null(detect_obs_times)) {
      seq(4L, nquar, by = 8L)
    } else {
      as.integer(detect_obs_times)
    },
    gamma_obs_grouped = if (length(gamma_obs_grouped) == 1L) {
      rep(gamma_obs_grouped, nquar)
    } else {
      as.logical(gamma_obs_grouped)
    },
    gamma_obs_sd = gamma_obs_sd,
    gamma_obs_scale = gamma_obs_scale,
    abroad_age_trend = abroad_age_trend,
    detect_step_quarter = if (is.null(detect_step_quarter)) {
      NULL
    } else {
      as.integer(detect_step_quarter)
    },
    detect_step_size = detect_step_size,
    gamma_obs_scale_trend = gamma_obs_scale_trend,
    gamma_obs_drift_sd = gamma_obs_drift_sd,
    recent_arrival_completeness = recent_arrival_completeness
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
#' - "p_increasing": The born-abroad:UK-born UK-acquired incidence ratio increases over time.
#' - "m_differing": Migration probabilities (mu) differ by diagnosed state.
#' - "m_varying": Migration probabilities (mu) vary over time.
#' - "combo_1": Combination of varying infections, step change migration, and varying diagnosis probabilities.
#' - "combo_2": Combination of varying infections, step change migration, varying diagnosis probabilities, and an increasing UK-acquired incidence ratio.
#' - "combo_3": Combination of early-peaking infections, increasing migration, varying diagnosis probabilities, an increasing UK-acquired incidence ratio, varying migration probabilities, and age effects.
#' - "combo_4": Late-peaking infections (varying_2) with a late upward step change in migration (step_change_2),
#'   representing a surge in arrivals from high-prevalence regions. Shares varying diagnosis probabilities,
#'   an increasing UK-acquired incidence ratio, varying migration probabilities, and age effects with combo_3.
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
      # State-invariant offset: the model's alpha_diag_mig is a single scalar
      # shared across CD4 states, so the truth must be too.
      config$mig_diag_offset <- rep(-0.25, if (config$rita) 5 else 4)
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
      # State-invariant offset: the model's alpha_diag_mig is a single scalar
      # shared across CD4 states, so the truth must be too.
      config$mig_diag_offset <- rep(-0.25, if (config$rita) 5 else 4)
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
#' \donttest{
#' q_age <- matrix(rep(c(0.2, 0.3, 0.4, 0.5), 52), nrow = 52, byrow = TRUE)
#' d_probs <- array(0.1, dim = c(10, 52, 4))
#' cd4backcalc:::Q_mat(t = 10, d_probs, q_age, nlatst = 4, nage = 52, rita = FALSE)
#' d_probs <- array(0.1, dim = c(10, 52, 5))
#' cd4backcalc:::Q_mat(t = 10, d_probs, q_age, nlatst = 7, nage = 52, rita = TRUE)
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
#' \donttest{
#' q_age <- matrix(rep(c(0.2, 0.3, 0.4, 0.5), 52), nrow = 52, byrow = TRUE)
#' d_probs <- array(0.1, dim = c(10, 52, 4))
#' cd4backcalc:::D_mat(t = 10, d_probs = d_probs, q_age,
#'  nlatst = 4, ndiagst = 5, nage = 52, rita = FALSE)
#' d_probs <- array(0.1, dim = c(10, 52, 5))
#' cd4backcalc:::D_mat(t = 10, d_probs = d_probs, q_age,
#'  nlatst = 7, ndiagst = 6, nage = 52, rita = TRUE)
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
#' \donttest{
#' config <- cd4backcalc:::define_sim_config(nquar = 10)
#' q_age <- matrix(rep(c(0.0976, 0.115, 0.116, 0.148), 52), nrow = 52, byrow = TRUE)
#' d_probs <- array(0.1, dim = c(10, 52, 4))
#' Q_l <- lapply(1:10, function(x) {
#'    cd4backcalc:::Q_mat(x, d_probs, q_age, nlatst = 4)
#' })
#' D_l <- lapply(1:10, function(x) {
#'    cd4backcalc:::D_mat(x, d_probs, q_age, nlatst = 4, ndiagst = 5)
#' })
#' lat_init <- array(0, dim = c(52, 4))
#' lat_init[, 1] <- c(40, 30, 20, 10)
#' h <- array(matrix(rep(10, 20), nrow = 10), dim = c(52, 10, 4))
#' cd4backcalc:::epi_mat(Q_l, D_l, h, lat_init, config)
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
#' @param p_pattern Pattern for the born-abroad:UK-born ratio of UK-acquired
#'   incidence ("constant", "increasing").
#' @param m_pattern Pattern type for mu ("constant", "differing", or
#'   "varying").
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
#' @param abroad_prop Proportion of infections acquired abroad at the outset.
#' @param cd4_completeness CD4 completeness.
#' @param rita_completeness RITA completeness.
#' @param rita_completeness_by_state Optional numeric vector of length 5
#'   giving state-specific completeness for RITA simulations. Elements
#'   correspond to RITA positive, CD4 500 or above, CD4 350--499,
#'   CD4 200--349, and CD4 below 200.
#'   Use lower values in lower CD4 states to simulate informative missingness.
#' @param age_diag_effect Numeric scalar controlling the strength of age variation in diagnosis probabilities. Applied as a logit-additive offset scaled linearly across ages.
#' @param mig_diag_offset Numeric vector of per-state logit-scale offsets applied
#'   to diagnosis probabilities for people born abroad (`y = 2`). Length must
#'   equal the number of diagnosis states (4 for non-RITA, 5 for RITA). `NULL`
#'   (default) applies no offset.
#' @param frr False recency rate for RITA misclassification.
#' @param eta_sd Standard deviation of the approximate melding model for the
#'   log ratio of UK-acquired incidence among people born abroad to UK-acquired
#'   incidence among UK-born people. Smaller values anchor the model more
#'   tightly to the external \code{ratio_abroad_uk} evidence. Set to 0 when
#'   that evidence is not independent of outcomes in the likelihood. Default
#'   is 0.4.
#' @param detect_obs_sd Standard deviation, on the logit scale, of independent
#'   external evidence about the probability that a previous diagnosis abroad
#'   is detected. Set to 0 to disable this auxiliary likelihood.
#' @param detect_obs_times Integer vector of quarter indices (1-based) at which external
#'   detection-probability observations are generated. Defaults to \code{NULL}, which uses
#'   every 20 quarters (approx. every 5 years). Supply a custom vector, e.g. \code{c(4, 10, 30)},
#'   to place observations at specific quarters.
#' @param gamma_obs_grouped Logical scalar or logical vector of length
#'   \code{nquar} controlling which quarters receive age-stratified external
#'   migration observations on the same 4 age groups as \code{ABROAD_group}.
#'   The scalar \code{gamma_obs} series is always simulated; quarters where
#'   \code{gamma_obs_grouped} is \code{FALSE} have \code{-1} rows in
#'   \code{gamma_obs_group}.
#' @param gamma_obs_sd Standard deviation, on the log scale, of an independent
#'   rolling-annual migration index used by migration models.
#'   This evidence informs temporal trend, not absolute migration level. Set
#'   to 0 to disable this auxiliary likelihood.
#' @param gamma_obs_scale Ratio of total external migration to HIV-positive
#'   migration volume used to construct the external migration index.
#' @param abroad_age_trend Total logit-scale change (start to end) in the
#'   15-24/25-34 age-group shares of born-abroad arrivals relative to the 45+
#'   reference group; 0 = time-constant age composition.
#' @param detect_step_quarter Quarter index (1-based) of an ascertainment
#'   step in detect_prob (e.g. a surveillance reporting change). \code{NULL}
#'   (default) means no step.
#' @param detect_step_size Logit-scale jump in detect_prob applied from
#'   \code{detect_step_quarter} onward.
#' @param gamma_obs_scale_trend Total log-scale change (start to end) in the
#'   HIV-positive share of the external migration index (gamma_obs_scale);
#'   0 = constant share.
#' @param gamma_obs_drift_sd Quarterly random-walk innovation SD on the log
#'   ratio of the external migration index to HIV-positive migration. This
#'   controls stochastic variation in migrant HIV prevalence and is passed to
#'   the matching Stan model. Set to 0 to disable.
#' @param recent_arrival_completeness Proportion of eligible diagnoses with
#'   avidity testing done and a recorded year of arrival. Combined
#'   multiplicatively with \code{cd4_completeness} (a recorded CD4 value is
#'   also required) to give the overall completeness of the recent-arrival
#'   CD4 auxiliary data. Set to 0 to disable.
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
  abroad_prop = 0.2, # proportion acquired abroad at outset
  cd4_completeness = 1, # CD4 completeness
  rita_completeness = 1, # RITA completeness
  rita_completeness_by_state = NULL, # optional state-specific RITA completeness
  age_diag_effect = 0, # age effect on diagnosis probabilities
  mig_diag_offset = NULL, # per-state logit offsets for migrant diagnosis probs
  frr = 0.019, # false recency rate for RITA misclassification
  eta_sd = 0.4, # sd of the melding prior for log(ratio_abroad_uk)
  detect_obs_sd = 0,
  detect_obs_times = NULL,
  gamma_obs_grouped = FALSE,
  gamma_obs_sd = 0, # sd of the external migration-index trend likelihood (0 disables)
  gamma_obs_scale = 500, # ratio of total external migration to HIV+ migration volume
  abroad_age_trend = 0, # logit-scale drift in born-abroad arrival age composition
  detect_step_quarter = NULL, # quarter of an ascertainment step in detect_prob
  detect_step_size = 0, # logit-scale size of that step
  gamma_obs_scale_trend = 0, # log-scale drift in gamma_obs_scale
  gamma_obs_drift_sd = 0.02, # quarterly log-scale random-walk innovation sd
  recent_arrival_completeness = 0
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
    abroad_prop,
    cd4_completeness,
    rita_completeness,
    rita_completeness_by_state,
    age_diag_effect,
    mig_diag_offset,
    frr,
    detect_obs_times,
    gamma_obs_grouped,
    gamma_obs_sd,
    gamma_obs_scale,
    abroad_age_trend,
    detect_step_quarter,
    detect_step_size,
    gamma_obs_scale_trend,
    gamma_obs_drift_sd,
    recent_arrival_completeness
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

  # Detection probability for previous diagnoses abroad, with an optional
  # ascertainment step (detect_step_quarter/detect_step_size), e.g. a
  # surveillance reporting change.
  detect_logit_trend <- logit(seq(0.7, 0.8, length.out = sim_config$nquar))
  if (!is.null(sim_config$detect_step_quarter)) {
    step_idx <- sim_config$detect_step_quarter:sim_config$nquar
    detect_logit_trend[step_idx] <- detect_logit_trend[step_idx] +
      sim_config$detect_step_size
  }
  detect_prob <- if (sim_config$dist) {
    invlogit(stats::rnorm(sim_config$nquar, detect_logit_trend, 0.1))
  } else {
    invlogit(detect_logit_trend)
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
  # have a systematically different baseline diagnosis rate than UK-born people
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

  # Derive equilibrium initial state proportions analytically from t=1 diagnosis
  # probabilities and progression rates, so init_props always match the
  # d_pattern/mig_diag_offset used in this simulation call.
  sim_config$init_props <- derive_init_props(
    d_probs_t1 = d_probs[1, 1, ],
    q = sim_config$q,
    rita = sim_config$rita
  )
  sim_config$init_props_abroad <- derive_init_props(
    d_probs_t1 = d_probs_mig[1, 1, ],
    q = sim_config$q,
    mu_lat_t1 = mu[1, 1:sim_config$nlatst],
    abroad = TRUE,
    rita = sim_config$rita
  )

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

  # Generate Q and D matrices for the UK-born group (y = 1)
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

  # Generate Q and D matrices for the born-abroad group (y = 2)
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
  age_group_lookup_r <- pmin(4L, 1L + (seq_len(sim_config$nage) - 1L) %/% 10L)
  agedist_abroad <- generate_age_dist(
    sim_config$nquar,
    sim_config$nage,
    sim_config$a_pattern,
    20,
    40
  )
  # Optional linear trend shifts share from the 15-24/25-34 groups 
  # to 35-44/45+ over time
  if (sim_config$abroad_age_trend != 0) {
    t_std <- (seq_len(sim_config$nquar) - 0.5 * (sim_config$nquar + 1)) /
      sim_config$nquar
    group_trend <- outer(t_std, c(-1, -1, 1, 1) * sim_config$abroad_age_trend)
    agedist_abroad <- agedist_abroad * exp(group_trend[, age_group_lookup_r])
    agedist_abroad <- agedist_abroad / rowSums(agedist_abroad)
  }
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
  # detected counts use the detection probability; remaining counts are
  # allocated using the arrival-age distribution.
  if (sim_config$dist) {
    total_diag <- as.integer(m[, sim_config$nlatst + 1])
    detected_total <- stats::rbinom(sim_config$nquar, total_diag, detect_prob)
    ABROAD_age_mat <- matrix(0L, sim_config$nage, sim_config$nquar)
    misclassified_age <- matrix(0L, sim_config$nage, sim_config$nquar)
    for (t in seq_len(sim_config$nquar)) {
      ABROAD_age_mat[, t] <- stats::rmultinom(
        1,
        detected_total[t],
        agedist_abroad[t, ]
      )
      misclassified_age[, t] <- stats::rmultinom(
        1,
        total_diag[t] - detected_total[t],
        agedist_abroad[t, ]
      )
    }
  } else {
    total_diag <- m[, sim_config$nlatst + 1]
    detected_total <- round(total_diag * detect_prob)
    ABROAD_age_mat <- vapply(seq_len(sim_config$nquar), function(t) {
      round_preserve(detected_total[t] * agedist_abroad[t, ])
    }, numeric(sim_config$nage))
    misclassified_age <- vapply(seq_len(sim_config$nquar), function(t) {
      round_preserve((total_diag[t] - detected_total[t]) * agedist_abroad[t, ])
    }, numeric(sim_config$nage))
  }
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

  recent_arrival_cd4 <- if (!sim_config$migration) {
    empty_recent_arrival_cd4 <- array(
      0,
      dim = c(sim_config$nquar %/% 4L, 4L, 4L)
    )
    list(
      observed = array(0L, dim = dim(empty_recent_arrival_cd4)),
      expected = empty_recent_arrival_cd4
    )
  } else {
    generate_recent_arrival_cd4(
      m_age = m_age[,, if (sim_config$rita) 4:7 else 1:4, drop = FALSE],
      d_probs_mig = d_probs_mig[,, if (sim_config$rita) 2:5 else 1:4, drop = FALSE],
      q_age = sim_config$q_age,
      completeness = 0,
      dist = FALSE
    )
  }

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
    "misclassified_age" = misclassified_age,
    "recent_arrival_cd4_expected" = recent_arrival_cd4$expected
  )

  # Perform validation checks before rounding
  validate_sim_data(sim_data)

  # Apply rounding which preserves sums
  epi_uk <- round_fields(
    sim_data$epi_uk,
    stochastic_diags_age = sim_config$dist
  )
  epi_uk_abroad <- round_fields(
    sim_data$epi_uk_abroad,
    stochastic_diags_age = sim_config$dist
  )
  epi_abroad <- round_fields(
    sim_data$epi_abroad,
    stochastic_diags_age = sim_config$dist
  )

  if (sim_config$dist) {
    sim_data$epi_uk <- epi_uk
    sim_data$epi_uk_abroad <- epi_uk_abroad
    sim_data$epi_abroad <- epi_abroad
  }

  # Prepare Stan data
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

  # Misclassified previous diagnoses abroad enter born-abroad HIV totals
  # immediately, without CD4 progression or diagnosis probabilities.
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
    # Sum over country-of-birth groups
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
  # detect_obs_times controls which quarters have observations.
  detect_obs_index <- if (detect_obs_sd > 0) {
    as.integer(sim_config$detect_obs_times)
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

  # Time-varying share of the external index that is HIV-positive migration
  t_std_gamma <- (seq_len(sim_config$nquar) - 0.5 * (sim_config$nquar + 1)) /
    sim_config$nquar
  gamma_obs_drift <- numeric(sim_config$nquar)
  if (sim_config$dist && sim_config$gamma_obs_drift_sd > 0) {
    gamma_obs_drift[-1L] <- cumsum(stats::rnorm(
      sim_config$nquar - 1L,
      0,
      sim_config$gamma_obs_drift_sd
    ))
    gamma_obs_drift <- gamma_obs_drift - mean(gamma_obs_drift)
  }
  gamma_obs_scale_t <- sim_config$gamma_obs_scale * exp(
    sim_config$gamma_obs_scale_trend * t_std_gamma + gamma_obs_drift
  )
  sim_data$migrant_hiv_prevalence <- 1 / gamma_obs_scale_t

  # External migration-index auxiliary data (e.g. ONS total immigration trend)
  # by the same 4 age groups used by ABROAD_group. gamma_obs_sd = 0 disables the
  # auxiliary likelihood in Stan
  gamma_obs <- rep(0, sim_config$nquar)
  gamma_obs_group <- matrix(-1, nrow = sim_config$nquar, ncol = 4L)
  if (gamma_obs_sd > 0) {
    omega_trend <- generate_omega(
      sim_config$nquar,
      sim_config$o_pattern,
      sim_config$target_mig,
      dist = FALSE
    )
    omega_start_slope <- if (length(omega_trend) > 1L) {
      omega_trend[2L] - omega_trend[1L]
    } else {
      0
    }
    omega_pre <- pmax(
      .Machine$double.eps,
      omega_trend[1L] - 3:1 * omega_start_slope
    )
    omega_ext <- c(omega_pre, omega_trend)
    rolling_vol <- vapply(
      seq_len(sim_config$nquar),
      function(t) sum(omega_ext[t:(t + 3L)]),
      numeric(1)
    )

    gamma_obs <- stats::rnorm(
      seq_len(sim_config$nquar),
      log(pmax(1, round(gamma_obs_scale_t * rolling_vol, -3))),
      gamma_obs_sd
    )

    if (any(sim_config$gamma_obs_grouped)) {
      omega_age_trend <- sweep(agedist_abroad, 1, omega_trend, "*")
      omega_group_trend <- matrix(0, nrow = sim_config$nquar, ncol = 4L)
      for (ag in 1:4) {
        idx <- which(age_group_lookup_r == ag)
        omega_group_trend[, ag] <- rowSums(omega_age_trend[, idx, drop = FALSE])
      }
      omega_group_start_slope <- if (sim_config$nquar > 1L) {
        omega_group_trend[2L, ] - omega_group_trend[1L, ]
      } else {
        rep(0, 4L)
      }
      omega_group_pre <- matrix(0, nrow = 3L, ncol = 4L)
      for (ag in 1:4) {
        omega_group_pre[, ag] <- pmax(
          .Machine$double.eps,
          omega_group_trend[1L, ag] - 3:1 * omega_group_start_slope[ag]
        )
      }
      omega_group_ext <- rbind(omega_group_pre, omega_group_trend)
      rolling_vol_group <- matrix(0, nrow = sim_config$nquar, ncol = 4L)
      for (ag in 1:4) {
        rolling_vol_group[, ag] <- vapply(
          seq_len(sim_config$nquar),
          function(t) sum(omega_group_ext[t:(t + 3L), ag]),
          numeric(1)
        )
      }
      gamma_obs_group_full <- matrix(
        stats::rnorm(
          sim_config$nquar * 4L,
          as.vector(log(pmax(
            1,
            round(gamma_obs_scale_t * rolling_vol_group, -3)
          ))),
          gamma_obs_sd
        ),
        nrow = sim_config$nquar,
        ncol = 4L
      )
      gamma_obs_group[sim_config$gamma_obs_grouped, ] <-
        gamma_obs_group_full[sim_config$gamma_obs_grouped, , drop = FALSE]
    }
  }

  # ABROAD data: detected portion of diagnosed migrants
  ABROAD <- colSums(ABROAD_age_mat)
  ABROAD_age <- t(ABROAD_age_mat) # [nquar, nage], kept for R-side use

  # Pre-aggregate ABROAD into the 4 diagnosis age groups used by the Stan model
  # (1-10 -> 1, 11-20 -> 2, 21-30 -> 3, 31+ -> 4), matching age_group_lookup.
  # This is the sufficient statistic for the multinomial age-composition likelihood;
  # the raw single-year ABROAD_age is retained in stan_data for diagnostics/plots.
  ABROAD_group <- matrix(0L, nrow = sim_config$nquar, ncol = 4L)
  for (ag in 1:4) {
    idx <- which(age_group_lookup_r == ag)
    ABROAD_group[, ag] <- as.integer(rowSums(ABROAD_age[, idx, drop = FALSE]))
  }

  # Make stan data list
  stan_data <- list(
    "HIV" = hiv_dx,
    "HIV_age" = hiv_dx_age,
    "AIDS" = aids_dx,
    "AIDS_age" = aids_dx_age,
    "CD4" = cd4_dx,
    "CD4_age" = cd4_dx_age,
    "ABROAD" = ABROAD,
    "ABROAD_group" = ABROAD_group,
    "ABROAD_age" = ABROAD_age, # retained for R-side diagnostics; not passed to Stan
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
    "detect_step_quarter" = if (is.null(sim_config$detect_step_quarter)) {
      integer(0)
    } else {
      as.integer(sim_config$detect_step_quarter)
    },
    "gamma_obs" = gamma_obs,
    "gamma_obs_group" = gamma_obs_group,
    "gamma_obs_sd" = gamma_obs_sd,
    "gamma_obs_drift_sd" = sim_config$gamma_obs_drift_sd,
    "sim_config" = sim_config
  )
  if (sim_config$migration) {
    stan_data$nyear <- sim_config$nquar %/% 4L
  }

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
  # Adjust RITA counts for completeness; optionally use state-specific values
  # to induce informative missingness (e.g. lower completeness at low CD4).
  if (sim_config$rita) {
    rita_state_probs <- if (is.null(sim_config$rita_completeness_by_state)) {
      rep(sim_config$rita_completeness, sim_config$ndiagst - 1)
    } else {
      as.numeric(sim_config$rita_completeness_by_state)
    }

    if (any(rita_state_probs < 1)) {
      apply_rita_completeness <- function(x, state_probs, dist) {
        dims <- dim(x)
        nstates <- dims[length(dims)]
        if (nstates != length(state_probs)) {
          stop(
            "RITA completeness vector length must match CD4 state dimension."
          )
        }

        n_per_state <- prod(dims[-length(dims)])
        probs <- rep(state_probs, each = n_per_state)
        counts <- pmax(0L, as.integer(round(as.vector(x))))
        out <- if (dist) {
          stats::rbinom(length(counts), counts, probs)
        } else {
          round(counts * probs)
        }
        array(out, dim = dims)
      }

      stan_data$CD4 <- apply_rita_completeness(
        stan_data$CD4,
        rita_state_probs,
        sim_config$dist
      )
      stan_data$CD4_age <- apply_rita_completeness(
        stan_data$CD4_age,
        rita_state_probs,
        sim_config$dist
      )
    }
  }

  if (sim_config$migration) {
    expected_total <- aggregate_migrant_cd4(cd4_dx_age, sim_config$rita)
    observed_total <- aggregate_migrant_cd4(
      stan_data$CD4_age,
      sim_config$rita
    )
    expected_recent <- recent_arrival_cd4$expected

    if (!sim_config$age) {
      expected_total <- apply(expected_total, c(1, 3), sum)
      observed_total <- apply(observed_total, c(1, 3), sum)
      expected_recent <- apply(expected_recent, c(1, 3), sum)
    }

    eligible_total <- if (sim_config$dist) {
      array(
        stats::rbinom(
          length(observed_total),
          as.integer(observed_total),
          sim_config$recent_arrival_completeness
        ),
        dim = dim(observed_total)
      )
    } else {
      round(observed_total * sim_config$recent_arrival_completeness)
    }
    recent_prob <- pmin(
      1,
      expected_recent / pmax(expected_total, .Machine$double.eps)
    )
    stan_data$recent_arrival_cd4_total <- eligible_total
    stan_data$recent_arrival_cd4 <- if (sim_config$dist) {
      array(
        stats::rbinom(
          length(eligible_total),
          as.integer(eligible_total),
          recent_prob
        ),
        dim = dim(eligible_total)
      )
    } else {
      round(eligible_total * recent_prob)
    }
  }

  return(list(
    "stan_data" = stan_data,
    "sim_data" = sim_data
  ))
}
