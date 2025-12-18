#' Function to get quantiles
#'
#' @param x Numeric vector
#' @return A named numeric vector with 2.5%, 50%, mean, and 97.5% quantiles
get_quantiles <- function(x) {
  x <- as.numeric(x)
  qs <- stats::quantile(
    x,
    probs = c(0.025, 0.5, 0.975),
    na.rm = TRUE,
    names = FALSE,
    type = 7
  )
  out <- c(qs[1], qs[2], mean(x, na.rm = TRUE), qs[3])
  names(out) <- c("2.5%", "50%", "mean", "97.5%")
  out
}

.validate_infs_dims <- function(infs, nt, nage) {
  if (!is.matrix(infs)) {
    stop("infs must be a matrix (draws x (nt*nage)).")
  }
  if (length(nt) != 1 || is.na(nt) || nt <= 0) {
    stop("nt must be a positive scalar.")
  }
  if (length(nage) != 1 || is.na(nage) || nage <= 0) {
    stop("nage must be a positive scalar.")
  }
  if (ncol(infs) != nt * nage) {
    stop(sprintf(
      "Expected ncol(infs) == nt*nage (%d), got %d.",
      nt * nage,
      ncol(infs)
    ))
  }
  invisible(TRUE)
}

.age_to_idx <- function(age.st, age.end, nage, offset = 14) {
  if (length(age.st) != 1 || length(age.end) != 1) {
    stop("age.st and age.end must be scalars.")
  }
  if (age.end < age.st) {
    stop("age.end must be >= age.st.")
  }
  a1 <- age.st - offset
  a2 <- age.end - offset
  if (a1 < 1 || a2 > nage) {
    stop(sprintf(
      "Age range [%d,%d] maps to indices [%d,%d] outside 1:%d.",
      age.st,
      age.end,
      a1,
      a2,
      nage
    ))
  }
  a1:a2
}

.aggregate_quarters_to_years <- function(x, nt) {
  if (nt %% 4 != 0) {
    stop("nt must be divisible by 4 when aggregating quarters to years.")
  }
  grp <- rep(seq_len(nt / 4), each = 4)
  as.numeric(rowsum(x, grp, reorder = FALSE))
}

#' Function to get posterior functions of incidence
#'
#' This function generates incidence by age group from log-infections (per each posterior iteration).
#' This function can aggregates incidence by year, depending on qt.flag.
#'
#' @param infs A matrix with quarterly log-infections for each posterior sample in each row.
#' @param age.st Starting age for incidence calculation.
#' @param age.end Ending age for incidence calculation.
#' @param nt Total number of time points (quarters).
#' @param type Type of incidence calculation, either "time" or "age".
#' @param qt.flag Logical flag indicating whether to aggregate quarterly infections into yearly.
#' @param nage Number of age groups.
#' @return A matrix with summary statistics (2.5%, 50%, mean, 97.5%) for each time point or age group.
infs.fct <- function(
  infs,
  age.st,
  age.end,
  nt,
  type = "time",
  qt.flag = FALSE,
  nage = 52
) {
  if (!type %in% c("time", "age")) {
    stop("type can only be 'time' or 'age'.")
  }
  .validate_infs_dims(infs, nt, nage)
  age_idx <- .age_to_idx(age.st, age.end, nage = nage)

  exp_infs <- exp(infs)

  if (type == "time") {
    x <- t(apply(exp_infs, 1, function(v) {
      mat <- matrix(v, nrow = nt, ncol = nage)
      ts <- rowSums(mat[, age_idx, drop = FALSE])
      if (qt.flag) .aggregate_quarters_to_years(ts, nt) else ts
    }))
  } else {
    x <- t(apply(exp_infs, 1, function(v) {
      mat <- matrix(v, nrow = nt, ncol = nage)
      colSums(mat[, age_idx, drop = FALSE])
    }))
  }

  apply(x, 2, get_quantiles)
}

#' Function to get posterior functions of incidence
#'
#' This function generates incidence by age group from log-infections (per each posterior iteration).
#' This function always aggregates quarterly incidence curve.
#'
#' @param infs A matrix with quarterly log-infections for each posterior sample in each row.
#' @param age.st Starting age for incidence calculation.
#' @param age.end Ending age for incidence calculation.
#' @param nt Total number of time points (quarters).
#' @param qt.flag Logical flag indicating whether to aggregate quarterly infections into yearly.
#' @param nage Number of age groups.
#' @return A matrix with summary statistics (2.5%, 50%, mean, 97.5%) for each time point.
infs.fct.qt <- function(infs, age.st, age.end, nt, qt.flag = FALSE, nage = 52) {
  age_idx <- .age_to_idx(age.st, age.end, nage = nage)

  if (!qt.flag) {
    # infs is on log scale; to split yearly total into 4 equal quarters:
    # log(I/4) = log(I) - log(4)
    infs <- infs[, rep(seq_len(ncol(infs)), each = 4), drop = FALSE] - log(4)
    nt <- nt * 4
  }

  .validate_infs_dims(infs, nt, nage)
  exp_infs <- exp(infs)

  x <- t(apply(exp_infs, 1, function(v) {
    mat <- matrix(v, nrow = nt, ncol = nage)
    rowSums(mat[, age_idx, drop = FALSE])
  }))

  apply(x, 2, get_quantiles)
}

#' Function to get posterior functions of incidence rate
#'
#' This function generates incidence rate by age class from log-infections (per each posterior iteration).
#' This function can aggregates incidence by year, depending on qt.flag.
#'
#' @param infs A matrix with quarterly log-infections for each posterior sample in each row.
#' @param age.cl Age class for incidence rate calculation. Options are "Tot", "15-34", "35-44", "45+".
#' @param last.yr Number of last years to consider for incidence rate calculation.
#' @param nt Total number of time points (quarters).
#' @param qt.flag Logical flag indicating whether to aggregate quarterly infections into yearly.
#' @param totals A list containing total population counts for different age classes.
#' @param nage Number of age groups.
#' @return A matrix with summary statistics (2.5%, 50%, mean, 97.5%) for each time point.
inc.fct <- function(
  infs,
  age.cl,
  last.yr,
  nt,
  qt.flag = FALSE,
  totals,
  nage = 52
) {
  ## infs is a matrix with quarterly log-infections for each posterior sample in each row
  ## qt.flag=TRUE aggregates quarterly infections and makes them yearly
  if (!age.cl %in% c("Tot", "15-34", "35-44", "45+")) {
    stop("age.cl can only be 'Tot','15-34', '35-44','45+' ")
  }
  .validate_infs_dims(infs, nt, nage)

  if (age.cl == "Tot") {
    a1 <- 1
    a2 <- 52
    den <- totals$N.msm.tot
  }
  if (age.cl == "15-34") {
    a1 <- 1
    a2 <- 20
    den <- totals$N.msm.1534
  }
  if (age.cl == "35-44") {
    a1 <- 21
    a2 <- 30
    den <- totals$N.msm.3544
  }
  if (age.cl == "45+") {
    a1 <- 31
    a2 <- 52
    den <- totals$N.msm.45
  }

  if (length(den) != 1 || is.na(den) || den <= 0) {
    stop("Population denominator must be a positive scalar.")
  }

  exp_infs <- exp(infs)

  if (!qt.flag) {
    if (last.yr > nt) {
      stop("last.yr cannot exceed nt when qt.flag=FALSE.")
    }
    ind <- (nt - last.yr + 1):nt
    x <- t(apply(exp_infs, 1, function(v) {
      mat <- matrix(v, nrow = nt, ncol = nage)
      rowSums(mat[ind, a1:a2, drop = FALSE]) / den
    }))
  } else {
    if (nt %% 4 != 0) {
      stop("nt must be divisible by 4 when qt.flag=TRUE.")
    }
    if (last.yr > (nt / 4)) {
      stop("Specified last.yr is too large for quarterly nt.")
    }
    ind <- (nt - 4 * last.yr + 1):nt
    grp <- rep(seq_len(last.yr), each = 4)

    x <- t(apply(exp_infs, 1, function(v) {
      mat <- matrix(v, nrow = nt, ncol = nage)
      ts <- rowSums(mat[ind, a1:a2, drop = FALSE])
      as.numeric(rowsum(ts, grp, reorder = FALSE)) / den
    }))
  }

  apply(x, 2, get_quantiles)
}

#' Function to identify the peak in infection in the last x years
#'
#' @param infs A matrix with quarterly log-infections for each posterior sample in each row.
#' @param age.st Starting age for peak identification.
#' @param age.end Ending age for peak identification.
#' @param nt Total number of time points (quarters).
#' @param last.yr Number of last years to consider for peak identification.
#' @param qt.flag Logical flag indicating whether to aggregate quarterly infections into yearly.
#' @param nage Number of age groups.
#' @return A numeric vector indicating the peak position for each posterior sample.
peak.fct <- function(
  infs,
  age.st,
  age.end,
  nt,
  last.yr = 10,
  qt.flag = FALSE,
  nage = 52
) {
  .validate_infs_dims(infs, nt, nage)
  age_idx <- .age_to_idx(age.st, age.end, nage = nage)

  exp_infs <- exp(infs)

  if (!qt.flag) {
    if (last.yr > nt) {
      stop("Specified last.yr is too large")
    }
    ind <- (nt - last.yr + 1):nt
    x <- t(apply(exp_infs, 1, function(v) {
      mat <- matrix(v, nrow = nt, ncol = nage)
      rowSums(mat[ind, age_idx, drop = FALSE])
    }))
  } else {
    if (nt %% 4 != 0) {
      stop("nt must be divisible by 4 when qt.flag=TRUE.")
    }
    if (last.yr > (nt / 4)) {
      stop("Specified last.yr is too large")
    }
    ind <- (nt - 4 * last.yr + 1):nt
    grp <- rep(seq_len(last.yr), each = 4)

    x <- t(apply(exp_infs, 1, function(v) {
      mat <- matrix(v, nrow = nt, ncol = nage)
      ts <- rowSums(mat[ind, age_idx, drop = FALSE])
      as.numeric(rowsum(ts, grp, reorder = FALSE))
    }))
  }

  apply(x, 1, which.max)
}

#' Function for the shape of incidence curve
#'
#' This function classifies the shape of the incidence curve in the last two years as either monotonically increasing, increasing (not monotonic), decreasing (not monotonic), or monotonically decreasing.
#' @param x A numeric vector representing the incidence curve.
#' @return A character string indicating the shape of the incidence curve.
shape.fct <- function(x) {
  ## x arbitrary vector...
  ## do not consider flat case here as not really the case...
  l <- length(x)
  o <- NA
  if (l <= 1) {
    stop("length of x must be greater than 1")
  } else if (l == 2) {
    if (x[1] < x[l]) o <- "inc" else o <- "dec"
  } else {
    if (x[1] < x[l]) {
      if (all(x == cummax(x))) {
        o <- "mon inc"
      } else {
        o <- "inc (not mon)"
      }
    } else {
      if (all(x == cummin(x))) {
        o <- "mon dec"
      } else {
        o <- "dec (not mon)"
      }
    }
  }
  return(o)
}

#' Function to identify the trend in infection in the last x years
#'
#' @param infs A matrix with quarterly log-infections for each posterior sample in each row.
#' @param age.st Starting age for trend identification.
#' @param age.end Ending age for trend identification.
#' @param nt Total number of time points (quarters).
#' @param last.yr Number of last years to consider for trend identification.
#' @param qt.flag Logical flag indicating whether to aggregate quarterly infections into yearly.
#' @param nage Number of age groups.
#' @return A character vector indicating the trend for each posterior sample.
trend.fct <- function(
  infs,
  age.st,
  age.end,
  nt,
  last.yr = 3,
  qt.flag = FALSE,
  nage = 52
) {
  .validate_infs_dims(infs, nt, nage)
  age_idx <- .age_to_idx(age.st, age.end, nage = nage)

  exp_infs <- exp(infs)

  if (!qt.flag) {
    if (last.yr > nt) {
      stop("Specified last.yr is too large")
    }
    ind <- (nt - last.yr + 1):nt
    x <- t(apply(exp_infs, 1, function(v) {
      mat <- matrix(v, nrow = nt, ncol = nage)
      rowSums(mat[ind, age_idx, drop = FALSE])
    }))
  } else {
    if (nt %% 4 != 0) {
      stop("nt must be divisible by 4 when qt.flag=TRUE.")
    }
    if (last.yr > (nt / 4)) {
      stop("Specified last.yr is too large")
    }
    ind <- (nt - 4 * last.yr + 1):nt
    grp <- rep(seq_len(last.yr), each = 4)

    x <- t(apply(exp_infs, 1, function(v) {
      mat <- matrix(v, nrow = nt, ncol = nage)
      ts <- rowSums(mat[ind, age_idx, drop = FALSE])
      as.numeric(rowsum(ts, grp, reorder = FALSE))
    }))
  }

  apply(x, 1, shape.fct)
}

#' Function to get posterior functions of goodness of fit
#'
#' This function generates goodness of fit statistics for HIV, AIDS, and CD4 diagnoses by age group from their respective posterior samples (per each posterior iteration).
#' This function can aggregates diagnoses by year, depending on qt.flag.
#' @param type A character string indicating the type of diagnosis ("HIV", "AIDS", "CD4.st1", "CD4.st2", "CD4.st3", "CD4.st4").
#' @param age.st Starting age for goodness of fit calculation.
#' @param age.end Ending age for goodness of fit calculation.
#' @param HIV A 3D array of posterior samples for HIV diagnoses.
#' @param AIDS A 3D array of posterior samples for AIDS diagnoses.
#' @param CD4.st1 A 3D array of posterior samples for CD4 stage 1 diagnoses.
#' @param CD4.st2 A 3D array of posterior samples for CD4 stage 2 diagnoses.
#' @param CD4.st3 A 3D array of posterior samples for CD4 stage 3 diagnoses.
#' @param CD4.st4 A 3D array of posterior samples for CD4 stage 4 diagnoses.
#' @param CD4.st5 A 3D array of posterior samples for CD4 stage 5 diagnoses (optional).
#' @param nt Total number of time points (quarters).
#' @param qt.flag Logical flag indicating whether to aggregate quarterly diagnoses into yearly.
#' @return A matrix with summary statistics (2.5%, 50%, mean, 97.5%) for each time point or age group.
gof.fct <- function(
  type = "HIV",
  age.st,
  age.end,
  HIV,
  AIDS,
  CD4.st1,
  CD4.st2,
  CD4.st3,
  CD4.st4,
  CD4.st5 = NULL,
  nt,
  qt.flag
) {
  if (
    !type %in%
      c("HIV", "AIDS", "CD4.st1", "CD4.st2", "CD4.st3", "CD4.st4", "CD4.st5")
  ) {
    stop(
      "Type can only be 'HIV','AIDS','CD4.st1','CD4.st2','CD4.st3','CD4.st4','CD4.st5'"
    )
  }

  arr <- switch(
    type,
    HIV = HIV,
    AIDS = AIDS,
    CD4.st1 = CD4.st1,
    CD4.st2 = CD4.st2,
    CD4.st3 = CD4.st3,
    CD4.st4 = CD4.st4,
    CD4.st5 = CD4.st5
  )
  if (is.null(arr)) {
    stop(sprintf("Requested type '%s' is NULL.", type))
  }

  # Expect arr dims: (t, age, draw)
  nage <- dim(arr)[2]
  age_idx <- .age_to_idx(age.st, age.end, nage = nage)

  arr_sub <- arr[, age_idx, , drop = FALSE]

  # Sum over age -> matrix (t, draw), then transpose -> (draw, t)
  y <- t(apply(arr_sub, c(1, 3), sum))

  if (qt.flag) {
    if (nt %% 4 != 0) {
      stop("nt must be divisible by 4 when qt.flag=TRUE.")
    }
    grp <- rep(seq_len(nt / 4), each = 4)
    y <- t(apply(y, 1, function(v) as.numeric(rowsum(v, grp, reorder = FALSE))))
  }

  apply(y, 2, get_quantiles)
}

#' Function to get posterior functions of prevalence
#'
#' This function generates prevalence by age group from prevalence arrays (per each posterior iteration).
#'
#' @param prev.arr A 4D array with prevalence for each posterior sample.
#' @param state State for prevalence calculation, either "all" or 1 to 7.
#' @param age.st Starting age for prevalence calculation.
#' @param age.end Ending age for prevalence calculation.
#' @param rita.flag Logical flag indicating if RITA diagnoses are considered.
#' @return A matrix with summary statistics (2.5%, 50%, mean, 97.5%) for each time point.
prev.fct <- function(prev.arr, state, age.st, age.end, rita.flag) {
  if (!state %in% c("all", 1:7)) {
    stop("State can only be 'all' or 1:7")
  }

  nage <- dim(prev.arr)[2]
  age_idx <- .age_to_idx(age.st, age.end, nage = nage)

  states <- if (identical(state, "all")) {
    if (rita.flag) 1:7 else 1:4
  } else {
    state
  }

  # prev.arr dims expected: (t, age, state, draw)
  prev_sub <- prev.arr[,, states, , drop = FALSE]
  prev <- apply(prev_sub, c(1, 2, 4), sum) # (t, age, draw)

  x <- t(apply(prev[, age_idx, , drop = FALSE], c(1, 3), sum)) # (draw, t)
  apply(x, 2, get_quantiles)
}

#' Function to get diagnosis matrix per iteration
#'
#' @param fit.mat A matrix of posterior samples from the fitted model.
#' @param iter The iteration number to extract data for.
#' @param age.dx.flag Logical flag indicating if age-specific diagnoses are considered.
#' @param age.dx.spl.flag Logical flag indicating if age-specific diagnoses with splines are considered.
#' @param rita.flag Logical flag indicating if RITA diagnoses are considered.
#' @param nt Total number of time points.
#' @param nage Number of age groups.
#' @param ind_list A list of indices for diagnosis data extraction.
#' @return A diagnosis matrix for the specified iteration.
d.mat.iter <- function(
  fit.mat,
  iter,
  age.dx.flag,
  age.dx.spl.flag,
  rita.flag,
  nt,
  nage,
  ind_list
) {
  if (!age.dx.flag) {
    d.mat <- matrix(NA, nrow = 4, ncol = nt)
    d.mat[1, ] <- fit.mat[iter, ind_list$d1.ind]
    d.mat[2, ] <- fit.mat[iter, ind_list$d2.ind]
    d.mat[3, ] <- fit.mat[iter, ind_list$d3.ind]
    d.mat[4, ] <- fit.mat[iter, ind_list$d4.ind]
  } else {
    if (!age.dx.spl.flag) {
      if (rita.flag) {
        d.mat <- array(NA, dim = c(5, nt, 4))
        d.mat[5, , 1] <- fit.mat[iter, ind_list$d5.ind.y]
        d.mat[5, , 2] <- fit.mat[iter, ind_list$d5.ind.my]
        d.mat[5, , 3] <- fit.mat[iter, ind_list$d5.ind.mo]
        d.mat[5, , 4] <- fit.mat[iter, ind_list$d5.ind.o]
      } else {
        d.mat <- array(NA, dim = c(4, nt, 4)) ## array to be consistent with gof_iter_fct in stancppfcts_simpleperv_agediag.cpp
      }
      d.mat[1, , 1] <- fit.mat[iter, ind_list$d1.ind.y]
      d.mat[2, , 1] <- fit.mat[iter, ind_list$d2.ind.y]
      d.mat[3, , 1] <- fit.mat[iter, ind_list$d3.ind.y]
      d.mat[4, , 1] <- fit.mat[iter, ind_list$d4.ind.y]
      d.mat[1, , 2] <- fit.mat[iter, ind_list$d1.ind.my]
      d.mat[2, , 2] <- fit.mat[iter, ind_list$d2.ind.my]
      d.mat[3, , 2] <- fit.mat[iter, ind_list$d3.ind.my]
      d.mat[4, , 2] <- fit.mat[iter, ind_list$d4.ind.my]
      d.mat[1, , 3] <- fit.mat[iter, ind_list$d1.ind.mo]
      d.mat[2, , 3] <- fit.mat[iter, ind_list$d2.ind.mo]
      d.mat[3, , 3] <- fit.mat[iter, ind_list$d3.ind.mo]
      d.mat[4, , 3] <- fit.mat[iter, ind_list$d4.ind.mo]
      d.mat[1, , 4] <- fit.mat[iter, ind_list$d1.ind.o]
      d.mat[2, , 4] <- fit.mat[iter, ind_list$d2.ind.o]
      d.mat[3, , 4] <- fit.mat[iter, ind_list$d3.ind.o]
      d.mat[4, , 4] <- fit.mat[iter, ind_list$d4.ind.o]
    } else {
      d.mat <- array(NA, dim = c(4, nt, nage))
      d.mat[1, , ] <- matrix(
        fit.mat[iter, ind_list$d1.ind],
        nrow = nt,
        ncol = nage
      )
      d.mat[2, , ] <- matrix(
        fit.mat[iter, ind_list$d2.ind],
        nrow = nt,
        ncol = nage
      )
      d.mat[3, , ] <- matrix(
        fit.mat[iter, ind_list$d3.ind],
        nrow = nt,
        ncol = nage
      )
      d.mat[4, , ] <- matrix(
        fit.mat[iter, ind_list$d4.ind],
        nrow = nt,
        ncol = nage
      )
    }
  }

  return(d.mat)
}

#' Function to deal with under-reporting of AIDS diagnoses
#'
#' @param urep Under-reporting factor.
#' @param AIDS A matrix of AIDS diagnoses.
#' @param urep.flag Logical flag indicating if under-reporting is considered.
#' @param qt.flag Logical flag indicating if quarterly data is considered.
#' @param nt Total number of time points.
#' @return A matrix of AIDS diagnoses adjusted for under-reporting if applicable.
urep.fct <- function(urep, AIDS, urep.flag, qt.flag, nt) {
  out <- matrix(NA_real_, nrow = nrow(AIDS), ncol = ncol(AIDS))

  # Under-reporting starts from year 2000:
  # yearly model: year index 5 (if start is 1996-like), quarterly model: quarter 20
  urep.start <- if (qt.flag) 20 else 5
  urep.start <- min(urep.start, nt)

  if (urep.flag) {
    out[1:urep.start, ] <- AIDS[1:urep.start, , drop = FALSE]
    if (urep.start < nt) {
      out[(urep.start + 1):nt, ] <- urep *
        AIDS[(urep.start + 1):nt, , drop = FALSE]
    }
  } else {
    out <- AIDS
  }

  out
}
