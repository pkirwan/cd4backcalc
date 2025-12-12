#' Function to get quantiles
#'
#' @param x Numeric vector
#' @return A named numeric vector with 2.5%, 50%, mean, and 97.5% quantiles
get_quantiles <- function(x) {
  ## gives 95% credible intervals
  ## median and mode
  out <- rep(NA, 4)
  out[1] <- stats::quantile(x, probs = 0.025)
  out[2] <- stats::quantile(x, probs = 0.5)
  out[3] <- mean(x)
  out[4] <- stats::quantile(x, probs = 0.975)
  names(out) <- c("2.5%", "50%", "mean", "97.5%")
  return(out)
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
#' @return A matrix with summary statistics (2.5%, 50%, mean, 97.5%) for each time point or age group.
infs.fct <- function(
  infs,
  age.st,
  age.end,
  nt,
  type = "time",
  qt.flag = FALSE
) {
  ## infs is a matrix with quarterly log-infections for each posterior sample in each row
  ## qt.flag=TRUE aggregates quarterly infections and makes them yearly
  if (!type %in% c("time", "age")) {
    stop("Type can only be 'time' or 'age'")
  }
  a1 <- age.st - 14
  a2 <- age.end - 14
  if (qt.flag) {
    qt.indx <- rep(1:(nt / 4), each = 4)
  }
  if (type == "time") {
    if (!qt.flag) {
      call <- "rowSums(matrix(x, nr=nt)[,a1:a2])"
    } else {
      call <- "tapply(rowSums(matrix(x, nr=nt)[,a1:a2]), qt.indx, sum)"
    }
  } else {
    call <- "colSums(matrix(x, nr=nt)[,a1:a2])"
  }
  x <- t(apply(exp(infs), 1, function(x) eval(parse(text = call))))
  return(apply(x, 2, get_quantiles))
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
#' @return A matrix with summary statistics (2.5%, 50%, mean, 97.5%) for each time point.
infs.fct.qt <- function(infs, age.st, age.end, nt, qt.flag = FALSE) {
  ## infs is a matrix with quarterly log-infections for each posterior sample in each row
  ## qt.flag=TRUE aggregates quarterly infections and makes them yearly
  a1 <- age.st - 14
  a2 <- age.end - 14
  if (!qt.flag) {
    infs <- infs[, rep(1:ncol(infs), each = 4)] / 4 ## making yearly infs constant in quarters
    nt <- nt * 4 ## nt is quarter also for yearly model
  }
  call <- "rowSums(matrix(x, nr=nt)[,a1:a2])"
  x <- t(apply(exp(infs), 1, function(x) eval(parse(text = call))))
  return(apply(x, 2, get_quantiles))
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
#' @return A matrix with summary statistics (2.5%, 50%, mean, 97.5%) for each time point.
inc.fct <- function(infs, age.cl, last.yr, nt, qt.flag = FALSE, totals) {
  ## infs is a matrix with quarterly log-infections for each posterior sample in each row
  ## qt.flag=TRUE aggregates quarterly infections and makes them yearly
  if (!age.cl %in% c("Tot", "15-34", "35-44", "45+")) {
    stop("age.cl can only be 'Tot','15-34', '35-44','45+' ")
  }
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

  if (!qt.flag) {
    ind <- (nt - last.yr + 1):nt
    call <- "rowSums(matrix(x, nr=nt)[ind,a1:a2]) / den "
  } else {
    qt.indx <- rep(1:last.yr, each = 4)
    ind <- (nt - 4 * last.yr + 1):nt
    call <- "tapply(rowSums(matrix(x, nr=nt)[ind,a1:a2]), qt.indx, sum) / den"
  }

  x <- t(apply(exp(infs), 1, function(x) eval(parse(text = call))))
  return(apply(x, 2, get_quantiles))
}

#' Function to identify the peak in infection in the last x years
#'
#' @param infs A matrix with quarterly log-infections for each posterior sample in each row.
#' @param age.st Starting age for peak identification.
#' @param age.end Ending age for peak identification.
#' @param nt Total number of time points (quarters).
#' @param last.yr Number of last years to consider for peak identification.
#' @param qt.flag Logical flag indicating whether to aggregate quarterly infections into yearly.
#' @return A numeric vector indicating the peak position for each posterior sample.
peak.fct <- function(infs, age.st, age.end, nt, last.yr = 10, qt.flag = FALSE) {
  ## infs is a matrix with quarterly log-infections for each posterior sample in each row
  ## qt.flag=TRUE aggregates quarterly infections and makes them yearly
  if (last.yr > (nt / 4)) {
    stop("Specified last year is too large")
  }
  a1 <- age.st - 14
  a2 <- age.end - 14
  if (!qt.flag) {
    ind <- (nt - last.yr + 1):nt
    call <- "rowSums(matrix(x, nr=nt)[ind,a1:a2])"
  } else {
    qt.indx <- rep(1:last.yr, each = 4)
    ind <- (nt - 4 * last.yr + 1):nt
    call <- "tapply(rowSums(matrix(x, nr=nt)[ind,a1:a2]), qt.indx, sum)"
  }
  x <- t(apply(exp(infs), 1, function(x) eval(parse(text = call))))
  ### yields maximum in the last years
  return(apply(x, 1, which.max))
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
#' @return A character vector indicating the trend for each posterior sample.
trend.fct <- function(infs, age.st, age.end, nt, last.yr = 3, qt.flag = FALSE) {
  if (last.yr > (nt / 4)) {
    stop("Specified last year is too large")
  }
  a1 <- age.st - 14
  a2 <- age.end - 14
  if (!qt.flag) {
    ind <- (nt - last.yr + 1):nt
    call <- "rowSums(matrix(x, nr=nt)[ind,a1:a2])"
  } else {
    qt.indx <- rep(1:last.yr, each = 4)
    ind <- (nt - 4 * last.yr + 1):nt
    call <- "tapply(rowSums(matrix(x, nr=nt)[ind,a1:a2]), qt.indx, sum)"
  }
  x <- t(apply(exp(infs), 1, function(x) eval(parse(text = call))))
  return(apply(x, 1, shape.fct))
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
  ## infs is a matrix with log-infections for each posterior sample in each row
  if (
    !type %in%
      c("HIV", "AIDS", "CD4.st1", "CD4.st2", "CD4.st3", "CD4.st4", "CD4.st5")
  ) {
    stop(
      "Type can only be 'HIV','AIDS','CD4.st1','CD4.st2','CD4.st3','CD4.st4','CD4.st5'"
    )
  }
  a1 <- age.st - 14
  a2 <- age.end - 14
  if (qt.flag) {
    qt.indx <- rep(1:(nt / 4), each = 4)
  }
  call <- paste(type, "[,a1:a2,]", sep = "")
  y <- t(apply(eval(parse(text = call)), 3, rowSums))
  if (qt.flag) {
    y <- t(apply(y, 1, function(x) tapply(x, qt.indx, sum)))
  }
  return(apply(y, 2, get_quantiles))
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
  ## infs is a matrix with log-infections for each posterior sample in each row
  if (!state %in% c("all", 1:7)) {
    stop("State can only be 'all' or 1:7")
  }
  a1 <- age.st - 14
  a2 <- age.end - 14

  if (state == "all") {
    if (rita.flag) {
      prev <- prev.arr[,, 1, ] +
        prev.arr[,, 2, ] +
        prev.arr[,, 3, ] +
        prev.arr[,, 4, ] +
        prev.arr[,, 5, ] +
        prev.arr[,, 6, ] +
        prev.arr[,, 7, ]
    } else {
      prev <- prev.arr[,, 1, ] +
        prev.arr[,, 2, ] +
        prev.arr[,, 3, ] +
        prev.arr[,, 4, ]
    }
  } else {
    prev <- prev.arr[,, state, ]
  }

  x <- t(apply(prev[, a1:a2, ], 3, rowSums))
  return(apply(x, 2, get_quantiles))
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
      d.mat[1, , ] <- matrix(fit.mat[iter, ind_list$d1.ind], nrow = nt)
      d.mat[2, , ] <- matrix(fit.mat[iter, ind_list$d2.ind], nrow = nt)
      d.mat[3, , ] <- matrix(fit.mat[iter, ind_list$d3.ind], nrow = nt)
      d.mat[4, , ] <- matrix(fit.mat[iter, ind_list$d4.ind], nrow = nt)
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
  ### including under-reporting in expected AIDS, for each iter
  out <- matrix(NA, nrow = nrow(AIDS), ncol = ncol(AIDS))
  ### Under-reporting starts from year 2000
  ### which is year 5 (qt.flag=FALSE)
  ### or quarter 20 (qt.flag=FALSE)
  urep.start <- ifelse(qt.flag == TRUE, 20, 5)

  if (urep.flag) {
    out[1:urep.start, ] <- AIDS[1:urep.start, ] ### non under-reporting years
    out[(urep.start + 1):nt, ] <- urep * AIDS[(urep.start + 1):nt, ] ### under-reporting years
  } else {
    out <- AIDS
  }
  return(out)
}
