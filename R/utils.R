#' Logit function
#' Computes the logit transformation of a numeric vector.
#' @param x Numeric vector input.
#' @return Logit transformed values.
logit <- function(x) {
  x <- pmin(pmax(x, .Machine$double.eps), 1 - .Machine$double.eps)
  stats::qlogis(x)
}

#' Inverse logit function
#' Computes the inverse logit transformation of a numeric vector.
#' @param x Numeric vector input.
#' @return Inverse logit transformed values.
invlogit <- function(x) stats::plogis(x)

#' Controlled rounding function to preserve totals
#'
#' @param x Numeric vector to be rounded.
#' @return Rounded numeric vector with preserved total.
round_preserve <- function(x) {
  # Get decimal parts and their ranking
  dec <- x - floor(x)
  n_round_up <- round(sum(dec)) # Number of values to round up

  # Determine which values to round up (those with highest decimals)
  round_up <- rank(-dec, ties.method = "random") <= n_round_up

  # Apply rounding
  out <- floor(x)
  out[round_up] <- ceiling(x[round_up])

  return(out)
}

#' Matrix rounding function to preserve totals
#'
#' @param mat Numeric matrix to be rounded.
#' @param colwise Logical; if TRUE, round each column separately.
#' @param rowwise Logical; if TRUE, round each row separately.
#' @return Rounded numeric matrix with preserved total.
round_matrix_preserve <- function(mat, colwise = FALSE, rowwise = FALSE) {
  # Convert matrix to vector, round, and convert back
  dims <- dim(mat)
  vec <- as.vector(mat)
  if (colwise) {
    # Round each column separately
    rounded_vec <- unlist(
      lapply(
        split(vec, rep(1:dims[2], each = dims[1])),
        round_preserve
      )
    )
    rounded <- rounded_vec
  } else if (rowwise) {
    # Round each row separately
    rounded_vec <- unlist(
      lapply(
        split(vec, rep(1:dims[1], times = dims[2])),
        round_preserve
      )
    )
    return(matrix(rounded_vec, nrow = dims[1], ncol = dims[2], byrow = TRUE))
  } else {
    rounded <- round_preserve(vec)
  }
  matrix(rounded, nrow = dims[1], ncol = dims[2])
}

#' Array rounding function to preserve totals
#'
#' @param arr Numeric array to be rounded.
#' @return Rounded numeric array with preserved total.
round_array_preserve <- function(arr) {
  # Convert array to vector, round, and convert back
  dims <- dim(arr)
  vec <- as.vector(arr)
  rounded <- round_preserve(vec)
  array(rounded, dim = dims)
}

#' Round age-stratified arrays to match pre-rounded quarter/state totals
#'
#' @param arr Numeric array with age in the first dimension.
#' @param totals Numeric array giving the target sums for each combination of
#'   the remaining dimensions.
#' @param stochastic Logical; if TRUE, allocate each total across ages with a
#'   multinomial draw instead of largest-remainder rounding.
#' @return Rounded array whose sums over age match totals exactly.
round_age_array_to_totals <- function(arr, totals, stochastic = FALSE) {
  dims <- dim(arr)

  allocate <- function(slice, target) {
    slice_total <- sum(slice)
    if (target <= 0 || slice_total <= 0) {
      return(rep(0, length(slice)))
    }
    if (stochastic) {
      return(as.vector(stats::rmultinom(1, target, slice / slice_total)))
    }
    round_preserve(target * slice / slice_total)
  }

  if (length(dims) < 2) {
    return(round_preserve(arr))
  }

  if (length(dims) == 2) {
    rounded <- matrix(0, nrow = dims[1], ncol = dims[2])
    for (j in seq_len(dims[2])) {
      slice <- arr[, j]
      target <- totals[j]
      rounded[, j] <- allocate(slice, target)
    }
    return(rounded)
  }

  if (length(dims) != 3) {
    stop("round_age_array_to_totals expects a 2D or 3D array with age in the first dimension")
  }

  rounded <- array(0, dim = dims)
  for (j in seq_len(dims[2])) {
    for (k in seq_len(dims[3])) {
      slice <- arr[, j, k]
      target <- totals[j, k]
      rounded[, j, k] <- allocate(slice, target)
    }
  }

  rounded
}

#' Round specified fields in an epidemic object
#'
#' @param epi_obj Epidemic object containing fields to round.
#' @param fields Character vector of field names to round.
#' @param stochastic_diags_age Logical; if TRUE, allocate diagnosis totals
#'   across ages with multinomial draws.
#' @return Epidemic object with specified fields rounded.
round_fields <- function(
  epi_obj,
  fields = c("lat", "diags", "lat_age", "diags_age"),
  stochastic_diags_age = FALSE
) {
  if ("lat_age" %in% fields && !is.null(epi_obj$lat_age)) {
    lat_totals <- if ("lat" %in% fields && !is.null(epi_obj$lat)) {
      round_matrix_preserve(epi_obj$lat)
    } else {
      round_matrix_preserve(apply(epi_obj$lat_age, c(2, 3), sum))
    }
    epi_obj$lat_age <- round_age_array_to_totals(epi_obj$lat_age, lat_totals)
    if ("lat" %in% fields) {
      epi_obj$lat <- lat_totals
    }
  } else if ("lat" %in% fields && !is.null(epi_obj$lat)) {
    epi_obj$lat <- round_matrix_preserve(epi_obj$lat)
  }

  if ("diags_age" %in% fields && !is.null(epi_obj$diags_age)) {
    diags_totals <- if ("diags" %in% fields && !is.null(epi_obj$diags)) {
      round_matrix_preserve(epi_obj$diags)
    } else {
      round_matrix_preserve(apply(epi_obj$diags_age, c(2, 3), sum))
    }
    epi_obj$diags_age <- round_age_array_to_totals(
      epi_obj$diags_age,
      diags_totals,
      stochastic = stochastic_diags_age
    )
    if ("diags" %in% fields) {
      epi_obj$diags <- diags_totals
    }
  } else if ("diags" %in% fields && !is.null(epi_obj$diags)) {
    epi_obj$diags <- round_matrix_preserve(epi_obj$diags)
  }

  remaining_fields <- setdiff(fields, c("lat", "diags", "lat_age", "diags_age"))
  for (field in remaining_fields) {
    if (!is.null(epi_obj[[field]])) {
      if (length(dim(epi_obj[[field]])) > 2) {
        epi_obj[[field]] <- round_array_preserve(epi_obj[[field]])
      } else {
        epi_obj[[field]] <- round_matrix_preserve(epi_obj[[field]])
      }
    }
  }
  return(epi_obj)
}
