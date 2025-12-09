#' Logit function
#' Computes the logit transformation of a numeric vector.
#' @param x Numeric vector input.
#' @return Logit transformed values.
logit <- function(x) log(x / (1 - x))

#' Inverse logit function
#' Computes the inverse logit transformation of a numeric vector.
#' @param x Numeric vector input.
#' @return Inverse logit transformed values.
invlogit <- function(x) exp(x) / (1 + exp(x))

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
#' @return Rounded numeric matrix with preserved total.
round_matrix_preserve <- function(mat) {
  # Convert matrix to vector, round, and convert back
  dims <- dim(mat)
  vec <- as.vector(mat)
  rounded <- round_preserve(vec)
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

#' Round specified fields in an epidemic object
#' 
#' @param epi_obj Epidemic object containing fields to round.
#' @param fields Character vector of field names to round.
#' @return Epidemic object with specified fields rounded.
round_fields <- function(epi_obj, fields = c("lat", "diags", "lat_age", "diags_age")) {
  for (field in fields) {
    if (!is.null(epi_obj[[field]])) {
      # Check if array (3+ dimensions) or matrix (2 dimensions)
      if (length(dim(epi_obj[[field]])) > 2) {
        epi_obj[[field]] <- round_array_preserve(epi_obj[[field]])
      } else {
        epi_obj[[field]] <- round_matrix_preserve(epi_obj[[field]])
      }
    }
  }
  return(epi_obj)
}