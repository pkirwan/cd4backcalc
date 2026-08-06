# Simulate epidemic over time

Simulates the epidemic over a specified number of time quarters using
the provided Q and D matrices, initial latent state, and new infections.

## Usage

``` r
epi_mat(Q_l, D_l, h, lat_init, config)
```

## Arguments

- Q_l:

  List of Q matrices for each time quarter.

- D_l:

  List of D matrices for each time quarter.

- h:

  Matrix of new infections over time.

- lat_init:

  Initial latent state vector.

- config:

  Configuration list containing simulation parameters.

## Value

A list containing matrices of latent and diagnosed states over time, as
well as detailed quarterly records.

## Examples

``` r
# \donttest{
config <- define_sim_config(nquar = 10)
#> Error in define_sim_config(nquar = 10): could not find function "define_sim_config"
q_age <- matrix(rep(c(0.0976, 0.115, 0.116, 0.148), 52), nrow = 52, byrow = TRUE)
d_probs <- matrix(rep(0.1, 40), nrow = 10)
Q_l <- lapply(1:10, function(x) {
   cd4backcalc:::Q_mat(x, d_probs, q_age, 4)
})
#> Error in d_probs[t, a, ]: incorrect number of dimensions
D_l <- lapply(1:10, function(x) {
   cd4backcalc:::D_mat(x, d_probs, q_age, 4, 5)
})
#> Error in d_probs[t, a, ]: incorrect number of dimensions
lat_init <- array(0, dim = c(52, 4))
lat_init[, 1] <- c(40, 30, 20, 10)
h <- array(matrix(rep(10, 20), nrow = 10), dim = c(52, 10, 4))
cd4backcalc:::epi_mat(Q_l, D_l, h, lat_init, config)
#> Error: object 'config' not found
# }
```
