# Make fit indices

Creates a list of indices for extracting parameters from the fit matrix.

## Usage

``` r
make_fit_indices(fit_mat, config, nquar)
```

## Arguments

- fit_mat:

  Matrix of posterior samples from stan model.

- config:

  Configuration list used in stan model fitting.

- nquar:

  Number of quarters in the model.

## Value

List of indices for various parameters in the fit matrix.
