# Prepare spline basis for back-calculation model

This function prepares the spline basis and augments the stan_data list
with the design matrices for the back-calculation model.

## Usage

``` r
model_type_prep(
  stan_data,
  migration = NULL,
  rita = NULL,
  age = NULL,
  inf_model = 1,
  diag_model = 1,
  model_seed = NULL,
  basis = "ts",
  knots = 10,
  penalty = 2
)
```

## Arguments

- stan_data:

  A list of data formatted for Stan.

- migration:

  Logical indicating whether to include migration in the model.

- rita:

  Logical indicating whether to include RITA data in the model.

- age:

  Logical indicating whether to include age-specific data in the model.

- inf_model:

  Infection model type (1 = spline, 2 = random walk, 3 = GP).

- diag_model:

  Diagnosis model type (1 = spline, 2 = AR(1)).

- model_seed:

  Seed for model fitting.

- basis:

  Basis type for the spline (default is "ts" for thin plate regression
  spline).

- knots:

  Number of knots for the spline (default is 10).

- penalty:

  Order of the penalty for the spline (default is 2).

## Value

The augmented stan_data list with spline design matrices and prior
information.
