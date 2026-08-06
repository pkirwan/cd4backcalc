# Define configuration for fitting back-calculation model

This function defines the configuration settings for fitting the CD4
back-calculation model.

## Usage

``` r
stan_model_config(
  migration = FALSE,
  rita = FALSE,
  age = FALSE,
  inf_model = 1,
  diag_model = 1,
  model_seed = NULL
)
```

## Arguments

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

## Value

A list containing configuration settings for model fitting.
