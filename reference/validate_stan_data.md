# Validate Stan data

Validates the consistency of the Stan data with the simulated epidemic
data.

## Usage

``` r
validate_stan_data(stan_data, sim_data, config)
```

## Arguments

- stan_data:

  List of data prepared for Stan model fitting.

- sim_data:

  List of simulated epidemic data.

- config:

  Configuration list defining simulation parameters.

## Value

NULL. Throws an error if validation checks fail.
