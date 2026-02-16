# Build posterior predictive quantities list

Constructs a list of posterior predictive quantities from the provided
distributions.

## Usage

``` r
build_ppq_list(stan_data, prev_dist, inf_dist, fit_slice, start_yr)
```

## Arguments

- stan_data:

  Data list used in stan model fitting.

- prev_dist:

  List of prevalence distributions.

- inf_dist:

  List of infection distributions.

- fit_slice:

  Matrix of posterior samples from stan model.

- start_yr:

  Starting year for posterior predictive quantities.

## Value

List of posterior predictive quantities.
