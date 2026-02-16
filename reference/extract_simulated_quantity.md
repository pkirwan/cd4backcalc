# Extract simulated quantities based on the provided configuration.

Extract simulated quantities based on the provided configuration.

## Usage

``` r
extract_simulated_quantity(sim_diags, quantity)
```

## Arguments

- sim_diags:

  The simulated diagnosis data containing stan_data and sim_data.

- quantity:

  The quantity to extract.

## Value

A list of model estimates including undiagnosed prevalence, incidence,
migration, diagnosis probabilities, diagnoses, migration probabilities,
proportion of UK-acquired infections that are migrants, transfers of
care, and CD4 strata.
