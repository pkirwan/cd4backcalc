# Remodel the simulated diagnosis data for fitting a different model type.

This function takes the output from
[`simulate_diagnoses()`](https://pkirwan.github.io/cd4backcalc/reference/simulate_diagnoses.md)
and remodels the data to be suitable for fitting a different model type
(e.g., from RITA to CD4-only model).

## Usage

``` r
sim_remodel(sim_diags, rita = TRUE)
```

## Arguments

- sim_diags:

  A list output from
  [`simulate_diagnoses()`](https://pkirwan.github.io/cd4backcalc/reference/simulate_diagnoses.md).

- rita:

  Logical; if `TRUE`, remodels the data for fitting an RITA model; if
  `FALSE`, remodels the data for fitting a CD4-only model.

## Value

A list containing the remodeled stan data and configuration.

## Examples

``` r
#' # Simulate diagnoses data using RITA model
sim_diags <- simulate_diagnoses(nquar = 60, rita = TRUE, sim_seed = 123)
# Remodel data for fitting CD4-only model
remodeled_data <- sim_remodel(sim_diags, rita = FALSE)
```
