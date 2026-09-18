# Extract real data for plotting

Extract real data for plotting

## Usage

``` r
extract_real_data(hiv_list, quantity)
```

## Arguments

- hiv_list:

  A list of HIV model results containing `stan_data`.

- quantity:

  A string specifying the quantity to extract.

## Value

A data frame of observed values for plotting, or `NULL` when the
requested quantity is not present in the real-data inputs.

## Examples

``` r
# See ?run_backcalc for a full end-to-end workflow.
invisible(NULL)
```
