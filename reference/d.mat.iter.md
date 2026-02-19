# Function to get diagnosis matrix per iteration

Function to get diagnosis matrix per iteration

## Usage

``` r
d.mat.iter(
  fit_mat,
  iter,
  age.dx.flag,
  age.dx.spl.flag,
  rita.flag,
  nt,
  nage,
  ind_list
)
```

## Arguments

- fit_mat:

  A matrix of posterior samples from the fitted model.

- iter:

  The iteration number to extract data for.

- age.dx.flag:

  Logical flag indicating if age-specific diagnoses are considered.

- age.dx.spl.flag:

  Logical flag indicating if age-specific diagnoses with splines are
  considered.

- rita.flag:

  Logical flag indicating if RITA diagnoses are considered.

- nt:

  Total number of time points.

- nage:

  Number of age groups.

- ind_list:

  A list of indices for diagnosis data extraction.

## Value

A diagnosis matrix for the specified iteration.
