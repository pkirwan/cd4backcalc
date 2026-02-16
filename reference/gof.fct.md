# Function to get posterior functions of goodness of fit

This function generates goodness of fit statistics for HIV, AIDS, and
CD4 diagnoses by age group from their respective posterior samples (per
each posterior iteration).

## Usage

``` r
gof.fct(
  type = "HIV",
  age.st,
  age.end,
  HIV,
  AIDS,
  CD4.st1,
  CD4.st2,
  CD4.st3,
  CD4.st4,
  CD4.st5 = NULL,
  nt
)
```

## Arguments

- type:

  A character string indicating the type of diagnosis ("HIV", "AIDS",
  "CD4.st1", "CD4.st2", "CD4.st3", "CD4.st4").

- age.st:

  Starting age for goodness of fit calculation.

- age.end:

  Ending age for goodness of fit calculation.

- HIV:

  A 3D array of posterior samples for HIV diagnoses.

- AIDS:

  A 3D array of posterior samples for AIDS diagnoses.

- CD4.st1:

  A 3D array of posterior samples for CD4 stage 1 diagnoses.

- CD4.st2:

  A 3D array of posterior samples for CD4 stage 2 diagnoses.

- CD4.st3:

  A 3D array of posterior samples for CD4 stage 3 diagnoses.

- CD4.st4:

  A 3D array of posterior samples for CD4 stage 4 diagnoses.

- CD4.st5:

  A 3D array of posterior samples for CD4 stage 5 diagnoses (optional).

- nt:

  Total number of time points (quarters).

## Value

A matrix with summary statistics (2.5%, 50%, mean, 97.5%) for each time
point or age group.
