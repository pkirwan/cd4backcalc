# Function to convert a matrix to a tibble

Function to convert a matrix to a tibble

## Usage

``` r
make_tibble(quantiles, nt, start_yr, qt.flag = TRUE)
```

## Arguments

- quantiles:

  A matrix of quantiles.

- nt:

  Total number of time points.

- start_yr:

  Starting year for the time series.

- qt.flag:

  Logical flag indicating if this is quarterly data.

## Value

A tibble with quartile, year, and value columns.
