# Check that cmdstanr and CmdStan are available

`cmdstanr` is not on CRAN, so it is a suggested package rather than a
hard dependency. Model fitting needs both the R package and a CmdStan
installation, so this stops with installation instructions if either is
missing.

## Usage

``` r
check_cmdstanr()
```

## Value

`NULL`, invisibly.
