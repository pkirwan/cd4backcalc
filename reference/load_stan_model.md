# Choose Stan model based on configuration

This function chooses the appropriate Stan model based on the provided
configuration. The model is compiled at runtime using cmdstanr, with the
compiled binary cached so subsequent calls are fast.

## Usage

``` r
load_stan_model(config)
```

## Arguments

- config:

  Configuration object containing model settings.

## Value

The Stan model object (a CmdStanModel).
