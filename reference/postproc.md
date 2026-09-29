# Post-processing model fit

This function decides which post-processing routine needs to be run and
passes the relevant data and arguments.

## Usage

``` r
postproc(fit, stan_data, start_yr = 0, nworkers = 1, chain_results = FALSE)
```

## Arguments

- fit:

  Fitted stan model object and stan_data list.

- stan_data:

  Data list used in stan model fitting.

- start_yr:

  Starting year for posterior predictive quantities.

- nworkers:

  Number of parallel workers to use.

- chain_results:

  Logical indicating whether to return results by chain.

## Value

List of posterior predictive quantities.

## Examples

``` r
# \donttest{
sim_diags <- simulate_diagnoses(sim_seed = 123)
stan_data <- cd4backcalc:::model_type_prep(sim_diags, model_seed = 123)
fit <- cd4backcalc:::fit_stan_model(
   stan_data,
   iter_warmup = 10,
   iter_sampling = 10,
   chains = 1,
   refresh = 0,
   messages = FALSE
)
#> CmdStan model cache dir: /tmp/Rtmp4ImGLk/cd4backcalc-cmdstan
#> R tempdir(): /tmp/Rtmp4ImGLk
#> Error: CmdStan path has not been set yet. See ?set_cmdstan_path.
hiv_list <- cd4backcalc:::postproc(fit, stan_data, start_yr = 0, nworkers = 1)
#> Error: object 'fit' not found
# }
```
