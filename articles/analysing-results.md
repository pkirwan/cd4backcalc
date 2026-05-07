# Analysing results

This article covers extracting and visualising results from fitted
back-calculation models, including the full set of available quantities,
model comparison, and working with real data.

``` r

library(cd4backcalc)
```

## Available quantities

After fitting a model with
[`run_backcalc()`](https://pkirwan.github.io/cd4backcalc/reference/run_backcalc.md),
several quantities can be plotted or extracted. The available quantities
depend on the model type:

**All models:**

| Quantity        | Description                                      |
|-----------------|--------------------------------------------------|
| `"incidence"`   | Estimated HIV incidence over time                |
| `"undiag_prev"` | Estimated undiagnosed HIV prevalence             |
| `"diagnoses"`   | Expected vs observed diagnoses                   |
| `"diag_prob"`   | Estimated diagnosis probabilities by CD4 stratum |

**Migration models** (additional):

| Quantity | Description |
|----|----|
| `"all_migration"` | Estimated overall migration with HIV |
| `"undiag_migration"` | Estimated undiagnosed migration |
| `"diag_prob_mig"` | Estimated diagnosis probabilities for migrants by CD4 stratum |
| `"migration_prob"` | Estimated migration probabilities |
| `"ratio_abroad_uk"` | Proportion of infections acquired in the UK |
| `"transfers"` | Estimated transfers of care |
| `"detect_prob"` | Estimated detection probabilities |

**Age-dependent models** (additional):

| Quantity            | Description                                    |
|---------------------|------------------------------------------------|
| `"incidence_age"`   | Incidence stratified by age group              |
| `"undiag_prev_age"` | Undiagnosed prevalence by age group            |
| `"diagnoses_age"`   | Diagnoses by age group                         |
| `"diag_prob_age"`   | Diagnosis probabilities by CD4 stratum and age |

**Age-dependent migration models** (additional):

| Quantity | Description |
|----|----|
| `"undiag_migration_age"` | Undiagnosed migration by age group |
| `"diag_prob_mig_age"` | Migrant diagnosis probabilities by CD4 stratum and age |

## Plotting simulated data

Before fitting, simulated quantities can be inspected with
[`plot_simulations()`](https://pkirwan.github.io/cd4backcalc/reference/plot_simulations.md):

``` r

sim_diags <- simulate_diagnoses(migration = TRUE, h_pattern = "increasing")

plot_simulations(sim_diags, quantity = "incidence")
plot_simulations(sim_diags, quantity = "undiag_prev")
plot_simulations(sim_diags, quantity = "diag_prob")
plot_simulations(sim_diags, quantity = "all_migration")
plot_simulations(sim_diags, quantity = "undiag_migration")
```

## Plotting model estimates

[`plot_estimates()`](https://pkirwan.github.io/cd4backcalc/reference/plot_estimates.md)
shows posterior median and 95% credible intervals. When simulated data
is present, the true values are shown for comparison:

``` r

hiv_list <- run_backcalc(sim_diags)
plot_estimates(hiv_list, quantity = "incidence")
```

By default, plotting and summary diagnostics compare fitted values with
the generative expectation from the simulation. To compare against the
sampled realisation from a particular simulation seed instead, set
`sim_comparison = "sampled"`:

``` r

plot_estimates(
  hiv_list,
  quantity = "incidence",
  sim_comparison = "sampled"
)
```

The same argument is available in
[`bias_plot()`](https://pkirwan.github.io/cd4backcalc/reference/bias_plot.md)
and
[`mse_plot()`](https://pkirwan.github.io/cd4backcalc/reference/mse_plot.md).

## Extracting estimates as data

Use
[`get_estimates()`](https://pkirwan.github.io/cd4backcalc/reference/get_estimates.md)
to obtain the posterior summary as a tibble for custom analysis:

``` r

inc_df <- get_estimates(hiv_list, quantity = "incidence")
head(inc_df)
```

## Extracting simulated quantities

Use
[`extract_simulated_quantity()`](https://pkirwan.github.io/cd4backcalc/reference/extract_simulated_quantity.md)
to extract the true simulated values for comparison:

``` r

sim_inc <- extract_simulated_quantity(sim_diags, quantity = "incidence")
head(sim_inc)
```

Use
[`extract_real_data()`](https://pkirwan.github.io/cd4backcalc/reference/extract_real_data.md)
to extract the observed data used for fitting:

``` r

obs <- extract_real_data(hiv_list, quantity = "diagnoses")
head(obs)
```

## Bias and MSE

When fitting to simulated data,
[`bias_plot()`](https://pkirwan.github.io/cd4backcalc/reference/bias_plot.md)
shows how much the model over- or under-estimates each quantity:

``` r

bias_plot(hiv_list, quantity = "incidence")
```

When fitting to simulated data,
[`mse_plot()`](https://pkirwan.github.io/cd4backcalc/reference/mse_plot.md)
summarises mean squared error across one or more simulated fits:

``` r

mse_plot(hiv_list, quantity = "incidence")
```

## Comparing two models

[`plot_estimates()`](https://pkirwan.github.io/cd4backcalc/reference/plot_estimates.md),
[`bias_plot()`](https://pkirwan.github.io/cd4backcalc/reference/bias_plot.md),
and
[`mse_plot()`](https://pkirwan.github.io/cd4backcalc/reference/mse_plot.md)
all accept a second model for side-by-side comparison:

``` r

sim_diags <- simulate_diagnoses(rita = TRUE, h_pattern = "increasing")
model_rita <- run_backcalc(sim_diags)

sim_cd4 <- sim_remodel(sim_diags, rita = FALSE)
model_cd4 <- run_backcalc(sim_cd4, rita = FALSE)

plot_estimates(
  model_rita,
  hiv_list_2 = model_cd4,
  quantity = "incidence",
  model_1_name = "CD4 + RITA",
  model_2_name = "CD4 only"
)

bias_plot(
  model_rita,
  hiv_list_2 = model_cd4,
  quantity = "incidence",
  model_1_name = "CD4 + RITA",
  model_2_name = "CD4 only"
)
```

## Comparing multiple models

Pass a list of model fits to compare more than two models
simultaneously:

``` r

sim_diags <- simulate_diagnoses(h_pattern = "increasing")
models <- list()
for (i in 1:3) {
  models[[i]] <- run_backcalc(sim_diags, model_seed = 100 + i)
}

plot_estimates(models, quantity = "incidence")
mse_plot(models, quantity = "incidence")
```

## Working with real data

When fitting to real (non-simulated) data, plotting and extraction
functions still work — the simulated truth overlay is simply omitted.
The `stan_data` object should contain the diagnosis arrays required for
the chosen model, with matching `migration`, `rita`, and `age` flags
supplied to
[`run_backcalc()`](https://pkirwan.github.io/cd4backcalc/reference/run_backcalc.md).

``` r

hiv_list <- run_backcalc(
  real_stan_data,
  migration = TRUE,
  rita = TRUE,
  age = FALSE
)

plot_estimates(hiv_list, quantity = "incidence")
plot_estimates(hiv_list, quantity = "diagnoses")
extract_real_data(hiv_list, quantity = "diagnoses")
```

## Per-chain results

For detailed convergence analysis, set `chain_results = TRUE` in
[`run_backcalc()`](https://pkirwan.github.io/cd4backcalc/reference/run_backcalc.md)
to obtain separate posterior summaries for each chain:

``` r

hiv_list <- run_backcalc(sim_diags, chain_results = TRUE)
plot_estimates(hiv_list, quantity = "incidence", chain_results = TRUE)
```
