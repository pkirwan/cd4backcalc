# Analysing results

This article covers extracting and visualising results from fitted
back-calculation models, including the full set of available quantities,
model comparison, and working with real data.

``` r

library(cd4backcalc)
library(patchwork)
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
| `"incidence_cob"` | UK-acquired incidence among UK-born and born-abroad people |
| `"all_migration"` | Estimated overall migration with HIV |
| `"undiag_migration"` | Estimated undiagnosed migration |
| `"diag_prob_mig"` | Estimated diagnosis probabilities for migrants by CD4 stratum |
| `"alpha_diag_mig"` | Born-abroad diagnosis offset odds ratio (vs UK-born), shared across CD4 strata, age and time |
| `"migration_prob"` | Estimated migration probabilities |
| `"ratio_abroad_uk"` | Ratio of UK-acquired incidence among people born abroad to that among UK-born people |
| `"transfers"` | Estimated transfers of care |
| `"transfers_detected"` | Estimated detected transfers of care |
| `"detect_prob"` | Estimated detection probabilities |
| `"migrant_hiv_prevalence"` | Estimated HIV-positive share of the external migration index |

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

sim_diags <- simulate_diagnoses(sim_type = "combo_3", migration = TRUE)

p1 <- plot_simulations(sim_diags, quantity = "incidence")
p2 <- plot_simulations(sim_diags, quantity = "undiag_prev")
p3 <- plot_simulations(sim_diags, quantity = "all_migration")
p4 <- plot_simulations(sim_diags, quantity = "undiag_migration")
p5 <- plot_simulations(sim_diags, quantity = "diag_prob")

(p1 + p2) / (p3 + p4) / p5
```

![Figure: Patchwork of simulated incidence, undiagnosed prevalence,
diagnosis probabilities, all migration, and undiagnosed
migration.](figures/analysing-results-sim-patchwork.png)

## Plotting model estimates

[`plot_estimates()`](https://pkirwan.github.io/cd4backcalc/reference/plot_estimates.md)
shows posterior median and 95% credible intervals. When simulated data
is present, the true values are shown for comparison:

``` r

hiv_list <- run_backcalc(sim_diags)

p1 <- plot_estimates(hiv_list, quantity = "incidence")
p2 <- plot_estimates(hiv_list, quantity = "undiag_prev")
p3 <- plot_estimates(hiv_list, quantity = "diag_prob")

(p1 + p2) / p3
```

![Figure: Patchwork of estimated incidence, undiagnosed prevalence, and
diagnosis
probabilities.](figures/analysing-results-estimates-patchwork.png)

By default, plotting and summary diagnostics compare fitted values with
the sampled realisation from the simulation. To compare against the
generative expectation from a particular simulation seed instead, set
`sim_comparison = "generative"`:

``` r

plot_estimates(
  hiv_list,
  quantity = "incidence",
  sim_comparison = "generative"
)
```

![Figure: Estimated incidence compared against the generative
expectation from the
simulation.](figures/analysing-results-est-incidence-generative.png)

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

    # A tibble: 6 × 7
       year `2.5%` `50%`  mean `97.5%` group   iteration
      <dbl>  <dbl> <dbl> <dbl>   <dbl> <chr>       <int>
    1  0        0     0     0        0 Model 1         1
    2  0.25   966  1073  1072.    1175 Model 1         1
    3  0.5   1006. 1101  1102.    1196 Model 1         1
    4  0.75  1052  1134. 1135.    1222 Model 1         1
    5  1     1082  1162  1161.    1244 Model 1         1
    6  1.25  1108  1186  1187.    1267 Model 1         1

## Extracting simulated quantities

Use
[`extract_simulated_quantity()`](https://pkirwan.github.io/cd4backcalc/reference/extract_simulated_quantity.md)
to extract the true simulated values for comparison:

``` r

sim_inc <- extract_simulated_quantity(sim_diags, quantity = "incidence")
head(sim_inc)
```

    # A tibble: 6 × 3
      value  year group            
      <dbl> <dbl> <chr>            
    1   978  0    Generative values
    2  1029  0.25 Generative values
    3  1039  0.5  Generative values
    4  1149  0.75 Generative values
    5  1195  1    Generative values
    6  1152  1.25 Generative values

Use
[`extract_real_data()`](https://pkirwan.github.io/cd4backcalc/reference/extract_real_data.md)
to extract the observed data used for fitting:

``` r

obs <- extract_real_data(hiv_list, quantity = "diagnoses")
head(obs)
```

    # A tibble: 6 × 3
      value  year group          
      <dbl> <dbl> <chr>          
    1   984  0    Observed values
    2   966  0.25 Observed values
    3  1008  0.5  Observed values
    4   929  0.75 Observed values
    5   990  1    Observed values
    6   939  1.25 Observed values

## Bias and MSE

When fitting to simulated data,
[`bias_plot()`](https://pkirwan.github.io/cd4backcalc/reference/bias_plot.md)
shows how much the model over- or under-estimates each quantity and
[`mse_plot()`](https://pkirwan.github.io/cd4backcalc/reference/mse_plot.md)
summarises mean squared error across one or more simulated fits:

``` r

p1 <- bias_plot(list(hiv_list), quantity = "incidence")
p2 <- mse_plot(list(hiv_list), quantity = "incidence")

p1 + p2
```

![Figure: Bias and MSE plots for
incidence.](figures/analysing-results-bias-mse-incidence.png)

## Comparing two models

[`bias_plot()`](https://pkirwan.github.io/cd4backcalc/reference/bias_plot.md)
and
[`mse_plot()`](https://pkirwan.github.io/cd4backcalc/reference/mse_plot.md)
accept a named list of models for side-by-side comparison:

``` r

sim_diags <- simulate_diagnoses(sim_type = "combo_3", rita = TRUE)
model_rita <- run_backcalc(sim_diags)

sim_cd4 <- sim_remodel(sim_diags, rita = FALSE)
model_cd4 <- run_backcalc(sim_cd4, rita = FALSE)

p1 <- plot_estimates(
  model_rita,
  hiv_list_2 = model_cd4,
  quantity = "incidence",
  model_1_name = "CD4 + RITA",
  model_2_name = "CD4 only"
)

p2 <- bias_plot(
  list("CD4 + RITA" = model_rita, "CD4 only" = model_cd4),
  quantity = "incidence"
)

p1 / p2
```

![Figure: Patchwork of incidence and bias comparisons for two
models.](figures/analysing-results-compare-two-patchwork.png)

## Comparing multiple models

Pass a list of model fits to compare more than two models
simultaneously:

``` r

sim_diags <- simulate_diagnoses(
  sim_type = "combo_3",
  sim_seed = 123
)

models <- list(
  "Spline incidence" = run_backcalc(
    sim_diags,
    inf_model = 1, diag_model = 1, model_seed = 123
  ),
  "Random-walk incidence" = run_backcalc(
    sim_diags,
    inf_model = 2, diag_model = 1, model_seed = 123
  ),
  "Gaussian-process incidence" = run_backcalc(
    sim_diags,
    inf_model = 3, diag_model = 1, model_seed = 123
  )
)

p1 <- plot_estimates(models, quantity = "incidence")
p2 <- bias_plot(models, quantity = "incidence")

p1 / p2
```

![Figure: Patchwork of incidence and bias comparisons for multiple
models.](figures/analysing-results-compare-smoothing-patchwork.png)

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

For chain-specific estimates, set `chain_results = TRUE` in
[`run_backcalc()`](https://pkirwan.github.io/cd4backcalc/reference/run_backcalc.md)
to obtain separate posterior summaries for each chain. A small number of
iterations are used in the example below to demonstrate poor convergence
behavior.

``` r

hiv_list <- run_backcalc(
  sim_diags,
  chain_results = TRUE,
  iter_warmup = 50,
  iter_sampling = 50
)
plot_estimates(hiv_list, quantity = "incidence", chain_results = TRUE)
```

![Figure: Per-chain incidence
estimates.](figures/analysing-results-chain-results-incidence.png)
