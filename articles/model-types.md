# Model types

The `cd4backcalc` package supports eight model configurations combining
CD4-only or CD4 + RITA data, with or without migration, and
age-independent or age-dependent structures. This article shows how to
simulate data and fit each variant.

## Supported models

| Model                  |             Age-independent              |              Age-dependent              |
|:-----------------------|:----------------------------------------:|:---------------------------------------:|
| CD4                    | `rita=FALSE, migration=FALSE, age=FALSE` | `rita=FALSE, migration=FALSE, age=TRUE` |
| CD4 + RITA             | `rita=TRUE, migration=FALSE, age=FALSE`  | `rita=TRUE, migration=FALSE, age=TRUE`  |
| CD4 + Migration        | `rita=FALSE, migration=TRUE, age=FALSE`  | `rita=FALSE, migration=TRUE, age=TRUE`  |
| CD4 + RITA + Migration |  `rita=TRUE, migration=TRUE, age=FALSE`  |  `rita=TRUE, migration=TRUE, age=TRUE`  |

The `rita`, `migration`, and `age` flags are passed to both
[`simulate_diagnoses()`](https://pkirwan.github.io/cd4backcalc/reference/simulate_diagnoses.md)
and
[`run_backcalc()`](https://pkirwan.github.io/cd4backcalc/reference/run_backcalc.md).
When fitting to real data these flags are passed to
[`run_backcalc()`](https://pkirwan.github.io/cd4backcalc/reference/run_backcalc.md)
only.

## CD4-only model

The CD4-only, age-independent model uses CD4 count at diagnosis and
known HIV progression probabilities to estimate incidence and
undiagnosed prevalence, see Birrell et al. (2012).

![Figure: CD4-staged back-calculation
model.](figures/backcalc-model-1.png)

Figure: CD4-staged back-calculation model. The model compartments
represent undiagnosed HIV (round boxes) and diagnosed HIV or AIDS
(square boxes); the $h_{j}$ arrow represents HIV acquired at each time
point ($j$); HIV progression probabilities ($q_{i}$) vary by CD4
stratum; HIV diagnosis probabilities ($d_{i,j}$) vary by CD4 stratum and
time.

This model can be fitted with:

``` r
library(cd4backcalc)

sim_cd4 <- simulate_diagnoses(h_pattern = "increasing")
fit_cd4 <- run_backcalc(sim_cd4)
plot_estimates(fit_cd4, quantity = "incidence")
```

## Age-dependent models

Age-dependent models stratify HIV incidence and prevalence by age group,
see Brizzi et al. (2019). These models take substantially longer to fit
due to the increased number of parameters. Checkpointing is available
for fitting age-dependent models on resource-constrained systems — see
[Checkpointing](https://pkirwan.github.io/cd4backcalc/articles/checkpointing.md).

![Figure: Age-dependent CD4 back-calculation
model.](figures/agedep-model-1.png)

Figure: Age-dependent CD4 back-calculation model with $k$ representing
the current age-category and $k_{0}$ representing the age category at
HIV acquisition.

``` r
sim_age <- simulate_diagnoses(age = TRUE, h_pattern = "increasing")
fit_age <- run_backcalc(sim_age)

# age-specific outputs
plot_estimates(fit_age, quantity = "incidence_age")
plot_estimates(fit_age, quantity = "undiag_prev_age")
plot_estimates(fit_age, quantity = "diagnoses_age")

# aggregate estimates are also available
plot_estimates(fit_age, quantity = "incidence")
```

## Adding RITA evidence

Evidence from RITA (Recent Infection Testing Algorithm) testing can
provide additional information on recent HIV acquisition, improving the
precision of incidence and undiagnosed prevalence estimates,
particularly for recent time periods and in populations with a high
proportion of recent diagnoses. Both the age-independent and
age-dependent models can be extended to include RITA evidence, see
Kirwan et al. (2026).

![Figure: Dual biomarker back-calculation model with recent incidence
assay and non-recent incidence assay CD4-staged
states.](figures/rita-model-1.png)

Figure: Dual biomarker back-calculation model with recent incidence
assay and non-recent incidence assay CD4-staged states. A/w = Acquired
within, mo. = months.

``` r
sim_rita <- simulate_diagnoses(rita = TRUE, h_pattern = "increasing")
fit_rita <- run_backcalc(sim_rita)
plot_estimates(fit_rita, quantity = "incidence")
```

## Migration-adjusted models

The migration-adjusted model uses information on country of birth to
distinguish between HIV acquired in the UK and abroad, and estimate
trends in migration. Both the age-independent and age-dependent models
can be extended to include migration transitions, as well as models with
RITA evidence, see Kirwan et al. (2026).

![Figure: Migration-adjusted CD4-staged back-calculation
model.](figures/migration-model-1.png)

Migration-adjusted CD4-staged back-calculation model. Solid lines
indicate HIV progression and transition from latent to diagnosed states;
dashed line is HIV acquisition among individuals born in the UK, dotted
lines are HIV acquisition or migration arrivals among individuals born
abroad.

``` r
sim_mig <- simulate_diagnoses(migration = TRUE, h_pattern = "increasing")
fit_mig <- run_backcalc(sim_mig)

# migration-specific outputs
plot_estimates(fit_mig, quantity = "undiag_migration")
plot_estimates(fit_mig, quantity = "migration_prob")
plot_estimates(fit_mig, quantity = "ratio_abroad_uk")
plot_estimates(fit_mig, quantity = "transfers")
plot_estimates(fit_mig, quantity = "detect_prob")
```

## Remodelling data

If data was simulated with RITA evidence but you want to fit a CD4-only
model (e.g. for comparison), use
[`sim_remodel()`](https://pkirwan.github.io/cd4backcalc/reference/sim_remodel.md)
to restructure the data:

``` r
# simulate with RITA
sim_rita <- simulate_diagnoses(rita = TRUE)

# remodel to CD4-only format
sim_cd4 <- sim_remodel(sim_rita, rita = FALSE)

# fit as CD4-only model
fit_cd4 <- run_backcalc(sim_cd4, rita = FALSE)
```

## Simulation types

For reproducible simulation studies, pre-defined simulation types
combine multiple parameter patterns. These are specified via the
`sim_type` argument:

``` r
# combo_3: varying acquisition + step-change migration + varying diagnosis +
#          increasing proportion UK + varying migration probabilities +
#          non-uniform age distribution
sim_c3 <- simulate_diagnoses(sim_type = "combo_3", migration = TRUE)

# combo_4: alternative incidence and migration patterns
sim_c4 <- simulate_diagnoses(sim_type = "combo_4", migration = TRUE)
```

Available simulation types include: `"constant"`, `"h_increasing"`,
`"h_varying"`, `"o_increasing"`, `"o_step_change"`, `"d_varying"`,
`"p_increasing"`, `"m_differing"`, `"m_varying"`, and `"combo_1"`
through `"combo_4"`.

## References

- Birrell PJ, Chadborn TR, Gill ON, Delpech VC, et al. (2012).
  Estimating trends in incidence, time-to-diagnosis and undiagnosed
  prevalence using a CD4-based Bayesian back-calculation. *Stat. Commun.
  Infect. Dis.* 4(1). doi:
  [10.1515/1948-4690.1055](https://doi.org/10.1515/1948-4690.1055).
- Brizzi F, Birrell PJ, Plummer MT, Kirwan PD, et al. (2019). Extending
  Bayesian back-calculation to estimate age and time specific HIV
  incidence. *Lifetime Data Anal.* 25(4), pp.757-780. doi:
  [10.1007/s10985-019-09465-1](https://doi.org/10.1007/s10985-019-09465-1).
- Kirwan PD, Presanis A, Birrell PJ, et al. (2026). Extending a Bayesian
  back-calculation model for HIV incidence to include biomarkers of
  recent acquisition. (in press).
- Kirwan PD, Presanis A, Birrell PJ, et al. (2026). HIV incidence among
  gay and bisexual men in England, Wales, and Northern Ireland:
  estimates from a migration-adjusted CD4-staged HIV back-calculation
  model. (in press).
