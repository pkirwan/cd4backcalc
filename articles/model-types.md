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

## CD4-only model (basic)

The simplest model uses CD4 count at diagnosis as the only evidence.

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

## Adding RITA evidence

RITA (Recent Infection Testing Algorithm) provides additional
information on recent infections, improving incidence estimation:

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

## Adding migration

When migration data is available, the model can distinguish between
infections acquired in the UK and abroad, and estimate undiagnosed
migration.

![Figure: CD4-staged back-calculation model including migration
transitions.](figures/migration-model-1.png)

CD4-staged back-calculation model including migration transitions. Solid
lines indicate HIV progression and transition from latent to diagnosed
states; dashed line is HIV acquisition among individuals born in the UK,
dotted lines are HIV acquisition or migration arrivals among individuals
born abroad.

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

## Age-dependent models

Age-dependent models stratify incidence and prevalence by age group.
These models are more complex and take substantially longer to fit:

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

For age-dependent models on resource-constrained systems, checkpointing
is recommended — see
[Checkpointing](https://pkirwan.github.io/cd4backcalc/articles/checkpointing.md).

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
# combo_3: varying infections + step-change migration + varying diagnosis +
#          increasing proportion UK + varying migration probabilities
sim_c3 <- simulate_diagnoses(sim_type = "combo_3", migration = TRUE)

# combo_4: alternative incidence and migration patterns
sim_c4 <- simulate_diagnoses(sim_type = "combo_4", migration = TRUE)
```

Available simulation types include: `"constant"`, `"h_increasing"`,
`"h_varying"`, `"o_increasing"`, `"o_step_change"`, `"d_varying"`,
`"p_increasing"`, `"m_differing"`, `"m_varying"`, `"combo_1"` through
`"combo_4"`, and `"combo_3_no_migration"`.
