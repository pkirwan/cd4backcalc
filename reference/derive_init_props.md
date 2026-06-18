# Derive equilibrium initial state proportions

Computes quasi-steady-state latent compartment proportions by solving
the linear no-aging system `(I - Q) lat = h`, where `Q` is the
single-age-group transition matrix built from the supplied diagnosis
probabilities and progression rates. The result gives good starting
values for `init_props` or `init_props_abroad` when the undiagnosed
stock is assumed to be near equilibrium before the simulation window.

## Usage

``` r
derive_init_props(
  d_probs_t1 = c(0.055, 0.065, 0.13, 0.35),
  q = c(0.0976, 0.115, 0.116, 0.148),
  mu_lat_t1 = NULL,
  abroad = FALSE,
  scale = if (!abroad) 7 else NULL,
  rita = FALSE
)
```

## Arguments

- d_probs_t1:

  Numeric vector of diagnosis probabilities at t = 1. Length 4 for
  CD4-only (`rita = FALSE`); length 5 for RITA (`rita = TRUE`): element
  1 is the RITA (recent-infection) diagnosis probability, elements 2–5
  are the CD4-stratified probabilities.

- q:

  Numeric vector of length 4: quarterly CD4 progression probabilities,
  matching the default `q = c(0.0976, 0.115, 0.116, 0.148)`.

- mu_lat_t1:

  Numeric vector of relative inflow weights across latent states at t
  = 1. Required when `abroad = TRUE`. For RITA, the first three elements
  should be 0 (migrants enter long-standing states only). Will be
  normalised internally to sum to 1.

- abroad:

  Logical. `FALSE` (default): solve for the native group with all inflow
  entering state 1. `TRUE`: solve for the migrant stock with inflow
  distributed via `mu_lat_t1`.

- scale:

  Numeric or `NULL`. When non-`NULL`, the output is scaled so that its
  first non-zero element equals `scale`. Defaults to 7 for the native
  group so the result is on the same scale as the built-in `init_props`
  default. Defaults to `NULL` for the migrant group: the unscaled solve
  already gives values in the correct units for direct use as
  `init_props_abroad` (latent stock per unit of quarterly migration into
  latent states).

- rita:

  Logical. Use the RITA compartment structure (7 latent states).

## Value

Numeric vector of length `nlatst` (4 for CD4-only, 7 for RITA).

## Details

The function mirrors the compartment logic in
[`Q_mat`](https://pkirwan.github.io/cd4backcalc/reference/Q_mat.md): for
the native group (`abroad = FALSE`) all inflow enters state 1; for the
migrant stock (`abroad = TRUE`) inflow is distributed across latent
states proportionally to `mu_lat_t1`.
