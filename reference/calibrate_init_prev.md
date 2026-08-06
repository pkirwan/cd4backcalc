# Calibrate initial prevalence from fitted draws

Back-calculates improved initial prevalence values from quarter-1
observed diagnoses and the fitted quarter-1 diagnosis probabilities. The
helper works with age-dependent and age-independent models, with or
without migration, and for both CD4-only and CD4-plus-RITA variants.

## Usage

``` r
calibrate_init_prev(
  object,
  stan_data = NULL,
  summary = c("median", "mean"),
  use_hiv_total = TRUE,
  smooth_age = TRUE,
  smooth_span = 0.5,
  recent_split = c("current", "equal"),
  min_prob = 1e-06,
  init_floor = 1e-08,
  return_draws = FALSE
)
```

## Arguments

- object:

  Either a list returned by
  [`run_backcalc()`](https://pkirwan.github.io/cd4backcalc/reference/run_backcalc.md)
  containing `fit` and `stan_data`, or a fitted model object with a
  `draws(format = "matrix")` method.

- stan_data:

  Optional Stan data list. Supply this when `object` is a fit rather
  than the full list returned by
  [`run_backcalc()`](https://pkirwan.github.io/cd4backcalc/reference/run_backcalc.md).

- summary:

  Posterior summary to return. One of `"median"` or `"mean"`.

- use_hiv_total:

  Logical; if `TRUE`, quarter-1 HIV totals are allocated across CD4/RITA
  categories using the observed category proportions. This is useful
  when CD4 completeness is below 100%.

- smooth_age:

  Logical; if `TRUE`, age-specific summaries are smoothed over age with
  loess after summarising across draws.

- smooth_span:

  Span passed to [`stats::loess()`](https://rdrr.io/r/stats/loess.html)
  when `smooth_age = TRUE`.

- recent_split:

  How to split the unidentified recent latent pool in RITA models.
  `"current"` preserves the relative weights in the current `init_prev`
  input; `"equal"` splits it equally across states 1:3.

- min_prob:

  Lower bound used when inverting fitted diagnosis probabilities.

- init_floor:

  Small positive floor applied to calibrated initial prevalence values
  to avoid non-finite, zero, or tiny values that can cause Stan
  initialisation failures.

- return_draws:

  Logical; if `TRUE`, include draw-level calibrated values in the
  output.

## Value

A list containing calibrated initial prevalence values, an updated
`stan_data` list ready for refitting, and optional draw-level arrays.

## Details

The returned values are posterior summaries of baseline latent
prevalence that would reproduce the first-quarter diagnoses under the
fitted diagnosis model. For RITA models, the first three latent "recent"
states are not separately identified from quarter-1 data; by default
they are split using the relative weights in the currently supplied
initial prevalence.

## Examples

``` r
# See ?run_backcalc for a full end-to-end workflow.
invisible(NULL)
```
