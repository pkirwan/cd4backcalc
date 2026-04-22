# Generate D Matrix

Generates the D matrix for the epidemic model. CD4-only diagnosed states
definition: Diagnosed state 1: CD4 \> 500 at diagnosis Diagnosed state
2: 500 \< CD4 \< 350 at diagnosis Diagnosed state 3: 350 \< CD4 \< 200
at diagnosis Diagnosed state 4: CD4 \< 200 at diagnosis

## Usage

``` r
D_mat(t, d_probs, q_age, nlatst, ndiagst, nage = 52, rita = FALSE)
```

## Arguments

- t:

  Time index.

- d_probs:

  Diagnosis probabilities, a 3D array with dimensions
  `[nquar, nage, nstates]`.

- q_age:

  Age-specific progression probabilities.

- nlatst:

  Number of latent states.

- ndiagst:

  Number of diagnosed states.

- nage:

  Number of age groups.

- rita:

  Logical indicating whether RITA-specific logic should be applied.

## Value

D matrix.

## Details

RITA diagnosed states definition: Diagnosed state 1: RITA positive
diagnosis Diagnosed state 2: CD4 \> 500 at diagnosis Diagnosed state 3:
500 \< CD4 \< 350 at diagnosis Diagnosed state 4: 350 \< CD4 \< 200 at
diagnosis Diagnosed state 5: CD4 \< 200 at diagnosis

## Examples

``` r
if (FALSE) { # \dontrun{
q_age <- matrix(rep(c(0.2, 0.3, 0.4, 0.5), 52), nrow = 52, byrow = TRUE)
d_probs <- array(0.1, dim = c(10, 52, 4))
D_mat(t = 10, d_probs = d_probs, q_age, nlatst = 4, ndiagst = 5, nage = 52, rita = FALSE)
d_probs <- array(0.1, dim = c(10, 52, 5))
D_mat(t = 10, d_probs = d_probs, q_age, nlatst = 7, ndiagst = 6, nage = 52, rita = TRUE)
} # }
```
