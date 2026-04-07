# Generate Q Matrix

Generates the Q matrix for the epidemic model. CD4-only latent states
definition: Latent state 1: CD4 \> 500 Latent state 2: 500 \< CD4 \< 350
Latent state 3: 200 \< CD4 \< 350 Latent state 4: CD4 \< 200

## Usage

``` r
Q_mat(t, d_probs, q_age, nlatst, nage = 52, rita = FALSE)
```

## Arguments

- t:

  Time index.

- d_probs:

  Diagnosis probabilities.

- q_age:

  Age-specific progression probabilities.

- nlatst:

  Number of latent states.

- nage:

  Number of age groups.

- rita:

  Logical indicating whether RITA-specific logic should be applied.

## Value

Q matrix.

## Details

RITA latent states definition: Latent state 1: CD4 \> 500, HIV acquired
in the last three months Latent state 2: CD4 \> 500, HIV acquired
between 3 and 6 months ago Latent state 3: 500 \< CD4 \< 350, HIV
acquired between 3 and 6 months ago Latent state 4: CD4 \> 500, HIV
acquired more than 6 months ago Latent state 5: 500 \< CD4 \< 350, HIV
acquired more than 6 months ago Latent state 6: 350 \< CD4 \< 200, HIV
acquired more than 6 months ago Latent state 7: CD4 \< 200, HIV acquired
more than 6 months ago

## Examples

``` r
if (FALSE) { # \dontrun{
q_age <- matrix(rep(c(0.2, 0.3, 0.4, 0.5), 52), nrow = 52, byrow = TRUE)
d_probs <- matrix(rep(0.1, 40), nrow = 10)
Q_mat(t = 10, d_probs, q_age, nlatst = 4, nage = 52, rita = FALSE)
d_probs <- matrix(rep(0.1, 50), nrow = 10)
Q_mat(t = 10, d_probs, q_age, nlatst = 7, nage = 52, rita = TRUE)
} # }
```
