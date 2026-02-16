# Define simulation configuration

Defines the patterns for various parameters in the epidemic simulation.
The number of latent and diagnosis states will be inferred from the
simulation specification.

## Usage

``` r
define_sim_config(
  h_pattern = "constant",
  o_pattern = "constant",
  d_pattern = "constant",
  p_pattern = "constant",
  m_pattern = "constant",
  r_pattern = "constant",
  a_pattern = "constant",
  sim_type = NULL,
  sim_seed = NULL,
  dist = TRUE,
  migration = FALSE,
  rita = FALSE,
  age = FALSE,
  nquar = 120,
  nage = 52,
  q = c(0.0976, 0.115, 0.116, 0.148),
  target_inf = 90000,
  target_mig = 12000,
  init_props = c(7, 3.5, 1.75, 0.5),
  init_props_abroad = c(1.75, 2.25, 2.25, 1.25),
  abroad_prop = 0.2,
  cd4_completeness = 1,
  rita_completeness = 1
)
```

## Arguments

- h_pattern:

  Pattern type for infections ("constant", "increasing", "varying").

- o_pattern:

  Pattern type for omega ("constant", "increasing", "step_change",
  "none").

- d_pattern:

  Pattern type for diagnosis probabilities ("constant", "varying").

- p_pattern:

  Pattern type for proportion UK ("constant", "increasing").

- m_pattern:

  Pattern type for mu ("constant", "differing", "varying").

- r_pattern:

  Pattern type for RITA data ("constant", "increasing", "none").

- a_pattern:

  Pattern type for age-based data ("constant", "wedge",
  "varying_wedge").

- sim_type:

  If provided, a preset type of simulation to define patterns for.

- sim_seed:

  Random seed for simulation. If NULL, a random seed is generated.

- dist:

  Logical indicating whether to add stochastic variation.

- migration:

  Logical indicating whether to include migration effects.

- rita:

  Logical indicating whether to simulate RITA (Recent Infection Testing
  Algorithm) data.

- age:

  Logical indicating whether to simulate age-stratified data.

- nquar:

  Number of time quarters.

- nage:

  Number of age groups.

- q:

  Vector of progression probabilities.

- target_inf:

  Target total number of HIV acquisitions.

- target_mig:

  Target total number of migrations.

- init_props:

  Initial state proportions.

- init_props_abroad:

  Initial state proportions for infections acquired abroad.

- abroad_prop:

  Proportion of infections acquired abroad at the outset.

- cd4_completeness:

  CD4 completeness.

- rita_completeness:

  RITA completeness.

## Value

Configuration list.

## Examples

``` r
define_sim_config(h_pattern = "constant")
#> $h_pattern
#> [1] "constant"
#> 
#> $o_pattern
#> [1] "none"
#> 
#> $d_pattern
#> [1] "constant"
#> 
#> $p_pattern
#> [1] "constant"
#> 
#> $m_pattern
#> [1] "constant"
#> 
#> $r_pattern
#> [1] "none"
#> 
#> $a_pattern
#> [1] "constant"
#> 
#> $sim_seed
#> [1] 2137934383
#> 
#> $dist
#> [1] TRUE
#> 
#> $migration
#> [1] FALSE
#> 
#> $rita
#> [1] FALSE
#> 
#> $age
#> [1] FALSE
#> 
#> $nquar
#> [1] 120
#> 
#> $nage
#> [1] 52
#> 
#> $nlatst
#> [1] 4
#> 
#> $ndiagst
#> [1] 5
#> 
#> $q
#> [1] 0.0976 0.1150 0.1160 0.1480
#> 
#> $q_age
#>         [,1]  [,2]  [,3]  [,4]
#>  [1,] 0.0976 0.115 0.116 0.148
#>  [2,] 0.0976 0.115 0.116 0.148
#>  [3,] 0.0976 0.115 0.116 0.148
#>  [4,] 0.0976 0.115 0.116 0.148
#>  [5,] 0.0976 0.115 0.116 0.148
#>  [6,] 0.0976 0.115 0.116 0.148
#>  [7,] 0.0976 0.115 0.116 0.148
#>  [8,] 0.0976 0.115 0.116 0.148
#>  [9,] 0.0976 0.115 0.116 0.148
#> [10,] 0.0976 0.115 0.116 0.148
#> [11,] 0.0976 0.115 0.116 0.148
#> [12,] 0.0976 0.115 0.116 0.148
#> [13,] 0.0976 0.115 0.116 0.148
#> [14,] 0.0976 0.115 0.116 0.148
#> [15,] 0.0976 0.115 0.116 0.148
#> [16,] 0.0976 0.115 0.116 0.148
#> [17,] 0.0976 0.115 0.116 0.148
#> [18,] 0.0976 0.115 0.116 0.148
#> [19,] 0.0976 0.115 0.116 0.148
#> [20,] 0.0976 0.115 0.116 0.148
#> [21,] 0.0976 0.115 0.116 0.148
#> [22,] 0.0976 0.115 0.116 0.148
#> [23,] 0.0976 0.115 0.116 0.148
#> [24,] 0.0976 0.115 0.116 0.148
#> [25,] 0.0976 0.115 0.116 0.148
#> [26,] 0.0976 0.115 0.116 0.148
#> [27,] 0.0976 0.115 0.116 0.148
#> [28,] 0.0976 0.115 0.116 0.148
#> [29,] 0.0976 0.115 0.116 0.148
#> [30,] 0.0976 0.115 0.116 0.148
#> [31,] 0.0976 0.115 0.116 0.148
#> [32,] 0.0976 0.115 0.116 0.148
#> [33,] 0.0976 0.115 0.116 0.148
#> [34,] 0.0976 0.115 0.116 0.148
#> [35,] 0.0976 0.115 0.116 0.148
#> [36,] 0.0976 0.115 0.116 0.148
#> [37,] 0.0976 0.115 0.116 0.148
#> [38,] 0.0976 0.115 0.116 0.148
#> [39,] 0.0976 0.115 0.116 0.148
#> [40,] 0.0976 0.115 0.116 0.148
#> [41,] 0.0976 0.115 0.116 0.148
#> [42,] 0.0976 0.115 0.116 0.148
#> [43,] 0.0976 0.115 0.116 0.148
#> [44,] 0.0976 0.115 0.116 0.148
#> [45,] 0.0976 0.115 0.116 0.148
#> [46,] 0.0976 0.115 0.116 0.148
#> [47,] 0.0976 0.115 0.116 0.148
#> [48,] 0.0976 0.115 0.116 0.148
#> [49,] 0.0976 0.115 0.116 0.148
#> [50,] 0.0976 0.115 0.116 0.148
#> [51,] 0.0976 0.115 0.116 0.148
#> [52,] 0.0976 0.115 0.116 0.148
#> 
#> $target_inf
#> [1] 90000
#> 
#> $target_mig
#> [1] 12000
#> 
#> $init_props
#> [1] 7.00 3.50 1.75 0.50
#> 
#> $init_props_abroad
#> [1] 1.75 2.25 2.25 1.25
#> 
#> $abroad_prop
#> [1] 0.2
#> 
#> $cd4_completeness
#> [1] 1
#> 
#> $rita_completeness
#> [1] 1
#> 
```
