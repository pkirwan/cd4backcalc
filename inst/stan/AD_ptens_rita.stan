// Age-dependent tensor product model with RITA state

data {
  int<lower=1> nquar; // number of years 
  int<lower=1> nage; // number of ages 
  array[nquar, nage] int<lower=0> HIV_age; // data with new HIV_age diag by year and age
  array[nquar, nage] int<lower=0> AIDS_age; // data with new AIDS_age diag by year and age
  array[nquar, nage, 5] int<lower=0> CD4_age; // CD4 data by number of years and ages and diagnosed CD4 states(5)
  matrix<lower=0, upper=1>[4, nage] q_age; // matrix of progressions for different ages.. fixed
  matrix<lower=0>[7, nage] init_prev_age; // Initial prevalence by age at inf and latent CD4 states(7)
  
  // Splines
  int<lower=1> npars_inf;
  matrix[nquar * nage, npars_inf] X_inf; // design matrix ... 50 is number of inf params
  matrix[npars_inf - 1, 2 * npars_inf - 2] S1; // NOT SURE IF BOUNDARY CONSISTENT
  
}
transformed data {
  row_vector[6] empty_vec = rep_row_vector(0, 6);
  vector[7] empty_vec1 = rep_vector(0, 7);
  vector[npars_inf - 1] zero_vec = rep_vector(0, npars_inf - 1);
  int nagroup = 4; // number of age groups
  int np1 = npars_inf - 1;
  
  // Pre-split penalty matrices (symmetrized to avoid fp noise)
  matrix[np1, np1] S1a = 0.5 * (S1[1 : np1, 1 : np1] + S1[1 : np1, 1 : np1]');
  matrix[np1, np1] S1b = 0.5 * (S1[1 : np1, npars_inf : (2 * npars_inf - 2)] + S1[1 : np1, npars_inf : (2 * npars_inf - 2)]');
  
  // Age-group mapping and fast quarter->year index (birthdays at t = 1, 5, 9, ...)
  array[nage] int age_group_lookup;
  array[nquar] int year_index;
  for (t in 1 : nquar) {
    year_index[t] = (t - 1) %/% 4;
  }
  
  // age 1-10 -> 1, 11-20 -> 2, 21-30 -> 3, 31+ -> 4
  for (a0 in 1 : nage) {
    age_group_lookup[a0] = min(4, 1 + (a0 - 1) %/% 10);
  }
  
  // Pre-computed current-age lookup: age_at_t[t0, t, a0] would be 3D and large,
  // so store only the year_diff component (independent of a0).
  // current age for cohort (t0, a0) at time t = min(a0 + year_diff[t0, t], nage)
  array[nquar, nquar] int year_diff;
  for (t0 in 1 : nquar) {
    for (t in 1 : nquar) {
      year_diff[t0, t] = year_index[t] - year_index[t0];
    }
  }
}
parameters {
  // Infection parameters
  vector[npars_inf] beta_inf; //beta are the infection parameters
  vector<lower=0>[2] sigma_inf; // smoothing parameters
  
  // Diagnosis probability parameters (random walk)
  // Non-centred increments for each state
  vector[nquar] delta_raw_1;
  vector[nquar] delta_raw_2;
  vector[nquar] delta_raw_3;
  vector[nquar] delta_raw_4;
  vector[nquar] delta_raw_5;
  // Step variances: non-centred log-normal avoids funnel geometry.
  // log_vardelta = mu_log_var + sigma_log_var * log_vardelta_raw,
  // so vardelta = exp(log_vardelta) is always strictly positive.
  vector[5] log_vardelta_raw; // non-centred increments
  real mu_log_var;          // shared hyperprior mean (log scale)
  real<lower=0> sigma_log_var; // shared hyperprior SD (log scale)
  // Initial level intercepts
  real alpha_diag_1;
  real alpha_diag_2;
  real alpha_diag_3;
  real alpha_diag_4;
  real alpha_diag_5;
  array[5] sum_to_zero_vector[4] gamma_age_raw;
  real<lower=0> sigma_age; // shared shrinkage for age effects
  
  // AR(1) coefficient shared across diagnosis-probability walks.
  // Constrained to (0, 1) so the process is mean-reverting but persistent.
  real<lower=0, upper=1> rho_diag;
  
  // False recency rate
  real logit_frr;
}
transformed parameters {
  // Infection parameters
  vector[nquar * nage] lambda = X_inf * beta_inf;
  vector<lower=0>[nquar * nage] h = exp(lambda); //vector of expected log responses
  
  // Step variances on natural scale (non-centred)
  vector[5] log_vardelta = mu_log_var + sigma_log_var * log_vardelta_raw;
  vector<lower=0>[5] vardelta = exp(log_vardelta);
  
  // False recency rate
  real<lower=0, upper=1> frr = inv_logit(logit_frr);
  
  // Diagnosis probabilities
  array[nquar, nage] row_vector[6] exp_diag; // expected number of diag at time t age a state k
  
  // Expected observations
  matrix[nquar, nage] exp_HIV_dx; // expected number of HIV_age diag at time t age a
  matrix[nquar, nage] exp_AIDS_dx; // expected number of HIV_age diag at time t age a
  array[nquar, nage] simplex[5] exp_p; // proportion of HIV_age diagnoses in state k (time t and age a)
  array[4] matrix[nquar, 5] d; // matrix of diagnoses for different ages and age-groups
  
  // Smoothing parameters
  vector<lower=0>[2] tau; // smoothing parameters
  
  tau = inv(square(sigma_inf));
  
  // DIAGNOSES
  // Diagnosis probabilities via random walk (non-centred parameterisation).
  // Initial value uses a fixed small step scale (0.2) so the starting
  // uncertainty is not driven by vardelta, which governs evolution.
  vector[nquar] delta_1;
  vector[nquar] delta_2;
  vector[nquar] delta_3;
  vector[nquar] delta_4;
  vector[nquar] delta_5;
  // AR(1) non-centred parameterisation.
  // t = 1: initial condition drawn from stationary distribution
  //   (variance = step_variance / (1 - rho^2)).
  // t > 1: mean reversion toward alpha_diag_k prevents endpoint drift.
  real inv_stationary = inv_sqrt(1 - square(rho_diag));
  delta_1[1] = alpha_diag_1 + sqrt(vardelta[1]) * inv_stationary * delta_raw_1[1];
  delta_2[1] = alpha_diag_2 + sqrt(vardelta[2]) * inv_stationary * delta_raw_2[1];
  delta_3[1] = alpha_diag_3 + sqrt(vardelta[3]) * inv_stationary * delta_raw_3[1];
  delta_4[1] = alpha_diag_4 + sqrt(vardelta[4]) * inv_stationary * delta_raw_4[1];
  delta_5[1] = alpha_diag_5 + sqrt(vardelta[5]) * inv_stationary * delta_raw_5[1];
  for (t in 2 : nquar) {
    delta_1[t] = alpha_diag_1 + rho_diag * (delta_1[t - 1] - alpha_diag_1) + sqrt(vardelta[1]) * delta_raw_1[t];
    delta_2[t] = alpha_diag_2 + rho_diag * (delta_2[t - 1] - alpha_diag_2) + sqrt(vardelta[2]) * delta_raw_2[t];
    delta_3[t] = alpha_diag_3 + rho_diag * (delta_3[t - 1] - alpha_diag_3) + sqrt(vardelta[3]) * delta_raw_3[t];
    delta_4[t] = alpha_diag_4 + rho_diag * (delta_4[t - 1] - alpha_diag_4) + sqrt(vardelta[4]) * delta_raw_4[t];
    delta_5[t] = alpha_diag_5 + rho_diag * (delta_5[t - 1] - alpha_diag_5) + sqrt(vardelta[5]) * delta_raw_5[t];
  }
  
  // Built-in zero-sum age effects by diagnosis state.
  // The sqrt(4 / 3) factor corrects the constrained-vector variance so each
  // age effect has marginal SD sigma_age rather than sigma_age * sqrt(3 / 4).
  matrix[4, 5] gamma_age;
  real age_effect_scale = sigma_age * sqrt(4.0 / 3.0);
  for (k in 1 : 5) {
    for (ag in 1 : 4) {
      gamma_age[ag, k] = age_effect_scale * gamma_age_raw[k][ag];
    }
  }
  
  // Diagnosis probabilities by age group (ag), state (k), and time (t).
  // Stored as matrix[nquar, 5] so the inner time-loop accesses column k,
  // which is contiguous in column-major memory (cache-friendly).
  for (ag in 1 : nagroup) {
    d[ag,  : , 1] = inv_logit(delta_1 + gamma_age[ag, 1]);
    d[ag,  : , 2] = inv_logit(delta_2 + gamma_age[ag, 2]);
    d[ag,  : , 3] = inv_logit(delta_3 + gamma_age[ag, 3]);
    d[ag,  : , 4] = inv_logit(delta_4 + gamma_age[ag, 4]);
    d[ag,  : , 5] = inv_logit(delta_5 + gamma_age[ag, 5]);
  }
  
  // NOW NEED TO CREATE THE TRANSITION MATRICES AND CALCULATE EXPECTED DIAGNOSES
  // Initialise exp_diag to zero
  exp_diag = rep_array(empty_vec, nquar, nage);
  
  for (t0 in 1 : nquar) {
    for (a0 in 1 : nage) {
      // Progression rates depend only on age at infection
      real q1 = q_age[1, a0];
      real q2 = q_age[2, a0];
      real q3 = q_age[3, a0];
      real q4 = q_age[4, a0];
      real nq1 = 1 - q1;
      real nq2 = 1 - q2;
      real nq3 = 1 - q3;
      real nq4 = 1 - q4;
      
      vector[7] lat = empty_vec1;
      vector[7] lat_next;
      
      // At t0 == 1: seed latent state with prevalent infections and record their diagnoses
      if (t0 == 1) {
        int ag0 = age_group_lookup[a0];
        real d1 = d[ag0, 1, 1];
        real d2 = d[ag0, 1, 2];
        real d3 = d[ag0, 1, 3];
        real d4 = d[ag0, 1, 4];
        real d5 = d[ag0, 1, 5];
        real nd1 = 1 - d1;
        real nd2 = 1 - d2;
        real nd3 = 1 - d3;
        real nd4 = 1 - d4;
        real nd5 = 1 - d5;
        
        real v1 = init_prev_age[1, a0];
        real v2 = init_prev_age[2, a0];
        real v3 = init_prev_age[3, a0];
        real v4 = init_prev_age[4, a0];
        real v5 = init_prev_age[5, a0];
        real v6 = init_prev_age[6, a0];
        real v7 = init_prev_age[7, a0];
        
        // Diagnoses from prevalent infections at t=1
        exp_diag[1, a0, 1] += (v1 + v2 + v3) * d1;
        exp_diag[1, a0, 2] += v4 * d2;
        exp_diag[1, a0, 3] += v5 * d3;
        exp_diag[1, a0, 4] += v6 * d4;
        exp_diag[1, a0, 5] += v7 * d5;
        exp_diag[1, a0, 6] += v7 * nd5 * q4;
        
        // Survivors of prevalent pool after diagnosis and progression
        lat[1] = 0;
        lat[2] = (nd1 * nq1) * v1;
        lat[3] = (nd1 * q1) * v1;
        lat[4] = (nd1 * nq1) * v2 + (nd2 * nq1) * v4;
        lat[5] = (nd1 * q1) * v2 + (nd1 * nq2) * v3 + (nd2 * q1) * v4
                 + (nd3 * nq2) * v5;
        lat[6] = (nd1 * q2) * v3 + (nd3 * q2) * v5 + (nd4 * nq3) * v6;
        lat[7] = (nd4 * q3) * v6 + (nd5 * nq4) * v7;
      }
      
      // New infections at t0 enter state 1 (diagnosed only from t0+1 onwards)
      lat[1] += h[(a0 - 1) * nquar + t0];
      
      // Propagate cohort (t0, a0) forward from t0+1
      for (t in (t0 + 1) : nquar) {
        int a = min(a0 + year_diff[t0, t], nage);
        int ag = age_group_lookup[a];
        
        real d1 = d[ag, t, 1];
        real d2 = d[ag, t, 2];
        real d3 = d[ag, t, 3];
        real d4 = d[ag, t, 4];
        real d5 = d[ag, t, 5];
        real nd1 = 1 - d1;
        real nd2 = 1 - d2;
        real nd3 = 1 - d3;
        real nd4 = 1 - d4;
        real nd5 = 1 - d5;
        
        // Diagnoses (RITA: states 1-3 share d1)
        exp_diag[t, a, 1] += (lat[1] + lat[2] + lat[3]) * d1;
        exp_diag[t, a, 2] += lat[4] * d2;
        exp_diag[t, a, 3] += lat[5] * d3;
        exp_diag[t, a, 4] += lat[6] * d4;
        exp_diag[t, a, 5] += lat[7] * d5;
        exp_diag[t, a, 6] += lat[7] * nd5 * q4;
        
        // Progression among undiagnosed
        lat_next[1] = 0;
        lat_next[2] = (nd1 * nq1) * lat[1];
        lat_next[3] = (nd1 * q1) * lat[1];
        lat_next[4] = (nd1 * nq1) * lat[2] + (nd2 * nq1) * lat[4];
        lat_next[5] = (nd1 * q1) * lat[2] + (nd1 * nq2) * lat[3]
                      + (nd2 * q1) * lat[4] + (nd3 * nq2) * lat[5];
        lat_next[6] = (nd1 * q2) * lat[3] + (nd3 * q2) * lat[5]
                      + (nd4 * nq3) * lat[6];
        lat_next[7] = (nd4 * q3) * lat[6] + (nd5 * nq4) * lat[7];
        
        lat = lat_next;
      }
    }
  }
  
  // Expected HIV and AIDS diagnoses by (t,a)
  for (t in 1 : nquar) {
    for (a in 1 : nage) {
      exp_HIV_dx[t, a] = sum(exp_diag[t, a, 1 : 5]);
      exp_AIDS_dx[t, a] = exp_diag[t, a, 6];
    }
  }
  
  // Proportions for multinomial CD4 distribution
  for (t in 1 : nquar) {
    exp_p[t, 1, 1] = 1;
    exp_p[t, 1, 2 : 5] = rep_vector(0, 4);
  }
  
  for (t in 1 : nquar) {
    for (a in 2 : nage) {
      real total = exp_HIV_dx[t, a];
      if (total > 1e-6) {
        // Add small epsilon to avoid exact zeros in the simplex,
        // which cause NaN via 0*log(0) in the multinomial likelihood
        vector[5] raw_p = to_vector(exp_diag[t, a, 1 : 5]) + 1e-12;
        vector[5] p1 = raw_p / sum(raw_p);
        exp_p[t, a, 1] = p1[1] + frr * (1 - p1[1]);
        exp_p[t, a, 2] = (1 - frr) * p1[2];
        exp_p[t, a, 3] = (1 - frr) * p1[3];
        exp_p[t, a, 4] = (1 - frr) * p1[4];
        exp_p[t, a, 5] = (1 - frr) * p1[5];
      } else {
        // safeguard; should be rare outside a==1, but avoids inv(0)
        exp_p[t, a] = rep_vector(0.2, 5);
      }
    }
  }
}
model {
  // Infection surface priors.
  // beta_inf[1] is the log-scale intercept (constant column of X_inf after
  // jagam diagonalization). Observed incidence ~100-300/quarter for younger
  // ages => log-mean ~ log(150) ~ 5.0. normal(5.0, 1.5) covers [7, 3000]
  // at ±2 SD, well-matched to the empirical range.
  sigma_inf ~ normal(0, 5);
  beta_inf[1] ~ normal(5.0, 1.5);
  
  matrix[np1, np1] K1_prec = S1a * tau[1] + S1b * tau[2];
  beta_inf[2 : npars_inf] ~ multi_normal_prec(zero_vec, K1_prec);
  
  // Diagnosis priors
  // Initial level intercepts — centred on midpoints of observed empirical ranges.
  // State 1 (recently acquired / RITA+): 5–30%   → logit midpoint ~ -1.75
  // State 2 (CD4 ≥ 500):                2.5–10%  → logit midpoint ~ -2.90
  // State 3 (CD4 350–500):              2.5–11%  → logit midpoint ~ -2.75
  // State 4 (CD4 200–349):              5–25%    → logit midpoint ~ -2.00
  // State 5 (CD4 < 200):                20–60%   → logit midpoint ~ -0.50
  alpha_diag_1 ~ normal(-1.75, 0.5);
  alpha_diag_2 ~ normal(-2.9, 0.5);
  alpha_diag_3 ~ normal(-2.75, 0.5);
  alpha_diag_4 ~ normal(-2.0, 0.5);
  alpha_diag_5 ~ normal(-0.5, 0.5);
  
  // Non-centred RW increments
  delta_raw_1 ~ std_normal();
  delta_raw_2 ~ std_normal();
  delta_raw_3 ~ std_normal();
  delta_raw_4 ~ std_normal();
  delta_raw_5 ~ std_normal();
  
  // Hierarchical prior on log step variances (non-centred).
  // mu_log_var: centres exp(-4.0) ~ 0.018 => step SD ~0.13 logit/quarter,
  // stationary SD ~0.26 logit. SD=0.5 avoids the heavy upper tail of normal(-3.5,1)
  // which allowed unrealistically large swings in diagnosis probability.
  // sigma_log_var ~ normal(0, 0.5): allows ~2-fold between-state differences.
  mu_log_var ~ normal(-4.0, 0.5);
  sigma_log_var ~ normal(0, 0.5);
  log_vardelta_raw ~ std_normal();
  
  // AR(1) coefficient: prior with mode ~ 0.86 (high persistence) but not a
  // pure random walk. beta(7, 2) has mode = 6/7 ~ 0.857 and mean = 7/9 ~ 0.78,
  // placing almost no mass below 0.5 and ruling out the RW1 boundary at 1.
  rho_diag ~ beta(7, 2);
  
  // False recency rate prior centred on logit(0.019)
  logit_frr ~ normal(-3.94, 0.32);
  
  for (k in 1 : 5) {
    gamma_age_raw[k] ~ std_normal();
  }
  sigma_age ~ normal(0, 0.5);
  
  // Likelihood of model
  for (t in 1 : nquar) {
    // a is in 2 ages, because for age = 1 -> at some t get Po(0)
    // Individuals age at times 5,9,13...
    // So if individuals inf at t=5 at age=1
    // Cannot be dx at t=5 at age age=1 (dx t > 6)
    // remove age = 1 from likelihood
    HIV_age[t, 2 : nage] ~ poisson(exp_HIV_dx[t, 2 : nage]);
    AIDS_age[t, 2 : nage] ~ poisson(exp_AIDS_dx[t, 2 : nage]);
    for (a in 2 : nage) {
      CD4_age[t, a,  : ] ~ multinomial(exp_p[t, a,  : ]);
    }
  }
}
generated quantities {
  // Log-lik calculations to get WAIC / LOO/ DIC
  vector[nquar * nage] log_lik; // has to be vector to be consistent with R pckg
  for (t in 1 : nquar) {
    log_lik[(1 - 1) * nquar + t] = 0;
    
    for (a in 2 : nage) {
      log_lik[(a - 1) * nquar + t] = poisson_lpmf(HIV_age[t, a] |
                                                  exp_HIV_dx[t, a])
                                     + poisson_lpmf(AIDS_age[t, a] |
                                                    exp_AIDS_dx[t, a])
                                     + multinomial_lpmf(CD4_age[t, a,  : ] |
                                                        exp_p[t, a,  : ]);
    }
  }
}
