// Age-independent model with spline for infection
// and spline for diagnosis probabilities

functions {
  // Function to create transition matrices
  matrix tmat_fct(int t, vector q, matrix d) {
    matrix[9, 9] Lambda = rep_matrix(0, 9, 9);
    
    // Setting diagonal entries to one
    for (i in 5 : 9) {
      Lambda[i, i] = 1;
    }
    
    // Cache (1 - d[i, t]) to avoid repeated computation
    real s1 = 1 - d[1, t];
    real s2 = 1 - d[2, t];
    real s3 = 1 - d[3, t];
    real s4 = 1 - d[4, t];
    
    Lambda[1, 1] = (1 - q[1]) * s1;
    Lambda[1, 2] = q[1] * s1;
    Lambda[1, 6] = d[1, t];
    
    Lambda[2, 2] = (1 - q[2]) * s2;
    Lambda[2, 3] = q[2] * s2;
    Lambda[2, 7] = d[2, t];
    
    Lambda[3, 3] = (1 - q[3]) * s3;
    Lambda[3, 4] = q[3] * s3;
    Lambda[3, 8] = d[3, t];
    
    Lambda[4, 4] = (1 - q[4]) * s4;
    Lambda[4, 5] = q[4] * s4;
    Lambda[4, 9] = d[4, t];
    
    return Lambda;
  }
}
data {
  int<lower=1> nquar; // number of quarters considered
  array[nquar] int<lower=0> HIV; // data with new HIV diag by year
  array[nquar] int<lower=0> AIDS; // data with new AIDS diag by year
  array[nquar, 4] int<lower=0> CD4; // CD4 data by number of years and CD4 states
  vector<lower=0, upper=1>[4] q; // Fixed progression probabilities
  vector<lower=0>[4] init_prev; // Initial prevalence
  
  // Splines
  int<lower=1> npars_inf;
  matrix[nquar, npars_inf] X_inf;
}
parameters {
  // Infection parameters
  vector[npars_inf] beta_inf;
  real log_sigma_inf;
  
  // Diagnosis probability parameters (random walk)
  vector[nquar] delta_raw_1;
  vector[nquar] delta_raw_2;
  vector[nquar] delta_raw_3;
  vector[nquar] delta_raw_4;
  vector[4] log_vardelta_raw;
  real mu_log_var;
  real<lower=0> sigma_log_var;
  real alpha_diag_1;
  real alpha_diag_2;
  real alpha_diag_3;
  real alpha_diag_4;
  
  // AR(1) coefficient shared across diagnosis-probability walks.
  // Constrained to (0, 1) so the process is mean-reverting but persistent.
  real<lower=0, upper=1> rho_diag;
}
transformed parameters {
  // Infection process
  real<lower=0> sigma_inf = exp(log_sigma_inf);
  vector[nquar] lambda = X_inf * beta_inf;
  vector<lower=0>[nquar] h = exp(lambda);
  
  // Step variances on natural scale (non-centred)
  vector[4] log_vardelta = mu_log_var + sigma_log_var * log_vardelta_raw;
  vector<lower=0>[4] vardelta = exp(log_vardelta);
  // Diagnosis probabilities via random walk
  vector[nquar] delta_1;
  vector[nquar] delta_2;
  vector[nquar] delta_3;
  vector[nquar] delta_4;
  // AR(1) non-centred parameterisation.
  // t = 1: initial condition drawn from stationary distribution
  //   (variance = step_variance / (1 - rho^2)).
  // t > 1: mean reversion toward alpha_diag_k prevents endpoint drift.
  real inv_stationary = inv_sqrt(1 - square(rho_diag));
  delta_1[1] = alpha_diag_1 + sqrt(vardelta[1]) * inv_stationary * delta_raw_1[1];
  delta_2[1] = alpha_diag_2 + sqrt(vardelta[2]) * inv_stationary * delta_raw_2[1];
  delta_3[1] = alpha_diag_3 + sqrt(vardelta[3]) * inv_stationary * delta_raw_3[1];
  delta_4[1] = alpha_diag_4 + sqrt(vardelta[4]) * inv_stationary * delta_raw_4[1];
  for (t in 2 : nquar) {
    delta_1[t] = alpha_diag_1 + rho_diag * (delta_1[t - 1] - alpha_diag_1) + sqrt(vardelta[1]) * delta_raw_1[t];
    delta_2[t] = alpha_diag_2 + rho_diag * (delta_2[t - 1] - alpha_diag_2) + sqrt(vardelta[2]) * delta_raw_2[t];
    delta_3[t] = alpha_diag_3 + rho_diag * (delta_3[t - 1] - alpha_diag_3) + sqrt(vardelta[3]) * delta_raw_3[t];
    delta_4[t] = alpha_diag_4 + rho_diag * (delta_4[t - 1] - alpha_diag_4) + sqrt(vardelta[4]) * delta_raw_4[t];
  }
  matrix<lower=0, upper=1>[4, nquar] d;
  
  // Transform diag probs using logit transformation
  for (t in 1 : nquar) {
    d[ : , t] = inv_logit([delta_1[t], delta_2[t], delta_3[t], delta_4[t]]');
  }
  
  // Epidemic state variables
  array[nquar] vector<lower=0>[9] cum_diags = rep_array(rep_vector(0, 9),
                                                        nquar);
  
  // Expected observations
  row_vector<lower=0>[nquar] exp_HIV_dx;
  row_vector<lower=0>[nquar] exp_AIDS_dx;
  array[nquar] simplex[4] exp_p = rep_array(rep_vector(0.25, 4), nquar);
  
  // Expected incidence
  row_vector<lower=0>[nquar] incidence_state_5;
  row_vector<lower=0>[nquar] incidence_state_6;
  row_vector<lower=0>[nquar] incidence_state_7;
  row_vector<lower=0>[nquar] incidence_state_8;
  row_vector<lower=0>[nquar] incidence_state_9;
  array[nquar] vector[4] prev_mat;
  
  // Initial time period (t=1, y=1)
  vector[9] init_vec = rep_vector(0, 9);
  init_vec[1] = h[1];
  
  vector[9] tmpvec = append_row(init_prev, rep_vector(0, 5));
  
  matrix[9, 9] Lambda1 = tmat_fct(1, q, d);
  cum_diags[1] = Lambda1' * tmpvec + init_vec;
  
  // Evolution over time
  for (t in 2 : nquar) {
    matrix[9, 9] Lambda_t = tmat_fct(t, q, d);
    cum_diags[t] = Lambda_t' * cum_diags[t - 1];
    cum_diags[t, 1] += h[t];
  }
  
  // Calculate incidences
  // First time period
  incidence_state_5[1] = cum_diags[1, 5];
  incidence_state_6[1] = cum_diags[1, 6];
  incidence_state_7[1] = cum_diags[1, 7];
  incidence_state_8[1] = cum_diags[1, 8];
  incidence_state_9[1] = cum_diags[1, 9];
  
  // Subsequent periods
  for (t in 2 : nquar) {
    incidence_state_5[t] = cum_diags[t, 5] - cum_diags[t - 1, 5];
    incidence_state_6[t] = cum_diags[t, 6] - cum_diags[t - 1, 6];
    incidence_state_7[t] = cum_diags[t, 7] - cum_diags[t - 1, 7];
    incidence_state_8[t] = cum_diags[t, 8] - cum_diags[t - 1, 8];
    incidence_state_9[t] = cum_diags[t, 9] - cum_diags[t - 1, 9];
  }
  
  // Prevalence matrix
  for (t in 1 : nquar) {
    prev_mat[t] = cum_diags[t, 1 : 4];
  }
  
  // Expected diagnoses
  exp_HIV_dx = incidence_state_6 + incidence_state_7 + incidence_state_8
               + incidence_state_9;
  exp_AIDS_dx = incidence_state_5;
  
  // Including CD4 counts
  for (t in 1 : nquar) {
    if (exp_HIV_dx[t] > 0) {
      exp_p[t, 1] = incidence_state_6[t] / exp_HIV_dx[t];
      exp_p[t, 2] = incidence_state_7[t] / exp_HIV_dx[t];
      exp_p[t, 3] = incidence_state_8[t] / exp_HIV_dx[t];
      exp_p[t, 4] = incidence_state_9[t] / exp_HIV_dx[t];
    }
  }
}
model {
  // PRIORS
  
  // Infection priors
  beta_inf[1] ~ normal(5, 1.5); // Vaguely informative priors: corresponds to quarter inf in (exp(2),exp(8)) = (8,4390)
  beta_inf[2 : npars_inf] ~ normal(0, sigma_inf);
  log_sigma_inf ~ normal(-1, 1);
  
  // Diagnosis probability priors — centred on midpoints of observed empirical ranges.
  // State 1 (CD4 ≥ 500):    2.5–10%  → logit midpoint ~ -2.90
  // State 2 (CD4 350–500):  2.5–11%  → logit midpoint ~ -2.75
  // State 3 (CD4 200–349):  5–25%    → logit midpoint ~ -2.00
  // State 4 (CD4 < 200):    20–60%   → logit midpoint ~ -0.50
  alpha_diag_1 ~ normal(-2.9, 0.5);
  alpha_diag_2 ~ normal(-2.75, 0.5);
  alpha_diag_3 ~ normal(-2.0, 0.5);
  alpha_diag_4 ~ normal(-0.5, 0.5);
  delta_raw_1 ~ std_normal();
  delta_raw_2 ~ std_normal();
  delta_raw_3 ~ std_normal();
  delta_raw_4 ~ std_normal();
  // mu_log_var: centres exp(-4.0) ~ 0.018 => step SD ~0.13 logit/quarter,
  // stationary SD ~0.26 logit. SD=0.5 avoids the heavy upper tail of normal(-3.5,1)
  // which allowed unrealistically large swings in diagnosis probability.
  mu_log_var ~ normal(-4.0, 0.5);
  sigma_log_var ~ normal(0, 0.5);
  log_vardelta_raw ~ std_normal();
  
  // AR(1) coefficient: prior with mode ~ 0.86 (high persistence) but not a
  // pure random walk. beta(7, 2) has mode = 6/7 ~ 0.857 and mean = 7/9 ~ 0.78,
  // placing almost no mass below 0.5 and ruling out the RW1 boundary at 1.
  rho_diag ~ beta(7, 2);
  
  // Likelihood
  HIV ~ poisson(exp_HIV_dx);
  AIDS ~ poisson(exp_AIDS_dx);
  
  for (t in 1 : nquar) {
    CD4[t,  : ] ~ multinomial(exp_p[t,  : ]);
  }
}
generated quantities {
  // Log-lik calculations to get WAIC / LOO/ DIC
  vector[nquar] log_lik;
  for (t in 1 : nquar) {
    log_lik[t] = poisson_lpmf(HIV[t] | exp_HIV_dx[t])
                 + poisson_lpmf(AIDS[t] | exp_AIDS_dx[t])
                 + multinomial_lpmf(CD4[t,  : ] | exp_p[t,  : ]);
  }
}
