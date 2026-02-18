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
  int<lower=1> npars_diag;
  matrix[nquar, npars_inf] X_inf;
  matrix[nquar, npars_diag] X_diag;
}
parameters {
  // Infection parameters
  vector[npars_inf] beta_inf;
  real log_sigma_inf;
  
  // Diagnosis probability parameters
  vector[npars_diag] beta_diag_1;
  vector[npars_diag] beta_diag_2;
  vector[npars_diag] beta_diag_3;
  vector[npars_diag] beta_diag_4;
  real<lower=0> sigma_diag_1; // smoothing parameter for delta1 spline
  real<lower=0> sigma_diag_2;
  real<lower=0> sigma_diag_3;
  real<lower=0> sigma_diag_4;
  real alpha_diag_1; // Initial value for delta1
  real alpha_diag_2;
  real alpha_diag_3;
  real alpha_diag_4;
}
transformed parameters {
  // Infection process
  real<lower=0> sigma_inf = exp(log_sigma_inf);
  vector[nquar] lambda = X_inf * beta_inf;
  vector<lower=0>[nquar] h = exp(lambda);
  
  // Diagnosis probabilities
  vector[nquar] delta_1 = alpha_diag_1 + 0.5 * (X_diag * beta_diag_1);
  vector[nquar] delta_2 = alpha_diag_2 + 0.5 * (X_diag * beta_diag_2);
  vector[nquar] delta_3 = alpha_diag_3 + 0.5 * (X_diag * beta_diag_3);
  vector[nquar] delta_4 = alpha_diag_4 + 0.7 * (X_diag * beta_diag_4);
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
  log_sigma_inf ~ normal(-2, 0.5);
  
  // Diagnosis probability priors
  alpha_diag_1 ~ normal(-3.2, 0.5);
  alpha_diag_2 ~ normal(-3.2, 0.5);
  alpha_diag_3 ~ normal(-3.0, 0.5);
  alpha_diag_4 ~ normal(-2.5, 0.5);
  beta_diag_1[1] ~ std_normal();
  beta_diag_2[1] ~ std_normal();
  beta_diag_3[1] ~ std_normal();
  beta_diag_4[1] ~ std_normal();
  beta_diag_1[2 : npars_diag] ~ normal(0, sigma_diag_1);
  beta_diag_2[2 : npars_diag] ~ normal(0, sigma_diag_2);
  beta_diag_3[2 : npars_diag] ~ normal(0, sigma_diag_3);
  beta_diag_4[2 : npars_diag] ~ normal(0, sigma_diag_4);
  sigma_diag_1 ~ gamma(2, 8);
  sigma_diag_2 ~ gamma(2, 8);
  sigma_diag_3 ~ gamma(2, 8);
  sigma_diag_4 ~ gamma(2, 8);
  
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
