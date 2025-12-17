// Age-independent spline model with RITA state

functions {
  // Function to create transition matrices
  matrix tmat_fct(int t, vector q, matrix d) {
    matrix[13, 13] Lambda = rep_matrix(0, 13, 13);
    
    // Setting diagonal entries to one
    for (i in 8 : 13) {
      Lambda[i, i] = 1;
    }
    
    Lambda[1, 1] = 0;
    Lambda[1, 2] = (1 - q[1]) * (1 - d[1, t]);
    Lambda[1, 3] = q[1] * (1 - d[1, t]);
    Lambda[1, 9] = d[1, t];
    
    Lambda[2, 2] = 0;
    Lambda[2, 4] = (1 - q[1]) * (1 - d[1, t]);
    Lambda[2, 5] = q[1] * (1 - d[1, t]);
    Lambda[2, 9] = d[1, t];
    
    Lambda[3, 3] = 0;
    Lambda[3, 5] = (1 - q[2]) * (1 - d[1, t]);
    Lambda[3, 6] = q[2] * (1 - d[1, t]);
    Lambda[3, 9] = d[1, t];
    
    Lambda[4, 4] = (1 - q[1]) * (1 - d[2, t]);
    Lambda[4, 5] = q[1] * (1 - d[2, t]);
    Lambda[4, 10] = d[2, t];
    
    Lambda[5, 5] = (1 - q[2]) * (1 - d[3, t]);
    Lambda[5, 6] = q[2] * (1 - d[3, t]);
    Lambda[5, 11] = d[3, t];
    
    Lambda[6, 6] = (1 - q[3]) * (1 - d[4, t]);
    Lambda[6, 7] = q[3] * (1 - d[4, t]);
    Lambda[6, 12] = d[4, t];
    
    Lambda[7, 7] = (1 - q[4]) * (1 - d[5, t]);
    Lambda[7, 8] = q[4] * (1 - d[5, t]);
    Lambda[7, 13] = d[5, t];
    
    return Lambda;
  }
}
data {
  int<lower=1> nquar; // number of quarters
  array[nquar] int<lower=0> HIV; //Number of HIV diagnoses
  array[nquar] int<lower=0> AIDS; // AIDS diag data
  array[nquar, 5] int<lower=0> CD4; // diagnosis data by states
  vector<lower=0, upper=1>[4] q; // Fixed progression probabilities
  vector<lower=0>[7] init_prev; // Initial prevalence
  
  // Splines
  int<lower=1> npars_inf;
  int<lower=1> npars_diag;
  matrix[nquar, npars_inf] X_inf;
  matrix[nquar, npars_diag] X_diag;
}
transformed data {
  vector[12] zeroes = rep_vector(0, 12);
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
  vector[npars_diag] beta_diag_5;
  real<lower=0> sigma_diag_1; // smoothing parameter for delta1 spline
  real<lower=0> sigma_diag_2;
  real<lower=0> sigma_diag_3;
  real<lower=0> sigma_diag_4;
  real<lower=0> sigma_diag_5;
  real alpha_diag_1; // Initial value for delta1
  real alpha_diag_2;
  real alpha_diag_3;
  real alpha_diag_4;
  real alpha_diag_5;
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
  vector[nquar] delta_5 = alpha_diag_5 + 0.7 * (X_diag * beta_diag_5);
  matrix<lower=0, upper=1>[5, nquar] d;
  
  // Transform diag probs using logit transformation
  for (t in 1 : nquar) {
    d[ : , t] = inv_logit(
                          [delta_1[t], delta_2[t], delta_3[t], delta_4[t],
                           delta_5[t]]');
  }
  
  // Epidemic state variables
  array[nquar] vector<lower=0>[13] cum_diags = rep_array(rep_vector(0, 13),
                                                         nquar);
  
  // Expected observations
  row_vector<lower=0>[nquar] exp_HIV_dx;
  row_vector<lower=0>[nquar] exp_AIDS_dx;
  array[nquar] simplex[5] exp_p = rep_array(rep_vector(0.20, 5), nquar);
  
  // Expected incidence
  row_vector<lower=0>[nquar] incidence_state_8;
  row_vector<lower=0>[nquar] incidence_state_9;
  row_vector<lower=0>[nquar] incidence_state_10;
  row_vector<lower=0>[nquar] incidence_state_11;
  row_vector<lower=0>[nquar] incidence_state_12;
  row_vector<lower=0>[nquar] incidence_state_13;
  array[nquar] vector<lower=0>[7] prev_mat;
  
  // Initial time period (t=1, y=1)
  vector[13] init_vec = rep_vector(0, 13);
  init_vec[1] = h[1];
  
  vector[13] tmpvec = append_row(init_prev, zeroes[1 : 6]);
  
  matrix[13, 13] Lambda1 = tmat_fct(1, q, d);
  cum_diags[1] = Lambda1' * tmpvec + init_vec;
  
  // Evolution over time
  for (t in 2 : nquar) {
    vector[13] inflow = rep_vector(0, 13);
    inflow[1] = h[t];
    matrix[13, 13] Lambda_t = tmat_fct(t, q, d);
    cum_diags[t] = Lambda_t' * cum_diags[t - 1] + inflow;
  }
  
  // Calculate incidences
  // First time period t=1
  incidence_state_8[1] = cum_diags[1, 8];
  incidence_state_9[1] = cum_diags[1, 9];
  incidence_state_10[1] = cum_diags[1, 10];
  incidence_state_11[1] = cum_diags[1, 11];
  incidence_state_12[1] = cum_diags[1, 12];
  incidence_state_13[1] = cum_diags[1, 13];
  
  // t = 2, ..., nquar
  for (t in 2 : nquar) {
    incidence_state_8[t] = fmax(0, cum_diags[t, 8] - cum_diags[t - 1, 8]);
    incidence_state_9[t] = fmax(0, cum_diags[t, 9] - cum_diags[t - 1, 9]);
    incidence_state_10[t] = fmax(0, cum_diags[t, 10] - cum_diags[t - 1, 10]);
    incidence_state_11[t] = fmax(0, cum_diags[t, 11] - cum_diags[t - 1, 11]);
    incidence_state_12[t] = fmax(0, cum_diags[t, 12] - cum_diags[t - 1, 12]);
    incidence_state_13[t] = fmax(0, cum_diags[t, 13] - cum_diags[t - 1, 13]);
  }
  
  // Prevalence matrix
  for (t in 1 : nquar) {
    prev_mat[t] = cum_diags[t, 1 : 7];
  }
  
  // Expected diagnoses
  exp_HIV_dx = incidence_state_9 + incidence_state_10 + incidence_state_11
               + incidence_state_12 + incidence_state_13;
  exp_AIDS_dx = incidence_state_8;
  
  // Including CD4 counts
  for (t in 1 : nquar) {
    if (exp_HIV_dx[t] > 0) {
      exp_p[t, 1] = incidence_state_9[t] / exp_HIV_dx[t];
      exp_p[t, 2] = incidence_state_10[t] / exp_HIV_dx[t];
      exp_p[t, 3] = incidence_state_11[t] / exp_HIV_dx[t];
      exp_p[t, 4] = incidence_state_12[t] / exp_HIV_dx[t];
      exp_p[t, 5] = incidence_state_13[t] / exp_HIV_dx[t];
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
  alpha_diag_5 ~ normal(-2.0, 0.5);
  beta_diag_1[1] ~ std_normal();
  beta_diag_2[1] ~ std_normal();
  beta_diag_3[1] ~ std_normal();
  beta_diag_4[1] ~ std_normal();
  beta_diag_5[1] ~ std_normal();
  beta_diag_1[2 : npars_diag] ~ normal(0, sigma_diag_1);
  beta_diag_2[2 : npars_diag] ~ normal(0, sigma_diag_2);
  beta_diag_3[2 : npars_diag] ~ normal(0, sigma_diag_3);
  beta_diag_4[2 : npars_diag] ~ normal(0, sigma_diag_4);
  beta_diag_5[2 : npars_diag] ~ normal(0, sigma_diag_5);
  sigma_diag_1 ~ gamma(2, 64);
  sigma_diag_2 ~ gamma(2, 64);
  sigma_diag_3 ~ gamma(2, 64);
  sigma_diag_4 ~ gamma(2, 64);
  sigma_diag_5 ~ gamma(2, 64);
  
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
