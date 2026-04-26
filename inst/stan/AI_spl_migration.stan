// Age-independent spline model including migration

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
  array[nquar, 2] int<lower=0> HIV; // data with new HIV diag by year and year of arrival
  array[nquar, 2] int<lower=0> AIDS; // data with new AIDS diag by year and year of arrival
  array[nquar, 2, 4] int<lower=0> CD4; // CD4 data by number of years and ages and CD4 states
  array[nquar] int<lower=0> ABROAD; // Prior diagnoses from abroad
  vector<lower=0, upper=1>[4] q; // Fixed progression probabilities
  vector<lower=0>[4] init_prev; // Initial prevalence (as for age-indept model)
  vector<lower=0>[4] init_prev_abroad; // Initial prevalence (as for age-indept model)
  
  // Splines
  int<lower=1> npars_inf;
  int<lower=1> npars_mig;
  int<lower=1> npars_diag;
  int<lower=1> npars_cob;
  int<lower=1> npars_mp;
  matrix[nquar, npars_inf] X_inf;
  matrix[nquar, npars_mig] X_mig;
  matrix[nquar, npars_diag] X_diag;
  matrix[nquar, npars_cob] X_cob;
  matrix[nquar, npars_mp] X_mp;
  
  // Prior information for evidence synthesis
  vector[nquar] eta_mean;
  real<lower=0> eta_sd;
}
transformed data {
  // No pre-computed constants needed
}
parameters {
  // Infection parameters
  vector[npars_inf] beta_inf; // infection parameters
  real log_sigma_inf;
  
  // Migration parameters
  vector[npars_mig] beta_mig;
  real log_sigma_mig;
  
  // Diagnosis probability parameters
  vector[npars_diag] beta_diag_1;
  vector[npars_diag] beta_diag_2;
  vector[npars_diag] beta_diag_3;
  vector[npars_diag] beta_diag_4;
  real<lower=0> sigma_diag_1; // smoothing parameter for delta1 spline
  real<lower=0> sigma_diag_2; // smoothing parameter for delta2 spline
  real<lower=0> sigma_diag_3; // smoothing parameter for delta3 spline
  real<lower=0> sigma_diag_4; // smoothing parameter for delta4 spline
  real alpha_diag_1; // Initial value for delta1
  real alpha_diag_2; // Initial value for delta2
  real alpha_diag_3; // Initial value for delta3
  real alpha_diag_4; // Initial value for delta4
  vector[4] alpha_diag_mig; // migrant-specific diagnosis intercept offset per state
  
  // Proportion UK-acquired parameters
  vector[npars_cob] beta_cob; // Proportion UK-acquired for each quarter
  real<lower=0> sigma_cob; // Random walk standard deviation
  real alpha_cob;
  real<lower=-1, upper=1> rho_eta;
  
  // Probability of migration into each compartment
  vector[npars_mp] beta_mp_1;
  vector[npars_mp] beta_mp_2;
  vector[npars_mp] beta_mp_3;
  vector[npars_mp] beta_mp_5;
  real<lower=0> sigma_mp_1;
  real<lower=0> sigma_mp_2;
  real<lower=0> sigma_mp_3;
  real<lower=0> sigma_mp_5;
  vector[4] alpha_mp;
  
  // Parameters for misclassification of abroad status
  real alpha_detect; // Intercept for detection probability
  vector[npars_mp] beta_detect; // Time-varying component 
  real<lower=0> sigma_detect; // Smoothing parameter
}
transformed parameters {
  // Infection process
  real<lower=0> sigma_inf = exp(log_sigma_inf);
  vector[nquar] lambda = X_inf * beta_inf;
  vector<lower=0>[nquar] h = exp(lambda);
  
  // Migration process
  real<lower=0> sigma_mig = exp(log_sigma_mig);
  vector[nquar] gamma = X_mig * beta_mig;
  vector<lower=0>[nquar] omega = exp(gamma);
  
  // Diagnosis probabilities
  vector[nquar] delta_1 = alpha_diag_1 + (X_diag * beta_diag_1);
  vector[nquar] delta_2 = alpha_diag_2 + (X_diag * beta_diag_2);
  vector[nquar] delta_3 = alpha_diag_3 + (X_diag * beta_diag_3);
  vector[nquar] delta_4 = alpha_diag_4 + (X_diag * beta_diag_4);
  matrix<lower=0, upper=1>[4, nquar] d_nat;
  matrix<lower=0, upper=1>[4, nquar] d_mig;
  
  // Transform diag probs using logit transformation
  for (t in 1 : nquar) {
    d_nat[ : , t] = inv_logit([delta_1[t], delta_2[t], delta_3[t], delta_4[t]]');
    d_mig[ : , t] = inv_logit([delta_1[t] + alpha_diag_mig[1],
                               delta_2[t] + alpha_diag_mig[2],
                               delta_3[t] + alpha_diag_mig[3],
                               delta_4[t] + alpha_diag_mig[4]]');
  }
  
  // Proportion UK-acquired
  vector[nquar] eta = alpha_cob + 0.5 * (X_cob * beta_cob);
  vector<lower=0>[nquar] p;
  
  // Migration probabilities using softmax, reference = compartment 4
  vector[nquar] theta_1 = alpha_mp[2] + X_mp * beta_mp_1;
  vector[nquar] theta_2 = alpha_mp[3] + X_mp * beta_mp_2;
  vector[nquar] theta_3 = alpha_mp[4] + X_mp * beta_mp_3;
  vector[nquar] theta_5 = alpha_mp[1] + X_mp * beta_mp_5;
  matrix<lower=0>[5, nquar] m;
  
  for (t in 1 : nquar) {
    vector[5] logits;
    logits[1] = theta_1[t];
    logits[2] = theta_2[t];
    logits[3] = theta_3[t];
    logits[4] = 0;  // reference category
    logits[5] = theta_5[t];
    m[ : , t] = softmax(logits);
  }
  
  // Epidemic state variables
  array[nquar, 2] vector<lower=0>[9] cum_diags = rep_array(rep_vector
                                                           (0, 9), nquar, 2);
  
  // Expected observations
  matrix<lower=0>[nquar, 2] exp_HIV_dx;
  matrix<lower=0>[nquar, 2] exp_AIDS_dx;
  vector<lower=0>[nquar] exp_ABROAD_dx;
  array[nquar, 2] simplex[4] exp_p = rep_array(rep_vector(0.25, 4), nquar, 2);
  
  // Detection probability for abroad diagnoses
  vector[nquar] detect_logit = alpha_detect + 0.5 * (X_mp * beta_detect);
  vector<lower=0, upper=1>[nquar] detect_prob;
  vector<lower=0>[nquar] true_abroad_dx;
  matrix<lower=0>[nquar, 2] misclassified_dx = rep_matrix(0, nquar, 2);
  
  for (t in 1 : nquar) {
    detect_prob[t] = inv_logit(detect_logit[t]);
    true_abroad_dx[t] = omega[t] * m[5, t];
    exp_ABROAD_dx[t] = true_abroad_dx[t] * detect_prob[t];
    misclassified_dx[t, 2] = true_abroad_dx[t] * (1 - detect_prob[t]);
  }
  
  // Expected incidence
  matrix<lower=0>[nquar, 2] incidence_state_5 = rep_matrix(0, nquar, 2);
  matrix<lower=0>[nquar, 2] incidence_state_6 = rep_matrix(0, nquar, 2);
  matrix<lower=0>[nquar, 2] incidence_state_7 = rep_matrix(0, nquar, 2);
  matrix<lower=0>[nquar, 2] incidence_state_8 = rep_matrix(0, nquar, 2);
  matrix<lower=0>[nquar, 2] incidence_state_9 = rep_matrix(0, nquar, 2);
  array[nquar, 2] vector[4] prev_mat;
  
  // Initial time period (t=1)
  matrix[9, 9] Lambda1_nat = tmat_fct(1, q, d_nat);
  matrix[9, 9] Lambda1_mig = tmat_fct(1, q, d_mig);
  
  cum_diags[1, 1,  : ] = Lambda1_nat' * append_row(init_prev, rep_vector(0, 5));
  cum_diags[1, 1, 1] += h[1];
  
  // Initial migrant group
  p[1] = exp(eta[1]);
  cum_diags[1, 2,  : ] = Lambda1_mig' * append_row(init_prev_abroad, rep_vector(0, 5));
  cum_diags[1, 2, 1] += h[1] * p[1];
  cum_diags[1, 2, 1 : 4] += omega[1] * m[1 : 4, 1];
  
  for (t in 2 : nquar) {
    matrix[9, 9] Lambda_t_nat = tmat_fct(t, q, d_nat);
    matrix[9, 9] Lambda_t_mig = tmat_fct(t, q, d_mig);
    
    // Evolution for native group
    cum_diags[t, 1,  : ] = Lambda_t_nat' * cum_diags[t - 1, 1,  : ];
    cum_diags[t, 1, 1] += h[t];
    
    // Evolution for migrant group
    p[t] = exp(eta[t]);
    cum_diags[t, 2,  : ] = Lambda_t_mig' * cum_diags[t - 1, 2,  : ];
    cum_diags[t, 2, 1] += h[t] * p[t];
    cum_diags[t, 2, 1 : 4] += omega[t] * m[1 : 4, t];
  }
  
  // Calculate incidences
  for (y in 1 : 2) {
    // First time period
    incidence_state_5[1, y] = cum_diags[1, y, 5];
    incidence_state_6[1, y] = cum_diags[1, y, 6];
    incidence_state_7[1, y] = cum_diags[1, y, 7];
    incidence_state_8[1, y] = cum_diags[1, y, 8];
    incidence_state_9[1, y] = cum_diags[1, y, 9];
    
    // Subsequent periods
    for (t in 2 : nquar) {
      incidence_state_5[t, y] = cum_diags[t, y, 5] - cum_diags[t - 1, y, 5];
      incidence_state_6[t, y] = cum_diags[t, y, 6] - cum_diags[t - 1, y, 6];
      incidence_state_7[t, y] = cum_diags[t, y, 7] - cum_diags[t - 1, y, 7];
      incidence_state_8[t, y] = cum_diags[t, y, 8] - cum_diags[t - 1, y, 8];
      incidence_state_9[t, y] = cum_diags[t, y, 9] - cum_diags[t - 1, y, 9];
    }
  }
  
  // Prevalence matrix
  for (t in 1 : nquar) {
    for (y in 1 : 2) {
      prev_mat[t, y] = cum_diags[t, y, 1 : 4];
    }
  }
  
  // Expected diagnoses
  exp_HIV_dx = incidence_state_6 + incidence_state_7 + incidence_state_8
               + incidence_state_9 + misclassified_dx;
  exp_AIDS_dx = incidence_state_5;
  
  // Including CD4 counts
  for (t in 1 : nquar) {
    // Native population (y=1)
    if (exp_HIV_dx[t, 1] > 0) {
      exp_p[t, 1, 1] = incidence_state_6[t, 1] / exp_HIV_dx[t, 1];
      exp_p[t, 1, 2] = incidence_state_7[t, 1] / exp_HIV_dx[t, 1];
      exp_p[t, 1, 3] = incidence_state_8[t, 1] / exp_HIV_dx[t, 1];
      exp_p[t, 1, 4] = incidence_state_9[t, 1] / exp_HIV_dx[t, 1];
    }
    // Migrant population (y=2): include misclassified arrivals in exp_p,
    // proportional to m[1:4,t], matching what the simulation puts in CD4 data.
    {
      vector[4] raw_p = [incidence_state_6[t, 2], incidence_state_7[t, 2],
                         incidence_state_8[t, 2], incidence_state_9[t, 2]]'
                        + 1e-12;
      if (misclassified_dx[t, 2] > 1e-10) {
        real total_lat = sum(m[1 : 4, t]);
        vector[4] m_lat = (total_lat > 1e-10)
                          ? m[1 : 4, t] / total_lat
                          : rep_vector(0.25, 4);
        raw_p += misclassified_dx[t, 2] * m_lat;
      }
      real total_p = sum(raw_p);
      exp_p[t, 2] = (total_p > 1e-10) ? raw_p / total_p : rep_vector(0.25, 4);
    }
  }
}
model {
  // PRIORS
  
  // Infection priors
  beta_inf[1] ~ normal(5, 1.5); // Vaguely informative priors: corresponds to quarter inf in (exp(2),exp(8)) = (8,4390)
  beta_inf[2 : npars_inf] ~ normal(0, sigma_inf);
  log_sigma_inf ~ normal(-2, 0.5);
  
  // Migration priors
  beta_mig[1] ~ normal(4, 1);
  beta_mig[2 : npars_mig] ~ normal(0, sigma_mig);
  log_sigma_mig ~ normal(-2, 0.5);
  
  // Diagnosis probability priors
  alpha_diag_1 ~ normal(-3.2, 0.5);
  alpha_diag_2 ~ normal(-3.0, 0.5);
  alpha_diag_3 ~ normal(-2.5, 0.5);
  alpha_diag_4 ~ normal(-2.0, 0.5);
  beta_diag_1[1] ~ std_normal();
  beta_diag_2[1] ~ std_normal();
  beta_diag_3[1] ~ std_normal();
  beta_diag_4[1] ~ std_normal();
  beta_diag_1[2 : npars_diag] ~ normal(0, sigma_diag_1);
  beta_diag_2[2 : npars_diag] ~ normal(0, sigma_diag_2);
  beta_diag_3[2 : npars_diag] ~ normal(0, sigma_diag_3);
  beta_diag_4[2 : npars_diag] ~ normal(0, sigma_diag_4);
  sigma_diag_1 ~ gamma(4, 32);
  sigma_diag_2 ~ gamma(4, 32);
  sigma_diag_3 ~ gamma(4, 32);
  sigma_diag_4 ~ gamma(4, 32);
  
  alpha_diag_mig ~ normal(0, 1);
  
  // Proportion UK-acquired priors
  alpha_cob ~ normal(0, 1.5);
  beta_cob[1] ~ std_normal();
  beta_cob[2 : npars_cob] ~ normal(0, sigma_cob);
  sigma_cob ~ gamma(2, 32);
  
  // Approximate melding for proportion UK-acquired - AR(1)
  eta_mean[1] ~ normal(eta[1], eta_sd);
  eta_mean[2 : nquar] ~ normal(
                               eta[2 : nquar]
                               + rho_eta
                                 * (eta_mean[1 : (nquar - 1)]
                                    - eta[1 : (nquar - 1)]),
                               eta_sd * sqrt(1 - rho_eta ^ 2));
  rho_eta ~ normal(0, 0.5);
  
  // Detection probability priors
  alpha_detect ~ normal(1, 1);
  beta_detect[1] ~ normal(0, 0.5);
  beta_detect[2 : npars_mp] ~ normal(0, sigma_detect);
  sigma_detect ~ gamma(2, 32);
  
  // Migration distribution priors
  alpha_mp[1] ~ normal(0, 1.5);
  alpha_mp[2 : 4] ~ std_normal();
  beta_mp_1[1] ~ std_normal();
  beta_mp_2[1] ~ std_normal();
  beta_mp_3[1] ~ std_normal();
  beta_mp_5[1] ~ std_normal();
  beta_mp_1[2 : npars_mp] ~ normal(0, sigma_mp_1);
  beta_mp_2[2 : npars_mp] ~ normal(0, sigma_mp_2);
  beta_mp_3[2 : npars_mp] ~ normal(0, sigma_mp_3);
  beta_mp_5[2 : npars_mp] ~ normal(0, sigma_mp_5);
  sigma_mp_1 ~ gamma(2, 32);
  sigma_mp_2 ~ gamma(2, 32);
  sigma_mp_3 ~ gamma(2, 32);
  sigma_mp_5 ~ gamma(2, 32);
  
  // Likelihood
  for (t in 1 : nquar) {
    for (y in 1 : 2) {
      HIV[t, y] ~ poisson(exp_HIV_dx[t, y]);
      AIDS[t, y] ~ poisson(exp_AIDS_dx[t, y]);
      CD4[t, y,  : ] ~ multinomial(exp_p[t, y,  : ]);
    }
  }
  ABROAD ~ poisson(exp_ABROAD_dx);
}
generated quantities {
  // Log-lik calculations to get WAIC / LOO/ DIC
  matrix[nquar, 2] log_lik;
  vector[nquar] log_lik_abroad;
  for (t in 1 : nquar) {
    log_lik_abroad[t] = poisson_lpmf(ABROAD[t] | exp_ABROAD_dx[t]);
    
    for (y in 1 : 2) {
      log_lik[t, y] = poisson_lpmf(HIV[t, y] | exp_HIV_dx[t, y])
                      + poisson_lpmf(AIDS[t, y] | exp_AIDS_dx[t, y])
                      + multinomial_lpmf(CD4[t, y,  : ] | exp_p[t, y,  : ]);
    }
  }
}
