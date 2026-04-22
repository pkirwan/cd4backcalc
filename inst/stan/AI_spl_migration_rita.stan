// Age-independent spline model including migration and RITA state

functions {
  // Function to create transition matrices
  matrix tmat_fct(int t, vector q, matrix d) {
    matrix[13, 13] Lambda = rep_matrix(0, 13, 13);
    
    // Setting diagonal entries to one
    for (i in 8 : 13) {
      Lambda[i, i] = 1;
    }
    
    // Cache (1 - d[i, t]) and repeated products to avoid redundant computation
    real s1 = 1 - d[1, t];
    real s2 = 1 - d[2, t];
    real s3 = 1 - d[3, t];
    real s4 = 1 - d[4, t];
    real s5 = 1 - d[5, t];
    // States 1,2,3 all share d[1,t]; cache the two products used across rows
    real q1s1  = q[1] * s1;
    real nq1s1 = (1 - q[1]) * s1;
    real q2s1  = q[2] * s1;
    real nq2s1 = (1 - q[2]) * s1;
    
    Lambda[1, 1] = 0;
    Lambda[1, 2] = nq1s1;
    Lambda[1, 3] = q1s1;
    Lambda[1, 9] = d[1, t];
    
    Lambda[2, 2] = 0;
    Lambda[2, 4] = nq1s1;
    Lambda[2, 5] = q1s1;
    Lambda[2, 9] = d[1, t];
    
    Lambda[3, 3] = 0;
    Lambda[3, 5] = nq2s1;
    Lambda[3, 6] = q2s1;
    Lambda[3, 9] = d[1, t];
    
    Lambda[4, 4] = (1 - q[1]) * s2;
    Lambda[4, 5] = q[1] * s2;
    Lambda[4, 10] = d[2, t];
    
    Lambda[5, 5] = (1 - q[2]) * s3;
    Lambda[5, 6] = q[2] * s3;
    Lambda[5, 11] = d[3, t];
    
    Lambda[6, 6] = (1 - q[3]) * s4;
    Lambda[6, 7] = q[3] * s4;
    Lambda[6, 12] = d[4, t];
    
    Lambda[7, 7] = (1 - q[4]) * s5;
    Lambda[7, 8] = q[4] * s5;
    Lambda[7, 13] = d[5, t];
    
    return Lambda;
  }
}
data {
  int<lower=1> nquar; // number of quarters considered
  array[nquar, 2] int<lower=0> HIV; // new HIV diagnoses by quarter and native/migrant
  array[nquar, 2] int<lower=0> AIDS; // new AIDS diagnoses by quarter and native/migrant
  array[nquar, 2, 5] int<lower=0> CD4; // CD4 distribution at HIV diagnosis (5 states)
  array[nquar] int<lower=0> ABROAD; // prior diagnoses from abroad (detected)
  vector<lower=0, upper=1>[4] q; // fixed progression probabilities
  vector<lower=0>[7] init_prev; // initial prevalence (latent states)
  vector<lower=0>[7] init_prev_abroad; // initial prevalence (latent states) for migrant group
  
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
  vector[npars_inf] beta_inf;
  real log_sigma_inf;
  
  // Migration parameters
  vector[npars_mig] beta_mig;
  real log_sigma_mig;
  
  // False recency rate
  real logit_frr;
  
  // Diagnosis probability parameters
  vector[npars_diag] beta_diag_1;
  vector[npars_diag] beta_diag_2;
  vector[npars_diag] beta_diag_3;
  vector[npars_diag] beta_diag_4;
  vector[npars_diag] beta_diag_5;
  real<lower=0> sigma_diag_1;
  real<lower=0> sigma_diag_2;
  real<lower=0> sigma_diag_3;
  real<lower=0> sigma_diag_4;
  real<lower=0> sigma_diag_5;
  real alpha_diag_1;
  real alpha_diag_2;
  real alpha_diag_3;
  real alpha_diag_4;
  real alpha_diag_5;
  
  // Proportion UK-acquired parameters
  vector[npars_cob] beta_cob;
  real<lower=0> sigma_cob;
  real alpha_cob;
  real<lower=-1, upper=1> rho_eta;
  
  // Probability of migration into each compartment (4 latent CD4 states + abroad-diagnosed)
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
  real alpha_detect;
  vector[npars_mp] beta_detect;
  real<lower=0> sigma_detect;
}
transformed parameters {
  // False recency rate
  real<lower=0, upper=1> frr = inv_logit(logit_frr);
  
  // Infection process
  real<lower=0> sigma_inf = exp(log_sigma_inf);
  vector[nquar] lambda = X_inf * beta_inf;
  vector<lower=0>[nquar] h = exp(lambda);
  
  // Migration process
  real<lower=0> sigma_mig = exp(log_sigma_mig);
  vector[nquar] gamma = X_mig * beta_mig;
  vector<lower=0>[nquar] omega = exp(gamma);
  
  // Diagnosis probabilities
  vector[nquar] delta_1 = alpha_diag_1 + 0.5 * (X_diag * beta_diag_1);
  vector[nquar] delta_2 = alpha_diag_2 + 0.5 * (X_diag * beta_diag_2);
  vector[nquar] delta_3 = alpha_diag_3 + 0.5 * (X_diag * beta_diag_3);
  vector[nquar] delta_4 = alpha_diag_4 + 0.7 * (X_diag * beta_diag_4);
  vector[nquar] delta_5 = alpha_diag_5 + 0.7 * (X_diag * beta_diag_5);
  matrix<lower=0, upper=1>[5, nquar] d;
  
  for (t in 1 : nquar) {
    d[ : , t] = inv_logit(
                          [delta_1[t], delta_2[t], delta_3[t], delta_4[t],
                           delta_5[t]]');
  }
  
  // Proportion UK-acquired (log-scale)
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
  array[nquar, 2] vector<lower=0>[13] cum_diags = rep_array(
                                                            rep_vector
                                                            (0, 13), nquar,
                                                            2);
  
  // Expected observations
  matrix<lower=0>[nquar, 2] exp_HIV_dx;
  matrix<lower=0>[nquar, 2] exp_AIDS_dx;
  vector<lower=0>[nquar] exp_ABROAD_dx;
  array[nquar, 2] simplex[5] exp_p = rep_array(rep_vector(0.20, 5), nquar, 2);
  
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
  
  // Expected incidence by diagnosed states (AIDS + 5 HIV CD4 states)
  matrix<lower=0>[nquar, 2] incidence_state_8 = rep_matrix(0, nquar, 2);
  matrix<lower=0>[nquar, 2] incidence_state_9 = rep_matrix(0, nquar, 2);
  matrix<lower=0>[nquar, 2] incidence_state_10 = rep_matrix(0, nquar, 2);
  matrix<lower=0>[nquar, 2] incidence_state_11 = rep_matrix(0, nquar, 2);
  matrix<lower=0>[nquar, 2] incidence_state_12 = rep_matrix(0, nquar, 2);
  matrix<lower=0>[nquar, 2] incidence_state_13 = rep_matrix(0, nquar, 2);
  array[nquar, 2] vector[7] prev_mat;
  
  // Initial time period (t=1)
  matrix[13, 13] Lambda1 = tmat_fct(1, q, d);
  
  cum_diags[1, 1,  : ] = Lambda1' * append_row(init_prev, rep_vector(0, 6));
  cum_diags[1, 1, 1] += h[1];
  
  // Initial migrant group
  p[1] = exp(eta[1]);
  cum_diags[1, 2,  : ] = Lambda1' * append_row(init_prev_abroad, rep_vector(0, 6));
  cum_diags[1, 2, 1] += h[1] * p[1];
  cum_diags[1, 2, 4 : 7] += omega[1] * m[1 : 4, 1];
  
  // Evolution over time
  for (t in 2 : nquar) {
    matrix[13, 13] Lambda_t = tmat_fct(t, q, d);
    
    // Native group
    cum_diags[t, 1,  : ] = Lambda_t' * cum_diags[t - 1, 1,  : ];
    cum_diags[t, 1, 1] += h[t];
    
    // Migrant group
    p[t] = exp(eta[t]);
    cum_diags[t, 2,  : ] = Lambda_t' * cum_diags[t - 1, 2,  : ];
    cum_diags[t, 2, 1] += h[t] * p[t];
    cum_diags[t, 2, 4 : 7] += omega[t] * m[1 : 4, t];
  }
  
  // Calculate incidences
  for (y in 1 : 2) {
    incidence_state_8[1, y] = cum_diags[1, y, 8];
    incidence_state_9[1, y] = cum_diags[1, y, 9];
    incidence_state_10[1, y] = cum_diags[1, y, 10];
    incidence_state_11[1, y] = cum_diags[1, y, 11];
    incidence_state_12[1, y] = cum_diags[1, y, 12];
    incidence_state_13[1, y] = cum_diags[1, y, 13];
    
    for (t in 2 : nquar) {
      incidence_state_8[t, y] = cum_diags[t, y, 8] - cum_diags[t - 1, y, 8];
      incidence_state_9[t, y] = cum_diags[t, y, 9] - cum_diags[t - 1, y, 9];
      incidence_state_10[t, y] = cum_diags[t, y, 10] - cum_diags[t - 1, y, 10];
      incidence_state_11[t, y] = cum_diags[t, y, 11] - cum_diags[t - 1, y, 11];
      incidence_state_12[t, y] = cum_diags[t, y, 12] - cum_diags[t - 1, y, 12];
      incidence_state_13[t, y] = cum_diags[t, y, 13] - cum_diags[t - 1, y, 13];
    }
  }
  
  // Prevalence matrix (latent states)
  for (t in 1 : nquar) {
    for (y in 1 : 2) {
      prev_mat[t, y] = cum_diags[t, y, 1 : 7];
    }
  }
  
  // Expected diagnoses
  exp_HIV_dx = incidence_state_9 + incidence_state_10 + incidence_state_11
               + incidence_state_12 + incidence_state_13 + misclassified_dx;
  exp_AIDS_dx = incidence_state_8;
  
  // Including CD4 / RITA counts, adjusted for false recency rate
  for (t in 1 : nquar) {
    // Native population (y=1)
    if (exp_HIV_dx[t, 1] > 0) {
      real p1 = incidence_state_9[t, 1] / exp_HIV_dx[t, 1];
      exp_p[t, 1, 1] = p1 + frr * (1 - p1);
      exp_p[t, 1, 2] = (1 - frr) * incidence_state_10[t, 1] / exp_HIV_dx[t, 1];
      exp_p[t, 1, 3] = (1 - frr) * incidence_state_11[t, 1] / exp_HIV_dx[t, 1];
      exp_p[t, 1, 4] = (1 - frr) * incidence_state_12[t, 1] / exp_HIV_dx[t, 1];
      exp_p[t, 1, 5] = (1 - frr) * incidence_state_13[t, 1] / exp_HIV_dx[t, 1];
    }
    
    // Migrant population (y=2)
    real denom = exp_HIV_dx[t, 2] - misclassified_dx[t, 2];
    if (denom > 0) {
      real p1 = incidence_state_9[t, 2] / denom;
      exp_p[t, 2, 1] = p1 + frr * (1 - p1);
      exp_p[t, 2, 2] = (1 - frr) * incidence_state_10[t, 2] / denom;
      exp_p[t, 2, 3] = (1 - frr) * incidence_state_11[t, 2] / denom;
      exp_p[t, 2, 4] = (1 - frr) * incidence_state_12[t, 2] / denom;
      exp_p[t, 2, 5] = (1 - frr) * incidence_state_13[t, 2] / denom;
    }
  }
}
model {
  // PRIORS
  
  // Infection priors
  beta_inf[1] ~ normal(5, 1.5);
  beta_inf[2 : npars_inf] ~ normal(0, sigma_inf);
  log_sigma_inf ~ normal(-2, 0.5);
  
  // Migration priors
  beta_mig[1] ~ normal(4, 1);
  beta_mig[2 : npars_mig] ~ normal(0, sigma_mig);
  log_sigma_mig ~ normal(-2, 0.5);
  
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
  sigma_diag_1 ~ gamma(2, 8);
  sigma_diag_2 ~ gamma(2, 8);
  sigma_diag_3 ~ gamma(2, 8);
  sigma_diag_4 ~ gamma(2, 8);
  sigma_diag_5 ~ gamma(2, 8);
  
  // False recency rate prior centred on logit(0.019)
  logit_frr ~ normal(-3.94, 0.32);
  
  // Proportion UK-acquired priors
  alpha_cob ~ std_normal();
  beta_cob[1] ~ std_normal();
  beta_cob[2 : npars_cob] ~ normal(0, sigma_cob);
  sigma_cob ~ gamma(2, 64);
  
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
  alpha_detect ~ normal(0, 0.2);
  beta_detect[1] ~ normal(0, 0.2);
  beta_detect[2 : npars_mp] ~ normal(0, sigma_detect);
  sigma_detect ~ gamma(8, 64);
  
  // Migration distribution priors
  alpha_mp[1] ~ normal(-2, 1);
  alpha_mp[2 : 4] ~ std_normal();
  beta_mp_1[1] ~ std_normal();
  beta_mp_2[1] ~ std_normal();
  beta_mp_3[1] ~ std_normal();
  beta_mp_5[1] ~ std_normal();
  beta_mp_1[2 : npars_mp] ~ normal(0, sigma_mp_1);
  beta_mp_2[2 : npars_mp] ~ normal(0, sigma_mp_2);
  beta_mp_3[2 : npars_mp] ~ normal(0, sigma_mp_3);
  beta_mp_5[2 : npars_mp] ~ normal(0, sigma_mp_5);
  sigma_mp_1 ~ gamma(2, 64);
  sigma_mp_2 ~ gamma(2, 64);
  sigma_mp_3 ~ gamma(2, 64);
  sigma_mp_5 ~ gamma(2, 64);
  
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
