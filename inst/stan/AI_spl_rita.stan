// Age-independent spline model with RITA state

functions {
  // Build transition matrix for quarter t
  matrix tmat_fct(int t, vector q, matrix d) {
    matrix[13, 13] Lambda = rep_matrix(0, 13, 13);
    
    // Diagnosed states are absorbing
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
    real q1s1 = q[1] * s1;
    real nq1s1 = (1 - q[1]) * s1;
    real q2s1 = q[2] * s1;
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
  int<lower=1> nquar; // quarters
  array[nquar] int<lower=0> HIV;
  array[nquar] int<lower=0> AIDS;
  array[nquar, 5] int<lower=0> CD4;
  vector<lower=0, upper=1>[4] q; // progression probabilities
  vector<lower=0>[7] init_prev;
  
  // Splines
  int<lower=1> npars_inf;
  matrix[nquar, npars_inf] X_inf;
}
parameters {
  // Infection parameters
  vector[npars_inf] beta_inf;
  real log_sigma_inf;
  
  // False recency rate
  real logit_frr;
  
  // Diagnosis probability parameters (AR(1))
  vector[nquar] delta_raw_1;
  vector[nquar] delta_raw_2;
  vector[nquar] delta_raw_3;
  vector[nquar] delta_raw_4;
  vector[nquar] delta_raw_5;
  vector[5] log_vardelta_raw;
  real mu_log_var;
  real<lower=0> sigma_log_var;
  real alpha_diag_1;
  real alpha_diag_2;
  real alpha_diag_3;
  real alpha_diag_4;
  real alpha_diag_5;
  
  // AR(1) coefficient for diagnosis-probability walks; (0,1) is mean-reverting
  real<lower=0, upper=1> rho_diag;
}
transformed parameters {
  // False recency rate
  real<lower=0, upper=1> frr = inv_logit(logit_frr);
  
  // Infection process
  real<lower=0> sigma_inf = exp(log_sigma_inf);
  vector[nquar] lambda = X_inf * beta_inf;
  vector<lower=0>[nquar] h = exp(lambda);
  
  // Step variances on natural scale (non-centred)
  vector[5] log_vardelta = mu_log_var + sigma_log_var * log_vardelta_raw;
  vector<lower=0>[5] vardelta = exp(log_vardelta);
  // Diagnosis probabilities via AR(1)
  vector[nquar] delta_1;
  vector[nquar] delta_2;
  vector[nquar] delta_3;
  vector[nquar] delta_4;
  vector[nquar] delta_5;
  // AR(1), non-centred: t=1 draws from stationary dist, t>1 mean-reverts to alpha_diag_k
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
  matrix<lower=0, upper=1>[5, nquar] d;
  
  // Transform diag probs using logit transformation
  for (t in 1 : nquar) {
    d[ : , t] = inv_logit(
                          [delta_1[t], delta_2[t], delta_3[t], delta_4[t],
                           delta_5[t]]');
  }
  
  // Expected observations
  row_vector<lower=0>[nquar] exp_HIV_dx;
  row_vector<lower=0>[nquar] exp_AIDS_dx;
  array[nquar] simplex[5] exp_p = rep_array(rep_vector(0.20, 5), nquar);
  
  {
  // Local: not saved as draws; only exp_HIV_dx/exp_AIDS_dx/exp_p needed downstream
  array[nquar] vector[13] cum_diags = rep_array(rep_vector(0, 13),
                                                         nquar);
  row_vector[nquar] incidence_state_8;
  row_vector[nquar] incidence_state_9;
  row_vector[nquar] incidence_state_10;
  row_vector[nquar] incidence_state_11;
  row_vector[nquar] incidence_state_12;
  row_vector[nquar] incidence_state_13;
  array[nquar] vector[7] prev_mat;
  
  // Initial time period (t=1, y=1)
  vector[13] init_vec = rep_vector(0, 13);
  init_vec[1] = h[1];
  
  vector[13] tmpvec = append_row(init_prev, rep_vector(0, 6));
  
  matrix[13, 13] Lambda1 = tmat_fct(1, q, d);
  cum_diags[1] = Lambda1' * tmpvec + init_vec;
  
  // Evolution over time
  for (t in 2 : nquar) {
    matrix[13, 13] Lambda_t = tmat_fct(t, q, d);
    cum_diags[t] = Lambda_t' * cum_diags[t - 1];
    cum_diags[t, 1] += h[t];
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
    incidence_state_8[t] = cum_diags[t, 8] - cum_diags[t - 1, 8];
    incidence_state_9[t] = cum_diags[t, 9] - cum_diags[t - 1, 9];
    incidence_state_10[t] = cum_diags[t, 10] - cum_diags[t - 1, 10];
    incidence_state_11[t] = cum_diags[t, 11] - cum_diags[t - 1, 11];
    incidence_state_12[t] = cum_diags[t, 12] - cum_diags[t - 1, 12];
    incidence_state_13[t] = cum_diags[t, 13] - cum_diags[t - 1, 13];
  }
  
  // Prevalence matrix
  for (t in 1 : nquar) {
    prev_mat[t] = cum_diags[t, 1 : 7];
  }
  
  // Expected diagnoses
  exp_HIV_dx = incidence_state_9 + incidence_state_10 + incidence_state_11
               + incidence_state_12 + incidence_state_13;
  exp_AIDS_dx = incidence_state_8;
  
  // Including CD4 / RITA counts, adjusted for false recency rate
  for (t in 1 : nquar) {
    if (exp_HIV_dx[t] > 0) {
      real p1 = incidence_state_9[t] / exp_HIV_dx[t];
      exp_p[t, 1] = p1 + frr * (1 - p1);
      exp_p[t, 2] = (1 - frr) * incidence_state_10[t] / exp_HIV_dx[t];
      exp_p[t, 3] = (1 - frr) * incidence_state_11[t] / exp_HIV_dx[t];
      exp_p[t, 4] = (1 - frr) * incidence_state_12[t] / exp_HIV_dx[t];
      exp_p[t, 5] = (1 - frr) * incidence_state_13[t] / exp_HIV_dx[t];
    }
  }
  }
}
model {
  // Infection priors
  // Note: unlike the AD model which penalises curvature via a precision matrix,
  // this model shrinks coefficient magnitude (ridge); normal(-1, 1) loosens the
  // prior so large peaks aren't suppressed, matching the AD model's permissiveness.
  beta_inf[1] ~ normal(5, 1.5);
  beta_inf[2 : npars_inf] ~ normal(0, sigma_inf);
  log_sigma_inf ~ normal(-1, 1);
  
  // Diagnosis probability priors
  // Initial level intercepts, centred on empirical range midpoints:
  // State 1 (recent/RITA+):   5–30%  → logit ~ -1.75
  // State 2 (CD4 ≥ 500):     2.5–10%  → logit ~ -2.90
  // State 3 (CD4 350–500):   2.5–11%  → logit ~ -2.75
  // State 4 (CD4 200–349):   5–25%    → logit ~ -2.00
  // State 5 (CD4 < 200):     20–60%   → logit ~ -0.50
  alpha_diag_1 ~ normal(-1.75, 0.5);
  alpha_diag_2 ~ normal(-2.9, 0.5);
  alpha_diag_3 ~ normal(-2.75, 0.5);
  alpha_diag_4 ~ normal(-2.0, 0.5);
  alpha_diag_5 ~ normal(-0.5, 0.5);
  delta_raw_1 ~ std_normal();
  delta_raw_2 ~ std_normal();
  delta_raw_3 ~ std_normal();
  delta_raw_4 ~ std_normal();
  delta_raw_5 ~ std_normal();
  // Step SD ~0.13 logit/quarter (stationary ~0.26); avoids heavy-tailed swings
  mu_log_var ~ normal(-4.0, 0.5);
  sigma_log_var ~ normal(0, 0.5);
  log_vardelta_raw ~ std_normal();
  
  // Mode ~0.86: persistent but rules out the random-walk boundary at rho=1
  rho_diag ~ beta(7, 2);
  
  // False recency rate prior centred on logit(0.019)
  logit_frr ~ normal(-3.94, 0.32);
  
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
