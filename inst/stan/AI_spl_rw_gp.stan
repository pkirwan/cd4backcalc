// Age-independent model with spline, random walk or GP for infection
// and spline or random walk for diagnosis probabilities

functions {
  // Function to create transition matrices
  matrix tmat_fct(int t, vector q, matrix d) {
    matrix[9, 9] Lambda = rep_matrix(0, 9, 9);
    
    // Setting diagonal entries to one
    for (i in 5 : 9) {
      Lambda[i, i] = 1;
    }
    
    Lambda[1, 1] = (1 - q[1]) * (1 - d[1, t]);
    Lambda[1, 2] = q[1] * (1 - d[1, t]);
    Lambda[1, 6] = d[1, t];
    
    Lambda[2, 2] = (1 - q[2]) * (1 - d[2, t]);
    Lambda[2, 3] = q[2] * (1 - d[2, t]);
    Lambda[2, 7] = d[2, t];
    
    Lambda[3, 3] = (1 - q[3]) * (1 - d[3, t]);
    Lambda[3, 4] = q[3] * (1 - d[3, t]);
    Lambda[3, 8] = d[3, t];
    
    Lambda[4, 4] = (1 - q[4]) * (1 - d[4, t]);
    Lambda[4, 5] = q[4] * (1 - d[4, t]);
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
  
  // Model type indicators
  int<lower=1, upper=3> inf_model; // 1 = spline, 2 = random walk, 3 = GP
  int<lower=1, upper=2> diag_model; // 1 = spline, 2 = random walk
  
  // Splines
  int<lower=1> npars_inf;
  int<lower=1> npars_diag;
  matrix[nquar, npars_inf] X_inf;
  matrix[nquar, npars_diag] X_diag;
  
  // GP parameters
  array[nquar] real x; // inputs of GP - 1:nquar
  real<lower=0> sigma_sq; // Arbitrary noise of GP
}
transformed data {
  vector[8] zeroes = rep_vector(0, 8);
  vector[nquar] mu = rep_vector(0, nquar);
}
parameters {
  // Infection parameters
  // Spline model
  vector[npars_inf] beta_inf;
  real log_sigma_inf;
  // Random walk model
  vector[nquar] gamma;
  real<lower=0> var1;
  // GP model
  real<lower=0> eta;
  real<lower=0> rho;
  
  // Diagnosis probability parameters
  // Spline model
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
  
  // Random walk model
  vector[nquar] delta_raw_1;
  vector[nquar] delta_raw_2;
  vector[nquar] delta_raw_3;
  vector[nquar] delta_raw_4;
  vector<lower=0>[4] vardelta; // variances of diagnoses
}
transformed parameters {
  // Infection process
  real<lower=0> sigma_inf = exp(log_sigma_inf);
  vector[nquar] lambda = X_inf * beta_inf;
  vector<lower=0>[nquar] h;
  
  // GP parameters
  matrix[nquar, nquar] Sigma;
  matrix[nquar, nquar] L;
  
  // Diagnosis probabilities
  vector[nquar] delta_1;
  vector[nquar] delta_2;
  vector[nquar] delta_3;
  vector[nquar] delta_4;
  matrix<lower=0, upper=1>[4, nquar] d;
  
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
  
  // Temporary variables
  vector[9] tmpvec;
  vector[9] init_vec;
  real<lower=0> inv_rho = inv(rho);
  
  // Initialize h based on infection model
  if (inf_model == 1) {
    h = exp(lambda);
  } else {
    h = exp(gamma);
  }
  
  if (inf_model == 3) {
    Sigma = gp_exp_quad_cov(x, eta, inv_rho);
    for (i in 1 : nquar) {
      Sigma[i, i] = Sigma[i, i] + sigma_sq; // ensures numerical PD
    }
    L = cholesky_decompose(Sigma);
  }
  
  // Initializing delta based on diagnosis model
  if (diag_model == 1) {
    // Spline model for diagnoses
    delta_1 = alpha_diag_1 + 0.5 * (X_diag * beta_diag_1);
    delta_2 = alpha_diag_2 + 0.5 * (X_diag * beta_diag_2);
    delta_3 = alpha_diag_3 + 0.5 * (X_diag * beta_diag_3);
    delta_4 = alpha_diag_4 + 0.7 * (X_diag * beta_diag_4);
  } else if (diag_model == 2) {
    // Random walk for diagnoses
    delta_1[1] = -3.2 + 0.2 * delta_raw_1[1];
    delta_2[1] = -3.2 + 0.2 * delta_raw_2[1];
    delta_3[1] = -3.0 + 0.2 * delta_raw_3[1];
    delta_4[1] = -2.5 + 0.3 * delta_raw_4[1];
    
    for (t in 2 : nquar) {
      delta_1[t] = delta_1[t - 1] + sqrt(vardelta[1]) * delta_raw_1[t];
      delta_2[t] = delta_2[t - 1] + sqrt(vardelta[2]) * delta_raw_2[t];
      delta_3[t] = delta_3[t - 1] + sqrt(vardelta[3]) * delta_raw_3[t];
      delta_4[t] = delta_4[t - 1] + sqrt(vardelta[4]) * delta_raw_4[t];
    }
  }
  
  // Transform diag probs using logit transformation
  for (t in 1 : nquar) {
    d[ : , t] = inv_logit([delta_1[t], delta_2[t], delta_3[t], delta_4[t]]');
  }
  
  // Initial time period (t=1, y=1)
  init_vec = rep_vector(0, 9);
  init_vec[1] = h[1];
  
  tmpvec = append_row(init_prev, zeroes[1 : 5]);
  
  matrix[9, 9] Lambda1 = tmat_fct(1, q, d);
  cum_diags[1] = Lambda1' * tmpvec + init_vec;
  
  // Evolution over time
  for (t in 2 : nquar) {
    vector[9] inflow = rep_vector(0, 9);
    inflow[1] = h[t];
    matrix[9, 9] Lambda_t = tmat_fct(t, q, d);
    cum_diags[t] = Lambda_t' * cum_diags[t - 1] + inflow;
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
  if (inf_model == 1) {
    beta_inf[1] ~ normal(5, 1.5); // Vaguely informative priors: corresponds to quarter inf in (exp(2),exp(8)) = (8,4390)
    beta_inf[2 : npars_inf] ~ normal(0, sigma_inf);
    log_sigma_inf ~ normal(-2, 0.5);
    
    // weak priors for inactive infection params
    gamma ~ normal(0, 5);
    var1 ~ gamma(2, 2);
    eta ~ normal(0, 2);
    rho ~ lognormal(0, 1);
  } else if (inf_model == 2) {
    gamma[1] ~ normal(5.521461, 0.2938935);
    var1 ~ gamma(1, 32);
    tail(gamma, nquar - 1) ~ normal(head(gamma, nquar - 1), sqrt(var1));
    
    // weak priors for inactive infection params
    beta_inf ~ normal(0, 5);
    log_sigma_inf ~ normal(-2, 1);
    eta ~ normal(0, 2);
    rho ~ lognormal(0, 1);
  } else if (inf_model == 3) {
    rho ~ student_t(4, 0, 1);
    eta ~ normal(4, 1);
    gamma ~ multi_normal_cholesky(mu, L);
    
    // weak priors for inactive infection params
    beta_inf ~ normal(0, 5);
    log_sigma_inf ~ normal(-2, 1);
    var1 ~ gamma(2, 2);
  }
  
  // Diagnosis probability priors
  if (diag_model == 1) {
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
    sigma_diag_1 ~ gamma(2, 64);
    sigma_diag_2 ~ gamma(2, 64);
    sigma_diag_3 ~ gamma(2, 64);
    sigma_diag_4 ~ gamma(2, 64);
    
    // weak priors for inactive RW diag params
    delta_raw_1 ~ normal(0, 1);
    delta_raw_2 ~ normal(0, 1);
    delta_raw_3 ~ normal(0, 1);
    delta_raw_4 ~ normal(0, 1);
    vardelta ~ gamma(2, 2);
  } else if (diag_model == 2) {
    delta_raw_1 ~ normal(0, 1);
    delta_raw_2 ~ normal(0, 1);
    delta_raw_3 ~ normal(0, 1);
    delta_raw_4 ~ normal(0, 1);
    vardelta ~ gamma(1, 32);
    
    // weak priors for inactive spline diag params
    alpha_diag_1 ~ normal(0, 2);
    alpha_diag_2 ~ normal(0, 2);
    alpha_diag_3 ~ normal(0, 2);
    alpha_diag_4 ~ normal(0, 2);
    beta_diag_1 ~ normal(0, 2);
    beta_diag_2 ~ normal(0, 2);
    beta_diag_3 ~ normal(0, 2);
    beta_diag_4 ~ normal(0, 2);
    sigma_diag_1 ~ gamma(2, 2);
    sigma_diag_2 ~ gamma(2, 2);
    sigma_diag_3 ~ gamma(2, 2);
    sigma_diag_4 ~ gamma(2, 2);
  }
  
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
