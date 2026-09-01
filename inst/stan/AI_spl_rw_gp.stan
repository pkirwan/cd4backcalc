// Age-independent model: spline/RW/GP infections, spline or AR(1) diagnosis

functions {
  // Build transition matrix for quarter t
  matrix tmat_fct(int t, vector q, matrix d) {
    matrix[9, 9] Lambda = rep_matrix(0, 9, 9);
    
    // Diagnosed states are absorbing
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
  int<lower=1> nquar; // quarters
  array[nquar] int<lower=0> HIV;
  array[nquar] int<lower=0> AIDS;
  array[nquar, 4] int<lower=0> CD4;
  vector<lower=0, upper=1>[4] q; // progression probabilities
  vector<lower=0>[4] init_prev;
  
  // Model type indicators
  int<lower=1, upper=3> inf_model; // 1 = spline, 2 = random walk, 3 = GP
  int<lower=1, upper=2> diag_model; // 1 = spline, 2 = AR(1)
  
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
  real alpha_diag_1;
  real alpha_diag_2;
  real alpha_diag_3;
  real alpha_diag_4;
  
  // AR(1) model
  vector[nquar] delta_raw_1;
  vector[nquar] delta_raw_2;
  vector[nquar] delta_raw_3;
  vector[nquar] delta_raw_4;
  vector<lower=0>[4] vardelta; // variances of diagnoses
  // AR(1) coefficient for diagnosis-probability walks (diag_model == 2); (0,1) is mean-reverting
  real<lower=0, upper=1> rho_diag;
}
transformed parameters {
  // Infection process
  real<lower=0> sigma_inf = exp(log_sigma_inf);
  vector[nquar] lambda = X_inf * beta_inf;
  vector<lower=0>[nquar] h;
  
  // Diagnosis probabilities
  vector[nquar] delta_1;
  vector[nquar] delta_2;
  vector[nquar] delta_3;
  vector[nquar] delta_4;
  matrix<lower=0, upper=1>[4, nquar] d;
  
  // Expected observations
  row_vector<lower=0>[nquar] exp_HIV_dx;
  row_vector<lower=0>[nquar] exp_AIDS_dx;
  array[nquar] simplex[4] exp_p = rep_array(rep_vector(0.25, 4), nquar);
  
  // Temporary variables
  real<lower=0> inv_rho = inv(rho);
  
  // Initialize h based on infection model
  if (inf_model == 1) {
    h = exp(lambda);
  } else {
    h = exp(gamma);
  }
  
  // Diagnosis process
  if (diag_model == 1) {
    // Spline model for diagnoses
    delta_1 = alpha_diag_1 + (X_diag * beta_diag_1);
    delta_2 = alpha_diag_2 + (X_diag * beta_diag_2);
    delta_3 = alpha_diag_3 + (X_diag * beta_diag_3);
    delta_4 = alpha_diag_4 + (X_diag * beta_diag_4);
  } else if (diag_model == 2) {
    // Mean-reverting AR(1), non-centred: t=1 draws from stationary dist, t>1 mean-reverts to alpha_diag_k
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
  }
  
  // Transform diag probs using logit transformation
  for (t in 1 : nquar) {
    d[ : , t] = inv_logit([delta_1[t], delta_2[t], delta_3[t], delta_4[t]]');
  }
  
  {
  // Local: not saved as draws; only exp_HIV_dx/exp_AIDS_dx/exp_p needed downstream
  array[nquar] vector[9] cum_diags = rep_array(rep_vector(0, 9),
                                                        nquar);
  row_vector[nquar] incidence_state_5;
  row_vector[nquar] incidence_state_6;
  row_vector[nquar] incidence_state_7;
  row_vector[nquar] incidence_state_8;
  row_vector[nquar] incidence_state_9;
  array[nquar] vector[4] prev_mat;
  
  // Initial time period (t=1, y=1)
  vector[9] tmpvec = append_row(init_prev, rep_vector(0, 5));
  matrix[9, 9] Lambda1 = tmat_fct(1, q, d);
  cum_diags[1] = Lambda1' * tmpvec;
  cum_diags[1, 1] += h[1];
  
  // Evolution over time
  for (t in 2 : nquar) {
    matrix[9, 9] Lambda_t = tmat_fct(t, q, d);
    cum_diags[t] = Lambda_t' * cum_diags[t - 1];
    cum_diags[t, 1] += h[t];
  }
  
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
}
model {
  // Infection priors
  if (inf_model == 1) {
    beta_inf[1] ~ normal(5, 1.5);
    beta_inf[2 : npars_inf] ~ normal(0, sigma_inf);
    // Loose prior avoids suppressing genuine peaks, as in AI_spl.stan
    log_sigma_inf ~ normal(-1, 1);
    
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
    // Sigma/L are local: not saved as (nquar x nquar) draws, only needed here
    matrix[nquar, nquar] Sigma = gp_exp_quad_cov(x, eta, inv_rho);
    for (i in 1 : nquar) {
      Sigma[i, i] += sigma_sq; // ensures numerical PD
    }
    matrix[nquar, nquar] L = cholesky_decompose(Sigma);
    
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
    // Empirical range midpoints (as elsewhere): State1 2.5-10%, State2 2.5-11%, State3 5-25%, State4 20-60%
    alpha_diag_1 ~ normal(-2.9, 0.5);
    alpha_diag_2 ~ normal(-2.75, 0.5);
    alpha_diag_3 ~ normal(-2.0, 0.5);
    alpha_diag_4 ~ normal(-0.5, 0.5);
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
    
    // Weak priors for inactive AR(1) diagnosis parameters
    delta_raw_1 ~ std_normal();
    delta_raw_2 ~ std_normal();
    delta_raw_3 ~ std_normal();
    delta_raw_4 ~ std_normal();
    vardelta ~ gamma(2, 2);
    rho_diag ~ beta(1, 1);
  } else if (diag_model == 2) {
    delta_raw_1 ~ std_normal();
    delta_raw_2 ~ std_normal();
    delta_raw_3 ~ std_normal();
    delta_raw_4 ~ std_normal();
    vardelta ~ gamma(1, 32);
    // Mode ~0.86: persistent but rules out the random-walk boundary at rho=1
    rho_diag ~ beta(7, 2);
    
    // Priors for AR(1) mean-reversion targets
    alpha_diag_1 ~ normal(-2.9, 0.5);
    alpha_diag_2 ~ normal(-2.75, 0.5);
    alpha_diag_3 ~ normal(-2.0, 0.5);
    alpha_diag_4 ~ normal(-0.5, 0.5);
    
    // Weak priors for inactive spline diagnosis parameters
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
