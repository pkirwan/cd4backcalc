// Age-independent spline model with migration
// Separate undiagnosed/transfer flows and a fixed CD4-composition simplex,
// matching AD_ptens_migration.stan.

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
  // y = 1: UK-born, y = 2: born abroad
  array[nquar, 2] int<lower=0> HIV;
  array[nquar, 2] int<lower=0> AIDS;
  array[nquar, 2, 4] int<lower=0> CD4;
  array[nquar] int<lower=0> ABROAD; // previous diagnoses abroad
  int<lower=1> nyear;
  // Numerator and auxiliary-observed risk set for recent-arrival membership
  array[nyear, 4] int<lower=0> recent_arrival_cd4;
  array[nyear, 4] int<lower=0> recent_arrival_cd4_total;
  vector<lower=0, upper=1>[4] q; // progression probabilities
  vector<lower=0>[4] init_prev;
  vector<lower=0>[4] init_prev_abroad;
  
  // Splines
  int<lower=1> npars_inf;
  int<lower=1> npars_mig;
  int<lower=1> npars_cob;
  matrix[nquar, npars_inf] X_inf;
  matrix[nquar, npars_mig] X_mig;
  matrix[nquar, npars_cob] X_cob;
  
  // External evidence for the log born-abroad:UK-born UK-acquired ratio
  vector[nquar] eta_mean;
  real<lower=0> eta_sd;
  
  // External data-linkage estimates of detection probability (logit scale, can be sparse)
  int<lower=0> n_detect_obs;
  array[n_detect_obs] int<lower=1, upper=nquar> detect_obs_index;
  vector[n_detect_obs] detect_logit_obs;
  real<lower=0> detect_obs_sd;
  
  // Quarters where abroad-diagnosis ascertainment steps up/down (e.g. reporting change)
  int<lower=0> n_detect_steps;
  array[n_detect_steps] int<lower=2, upper=nquar> detect_step_quarter;

  // Independent rolling-annual migration index; sd=0 disables the likelihood
  vector[nquar] gamma_obs;
  real<lower=0> gamma_obs_sd;
  real<lower=0> gamma_obs_drift_sd;
  real<lower=0> gamma_obs_slope_sd;
}
transformed data {
  vector[nquar] t_std;
  array[nquar] int year_index;
  for (t in 1 : nquar) {
    t_std[t] = (t - 0.5 * (nquar + 1)) / nquar;
    year_index[t] = (t - 1) %/% 4;
  }
  if (nyear != nquar %/% 4) {
    reject("nyear must equal nquar / 4");
  }
  
  // Step indicators for ascertainment changes in the abroad-diagnosis field
  matrix[nquar, n_detect_steps] detect_step_design = rep_matrix(0, nquar,
                                                                n_detect_steps);
  for (s in 1 : n_detect_steps) {
    for (t in detect_step_quarter[s] : nquar) {
      detect_step_design[t, s] = 1;
    }
  }
}
parameters {
  // Infection parameters
  vector[npars_inf] beta_inf; // infection parameters
  real log_sigma_inf;
  
  // Total migration (omega) and the transfers-of-care arrivals are
  // parameterised directly
  real beta_omega_raw_1;
  vector[npars_mig - 1] beta_omega_wiggle_raw;
  real log_sigma_omega;
  real alpha_share;
  vector[npars_mig - 1] beta_share_raw;
  real log_sigma_share;
  real gamma_obs_offset;
  real gamma_obs_slope;
  vector[nquar - 1] gamma_obs_drift_raw;
  
  // Diagnosis probability parameters (AR(1))
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
  real alpha_diag_mig; // shared born-abroad diagnosis offset

  // Born-abroad:UK-born ratio of UK-acquired incidence
  vector[npars_cob] beta_cob_raw;
  real<lower=0> sigma_cob; // spline smoothing parameter
  real alpha_cob;
  real<lower=0, upper=0.9> rho_eta;
  
  // Time-constant CD4 composition among undiagnosed arrivals
  simplex[4] undiag_cd4_prop;
  
  // Misclassification of abroad diagnoses: logit trend plus ascertainment change-points
  real alpha_detect;
  real beta_detect_slope;
  vector[n_detect_steps] delta_detect;
  
  // AR(1) coefficient shared across diagnosis-probability walks; (0,1) is mean-reverting
  real<lower=0, upper=1> rho_diag;
}
transformed parameters {
  // Infection process
  real<lower=0> sigma_inf = exp(log_sigma_inf);
  vector[nquar] lambda = X_inf * beta_inf;
  vector<lower=0>[nquar] h = exp(lambda);
  
  // omega (total) and the transfer share are separate splines
  real<lower=0> sigma_omega = exp(log_sigma_omega);
  real<lower=0> sigma_share = exp(log_sigma_share);
  vector[npars_mig] beta_omega;
  beta_omega[1] = beta_omega_raw_1;
  beta_omega[2 : npars_mig] = beta_omega_wiggle_raw * sigma_omega;
  vector[nquar] gamma_omega = X_mig * beta_omega;
  vector<lower=0>[nquar] omega = exp(gamma_omega);
  vector[nquar] logit_transfer_share = alpha_share
                                       + X_mig[ : , 2 : npars_mig]
                                         * (beta_share_raw * sigma_share);
  vector<lower=0, upper=1>[nquar] transfer_share = inv_logit(
                                                    logit_transfer_share);
  vector<lower=0>[nquar] transfers = omega .* transfer_share;
  vector<lower=0>[nquar] undiag_migration = omega .* (1 - transfer_share);
  vector[nquar] gamma = log(omega);
  vector<lower=0>[3] omega_pre = rep_vector(omega[1], 3);
  vector[nquar] gamma_agg;
  {
    real running = sum(omega_pre);
    for (t in 1 : nquar) {
      running += omega[t];
      if (t > 1) {
        running -= t <= 4 ? omega_pre[t - 1] : omega[t - 4];
      }
      gamma_agg[t] = log(running);
    }
  }

  vector[nquar] gamma_obs_offset_t;
  {
    real drift = 0;
    for (t in 1 : nquar) {
      if (t > 1) {
        drift += gamma_obs_drift_sd * gamma_obs_drift_raw[t - 1];
      }
      gamma_obs_offset_t[t] = gamma_obs_offset + gamma_obs_slope * t_std[t]
                              + drift;
    }
  }
  
  // Step variances on natural scale (non-centred)
  vector[4] log_vardelta = mu_log_var + sigma_log_var * log_vardelta_raw;
  vector<lower=0>[4] vardelta = exp(log_vardelta);
  // Diagnosis probabilities via AR(1)
  vector[nquar] delta_1;
  vector[nquar] delta_2;
  vector[nquar] delta_3;
  vector[nquar] delta_4;
  // AR(1), non-centred: t=1 draws from stationary dist, t>1 mean-reverts to alpha_diag_k
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
  matrix<lower=0, upper=1>[4, nquar] d_nat;
  matrix<lower=0, upper=1>[4, nquar] d_mig;
  
  // d_mig adds a shared offset to d_nat, applied uniformly across states
  for (t in 1 : nquar) {
    d_nat[ : , t] = inv_logit([delta_1[t], delta_2[t], delta_3[t], delta_4[t]]');
    d_mig[ : , t] = inv_logit([delta_1[t] + alpha_diag_mig,
                               delta_2[t] + alpha_diag_mig,
                               delta_3[t] + alpha_diag_mig,
                               delta_4[t] + alpha_diag_mig]');
  }
  
  // Born-abroad:UK-born ratio of UK-acquired incidence (non-centred spline)
  vector[npars_cob] beta_cob;
  beta_cob[1] = beta_cob_raw[1];
  beta_cob[2 : npars_cob] = beta_cob_raw[2 : npars_cob] * sigma_cob;
  vector[nquar] eta = alpha_cob + 0.5 * (X_cob * beta_cob);
  vector<lower=0>[nquar] p;
  
  // Derived compatibility outputs for existing post-processing.
  matrix<lower=0>[5, nquar] m;
  matrix[5, nquar] omega_m;
  for (t in 1 : nquar) {
    omega_m[1 : 4, t] = undiag_cd4_prop * undiag_migration[t];
    omega_m[5, t] = transfers[t];
    m[ : , t] = omega_m[ : , t] / omega[t];
  }
  
  // Misclassification model for abroad diagnoses: logit trend in time, plus
  // level shifts at any supplied ascertainment change-points.
  vector[nquar] detect_logit = alpha_detect + beta_detect_slope * t_std
                               + detect_step_design * delta_detect;
  vector<lower=0, upper=1>[nquar] detect_prob;
  vector<lower=0>[nquar] exp_ABROAD_dx;
  matrix<lower=0>[nquar, 2] misclassified_dx = rep_matrix(0, nquar, 2);
  
  for (t in 1 : nquar) {
    detect_prob[t] = inv_logit(detect_logit[t]);
    exp_ABROAD_dx[t] = transfers[t] * detect_prob[t];
    misclassified_dx[t, 2] = transfers[t] * (1 - detect_prob[t]);
  }
  
  // Expected observations
  matrix<lower=0>[nquar, 2] exp_HIV_dx;
  matrix<lower=0>[nquar, 2] exp_AIDS_dx;
  array[nquar, 2] simplex[4] exp_p = rep_array(rep_vector(0.25, 4), nquar, 2);
  array[nyear] vector<lower=0, upper=1>[4] recent_arrival_cd4_prob;
  
  {
  // Local: not saved as draws; only exp_HIV_dx/exp_AIDS_dx/exp_p needed downstream
  array[nquar, 2] vector[9] cum_diags = rep_array(rep_vector
                                                           (0, 9), nquar, 2);
  matrix[nquar, 2] incidence_state_5 = rep_matrix(0, nquar, 2);
  matrix[nquar, 2] incidence_state_6 = rep_matrix(0, nquar, 2);
  matrix[nquar, 2] incidence_state_7 = rep_matrix(0, nquar, 2);
  matrix[nquar, 2] incidence_state_8 = rep_matrix(0, nquar, 2);
  matrix[nquar, 2] incidence_state_9 = rep_matrix(0, nquar, 2);
  array[nquar, 2] vector[4] prev_mat;
  
  // Initial time period (t=1)
  matrix[9, 9] Lambda1_nat = tmat_fct(1, q, d_nat);
  matrix[9, 9] Lambda1_mig = tmat_fct(1, q, d_mig);
  
  cum_diags[1, 1,  : ] = Lambda1_nat' * append_row(init_prev, rep_vector(0, 5));
  cum_diags[1, 1, 1] += h[1];
  
  // Initial born-abroad group
  p[1] = exp(eta[1]);
  cum_diags[1, 2,  : ] = Lambda1_mig' * append_row(init_prev_abroad, rep_vector(0, 5));
  cum_diags[1, 2, 1] += h[1] * p[1];
  cum_diags[1, 2, 1 : 4] += omega_m[1 : 4, 1];
  
  for (t in 2 : nquar) {
    matrix[9, 9] Lambda_t_nat = tmat_fct(t, q, d_nat);
    matrix[9, 9] Lambda_t_mig = tmat_fct(t, q, d_mig);
    
    // Evolution for UK-born group
    cum_diags[t, 1,  : ] = Lambda_t_nat' * cum_diags[t - 1, 1,  : ];
    cum_diags[t, 1, 1] += h[t];
    
    // Evolution for born-abroad group
    p[t] = exp(eta[t]);
    cum_diags[t, 2,  : ] = Lambda_t_mig' * cum_diags[t - 1, 2,  : ];
    cum_diags[t, 2, 1] += h[t] * p[t];
    cum_diags[t, 2, 1 : 4] += omega_m[1 : 4, t];
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
    // UK-born population (y = 1)
    if (exp_HIV_dx[t, 1] > 0) {
      exp_p[t, 1, 1] = incidence_state_6[t, 1] / exp_HIV_dx[t, 1];
      exp_p[t, 1, 2] = incidence_state_7[t, 1] / exp_HIV_dx[t, 1];
      exp_p[t, 1, 3] = incidence_state_8[t, 1] / exp_HIV_dx[t, 1];
      exp_p[t, 1, 4] = incidence_state_9[t, 1] / exp_HIV_dx[t, 1];
    }
    // Born-abroad (y=2): include misclassified arrivals, proportional to
    // m[1:4,t], matching what the simulation puts in CD4 data.
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

  {
  array[nyear] vector[4] exp_recent_arrival_cd4
    = rep_array(rep_vector(0, 4), nyear);
  array[nyear] vector[4] exp_migrant_cd4
    = rep_array(rep_vector(0, 4), nyear);
  for (t in 1 : nquar) {
    exp_migrant_cd4[year_index[t] + 1]
      += exp_HIV_dx[t, 2] * exp_p[t, 2];
  }
  for (t0 in 1 : nquar) {
    int arrival_year = year_index[t0] + 1;
    int last_t = min(nquar, 4 * (arrival_year + 1));
    if (t0 < last_t) {
      vector[4] lat = omega_m[1 : 4, t0];
      for (t in (t0 + 1) : last_t) {
        int diagnosis_year = year_index[t] + 1;
        vector[4] diagnosis = lat .* d_mig[ : , t];
        vector[4] not_diagnosed = lat - diagnosis;
        vector[4] lat_next;
        exp_recent_arrival_cd4[diagnosis_year] += diagnosis;
        lat_next[1] = not_diagnosed[1] * (1 - q[1]);
        lat_next[2] = not_diagnosed[1] * q[1]
                      + not_diagnosed[2] * (1 - q[2]);
        lat_next[3] = not_diagnosed[2] * q[2]
                      + not_diagnosed[3] * (1 - q[3]);
        lat_next[4] = not_diagnosed[3] * q[3]
                      + not_diagnosed[4] * (1 - q[4]);
        lat = lat_next;
      }
    }
  }
  for (year in 1 : nyear) {
    for (state in 1 : 4) {
      recent_arrival_cd4_prob[year, state] = fmin(
        1,
        exp_recent_arrival_cd4[year, state]
        / (exp_migrant_cd4[year, state] + 1e-12)
      );
    }
  }
  }
}
model {
  // PRIORS
  
  // Infection priors
  beta_inf[1] ~ normal(5, 1.5); // quarterly infections in (exp(2),exp(8)) = (8,4390) at +-2 SD
  beta_inf[2 : npars_inf] ~ normal(0, sigma_inf);
  log_sigma_inf ~ normal(-1, 1);
  
  // omega's own level/trend prior; transfer_share's own baseline + smoothness
  beta_omega_raw_1 ~ normal(4, 1);
  beta_omega_wiggle_raw ~ std_normal();
  log_sigma_omega ~ normal(-0.5, 0.5);
  alpha_share ~ normal(logit(0.30), 1);
  beta_share_raw ~ std_normal();
  log_sigma_share ~ normal(-1, 0.5);
  gamma_obs_offset ~ normal(0, 10);
  gamma_obs_slope ~ normal(0, gamma_obs_slope_sd);
  gamma_obs_drift_raw ~ std_normal();

  if (gamma_obs_sd > 0) {
    gamma_obs ~ normal(gamma_agg + gamma_obs_offset_t, gamma_obs_sd);
  }

  // Priors centred on empirical range midpoints:
  // State 1 (CD4 ≥ 500):   2.5–10%  → logit ~ -2.90
  // State 2 (CD4 350–500): 2.5–11%  → logit ~ -2.75
  // State 3 (CD4 200–349): 5–25%    → logit ~ -2.00
  // State 4 (CD4 < 200):   20–60%   → logit ~ -0.50
  alpha_diag_1 ~ normal(-2.9, 0.5);
  alpha_diag_2 ~ normal(-2.75, 0.5);
  alpha_diag_3 ~ normal(-2.0, 0.5);
  alpha_diag_4 ~ normal(-0.5, 0.5);
  delta_raw_1 ~ std_normal();
  delta_raw_2 ~ std_normal();
  delta_raw_3 ~ std_normal();
  delta_raw_4 ~ std_normal();
  // Step SD ~0.13 logit/quarter (stationary ~0.26); avoids heavy-tailed swings
  mu_log_var ~ normal(-4.0, 0.5);
  sigma_log_var ~ normal(0, 0.5);
  log_vardelta_raw ~ std_normal();
  
  // Mode ~0.86: persistent but rules out the random-walk boundary at rho=1
  rho_diag ~ beta(7, 2);
  
  alpha_diag_mig ~ normal(0, 0.5);
  
  // UK-acquired incidence-ratio priors
  alpha_cob ~ normal(0, 1.5);
  beta_cob_raw ~ std_normal();
  sigma_cob ~ gamma(2, 32);
  
  // Approximate melding with independent evidence; eta_sd = 0 disables it.
  if (eta_sd > 0) {
    mean(eta_mean - eta) ~ normal(0, eta_sd);
    eta_mean[2 : nquar] ~ normal(
                                 eta[2 : nquar]
                                 + rho_eta
                                   * (eta_mean[1 : (nquar - 1)]
                                      - eta[1 : (nquar - 1)]),
                                 eta_sd * sqrt(1 - rho_eta ^ 2));
  }
  
  // Detection probability priors
  alpha_detect ~ normal(1.1, 0.5);
  beta_detect_slope ~ std_normal();
  delta_detect ~ std_normal();
  
  // Mean proportions (2/7, 2/7, 2/7, 1/7), with moderate concentration.
  undiag_cd4_prop ~ dirichlet([4, 4, 4, 2]');
  
  // Auxiliary likelihood for detect_prob from a data-linkage study, when observed
  if (detect_obs_sd > 0 && n_detect_obs > 0) {
    for (i in 1 : n_detect_obs) {
      detect_logit_obs[i] ~ normal(detect_logit[detect_obs_index[i]], detect_obs_sd);
    }
  }
  
  // Likelihood
  for (t in 1 : nquar) {
    for (y in 1 : 2) {
      HIV[t, y] ~ poisson(exp_HIV_dx[t, y]);
      AIDS[t, y] ~ poisson(exp_AIDS_dx[t, y]);
      CD4[t, y,  : ] ~ multinomial(exp_p[t, y,  : ]);
    }
  }
  ABROAD ~ poisson(exp_ABROAD_dx);
  for (year in 1 : nyear) {
    for (state in 1 : 4) {
      recent_arrival_cd4[year, state]
        ~ binomial(recent_arrival_cd4_total[year, state],
                   recent_arrival_cd4_prob[year, state]);
    }
  }
}
generated quantities {
  // Pointwise log likelihood
  matrix[nquar, 2] log_lik;
  vector[nquar] log_lik_abroad;
  array[nyear] real log_lik_recent_arrival_cd4;
  for (year in 1 : nyear) {
    log_lik_recent_arrival_cd4[year] = 0;
    for (state in 1 : 4) {
      log_lik_recent_arrival_cd4[year]
        += binomial_lpmf(recent_arrival_cd4[year, state]
                         | recent_arrival_cd4_total[year, state],
                           recent_arrival_cd4_prob[year, state]);
    }
  }
  for (t in 1 : nquar) {
    log_lik_abroad[t] = poisson_lpmf(ABROAD[t] | exp_ABROAD_dx[t]);
    
    for (y in 1 : 2) {
      log_lik[t, y] = poisson_lpmf(HIV[t, y] | exp_HIV_dx[t, y])
                      + poisson_lpmf(AIDS[t, y] | exp_AIDS_dx[t, y])
                      + multinomial_lpmf(CD4[t, y,  : ] | exp_p[t, y,  : ]);
    }
  }
}
