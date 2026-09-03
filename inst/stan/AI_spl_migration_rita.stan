// Age-independent spline model with migration and RITA state
// Separate undiagnosed/transfer flows and a fixed CD4-composition simplex,
// matching AD_ptens_migration_rita.stan.

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
  int<lower=1> nquar; // quarters
  // y = 1: UK-born, y = 2: born abroad
  array[nquar, 2] int<lower=0> HIV;
  array[nquar, 2] int<lower=0> AIDS;
  array[nquar, 2, 5] int<lower=0> CD4;
  array[nquar] int<lower=0> ABROAD; // previous diagnoses abroad
  vector<lower=0, upper=1>[4] q; // progression probabilities
  vector<lower=0>[7] init_prev;
  vector<lower=0>[7] init_prev_abroad;
  
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
}
transformed data {
  vector[nquar] t_std;
  for (t in 1 : nquar) {
    t_std[t] = (t - 0.5 * (nquar + 1)) / nquar;
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
  vector[npars_inf] beta_inf;
  real log_sigma_inf;
  
  // Separate undiagnosed and previously diagnosed migration flows
  vector[npars_mig] beta_undiag_mig_raw;
  vector[npars_mig] beta_transfer_raw;
  real log_sigma_undiag_mig;
  real log_sigma_transfer;
  // Extra smoothness on transfers of care
  real log_sigma_transfer_share;
  
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
  real alpha_diag_mig; // shared born-abroad diagnosis offset
  
  // Born-abroad:UK-born ratio of UK-acquired incidence
  vector[npars_cob] beta_cob_raw;
  real<lower=0> sigma_cob;
  real alpha_cob;
  // Capped below 1: rho_eta = 1 degenerates to a random walk, uninformative about eta's level
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
  // False recency rate
  real<lower=0, upper=1> frr = inv_logit(logit_frr);
  
  // Infection process
  real<lower=0> sigma_inf = exp(log_sigma_inf);
  vector[nquar] lambda = X_inf * beta_inf;
  vector<lower=0>[nquar] h = exp(lambda);
  
  // Separate non-centred splines for undiagnosed arrivals vs transfers;
  // omega and m are retained as derived compatibility outputs.
  real<lower=0> sigma_undiag_mig = exp(log_sigma_undiag_mig);
  real<lower=0> sigma_transfer = exp(log_sigma_transfer);
  real<lower=0> sigma_transfer_share = exp(log_sigma_transfer_share);
  vector[npars_mig] beta_undiag_mig;
  vector[npars_mig] beta_transfer;
  beta_undiag_mig[1] = beta_undiag_mig_raw[1];
  beta_transfer[1] = beta_transfer_raw[1];
  beta_undiag_mig[2 : npars_mig]
    = beta_undiag_mig_raw[2 : npars_mig] * sigma_undiag_mig;
  beta_transfer[2 : npars_mig]
    = beta_transfer_raw[2 : npars_mig] * sigma_transfer;
  vector[nquar] gamma_undiag_mig = X_mig * beta_undiag_mig;
  vector[nquar] gamma_transfer = X_mig * beta_transfer;
  vector<lower=0>[nquar] undiag_migration = exp(gamma_undiag_mig);
  vector<lower=0>[nquar] transfers = exp(gamma_transfer);
  vector<lower=0>[nquar] omega = undiag_migration + transfers;
  vector[nquar] gamma = log(omega);
  
  // Step variances on natural scale (non-centred)
  vector[5] log_vardelta = mu_log_var + sigma_log_var * log_vardelta_raw;
  vector<lower=0>[5] vardelta = exp(log_vardelta);
  // Diagnosis probabilities via AR(1), non-centred: t=1 stationary draw, t>1 mean-reverts to alpha_diag_k
  vector[nquar] delta_1;
  vector[nquar] delta_2;
  vector[nquar] delta_3;
  vector[nquar] delta_4;
  vector[nquar] delta_5;
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
  matrix<lower=0, upper=1>[5, nquar] d_nat;
  matrix<lower=0, upper=1>[5, nquar] d_mig;
  
  // d_mig adds a single shared offset applied uniformly across states.
  for (t in 1 : nquar) {
    d_nat[ : , t] = inv_logit(
                              [delta_1[t], delta_2[t], delta_3[t], delta_4[t],
                               delta_5[t]]');
    d_mig[ : , t] = inv_logit(
                              [delta_1[t] + alpha_diag_mig,
                               delta_2[t] + alpha_diag_mig,
                               delta_3[t] + alpha_diag_mig,
                               delta_4[t] + alpha_diag_mig,
                               delta_5[t] + alpha_diag_mig]');
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
  array[nquar, 2] simplex[5] exp_p = rep_array(rep_vector(0.20, 5), nquar, 2);
  
  {
  // Local: not saved as draws; only exp_HIV_dx/exp_AIDS_dx/exp_p needed downstream
  array[nquar, 2] vector[13] cum_diags = rep_array(
                                                            rep_vector
                                                            (0, 13), nquar,
                                                            2);
  matrix[nquar, 2] incidence_state_8 = rep_matrix(0, nquar, 2);
  matrix[nquar, 2] incidence_state_9 = rep_matrix(0, nquar, 2);
  matrix[nquar, 2] incidence_state_10 = rep_matrix(0, nquar, 2);
  matrix[nquar, 2] incidence_state_11 = rep_matrix(0, nquar, 2);
  matrix[nquar, 2] incidence_state_12 = rep_matrix(0, nquar, 2);
  matrix[nquar, 2] incidence_state_13 = rep_matrix(0, nquar, 2);
  array[nquar, 2] vector[7] prev_mat;
  
  // Initial time period (t=1)
  matrix[13, 13] Lambda1_nat = tmat_fct(1, q, d_nat);
  matrix[13, 13] Lambda1_mig = tmat_fct(1, q, d_mig);
  
  cum_diags[1, 1,  : ] = Lambda1_nat' * append_row(init_prev, rep_vector(0, 6));
  cum_diags[1, 1, 1] += h[1];
  
  // Initial born-abroad group
  p[1] = exp(eta[1]);
  cum_diags[1, 2,  : ] = Lambda1_mig' * append_row(init_prev_abroad, rep_vector(0, 6));
  cum_diags[1, 2, 1] += h[1] * p[1];
  cum_diags[1, 2, 4 : 7] += omega_m[1 : 4, 1];
  
  // Evolution over time
  for (t in 2 : nquar) {
    matrix[13, 13] Lambda_t_nat = tmat_fct(t, q, d_nat);
    matrix[13, 13] Lambda_t_mig = tmat_fct(t, q, d_mig);
    
    // UK-born group
    cum_diags[t, 1,  : ] = Lambda_t_nat' * cum_diags[t - 1, 1,  : ];
    cum_diags[t, 1, 1] += h[t];
    
    // Born-abroad group
    p[t] = exp(eta[t]);
    cum_diags[t, 2,  : ] = Lambda_t_mig' * cum_diags[t - 1, 2,  : ];
    cum_diags[t, 2, 1] += h[t] * p[t];
    cum_diags[t, 2, 4 : 7] += omega_m[1 : 4, t];
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
    // UK-born population (y = 1)
    if (exp_HIV_dx[t, 1] > 0) {
      real p1 = incidence_state_9[t, 1] / exp_HIV_dx[t, 1];
      exp_p[t, 1, 1] = p1 + frr * (1 - p1);
      exp_p[t, 1, 2] = (1 - frr) * incidence_state_10[t, 1] / exp_HIV_dx[t, 1];
      exp_p[t, 1, 3] = (1 - frr) * incidence_state_11[t, 1] / exp_HIV_dx[t, 1];
      exp_p[t, 1, 4] = (1 - frr) * incidence_state_12[t, 1] / exp_HIV_dx[t, 1];
      exp_p[t, 1, 5] = (1 - frr) * incidence_state_13[t, 1] / exp_HIV_dx[t, 1];
    }
    
    // Born-abroad (y=2): include misclassified arrivals in exp_p. They enter
    // non-RITA states (diagnosed 2-5), proportional to m[1:4,t]; FRR applied as usual.
    {
      vector[5] raw_p = [incidence_state_9[t, 2], incidence_state_10[t, 2],
                         incidence_state_11[t, 2], incidence_state_12[t, 2],
                         incidence_state_13[t, 2]]'
                        + 1e-12;
      if (misclassified_dx[t, 2] > 1e-10) {
        real total_lat = sum(m[1 : 4, t]);
        vector[4] m_lat = (total_lat > 1e-10)
                          ? m[1 : 4, t] / total_lat
                          : rep_vector(0.25, 4);
        raw_p[2] += misclassified_dx[t, 2] * m_lat[1];
        raw_p[3] += misclassified_dx[t, 2] * m_lat[2];
        raw_p[4] += misclassified_dx[t, 2] * m_lat[3];
        raw_p[5] += misclassified_dx[t, 2] * m_lat[4];
      }
      real total_p_sum = sum(raw_p);
      if (total_p_sum > 1e-10) {
        vector[5] p1 = raw_p / total_p_sum;
        exp_p[t, 2, 1] = p1[1] + frr * (1 - p1[1]);
        exp_p[t, 2, 2] = (1 - frr) * p1[2];
        exp_p[t, 2, 3] = (1 - frr) * p1[3];
        exp_p[t, 2, 4] = (1 - frr) * p1[4];
        exp_p[t, 2, 5] = (1 - frr) * p1[5];
      } else {
        exp_p[t, 2] = rep_vector(0.2, 5);
      }
    }
  }
  }
}
model {
  // PRIORS
  
  // Infection priors
  beta_inf[1] ~ normal(5, 1.5);
  beta_inf[2 : npars_inf] ~ normal(0, sigma_inf);
  log_sigma_inf ~ normal(-1, 1);
  
  // Separate flow priors
  beta_undiag_mig_raw[1] ~ normal(4 + log(0.70), 1);
  beta_transfer_raw[1] ~ normal(4 + log(0.30), 1);
  beta_undiag_mig_raw[2 : npars_mig] ~ std_normal();
  beta_transfer_raw[2 : npars_mig] ~ std_normal();
  log_sigma_undiag_mig ~ normal(-0.5, 0.5);
  log_sigma_transfer ~ normal(-0.5, 0.5);

  // Penalize quarter-to-quarter change in logit(transfers share of arrivals)
  vector[nquar - 1] share_logit_diff = (gamma_transfer[2 : nquar]
                                        - gamma_undiag_mig[2 : nquar])
                                       - (gamma_transfer[1 : (nquar - 1)]
                                          - gamma_undiag_mig[1 : (nquar - 1)]);
  share_logit_diff ~ normal(0, sigma_transfer_share);
  log_sigma_transfer_share ~ normal(-1.5, 0.5);

  // Priors centred on empirical range midpoints:
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
