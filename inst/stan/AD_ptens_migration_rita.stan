// Age-dependent tensor product model with migration flows and RITA state.
// Same data contract and migration parameterisation as AD_ptens_migration.stan,
// plus the 7-state recent/non-recent structure and false recency rate from
// AD_ptens_rita.stan.

data {
  int<lower=1> nquar; // number of quarters
  int<lower=1> nage; // number of ages
  
  // Observations by quarter, age, and country-of-birth group
  // y = 1: UK-born, y = 2: born abroad
  array[nquar, nage, 2] int<lower=0> HIV_age; // new HIV diagnoses
  array[nquar, nage, 2] int<lower=0> AIDS_age; // new AIDS diagnoses
  array[nquar, nage, 2, 5] int<lower=0> CD4_age; // recent + CD4 categories among HIV diagnoses
  
  // Previous diagnoses abroad, as quarterly totals and four age groups
  array[nquar] int<lower=0> ABROAD;           // quarterly total detected abroad diagnoses
  array[nquar, 4] int<lower=0> ABROAD_group;  // grouped by diagnosis age group

  // Replication weight on ABROAD_group's multinomial likelihood, offsetting its
  // otherwise small influence vs. the much larger migrant HIV_age/CD4_age
  // likelihood on the same age-composition parameter. Default (set in
  // model_type_prep()): sum(HIV_age[,,2]) / sum(ABROAD_group). 1 = no reweighting.
  real<lower=0> abroad_group_weight;

  // Annual numerator and auxiliary-observed non-recent risk set by age/state
  int<lower=1> nyear;
  array[nyear, 4, 4] int<lower=0> recent_arrival_cd4;
  array[nyear, 4, 4] int<lower=0> recent_arrival_cd4_total;
  
  // Progression by age at acquisition or arrival
  matrix<lower=0, upper=1>[4, nage] q_age;
  
  // Initial prevalence at t = 1 by age and latent RITA state (1..7)
  matrix<lower=0>[7, nage] init_prev_age;
  matrix<lower=0>[7, nage] init_prev_abroad_age;
  
  // Infection spline (tensor product, as in AD_ptens)
  int<lower=1> npars_inf;
  matrix[nquar * nage, npars_inf] X_inf;
  matrix[npars_inf - 1, 2 * npars_inf - 2] S1;
  
  // Diagnosis process
  
  // Migration design matrices
  int<lower=1> npars_mig;
  int<lower=1> npars_cob;
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
  
  // Independent rolling-annual migration index (log scale). Complete
  // gamma_obs_group rows use grouped data, else gamma_obs; sd=0 disables this likelihood
  vector[nquar] gamma_obs;
  matrix[nquar, 4] gamma_obs_group;
  real<lower=0> gamma_obs_sd;
  // Quarter-to-quarter drift allowed in log(HIV-positive migration)
  real<lower=0> gamma_obs_drift_sd;
  // Prior sd for the linear trend in log(HIV-positive migration); tighten toward
  // 0 if little drift is expected over the study period
  real<lower=0> gamma_obs_slope_sd;
}
transformed data {
  row_vector[6] empty_vec = rep_row_vector(0, 6);
  vector[7] empty_vec1 = rep_vector(0, 7);
  vector[npars_inf - 1] zero_vec = rep_vector(0, npars_inf - 1);
  
  int np1 = npars_inf - 1;
  // Pre-split penalty matrices (symmetrized to avoid fp noise)
  matrix[np1, np1] S1a = 0.5 * (S1[1 : np1, 1 : np1] + S1[1 : np1, 1 : np1]');
  matrix[np1, np1] S1b = 0.5
                         * (S1[1 : np1, npars_inf : (2 * npars_inf - 2)]
                            + S1[1 : np1, npars_inf : (2 * npars_inf - 2)]');
  
  // Age-group mapping, centred time trend, and fast quarter->year index
  // (birthdays at t = 1, 5, 9, ...)
  array[nage] int age_group_lookup;
  array[nquar] int year_index;
  vector[nquar] t_std;
  for (t in 1 : nquar) {
    year_index[t] = (t - 1) %/% 4;
    t_std[t] = (t - 0.5 * (nquar + 1)) / nquar;
  }
  if (nyear != nquar %/% 4) {
    reject("nyear must equal nquar / 4");
  }
  
  // age 1-10 -> 1, 11-20 -> 2, 21-30 -> 3, 31+ -> 4
  for (a0 in 1 : nage) {
    age_group_lookup[a0] = min(4, 1 + (a0 - 1) %/% 10);
  }
  
  // Pre-computed year difference: current age for cohort (t0, a0) at time t
  // = min(a0 + year_diff[t0, t], nage)
  array[nquar, nquar] int year_diff;
  for (t0 in 1 : nquar) {
    for (t in 1 : nquar) {
      year_diff[t0, t] = year_index[t] - year_index[t0];
    }
  }
  
  // Number of single-year ages within each diagnosis age group
  array[4] int<lower=0> age_group_size;
  for (ag in 1 : 4) {
    age_group_size[ag] = 0;
  }
  for (a in 1 : nage) {
    age_group_size[age_group_lookup[a]] += 1;
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
  // Infection parameters (tensor product)
  vector[npars_inf] beta_inf;
  vector<lower=0>[2] sigma_inf;
  
  // Total migration (omega) and the transfers-of-care arrivals are
  // parameterised directly
  real beta_omega_raw_1;
  vector[npars_mig - 1] beta_omega_wiggle_raw;
  real log_sigma_omega;
  real alpha_share;
  vector[npars_mig - 1] beta_share_raw;
  real log_sigma_share;
  // Offset leaves the level of omega unidentified by the external index
  real gamma_obs_offset;
  real gamma_obs_slope;
  vector[nquar - 1] gamma_obs_drift_raw;

  // Age composition of the born-abroad inflow on the 4 diagnosis age groups,
  // as a softmax with a linear time trend (reference group = 45+).
  vector[3] alpha_abroad_age;
  vector[3] beta_abroad_age;
  // Diagnosis probability parameters (AR(1)): 1 recent + 4 non-recent CD4 states
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
  // Age effect on diagnosis probability: a profile shared across the 5
  // diagnosis states, plus small state-specific deviations
  sum_to_zero_vector[4] gamma_age_shared_raw;
  real<lower=0> sigma_age_shared;
  array[5] sum_to_zero_vector[4] gamma_age_dev_raw;
  real<lower=0> sigma_age_dev;
  real alpha_diag_mig; // shared born-abroad diagnosis offset
  
  // Born-abroad:UK-born ratio of UK-acquired incidence
  vector[npars_cob] beta_cob_raw;
  real<lower=0> sigma_cob;
  real alpha_cob;
  real<lower=0, upper=0.9> rho_eta;
  
  // Time-constant CD4 composition among undiagnosed arrivals
  simplex[4] undiag_cd4_prop;
  
  // Misclassification of abroad diagnoses
  real alpha_detect;
  real beta_detect_slope;
  vector[n_detect_steps] delta_detect;
  
  // AR(1) coefficient shared across diagnosis-probability walks; (0,1) is mean-reverting
  real<lower=0, upper=1> rho_diag;
  
  // False recency rate
  real logit_frr;
}
transformed parameters {
  // Infection intensity by (t, a)
  vector[nquar * nage] lambda = X_inf * beta_inf;
  vector<lower=0>[nquar * nage] h = exp(lambda);
  
  // Tensor product prior precision
  vector<lower=0>[2] tau = inv(square(sigma_inf));
  
  // False recency rate
  real<lower=0, upper=1> frr = inv_logit(logit_frr);
  
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
  
  // Rolling 4-quarter aggregation of omega (log scale), matching the quarterly
  // rolling-annual external migration index. For t=1:3, the window is completed
  // with omega_pre (the three quarters preceding the surveillance window).
  vector[nquar] gamma_agg;
  {
    real running = sum(omega_pre);
    for (t in 1 : nquar) {
      running += omega[t];
      if (t > 1) {
        if (t <= 4) {
          running -= omega_pre[t - 1];
        } else {
          running -= omega[t - 4];
        }
      }
      gamma_agg[t] = log(running);
    }
  }
  
  // Offset between the external all-cause migration index and log(omega), i.e.
  // the log HIV prevalence of the inflow: level, trend, and drift.
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
  
  // Age distribution for allocating born-abroad flows by age; group-level
  // composition is expanded uniformly within each decade band (divide by
  // age_group_size) so per-age values sum to 1.
  matrix<lower=0, upper=1>[nquar, 4] abroad_age_group_prop;
  matrix<lower=0, upper=1>[nquar, nage] abroad_age_prop;
  matrix[nquar, 4] gamma_agg_group;
  for (t in 1 : nquar) {
    vector[4] logits_age;
    logits_age[1 : 3] = alpha_abroad_age + beta_abroad_age * t_std[t];
    logits_age[4] = 0; // reference category
    abroad_age_group_prop[t,  : ] = softmax(logits_age)';
  }
  for (t in 1 : nquar) {
    for (a in 1 : nage) {
      abroad_age_prop[t, a] = abroad_age_group_prop[t, age_group_lookup[a]]
                              / age_group_size[age_group_lookup[a]];
    }
  }
  for (ag in 1 : 4) {
    gamma_agg_group[ : , ag] = gamma_agg + log(abroad_age_group_prop[ : , ag]);
  }
  
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
  delta_1[1] = alpha_diag_1
               + sqrt(vardelta[1]) * inv_stationary * delta_raw_1[1];
  delta_2[1] = alpha_diag_2
               + sqrt(vardelta[2]) * inv_stationary * delta_raw_2[1];
  delta_3[1] = alpha_diag_3
               + sqrt(vardelta[3]) * inv_stationary * delta_raw_3[1];
  delta_4[1] = alpha_diag_4
               + sqrt(vardelta[4]) * inv_stationary * delta_raw_4[1];
  delta_5[1] = alpha_diag_5
               + sqrt(vardelta[5]) * inv_stationary * delta_raw_5[1];
  for (t in 2 : nquar) {
    delta_1[t] = alpha_diag_1 + rho_diag * (delta_1[t - 1] - alpha_diag_1)
                 + sqrt(vardelta[1]) * delta_raw_1[t];
    delta_2[t] = alpha_diag_2 + rho_diag * (delta_2[t - 1] - alpha_diag_2)
                 + sqrt(vardelta[2]) * delta_raw_2[t];
    delta_3[t] = alpha_diag_3 + rho_diag * (delta_3[t - 1] - alpha_diag_3)
                 + sqrt(vardelta[3]) * delta_raw_3[t];
    delta_4[t] = alpha_diag_4 + rho_diag * (delta_4[t - 1] - alpha_diag_4)
                 + sqrt(vardelta[4]) * delta_raw_4[t];
    delta_5[t] = alpha_diag_5 + rho_diag * (delta_5[t - 1] - alpha_diag_5)
                 + sqrt(vardelta[5]) * delta_raw_5[t];
  }
  
  // Age effect on diagnosis probability = shared[ag] + deviation[ag,k]
  // (partial pooling across the 5 diagnosis states). sqrt(4/3) corrects
  // sum-to-zero variance so each component's marginal SD equals its sigma.
  matrix[4, 5] gamma_age;
  real age_effect_scale_shared = sigma_age_shared * sqrt(4.0 / 3.0);
  real age_effect_scale_dev = sigma_age_dev * sqrt(4.0 / 3.0);
  vector[4] gamma_age_shared = age_effect_scale_shared * gamma_age_shared_raw;
  for (k in 1 : 5) {
    for (ag in 1 : 4) {
      gamma_age[ag, k] = gamma_age_shared[ag]
                         + age_effect_scale_dev * gamma_age_dev_raw[k][ag];
    }
  }
  
  // Diagnosis probabilities by age group, state, and time; stored as
  // matrix[nquar, 5] so column k is contiguous (cache-friendly).
  // d_nat: UK-born; d_mig: born abroad (adds a state-specific offset).
  array[4] matrix<lower=0, upper=1>[nquar, 5] d_nat;
  array[4] matrix<lower=0, upper=1>[nquar, 5] d_mig;
  for (ag in 1 : 4) {
    d_nat[ag,  : , 1] = inv_logit(delta_1 + gamma_age[ag, 1]);
    d_nat[ag,  : , 2] = inv_logit(delta_2 + gamma_age[ag, 2]);
    d_nat[ag,  : , 3] = inv_logit(delta_3 + gamma_age[ag, 3]);
    d_nat[ag,  : , 4] = inv_logit(delta_4 + gamma_age[ag, 4]);
    d_nat[ag,  : , 5] = inv_logit(delta_5 + gamma_age[ag, 5]);
    d_mig[ag,  : , 1] = inv_logit(
                                  delta_1 + gamma_age[ag, 1]
                                  + alpha_diag_mig);
    d_mig[ag,  : , 2] = inv_logit(
                                  delta_2 + gamma_age[ag, 2]
                                  + alpha_diag_mig);
    d_mig[ag,  : , 3] = inv_logit(
                                  delta_3 + gamma_age[ag, 3]
                                  + alpha_diag_mig);
    d_mig[ag,  : , 4] = inv_logit(
                                  delta_4 + gamma_age[ag, 4]
                                  + alpha_diag_mig);
    d_mig[ag,  : , 5] = inv_logit(
                                  delta_5 + gamma_age[ag, 5]
                                  + alpha_diag_mig);
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
  matrix<lower=0>[nquar, nage] misclassified_dx = rep_matrix(0, nquar, nage);
  
  for (t in 1 : nquar) {
    detect_prob[t] = inv_logit(detect_logit[t]);
    exp_ABROAD_dx[t] = transfers[t] * detect_prob[t];
    misclassified_dx[t,  : ] = transfers[t] * (1 - detect_prob[t])
                               * abroad_age_prop[t,  : ];
  }
  
  // Expected diagnoses by (t, a, y) and diagnosed state: 1 recent HIV,
  // 2-5 non-recent HIV by CD4 stage, 6 AIDS.
  array[nquar, nage, 2] real<lower=0> exp_HIV_dx;
  array[nquar, nage, 2] real<lower=0> exp_AIDS_dx;
  array[nquar, nage, 2] simplex[5] exp_p;
  array[nyear, 4] vector<lower=0, upper=1>[4] recent_arrival_cd4_prob;
  
  {
  array[nquar, nage, 2] row_vector[6] exp_diag
    = rep_array(empty_vec, nquar, nage, 2);
  
  // For each entry time (t0) and age (a0), propagate diagnosis and
  // progression forward, with ageing in diagnosis probabilities.
  for (t0 in 1 : nquar) {
    // Time-varying born-abroad:UK-born UK-acquired incidence ratio
    p[t0] = exp(eta[t0]);
    
    for (a0 in 1 : nage) {
      // UK-born (y = 1) and born abroad (y = 2)
      for (y in 1 : 2) {
        vector[7] lat;
        vector[7] lat_next;
        vector[7] lat_p = empty_vec1;
        
        // Progression depends on age at entry a0
        real q1 = q_age[1, a0];
        real q2 = q_age[2, a0];
        real q3 = q_age[3, a0];
        real q4 = q_age[4, a0];
        real nq1 = 1 - q1;
        real nq2 = 1 - q2;
        real nq3 = 1 - q3;
        real nq4 = 1 - q4;
        
        if (t0 == 1) {
          int ag0 = age_group_lookup[a0];
          
          // d[ag, t, k] layout: row = time, col = state
          real d1 = (y == 1) ? d_nat[ag0, 1, 1] : d_mig[ag0, 1, 1];
          real d2 = (y == 1) ? d_nat[ag0, 1, 2] : d_mig[ag0, 1, 2];
          real d3 = (y == 1) ? d_nat[ag0, 1, 3] : d_mig[ag0, 1, 3];
          real d4 = (y == 1) ? d_nat[ag0, 1, 4] : d_mig[ag0, 1, 4];
          real d5 = (y == 1) ? d_nat[ag0, 1, 5] : d_mig[ag0, 1, 5];
          real nd1 = 1 - d1;
          real nd2 = 1 - d2;
          real nd3 = 1 - d3;
          real nd4 = 1 - d4;
          real nd5 = 1 - d5;
          
          real v1;
          real v2;
          real v3;
          real v4;
          real v5;
          real v6;
          real v7;
          if (y == 1) {
            v1 = init_prev_age[1, a0];
            v2 = init_prev_age[2, a0];
            v3 = init_prev_age[3, a0];
            v4 = init_prev_age[4, a0];
            v5 = init_prev_age[5, a0];
            v6 = init_prev_age[6, a0];
            v7 = init_prev_age[7, a0];
          } else {
            v1 = init_prev_abroad_age[1, a0];
            v2 = init_prev_abroad_age[2, a0];
            v3 = init_prev_abroad_age[3, a0];
            v4 = init_prev_abroad_age[4, a0];
            v5 = init_prev_abroad_age[5, a0];
            v6 = init_prev_abroad_age[6, a0];
            v7 = init_prev_abroad_age[7, a0];
          }
          
          // diagnoses from prevalent infections at quarter 1 (RITA: states
          // 1-3 share d1)
          exp_diag[t0, a0, y, 1] += (v1 + v2 + v3) * d1;
          exp_diag[t0, a0, y, 2] += v4 * d2;
          exp_diag[t0, a0, y, 3] += v5 * d3;
          exp_diag[t0, a0, y, 4] += v6 * d4;
          exp_diag[t0, a0, y, 5] += v7 * d5;
          exp_diag[t0, a0, y, 6] += v7 * nd5 * q4;
          
          // survivors after diagnosis and progression
          lat_p[1] = 0;
          lat_p[2] = (nd1 * nq1) * v1;
          lat_p[3] = (nd1 * q1) * v1;
          lat_p[4] = (nd1 * nq1) * v2 + (nd2 * nq1) * v4;
          lat_p[5] = (nd1 * q1) * v2 + (nd1 * nq2) * v3 + (nd2 * q1) * v4
                     + (nd3 * nq2) * v5;
          lat_p[6] = (nd1 * q2) * v3 + (nd3 * q2) * v5 + (nd4 * nq3) * v6;
          lat_p[7] = (nd4 * q3) * v6 + (nd5 * nq4) * v7;
        }
        
        // Inflows at entry time t0 and age a0
        real inflow_native = h[(a0 - 1) * nquar + t0];
        
        // For migrants: UK-acquired inflow uses same age/time infection surface, scaled by p[t0]
        real inflow_mig_uk = inflow_native * p[t0];
        
        // Abroad-acquired inflow arrives directly in the non-recent CD4
        // states 4..7 (no recent states 1..3).
        real af = abroad_age_prop[t0, a0];
        vector[4] inflow_mig_ab;
        inflow_mig_ab[1] = omega_m[1, t0] * af;
        inflow_mig_ab[2] = omega_m[2, t0] * af;
        inflow_mig_ab[3] = omega_m[3, t0] * af;
        inflow_mig_ab[4] = omega_m[4, t0] * af;
        
        if (y == 1) {
          lat[1] = inflow_native + lat_p[1];
          lat[2] = lat_p[2];
          lat[3] = lat_p[3];
          lat[4] = lat_p[4];
          lat[5] = lat_p[5];
          lat[6] = lat_p[6];
          lat[7] = lat_p[7];
        } else {
          lat[1] = inflow_mig_uk + lat_p[1];
          lat[2] = lat_p[2];
          lat[3] = lat_p[3];
          lat[4] = inflow_mig_ab[1] + lat_p[4];
          lat[5] = inflow_mig_ab[2] + lat_p[5];
          lat[6] = inflow_mig_ab[3] + lat_p[6];
          lat[7] = inflow_mig_ab[4] + lat_p[7];
        }
        
        if (t0 < nquar) {
          for (t in (t0 + 1) : nquar) {
            int a = min(a0 + year_diff[t0, t], nage);
            int ag = age_group_lookup[a];
            
            // d[ag, t, k]: column k is contiguous in memory, cache-friendly for varying t
            real d1 = (y == 1) ? d_nat[ag, t, 1] : d_mig[ag, t, 1];
            real d2 = (y == 1) ? d_nat[ag, t, 2] : d_mig[ag, t, 2];
            real d3 = (y == 1) ? d_nat[ag, t, 3] : d_mig[ag, t, 3];
            real d4 = (y == 1) ? d_nat[ag, t, 4] : d_mig[ag, t, 4];
            real d5 = (y == 1) ? d_nat[ag, t, 5] : d_mig[ag, t, 5];
            real nd1 = 1 - d1;
            real nd2 = 1 - d2;
            real nd3 = 1 - d3;
            real nd4 = 1 - d4;
            real nd5 = 1 - d5;
            
            // Diagnosis occurs before progression (RITA: states 1-3 share d1)
            exp_diag[t, a, y, 1] += (lat[1] + lat[2] + lat[3]) * d1;
            exp_diag[t, a, y, 2] += lat[4] * d2;
            exp_diag[t, a, y, 3] += lat[5] * d3;
            exp_diag[t, a, y, 4] += lat[6] * d4;
            exp_diag[t, a, y, 5] += lat[7] * d5;
            exp_diag[t, a, y, 6] += lat[7] * nd5 * q4;
            
            // Progression among undiagnosed (states 1-3 always advance)
            lat_next[1] = 0;
            lat_next[2] = (nd1 * nq1) * lat[1];
            lat_next[3] = (nd1 * q1) * lat[1];
            lat_next[4] = (nd1 * nq1) * lat[2] + (nd2 * nq1) * lat[4];
            lat_next[5] = (nd1 * q1) * lat[2] + (nd1 * nq2) * lat[3]
                          + (nd2 * q1) * lat[4] + (nd3 * nq2) * lat[5];
            lat_next[6] = (nd1 * q2) * lat[3] + (nd3 * q2) * lat[5]
                          + (nd4 * nq3) * lat[6];
            lat_next[7] = (nd4 * q3) * lat[6] + (nd5 * nq4) * lat[7];
            
            lat = lat_next;
          }
        }
      }
    }
  }

  // Expected HIV and AIDS diagnoses and CD4/RITA proportions
  for (t in 1 : nquar) {
    for (a in 1 : nage) {
      for (y in 1 : 2) {
        real hiv_dx = sum(exp_diag[t, a, y, 1 : 5]);
        real aids_dx = exp_diag[t, a, y, 6];
        
        // Add misclassified abroad diagnoses to migrant HIV totals only
        if (y == 2) {
          hiv_dx += misclassified_dx[t, a];
        }
        
        exp_HIV_dx[t, a, y] = hiv_dx;
        exp_AIDS_dx[t, a, y] = aids_dx;
        
        // CD4/RITA composition: include misclassified arrivals for migrants,
        // added to the non-recent categories (2-5) before the FRR adjustment;
        // otherwise each migrant cell's likelihood independently pushes
        // detect_prob -> 99% to eliminate the mismatch.
        {
          vector[5] raw_p = to_vector(exp_diag[t, a, y, 1 : 5]) + 1e-12;
          if (y == 2 && misclassified_dx[t, a] > 1e-10) {
            real total_lat = sum(m[1 : 4, t]);
            vector[4] m_lat = (total_lat > 1e-10) ? m[1 : 4, t] / total_lat
                              : rep_vector(0.25, 4);
            raw_p[2 : 5] += misclassified_dx[t, a] * m_lat;
          }
          real total_p = sum(raw_p);
          vector[5] p1 = (total_p > 1e-10) ? raw_p / total_p
                         : rep_vector(0.2, 5);
          exp_p[t, a, y, 1] = p1[1] + frr * (1 - p1[1]);
          exp_p[t, a, y, 2] = (1 - frr) * p1[2];
          exp_p[t, a, y, 3] = (1 - frr) * p1[3];
          exp_p[t, a, y, 4] = (1 - frr) * p1[4];
          exp_p[t, a, y, 5] = (1 - frr) * p1[5];
        }
      }
    }
  }
  }
  
  // Diagnoses from infections present on arrival, limited to the diagnosis
  // calendar year and the preceding arrival year. Arrivals enter the
  // non-recent CD4 states directly, so this uses d_mig[.,.,2:5] rather than
  // the recent-diagnosis category.
  {
  array[nyear, 4] vector[4] exp_recent_arrival_cd4
    = rep_array(rep_vector(0, 4), nyear, 4);
  array[nyear, 4] vector[4] exp_migrant_cd4
    = rep_array(rep_vector(0, 4), nyear, 4);
  for (t in 1 : nquar) {
    int diagnosis_year = year_index[t] + 1;
    for (a in 1 : nage) {
      int ag = age_group_lookup[a];
      exp_migrant_cd4[diagnosis_year, ag]
        += exp_HIV_dx[t, a, 2] * exp_p[t, a, 2, 2 : 5];
    }
  }
  for (t0 in 1 : nquar) {
    int arrival_year = year_index[t0] + 1;
    int last_t = min(nquar, 4 * (arrival_year + 1));
    if (t0 < last_t) {
      for (a0 in 1 : nage) {
        vector[4] lat = omega_m[1 : 4, t0] * abroad_age_prop[t0, a0];
        real q1 = q_age[1, a0];
        real q2 = q_age[2, a0];
        real q3 = q_age[3, a0];
        real q4 = q_age[4, a0];
        for (t in (t0 + 1) : last_t) {
          int diagnosis_year = year_index[t] + 1;
          int a = min(a0 + year_diff[t0, t], nage);
          int ag = age_group_lookup[a];
          vector[4] diagnosis;
          vector[4] not_diagnosed;
          vector[4] lat_next;
          diagnosis[1] = lat[1] * d_mig[ag, t, 2];
          diagnosis[2] = lat[2] * d_mig[ag, t, 3];
          diagnosis[3] = lat[3] * d_mig[ag, t, 4];
          diagnosis[4] = lat[4] * d_mig[ag, t, 5];
          not_diagnosed = lat - diagnosis;
          exp_recent_arrival_cd4[diagnosis_year, ag] += diagnosis;
          lat_next[1] = not_diagnosed[1] * (1 - q1);
          lat_next[2] = not_diagnosed[1] * q1
                        + not_diagnosed[2] * (1 - q2);
          lat_next[3] = not_diagnosed[2] * q2
                        + not_diagnosed[3] * (1 - q3);
          lat_next[4] = not_diagnosed[3] * q3
                        + not_diagnosed[4] * (1 - q4);
          lat = lat_next;
        }
      }
    }
  }
  for (year in 1 : nyear) {
    for (ag in 1 : 4) {
      for (state in 1 : 4) {
        recent_arrival_cd4_prob[year, ag, state] = fmin(
          1,
          (1 - frr) * exp_recent_arrival_cd4[year, ag, state]
          / (exp_migrant_cd4[year, ag, state] + 1e-12)
        );
      }
    }
  }
  }
}
model {
  // Infection surface priors
  // beta_inf[1]: log-scale intercept; normal(5.0, 1.5) covers [7, 3000]
  // infections/quarter at ±2 SD, matching observed incidence scale.
  sigma_inf ~ std_normal();
  beta_inf[1] ~ normal(5.0, 1.5);
  
  matrix[np1, np1] K1_prec = S1a * tau[1] + S1b * tau[2];
  beta_inf[2 : npars_inf] ~ multi_normal_prec(zero_vec, K1_prec);
  
  // omega's own level/trend prior; transfer_share's own baseline + smoothness
  beta_omega_raw_1 ~ normal(4, 1);
  beta_omega_wiggle_raw ~ std_normal();
  log_sigma_omega ~ normal(-0.5, 0.5);
  alpha_share ~ normal(logit(0.30), 1);
  beta_share_raw ~ std_normal();
  log_sigma_share ~ normal(-1, 0.5);

  // Age composition of previous diagnoses abroad, with a linear time trend
  alpha_abroad_age ~ normal(0, 1.5);
  beta_abroad_age ~ std_normal();
  for (t in 1 : nquar) {
    target += abroad_group_weight
              * multinomial_lpmf(ABROAD_group[t]
                                 | to_vector(abroad_age_group_prop[t,  : ]));
  }
  
  // Diagnosis priors centred on empirical range midpoints:
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
  
  gamma_age_shared_raw ~ std_normal();
  sigma_age_shared ~ normal(0, 0.5);
  for (k in 1 : 5) {
    gamma_age_dev_raw[k] ~ std_normal();
  }
  // Tighter than sigma_age_shared: state-specific deviations expected to be small
  sigma_age_dev ~ normal(0, 0.2);

  alpha_diag_mig ~ normal(0, 0.5);
  
  // UK-acquired incidence-ratio priors
  alpha_cob ~ normal(0, 1.5);
  beta_cob_raw ~ std_normal();
  sigma_cob ~ gamma(2, 32);
  
  // Approximate melding with independent evidence; eta_sd = 0 disables it
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
  
  // Independent migration index informs shape, not level or long-run growth
  gamma_obs_offset ~ normal(0, 10);
  gamma_obs_slope ~ normal(0, gamma_obs_slope_sd);
  gamma_obs_drift_raw ~ std_normal();
  if (gamma_obs_sd > 0) {
    for (t in 1 : nquar) {
      int has_group = 1;
      for (ag in 1 : 4) {
        if (gamma_obs_group[t, ag] < 0) {
          has_group = 0;
        }
      }
      if (has_group == 1) {
        gamma_obs_group[t, ]
          ~ normal(gamma_agg_group[t, ]
                   + rep_row_vector(gamma_obs_offset_t[t], 4),
                   gamma_obs_sd);
      } else {
        gamma_obs[t] ~ normal(gamma_agg[t] + gamma_obs_offset_t[t],
                              gamma_obs_sd);
      }
    }
  }
  
  // Likelihood
  // As in AD_ptens, omit age == 1 from the likelihood (often structurally zero)
  for (t in 1 : nquar) {
    for (y in 1 : 2) {
      HIV_age[t, 2 : nage, y] ~ poisson(exp_HIV_dx[t, 2 : nage, y]);
      AIDS_age[t, 2 : nage, y] ~ poisson(exp_AIDS_dx[t, 2 : nage, y]);
      for (a in 2 : nage) {
        CD4_age[t, a, y,  : ] ~ multinomial(exp_p[t, a, y]);
      }
    }
  }
  
  ABROAD ~ poisson(exp_ABROAD_dx);
  for (year in 1 : nyear) {
    for (ag in 1 : 4) {
      for (state in 1 : 4) {
        recent_arrival_cd4[year, ag, state]
          ~ binomial(recent_arrival_cd4_total[year, ag, state],
                     recent_arrival_cd4_prob[year, ag, state]);
      }
    }
  }
}
generated quantities {
  // Log-likelihood contributions for model comparison
  array[nquar, nage, 2] real log_lik;
  vector[nquar] log_lik_abroad;
  array[nyear, 4] real log_lik_recent_arrival_cd4;

  for (year in 1 : nyear) {
    for (ag in 1 : 4) {
      log_lik_recent_arrival_cd4[year, ag] = 0;
      for (state in 1 : 4) {
        log_lik_recent_arrival_cd4[year, ag]
          += binomial_lpmf(recent_arrival_cd4[year, ag, state]
                           | recent_arrival_cd4_total[year, ag, state],
                             recent_arrival_cd4_prob[year, ag, state]);
      }
    }
  }
  
  for (t in 1 : nquar) {
    log_lik_abroad[t] = poisson_lpmf(ABROAD[t] | exp_ABROAD_dx[t])
                        + abroad_group_weight
                          * multinomial_lpmf(ABROAD_group[t]
                                             | to_vector(
                                                abroad_age_group_prop[t,  : ]));
    
    for (y in 1 : 2) {
      // a == 1 is omitted from the sampling likelihood, so its log-lik
      // contribution is set to 0 (not computed)
      log_lik[t, 1, y] = 0;
      for (a in 2 : nage) {
        log_lik[t, a, y] = poisson_lpmf(HIV_age[t, a, y] |
                                        exp_HIV_dx[t, a, y])
                           + poisson_lpmf(AIDS_age[t, a, y] |
                                          exp_AIDS_dx[t, a, y])
                           + multinomial_lpmf(CD4_age[t, a, y,  : ] |
                                              exp_p[t, a, y]);
      }
    }
  }
}
