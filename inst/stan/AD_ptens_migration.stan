// Age-dependent tensor product model including migration

data {
  int<lower=1> nquar; // number of quarters
  int<lower=1> nage; // number of ages
  
  // Observations by quarter (t), age (a), and group (y)
  // y = 1: native, y = 2: migrant
  array[nquar, nage, 2] int<lower=0> HIV_age; // new HIV diagnoses
  array[nquar, nage, 2] int<lower=0> AIDS_age; // new AIDS diagnoses
  array[nquar, nage, 2, 4] int<lower=0> CD4_age; // CD4 categories among HIV diagnoses
  
  // Age distribution of prior diagnoses from abroad by quarter.
  // The model derives quarterly totals as ABROAD[t] = sum_a ABROAD_age[t,a].
  array[nquar, nage] int<lower=0> ABROAD_age;
  
  // Progression probabilities by age at entry (age at infection for natives/UK-acquired,
  // age at arrival for abroad-acquired migrants)
  matrix<lower=0, upper=1>[4, nage] q_age;
  
  // Initial prevalence at t = 1 by age and latent CD4 state (1..4)
  matrix<lower=0>[4, nage] init_prev_age;
  matrix<lower=0>[4, nage] init_prev_abroad_age;
  
  // Infection spline (tensor product, as in AD_ptens)
  int<lower=1> npars_inf;
  matrix[nquar * nage, npars_inf] X_inf;
  matrix[npars_inf - 1, 2 * npars_inf - 2] S1;
  
  // Time-varying diagnosis random walk (plus age-group offsets)
  // (npars_diag and X_diag removed; diagnosis now uses random walk)
  
  // Migration design matrices
  int<lower=1> npars_mig;
  int<lower=1> npars_cob;
  int<lower=1> npars_mp;
  matrix[nquar, npars_mig] X_mig;
  matrix[nquar, npars_cob] X_cob;
  // Basis for the undiagnosed-mix softmax and the detection-probability spline.
  matrix[nquar, npars_mp] X_mp;
  
  // Prior information for evidence synthesis (log proportion UK-acquired)
  vector[nquar] eta_mean;
  real<lower=0> eta_sd;
  
  // External data-linkage estimates of detection probability (on logit scale).
  // Observations can be sparse, e.g. occasional linkage studies.
  int<lower=0> n_detect_obs;
  array[n_detect_obs] int<lower=1, upper=nquar> detect_obs_index;
  vector[n_detect_obs] detect_logit_obs;
  real<lower=0> detect_obs_sd;
  
  // External auxiliary observations informing the TREND of log migration volume
  // (gamma = log omega). Intended for an INDEPENDENT external migration index on
  // an arbitrary (typically much larger) scale, e.g. ONS total long-term
  // immigration to England, reported as QUARTERLY ROLLING ANNUAL (year-ending)
  // counts. Each observation is matched to gamma_agg[idx] = log of the sum of
  // omega over the four quarters ending at idx, so the rolling-annual series is
  // used directly (no unstable deconvolution to quarterly values). The series
  // enters as
  //   gamma_obs_mean[i] ~ normal(gamma_agg[idx] + gamma_obs_offset, gamma_obs_sd).
  // The free parameter gamma_obs_offset absorbs the unknown level difference
  // (log ratio of total migration to HIV-positive migration), so only the
  // SHAPE/trend of the series constrains gamma; its absolute level (and so the
  // level of omega) stays identified by the HIV/AIDS/CD4/ABROAD counts. Provide
  // log-scale values (e.g. log of ONS counts), indexed by the END quarter of
  // each annual window (index >= 4). Do NOT feed in the posterior of another
  // model fit to these same counts: that double-uses the data. Set
  // gamma_obs_sd = 0 to disable.
  int<lower=0> n_gamma_obs;
  array[n_gamma_obs] int<lower=4, upper=nquar> gamma_obs_index;
  vector[n_gamma_obs] gamma_obs_mean;
  real<lower=0> gamma_obs_sd;
}
transformed data {
  row_vector[5] empty_vec = rep_row_vector(0, 5);
  vector[4] empty_vec1 = rep_vector(0, 4);
  vector[npars_inf - 1] zero_vec = rep_vector(0, npars_inf - 1);
  
  int nagroup = 4; // number of age groups
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
  
  // Derived totals for likelihood and log_lik
  array[nquar] int<lower=0> ABROAD;
  for (t in 1 : nquar) {
    ABROAD[t] = sum(ABROAD_age[t]);
  }
}
parameters {
  // Infection parameters (tensor product)
  vector[npars_inf] beta_inf;
  vector<lower=0>[2] sigma_inf;
  
  // Migration process
  vector[npars_mig] beta_mig;
  real log_sigma_mig;
  // Level offset between the external migration index (e.g. ONS total migration
  // counts) and log HIV-positive migration volume gamma = log(omega). The
  // external series informs only the SHAPE/trend of gamma; this free offset
  // absorbs the (large, unknown) scale difference so the level of omega stays
  // identified by the HIV/AIDS/CD4/ABROAD likelihood. When the gamma_obs hook is
  // disabled it is simply sampled from its prior (harmless).
  real gamma_obs_offset;
  // Diagnosis probability parameters (random walk)
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
  array[4] sum_to_zero_vector[4] gamma_age_raw;
  real<lower=0> sigma_age; // shared shrinkage for age effects
  vector[4] alpha_diag_mig; // migrant-specific diagnosis intercept offset per state
  
  // Proportion UK-acquired among migrants
  vector[npars_cob] beta_cob;
  real<lower=0> sigma_cob;
  real alpha_cob;
  real<lower=-1, upper=1> rho_eta;
  
  // Hybrid transfer-share parameterization: m[5,t] follows a linear logit
  // trend (intercept + slope), reflecting that the diagnosed fraction of
  // arriving HIV-positive migrants changes on decade scales, not quarterly.
  // The undiagnosed CD4 mix m[1:4,t] retains a time-varying softmax on X_mp.
  vector[4] alpha_mp;
  real beta_mp_5;
  vector[npars_mp] beta_mp_1_raw;
  vector[npars_mp] beta_mp_2_raw;
  vector[npars_mp] beta_mp_3_raw;
  real<lower=0> sigma_mp_1;
  real<lower=0> sigma_mp_2;
  real<lower=0> sigma_mp_3;
  
  // Misclassification of abroad diagnoses.
  // detect_prob follows a linear logit trend (intercept + slope) rather than a
  // spline: ABROAD only constrains omega * m5 * detect_prob, so a flexible
  // detect_prob on the shared X_mp basis competes with omega and m[1:4] and
  // leaves the product non-identifiable. A low-dimensional trend converts
  // ABROAD into a tight constraint on omega * m5 (and hence the misclassified
  // inflow), letting the migrant diagnosis flow identify omega * (1 - m5).
  real alpha_detect;
  real beta_detect_slope;
  
  // AR(1) coefficient shared across diagnosis-probability walks.
  // Constrained to (0, 1) so the process is mean-reverting but persistent.
  real<lower=0, upper=1> rho_diag;
}
transformed parameters {
  // Infection intensity by (t, a)
  vector[nquar * nage] lambda = X_inf * beta_inf;
  vector<lower=0>[nquar * nage] h = exp(lambda);
  
  // Tensor product prior precision
  vector<lower=0>[2] tau = inv(square(sigma_inf));
  
  // Migration process
  real<lower=0> sigma_mig = exp(log_sigma_mig);
  vector[nquar] gamma = X_mig * beta_mig;
  vector<lower=0>[nquar] omega = exp(gamma);
  
  // Rolling 4-quarter (annual) aggregation of omega on the log scale, matching
  // the quarterly rolling-annual external migration index. gamma_agg[t] is the
  // log of the migration volume over the four quarters ending at t; only
  // t >= 4 carries a full annual window and is supplied as an observation.
  vector[nquar] gamma_agg;
  {
    real running = 0;
    for (t in 1 : nquar) {
      running += omega[t];
      if (t > 4) {
        running -= omega[t - 4];
      }
      gamma_agg[t] = log(running);
    }
  }
  
  // Age distribution used for allocating migrant flows by age
  matrix<lower=0, upper=1>[nquar, nage] abroad_age_prop;
  for (t in 1 : nquar) {
    real total = ABROAD[t]; // ABROAD is int array; promoted to real here
    if (total > 1e-6) {
      for (a in 1 : nage) {
        abroad_age_prop[t, a] = ABROAD_age[t, a] / total;
      }
    } else {
      abroad_age_prop[t,  : ] = rep_row_vector(1.0 / nage, nage);
    }
  }
  
  // Step variances on natural scale (non-centred)
  vector[4] log_vardelta = mu_log_var + sigma_log_var * log_vardelta_raw;
  vector<lower=0>[4] vardelta = exp(log_vardelta);
  // Diagnosis probabilities via random walk
  vector[nquar] delta_1;
  vector[nquar] delta_2;
  vector[nquar] delta_3;
  vector[nquar] delta_4;
  // AR(1) non-centred parameterisation.
  // t = 1: initial condition drawn from stationary distribution
  //   (variance = step_variance / (1 - rho^2)).
  // t > 1: mean reversion toward alpha_diag_k prevents endpoint drift.
  real inv_stationary = inv_sqrt(1 - square(rho_diag));
  delta_1[1] = alpha_diag_1
               + sqrt(vardelta[1]) * inv_stationary * delta_raw_1[1];
  delta_2[1] = alpha_diag_2
               + sqrt(vardelta[2]) * inv_stationary * delta_raw_2[1];
  delta_3[1] = alpha_diag_3
               + sqrt(vardelta[3]) * inv_stationary * delta_raw_3[1];
  delta_4[1] = alpha_diag_4
               + sqrt(vardelta[4]) * inv_stationary * delta_raw_4[1];
  for (t in 2 : nquar) {
    delta_1[t] = alpha_diag_1 + rho_diag * (delta_1[t - 1] - alpha_diag_1)
                 + sqrt(vardelta[1]) * delta_raw_1[t];
    delta_2[t] = alpha_diag_2 + rho_diag * (delta_2[t - 1] - alpha_diag_2)
                 + sqrt(vardelta[2]) * delta_raw_2[t];
    delta_3[t] = alpha_diag_3 + rho_diag * (delta_3[t - 1] - alpha_diag_3)
                 + sqrt(vardelta[3]) * delta_raw_3[t];
    delta_4[t] = alpha_diag_4 + rho_diag * (delta_4[t - 1] - alpha_diag_4)
                 + sqrt(vardelta[4]) * delta_raw_4[t];
  }
  
  // State-specific zero-sum age effects: each CD4 state has its own
  // age-group profile, reflecting real-data evidence that the age disparity
  // in diagnosis probability differs across CD4 strata.
  // The sqrt(4 / 3) factor corrects the constrained-vector variance so each
  // age effect has marginal SD sigma_age rather than sigma_age * sqrt(3 / 4).
  matrix[4, 4] gamma_age;
  real age_effect_scale = sigma_age * sqrt(4.0 / 3.0);
  for (k in 1 : 4) {
    for (ag in 1 : 4) {
      gamma_age[ag, k] = age_effect_scale * gamma_age_raw[k][ag];
    }
  }
  
  // Diagnosis probabilities by age group (ag), state (k), and time (t).
  // Stored as matrix[nquar, 4] so the inner time-loop accesses column k,
  // which is contiguous in column-major memory (cache-friendly).
  // d_nat: natives; d_mig: migrants (adds alpha_diag_mig offset per state).
  array[4] matrix<lower=0, upper=1>[nquar, 4] d_nat;
  array[4] matrix<lower=0, upper=1>[nquar, 4] d_mig;
  for (ag in 1 : 4) {
    d_nat[ag,  : , 1] = inv_logit(delta_1 + gamma_age[ag, 1]);
    d_nat[ag,  : , 2] = inv_logit(delta_2 + gamma_age[ag, 2]);
    d_nat[ag,  : , 3] = inv_logit(delta_3 + gamma_age[ag, 3]);
    d_nat[ag,  : , 4] = inv_logit(delta_4 + gamma_age[ag, 4]);
    d_mig[ag,  : , 1] = inv_logit(
                                  delta_1 + gamma_age[ag, 1]
                                  + alpha_diag_mig[1]);
    d_mig[ag,  : , 2] = inv_logit(
                                  delta_2 + gamma_age[ag, 2]
                                  + alpha_diag_mig[2]);
    d_mig[ag,  : , 3] = inv_logit(
                                  delta_3 + gamma_age[ag, 3]
                                  + alpha_diag_mig[3]);
    d_mig[ag,  : , 4] = inv_logit(
                                  delta_4 + gamma_age[ag, 4]
                                  + alpha_diag_mig[4]);
  }
  
  // Proportion UK-acquired among migrants
  vector[nquar] eta = alpha_cob + 0.5 * (X_cob * beta_cob);
  vector<lower=0>[nquar] p;
  
  // Non-centred reconstruction of the softmax components for the undiagnosed
  // migrant CD4 mix (states 1:4, reference = state 4).
  vector[npars_mp] beta_mp_1;
  vector[npars_mp] beta_mp_2;
  vector[npars_mp] beta_mp_3;
  beta_mp_1[1] = beta_mp_1_raw[1];
  beta_mp_2[1] = beta_mp_2_raw[1];
  beta_mp_3[1] = beta_mp_3_raw[1];
  beta_mp_1[2 : npars_mp] = beta_mp_1_raw[2 : npars_mp] * sigma_mp_1;
  beta_mp_2[2 : npars_mp] = beta_mp_2_raw[2 : npars_mp] * sigma_mp_2;
  beta_mp_3[2 : npars_mp] = beta_mp_3_raw[2 : npars_mp] * sigma_mp_3;
  
  vector[nquar] theta_1 = alpha_mp[2] + X_mp * beta_mp_1;
  vector[nquar] theta_2 = alpha_mp[3] + X_mp * beta_mp_2;
  vector[nquar] theta_3 = alpha_mp[4] + X_mp * beta_mp_3;
  
  // Transfer share: linear trend on logit scale. Source-country testing
  // coverage evolves on decade scales so a slow drift is all that is
  // identifiable from the ABROAD signal.
  vector[nquar] theta_5 = alpha_mp[1] + beta_mp_5 * t_std;
  vector<lower=0, upper=1>[nquar] m5 = inv_logit(theta_5);
  matrix<lower=0>[5, nquar] m;
  
  for (t in 1 : nquar) {
    vector[4] logits_u;
    logits_u[1] = theta_1[t];
    logits_u[2] = theta_2[t];
    logits_u[3] = theta_3[t];
    logits_u[4] = 0; // reference category
    m[1 : 4, t] = (1 - m5[t]) * softmax(logits_u);
    m[5, t] = m5[t];
  }
  
  // Precompute omega_m[k, t] = omega[t] * m[k, t] (5 compartments x nquar).
  // Avoids recomputing omega[t0]*m[k,t0] for every a0.
  matrix[5, nquar] omega_m;
  for (t in 1 : nquar) {
    omega_m[ : , t] = m[ : , t] * omega[t];
  }
  
  // Misclassification model for abroad diagnoses: linear logit trend in time.
  vector[nquar] detect_logit = alpha_detect + beta_detect_slope * t_std;
  vector<lower=0, upper=1>[nquar] detect_prob;
  vector<lower=0>[nquar] true_abroad_dx;
  vector<lower=0>[nquar] exp_ABROAD_dx;
  matrix<lower=0>[nquar, nage] misclassified_dx = rep_matrix(0, nquar, nage);
  
  for (t in 1 : nquar) {
    detect_prob[t] = inv_logit(detect_logit[t]);
    true_abroad_dx[t] = omega_m[5, t]; // precomputed omega[t] * m[5, t]
    exp_ABROAD_dx[t] = true_abroad_dx[t] * detect_prob[t];
    // Vectorised row assignment replaces inner a loop
    misclassified_dx[t,  : ] = true_abroad_dx[t] * (1 - detect_prob[t])
                               * abroad_age_prop[t,  : ];
  }
  
  // Expected diagnoses by (t, a, y) and latent state (1..4) plus AIDS (5)
  array[nquar, nage, 2] row_vector<lower=0>[5] exp_diag;
  
  // Expected observations
  array[nquar, nage, 2] real<lower=0> exp_HIV_dx;
  array[nquar, nage, 2] real<lower=0> exp_AIDS_dx;
  array[nquar, nage, 2] simplex[4] exp_p;
  
  // Initialize expected diagnoses
  exp_diag = rep_array(empty_vec, nquar, nage, 2);
  for (t in 1 : nquar) {
    for (a in 1 : nage) {
      for (y in 1 : 2) {
        exp_HIV_dx[t, a, y] = 0;
        exp_AIDS_dx[t, a, y] = 0;
        exp_p[t, a, y] = rep_vector(0.25, 4);
      }
    }
  }
  
  // Main forward calculation: for each entry time (t0) and entry age (a0),
  // propagate diagnosis and progression forward with ageing in diagnosis probabilities.
  for (t0 in 1 : nquar) {
    // time-varying UK-acquired proportion
    p[t0] = exp(eta[t0]);
    
    for (a0 in 1 : nage) {
      // Native (y=1) and Migrant (y=2)
      for (y in 1 : 2) {
        vector[4] lat;
        vector[4] lat_next;
        vector[4] lat_p = empty_vec1;
        
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
          real nd1 = 1 - d1;
          real nd2 = 1 - d2;
          real nd3 = 1 - d3;
          real nd4 = 1 - d4;
          
          real v1;
          real v2;
          real v3;
          real v4;
          if (y == 1) {
            v1 = init_prev_age[1, a0];
            v2 = init_prev_age[2, a0];
            v3 = init_prev_age[3, a0];
            v4 = init_prev_age[4, a0];
          } else {
            v1 = init_prev_abroad_age[1, a0];
            v2 = init_prev_abroad_age[2, a0];
            v3 = init_prev_abroad_age[3, a0];
            v4 = init_prev_abroad_age[4, a0];
          }
          
          // diagnoses from prevalent infections at quarter 1
          exp_diag[t0, a0, y, 1] += v1 * d1;
          exp_diag[t0, a0, y, 2] += v2 * d2;
          exp_diag[t0, a0, y, 3] += v3 * d3;
          exp_diag[t0, a0, y, 4] += v4 * d4;
          exp_diag[t0, a0, y, 5] += v4 * nd4 * q4;
          
          // survivors after diagnosis and progression
          lat_p[1] = (nd1 * nq1) * v1;
          lat_p[2] = (nd1 * q1) * v1 + (nd2 * nq2) * v2;
          lat_p[3] = (nd2 * q2) * v2 + (nd3 * nq3) * v3;
          lat_p[4] = (nd3 * q3) * v3 + (nd4 * nq4) * v4;
        }
        
        // Inflows at entry time t0 and age a0
        real inflow_native = h[(a0 - 1) * nquar + t0];
        
        // For migrants: UK-acquired inflow uses same age/time infection surface, scaled by p[t0]
        real inflow_mig_uk = inflow_native * p[t0];
        
        // For migrants: abroad-acquired inflow arrives already in latent states 1..4.
        // omega_m[k,t0] = omega[t0]*m[k,t0] precomputed; multiply by age proportion
        // to get omega[t0]*abroad_age_prop[t0,a0]*m[k,t0] at the cost of 4 mults.
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
        } else {
          lat[1] = inflow_mig_uk + inflow_mig_ab[1] + lat_p[1];
          lat[2] = inflow_mig_ab[2] + lat_p[2];
          lat[3] = inflow_mig_ab[3] + lat_p[3];
          lat[4] = inflow_mig_ab[4] + lat_p[4];
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
            real nd1 = 1 - d1;
            real nd2 = 1 - d2;
            real nd3 = 1 - d3;
            real nd4 = 1 - d4;
            
            // Diagnosis occurs before progression within the quarter
            exp_diag[t, a, y, 1] += lat[1] * d1;
            exp_diag[t, a, y, 2] += lat[2] * d2;
            exp_diag[t, a, y, 3] += lat[3] * d3;
            exp_diag[t, a, y, 4] += lat[4] * d4;
            exp_diag[t, a, y, 5] += lat[4] * nd4 * q4;
            
            // Progression among undiagnosed
            lat_next[1] = (nd1 * nq1) * lat[1];
            lat_next[2] = (nd1 * q1) * lat[1] + (nd2 * nq2) * lat[2];
            lat_next[3] = (nd2 * q2) * lat[2] + (nd3 * nq3) * lat[3];
            lat_next[4] = (nd3 * q3) * lat[3] + (nd4 * nq4) * lat[4];
            
            lat = lat_next;
          }
        }
      }
    }
  }
  
  // Expected HIV and AIDS diagnoses and CD4 proportions
  for (t in 1 : nquar) {
    for (a in 1 : nage) {
      for (y in 1 : 2) {
        real hiv_dx = sum(exp_diag[t, a, y, 1 : 4]);
        real aids_dx = exp_diag[t, a, y, 5];
        
        // Add misclassified abroad diagnoses to migrant HIV totals only
        if (y == 2) {
          hiv_dx += misclassified_dx[t, a];
        }
        
        exp_HIV_dx[t, a, y] = hiv_dx;
        exp_AIDS_dx[t, a, y] = aids_dx;
        
        // CD4 composition: include misclassified arrivals for migrants,
        // proportional to m[1:4] (their latent-state distribution at arrival).
        // Without this, the multinomial likelihood for each of nquar*nage migrant
        // cells independently pushes detect_prob->99% to eliminate the mismatch.
        {
          vector[4] raw_p = to_vector(exp_diag[t, a, y, 1 : 4]) + 1e-12;
          if (y == 2 && misclassified_dx[t, a] > 1e-10) {
            real total_lat = sum(m[1 : 4, t]);
            vector[4] m_lat = (total_lat > 1e-10) ? m[1 : 4, t] / total_lat
                              : rep_vector(0.25, 4);
            raw_p += misclassified_dx[t, a] * m_lat;
          }
          real total_p = sum(raw_p);
          exp_p[t, a, y] = (total_p > 1e-10) ? raw_p / total_p
                           : rep_vector(0.25, 4);
        }
      }
    }
  }
}
model {
  // Priors for infection surface
  // beta_inf[1] is the log-scale intercept; normal(5.0, 1.5) covers [7, 3000]
  // infections/quarter at ±2 SD, matching observed incidence scale.
  sigma_inf ~ normal(0, 5);
  beta_inf[1] ~ normal(5.0, 1.5);
  
  matrix[np1, np1] K1_prec = S1a * tau[1] + S1b * tau[2];
  beta_inf[2 : npars_inf] ~ multi_normal_prec(zero_vec, K1_prec);
  
  // Migration priors
  // log_sigma_mig ~ normal(-0.5, 0.5): sigma_mig mode ~0.6, range [0.22, 1.65].
  // The previous normal(-2, 0.5) gave sigma_mig ~0.14, which is too tight for
  // increasing/step_change omega patterns where log(omega) shifts by ~1 unit
  // over the study; spline coefficients needed 7+ SD from prior — essentially
  // impossible, causing omega to be estimated with spurious precision.
  beta_mig[1] ~ normal(4, 1);
  beta_mig[2 : npars_mig] ~ normal(0, sigma_mig);
  log_sigma_mig ~ normal(-0.5, 0.5);
  
  // Diagnosis priors — centred on midpoints of observed empirical ranges.
  // State 1 (CD4 ≥ 500):    2.5–10%  → logit midpoint ~ -2.90
  // State 2 (CD4 350–500):  2.5–11%  → logit midpoint ~ -2.75
  // State 3 (CD4 200–349):  5–25%    → logit midpoint ~ -2.00
  // State 4 (CD4 < 200):    20–60%   → logit midpoint ~ -0.50
  alpha_diag_1 ~ normal(-2.9, 0.5);
  alpha_diag_2 ~ normal(-2.75, 0.5);
  alpha_diag_3 ~ normal(-2.0, 0.5);
  alpha_diag_4 ~ normal(-0.5, 0.5);
  delta_raw_1 ~ std_normal();
  delta_raw_2 ~ std_normal();
  delta_raw_3 ~ std_normal();
  delta_raw_4 ~ std_normal();
  // mu_log_var: centres exp(-4.0) ~ 0.018 => step SD ~0.13 logit/quarter,
  // stationary SD ~0.26 logit. SD=0.5 avoids the heavy upper tail of normal(-3.5,1)
  // which allowed unrealistically large swings in diagnosis probability.
  mu_log_var ~ normal(-4.0, 0.5);
  sigma_log_var ~ normal(0, 0.5);
  log_vardelta_raw ~ std_normal();
  
  // AR(1) coefficient: prior with mode ~ 0.86 (high persistence) but not a
  // pure random walk. beta(7, 2) has mode = 6/7 ~ 0.857 and mean = 7/9 ~ 0.78,
  // placing almost no mass below 0.5 and ruling out the RW1 boundary at 1.
  rho_diag ~ beta(7, 2);
  
  for (k in 1 : 4) {
    gamma_age_raw[k] ~ std_normal();
  }
  sigma_age ~ normal(0, 0.5);
  // States 1–3: normal(0, 1) allows up to ~2 logit-unit migrant offset.
  // State 4 (CD4 < 200): tighter normal(0, 0.25) — late-stage diagnosis is
  // largely symptom-driven, so native/migrant rates should be similar.
  // This prevents alpha_diag_mig[4] absorbing age-effect confounding.
  alpha_diag_mig[1 : 3] ~ normal(0, 1);
  alpha_diag_mig[4] ~ normal(0, 0.25);
  
  // Proportion UK-acquired priors
  alpha_cob ~ normal(0, 1.5);
  beta_cob[1] ~ std_normal();
  beta_cob[2 : npars_cob] ~ normal(0, sigma_cob);
  sigma_cob ~ gamma(2, 32);
  
  // Approximate melding for proportion UK-acquired - AR(1).
  // Guard: only active when eta_sd > 0 (i.e. independent external evidence is
  // available). Set eta_sd = 0 to let the likelihood alone identify eta — this
  // is the correct choice when the RITA-positive COB proportions used to derive
  // eta_mean are already in the diagnosis likelihood.
  if (eta_sd > 0) {
    eta_mean[1] ~ normal(eta[1], eta_sd);
    eta_mean[2 : nquar] ~ normal(
                                 eta[2 : nquar]
                                 + rho_eta
                                   * (eta_mean[1 : (nquar - 1)]
                                      - eta[1 : (nquar - 1)]),
                                 eta_sd * sqrt(1 - rho_eta ^ 2));
  }
  rho_eta ~ normal(0, 0.5);
  
  // Detection probability priors
  // detect_prob is not separately identified from ABROAD alone: the ABROAD
  // likelihood constrains omega * m[5] * detect_prob jointly. Restricting
  // detect_prob to a linear logit trend (intercept + slope) makes ABROAD a
  // tight constraint on omega * m5, which in turn pins the misclassified inflow
  // omega * m5 * (1 - detect_prob) and frees the migrant diagnosis flow to
  // identify omega * (1 - m5). The slope prior allows a ~13pp drift over the
  // full study period at 2 SD without forcing a trend; an external
  // data-linkage detection rate (detect_logit_obs) can still anchor the level.
  alpha_detect ~ normal(1.1, 0.5);
  beta_detect_slope ~ normal(0, 1);
  
  // Transfer share priors: vague intercept so m5 can sit well away from 25%
  // when data support it; unit-scale slope allows a ~13pp drift over the full
  // study period at 2 SD without forcing a trend. The undiagnosed CD4 mix
  // keeps a time-varying softmax.
  alpha_mp[1] ~ normal(logit(0.25), 2);
  alpha_mp[2 : 4] ~ std_normal();
  beta_mp_5 ~ normal(0, 1);
  beta_mp_1_raw ~ std_normal();
  beta_mp_2_raw ~ std_normal();
  beta_mp_3_raw ~ std_normal();
  sigma_mp_1 ~ gamma(2, 8);
  sigma_mp_2 ~ gamma(2, 8);
  sigma_mp_3 ~ gamma(2, 8);
  
  // Auxiliary likelihood for detect_prob from data-linkage study.
  // Only included when detect_obs_sd > 0 and a linkage study is observed.
  if (detect_obs_sd > 0 && n_detect_obs > 0) {
    for (i in 1 : n_detect_obs) {
      detect_logit_obs[i] ~ normal(detect_logit[detect_obs_index[i]], detect_obs_sd);
    }
  }
  
  // Auxiliary likelihood for the TREND of log migration volume from an EXTERNAL
  // migration index (e.g. ONS total immigration). The series sits on an
  // arbitrary, typically much larger scale; gamma_obs_offset absorbs that level
  // difference so only its trend constrains gamma = log(omega). Each observation
  // is matched to gamma_agg, the log rolling 4-quarter (annual) sum of omega, so
  // quarterly rolling-annual counts are used directly. The level of omega
  // remains identified by the HIV/AIDS/CD4/ABROAD counts below. Legitimate only
  // when the index is independent of those counts. Set gamma_obs_sd = 0 to
  // disable.
  gamma_obs_offset ~ normal(0, 10);
  if (gamma_obs_sd > 0 && n_gamma_obs > 0) {
    for (i in 1 : n_gamma_obs) {
      gamma_obs_mean[i] ~ normal(gamma_agg[gamma_obs_index[i]] + gamma_obs_offset,
                                 gamma_obs_sd);
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
}
generated quantities {
  // Log-likelihood contributions for model comparison
  array[nquar, nage, 2] real log_lik;
  vector[nquar] log_lik_abroad;
  
  for (t in 1 : nquar) {
    log_lik_abroad[t] = poisson_lpmf(ABROAD[t] | exp_ABROAD_dx[t]);
    
    for (y in 1 : 2) {
      for (a in 1 : nage) {
        // keep a==1 defined, even though it is omitted from sampling likelihood
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
