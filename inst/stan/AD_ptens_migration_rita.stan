// Age-dependent tensor product model including migration and RITA state

data {
  int<lower=1> nquar; // number of quarters
  int<lower=1> nage; // number of ages
  
  // Observations by quarter (t), age (a), and group (y)
  // y = 1: native, y = 2: migrant
  array[nquar, nage, 2] int<lower=0> HIV_age;
  array[nquar, nage, 2] int<lower=0> AIDS_age;
  array[nquar, nage, 2, 5] int<lower=0> CD4_age;
  
  // Age distribution of prior diagnoses from abroad by quarter.
  array[nquar, nage] int<lower=0> ABROAD_age;
  
  // Progression probabilities by age at entry
  matrix<lower=0, upper=1>[4, nage] q_age;
  
  // Initial prevalence at t = 1 by age and latent RITA state (1..7)
  matrix<lower=0>[7, nage] init_prev_age;
  matrix<lower=0>[7, nage] init_prev_abroad_age;
  
  // Infection spline (tensor product)
  int<lower=1> npars_inf;
  matrix[nquar * nage, npars_inf] X_inf;
  matrix[npars_inf - 1, 2 * npars_inf - 2] S1;
  
  // Time-varying diagnosis random walk (plus age-group offsets)
  // (npars_diag and X_diag removed; diagnosis now uses random walk)
  
  // Migration splines
  int<lower=1> npars_mig;
  int<lower=1> npars_cob;
  int<lower=1> npars_mp;
  matrix[nquar, npars_mig] X_mig;
  matrix[nquar, npars_cob] X_cob;
  matrix[nquar, npars_mp] X_mp;
  
  // Prior information for evidence synthesis (log proportion UK-acquired)
  vector[nquar] eta_mean;
  real<lower=0> eta_sd;
}
transformed data {
  row_vector[6] empty_vec = rep_row_vector(0, 6);
  vector[7] empty_vec1 = rep_vector(0, 7);
  vector[npars_inf - 1] zero_vec = rep_vector(0, npars_inf - 1);
  
  int nagroup = 4;
  int np1 = npars_inf - 1;
  
  // Pre-split penalty matrices (symmetrized to avoid fp noise)
  matrix[np1, np1] S1a = 0.5 * (S1[1 : np1, 1 : np1] + S1[1 : np1, 1 : np1]');
  matrix[np1, np1] S1b = 0.5
                         * (S1[1 : np1, npars_inf : (2 * npars_inf - 2)]
                            + S1[1 : np1, npars_inf : (2 * npars_inf - 2)]');
  
  // Age-group mapping and fast quarter->year index
  array[nage] int age_group_lookup;
  array[nquar] int year_index;
  for (t in 1 : nquar) {
    year_index[t] = (t - 1) %/% 4;
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
  
  // Derived totals
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
  
  // Diagnosis probability parameters (random walk)
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
  vector[4 * 5] alpha_age_diag; // 4 age groups x 5 states
  vector[5] alpha_diag_mig; // migrant-specific diagnosis intercept offset per state
  
  // Proportion UK-acquired among migrants
  vector[npars_cob] beta_cob;
  real<lower=0> sigma_cob;
  real alpha_cob;
  real<lower=-1, upper=1> rho_eta;
  
  // Probability of migration into each compartment (1..5)
  vector[npars_mp] beta_mp_1_raw;
  vector[npars_mp] beta_mp_2_raw;
  vector[npars_mp] beta_mp_3_raw;
  vector[npars_mp] beta_mp_5_raw;
  real<lower=0> sigma_mp_1;
  real<lower=0> sigma_mp_2;
  real<lower=0> sigma_mp_3;
  real<lower=0> sigma_mp_5;
  vector[4] alpha_mp;
  
  // Misclassification of abroad diagnoses
  real alpha_detect;
  vector[npars_mp] beta_detect;
  real<lower=0> sigma_detect;
  
  // AR(1) coefficient shared across diagnosis-probability walks.
  // Constrained to (0, 1) so the process is mean-reverting but persistent.
  real<lower=0, upper=1> rho_diag;
  
  // False recency rate
  real logit_frr;
}
transformed parameters {
  // Infection intensity by (t, a)
  vector[nquar * nage] lambda = X_inf * beta_inf;
  vector<lower=0>[nquar * nage] h = exp(lambda);
  
  // False recency rate
  real<lower=0, upper=1> frr = inv_logit(logit_frr);
  
  // Tensor product prior precision
  vector<lower=0>[2] tau = inv(square(sigma_inf));
  
  // Migration process
  real<lower=0> sigma_mig = exp(log_sigma_mig);
  vector[nquar] gamma = X_mig * beta_mig;
  vector<lower=0>[nquar] omega = exp(gamma);
  
  // Age distribution used for allocating migrant flows by age
  matrix<lower=0, upper=1>[nquar, nage] abroad_age_prop;
  for (t in 1 : nquar) {
    real total = ABROAD[t];
    if (total > 1e-6) {
      for (a in 1 : nage) {
        abroad_age_prop[t, a] = ABROAD_age[t, a] / total;
      }
    } else {
      abroad_age_prop[t,  : ] = rep_row_vector(1.0 / nage, nage);
    }
  }
  
  // Step variances on natural scale (non-centred)
  vector[5] log_vardelta = mu_log_var + sigma_log_var * log_vardelta_raw;
  vector<lower=0>[5] vardelta = exp(log_vardelta);
  // Diagnosis probabilities via random walk
  vector[nquar] delta_1;
  vector[nquar] delta_2;
  vector[nquar] delta_3;
  vector[nquar] delta_4;
  vector[nquar] delta_5;
  // AR(1) non-centred parameterisation.
  // t = 1: initial condition drawn from stationary distribution
  //   (variance = step_variance / (1 - rho^2)).
  // t > 1: mean reversion toward alpha_diag_k prevents endpoint drift.
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
  
  // Diagnosis probabilities by age group (ag), state (k), and time (t).
  // Stored as matrix[nquar, 5] so the inner time-loop accesses column k,
  // which is contiguous in column-major memory (cache-friendly).
  // d_nat: natives; d_mig: migrants (adds alpha_diag_mig offset per state).
  array[4] matrix<lower=0, upper=1>[nquar, 5] d_nat;
  array[4] matrix<lower=0, upper=1>[nquar, 5] d_mig;
  for (ag in 1 : 4) {
    int ag_offset = (ag - 1) * 5;
    d_nat[ag,  : , 1] = inv_logit(delta_1 + alpha_age_diag[ag_offset + 1]);
    d_nat[ag,  : , 2] = inv_logit(delta_2 + alpha_age_diag[ag_offset + 2]);
    d_nat[ag,  : , 3] = inv_logit(delta_3 + alpha_age_diag[ag_offset + 3]);
    d_nat[ag,  : , 4] = inv_logit(delta_4 + alpha_age_diag[ag_offset + 4]);
    d_nat[ag,  : , 5] = inv_logit(delta_5 + alpha_age_diag[ag_offset + 5]);
    d_mig[ag,  : , 1] = inv_logit(delta_1 + alpha_age_diag[ag_offset + 1] + alpha_diag_mig[1]);
    d_mig[ag,  : , 2] = inv_logit(delta_2 + alpha_age_diag[ag_offset + 2] + alpha_diag_mig[2]);
    d_mig[ag,  : , 3] = inv_logit(delta_3 + alpha_age_diag[ag_offset + 3] + alpha_diag_mig[3]);
    d_mig[ag,  : , 4] = inv_logit(delta_4 + alpha_age_diag[ag_offset + 4] + alpha_diag_mig[4]);
    d_mig[ag,  : , 5] = inv_logit(delta_5 + alpha_age_diag[ag_offset + 5] + alpha_diag_mig[5]);
  }
  
  // Proportion UK-acquired among migrants
  vector[nquar] eta = alpha_cob + 0.5 * (X_cob * beta_cob);
  vector<lower=0>[nquar] p;
  
  // Non-centred reconstruction of migration compartment splines
  vector[npars_mp] beta_mp_1;
  vector[npars_mp] beta_mp_2;
  vector[npars_mp] beta_mp_3;
  vector[npars_mp] beta_mp_5;
  beta_mp_1[1] = beta_mp_1_raw[1];
  beta_mp_2[1] = beta_mp_2_raw[1];
  beta_mp_3[1] = beta_mp_3_raw[1];
  beta_mp_5[1] = beta_mp_5_raw[1];
  beta_mp_1[2 : npars_mp] = beta_mp_1_raw[2 : npars_mp] * sigma_mp_1;
  beta_mp_2[2 : npars_mp] = beta_mp_2_raw[2 : npars_mp] * sigma_mp_2;
  beta_mp_3[2 : npars_mp] = beta_mp_3_raw[2 : npars_mp] * sigma_mp_3;
  beta_mp_5[2 : npars_mp] = beta_mp_5_raw[2 : npars_mp] * sigma_mp_5;
  
  // Migration probabilities into compartments (softmax, reference = compartment 4)
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
  
  // Precompute omega_m[k, t] = omega[t] * m[k, t] (5 compartments x nquar).
  // Avoids recomputing omega[t0]*m[k,t0] for every a0.
  matrix[5, nquar] omega_m;
  for (t in 1 : nquar) {
    omega_m[ : , t] = m[ : , t] * omega[t];
  }
  
  // Misclassification model for abroad diagnoses
  vector[nquar] detect_logit = alpha_detect + 0.5 * (X_mp * beta_detect);
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
  
  // Expected diagnoses by (t, a, y) and diagnosed state (1..5 HIV, 6 AIDS)
  array[nquar, nage, 2] row_vector<lower=0>[6] exp_diag;
  array[nquar, nage, 2] real<lower=0> exp_HIV_dx;
  array[nquar, nage, 2] real<lower=0> exp_AIDS_dx;
  array[nquar, nage, 2] simplex[5] exp_p;
  
  exp_diag = rep_array(empty_vec, nquar, nage, 2);
  for (t in 1 : nquar) {
    for (a in 1 : nage) {
      for (y in 1 : 2) {
        exp_HIV_dx[t, a, y] = 0;
        exp_AIDS_dx[t, a, y] = 0;
        exp_p[t, a, y] = rep_vector(0.2, 5);
      }
    }
  }
  
  // Forward calculation
  for (t0 in 1 : nquar) {
    p[t0] = exp(eta[t0]);
    
    for (a0 in 1 : nage) {
      for (y in 1 : 2) {
        vector[7] lat;
        vector[7] lat_next;
        vector[7] lat_p = empty_vec1;
        
        // Progression rates (depend only on age at infection a0, not current time)
        real q1 = q_age[1, a0];
        real q2 = q_age[2, a0];
        real q3 = q_age[3, a0];
        real q4 = q_age[4, a0];
        real nq1 = 1 - q1;
        real nq2 = 1 - q2;
        real nq3 = 1 - q3;
        real nq4 = 1 - q4;
        
        // Prevalent infections at t0 == 1
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
          
          // diagnoses from prevalent infections at quarter 1
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
        real inflow_mig_uk = inflow_native * p[t0];
        
        // Migrants acquired abroad enter states 4..7 only (no recent states 1..3).
        // omega_m[k,t0] = omega[t0]*m[k,t0] precomputed; multiply by age proportion.
        real af = abroad_age_prop[t0, a0];
        vector[7] inflow_mig_ab = rep_vector(0, 7);
        inflow_mig_ab[4] = omega_m[1, t0] * af;
        inflow_mig_ab[5] = omega_m[2, t0] * af;
        inflow_mig_ab[6] = omega_m[3, t0] * af;
        inflow_mig_ab[7] = omega_m[4, t0] * af;
        
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
          lat[4] = inflow_mig_ab[4] + lat_p[4];
          lat[5] = inflow_mig_ab[5] + lat_p[5];
          lat[6] = inflow_mig_ab[6] + lat_p[6];
          lat[7] = inflow_mig_ab[7] + lat_p[7];
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
  
  // Expected HIV and AIDS diagnoses and CD4 proportions
  for (t in 1 : nquar) {
    for (a in 1 : nage) {
      for (y in 1 : 2) {
        real hiv_dx = sum(exp_diag[t, a, y, 1 : 5]);
        real aids_dx = exp_diag[t, a, y, 6];
        
        // Add misclassified abroad diagnoses to migrant HIV totals only (excluded from CD4 composition)
        if (y == 2) {
          hiv_dx += misclassified_dx[t, a];
        }
        
        exp_HIV_dx[t, a, y] = hiv_dx;
        exp_AIDS_dx[t, a, y] = aids_dx;
        
        // CD4 composition: include misclassified arrivals for migrants.
        // Misclassified migrants arrive in non-RITA latent states 4-7,
        // corresponding to diagnosed states 2-5. Their distribution at arrival
        // is proportional to m[1:4,t]. FRR is applied afterwards as usual.
        {
          vector[5] raw_p = to_vector(exp_diag[t, a, y, 1 : 5]) + 1e-12;
          if (y == 2 && misclassified_dx[t, a] > 1e-10) {
            real total_lat = sum(m[1 : 4, t]);
            vector[4] m_lat = (total_lat > 1e-10)
                              ? m[1 : 4, t] / total_lat
                              : rep_vector(0.25, 4);
            raw_p[2] += misclassified_dx[t, a] * m_lat[1];
            raw_p[3] += misclassified_dx[t, a] * m_lat[2];
            raw_p[4] += misclassified_dx[t, a] * m_lat[3];
            raw_p[5] += misclassified_dx[t, a] * m_lat[4];
          }
          real total_p_sum = sum(raw_p);
          if (total_p_sum > 1e-10) {
            vector[5] p1 = raw_p / total_p_sum;
            exp_p[t, a, y, 1] = p1[1] + frr * (1 - p1[1]);
            exp_p[t, a, y, 2] = (1 - frr) * p1[2];
            exp_p[t, a, y, 3] = (1 - frr) * p1[3];
            exp_p[t, a, y, 4] = (1 - frr) * p1[4];
            exp_p[t, a, y, 5] = (1 - frr) * p1[5];
          } else {
            exp_p[t, a, y] = rep_vector(0.2, 5);
          }
        }
      }
    }
  }
}
model {
  // Infection surface priors.
  // half-normal(0, 5) is a middle ground: loose enough that genuine tail
  // spikes supported by the RITA signal are not suppressed, but tight enough
  // beta_inf[1] is the log-scale intercept; normal(5.0, 1.5) covers [7, 3000]
  // infections/quarter at ±2 SD, matching observed incidence scale.
  sigma_inf ~ normal(0, 5);
  beta_inf[1] ~ normal(5.0, 1.5);
  
  matrix[np1, np1] K1_prec = S1a * tau[1] + S1b * tau[2];
  beta_inf[2 : npars_inf] ~ multi_normal_prec(zero_vec, K1_prec);
  
  // Migration priors
  beta_mig[1] ~ normal(4, 1);
  beta_mig[2 : npars_mig] ~ normal(0, sigma_mig);
  log_sigma_mig ~ normal(-2, 0.5);
  
  // Diagnosis priors — centred on midpoints of observed empirical ranges.
  // State 1 (recently acquired / RITA+): 5–30%   → logit midpoint ~ -1.75
  // State 2 (CD4 ≥ 500):                2.5–10%  → logit midpoint ~ -2.90
  // State 3 (CD4 350–500):              2.5–11%  → logit midpoint ~ -2.75
  // State 4 (CD4 200–349):              5–25%    → logit midpoint ~ -2.00
  // State 5 (CD4 < 200):                20–60%   → logit midpoint ~ -0.50
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
  
  // False recency rate prior centred on logit(0.019)
  logit_frr ~ normal(-3.94, 0.32);
  
  alpha_age_diag ~ normal(0, 1);
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
  alpha_detect ~ normal(1.1, 0.5);
  beta_detect[1] ~ normal(0, 0.5);
  beta_detect[2 : npars_mp] ~ normal(0, sigma_detect);
  sigma_detect ~ gamma(2, 32);
  
  // Migration distribution priors
  alpha_mp[1] ~ normal(0, 1.5);
  alpha_mp[2 : 4] ~ std_normal();
  
  beta_mp_1_raw ~ std_normal();
  beta_mp_2_raw ~ std_normal();
  beta_mp_3_raw ~ std_normal();
  beta_mp_5_raw ~ std_normal();
  
  sigma_mp_1 ~ gamma(2, 32);
  sigma_mp_2 ~ gamma(2, 32);
  sigma_mp_3 ~ gamma(2, 32);
  sigma_mp_5 ~ gamma(2, 32);
  
  // Likelihood (omit age == 1)
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
  array[nquar, nage, 2] real log_lik;
  vector[nquar] log_lik_abroad;
  for (t in 1 : nquar) {
    log_lik_abroad[t] = poisson_lpmf(ABROAD[t] | exp_ABROAD_dx[t]);
    for (y in 1 : 2) {
      for (a in 1 : nage) {
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
