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
  
  // Time-varying diagnosis splines (as in AI_spl_migration), plus age-group offsets
  int<lower=1> npars_diag;
  matrix[nquar, npars_diag] X_diag;
  
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
  row_vector[5] empty_vec = rep_row_vector(0, 5);
  vector[4] empty_vec1 = rep_vector(0, 4);
  vector[npars_inf - 1] zero_vec = rep_vector(0, npars_inf - 1);
  
  int nagroup = 4; // number of age groups
  int np1 = npars_inf - 1;
  
  // Pre-split penalty matrices (symmetrized to avoid fp noise)
  matrix[np1, np1] S1a = 0.5 * (S1[1 : np1, 1 : np1] + S1[1 : np1, 1 : np1]');
  matrix[np1, np1] S1b = 0.5 * (S1[1 : np1, npars_inf : (2 * npars_inf - 2)] + S1[1 : np1, npars_inf : (2 * npars_inf - 2)]');
  
  // Age-group mapping and fast quarter->year index (birthdays at t = 1, 5, 9, ...)
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
  
  // Diagnosis probability parameters (time spline + age-group offset)
  vector[npars_diag] beta_diag_1;
  vector[npars_diag] beta_diag_2;
  vector[npars_diag] beta_diag_3;
  vector[npars_diag] beta_diag_4;
  real<lower=0> sigma_diag_1;
  real<lower=0> sigma_diag_2;
  real<lower=0> sigma_diag_3;
  real<lower=0> sigma_diag_4;
  real alpha_diag_1;
  real alpha_diag_2;
  real alpha_diag_3;
  real alpha_diag_4;
  vector[4 * 4] alpha_age_diag; // 4 age groups x 4 states
  
  // Proportion UK-acquired among migrants
  vector[npars_cob] beta_cob;
  real<lower=0> sigma_cob;
  real alpha_cob;
  real<lower=-1, upper=1> rho_eta;
  
  // Probability of migration into each compartment (1..5)
  vector[npars_mp] beta_mp_1;
  vector[npars_mp] beta_mp_2;
  vector[npars_mp] beta_mp_3;
  vector[npars_mp] beta_mp_5;
  real<lower=0> sigma_mp_1;
  real<lower=0> sigma_mp_2;
  real<lower=0> sigma_mp_3;
  real<lower=0> sigma_mp_5;
  vector[4] alpha_mp;
  
  // Misclassification of abroad diagnoses
  real alpha_detect;
  vector[npars_mp] beta_detect;
  real<lower=0> sigma_detect;
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
  
  // Diagnosis probabilities by time
  vector[nquar] delta_1 = alpha_diag_1 + 0.5 * (X_diag * beta_diag_1);
  vector[nquar] delta_2 = alpha_diag_2 + 0.5 * (X_diag * beta_diag_2);
  vector[nquar] delta_3 = alpha_diag_3 + 0.5 * (X_diag * beta_diag_3);
  vector[nquar] delta_4 = alpha_diag_4 + 0.7 * (X_diag * beta_diag_4);
  
  // Diagnosis probabilities by age group (ag), state (k), and time (t).
  // Stored as matrix[nquar, 4] so the inner time-loop accesses column k,
  // which is contiguous in column-major memory (cache-friendly).
  array[4] matrix<lower=0, upper=1>[nquar, 4] d;
  for (ag in 1 : 4) {
    int ag_offset = (ag - 1) * 4;
    d[ag,  : , 1] = inv_logit(delta_1 + alpha_age_diag[ag_offset + 1]);
    d[ag,  : , 2] = inv_logit(delta_2 + alpha_age_diag[ag_offset + 2]);
    d[ag,  : , 3] = inv_logit(delta_3 + alpha_age_diag[ag_offset + 3]);
    d[ag,  : , 4] = inv_logit(delta_4 + alpha_age_diag[ag_offset + 4]);
  }
  
  // Proportion UK-acquired among migrants
  vector[nquar] eta = alpha_cob + 0.5 * (X_cob * beta_cob);
  vector<lower=0>[nquar] p;
  
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
          real d1 = d[ag0, 1, 1];
          real d2 = d[ag0, 1, 2];
          real d3 = d[ag0, 1, 3];
          real d4 = d[ag0, 1, 4];
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
            real d1 = d[ag, t, 1];
            real d2 = d[ag, t, 2];
            real d3 = d[ag, t, 3];
            real d4 = d[ag, t, 4];
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
        
        // Add misclassified abroad diagnoses to migrant HIV totals only (excluded from CD4 composition)
        if (y == 2) {
          hiv_dx += misclassified_dx[t, a];
        }
        
        exp_HIV_dx[t, a, y] = hiv_dx;
        exp_AIDS_dx[t, a, y] = aids_dx;
        
        // CD4 composition excludes misclassified diagnoses for migrants
        {
          real denom = hiv_dx;
          if (y == 2) {
            denom -= misclassified_dx[t, a];
          }
          
          if (denom > 0) {
            exp_p[t, a, y] = to_vector(exp_diag[t, a, y, 1 : 4]) / denom;
          } else {
            exp_p[t, a, y] = rep_vector(0.25, 4);
          }
        }
      }
    }
  }
}
model {
  // Priors for infection surface
  sigma_inf ~ student_t(4, 0, 200);
  
  beta_inf[1] ~ normal(0, 30);
  matrix[np1, np1] K1_prec = S1a * tau[1] + S1b * tau[2];
  beta_inf[2 : npars_inf] ~ multi_normal_prec(zero_vec, K1_prec);
  
  // Migration priors
  beta_mig[1] ~ normal(4, 1);
  beta_mig[2 : npars_mig] ~ normal(0, sigma_mig);
  log_sigma_mig ~ normal(-2, 1);
  
  // Diagnosis priors (time spline + age group offsets)
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
  
  sigma_diag_1 ~ gamma(2, 16);
  sigma_diag_2 ~ gamma(2, 16);
  sigma_diag_3 ~ gamma(2, 16);
  sigma_diag_4 ~ gamma(2, 16);
  
  alpha_age_diag ~ normal(0, 1);
  
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
