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
  matrix[nquar * nage, npars_inf] X;
  matrix[npars_inf - 1, 2 * npars_inf - 2] S1;
  
  // Time-varying diagnosis splines (plus age-group offsets)
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
  row_vector[6] empty_vec = rep_row_vector(0, 6);
  vector[7] empty_vec1 = rep_vector(0, 7);
  vector[npars_inf - 1] zero_vec = rep_vector(0, npars_inf - 1);
  
  int nagroup = 4;
  
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
  
  // Derived totals
  array[nquar] int<lower=0> ABROAD;
  for (t in 1 : nquar) {
    ABROAD[t] = sum(ABROAD_age[t]);
  }
}
parameters {
  // Infection parameters (tensor product)
  vector[npars_inf] beta;
  vector<lower=0>[2] sigma;
  
  // Migration process
  vector[npars_mig] beta_mig;
  real log_sigma_mig;
  
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
  vector[4 * 5] alpha_age_diag; // 4 age groups x 5 states
  
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
  vector<lower=0>[nquar * nage] H = exp(X * beta);
  
  // Tensor product prior precision
  vector<lower=0>[2] lambda = inv(square(sigma));
  matrix[npars_inf - 1, npars_inf - 1] K1_inv;
  K1_inv = S1[1 : (npars_inf - 1), 1 : (npars_inf - 1)] * lambda[1]
           + S1[1 : (npars_inf - 1), npars_inf : (2 * npars_inf - 2)]
             * lambda[2];
  
  // Migration process
  real<lower=0> sigma_mig = exp(log_sigma_mig);
  vector[nquar] gamma = X_mig * beta_mig;
  vector<lower=0>[nquar] omega = exp(gamma);
  
  // Age distribution used for allocating migrant flows by age
  matrix<lower=0, upper=1>[nquar, nage] abroad_age_prop;
  for (t in 1 : nquar) {
    real total = ABROAD[t];
    if (total > 0) {
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
  vector[nquar] delta_5 = alpha_diag_5 + 0.7 * (X_diag * beta_diag_5);
  
  // Diagnosis probabilities by age group (ag), state (k), and time (t)
  array[4] matrix<lower=0, upper=1>[5, nquar] d;
  for (ag in 1 : 4) {
    int ag_offset = (ag - 1) * 5;
    d[ag, 1,  : ] = to_row_vector(
                                  inv_logit(
                                            delta_1
                                            + alpha_age_diag[ag_offset + 1]));
    d[ag, 2,  : ] = to_row_vector(
                                  inv_logit(
                                            delta_2
                                            + alpha_age_diag[ag_offset + 2]));
    d[ag, 3,  : ] = to_row_vector(
                                  inv_logit(
                                            delta_3
                                            + alpha_age_diag[ag_offset + 3]));
    d[ag, 4,  : ] = to_row_vector(
                                  inv_logit(
                                            delta_4
                                            + alpha_age_diag[ag_offset + 4]));
    d[ag, 5,  : ] = to_row_vector(
                                  inv_logit(
                                            delta_5
                                            + alpha_age_diag[ag_offset + 5]));
  }
  
  // Proportion UK-acquired among migrants
  vector[nquar] eta = alpha_cob + 0.5 * (X_cob * beta_cob);
  vector<lower=0>[nquar] p;
  
  // Migration probabilities into compartments
  vector[nquar] theta_1 = alpha_mp[2] + X_mp * beta_mp_1;
  vector[nquar] theta_2 = alpha_mp[3] + X_mp * beta_mp_2;
  vector[nquar] theta_3 = alpha_mp[4] + X_mp * beta_mp_3;
  vector[nquar] theta_5 = alpha_mp[1] + X_mp * beta_mp_5;
  matrix<lower=0>[5, nquar] m;
  
  for (t in 1 : nquar) {
    real remainder = 1.0;
    m[5, t] = inv_logit(theta_5[t]) * remainder;
    remainder -= m[5, t];
    m[1, t] = inv_logit(theta_1[t]) * remainder;
    remainder -= m[1, t];
    m[2, t] = inv_logit(theta_2[t]) * remainder;
    remainder -= m[2, t];
    m[3, t] = inv_logit(theta_3[t]) * remainder;
    remainder -= m[3, t];
    m[4, t] = remainder;
  }
  
  // Misclassification model for abroad diagnoses
  vector[nquar] detect_logit = alpha_detect + 0.5 * (X_mp * beta_detect);
  vector<lower=0, upper=1>[nquar] detect_prob;
  vector<lower=0>[nquar] true_abroad_dx;
  vector<lower=0>[nquar] exp_ABROAD_dx;
  matrix<lower=0>[nquar, nage] misclassified_dx_age = rep_matrix(0, nquar,
                                                                 nage);
  
  for (t in 1 : nquar) {
    detect_prob[t] = inv_logit(detect_logit[t]);
    true_abroad_dx[t] = omega[t] * m[5, t];
    exp_ABROAD_dx[t] = true_abroad_dx[t] * detect_prob[t];
    for (a in 1 : nage) {
      misclassified_dx_age[t, a] = true_abroad_dx[t] * (1 - detect_prob[t])
                                   * abroad_age_prop[t, a];
    }
  }
  
  // Expected diagnoses by (t, a, y) and diagnosed state (1..5 HIV, 6 AIDS)
  array[nquar, nage, 2] row_vector<lower=0>[6] exp_diag;
  array[nquar, nage, 2] real<lower=0> exp_HIV_dx;
  array[nquar, nage, 2] real<lower=0> exp_AIDS_dx;
  array[nquar, nage, 2] simplex[5] exp_p;
  
  for (t in 1 : nquar) {
    for (a in 1 : nage) {
      for (y in 1 : 2) {
        exp_diag[t, a, y] = empty_vec;
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
        
        // Prevalent infections at t0 == 1
        if (t0 == 1) {
          int ag0 = age_group_lookup[a0];
          
          real d1 = d[ag0, 1, 1];
          real d2 = d[ag0, 2, 1];
          real d3 = d[ag0, 3, 1];
          real d4 = d[ag0, 4, 1];
          real d5 = d[ag0, 5, 1];
          real nd1 = 1 - d1;
          real nd2 = 1 - d2;
          real nd3 = 1 - d3;
          real nd4 = 1 - d4;
          real nd5 = 1 - d5;
          
          real q1 = q_age[1, a0];
          real q2 = q_age[2, a0];
          real q3 = q_age[3, a0];
          real q4 = q_age[4, a0];
          real nq1 = 1 - q1;
          real nq2 = 1 - q2;
          real nq3 = 1 - q3;
          real nq4 = 1 - q4;
          
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
        real inflow_native = H[(a0 - 1) * nquar + t0];
        real inflow_mig_uk = inflow_native * p[t0];
        
        // Migrants acquired abroad enter states 4..7 only (no recent states 1..3)
        vector[7] inflow_mig_ab = rep_vector(0, 7);
        inflow_mig_ab[4] = omega[t0] * abroad_age_prop[t0, a0] * m[1, t0];
        inflow_mig_ab[5] = omega[t0] * abroad_age_prop[t0, a0] * m[2, t0];
        inflow_mig_ab[6] = omega[t0] * abroad_age_prop[t0, a0] * m[3, t0];
        inflow_mig_ab[7] = omega[t0] * abroad_age_prop[t0, a0] * m[4, t0];
        
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
            int a = min(a0 + year_index[t] - year_index[t0], nage);
            int ag = age_group_lookup[a];
            
            real d1 = d[ag, 1, t];
            real d2 = d[ag, 2, t];
            real d3 = d[ag, 3, t];
            real d4 = d[ag, 4, t];
            real d5 = d[ag, 5, t];
            real nd1 = 1 - d1;
            real nd2 = 1 - d2;
            real nd3 = 1 - d3;
            real nd4 = 1 - d4;
            real nd5 = 1 - d5;
            
            // Progression depends on age at entry a0
            real q1 = q_age[1, a0];
            real q2 = q_age[2, a0];
            real q3 = q_age[3, a0];
            real q4 = q_age[4, a0];
            real nq1 = 1 - q1;
            real nq2 = 1 - q2;
            real nq3 = 1 - q3;
            real nq4 = 1 - q4;
            
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
          hiv_dx += misclassified_dx_age[t, a];
        }
        
        exp_HIV_dx[t, a, y] = hiv_dx;
        exp_AIDS_dx[t, a, y] = aids_dx;
        
        {
          real denom = hiv_dx;
          if (y == 2) {
            denom -= misclassified_dx_age[t, a];
          }
          if (denom > 0) {
            exp_p[t, a, y] = to_vector(exp_diag[t, a, y, 1 : 5]) / denom;
          } else {
            exp_p[t, a, y] = rep_vector(0.2, 5);
          }
        }
      }
    }
  }
}
model {
  // Infection surface priors
  sigma ~ student_t(4, 0, 200);
  beta[1] ~ normal(0, 30);
  beta[2 : npars_inf] ~ multi_normal_prec(zero_vec, K1_inv);
  
  // Migration priors
  beta_mig[1] ~ normal(4, 1);
  beta_mig[2 : npars_mig] ~ normal(0, sigma_mig);
  log_sigma_mig ~ normal(-2, 1);
  
  // Diagnosis priors
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
  
  sigma_diag_1 ~ gamma(2, 16);
  sigma_diag_2 ~ gamma(2, 16);
  sigma_diag_3 ~ gamma(2, 16);
  sigma_diag_4 ~ gamma(2, 16);
  sigma_diag_5 ~ gamma(2, 16);
  
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
