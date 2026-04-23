// Age-dependent tensor product model

data {
  int<lower=1> nquar; // number of years 
  int<lower=1> nage; // number of ages 
  array[nquar, nage] int<lower=0> HIV_age; // data with new HIV_age diag by year and age
  array[nquar, nage] int<lower=0> AIDS_age; // data with new AIDS_age diag by year and age
  array[nquar, nage, 4] int<lower=0> CD4_age; // CD4 data by number of years and ages and CD4 states(4)
  matrix<lower=0, upper=1>[4, nage] q_age; // matrix of progressions for different ages.. fixed
  matrix<lower=0>[4, nage] init_prev_age; // Initial prevalence by age at inf and CD4 states(4) 
  
  // Splines
  int<lower=1> npars_inf;
  matrix[nquar * nage, npars_inf] X_inf; // design matrix ... 50 is number of inf params
  matrix[npars_inf - 1, 2 * npars_inf - 2] S1; // NOT SURE IF BOUNDARY CONSISTENT
  
  // Time-varying diagnosis splines, plus age-group offsets
  int<lower=1> npars_diag;
  matrix[nquar, npars_diag] X_diag;
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
}
parameters {
  // Infection parameters
  vector[npars_inf] beta_inf; //beta_inf are the infection parameters
  vector<lower=0>[2] sigma_inf; // smoothing parameters
  
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
  vector[4 * nagroup] alpha_age_diag; // 4 age groups x 4 states
}
transformed parameters {
  // Infection parameters
  vector[nquar * nage] lambda = X_inf * beta_inf; // expected infections by time and age
  vector<lower=0>[nquar * nage] h = exp(lambda); //vector of expected log responses
  
  // Tensor product prior precision
  vector<lower=0>[2] tau = inv(square(sigma_inf));
  
  // Diagnosis probabilities by time
  vector[nquar] delta_1 = alpha_diag_1 + 0.5 * (X_diag * beta_diag_1);
  vector[nquar] delta_2 = alpha_diag_2 + 0.5 * (X_diag * beta_diag_2);
  vector[nquar] delta_3 = alpha_diag_3 + 0.5 * (X_diag * beta_diag_3);
  vector[nquar] delta_4 = alpha_diag_4 + 0.7 * (X_diag * beta_diag_4);
  
  // Diagnosis probabilities by age group (ag), state (k), and time (t).
  // Stored as matrix[nquar, 4] so the inner time-loop accesses column k,
  // which is contiguous in column-major memory (cache-friendly).
  array[4] matrix<lower=0, upper=1>[nquar, 4] d;
  for (ag in 1 : nagroup) {
    int ag_offset = (ag - 1) * 4;
    d[ag,  : , 1] = inv_logit(delta_1 + alpha_age_diag[ag_offset + 1]);
    d[ag,  : , 2] = inv_logit(delta_2 + alpha_age_diag[ag_offset + 2]);
    d[ag,  : , 3] = inv_logit(delta_3 + alpha_age_diag[ag_offset + 3]);
    d[ag,  : , 4] = inv_logit(delta_4 + alpha_age_diag[ag_offset + 4]);
  }
  
  // Expected diagnoses
  array[nquar, nage] row_vector<lower=0>[5] exp_diag; // expected number of diag at time t age a state k
  
  // Expected observations
  matrix[nquar, nage] exp_HIV_dx; // expected number of HIV_age diag at time t age a
  matrix[nquar, nage] exp_AIDS_dx; // expected number of HIV_age diag at time t age a
  array[nquar, nage] simplex[4] exp_p; // proportion of HIV_age diagnoses in state k (time t and age a)
  
  // NOW NEED TO CREATE THE TRANSITION MATRICES AND CALCULATE EXPECTED DIAGNOSES
  
  exp_diag = rep_array(empty_vec, nquar, nage);
  
  for (t0 in 1 : nquar) {
    for (a0 in 1 : nage) {
      // Progression rates (depend only on age at infection a0, not current time)
      real q1 = q_age[1, a0];
      real q2 = q_age[2, a0];
      real q3 = q_age[3, a0];
      real q4 = q_age[4, a0];
      real nq1 = 1 - q1;
      real nq2 = 1 - q2;
      real nq3 = 1 - q3;
      real nq4 = 1 - q4;
      
      vector[4] lat = empty_vec1;
      vector[4] lat_next;
      
      if (t0 == 1) {
        int ag0 = age_group_lookup[a0];
        
        real d1 = d[ag0, 1, 1];
          real d2 = d[ag0, 1, 2];
          real d3 = d[ag0, 1, 3];
          real d4 = d[ag0, 1, 4];
        real nd1 = 1 - d1;
        real nd2 = 1 - d2;
        real nd3 = 1 - d3;
        real nd4 = 1 - d4;
        
        real v1 = init_prev_age[1, a0];
        real v2 = init_prev_age[2, a0];
        real v3 = init_prev_age[3, a0];
        real v4 = init_prev_age[4, a0];
        
        // diagnoses from prevalent infections at quarter 1
        exp_diag[1, a0, 1] = v1 * d1;
        exp_diag[1, a0, 2] = v2 * d2;
        exp_diag[1, a0, 3] = v3 * d3;
        exp_diag[1, a0, 4] = v4 * d4;
        exp_diag[1, a0, 5] = v4 * nd4 * q4;
        
        // survivors after diagnosis and progression
        lat[1] = (nd1 * nq1) * v1;
        lat[2] = (nd1 * q1) * v1 + (nd2 * nq2) * v2;
        lat[3] = (nd2 * q2) * v2 + (nd3 * nq3) * v3;
        lat[4] = (nd3 * q3) * v3 + (nd4 * nq4) * v4;
      }
      
      // latent state immediately after t0 inflow (+ any t0==1 prevalence survivors)
      lat[1] += h[(a0 - 1) * nquar + t0];
      
      for (t in (t0 + 1) : nquar) {
        int a = min(a0 + year_diff[t0, t], nage);
        int ag = age_group_lookup[a];
        
        // Diagnosis probabilities based on current age group
        real d1 = d[ag, t, 1];
        real d2 = d[ag, t, 2];
        real d3 = d[ag, t, 3];
        real d4 = d[ag, t, 4];
        
        real nd1 = 1 - d1;
        real nd2 = 1 - d2;
        real nd3 = 1 - d3;
        real nd4 = 1 - d4;
        
        // diagnosis uses current latent state (before progression)
        exp_diag[t, a, 1] += lat[1] * d1;
        exp_diag[t, a, 2] += lat[2] * d2;
        exp_diag[t, a, 3] += lat[3] * d3;
        exp_diag[t, a, 4] += lat[4] * d4;
        exp_diag[t, a, 5] += lat[4] * nd4 * q4;
        
        // progression among undiagnosed
        lat_next[1] = (nd1 * nq1) * lat[1];
        lat_next[2] = (nd1 * q1) * lat[1] + (nd2 * nq2) * lat[2];
        lat_next[3] = (nd2 * q2) * lat[2] + (nd3 * nq3) * lat[3];
        lat_next[4] = (nd3 * q3) * lat[3] + (nd4 * nq4) * lat[4];
        
        lat = lat_next;
      }
    }
  }
  
  // Expected HIV and AIDS diagnoses by (t,a)
  for (t in 1 : nquar) {
    for (a in 1 : nage) {
      exp_HIV_dx[t, a] = sum(exp_diag[t, a, 1 : 4]);
      exp_AIDS_dx[t, a] = exp_diag[t, a, 5];
    }
  }
  
  // Proportions for multinomial CD4 distribution
  for (t in 1 : nquar) {
    exp_p[t, 1, 1] = 1;
    exp_p[t, 1, 2 : 4] = rep_vector(0, 3);
  }
  
  for (t in 1 : nquar) {
    for (a in 2 : nage) {
      real total = exp_HIV_dx[t, a];
      if (total > 1e-6) {
        // Add small epsilon to avoid exact zeros in the simplex,
        // which cause NaN via 0*log(0) in the multinomial likelihood
        vector[4] raw_p = to_vector(exp_diag[t, a, 1 : 4]) + 1e-12;
        exp_p[t, a] = raw_p / sum(raw_p);
      } else {
        // safeguard; should be rare outside a==1, but avoids inv(0)
        exp_p[t, a] = rep_vector(0.25, 4);
      }
    }
  }
}
model {
  // Infection surface priors
  sigma_inf ~ normal(0, 50);
  beta_inf[1] ~ normal(0, 5);
  
  matrix[np1, np1] K1_prec = S1a * tau[1] + S1b * tau[2];
  beta_inf[2 : npars_inf] ~ multi_normal_prec(zero_vec, K1_prec);
  
  // Diagnosis priors (time spline + age group offsets)
  alpha_diag_1 ~ normal(-3.2, 0.5);
  alpha_diag_2 ~ normal(-3.0, 0.5);
  alpha_diag_3 ~ normal(-2.5, 0.5);
  alpha_diag_4 ~ normal(-2.0, 0.5);
  
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
  
  // Likelihood of model
  for (t in 1 : nquar) {
    // a is in 2 ages, because for age = 1 -> at some t get Po(0)
    // Individuals age at times 5,9,13...
    // So if individuals inf at t=5 at age=1
    // Cannot be dx at t=5 at age age=1 (dx t > 6)
    // remove age = 1 from likelihood
    HIV_age[t, 2 : nage] ~ poisson(exp_HIV_dx[t, 2 : nage]);
    AIDS_age[t, 2 : nage] ~ poisson(exp_AIDS_dx[t, 2 : nage]);
    for (a in 2 : nage) {
      CD4_age[t, a,  : ] ~ multinomial(exp_p[t, a,  : ]);
    }
  }
}
generated quantities {
  // Log-lik calculations to get WAIC / LOO/ DIC
  vector[nquar * nage] log_lik; // has to be vector to be consistent with R pckg
  for (t in 1 : nquar) {
    log_lik[(1 - 1) * nquar + t] = 0;
    
    for (a in 2 : nage) {
      log_lik[(a - 1) * nquar + t] = poisson_lpmf(HIV_age[t, a] |
                                                  exp_HIV_dx[t, a])
                                     + poisson_lpmf(AIDS_age[t, a] |
                                                    exp_AIDS_dx[t, a])
                                     + multinomial_lpmf(CD4_age[t, a,  : ] |
                                                        exp_p[t, a,  : ]);
    }
  }
}
