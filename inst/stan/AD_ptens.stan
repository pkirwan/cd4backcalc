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
  matrix[nquar * nage, npars_inf] X; // design matrix ... 50 is number of inf params
  matrix[npars_inf - 1, 2 * npars_inf - 2] S1; // NOT SURE IF BOUNDARY CONSISTENT
}
transformed data {
  row_vector[5] empty_vec = rep_row_vector(0, 5);
  vector[4] empty_vec1 = rep_vector(0, 4);
  vector[npars_inf - 1] zero_vec = rep_vector(0, npars_inf - 1);
  int nagroup = 4; // number of age groups
  
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
}
parameters {
  // Infection parameters
  vector[npars_inf] beta; //beta are the infection parameters
  vector<lower=0>[2] sigma; // smoothing parameters
  vector<lower=0>[4] vardelta;
  
  // Diagnosis probability parameters
  vector[nquar] delta_raw_1;
  vector[nquar] delta_raw_2;
  vector[nquar] delta_raw_3;
  vector[nquar] delta_raw_4;
  vector[4 * nagroup] alpha; // age group specific (ag) and state (k) specific linear time trend for dx probs --- 
}
transformed parameters {
  // Infection parameters
  vector<lower=0>[nquar * nage] H = exp(X * beta); //vector of expected log responses
  vector<lower=0>[4] lat_p;
  
  // Diagnosis probabilities
  row_vector<lower=0>[5] diag_p;
  array[nquar, nage] row_vector<lower=0>[5] exp_diag; // expected number of diag at time t age a state k
  
  // Expected observations
  matrix[nquar, nage] exp_HIV_dx; // expected number of HIV_age diag at time t age a
  matrix[nquar, nage] exp_AIDS_dx; // expected number of HIV_age diag at time t age a
  array[nquar, nage] simplex[4] exp_p; // proportion of HIV_age diagnoses in state k (time t and age a)
  array[4] matrix[4, nquar] d; // matrix of diagnoses for different ages and age-groups
  matrix[npars_inf - 1, npars_inf - 1] K1_inv;
  
  // Expected incidence
  vector[nquar] delta_1;
  vector[nquar] delta_2;
  vector[nquar] delta_3;
  vector[nquar] delta_4;
  vector<lower=0>[2] lambda; // smoothing parameters
  
  lambda = inv(square(sigma));
  
  // PRIORS for the tensor product
  K1_inv = S1[1 : (npars_inf - 1), 1 : (npars_inf - 1)] * lambda[1]
           + S1[1 : (npars_inf - 1), npars_inf : (2 * npars_inf - 2)]
             * lambda[2];
  
  // DIAGNOSES
  
  // Initializing delta
  delta_1[1] = -3.2 + 0.2 * delta_raw_1[1];
  delta_2[1] = -3.2 + 0.2 * delta_raw_2[1];
  delta_3[1] = -3.0 + 0.2 * delta_raw_3[1];
  delta_4[1] = -2.5 + 0.3 * delta_raw_4[1];
  
  // Random walk for diagnoses
  for (t in 2 : nquar) {
    delta_1[t] = delta_1[t - 1] + vardelta[1] * delta_raw_1[t];
    delta_2[t] = delta_2[t - 1] + vardelta[2] * delta_raw_2[t];
    delta_3[t] = delta_3[t - 1] + vardelta[3] * delta_raw_3[t];
    delta_4[t] = delta_4[t - 1] + vardelta[4] * delta_raw_4[t];
  }
  
  // diag probs are inv logit of delta
  // so that are normalised btwn 0 and 1
  // Age-specific additive term, interpreted as a linear time trend (for a given state)
  // which has an age and time specific gradient
  for (ag in 1 : nagroup) {
    int ag_offset = (ag - 1) * 4;
    d[ag, 1,  : ] = to_row_vector(inv_logit(delta_1 + alpha[ag_offset + 1]));
    d[ag, 2,  : ] = to_row_vector(inv_logit(delta_2 + alpha[ag_offset + 2]));
    d[ag, 3,  : ] = to_row_vector(inv_logit(delta_3 + alpha[ag_offset + 3]));
    d[ag, 4,  : ] = to_row_vector(inv_logit(delta_4 + alpha[ag_offset + 4]));
  }
  
  // NOW NEED TO CREATE THE TRANSITION MATRICES AND CALCULATE EXPECTED DIAGNOSES
  
  for (t in 1 : nquar) {
    for (a in 1 : nage) {
      exp_diag[t, a] = empty_vec;
    }
  }
  
  for (t0 in 1 : nquar) {
    for (a0 in 1 : nage) {
      vector[4] lat;
      vector[4] lat_next;

      lat_p = empty_vec1;
      diag_p = empty_vec;
      
      if (t0 == 1) {
        int ag0 = age_group_lookup[a0];
        
        real d1 = d[ag0, 1, 1];
        real d2 = d[ag0, 2, 1];
        real d3 = d[ag0, 3, 1];
        real d4 = d[ag0, 4, 1];
        real nd1 = 1 - d1;
        real nd2 = 1 - d2;
        real nd3 = 1 - d3;
        real nd4 = 1 - d4;
        real q1 = q_age[1, a0];
        real q2 = q_age[2, a0];
        real q3 = q_age[3, a0];
        real q4 = q_age[4, a0];
        real nq1 = 1 - q1;
        real nq2 = 1 - q2;
        real nq3 = 1 - q3;
        real nq4 = 1 - q4;
        
        real v1 = init_prev_age[1, a0];
        real v2 = init_prev_age[2, a0];
        real v3 = init_prev_age[3, a0];
        real v4 = init_prev_age[4, a0];
        
        // diagnoses from prevalent infections at quarter 1
        diag_p[1] = v1 * d1;
        diag_p[2] = v2 * d2;
        diag_p[3] = v3 * d3;
        diag_p[4] = v4 * d4;
        diag_p[5] = v4 * nd4 * q4;
        
        // survivors after diagnosis and progression
        lat_p[1] = (nd1 * nq1) * v1;
        lat_p[2] = (nd1 * q1) * v1 + (nd2 * nq2) * v2;
        lat_p[3] = (nd2 * q2) * v2 + (nd3 * nq3) * v3;
        lat_p[4] = (nd3 * q3) * v3 + (nd4 * nq4) * v4;
        
        // diagnoses at t0 occur at age a0
        exp_diag[t0, a0] += diag_p;
      }
      
      // latent state immediately after t0 inflow (+ any t0==1 prevalence survivors)
      lat[1] = H[(a0 - 1) * nquar + t0] + lat_p[1];
      lat[2] = lat_p[2];
      lat[3] = lat_p[3];
      lat[4] = lat_p[4];
      
      if (t0 < nquar) {
        for (t in (t0 + 1) : nquar) {
          
          int a = min(a0 + year_index[t] - year_index[t0], nage);
          int ag = age_group_lookup[a];
          
          // Diagnosis probabilities based on current age group
          real d1 = d[ag, 1, t];
          real d2 = d[ag, 2, t];
          real d3 = d[ag, 3, t];
          real d4 = d[ag, 4, t];
          real nd1 = 1 - d1;
          real nd2 = 1 - d2;
          real nd3 = 1 - d3;
          real nd4 = 1 - d4;
          
          // Progression is based on age at infection
          real q1 = q_age[1, a0];
          real q2 = q_age[2, a0];
          real q3 = q_age[3, a0];
          real q4 = q_age[4, a0];
          real nq1 = 1 - q1;
          real nq2 = 1 - q2;
          real nq3 = 1 - q3;
          real nq4 = 1 - q4;
          
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
      if (total > 0) {
        exp_p[t, a] = to_vector(exp_diag[t, a, 1 : 4]) / total;
      } else {
        // safeguard; should be rare outside a==1, but avoids inv(0)
        exp_p[t, a] = rep_vector(0.25, 4);
      }
    }
  }
}
model {
  // Priors
  // Splines
  // sigma ~ inv_gamma(.05,.005); // s.d of observations. Manual recommends using flat priors, following Gelman's paper
  // If commented default un-normslized uniform prior
  
  sigma ~ student_t(4, 0, 200); // Weakly informative priors
  
  beta[1] ~ normal(0, 30); // 30 s.d corresponding to precision given in JAGS
  beta[2 : npars_inf] ~ multi_normal_prec(zero_vec, K1_inv);
  
  vardelta ~ gamma(1, 8); // prior for logit random walk
  
  alpha ~ normal(0, 1); // Priors for age specific intercepts
  
  // Non centered parametrization for diagnoses
  delta_raw_1 ~ normal(0, 1);
  delta_raw_2 ~ normal(0, 1);
  delta_raw_3 ~ normal(0, 1);
  delta_raw_4 ~ normal(0, 1);
  
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
