// Age-dependent tensor product model with RITA state

functions {
  vector row_sums(matrix X) {
    vector[rows(X)] s;
    s = X * rep_vector(1, cols(X));
    return s;
  }
  
  row_vector col_sums(matrix X) {
    row_vector[cols(X)] s;
    s = rep_row_vector(1, rows(X)) * X;
    return s;
  }
  
  // Q_fct creates quarterly transition (only, no diag) matrices 
  // Assume diagnosis occurs first and if not then progression
  matrix Q_fct(matrix d, matrix q, int t, int a) {
    matrix[7, 7] Q = rep_matrix(0, 7, 7);
    real d1t = 1 - d[1, t];
    real d2t = 1 - d[2, t];
    real d3t = 1 - d[3, t];
    real d4t = 1 - d[4, t];
    real d5t = 1 - d[5, t];
    real q1a = 1 - q[1, a];
    real q2a = 1 - q[2, a];
    real q3a = 1 - q[3, a];
    real q4a = 1 - q[4, a];
    
    Q[1, 2] = d1t * q1a;
    Q[1, 3] = d1t * q[1, a];
    Q[2, 4] = d1t * q1a;
    Q[2, 5] = d1t * q[1, a];
    Q[3, 5] = d1t * q2a;
    Q[3, 6] = d1t * q[2, a];
    Q[4, 4] = d2t * q1a;
    Q[4, 5] = d2t * q[1, a];
    Q[5, 5] = d3t * q2a;
    Q[5, 6] = d3t * q[2, a];
    Q[6, 6] = d4t * q3a;
    Q[6, 7] = d4t * q[3, a];
    Q[7, 7] = d5t * q4a;
    
    return Q;
  }
  
  // Links age to age group
  // age 1-10 -> ag1
  // age 11-10 -> ag2
  // age 21-30 -> ag3
  // age 31-52 -> ag4
  int age_to_agroup(int a) {
    int ag;
    if (a <= 10) {
      ag = 1;
    } else if (a > 10 && a <= 20) {
      ag = 2;
    } else if (a > 20 && a <= 30) {
      ag = 3;
    } else {
      ag = 4;
    }
    return ag;
  }
  
  array[,] matrix diag_arr_fct(array[,] matrix diag_arr, matrix e_m,
                               int nquar, int nage, array[,,] int age_lookup) {
    array[nquar, nage] matrix[nquar, 5] diag_arr1;
    int a;
    
    // Initialize all entries to 0
    for (t0 in 1 : nquar) {
      for (a0 in 1 : nage) {
        diag_arr1[t0, a0] = e_m;
      }
    }
    
    // Accumulate diagnoses at correct ages - using lookup
    for (t0 in 1 : nquar) {
      for (a0 in 1 : nage) {
        for (t in t0 : nquar) {
          a = age_lookup[a0, t0, t];
          diag_arr1[t, a, t0,  : ] += diag_arr[t, a0, t0,  : ];
        }
      }
    }
    return diag_arr1;
  }
}
data {
  int<lower=1> nquar; // number of years 
  int<lower=1> nage; // number of ages 
  array[nquar, nage] int<lower=0> HIV_age; // data with new HIV_age diag by year and age
  array[nquar, nage] int<lower=0> AIDS_age; // data with new AIDS_age diag by year and age
  array[nquar, nage, 5] int<lower=0> CD4_age; // CD4_age data by number of years and ages and CD4_age states(3)
  matrix<lower=0, upper=1>[4, nage] q_age; // matrix of progressions for different ages.. fixed
  matrix<lower=0>[7, nage] init_prev_age; // Initial prevalence by age at inf and CD4_age states(3) 
  
  // Splines
  int<lower=1> npars_inf;
  matrix[nquar * nage, npars_inf] X; // design matrix ... 50 is number of inf params
  matrix[npars_inf - 1, 2 * npars_inf - 2] S1; // NOT SURE IF BOUNDARY CONSISTENT
}
transformed data {
  vector[6] empty_vec = rep_vector(0, 6);
  row_vector[6] empty_vec1 = rep_row_vector(0, 6);
  vector[7] empty_vec2 = rep_vector(0, 7);
  matrix[nquar, 6] empty_mat = rep_matrix(0, nquar, 6);
  vector[npars_inf - 1] zero_vec = rep_vector(0, npars_inf - 1);
  int nagroup = 4; // number of age groups
  
  // NEW: PRE-COMPUTE AGE MAPPINGS
  array[nage, nquar, nquar] int age_lookup;
  array[nage] int age_group_lookup;
  
  // Pre-compute all age transitions
  for (a0 in 1 : nage) {
    age_group_lookup[a0] = age_to_agroup(a0);
    
    for (t0 in 1 : nquar) {
      for (t in 1 : nquar) {
        if (t < t0) {
          age_lookup[a0, t0, t] = a0; // dummy value
        } else {
          // Compute age at time t given infection at t0 at age a0
          int a = a0;
          for (t1 in (t0 + 1) : t) {
            if (t1 % 4 == 1) {
              a = a + 1;
            }
          }
          age_lookup[a0, t0, t] = min(a, nage);
        }
      }
    }
  }
}
parameters {
  // Infection parameters
  vector[npars_inf] beta; //beta are the infection parameters
  vector<lower=0>[2] sigma; // smoothing parameters
  vector<lower=0>[5] vardelta;
  
  // Diagnosis probability parameters
  vector[nquar] delta_raw_1;
  vector[nquar] delta_raw_2;
  vector[nquar] delta_raw_3;
  vector[nquar] delta_raw_4;
  vector[nquar] delta_raw_5;
  vector[5 * nagroup] alpha; // age group specific (ag) and state (k) specific linear time trend for dx probs --- 
}
transformed parameters {
  // Infection parameters
  vector<lower=0>[nquar * nage] H = exp(X * beta); //vector of expected log responses
  array[nquar] vector<lower=0>[7] lat_arr; // expected number of diag at time t age a state k
  vector<lower=0>[7] lat_p;
  
  // Diagnosis probabilities
  row_vector<lower=0>[6] diag_p;
  array[nquar, nage] matrix[nquar, 6] diag_arr; // expected number of diag at time t age a state k
  array[nquar, nage] matrix<lower=0>[nquar, 6] diag_arr1; // expected number of diag at time t age a state k
  vector<lower=0>[7] temp; // expected number of diag at time t age a state k
  array[nquar, nage] row_vector<lower=0>[6] exp_diag; // expected number of diag at time t age a state k
  
  // Expected observations
  matrix[nquar, nage] exp_HIV_dx; // expected number of HIV_age diag at time t age a
  matrix[nquar, nage] exp_AIDS_dx; // expected number of HIV_age diag at time t age a
  array[nquar, nage] simplex[5] exp_p; // proportion of HIV_age diagnoses in state k (time t and age a)
  array[4] matrix[5, nquar] d; // matrix of diagnoses for different ages and age-groups
  matrix[npars_inf - 1, npars_inf - 1] K1_inv;
  
  // Expected incidence
  vector[nquar] delta_1;
  vector[nquar] delta_2;
  vector[nquar] delta_3;
  vector[nquar] delta_4;
  vector[nquar] delta_5;
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
  delta_3[1] = -3.2 + 0.2 * delta_raw_3[1];
  delta_4[1] = -3.0 + 0.2 * delta_raw_4[1];
  delta_5[1] = -2.5 + 0.3 * delta_raw_5[1];
  
  // Random walk for diagnoses
  for (t in 2 : nquar) {
    delta_1[t] = delta_1[t - 1] + vardelta[1] * delta_raw_1[t];
    delta_2[t] = delta_2[t - 1] + vardelta[2] * delta_raw_2[t];
    delta_3[t] = delta_3[t - 1] + vardelta[3] * delta_raw_3[t];
    delta_4[t] = delta_4[t - 1] + vardelta[4] * delta_raw_4[t];
    delta_5[t] = delta_5[t - 1] + vardelta[5] * delta_raw_5[t];
  }
  
  // diag probs are inv logit of delta
  // so that are normalised btwn 0 and 1
  // Age-specific additive term, interpreted as a linear time trend (for a given state)
  // which has an age and time specific gradient
  for (ag in 1 : nagroup) {
    int ag_offset = (ag - 1) * 5;
    d[ag, 1,  : ] = to_row_vector(inv_logit(delta_1 + alpha[ag_offset + 1]));
    d[ag, 2,  : ] = to_row_vector(inv_logit(delta_2 + alpha[ag_offset + 2]));
    d[ag, 3,  : ] = to_row_vector(inv_logit(delta_3 + alpha[ag_offset + 3]));
    d[ag, 4,  : ] = to_row_vector(inv_logit(delta_4 + alpha[ag_offset + 4]));
    d[ag, 5,  : ] = to_row_vector(inv_logit(delta_5 + alpha[ag_offset + 5]));
  }
  
  // NOW NEED TO CREATE THE TRANSITION MATRICES AND CALCULATE EXPECTED DIAGNOSES
  
  for (t0 in 1 : nquar) {
    for (a0 in 1 : nage) {
      // At time = time of inf
      // New infs enter model in state 1
      
      // At time 1 we need to consider initial prevalence
      // Prevalence is latent in 1994Q4 (t = 0)
      // So consider diagnosis and progressions
      lat_p = empty_vec2;
      diag_p = empty_vec1;
      
      if (t0 == 1) {
        int ag0 = age_group_lookup[a0];
        matrix[7, 7] PA_current = Q_fct(d[ag0,  : ,  : ], q_age, t0, a0);
        lat_p = (PA_current' * init_prev_age[ : , a0]);
        diag_p[1] = (init_prev_age[1, a0] * d[ag0, 1, 1])
                    + (init_prev_age[2, a0] * d[ag0, 1, 1])
                    + (init_prev_age[3, a0] * d[ag0, 1, 1]);
        diag_p[2] = init_prev_age[4, a0] * d[ag0, 2, 1];
        diag_p[3] = init_prev_age[5, a0] * d[ag0, 3, 1];
        diag_p[4] = init_prev_age[6, a0] * d[ag0, 4, 1];
        diag_p[5] = init_prev_age[7, a0] * d[ag0, 5, 1];
        diag_p[6] = init_prev_age[7, a0] * (1 - d[ag0, 5, 1]) * q_age[4, a0];
      }
      lat_arr[t0, 1] = H[(a0 - 1) * nquar + t0] + lat_p[1];
      lat_arr[t0, 2 : 7] = empty_vec + lat_p[2 : 7];
      diag_arr[t0, a0, t0,  : ] = empty_vec1 + diag_p;
      
      for (t in (t0 + 1) : nquar) {
        // might worry if t = nquar, but in STAN loops in (nquar+1):nquar are not evaluated (unlike R, like JAGS) --- saving an if loop
        int a = age_lookup[a0, t0, t];
        int ag = age_group_lookup[a];
        
        matrix[7, 7] PA_current = Q_fct(d[ag,  : ,  : ], q_age, t, a);
        temp = lat_arr[t - 1,  : ];
        
        lat_arr[t,  : ] = (PA_current' * temp);
        // mapping a0,t and t0 to corresponding current age via age_lookup
        diag_arr[t, a0, t0, 1] = (temp[1] * d[ag, 1, t])
                                 + (temp[2] * d[ag, 1, t])
                                 + (temp[3] * d[ag, 1, t]);
        diag_arr[t, a0, t0, 2] = temp[4] * d[ag, 2, t];
        diag_arr[t, a0, t0, 3] = temp[5] * d[ag, 3, t];
        diag_arr[t, a0, t0, 4] = temp[6] * d[ag, 4, t];
        diag_arr[t, a0, t0, 5] = temp[7] * d[ag, 5, t];
        diag_arr[t, a0, t0, 6] = temp[7] * (1 - d[ag, 5, t]) * q_age[4, a0];
      }
    }
  }
  
  diag_arr1 = diag_arr_fct(diag_arr, empty_mat, nquar, nage, age_lookup);
  
  for (t in 1 : nquar) {
    for (a in 1 : nage) {
      exp_diag[t, a,  : ] = col_sums(diag_arr1[t, a,  : ,  : ]);
      // expected HIV diagnoses
      exp_HIV_dx[t, a] = sum(exp_diag[t, a, 1 : 5]);
      // Expected AIDS diagnoses
      exp_AIDS_dx[t, a] = exp_diag[t, a, 6];
    }
  }
  
  /// For a = 1 and some times (t=5, t=9 ...)
  // there are 0 diagnes arrivals
  // This is due to model assumptions
  // So removing age 1 from likelihood
  // Fixing a=1 to get a valid simplex
  for (t in 1 : nquar) {
    exp_p[t, 1, 1] = 1;
    exp_p[t, 1, 2 : 5] = rep_vector(0, 4);
  }
  
  for (t in 1 : nquar) {
    for (a in 2 : nage) {
      // Now we want the proportion of CD4_age that
      // is diangosed with HIV_age in state k at time t and age a (p_kta)
      real inv_total = inv(exp_HIV_dx[t, a]);
      exp_p[t, a, 1] = exp_diag[t, a, 1] * inv_total;
      exp_p[t, a, 2] = exp_diag[t, a, 2] * inv_total;
      exp_p[t, a, 3] = exp_diag[t, a, 3] * inv_total;
      exp_p[t, a, 4] = exp_diag[t, a, 4] * inv_total;
      exp_p[t, a, 5] = exp_diag[t, a, 5] * inv_total;
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
  delta_raw_5 ~ normal(0, 1);
  
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
    log_lik[(1 - 1) * nquar + t] = poisson_lpmf(HIV_age[t, 1] |
                                                exp_HIV_dx[t, 1]);
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
