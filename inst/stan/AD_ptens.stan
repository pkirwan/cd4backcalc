// Age-dependent tensor product model

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
    matrix[4, 4] Q = rep_matrix(0, 4, 4);
    Q[1, 1] = (1 - d[1, t]) * (1 - q[1, a]);
    Q[1, 2] = (1 - d[1, t]) * q[1, a];
    
    Q[2, 2] = (1 - d[2, t]) * (1 - q[2, a]);
    Q[2, 3] = (1 - d[2, t]) * q[2, a];
    
    Q[3, 3] = (1 - d[3, t]) * (1 - q[3, a]);
    Q[3, 4] = (1 - d[3, t]) * q[3, a];
    
    Q[4, 4] = (1 - d[4, t]) * (1 - q[4, a]);
    return Q;
  }
  
  int age_fct(int a0, int t, int t0, int nquar, int nage) {
    // Assume individuals age at quarters 5,9,13...
    // So when modulus(t,4) == 1
    real rem;
    int a; // current age
    int addyr;
    
    a = a0; // at t0, current age =a0
    for (t1 in min(t0 + 1, nquar) : t) {
      rem = t1 % 4;
      if (rem == 1) 
        addyr = 1;
      else 
        addyr = 0;
      a = a + addyr;
    }
    
    return min(a, nage);
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
                               int nquar, int nage) {
    array[nquar, nage] matrix[nquar, 5] diag_arr1;
    array[nquar, nage] row_vector[5] out;
    row_vector[5] x;
    int a;
    
    // Setting all entries of diag_arr1 to 0
    // (default: no new arrival)
    // Needed: if not get nan
    for (t0 in 1 : nquar) {
      for (a0 in 1 : nage) {
        diag_arr1[t0, a0] = e_m;
      }
    }
    
    // Specific t0,a0,t for which new arrivals
    // are non zero
    for (t0 in 1 : nquar) {
      for (a0 in 1 : nage) {
        for (t in t0 : nquar) {
          a = age_fct(a0, t, t0, nquar, nage);
          x = diag_arr1[t, a, t0,  : ]; // To deal with last cumulative age
          diag_arr1[t, a, t0,  : ] = x + diag_arr[t, a0, t0,  : ];
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
  array[nquar, nage, 4] int<lower=0> CD4_age; // CD4_age data by number of years and ages and CD4_age states(3)
  matrix<lower=0, upper=1>[4, nage] q_age; // matrix of progressions for different ages.. fixed
  matrix<lower=0>[4, nage] init_prev_age; // Initial prevalence by age at inf and CD4_age states(3) 
  
  // Splines
  int<lower=1> npars_inf; // number of infection parameters
  matrix[nquar * nage, npars_inf] X; // design matrix ... 50 is number of inf params
  matrix[npars_inf - 1, 2 * npars_inf - 2] S1; // NOT SURE IF BOUNDARY CONSISTENT
}
transformed data {
  vector[3] empty_vec = rep_vector(0, 3);
  row_vector[5] empty_vec1 = rep_row_vector(0, 5);
  vector[4] empty_vec2 = rep_vector(0, 4);
  matrix[nquar, 5] empty_mat = rep_matrix(0, nquar, 5);
  vector[npars_inf - 1] zero_vec = rep_vector(0, npars_inf - 1);
  int nagroup = 4; // number of age groups
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
  array[nquar] vector<lower=0>[4] lat_arr; // expected number of diag at time t age a state k
  vector<lower=0>[4] lat_p;
  
  // Diagnosis probabilities
  row_vector<lower=0>[5] diag_p;
  array[nquar, nage] matrix[nquar, 5] diag_arr; // expected number of diag at time t age a state k
  array[nquar, nage] matrix<lower=0>[nquar, 5] diag_arr1; // expected number of diag at time t age a state k
  vector<lower=0>[4] temp; // expected number of diag at time t age a state k
  array[nquar, nage] row_vector<lower=0>[5] exp_diag; // expected number of diag at time t age a state k
  
  // Epidemic state variables
  array[nquar, nage, nagroup] matrix<lower=0, upper=1>[4, 4] PA; // transition matrices
  
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
  
  // INCIDENCE
  H = exp(X * beta); // These are log infection...
  
  // DIAGNOSES
  
  // Initializing delta
  delta_1[1] = -3.2 + 0.2 * delta_raw_1[1]; // worth reconsidering? Originally for diag in 85
  delta_2[1] = -3.2 + 0.2 * delta_raw_2[1]; // worth reconsidering? Originally for diag in 85
  delta_3[1] = -3.0 + 0.2 * delta_raw_3[1]; // worth reconsidering? Originally for diag in 85
  delta_4[1] = -2.5 + 0.3 * delta_raw_4[1]; // worth reconsidering? Originally for diag in 85
  
  // Random walk for diagnoses
  for (t in 2 : nquar) {
    delta_1[t] = delta_1[t - 1] + vardelta[1] * delta_raw_1[t];
    delta_2[t] = delta_2[t - 1] + vardelta[2] * delta_raw_2[t];
    delta_3[t] = delta_3[t - 1] + vardelta[3] * delta_raw_3[t];
    delta_4[t] = delta_4[t - 1] + vardelta[4] * delta_raw_4[t];
  }
  
  // diag probs are inv logit of delta
  // so that are normalised btwn 0 and 1
  // inv logit not vectorized -> use loop
  // Age-specific additive term, interpreted as a linear time trend (for a given state)
  // which has an age and time specific gradient
  for (ag in 1 : nagroup) {
    for (i in 1 : nquar) {
      d[ag, 1, i] = inv_logit(delta_1[i] + alpha[(ag - 1) * nagroup + 1]);
      d[ag, 2, i] = inv_logit(delta_2[i] + alpha[(ag - 1) * nagroup + 2]);
      d[ag, 3, i] = inv_logit(delta_3[i] + alpha[(ag - 1) * nagroup + 3]);
      d[ag, 4, i] = inv_logit(delta_4[i] + alpha[(ag - 1) * nagroup + 4]);
    }
  }
  
  // NOW NEED TO CREATE THE TRANSITION MATRICES
  
  // Useful to define a zero-matrix
  // For later calculations
  
  // Looping over t and a to get trans mat
  // for all ages and times
  for (t in 1 : nquar) {
    for (a in 1 : nage) {
      for (ag in 1 : nagroup) {
        //QUARTERLY progression probabilities to observed states
        PA[t, a, ag,  : ,  : ] = Q_fct(d[ag,  : ,  : ], q_age, t, a);
      }
    }
  }
  
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
        int ag0 = age_to_agroup(a0);
        lat_p = (PA[1, a0, ag0,  : ]' * init_prev_age[ : , a0]);
        diag_p[1] = init_prev_age[1, a0] * d[ag0, 1, 1];
        diag_p[2] = init_prev_age[2, a0] * d[ag0, 2, 1];
        diag_p[3] = init_prev_age[3, a0] * d[ag0, 3, 1];
        diag_p[4] = init_prev_age[4, a0] * d[ag0, 4, 1];
        diag_p[5] = init_prev_age[4, a0] * (1 - d[ag0, 4, 1]) * q_age[4, a0];
      }
      lat_arr[t0, 1] = H[(a0 - 1) * nquar + t0] + lat_p[1];
      lat_arr[t0, 2 : 4] = empty_vec + lat_p[2 : 4];
      diag_arr[t0, a0, t0,  : ] = empty_vec1 + diag_p;
      
      for (t in (t0 + 1) : nquar) {
        // might worry if t = nquar, but in STAN loops in (nquar+1):nquar are not evaluated (unlike R, like JAGS) --- saving an if loop
        int a = age_fct(a0, t, t0, nquar, nage);
        int ag = age_to_agroup(a);
        
        temp = lat_arr[t - 1,  : ];
        lat_arr[t,  : ] = (PA[t, a0, ag,  : ]' * temp);
        // mapping a0,t and t0 to corresponding current age via age_fct
        diag_arr[t, a0, t0, 1] = temp[1] * d[ag, 1, t];
        diag_arr[t, a0, t0, 2] = temp[2] * d[ag, 2, t];
        diag_arr[t, a0, t0, 3] = temp[3] * d[ag, 3, t];
        diag_arr[t, a0, t0, 4] = temp[4] * d[ag, 4, t];
        diag_arr[t, a0, t0, 5] = temp[4] * (1 - d[ag, 4, t]) * q_age[4, a0];
      }
    }
  }
  
  diag_arr1 = diag_arr_fct(diag_arr, empty_mat, nquar, nage);
  
  for (t in 1 : nquar) {
    for (a in 1 : nage) {
      exp_diag[t, a,  : ] = col_sums(diag_arr1[t, a,  : ,  : ]);
      // expected HIV_age diagnoses
      exp_HIV_dx[t, a] = exp_diag[t, a, 1] + exp_diag[t, a, 2]
                          + exp_diag[t, a, 3] + exp_diag[t, a, 4];
      // Expected AIDS_age diagnoses
      exp_AIDS_dx[t, a] = exp_diag[t, a, 5];
    }
  }
  
  /// For a = 1 and some times (t=5, t=9 ...)
  // there are 0 diagnes arrivals
  // This is due to model assumptions
  // So removing age 1 from likelihood
  // Fixing a=1 to get a valid simplex
  for (t in 1 : nquar) {
    exp_p[t, 1, 1] = 1;
    exp_p[t, 1, 2] = 0;
    exp_p[t, 1, 3] = 0;
    exp_p[t, 1, 4] = 0;
  }
  
  for (t in 1 : nquar) {
    for (a in 2 : nage) {
      // Now we want the proportion of CD4_age that
      // is diangosed with HIV_age in state k at time t and age a (p_kta)
      exp_p[t, a, 1] = exp_diag[t, a, 1] / exp_HIV_dx[t, a];
      exp_p[t, a, 2] = exp_diag[t, a, 2] / exp_HIV_dx[t, a];
      exp_p[t, a, 3] = exp_diag[t, a, 3] / exp_HIV_dx[t, a];
      exp_p[t, a, 4] = exp_diag[t, a, 4] / exp_HIV_dx[t, a];
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
    for (a in 2 : nage) {
      HIV_age[t, a] ~ poisson(exp_HIV_dx[t, a]);
      AIDS_age[t, a] ~ poisson(exp_AIDS_dx[t, a]);
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
