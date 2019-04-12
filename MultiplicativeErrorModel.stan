functions {
  real[] LVC1(real[] t, real r, real alpha, real x0){
    real k = r / alpha;
    real c = x0 / (k - x0);
    int N = size(t);
    vector[N] x = rep_vector(k, N) ./ (1 + exp(- r * to_vector(t)) / c);
    
    return(to_array_1d(x));
  }
  
  real[] LVC(real t,
             real[] p,
             real[] theta,
             real[] x_r,
             int[] x_i){
    int d = x_i[1];
    real r[d] = theta[1:d];
    matrix[d, d] alpha = to_matrix(to_vector(segment(theta, d + 1, d * d)), d, d); // Matrix filled in column-major order
    real dpdt[d] = to_array_1d(to_vector(p) .*  (to_vector(r) - alpha * to_vector(p)));
    return(dpdt);
  }
}

data {
  int<lower = 1> D; // Number of species (2 here)
  
  // Species 1, Experiment A
  int<lower = 0> N1_A; // Number of observations
  real<lower = 0> t1_A[N1_A]; // Measurement times
  real<lower = 0> y1_A[N1_A]; // Measured population
  
  // Species 1, Experiment B
  int<lower = 0> N1_B; // Number of observations
  real<lower = 0> t1_B[N1_B]; // Measurement times
  real<lower = 0> y1_B[N1_B]; // Measured population
  
  // Species 2, Experiment A
  int<lower = 0> N2_A; // Number of observations
  real<lower = 0> t2_A[N2_A]; // Measurement times
  real<lower = 0> y2_A[N2_A]; // Measured population
  
  // Species 2, Experiment B
  int<lower = 0> N2_B; // Number of observations
  real<lower = 0> t2_B[N2_B]; // Measurement times
  real<lower = 0> y2_B[N2_B]; // Measured population

  // Species 1 and 2, Experiment A
  int<lower = 0> N12_A; // Number of observations
  real<lower = 0> t12_A[N12_A]; // Measurement times
  real<lower = 0> y12_A[N12_A, D]; // Measured populations
  
  // Species 1 and 2, Experiment B
  int<lower = 0> N12_B; // Number of observations
  real<lower = 0> t12_B[N12_B]; // Measurement times
  real<lower = 0> y12_B[N12_B, D]; // Measured populations
  
  // Replications
  int<lower = 0> N_rep; // Number of timepoints for replications
  real<lower = 0> t_rep[N_rep]; // Measurement times for replications
}

transformed data {
  int N = N1_A + N1_B + N2_A + N2_B + 2 * (N12_A + N12_B);
}

parameters {
  real<lower = 0> r[D]; // growth parameters
  matrix<lower = 0>[D, D] alpha; // competition rates
  real<lower = 0> sigma; // measurement noise standard deviation
  real<lower = 0> sigma_alpha; // scale for alpha
  
  real<lower = 0> f0[8]; // initial populations
  real<lower = 0> sigma_f0; // population std for f0
}

transformed parameters {
  vector[2] k_uni = to_vector(r) ./ diagonal(alpha); // Steady states for single species
  vector[2] k_multi = inverse(alpha) * to_vector(r); // Steady states for multiple species
  
  real f1_A[N1_A] = LVC1(t1_A, r[1], alpha[1, 1], f0[1]);
  real f1_B[N1_B] = LVC1(t1_B, r[1], alpha[1, 1], f0[2]);
  real f2_A[N2_A] = LVC1(t2_A, r[2], alpha[2, 2], f0[3]);
  real f2_B[N2_B] = LVC1(t2_B, r[2], alpha[2, 2], f0[4]);
  
  real theta_12[D + D * D];
  real<lower = 0> f12_A[N12_A, 2];
  real<lower = 0> f12_B[N12_B, 2];
  
  theta_12[1:D] = r;
  theta_12[(D+1):(D + D * D)] = to_array_1d(to_vector(alpha));
  f12_A = integrate_ode_bdf(LVC, { f0[5], f0[6] }, 0, t12_A, theta_12, rep_array(0.0, 0), { D });
  f12_B = integrate_ode_bdf(LVC, { f0[7], f0[8] }, 0, t12_B, theta_12, rep_array(0.0, 0), { D });
}

model {
  r ~ cauchy(0, 1);
  sigma ~ cauchy(0, 0.1);
  sigma_alpha ~ cauchy(0, 0.5);
  
  f0 ~ normal(0, sigma_f0);
  sigma_f0 ~ cauchy(0, 0.1);
  
  y1_A ~ lognormal(log(f1_A), sigma);
  y1_B ~ lognormal(log(f1_B), sigma);
  y2_A ~ lognormal(log(f2_A), sigma);
  y2_B ~ lognormal(log(f2_B), sigma);
  
  for (i in 1:D){
    for (j in 1:D){
      alpha[i, j] ~ cauchy(0, sigma_alpha);
    }
    y12_A[, i] ~ lognormal(log(f12_A[, i]), sigma);
    y12_B[, i] ~ lognormal(log(f12_B[, i]), sigma);
  }
}

generated quantities {
  real log_lik[N];
  real y1_rep[N_rep + 1, 2]; 
  real y2_rep[N_rep + 1, 2];
  real y12_rep[N_rep + 1, D, 2];
  
  // Posterior predictive checks
  {
    real f1_rep[N_rep, 2];
    real f2_rep[N_rep, 2];
    real f12_rep[N_rep, D, 2];
    
    f1_rep[, 1] = LVC1(t_rep, r[1], alpha[1, 1], f0[1]);
    f1_rep[, 2] = LVC1(t_rep, r[1], alpha[1, 1], f0[2]);
    
    f2_rep[, 1] =  LVC1(t_rep, r[2], alpha[2, 2], f0[3]);
    f2_rep[, 2] =  LVC1(t_rep, r[2], alpha[2, 2], f0[4]);
    
    f12_rep[, , 1] = integrate_ode_bdf(LVC, { f0[5], f0[6] }, 0, t_rep, theta_12, rep_array(0.0, 0), { D });
    f12_rep[, , 2] = integrate_ode_bdf(LVC, { f0[7], f0[8] }, 0, t_rep, theta_12, rep_array(0.0, 0), { D });
  
    y1_rep[1, 1:2] = lognormal_rng(log(f0[1:2]), sigma);
    y2_rep[1, 1:2] = lognormal_rng(log(f0[3:4]), sigma);
    y12_rep[1, 1:2, 1] = lognormal_rng(log(f0[5:6]), sigma);
    y12_rep[1, 1:2, 2] = lognormal_rng(log(f0[7:8]), sigma);
  
    for (ex in 1:2){
      y1_rep[2:(N_rep + 1), ex] = lognormal_rng(log(f1_rep[, ex]), sigma);
      y2_rep[2:(N_rep + 1), ex] = lognormal_rng(log(f2_rep[, ex]), sigma);
      for (s in 1:2){
        y12_rep[2:(N_rep + 1), s, ex] = lognormal_rng(log(f12_rep[, s, ex]), sigma);
      }
    }
  }

  // Log Likelihood
  {
    int ct = 1;
    for (i in 1:N1_A){
      log_lik[ct] = lognormal_lpdf(y1_A[i] | log(f1_A[i]), sigma);
      ct = ct + 1;
    }
    for (i in 1:N1_B){
      log_lik[ct] = lognormal_lpdf(y1_B[i] | log(f1_B[i]), sigma);
      ct = ct + 1;
    }
    for (i in 1:N2_A){
      log_lik[ct] = lognormal_lpdf(y2_A[i] | log(f2_A[i]), sigma);
      ct = ct + 1;
    }
    for (i in 1:N2_B){
      log_lik[ct] = lognormal_lpdf(y2_B[i] | log(f2_B[i]), sigma);
      ct = ct + 1;
    }
    for (i in 1:N12_A){
      for (s in 1:2){
        log_lik[ct] = lognormal_lpdf(y12_A[i, s] | log(f12_A[i, s]), sigma);
        ct = ct + 1;
      }
    }
    for (i in 1:N12_B){
      for (s in 1:2){
        log_lik[ct] = lognormal_lpdf(y12_B[i, s] | log(f12_B[i, s]), sigma);
        ct = ct + 1;
      }
    }
  }
  
}