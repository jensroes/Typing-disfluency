/*
  Mixture model for disfluencies
  With by-ppt mixing proportion
  Random intercepts for ppts
  Random intercepts for bigrams
  By-ppts slope adjustments
*/

data {
  int<lower=1> N;                    // number of observations
  int<lower=1> nS;                   // number of ppts
  int<lower=1, upper=nS> ppts[N];    // ppts identifier
  int<lower=1> nB;                   // number of bigrams
  int<lower=1, upper=nB> bigrams[N]; // ppts identifier
  int<lower=1> K;                    // number of conditions
  int condition[N];                  // vector with condition levels
  vector[N] y;                       // outcome
}


parameters {
  real beta_mu;                   // fluent typing
  vector[K] beta_raw;
  real<lower=0> beta_sigma;

  vector<lower = 0>[K] delta;    // slowdown for hesitations

  vector[K] alpha;               // mixing proportion
  real tau;
  matrix[K,nS] alpha_s;          // by-ppt mixing proportion

  vector<lower = 0>[K] sigma;    // residual error
  vector<lower=0>[K] sigma_diff; // difference between mixture components

  // For random effects
  vector<lower=0>[K] sigma_u;   // ppts slopes
  cholesky_factor_corr[K] L_u;  // variance-covariance matrix
  matrix[K,nS] z_u;             // random effects for participants
  
  vector[nB] w;                 // bigram intercepts
  real<lower=0> sigma_w;        // bigram sd
}

transformed parameters{
  vector[K] theta = 1 - inv_logit(alpha);         // convert mixing proportion to proportion of second component
  matrix[K,nS] theta_s = 1 - inv_logit(alpha_s);
  matrix[K,nS] log_alpha_s_1 = log_inv_logit(alpha_s);
  matrix[K,nS] log_alpha_s_2 = log1m_inv_logit(alpha_s);
  
  vector[K] beta = beta_mu + beta_sigma * beta_raw;

  vector<lower=0>[K] sigmap_e = sigma + sigma_diff; // variance for second mixture component is larger
  vector<lower=0>[K] sigma_e = sigma - sigma_diff;  // variance for first mixture component is smaller      

  matrix[K,nS] u = diag_pre_multiply(sigma_u, L_u) * z_u; // ppts random effects

  vector[N] mu;
  
  for(n in 1:N){
    mu[n] = beta[condition[n]] + u[condition[n], ppts[n]] + w[bigrams[n]];
  }
}

model {
  vector[2] lp_parts;

  // Priors
  beta_mu ~ normal(5, 2.5);
  beta_sigma ~ cauchy(0, 2.5);
  beta_raw ~ normal(0, 1);

  delta ~ normal(0, 1);

  for(s in 1:nS){
    for(k in 1:K){
     alpha_s[k,s] ~ normal(alpha[k], tau);  
    }
  }

  alpha ~ normal(0, 1);
  tau ~ cauchy(0, 1);

  sigma ~ cauchy(0, 2.5);
  sigma_diff ~ normal(0, 1);

  // REs priors
  L_u ~ lkj_corr_cholesky(2.0);
  to_vector(z_u) ~ normal(0,1);
  sigma_u ~ normal(0,2.5);      // ppts random effects
  
  sigma_w ~ normal(0,2.5);
  w ~ normal(0, sigma_w);       // bigram random effects

  // likelihood
  for(n in 1:N){
    lp_parts[1] = log_alpha_s_1[condition[n], ppts[n]] + lognormal_lpdf(y[n] | mu[n], sigma_e[condition[n]]); // fluent typing
    lp_parts[2] = log_alpha_s_2[condition[n], ppts[n]] + lognormal_lpdf(y[n] | mu[n] + delta[condition[n]], sigmap_e[condition[n]]); // hesitatant typing
    target += log_sum_exp(lp_parts);
  }
}

generated quantities{
  vector[2] lp_parts;
  real<lower=0,upper=1> prob_tilde; 

  vector[N] log_lik; // log likelihood
  vector[N] y_tilde; // predicted data
  vector[K] beta2 = beta + delta;    // second mixture component
  matrix[K,nS] beta_s;               // by-ppt fluent typing estimates

  for(s in 1:nS){
    for(k in 1:K){
      beta_s[k, s] = beta[k] + u[k, s];
    }
  }

  for(n in 1:N){
    // calculate log likelihood
    lp_parts[1] = log_alpha_s_1[condition[n], ppts[n]] + lognormal_lpdf(y[n] | mu[n], sigma_e[condition[n]]); // fluent typing
    lp_parts[2] = log_alpha_s_2[condition[n], ppts[n]] + lognormal_lpdf(y[n] | mu[n] + delta[condition[n]], sigmap_e[condition[n]]); // hesitant typing
    log_lik[n] = log_sum_exp(lp_parts);
   	prob_tilde = bernoulli_rng(theta_s[condition[n], ppts[n]]); 
   	// predict data
    if(prob_tilde) { 
      y_tilde[n] = lognormal_rng(mu[n] + delta[condition[n]], sigmap_e[condition[n]]);
    }
    else{
      y_tilde[n] = lognormal_rng(mu[n], sigma_e[condition[n]]);
    }
  }
}
