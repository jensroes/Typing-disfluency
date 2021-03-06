// LMM with REs for subject and bigrams 

data {
	int<lower=1> N;                    // Number of observations
  int<lower=1> nS;            //number of subjects
  int<lower=1> nB[nS];        // total number of bigrams produced by ppt
  int<lower=1> maxB;          // Max number of bigrams for matrix
  matrix[nS,maxB] y;            //outcome: for each subject one IKI per column/bigram (in 
}


parameters {
	real<lower=0> delta;
	real theta;

  real beta_mu;
  real beta_raw;
  real<lower=0> beta_sigma;

	real<lower=0> sigma;		// residual sd
  real<lower=0> sigma_diff;

   // For random effects
	vector[nS] u; //subj intercepts
  real<lower=0> sigma_u;//subj sd
  
  vector[maxB-1] w; //bigram intercepts
  real<lower=0> sigma_w;//bigram sd

}

transformed parameters{
  vector[2] log_theta;
  real prob = 1 - inv_logit(theta); 
  real beta = beta_mu + beta_sigma * beta_raw;
  real<lower=0> sigmap_e = sigma + sigma_diff;
  real<lower=0> sigma_e = sigma - sigma_diff;

  log_theta[1] = log_inv_logit(theta);
  log_theta[2] = log1m_inv_logit(theta);
}

model {
  vector[2] lp_parts;

  // Priors
  beta_mu ~ normal(5, 5);
  beta_sigma ~ cauchy(0, 2.5);
  beta_raw ~ normal(0, 1);
  sigma ~ cauchy(0, 2.5);
  delta ~ normal(0, 1);
  theta ~ normal(0, 1);

	// REs priors
  sigma_u ~ normal(0,2.5);
  u ~ normal(0, sigma_u); //subj random effects

  sigma_w ~ normal(0,2.5);
  w ~ normal(0, sigma_w); //bigram random effects

  // likelihood
  for(s in 1:nS){
    int nBS = nB[s];
    for(b in 2:nBS){
      real mu = beta + u[s] + w[b-1];
      lp_parts[1] = log_theta[1] + lognormal_lpdf(y[s,b] | mu, sigma_e); // typing
      lp_parts[2] = log_theta[2] + lognormal_lpdf(y[s,b] | mu + delta, sigmap_e); // looking at text
      target += log_sum_exp(lp_parts);
    }
  }
}

generated quantities{
	vector[N-nS] log_lik;
	vector[N-nS] y_tilde;
  vector[2] lp_parts;
  real<lower=0,upper=1> prob_tilde; 
  real beta2 = beta + delta;
  int n = 0;
 
  for(s in 1:nS){
    int nBS = nB[s];
    for(b in 2:nBS){
      real mu = beta + u[s] + w[b-1];
      n += 1;
      lp_parts[1] = log_theta[1] + lognormal_lpdf(y[s,b] | mu, sigma_e); // typing
      lp_parts[2] = log_theta[2] + lognormal_lpdf(y[s,b] | mu + delta, sigmap_e); // looking at text
      log_lik[n] = log_sum_exp(lp_parts);
   		prob_tilde = bernoulli_rng(prob); 
      if(prob_tilde) { 
        y_tilde[n] = lognormal_rng(mu + delta, sigmap_e);
      }
      else{
        y_tilde[n] = lognormal_rng(mu, sigma_e);
      }
    }
  }
}
