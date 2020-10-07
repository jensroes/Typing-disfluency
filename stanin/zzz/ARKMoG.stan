/*
  Serial order of bigrams auto regression model
  Random intercepts for subject 
*/
  
data {
  int<lower=0> N;
  int<lower=1> nS;            //number of subjects
  int<lower=1> nB[nS];        // total number of bigrams produced by ppt
  int<lower=1> maxB;          // Max number of bigrams for matrix
  matrix[nS,maxB] y;            //outcome: for each subject one IKI per column/bigram (in order as produced)
  int<lower=0> K; // order of autocorrelation
}

transformed data{
  matrix[nS, maxB] logy;
  
  for(s in 1:nS){
    int nBS = nB[s];
    for(b in 1:nBS){ 
      logy[s,b] = log(y[s,b]);
    }
  }
}


parameters {
  real<lower=0> delta; // beta + delta for mixture
  real theta; // mixing proportion
  real phi[K]; // for autoregression
  
	// Parameters for non-centering
	real beta_mu;
  real<lower =0> beta_sigma;	
	real beta_raw;			// distributions

  // Residual error
  real<lower=0> sigma;
  real<lower=0> sigma_diff;

  // For random effects
  vector[nS] u; //subject intercepts
  real<lower=0> sigma_u;	// subj sd

  vector[maxB-1] w; //bigram intercepts
  real<lower=0> sigma_w;//bigram sd

}

transformed parameters{
  vector[2] log_theta;
  real<lower=0> sigmap_e = sigma + sigma_diff;
  real<lower=0> sigma_e = sigma - sigma_diff;
  real beta = beta_mu + beta_sigma * beta_raw;
  real prob = 1 - inv_logit(theta);

  log_theta[1] = log_inv_logit(theta);
  log_theta[2] = log1m_inv_logit(theta);
}



model {
  vector[2] lp_parts;
  
  // Priors
  beta_mu ~ normal(5, 5);
  beta_sigma ~ cauchy(0, 2.5);
  beta_raw ~ normal(0, 1);
  delta ~ normal(0, 1);
  theta ~ normal(0, 2);
  phi ~ normal(0, 1);
  
  sigma ~ cauchy(0, 2.5);
  sigma_diff ~ normal(0, 1);

  // REs priors
  sigma_u ~ normal(0,2.5);
  u ~ normal(0, sigma_u); //subj random effects
  
  sigma_w ~ normal(0,2.5);
  w ~ normal(0, sigma_w); //bigram random effects
  
  // Likelihood	
  for(s in 1:nS){
    int nBS = nB[s];
    for(b in (K+1):nBS){
      real mu = u[s] + beta;
      for (k in 1:K){
        mu += phi[k] * logy[s,b-k];
        lp_parts[1] = log_theta[1] + lognormal_lpdf(y[s,b] | mu, sigma_e); // typing
        lp_parts[2] = log_theta[2] + lognormal_lpdf(y[s,b] | beta + delta + u[s] + w[b-1], sigmap_e); // looking at text
        target += log_sum_exp(lp_parts);
      }
    }
  }
}

generated quantities{
  vector[N-(nS*K)] log_lik;
  vector[N-(nS*K)] y_tilde;
  vector[2] lp_parts;
  real<lower=0,upper=1> theta_tilde; 
  int n = 0;
  real beta2 = beta + delta;
  vector[nS] beta_s = beta + u;
  vector[nS] beta2_s = beta + delta + u;

  for(s in 1:nS){
    int nBS = nB[s];
    for(b in (K+1):nBS){
      real mu = beta + u[s];
      for(k in 1:K){
        n += 1;
        mu += phi[k] * logy[s,b-k];
        lp_parts[1] = log_theta[1] + lognormal_lpdf(y[s,b] | mu, sigma_e); 
        lp_parts[2] = log_theta[2] + lognormal_lpdf(y[s,b] | beta + delta + u[s] + w[b-1], sigmap_e); 
        log_lik[n] = log_sum_exp(lp_parts); 
     		theta_tilde = bernoulli_rng(prob); 
        if(theta_tilde) { 
          y_tilde[n] = lognormal_rng(beta + delta + u[s] + w[b-1], sigmap_e);
        }
        else{
          y_tilde[n] = lognormal_rng(mu, sigma_e);
        }
      }
    }
  }
}
