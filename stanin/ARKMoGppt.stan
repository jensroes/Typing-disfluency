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
  vector[nS] theta_s; // by-ppt mixing proportion
  real<lower=0> tau_theta;

  real phi[K]; // for autoregression
  vector[K] phi_s[nS];
  real<lower=0> tau_phi;

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
  vector[nS] log_theta_s[2];
  real prob = 1 - inv_logit(theta); 
  vector[nS] prob_s = 1 - inv_logit(theta_s);
  real beta = beta_mu + beta_sigma * beta_raw;
  real<lower=0> sigmap_e = sigma + sigma_diff;
  real<lower=0> sigma_e = sigma - sigma_diff;

  log_theta_s[1] = log_inv_logit(theta_s);
  log_theta_s[2] = log1m_inv_logit(theta_s);
}



model {
  vector[2] lp_parts;
  
  // Priors
  beta_mu ~ normal(5, 1);
  beta_sigma ~ normal(0, 1);
  beta_raw ~ normal(0, 1);
  delta ~ normal(0, 1);
  sigma ~ cauchy(0, 2.5);
  sigma_diff ~ normal(0, 1);


  for(s in 1:nS){
    phi_s[s] ~ normal(phi, tau_phi);
  }

  phi ~ normal(0, 1);
  tau_phi ~ normal(0, 3);
  
  theta_s ~ normal(theta, tau_theta);
  theta ~ normal(0, 1);
  tau_theta ~ normal(0, 2.5);

  
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
        mu += phi_s[s,k] * logy[s,b-k];
        lp_parts[1] = log_theta_s[1,s] + lognormal_lpdf(y[s,b] | mu, sigma_e); // fluent
        lp_parts[2] = log_theta_s[2,s] + lognormal_lpdf(y[s,b] | beta + delta + u[s] + w[b-1], sigmap_e); // disfluent
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

  for(s in 1:nS){
    int nBS = nB[s];
    for(b in (K+1):nBS){
      real mu = beta + u[s];
      for(k in 1:K){
        n += 1;
        mu += phi_s[s,k] * logy[s,b-k];
        lp_parts[1] = log_theta_s[1,s] + lognormal_lpdf(y[s,b] | mu, sigma_e); 
        lp_parts[2] = log_theta_s[2,s] + lognormal_lpdf(y[s,b] | beta + delta + u[s] + w[b-1], sigmap_e); 
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
