/*
  Serial order of bigrams auto regression model (first order)
  Different typing speed for correct and incorrect responses
  Random intercepts for subject 
*/
  
data {
    int<lower=0> N;
    int<lower=1> nS;            //number of subjects
    int<lower=1> nB[nS];        // total number of bigrams produced by ppt
    int<lower=1> maxB;          // Max number of bigrams for matrix
    matrix[nS,maxB] y;            //outcome: for each subject one IKI per column/bigram (in order as produced)
    matrix[nS,maxB] correct; // bigram correct 1 or not 0 
    int K;
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
  real delta;
  real phi[K]; // for autoregression

  // Residual error
  real<lower=0> sigma;
  real<lower=0> sigma_diff;

	// Parameters for non-centering
	real beta_mu;
  real<lower =0> beta_sigma;	
	real beta_raw;			// distributions

  // For random effects
  vector[nS] u; //subject intercepts
  real<lower=0> sigma_u;	// subj sd
}

transformed parameters{
  real<lower=0> sigmap_e;
  real<lower=0> sigma_e;
  real beta = beta_sigma * beta_raw + beta_mu;
  
  sigmap_e = sigma + sigma_diff; // larger variance for incorrect responses
  sigma_e = sigma - sigma_diff;
}

model {
  int n = 0;
  vector[2] lp_parts;

  // Priors
  beta_mu ~ cauchy(5, 2.5);
  beta_sigma ~ normal(0, 2);
  beta_raw ~ normal(0, 1);
  phi ~ normal(0, 1);
  delta ~ normal(0, 1);

  // Residual error  
  sigma ~ cauchy(0, 1);
  sigma_diff ~ normal(0, 2);

  // REs priors
  sigma_u ~ normal(0,2.5);
  u ~ normal(0, sigma_u); //subj random effects
  
  // Likelihood	
  for(s in 1:nS){
    int nBS = nB[s];
    for(b in (K+1):nBS){
      real mu = beta + u[s];
      for(k in 1:K){
        mu += phi[k] * logy[s,b-k];
        if(correct[s, b] == 1){// Correct responses 
          y[s, b] ~ lognormal(mu, sigma_e);
        }
        if(correct[s, b] == 0){// Incorrect responses 
          y[s, b] ~ lognormal(mu + delta, sigmap_e);
        }
      }
    }
  }
}


generated quantities{
  vector[N-(nS*K)] log_lik;
  vector[N-(nS*K)] y_tilde;
  int n = 0;

  for(s in 1:nS){
    int nBS = nB[s];
    for(b in (K+1):nBS){
      real mu = beta + u[s];
      for (k in 1:K){
        n += 1;
        mu += phi[k] * logy[s,b-k];
        if(correct[s, b] == 1){// Correct responses 
          log_lik[n] = lognormal_lpdf(y[s, b] | mu, sigma_e); 
          y_tilde[n] = lognormal_rng(mu, sigma_e);
        }
        if(correct[s, b] == 0){// Correct responses 
          log_lik[n] = lognormal_lpdf(y[s, b] | mu + delta, sigmap_e); 
          y_tilde[n] = lognormal_rng(mu + delta, sigmap_e);
        }
      }
    }
  }
}
