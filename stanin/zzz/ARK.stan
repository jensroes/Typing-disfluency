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
  real phi[K];
  real<lower=0> sigma;

	// Parameters for non-centering
	real beta_mu;
  real<lower =0> beta_sigma;	
	real beta_raw;			// distributions

  // For random effects
  vector[nS] u; //subject intercepts
  real<lower=0> sigma_u;	// subj sd
}

transformed parameters{
  real beta = beta_sigma * beta_raw + beta_mu;
}

model {
  // Priors
  beta_mu ~ normal(5, 4);
  beta_sigma ~ normal(0, 10);
  beta_raw ~ normal(0, 1);
  sigma ~ cauchy(0, 2.5);
  phi ~ normal(0, 1);
  
  // REs priors
  sigma_u ~ normal(0,2.5);
  u ~ normal(0, sigma_u); //subj random effects
  
  // Likelihood	
  for(s in 1:nS){
    int nBS = nB[s];
    for(b in (K+1):nBS){
      real mu = beta + u[s];
      for (k in 1:K){
        mu += phi[k] * logy[s,b-k];
        y[s, b] ~ lognormal(mu, sigma);
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
        log_lik[n] = lognormal_lpdf(y[s, b] | mu, sigma); 
        y_tilde[n] = lognormal_rng(mu, sigma);
      }
    }
  }
}


