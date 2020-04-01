// LMM with REs for subject and bigrams 

data {
  int<lower=1> N;                    // Number of observations
  int<lower=1> nS;            //number of subjects
  int<lower=1> nB[nS];        // total number of bigrams produced by ppt
  int<lower=1> maxB;          // Max number of bigrams for matrix
  matrix[nS,maxB] y;            //outcome: for each subject one IKI per column/bigram
  matrix[nS,maxB] correct; // bigram correct 1 or not 0 
}


parameters {
  real<lower=0> sigma;		// residual sd
  real delta;
  
  // Parameters for non-centering
  real beta_mu;
  real<lower =0> beta_sigma;	
  real beta_raw;			// distributions
  
  // For random effects
  vector[nS] u; //subj intercepts
  real<lower=0> sigma_u;//subj sd
  
  vector[maxB-1] w; //bigram intercepts
  real<lower=0> sigma_w;//bigram sd
  
}

transformed parameters{
  real beta = beta_sigma * beta_raw + beta_mu;
}

model {
  // Priors
  beta_mu ~ cauchy(5, 2.5);
  beta_sigma ~ normal(0, 2);
  beta_raw ~ normal(0, 1);
  delta ~ normal(0, 5);
  sigma ~ cauchy(0, 2.5);
  
  // REs priors
  sigma_u ~ normal(0,2.5);
  u ~ normal(0, sigma_u); //subj random effects
  
  sigma_w ~ normal(0,2.5);
  w ~ normal(0, sigma_w); //bigram random effects
  
  // likelihood
  for(s in 1:nS){
    int nBS = nB[s];
    for(b in 2:nBS){
      real mu = beta + delta*correct[s, b] + u[s] + w[b-1];
      y[s, b] ~ lognormal(mu, sigma); 
    }
  }
}

generated quantities{
  vector[N-nS] log_lik;
  vector[N-nS] y_tilde;
  int n = 0;
  
  for(s in 1:nS){
    int nBS = nB[s];
    for(b in 2:nBS){
      real mu = beta + delta*correct[s, b] + u[s] + w[b-1];
      n += 1;
      log_lik[n] = lognormal_lpdf(y[s, b] | mu, sigma); 
      y_tilde[n] = lognormal_rng(mu, sigma);
    }
  }
}
