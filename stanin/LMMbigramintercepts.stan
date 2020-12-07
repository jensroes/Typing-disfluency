// LMM with REs for subject and bigrams 

data {
	int<lower=1> N;                    // Number of observations
  int<lower=1> nS;            //number of subjects
  int<lower=1> nB[nS];        // total number of bigrams produced by ppt
  int<lower=1> maxB;          // Max number of bigrams for matrix
  matrix[nS,maxB] y;            //outcome: for each subject one IKI per column/bigram (in 
}


parameters {
	real<lower=0> sigma;		// residual sd

	// Parameters for non-centering
	real beta_mu;
  real<lower =0> beta_sigma;	
	real beta_raw;			// distributions

   // For random effects
	vector[nS] u; //subj intercepts
  real<lower=0> sigma_u;//subj sd
  
  vector[maxB] w; //bigram intercepts
  real<lower=0> sigma_w;//bigram sd

}

transformed parameters{
  real beta = beta_mu + beta_sigma * beta_raw;
}

model {
  // Priors
  beta_mu ~ normal(5, 4);
  beta_sigma ~ normal(0, 10);
  beta_raw ~ normal(0, 1);
  sigma ~ cauchy(0, 2.5);
	
	// REs priors
  sigma_u ~ normal(0,2.5);
  u ~ normal(0, sigma_u); //subj random effects

  sigma_w ~ normal(0,2.5);
  w ~ normal(0, sigma_w); //bigram random effects

  // likelihood
  for(s in 1:nS){
    int nBS = nB[s];
    for(b in 1:nBS){
      real mu = beta + u[s] + w[b];
      y[s, b] ~ lognormal(mu, sigma); 
    }
  }
}

generated quantities{
	vector[N] log_lik;
	vector[N] y_tilde;
  int n = 0;

  for(s in 1:nS){
    int nBS = nB[s];
    for(b in 1:nBS){
      real mu = beta + u[s] + w[b];
      n += 1;
      log_lik[n] = lognormal_lpdf(y[s, b] | mu, sigma); 
      y_tilde[n] = lognormal_rng(mu, sigma);
    }
  }
}
