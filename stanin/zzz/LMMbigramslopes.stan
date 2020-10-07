// LMM with REs for subject and by-subjects slopes for bigram order

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
	real beta_raw;		

  // For random effects
  vector<lower=0>[maxB-1] sigma_u;	// subj sd
  cholesky_factor_corr[maxB-1] L_u;
  matrix[maxB-1,nS] z_u; // tmp for subject intercepts and slopes

  vector[maxB-1] w; //bigram intercepts
  real<lower=0> sigma_w;//bigram sd
}

transformed parameters{
  real beta = beta_sigma * beta_raw + beta_mu;
  matrix[nS, maxB-1] u;
  u = (diag_pre_multiply(sigma_u, L_u) * z_u)';	// subj random effects
}

model {
   // Priors
  beta_mu ~ normal(5, 4);
  beta_sigma ~ normal(0, 10);
  beta_raw ~ normal(0, 1);

	sigma ~ cauchy(0, 2.5);

	// REs priors
	sigma_u ~ normal(0, 2.5);
	L_u ~ lkj_corr_cholesky(2.0);
	to_vector(z_u) ~ normal(0,1);	

  sigma_w ~ normal(0,2.5);
  w ~ normal(0, sigma_w); //bigram random effects

  // likelihood
  for(s in 1:nS){
    int nBS = nB[s];
    for(b in 2:nBS){
      real mu = beta + u[s,b-1] + w[b-1];
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
      real mu = beta + u[s,b-1] + w[b-1];
      n += 1;
      log_lik[n] = lognormal_lpdf(y[s, b] | mu, sigma); 
      y_tilde[n] = lognormal_rng(mu, sigma);
    }
  }
}
