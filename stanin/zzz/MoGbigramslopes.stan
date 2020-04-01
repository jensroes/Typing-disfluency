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
	real<lower=.001, upper=.999> theta;

	real<lower=0> sigma;		// residual sd
  real<lower=0> sigma_diff;

	// Parameters for non-centering
	real beta_mu;
  real<lower =0> beta_sigma;	
	real beta_raw;			// distributions

   // For random effects
  vector<lower=0>[maxB-1] sigma_u;	// subj sd
  cholesky_factor_corr[maxB-1] L_u;
  matrix[maxB-1,nS] z_u; // tmp for subject intercepts and slopes

  
  vector[maxB-1] w; //bigram intercepts
  real<lower=0> sigma_w;//bigram sd

}

transformed parameters{
  vector[2] log_theta;
  real<lower=0> sigmap_e;
  real<lower=0> sigma_e;
  matrix[nS, maxB-1] u;
  real beta = beta_sigma * beta_raw + beta_mu;

  u = (diag_pre_multiply(sigma_u, L_u) * z_u)';	// subj random effects
  
  sigmap_e = sigma + sigma_diff;
  sigma_e = sigma - sigma_diff;
  
  log_theta[1] = log(theta);
  log_theta[2] = log1m(theta);
}

model {
  vector[2] lp_parts;

  // Priors
  beta_mu ~ cauchy(5, 2.5);
  beta_sigma ~ normal(0, 2);
  beta_raw ~ normal(0, 1);
  sigma ~ cauchy(0, 2.5);
	
  delta ~ normal(0, 1);
  theta ~ beta(2, 2);

	// REs priors
	sigma_u ~ normal(0, 2.5);
	L_u ~ lkj_corr_cholesky(3);
	to_vector(z_u) ~ normal(0,1);	


  sigma_w ~ normal(0,2.5);
  w ~ normal(0, sigma_w); //bigram random effects

  // likelihood
  for(s in 1:nS){
    int nBS = nB[s];
    for(b in 2:nBS){
      real mu = beta + u[s,b-1] + w[b-1];
      lp_parts[1] = log_theta[1] + lognormal_lpdf(y[s,b] | mu, sigma_e); // typing
      lp_parts[2] = log_theta[2] + lognormal_lpdf(y[s,b] | mu + delta, sigmap_e); // looking at text
      target += log_sum_exp(lp_parts);
    }
  }
}

generated quantities{
	vector[N-nS] log_lik;
	vector[N-nS] y_tilde;
  int n = 0;
  vector[2] lp_parts;
  real<lower=0,upper=1> theta_tilde; 

  for(s in 1:nS){
    int nBS = nB[s];
    for(b in 2:nBS){
      real mu = beta + u[s,b-1] + w[b-1];
      n += 1;
      lp_parts[1] = log_theta[1] + lognormal_lpdf(y[s,b] | mu, sigma_e); // typing
      lp_parts[2] = log_theta[2] + lognormal_lpdf(y[s,b] | mu + delta, sigmap_e); // looking at text
      log_lik[n] = log_sum_exp(lp_parts);
   		theta_tilde = bernoulli_rng(theta); 
      if(theta_tilde) { 
        y_tilde[n] = lognormal_rng(mu, sigma_e);
      }
      else{
        y_tilde[n] = lognormal_rng(mu + delta, sigmap_e);
      }
    }
  }
}
