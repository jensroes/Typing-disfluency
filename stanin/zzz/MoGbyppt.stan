// LMM with REs for subject and bigrams and by-subjects slopes

data {
  int<lower=1> N;                    // Number of observations
  int<lower=1> nS;            //number of subjects
  int<lower=1> nB[nS];        // total number of bigrams produced by ppt
  int<lower=1> maxB;          // Max number of bigrams for matrix
  matrix[nS,maxB] y;            //outcome: for each subject one IKI per column/bigram (in 
  int<lower=1, upper=nS> subj[N];   //subject id
  int<lower=1, upper=maxB> bigram[N];   //bigramid
}


parameters {
  //  real<lower=0> beta;
  real<lower = 0.001, upper = .999> theta; // mixing proportion

  // Parameters for non-centering
  real beta;


  real<lower=0> sigma;		// residual sd
	real<lower=0> sigma_diff;

  real<lower = 0> delta;
  
	vector<lower=0>[2] sigma_u;	// subj sd
	cholesky_factor_corr[2] L_u;
	matrix[2,nS] z_u; // tmp for subject intercepts and slopes
  
  vector[maxB] w; //bigram intercepts
  real<lower=0> sigma_w;//bigram sd
  
}

transformed parameters{
	matrix[nS,2] u;
	real sigmap_e;
	real sigma_e;
  vector[2] log_theta;

  sigmap_e = sigma + sigma_diff;
  sigma_e = sigma - sigma_diff;

  
	u = (diag_pre_multiply(sigma_u, L_u) * z_u)';	// subj random effects
  log_theta[1] = log(theta);
  log_theta[2] = log1m(theta);

}

model {
  int n = 0;
  vector[2] lp_parts;
    
  // Priors
  beta ~ normal(5, 2.5);

  sigma ~ cauchy(0, 1);
  sigma_diff ~ normal(0, 1);
  
  theta ~ beta(2, 2);
  delta ~ normal(0, 1);
  
  // REs priors
	sigma_u ~ normal(0, 2.5);
	L_u ~ lkj_corr_cholesky(2);
	to_vector(z_u) ~ normal(0,1);	
  
  sigma_w ~ normal(0,2.5);
  w ~ normal(0, sigma_w); //bigram random effects
  
  // likelihood
  for(s in 1:nS){
    int nBS = nB[s];
    for(b in 1:nBS){
      n += 1;
      lp_parts[1] = log_theta[1] + lognormal_lpdf(y[s,b] | beta + delta + u[subj[n],1] + w[bigram[n]], sigmap_e);
      lp_parts[2] = log_theta[2] + lognormal_lpdf(y[s,b] | beta + u[subj[n],2] + w[bigram[n]], sigma_e);
      target += log_sum_exp(lp_parts); 
    }
  }
}

generated quantities{
  real log_lik[N];
  vector[N] y_tilde;
  real<lower=0,upper=1> theta_tilde; 
  int n = 0;
  
  for(s in 1:nS){
    int nBS = nB[s];
    for(b in 1:nBS){
      n += 1;
      log_lik[n] = log_sum_exp(
 		      log_theta[1] + lognormal_lpdf(y[s,b] | beta + delta + u[subj[n],1] + w[bigram[n]], sigmap_e), 
   		    log_theta[2] + lognormal_lpdf(y[s,b] | beta + u[subj[n],2] + w[bigram[n]], sigma_e));
   		theta_tilde = bernoulli_rng(theta); 
      if(theta_tilde) { 
        y_tilde[n] = lognormal_rng(beta + delta + u[subj[n],1] + w[bigram[n]], sigmap_e);
      }
      else{
        y_tilde[n] = lognormal_rng(beta + u[subj[n],2] + w[bigram[n]], sigma_e);
      }
    }
  }
}
