/*
  Serial order of bigrams autoregression model
  Mixture model of first order and second order autoregression model
  Random intercepts for subject 
*/
  
data {
  int<lower=0> N;             // total number of observations
  int<lower=1> nS;              //number of subjects
  int nB[nS];                 // total number of bigrams produced by ppt
  int<lower=1> maxB;          // Max number of bigrams for matrix
  matrix[nS,maxB] y;            //outcome: for each subject one IKI per column (in order as produced)
	int<lower=1, upper=nS> subj[N];   //subject id
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
  real phi;
  real gamma;

  vector<lower = 0>[2] delta; 
  real<lower = 0> mu_delta;
  real<lower = 0> sigma_delta;
  real<lower = 0.001, upper = .999> theta[2]; // mixing proportion
	
  
  // Parameters for non-centering
	real<lower = 0> alpha_mu;
  real<lower = 0> alpha_sigma;	
	vector[2] alpha_raw;			// distributions
  
  real<lower=0> sigma;		// residual sd
	real<lower=0> sigma_diff;

  // For random effects
  vector[nS] u; //subject intercepts
  real<lower=0> sigma_u;	// subj sd
  
}


transformed parameters{
	vector[2] alpha = alpha_mu + alpha_sigma * alpha_raw;
	matrix[2,2] log_theta;
  vector[N] RE; // Random effects
	real sigmap_e;
	real sigma_e;

  log_theta[1,1] = log(theta[1]);
  log_theta[2,1] = log1m(theta[1]);
  log_theta[1,2] = log(theta[2]);
  log_theta[2,2] = log1m(theta[2]);

  sigmap_e = sigma + sigma_diff;
  sigma_e = sigma - sigma_diff;
  
  for(n in 1:N){
    RE[n] = u[subj[n]];
  }
}

model {
  int n = 0;
  vector[2] lp_parts;

  // Priors
  alpha_mu ~ cauchy(6, 2.5);
  alpha_sigma ~ cauchy(0, 2.5);
  alpha_raw ~ normal(0, 2);

  phi ~ normal(0,1);
  gamma ~ normal(0, 1);
  
  mu_delta ~ cauchy(0, 1.5);
  sigma_delta ~ cauchy(0, 2.5);
  delta ~ cauchy(mu_delta, sigma_delta);
  
  theta ~ beta(2, 2);
  sigma ~ cauchy(0, 2.5);
  sigma_diff ~ normal(0, 1);

  // REs priors
  sigma_u ~ normal(0,2.5);
  u ~ normal(0, sigma_u); //subj random effects
  

  // Likelihood	
  for(s in 1:nS){
    int nBS = nB[s];
    for(b in 1:nBS){ 
      n += 1;
      if(b == 1){
        lp_parts[1] = log_theta[1,1] + normal_lpdf(logy[s,b] | alpha[1] + delta[1] + RE[n], sigmap_e);
        lp_parts[2] = log_theta[2,1] + normal_lpdf(logy[s,b] | alpha[1] + RE[n], sigma_e);
      }
      if(b == 2){
        lp_parts[1] = log_theta[1,2] + normal_lpdf(logy[s,b] | alpha[2] + delta[2] + phi*logy[s, b-1] + RE[n], sigmap_e);
        lp_parts[2] = log_theta[2,2] + normal_lpdf(logy[s,b] | alpha[2] + phi*logy[s, b-1] + RE[n], sigma_e);
      }
      if(b > 2){
        lp_parts[1] = log_theta[1,2] + normal_lpdf(logy[s,b] | alpha[2] + delta[2] + phi*logy[s, b-1] + gamma*logy[s, b-2] + RE[n], sigmap_e);
        lp_parts[2] = log_theta[2,2] + normal_lpdf(logy[s,b] | alpha[2] + phi*logy[s, b-1] + gamma*logy[s, b-2] + RE[n], sigma_e);
      }
      target += log_sum_exp(lp_parts); 
    }
  }
}

generated quantities{
  vector[N] log_lik;
  vector[N] y_tilde;
  real<lower=0,upper=1> theta_tilde; 
  int n = 0;
  
  // likelihood: 
  for(s in 1:nS){
    int nBS = nB[s];
    for(b in 1:nBS){
      n += 1;
      if(b == 1){
        log_lik[n] = log_sum_exp(
 		      log_theta[1,1] + normal_lpdf(logy[s,b] | alpha[1] + delta[1] + RE[n], sigmap_e), 
   		    log_theta[2,1] + normal_lpdf(logy[s,b] | alpha[1] + RE[n], sigma_e));
   		  theta_tilde = bernoulli_rng(theta[1]); 
        if(theta_tilde) { 
           y_tilde[n] = normal_rng(alpha[1] + delta[1] + RE[n], sigmap_e);
        }
        else{
            y_tilde[n] = normal_rng(alpha[1] + RE[n], sigma_e);
        }
      }
      if(b == 2){      
	      log_lik[n] = log_sum_exp(
 		      log_theta[1,2] + normal_lpdf(logy[s,b] | alpha[2] + delta[2] + phi*logy[s, b-1] + RE[n], sigmap_e), 
   		    log_theta[2,2] + normal_lpdf(logy[s,b] | alpha[2] + phi*logy[s, b-1] + RE[n], sigma_e));
   		  theta_tilde = bernoulli_rng(theta[2]); 
        if(theta_tilde) { 
          y_tilde[n] = normal_rng(alpha[2] + delta[2] + phi*logy[s, b-1] + RE[n], sigmap_e);
        }
        else{
          y_tilde[n] = normal_rng(alpha[2] + phi*logy[s, b-1] + RE[n], sigma_e);
        }
      }
      if(b > 2){
        log_lik[n] = log_sum_exp(
   		    log_theta[1,2] + normal_lpdf(logy[s,b] | alpha[2] + delta[2] + phi*logy[s, b-1] + gamma*logy[s, b-2] + RE[n], sigmap_e), 
     		  log_theta[2,2] + normal_lpdf(logy[s,b] | alpha[2] + phi*logy[s, b-1] + gamma*logy[s, b-2] + RE[n], sigma_e));
        theta_tilde = bernoulli_rng(theta[2]); 
        if(theta_tilde) { 
          y_tilde[n] = normal_rng(alpha[2] + delta[2] + phi*logy[s, b-1] + gamma*logy[s, b-2] + RE[n], sigmap_e);
        }
        else{
          y_tilde[n] = normal_rng(alpha[2] + phi*logy[s, b-1] + gamma*logy[s, b-2] + RE[n], sigma_e);
        }
      }
    }
  }
}

