/*
  Serial order of bigrams auto regression model (first order)
  with mixture for correct responses and errors
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
  // Parameters for non-centering
	real beta_mu;
  real<lower =0> beta_sigma;	
	real beta_raw;			// distributions
	
  real<lower=0> delta_corr; // beta + delta for mixture in correct responses
  real<lower=0> delta_incorr; // beta + delta for mixture in incorrect responses
  real<lower = .001, upper = .999> theta_corr; // mixing proportion for correct responses
  real<lower = .001, upper = .999> theta_incorr; // mixing proportion for incorrect responses
  real phi[K]; // for autoregression

  // Residual error
  real<lower=0> sigma;
	real<lower=0> sigma_diff;
	real<lower=0> sigma_diff_corr;
	real<lower=0> sigma_diff_incorr;
  
  // For random effects
  vector[nS] u; //subject intercepts
  real<lower=0> sigma_u;	// subj sd
}

transformed parameters{
  vector[2] log_theta_corr;
  vector[2] log_theta_incorr;
	real<lower=0> sigmap_e_corr;
	real<lower=0> sigmap_e_incorr;
	real<lower=0> sigma_e;
  real beta = beta_mu + beta_sigma * beta_raw;

  sigmap_e_incorr = sigma + sigma_diff_incorr;
  sigmap_e_corr = sigma + sigma_diff_corr;
  sigma_e = sigma - sigma_diff;

  log_theta_corr[1] = log(theta_corr);
  log_theta_corr[2] = log1m(theta_corr);
  log_theta_incorr[1] = log(theta_incorr);
  log_theta_incorr[2] = log1m(theta_incorr);

}



model {
  int n = 0;
  vector[2] lp_parts;

  // Priors
  beta_mu ~ cauchy(5, 2.5);
  beta_sigma ~ normal(0, 2);
  beta_raw ~ normal(0, 1);
  delta_corr ~ normal(0, 1);
  delta_incorr ~ normal(0, 1);
  theta_corr ~ beta(2, 2);
  theta_incorr ~ beta(2, 2);
  phi ~ normal(0, 1);
  
  sigma ~ cauchy(0, 1);
  sigma_diff_corr ~ normal(0, 1);
  sigma_diff_incorr ~ normal(0, 1);

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
        if(correct[s, b] == 1){// Correct responses are a mixture of normal typing and pausing (e.g. stimulus encoding)
          lp_parts[1] = log_theta_corr[1] + lognormal_lpdf(y[s,b] | mu, sigma_e); // typing
          lp_parts[2] = log_theta_corr[2] + lognormal_lpdf(y[s,b] | u[s] + beta + delta_corr, sigmap_e_corr); // looking at text
          target += log_sum_exp(lp_parts);
        }
        if(correct[s, b] == 0){// Incorrect responses are a mixture of unnoticed mistakes and other mistakes. 
          lp_parts[1] = log_theta_incorr[1] + lognormal_lpdf(y[s,b] | mu, sigma_e); // unnoticed errors
          lp_parts[2] = log_theta_incorr[2] + lognormal_lpdf(y[s,b] | u[s] + beta + delta_incorr, sigmap_e_incorr); // other errors
          target += log_sum_exp(lp_parts); 
        }
      }
    }
  }
}


generated quantities{
  vector[N-(nS*K)] log_lik;
  vector[N-(nS*K)] y_tilde;
  int n = 0;
  vector[2] lp_parts;
  real<lower=0,upper=1> theta_tilde; 

  
  for(s in 1:nS){
    int nBS = nB[s];
    for(b in (K+1):nBS){
      real mu = beta + u[s];
      for(k in 1:K){
        n += 1;
        mu += phi[k] * logy[s,b-k];
        if(correct[s, b] == 1){// Correct responses are a mixture of normal typing and pausing (e.g. stimulus encoding)
          lp_parts[1] = log_theta_corr[1] + lognormal_lpdf(y[s,b] | mu, sigma_e); 
          lp_parts[2] = log_theta_corr[2] + lognormal_lpdf(y[s,b] | u[s] + beta + delta_corr, sigmap_e_corr); 
          log_lik[n] = log_sum_exp(lp_parts); 
       		theta_tilde = bernoulli_rng(theta_corr); 
          if(theta_tilde) { 
            y_tilde[n] = lognormal_rng(mu, sigma_e);
          }
          else{
            y_tilde[n] = lognormal_rng(u[s] + beta + delta_corr, sigmap_e_corr);
          }
        }
        if(correct[s, b] == 0){
          lp_parts[1] = log_theta_incorr[1] + lognormal_lpdf(y[s,b] | mu, sigma_e); 
          lp_parts[2] = log_theta_incorr[2] + lognormal_lpdf(y[s,b] | u[s] + beta + delta_incorr, sigmap_e_incorr); 
          log_lik[n] = log_sum_exp(lp_parts); 
       		theta_tilde = bernoulli_rng(theta_incorr); 
          if(theta_tilde) { 
            y_tilde[n] = lognormal_rng(mu, sigma_e);
          }
          else{
            y_tilde[n] = lognormal_rng(u[s] + beta + delta_incorr, sigmap_e_incorr);
          }
        }
      }
    }
  }
}





