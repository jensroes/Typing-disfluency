// LMM with REs for subject and bigrams 

data {
	int<lower=1> N;                    // Number of observations
  int<lower=1> nS;            //number of subjects
  int<lower=1> nB[nS];        // total number of bigrams produced by ppt
  int<lower=1> maxB;          // Max number of bigrams for matrix
  matrix[nS,maxB] y;            //outcome: for each subject one IKI per column/bigram (in 
  matrix[nS,maxB] correct; // bigram correct 1 or not 0 
}


parameters {
	real<lower=0> delta;
  real<lower=0> tau;

	real<lower=0> delta_corr;
	real<lower=0> delta_incorr;

  real theta;
  real eta;
  
  real theta_corr;         // Hyperparameter for mixing proportion
  real theta_incorr;

	real<lower=0> sigma;		// residual sd
  real<lower=0> sigma_diff;
  real<lower=0> sigma_diff_corr;
  real<lower=0> sigma_diff_incorr;

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
  vector[2] log_theta_corr;
  vector[2] log_theta_incorr;
	real<lower=0> sigmap_e_corr = sigma + sigma_diff_corr;
	real<lower=0> sigmap_e_incorr = sigma + sigma_diff_incorr;
  real<lower=0> sigma_e = sigma - sigma_diff;
  real prob = 1 - inv_logit(theta); 
  real prob_corr = 1 - inv_logit(theta_corr); //probability of extreme values
  real prob_incorr = 1 - inv_logit(theta_incorr); 
  real beta = beta_mu + beta_sigma * beta_raw;
  
  log_theta_corr[1] = log_inv_logit(theta_corr);
  log_theta_corr[2] = log1m_inv_logit(theta_corr);
  log_theta_incorr[1] = log_inv_logit(theta_incorr);
  log_theta_incorr[2] = log1m_inv_logit(theta_incorr);
}

model {
  vector[2] lp_parts;

  // Priors
  beta_mu ~ normal(5, 5);
  beta_sigma ~ cauchy(0, 2.5);
  beta_raw ~ normal(0, 1);
  sigma ~ cauchy(0, 2.5);
  sigma_diff_corr ~ normal(0, 1);
  sigma_diff_incorr ~ normal(0, 1);

  delta_corr ~ normal(delta, tau);
  delta_incorr ~ normal(delta, tau);

  theta_corr ~ normal(theta, eta);
  theta_incorr ~ normal(theta, eta);

  delta ~ normal(0, 1);
  tau ~ cauchy(0, 2.5);

	theta ~ normal(0, 1);
	eta ~ cauchy(0, 2.5);
	
	// REs priors
  sigma_u ~ normal(0,2.5);
  u ~ normal(0, sigma_u); //subj random effects

  sigma_w ~ normal(0,2.5);
  w ~ normal(0, sigma_w); //bigram random effects

  // likelihood
  for(s in 1:nS){
    int nBS = nB[s];
    for(b in 2:nBS){
      real mu = beta + u[s] + w[b-1];
      if(correct[s, b] == 1){
        lp_parts[1] = log_theta_corr[1] + lognormal_lpdf(y[s,b] | mu, sigma_e); // typing
        lp_parts[2] = log_theta_corr[2] + lognormal_lpdf(y[s,b] | mu + delta_corr, sigmap_e_corr); // looking at text
        target += log_sum_exp(lp_parts);
      }
      if(correct[s, b] == 0){
        lp_parts[1] = log_theta_incorr[1] + lognormal_lpdf(y[s,b] | mu, sigma_e); 
        lp_parts[2] = log_theta_incorr[2] + lognormal_lpdf(y[s,b] | mu + delta_incorr, sigmap_e_incorr); 
        target += log_sum_exp(lp_parts);
      }
    }
  }
}

generated quantities{
	vector[N-nS] log_lik;
	vector[N-nS] y_tilde;
  int n = 0;
  vector[2] lp_parts;
  real<lower=0,upper=1> prob_tilde; 
  real beta2 = beta + delta;
  real beta_corr = beta + delta_corr;
  real beta_incorr = beta + delta_incorr;
  real prob_diff = prob_incorr - prob_corr;


  for(s in 1:nS){
    int nBS = nB[s];
    for(b in 2:nBS){
      real mu = beta + u[s] + w[b-1];
      n += 1;
      if(correct[s, b] == 1){
        lp_parts[1] = log_theta_corr[1] + lognormal_lpdf(y[s,b] | mu, sigma_e); // typing
        lp_parts[2] = log_theta_corr[2] + lognormal_lpdf(y[s,b] | mu + delta_corr, sigmap_e_corr); // looking at text
        log_lik[n] = log_sum_exp(lp_parts);
     		prob_tilde = bernoulli_rng(prob_corr); 
        if(prob_tilde) { 
          y_tilde[n] = lognormal_rng(mu + delta_corr, sigmap_e_corr);
        }
        else{
          y_tilde[n] = lognormal_rng(mu, sigma_e);
        }
      }
      if(correct[s, b] == 0){
        lp_parts[1] = log_theta_incorr[1] + lognormal_lpdf(y[s,b] | mu, sigma_e); // typing
        lp_parts[2] = log_theta_incorr[2] + lognormal_lpdf(y[s,b] | mu + delta_incorr, sigmap_e_incorr); // looking at text
        log_lik[n] = log_sum_exp(lp_parts);
     		prob_tilde = bernoulli_rng(prob_incorr); 
        if(prob_tilde) { 
          y_tilde[n] = lognormal_rng(mu + delta_incorr, sigmap_e_incorr);
        }
        else{
          y_tilde[n] = lognormal_rng(mu, sigma_e);
        }
      }
    }
  }
}
