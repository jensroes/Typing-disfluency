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

//transformed data{
//  matrix[nS, maxB] logy;
//  
//  for(s in 1:nS){
//    int nBS = nB[s];
//    for(b in 1:nBS){ 
//      logy[s,b] = log(y[s,b]);
//    }
//  }
//}


parameters {
//  real<lower=0> beta;

	// Parameters for non-centering
	real beta_mu;
  real beta_sigma;	
	vector[2] beta_raw;			// distributions

	real sigma_mu;
  real sigma_sigma;	
	real sigma_raw;			// distributions

   // For random effects
	vector[nS] u; //subj intercepts
  real<lower=0> sigma_u;//subj sd
  
  vector[maxB] w; //bigram intercepts
  real<lower=0> sigma_w;//bigram sd

}

transformed parameters{
  vector[N] RE;
	real<lower=0> sigma;		// residual sd
  vector[2] beta;

  beta = beta_sigma * beta_raw + beta_mu;
  sigma = sigma_mu + sigma_sigma*sigma_raw;
  for(n in 1:N){
    RE[n] = u[subj[n]] + w[bigram[n]];
  }
}

model {
  int n = 0;

  // Priors
  beta_mu ~ normal(5, 2.5);
  beta_raw ~ normal(0, 1);
  beta_sigma ~ normal(0, 2);

  sigma_mu ~ cauchy(0, 2.5);
  sigma_sigma ~ cauchy(0, 2.5);
  sigma_raw ~ normal(0, 1);

	// REs priors
  sigma_u ~ normal(0,2.5);
  u ~ normal(0, sigma_u); //subj random effects

  sigma_w ~ normal(0,2.5);
  w ~ normal(0, sigma_w); //bigram random effects

  // likelihood
  for(s in 1:nS){
  int nBS = nB[s];
    for(b in 1:nBS){
      n += 1;
      if(b == 1){
        y[s, b] ~ lognormal(beta[1] + RE[n], sigma); 
      }
      if(b > 1){
        y[s, b] ~ lognormal(beta[2] + RE[n], sigma); 
      }
    }
  }
}

generated quantities{
	real log_lik[N];
	vector[N] y_tilde;
  int n = 0;


  for(s in 1:nS){
  int nBS = nB[s];
    for(b in 1:nBS){
      n += 1;
      if(b == 1){
        log_lik[n] = lognormal_lpdf(y[s, b] | beta[1] + RE[n], sigma); 
        y_tilde[n] = lognormal_rng(beta[1] + RE[n], sigma);
      }
      if(b > 1){
        log_lik[n] = lognormal_lpdf(y[s, b] | beta[2] + RE[n], sigma); 
        y_tilde[n] = lognormal_rng(beta[2] + RE[n], sigma);
      }
    }
  }
}
