/*
  Serial order of bigrams auto regression model (first order)
Random intercepts for subject 
*/
  
  data {
    int<lower=0> N;
    int<lower=1> nS;            //number of subjects
    int<lower=1> nB[nS];        // total number of bigrams produced by ppt
    int<lower=1> maxB;          // Max number of bigrams for matrix
    matrix[nS,maxB] y;            //outcome: for each subject one IKI per column/bigram (in order as produced)
    int<lower=1, upper=nS> subj[N];   //subject id
    //  int<lower=0> nBre;
    //	int<lower=1, upper=nBre> bigram[N];   //bigram id
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
  
  // Parameters for non-centering
  real beta_mu;
  real<lower = 0> beta_sigma;	
  real beta_raw;			// distributions
  
  real sigma_mu;
  real sigma_sigma;	
  real sigma_raw;			// distributions
  
  
  // For random effects
  vector[nS] u; //subject intercepts
  real<lower=0> sigma_u;	// subj sd
  //  vector[nBre] w; //bigram intercepts
  //  real<lower=0> sigma_w;	// bigram sd
}

transformed parameters{
  vector[N] RE; // Random effects
  real<lower=0> sigma;		// residual sd
  real beta;
  
  sigma = sigma_mu + sigma_sigma*sigma_raw;
  beta = beta_sigma * beta_raw + beta_mu;
  
  for(n in 1:N){
    RE[n] = u[subj[n]]; // + w[bigram[n]];
  } 
}



model {
  int n = 0;
  
  // Priors
  beta_mu ~ cauchy(0, 2.5);
  beta_sigma ~ cauchy(0, 2.5);
  beta_raw ~ normal(0, 2);
  
  sigma_mu ~ cauchy(0, 2.5);
  sigma_sigma ~ cauchy(0, 2.5);
  sigma_raw ~ normal(0, 2);
  
  phi ~ normal(0, 1);
  
  // REs priors
  sigma_u ~ normal(0,2.5);
  u ~ normal(0, sigma_u); //subj random effects
  
  // Likelihood	
  for(s in 1:nS){
    int nBS = nB[s];
    for(b in 1:nBS){
      n += 1;
      if(b == 1){
        y[s, b] ~ lognormal(beta + RE[n], sigma); 
      }
      if(b > 1){
        y[s, b] ~ lognormal(beta + RE[n] + phi * logy[s, b-1], sigma); 
      }
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
      n += 1;
      if(b == 1){
        log_lik[n] = lognormal_lpdf(y[s, b] | beta + RE[n], sigma); 
        y_tilde[n] = lognormal_rng(beta + RE[n], sigma);
      }
      if(b > 1){
        log_lik[n] = lognormal_lpdf(y[s, b] | beta + RE[n] + phi * logy[s, b-1], sigma); 
        y_tilde[n] = lognormal_rng(beta + RE[n] + phi * logy[s, b-1], sigma);
      }
    }
  }
}


