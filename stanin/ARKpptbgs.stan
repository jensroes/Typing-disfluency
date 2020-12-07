/*
  Serial order of bigrams auto regression model
  Random intercepts for subject 
*/
  
data {
  int<lower=0> N;
  int<lower=1> nS;            //number of subjects
  int<lower=1> nB[nS];        // total number of bigrams produced by ppt
  int<lower=1> maxB;          // Max number of bigrams for matrix
  matrix[nS,maxB] y;            //outcome: for each subject one IKI per column/bigram (in order as produced)
  int<lower=0> K; // order of autocorrelation
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
  vector[K] phi_s[nS];
  real phi[K];
  real<lower=0> tau; // between ppts variation in autoregressor
  
  real<lower=0> sigma;

	// Parameters for non-centering
	real beta_mu;
  real<lower =0> beta_sigma;	
	real beta_raw;			// distributions

  // For random effects
  vector[nS] u; //subject intercepts
  real<lower=0> sigma_u;	// subj sd
  
  vector[maxB] w; // bigram intercepts
  real<lower=0> sigma_w;	// bigram sd
  
}

transformed parameters{
  real beta = beta_sigma * beta_raw + beta_mu;
}

model {
  // Priors
  beta_mu ~ normal(5, 4);
  beta_sigma ~ normal(0, 10);
  beta_raw ~ normal(0, 1);
  sigma ~ cauchy(0, 2.5);
  
  for(s in 1:nS){
    phi_s[s] ~ normal(phi, tau);
  }
  
  phi ~ normal(0, 1);
  tau ~ cauchy(0, 1);
  
  // REs priors
  sigma_u ~ normal(0,2.5);
  u ~ normal(0, sigma_u); //subj random effects

  sigma_w ~ normal(0,2.5);
  w ~ normal(0, sigma_w); //bigram random effects

  // Likelihood	
  for(s in 1:nS){
    int nBS = nB[s];
    for(b in 1:K){
      real mu = beta + u[s] + w[b];
      if(b == 1){      // Model of the first keystrokes
        y[s, b] ~ lognormal(mu, sigma);
      }
      if(b > 1){       // Model of the first k keystrokes after the k=1 
       for(k in 1:(b-1)){
         mu += phi_s[s,k] * logy[s,b-k];
       }
       y[s, b] ~ lognormal(mu, sigma);
      }
    } // Model of the keystrokes with full k-degree autoregression
    for(b in (K+1):nBS){
      real mu = beta + u[s] + w[b];
      for(k in 1:K){
        mu += phi_s[s,k] * logy[s,b-k];
      }
      y[s, b] ~ lognormal(mu, sigma);
    }
  }
}

generated quantities{
//  vector[N-(nS*K)] log_lik;
//  vector[N-(nS*K)] y_tilde;
  vector[N] log_lik;
  vector[N] y_tilde;
  int n = 0;
  
  for(s in 1:nS){
    int nBS = nB[s];
    for(b in 1:K){
      real mu = beta + u[s] + w[b];
      n += 1;
      if(b > 1){       // Model of the first k keystrokes after the k=1 
        for(k in 1:(b-1)){
          mu += phi_s[s,k] * logy[s,b-k];
        }
      }
      log_lik[n] = lognormal_lpdf(y[s, b] | mu, sigma);
      y_tilde[n] = lognormal_rng(mu, sigma);
    } // Model of the keystrokes with full k-degree autoregression
    for(b in (K+1):nBS){
      real mu = beta + u[s] + w[b];
      n += 1;
      for(k in 1:K){
        mu += phi_s[s,k] * logy[s,b-k];
      }
      log_lik[n] = lognormal_lpdf(y[s, b] | mu, sigma); 
      y_tilde[n] = lognormal_rng(mu, sigma);
    }
  }
}


