/*
  Serial order of bigrams auto regression model
  Random intercepts for subject 
*/
  
data {
  int<lower=1> nS;              //number of subjects
  int nB[nS];                 // total number of bigrams produced by ppt
  int<lower=1> maxB;          // Max number of bigrams for matrix
  matrix[nS,maxB] y;            //outcome: for each subject one IKI per column/bigram (in order as produced)
}

transformed data{
  int N = sum(nB) - (nS*K); // This is for avoiding matrix for y_tilde and log_log
}

parameters {
  real alpha;
  vector beta;
  real<lower=0> sigma;		// residual sd

  // For random effects
  vector[nS] u; //subject intercepts
  real<lower=0> sigma_u;	// subj sd
}


model {
  alpha ~ cauchy(0, 100);

  // Priors
  beta ~ cauchy(0, 10);
  sigma ~ cauchy(0, 5);
  
  // REs priors
  sigma_u ~ normal(0,2.5);
  u ~ normal(0, sigma_u); //subj random effects
  
  // Likelihood	
  for(s in 1:nS){
    int nBS = nB[s];
    for(b in 2:nBS){
        y[s,b] ~ normal(alpha + beta * y[s, b-1] + u[s], sigma); 
    }
  }
}

generated quantities{
  vector[N] log_lik;
  vector[N] y_tilde;
  int n2 = 0;
  
  for(s in 1:nS){
    int nBS = nB[s];
    for(b in (K+1):nBS){
        log_lik[n2] = normal_lpdf(y[s,b] | alpha + beta * y[s, b-1] + u[s], sigma);
        y_tilde[n2] = normal_rng(alpha + beta * y[s, b-1] + u[s], sigma);
      }
    }
  }
}

