/*
  Serial order of bigrams autoregression model
  Mixture model of first order and second order autoregression model
  Random intercepts for subject 
*/
  
data {
  int<lower=1> nS;              //number of subjects
  int nB[nS];                 // total number of bigrams produced by ppt
  int<lower=1> maxB;          // Max number of bigrams for matrix
  matrix[nS,maxB] y;            //outcome: for each subject one IKI per column (in order as produced)
}

transformed data{
  int N = sum(nB) - nS;
}

parameters {
  vector[2] alpha; // First and second order parameters
  vector[2] beta;
  real gamma;
	simplex[2] theta; // mixing proportion
  
  real<lower=0> sigma;		// residual sd

  // For random effects
  vector[nS] u; //subject intercepts
  real<lower=0> sigma_u;	// subj sd
}


transformed parameters{
	vector[2] log_theta;

	log_theta[1] = log(theta[1]);
  log_theta[2] = log1m(theta[1]);

}

model {
  int n = 0;
  vector[2] lp_parts;

  // Priors
  alpha ~ cauchy(0, 10);
  beta ~ cauchy(0, 5);
  gamma ~ cauchy(0, 5);
  sigma ~ cauchy(0, 2.5);
  theta ~ beta(2, 2);
  // REs priors
  sigma_u ~ normal(0,2.5);
  u ~ normal(0, sigma_u); //subj random effects
  
  // Likelihood	
  for(s in 1:nS){
    int nBS = nB[s];
    for(b in 2:nBS){ 
      n += 1;
      if(b >= 3){
        lp_parts[1] = log_theta[2] + normal_lpdf(y[s,b] | alpha[2] + beta[2]*y[s, b-1] + gamma*y[s, b-2] + u[s], sigma);
      }
      else{
        lp_parts[1] = log_theta[1] + normal_lpdf(y[s,b] | alpha[1] + beta[1]*y[s, b-1] + u[s], sigma);
      }
      lp_parts[2] = log_theta[1] + normal_lpdf(y[s,b] | alpha[1] + beta[1]*y[s, b-1] + u[s], sigma);
      target += log_sum_exp(lp_parts); // play with the indecies in lpparts to get this working 
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
    for(b in 2:nBS){
      n += 1;
      if(b >= 3){
	    log_lik[n] = log_sum_exp(
 		    log_theta[2] + normal_lpdf(y[s,b] | alpha[2] + beta[2]*y[s, b-1] + gamma*y[s, b-2] + u[s], sigma), 
   		  log_theta[1] + normal_lpdf(y[s,b] | alpha[1] + beta[1]*y[s, b-1] + u[s], sigma)
	    );
       theta_tilde = bernoulli_rng(theta[2]); 
          if(theta_tilde) { 
              y_tilde[n] = normal_rng(alpha[2] + beta[2]*y[s, b-1] + gamma*y[s, b-2] + u[s], sigma);
          }
          else{
              y_tilde[n] = normal_rng(alpha[1] + beta[1]*y[s, b-1] + u[s], sigma);
          }
      }
      else{
        log_lik[n] = log_sum_exp(
 		      log_theta[1] + normal_lpdf(y[s,b] | alpha[1] + beta[1]*y[s, b-1] + u[s], sigma), 
   		    log_theta[1] + normal_lpdf(y[s,b] | alpha[1] + beta[1]*y[s, b-1] + u[s], sigma)
	      );
        y_tilde[n] = normal_rng(alpha[1] + beta[1]*y[s, b-1] + u[s], sigma);
      }
    } 
  }
}

