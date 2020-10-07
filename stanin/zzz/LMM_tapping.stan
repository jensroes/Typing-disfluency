// LMM with REs for subject and bigrams and by-subjects slopes

data {
	int<lower=1> N;                    // Number of observations
	vector[N] y;  		            //outcome
  int<lower =2> max_order;
  int<lower=1, upper=max_order> order[N];
	int<lower=1> S;                  //number of subjects
	int<lower=1, upper=S> subj[N];   //subject id
	int<lower=1> B;                  //number of bigrams
	int<lower=1, upper=B> bigram[N];   //bigramid
}

transformed data{
  vector[N] logy = log(y);
}


parameters {
//  real<lower=0> beta;
  real delta;

	real<lower=0> sigma;		// residual sd
	
  // For random effects
	vector[S] u; //subj intercepts
  real<lower=0> sigma_u;//subj sd
  
  vector[B] w; //bigram intercepts
  real<lower=0> sigma_w;//bigram sd

	// Parameters for non-centering
	real<lower = 0> mu_beta;
  real<lower = 0> sigma_beta;	
	real beta_raw;			// distributions
  
}

transformed parameters{
  vector[N] mu;
	real beta = mu_beta + sigma_beta * beta_raw;

  for(n in 1:N){
    mu[n] = beta + order[n]*delta + u[subj[n]] + w[bigram[n]];
  }
}

model {
  // Priors
  mu_beta ~ cauchy(5, 1.5);
  sigma_beta ~ cauchy(0, 1);
  beta_raw ~ normal(0, 2);
  //beta ~ cauchy(5, 1.5);
	sigma ~ cauchy(0, 2.5);
  delta ~ cauchy(0, 1.5);
	
	// REs priors
  sigma_u ~ normal(0,2.5);
  u ~ normal(0, sigma_u); //subj random effects

  sigma_w ~ normal(0,2.5);
  w ~ normal(0, sigma_w); //bigram random effects

  // likelihood
  logy ~ normal(mu, sigma);
}

generated quantities{
	real log_lik[N];
	vector[N] y_tilde;
	for (n in 1:N){
		log_lik[n] = normal_lpdf(logy[n] | mu[n], sigma);
		y_tilde[n] = normal_rng(mu[n], sigma);
	}
}
