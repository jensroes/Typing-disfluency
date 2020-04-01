# Load packages
library(tidyverse)
library(rstan)
library(loo)
library(plyr)
library(magrittr)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

# Sampling parameters
set.seed(125)
n_cores = 3
n_chain = 3
iterations = 6000

# Load df
d <- read_csv("data/ct.csv") %>% 
  filter(component == "Consonants", 
         !is.na(age), 
         age >= 50, 
         session == 1,
         !is.na(IKI), 
         IKI > 0, 
         target == 1) %>%
  select(bigram, component, IKI, subj) %>%
  mutate(subj = as.integer(factor(subj)),
         mutate(component = factor(component, levels = c("Tapping", "Sentence", "HF", "LF", "Consonants"), ordered = TRUE))) ; d

# Data as list for stan input
components <- d %$% as.numeric(mapvalues(component, from = unique(component), to= 1:length(unique(component))))
dat <-  within( list(), {
  N <- nrow(d)
  y <- d$IKI 
  subj <- d$subj
  S <- length(unique(d$subj))
  bigram <- as.integer(factor(d$bigram))
  B <- length(unique(d$bigram))
} );str(dat)


# Initialise start values
start <- 
  function(chain_id = 1){
    list(beta = mean(log(dat$y))
         , delta = 1
         , sigma = sd(log(dat$y))
         , sigma_diff = .01 
         , theta = rep(.5, 2)
         , u = rep(0.1, dat$S)
         , w = rep(0.1, dat$B)
         , sigma_u = 1
         , sigma_w = 1
         , alpha = chain_id
    )
  }

start_ll <- lapply(1:n_chain, function(id) start(chain_id = id) )

# --------------
# Stan models ##
# --------------
#---- 
# Mixture of two gaussians with unequal variance
#---- 
# Load model
mog <- stan_model(file = "stanin/MoG.stan")

# Parameters to omit in output
omit <- c("mu1", "mu2")

# Fit model
m <- sampling(mog, 
              data = dat,
              init = start_ll,
              iter = iterations,
              warmup = iterations/2,
              chains = n_chain, 
              cores = n_cores,
	            refresh = 250,
              save_warmup = FALSE, # Don't save the warmup
              include = FALSE, # Don't include the following parameters in the output
              pars = omit,
              thin = 1,
              seed = 81,
              control = list(max_treedepth = 16,
                             adapt_delta = 0.99)
)


# Save model
saveRDS(m, 
        file = "stanout/MoG.rda",
        compress = "xz")

# Traceplots
#param <- c("beta", "delta", "theta") 
#summary(print(m, pars = param))
#traceplot(m, param, inc_warmup = TRUE)
#traceplot(m, param, inc_warmup = F)

# Extract and save posterior and log likelihood seperately
# Get log likelihood
log_lik <- extract_log_lik(m, merge_chains = F) 
saveRDS(log_lik, 
        file = "stanout/MoG_loglik.rda",
        compress = "xz")


# Get parameter posterior
param <- c("beta","beta2", "delta", "theta", "sigma", "sigma_diff", "sigmap_e", "sigma_e")
samps <- as.data.frame(m, pars = param) 
saveRDS(samps, 
        file = "stanout/MoG_posterior.rda",
        compress = "xz")


# Get random effects
param <- c("u", "w", "sigma_u", "sigma_w")
re <- as.data.frame(m, pars = param) 
saveRDS(re, 
        file = "stanout/MoG_re.rda",
        compress = "xz")


# Get posterior predicted values
param <- c("y_tilde")
y_tilde <- as.data.frame(m, pars = param) 
saveRDS(y_tilde, 
        file = "stanout/MoG_ytilde.rda",
        compress = "xz")



