# Load packages
rm(list=ls())
library(tidyverse)
library(rstan)
library(loo)
library(magrittr)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
source("./functions/get_data.R")

# Sampling parameters
n_cores = 3
n_chain = 3
iterations = 30000

# Load df
path <- "./data/"
d <- get_data(path = path) %>% 
  filter(component == "Consonants") %>% 
  select(-component)

(maxB <- max(d$bigram))
(nS <- length(unique(d$subj)))
subj <- d$subj
bigram <- d$bigram
N <- nrow(d)
y <- array(-1, c(nS, maxB))
correct <- array(-1, c(nS, maxB))
d %>% group_by(subj) %>% dplyr::count() %>% ungroup() -> nB

for (i in 1:nS) {
  subj_data <- filter(d, subj == i)
  y[i, 1:nB[nB$subj == i,]$n] <- subj_data$IKI  
  correct[i, 1:nB[nB$subj == i,]$n] <- subj_data$target
}

dat <- within( list(), {
  nS <- nS
  nB <- nB$n
  #  nBre <- max(as.integer(factor(d$bigram)))
  maxB <- maxB
  y <- y
  correct <- correct
  N <- N
  K <- 1
} );str(dat)


# Initialise start values
start <- 
  function(chain_id = 1){
    list( phi = as.array(0)
          , delta_corr = 0.1
          , delta_incorr = 0.1
          , beta_mu = 5
          , beta_raw = 0
          , beta_sigma =.1
          , sigma = 2
          , theta_corr = .5
          , theta_incorr = .5
          , sigma_diff_corr = .1
          , sigma_diff_incorr = .1
          , u = rep(0.1, dat$nS)
          , sigma_u = 1
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
ark <- stan_model(file = "stanin/ARKtargetMoG.stan")

# Parameters to omit in output
#omit <- c("RE")

# Fit model
m <- sampling(ark, 
              data = dat,
              init = start_ll,
              iter = iterations,
              warmup = iterations/2,
              chains = n_chain, 
              cores = n_cores,
              refresh = 2000,
              save_warmup = FALSE, # Don't save the warmup
#              include = FALSE, # Don't include the following parameters in the output
#              pars = omit,
              thin = 1,
              seed = 81,
              control = list(max_treedepth = 16,
                             adapt_delta = 0.99,
                             stepsize = 2)
)


# Save model
saveRDS(m, 
        file = "stanout/consonants/ARKtargetMoG.rda",
        compress = "xz")

# Traceplots
#names(m)
param <- c("beta",  "phi", "theta_corr", "theta_incorr", "delta_corr", 
           "delta_incorr", "sigma", "sigma_diff_corr", "sigma_diff_incorr", 
           "sigma_e", "sigmap_e_corr", "sigmap_e_incorr") 
summary(print(m, pars = param, probs = c(.025,.975)))
traceplot(m, param, inc_warmup = F)
#traceplot(m, "u", inc_warmup = F)

