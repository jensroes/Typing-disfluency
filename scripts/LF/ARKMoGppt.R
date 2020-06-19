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
  filter(component == "LF", rep == 1) %>% 
  select(subj, bg, bigram, IKI)

(maxB <- max(d$bigram))
(nS <- length(unique(d$subj)))
subj <- d$subj
bigram <- d$bigram
N <- nrow(d)
y <- array(-1, c(nS, maxB))
d %>% group_by(subj) %>% dplyr::count() %>% ungroup() -> nB

for (i in 1:nS) {
  subj_data <- filter(d, subj == i)
  y[i, 1:nB[nB$subj == i,]$n] <- subj_data$IKI  
}

dat <- within( list(), {
  nS <- nS
  nB <- nB$n
  maxB <- maxB
  y <- y
  N <- N
  K <- 1
} );str(dat)


# Initialise start values
start <- 
  function(chain_id = 1){
    list( phi = as.array(0)
          , phi_s = matrix(0, nrow = dat$nS, ncol = dat$K)
          , delta = 0.1
          , beta_mu = 5
          , beta_sigma = .1
          , beta_raw = 0
          , theta = 0
          , theta_s = rep(0, dat$nS)
          , tau_theta = .1
          , tau_phi = .1
          , sigma = 1
          , sigma_diff = .1
          , u = rep(0, dat$nS)
          , w = rep(0, dat$maxB-1)
          , sigma_u = 0.1
          , sigma_w = 0.1
    )  }

start_ll <- lapply(1:n_chain, function(id) start(chain_id = id) )

# --------------
# Stan models ##
# --------------
#---- 
# Mixture of two gaussians with unequal variance
#---- 
# Load model
ark <- stan_model(file = "stanin/ARKMoGppt.stan")

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
              refresh = 1000,
              save_warmup = FALSE, # Don't save the warmup
              thin = 1,
              seed = 81,
              control = list(max_treedepth = 16,
                             adapt_delta = 0.99,
                             stepsize = 2)
)


# Save model
saveRDS(m, 
        file = "stanout/LF/ARKMoGppt.rda",
        compress = "xz")

# Traceplots
(param <- names(m)[!grepl("u|w|log_|n|lp_|_tilde|_s", names(m))])
#param <- c("beta", "beta_mu", "beta_sigma", "beta_raw", "phi", "theta", "delta",
#           "sigma", "sigma_diff", "sigma_e", "sigmap_e") 
summary(print(m, pars = param, probs = c(.025,.975)))
traceplot(m, param, inc_warmup = F)
#traceplot(m, "u", inc_warmup = F)

