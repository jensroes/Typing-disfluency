# Load packages
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
  filter(component == "Tapping") %>% 
  select(-component)

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
  #  nBre <- max(as.integer(factor(d$bigram)))
  maxB <- maxB
  y <- y
  K <- 1
  N <- N
} );str(dat)


# Initialise start values
start <- 
  function(chain_id = 1){
    list( phi = array(0)
          , beta_mu = 5
          , beta_sigma = .1
          , beta_raw = 0
          , sigma = 1
          , u = rep(0.1, dat$nS)
          , sigma_u = 1
    )
  }

start_ll <- lapply(1:n_chain, function(id) start(chain_id = id) )

# --------------
# Stan models ##
# --------------
# Load model
ar1 <- stan_model(file = "stanin/ARK.stan")
# Parameters to omit in output
#omit <- c("RE")

# Fit model
m <- sampling(ar1, 
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
              seed = 81, control = list(adapt_delta = 0.99,
                                        max_treedepth = 16))


# Save model
saveRDS(m, 
        file = "stanout/tapping/ARK1.rda",
        compress = "xz")

# Traceplots
#names(m)
param <- c("beta",  "phi", "sigma", "beta_mu", "beta_raw", "beta_sigma") 
summary(print(m, pars = param, probs = c(.025,.975)))
traceplot(m, param, inc_warmup = F)
#traceplot(m, "u", inc_warmup = F)

