# Load packages
library(tidyverse)
library(rstan)
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
  select(subj, bg, bigram, IKI)

# Prepare data and variables
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
} );str(dat)


# Initialise start values
start <- 
  function(chain_id = 1){
    list(   beta_mu = 5
          , beta_sigma = .1
          , beta_raw = 0
          , delta_s = rep(.1, dat$nS)
          , theta_s = rep(0, dat$nS)
          , theta = 0
          , delta = .1
          , tau_theta = .1
          , tau_delta = .1
          , sigma = 1
          , sigma_diff = .1
          , u = rep(0, dat$nS)
          , w = rep(0, dat$maxB)
          , sigma_u = 0.1
          , sigma_w = 0.1
    )
  }

start_ll <- lapply(1:n_chain, function(id) start(chain_id = id) )

# --------------
# Stan models ##
# --------------
#---- 
# Load model
mog <- stan_model(file = "stanin/MoGpptsbigraminterceptsdelta.stan")

# Parameters to omit in output
omit <- c("theta", "theta_s", "log_theta_s", "u")

# Fit model
m <- sampling(mog, 
              data = dat,
              init = start_ll,
              iter = iterations,
              warmup = iterations/2,
              chains = n_chain, 
              cores = n_cores,
              refresh = 2000,
              save_warmup = FALSE, # Don't save the warmup
              include = FALSE, # Don't include the following parameters in the output
              pars = omit,
              thin = 1,
              seed = 81,
              control = list(max_treedepth = 16,
                             adapt_delta = 0.99,
                             stepsize = 2)
)

# Save model
saveRDS(m, 
        file = "stanout/consonants/MoGpptsbigraminterceptsdelta.rda",
        compress = "xz")

# Traceplots
(param <- names(m)[!grepl("log_|y_|n|lp_|_s|w", names(m))])
#param <- c("beta", "delta", "theta", "sigma", "sigma_u", "sigma_w") 
summary(print(m, pars = param, probs = c(.025,.975)))
traceplot(m, param, inc_warmup = F)
#traceplot(m, "u", inc_warmup = F)

