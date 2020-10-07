# Load packages
rm(list=ls());gc()
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
  K <- 1
  N <- N
} );str(dat)


# Initialise start values
start <- 
  function(chain_id = 1){
    list( phi = array(0)
          , phi_s = matrix(0, nrow = dat$nS, ncol = dat$K)
          , beta_mu = 5
          , beta_sigma = .1
          , beta_raw = 0
          , sigma = .1
          , tau = .1
          , u = rep(0.1, dat$nS)
          , sigma_u = 1
    )
  }

start_ll <- lapply(1:n_chain, function(id) start(chain_id = id) )

# --------------
# Stan models ##
# --------------
# Load model
ark <- stan_model(file = "stanin/ARKppt.stan")

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
              thin = 1,
              seed = 81,
              control = list(max_treedepth = 16,
                             adapt_delta = 0.99,
                             stepsize = 2)
)

# Save model
saveRDS(m, 
        file = "stanout/consonants/ARK1ppt.rda",
        compress = "xz")

# Traceplots
(param <- names(m)[!grepl("log_lik|y_tilde|phi_s|u\\[",names(m))])
summary(print(m, pars = param, probs = c(.025,.975)))
traceplot(m, param, inc_warmup = F)
#traceplot(m, "u", inc_warmup = F)

