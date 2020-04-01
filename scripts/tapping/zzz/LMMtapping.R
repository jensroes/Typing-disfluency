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
iterations = 10000

setwd("Frontline/")

# Load df
d <- read_csv("data/ct.csv") %>% 
  filter(component == "Tapping", 
         !is.na(age), 
         age >= 50, 
         session == 1,
         !is.na(IKI), 
         IKI > 0, 
         target == 1) %>%
  mutate(logfile = paste0(subj, logfile),
         subj = as.integer(factor(logfile)),
         bigram = as.integer(factor(bigram))) %>% 
  group_by(subj) %>%
  dplyr::mutate(order = 1:n()) %>%
  ungroup() %>%
  select(-logfile) %>%
  select(bigram, IKI, subj, order)

S <- length(unique(d$subj))
B <- length(unique(d$bigram))
subj <- d$subj
bigram <- d$bigram
N <- nrow(d)
order <- d$order

dat <- within( list(), {
  S <- S
  B <- B
  y <- d$IKI
  subj <- subj
  bigram <- bigram
  N <- N
  order <- order
  max_order <- max(order)
} );str(dat)


# Initialise start values
start <- 
  function(chain_id = 1){
    list(sigma = .1
         , beta_raw = 0
         , mu_beta = 5
         , sigma_beta = .1
         , delta = 0
         , u = rep(0.1, dat$S)
         , sigma_u = 1
         , w = rep(0.1, dat$B)
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
lmm <- stan_model(file = "stanin/LMMtapping.stan")

# Parameters to omit in output
omit <- c("mu", "logy")

# Fit model
m <- sampling(lmm, 
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
                             adapt_delta = 0.99,
                             stepsize = .1)
)


# Save model
saveRDS(m, 
        file = "stanout/tapping/LMMtapping.rda",
        compress = "xz")

# Traceplots
#names(m)
param <- c("beta", "delta", "sigma") 
summary(print(m, pars = param, probs = c(.025, .975)))
#traceplot(m, param, inc_warmup = F)



