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
iterations = 4000

# Load df
d <- read_csv("data/ct.csv") %>% 
  filter(component == "Consonants", 
        # component == "Tapping", 
         !is.na(age), 
         age >= 50, 
         session == 1,
         !is.na(IKI), 
         IKI > 0, 
         target == 1) %>%
  select(bigram, IKI, subj) %>%
  mutate(subj = as.integer(factor(subj))) %>%
  filter(subj %in% 1:10); d

d %<>% group_by(subj) %>% dplyr::mutate(bid = 1:n())
d %>% group_by(subj) %>% dplyr::count() -> nB
nB %>% summarise(max = max(n)) -> maxB
nS <- length(unique(d$subj))
y <- array(-1, c(nS, maxB))

for (i in 1:nS) {
  subj_data <- filter(d, subj == i)
  y[i, 1:nB[nB$subj == i,]$n] <- subj_data$IKI  
}

dat <- within( list(), {
  nS <- nS
  nB <- nB$n
  maxB <- maxB$max
  y <- y
  K <- 1
} );str(dat)


# Initialise start values
start <- 
  function(chain_id = 1){
    list(sigma = .1
#         , beta = rep(1,1)
         , alpha = 0
         , u = rep(0.1, dat$nS)
         , sigma_u = 1
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
ar1 <- stan_model(file = "stanin/AR1.stan")

# Parameters to omit in output
#omit <- c("mu1", "mu2")

# Fit model
m <- sampling(ar1, 
              data = dat,
              init = start_ll,
              iter = iterations,
              warmup = iterations/2,
              chains = n_chain, 
              cores = n_cores,
	            refresh = 250,
 #             save_warmup = FALSE, # Don't save the warmup
#              include = FALSE, # Don't include the following parameters in the output
 #             pars = omit,
              thin = 1,
              seed = 81,
              control = list(max_treedepth = 16,
                             adapt_delta = 0.99)
)


# Save model
saveRDS(m, 
        file = "stanout/AR1.rda",
        compress = "xz")

# Traceplots
names(m)
param <- c("alpha","beta", "sigma") 
summary(print(m, pars = param))

traceplot(m, param, inc_warmup = TRUE)
traceplot(m, param, inc_warmup = F)


# Extract and save posterior and log likelihood seperately
# Get log likelihood
log_lik <- extract_log_lik(m, merge_chains = F) 

saveRDS(log_lik, 
        file = "stanout/AR1_loglik.rda",
        compress = "xz")


# Get parameter posterior
param <- c("beta","beta2", "delta", "theta", "sigma", "sigma_diff", "sigmap_e", "sigma_e")
samps <- as.data.frame(m, pars = param) 
saveRDS(samps, 
        file = "stanout/AR1_posterior.rda",
        compress = "xz")


# Get random effects
param <- c("u", "sigma_u")
re <- as.data.frame(m, pars = param) 
saveRDS(re, 
        file = "stanout/AR1_re.rda",
        compress = "xz")


# Get posterior predicted values
param <- c("y_tilde")
y_tilde <- as.data.frame(m, pars = param) 
saveRDS(y_tilde, 
        file = "stanout/AR1_ytilde.rda",
        compress = "xz")



