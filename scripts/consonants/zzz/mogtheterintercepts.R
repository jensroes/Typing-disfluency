# Load packages
library(brms)
library(tidyverse)

# Load df
d <- read_csv("data/ct.csv") %>% 
  select(bigram, component, subj, IKI, target, sex, age, session) %>% 
  filter(component == "Consonants", !is.na(age), age <=16, session == 1) %>%
  select(-component, -session) %>%
  filter(!is.na(IKI), IKI > 0, target == 1) %>%
  select(-target, -sex, -age) %>%
  mutate(subj = as.integer(factor(subj))) ; d


formula <- bf(IKI ~ 1 + (1 | subj) + (1 | bigram), 
              theta2 ~ 1 + (1 | subj) + (1 | bigram)) 

# Mixture model
mix <- mixture(lognormal, lognormal)

# Priors
prior <- c(
  prior(cauchy(0, 1.5), Intercept, dpar = mu1),
  prior(cauchy(0, 1.5), Intercept, dpar = mu2),
  prior(beta(2, 2), Intercept, dpar = theta2)
  )

m <- brm(formula,
         data = d, 
         chains = 3, 
         cores = 3, 
         iter = 4000,
         family = mix
         #prior = prior
         )

plot(m, pars = "^b_")
summary(m)
pp_check(m)
loo(m)

# Save posterior samples
saveRDS(m,
        file="stanout/mogthetaintercepts.rda",
        compress="xz")
