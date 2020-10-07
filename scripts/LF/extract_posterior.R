library(rstan)
library(loo)
library(tidyverse)
library(magrittr)

path <- "stanout/LF/"

(files <- dir(path, pattern = ".rda"))
m <- readRDS(paste0(path, files[5] ))
params <- names(m)[grepl("beta|prob|delta", names(m))]
params <- params[!grepl("beta2|tilde|_raw|_sigma|_mu",params)]
ps <- rstan::extract(m, params) %>% as_tibble()

write_csv(ps, "results/LF_posterior_MoG.csv")

#ytilde <- rstan::extract(m, "y_tilde") %>% as.data.frame() %>% as_tibble()
#write_csv(ytilde, "results/consonants_posterior_MoG_ytilde.csv")