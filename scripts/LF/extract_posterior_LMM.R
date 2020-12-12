library(tidyverse)

path <- "stanout/LF/"
(files <- dir(path, pattern = ".rda"))

m <- readRDS(paste0(path, files[3]))
(params <- names(m)[grepl("beta|u", names(m))])
(params <- params[!grepl("_", params)])
ps <- rstan::extract(m, params) %>% as_tibble()

names(ps) <- gsub("u", "beta_s", names(ps))

write_csv(ps, "results/LF_posterior_LMM.csv")
