library(loo)
library(tidyverse)
library(magrittr)

path <- "stanout/tapping/"
#path <- "stanout/"

(files <- dir(path, pattern = ".rda"))
ms <- gsub(files, pattern = ".rda", replacement = "")

for(i in 1:length(files)){
  m <- readRDS(paste0(path, files[i] ))
  log_lik <- extract_log_lik(m, merge_chains = F) 
  r_eff <- relative_eff(exp(log_lik)) 
  assign(paste0("loo_",ms[i]), loo(log_lik, r_eff = r_eff, cores = 2))#,  
  print(ms[i]); rm(list = "m"); gc()
}

(loos <- ls(pattern = "loo_"))
mc <- do.call(loo_compare, lapply(loos, as.name))

loos %>%
  tibble(name = ., model = paste0("model",1:length(.))) -> names

mc %<>% as.data.frame() %>%
  round(2) %>%
  mutate(model=row.names(.)) %>%
  select(model, elpd_diff:se_elpd_loo) %>%
  remove_rownames() %>%
  left_join(names) %>%
  select(-model) %>%
  mutate(name = gsub("loo_", "", name)) %>%
  rename(model = name);mc

file_out <- paste0("results/loo_results_tapping.csv")
write_csv(mc, file_out)
