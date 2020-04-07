library(loo)
library(tidyverse)
library(magrittr)

path <- "stanout/consonants/"
#path <- "stanout/"

(files <- dir(path, pattern = ".rda"))
(files <- files[grepl("MoGppts", files)])
ms <- gsub(files, pattern = ".rda", replacement = "")

for(i in 1:length(files)){
  m <- readRDS(paste0(path, files[i] ))
  log_lik <- extract_log_lik(m, merge_chains = F) 
  r_eff <- relative_eff(exp(log_lik)) 
  assign(paste0("loo_",ms[i]), loo(log_lik, r_eff = r_eff, cores = 2))  
  print(ms[i]); rm(list = "m"); gc()
}

(loos <- ls(pattern = "loo_MoGppts"))
mc <- do.call(what = loo_compare, args = lapply(loos, as.name))
tibble(name = loos, model = paste0("model",1:length(loos))) -> names

mc %<>% as.data.frame() %>%
  round(2) %>%
  mutate(model=row.names(.)) %>%
  select(model, elpd_diff:se_elpd_loo) %>%
  left_join(names) %>%
  select(-model) %>%
  mutate(name = gsub("loo_", "", name)) %>%
  select(name, elpd_diff:se_elpd_loo) %>%
  rename(Model = name);mc

mc %<>%
  filter(elpd_diff != 0) %>%
  mutate(Model = paste0(Model," - MoGpptsbigraminterceptsdelta"))


loo_compare(loo_MoGpptsbigramintercepts, loo_MoGpptsbigraminterceptsdelta2) %>%
  as.data.frame() %>%
  mutate(Model=row.names(.)) %>%
  select(Model, elpd_diff:se_elpd_loo) %>%
  filter(elpd_diff != 0) %>%
  mutate(Model = "MoGpptsbigraminterceptsdelta2 - MoGpptsbigramintercepts") %>%
  bind_rows(mc) -> mc;mc

file_out <- paste0("results/loo_results_consonants_mog.csv")
write_csv(mc, file_out)
