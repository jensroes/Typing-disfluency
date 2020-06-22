library(loo)
library(tidyverse)
library(magrittr)

path <- "stanout/LF/"
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

delta1 <- loo_compare(loo_MoGpptsbigramintercepts, loo_MoGpptsbigraminterceptsdelta)
delta2 <- loo_compare(loo_MoGpptsbigramintercepts, loo_MoGpptsbigraminterceptsdelta2)

delta1 %<>% as.data.frame() %>%
  rownames_to_column(var = "best_model") %>%
  mutate(model = "MoGpptsbigraminterceptsdelta")

delta2 %<>% as.data.frame() %>%
  rownames_to_column(var = "best_model") %>%
  mutate(model = "MoGpptsbigraminterceptsdelta2")

bind_rows(delta1, delta2) %>%
  filter(elpd_diff != 0) %>%
  mutate_if(is.numeric, round, 2) %>%
  select(best_model, model, elpd_diff:se_elpd_loo) -> mc

file_out <- paste0("results/loo_results_LF_mog.csv")
write_csv(mc, file_out)
