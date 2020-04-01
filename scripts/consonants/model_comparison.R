library(loo)
library(tidyverse)
library(magrittr)

path <- "stanout/consonants/"
#path <- "stanout/"

(files <- dir(path, pattern = ".rda"))
ms <- gsub(files, pattern = ".rda", replacement = "")

for(i in 1:length(files)){
  m <- readRDS(paste0(path, files[i] ))
  log_lik <- extract_log_lik(m, merge_chains = F) 
  r_eff <- relative_eff(exp(log_lik)) 
  assign(paste0("loo_",ms[i]), loo(log_lik, r_eff = r_eff, cores = 2))  
  print(ms[i]); rm(list = "m"); gc()
}

(loos <- ls(pattern = "loo_")[-c(1,3,7,9,10)])
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


loo_model_weights(list(loo_MoGpptsbigramintercepts, loo_MoGpptsbigraminterceptsdelta, loo_ARKMoGppt))

#print(loo_ARKtarget)
#plot(loo_ARKtarget)

# Compare all models individually
#Loos <- combn(loos, 2)
#mc2 <- tibble()
#for(i in 1:ncol(Loos)){
#  do.call(what = loo_compare, 
#          args = lapply(Loos[c(1,2),i], as.name)) %>%
#    as.data.frame() %>%
 #   rownames_to_column() %>%
#    rename(worse_model = rowname) %>%
#    filter(elpd_diff != 0)  %>%
#    mutate(model1 = Loos[1,i],
#           model2 = Loos[2,i]) %>%
#    bind_rows(mc2) -> mc2
#}

#mc2 %<>% 
#  mutate_at(vars(model1,model2), ~gsub(pattern = "loo_", replacement = "", .)) %>%
#  mutate(Model = ifelse(worse_model == "model2", 
#                        paste(model1, model2, sep = " - "),
#                        paste(model2, model1, sep = " - "))) %>%
#  select(Model, elpd_diff, se_diff) %>%
#  arrange(desc(elpd_diff));mc2

file_out <- paste0("results/loo_results_consonants.csv")
write_csv(mc, file_out)
#file_out <- paste0("results/loo_comparisons_consonants.csv")
#write_csv(mc2, file_out)
