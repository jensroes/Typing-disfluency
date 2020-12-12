library(tidyverse)
library(magrittr)

# Load posterior
path <- "stanout/consonants/"
(files <- dir(path))
m <- readRDS(paste0(path, files[4]))
y_pred <- rstan::extract(m, 'y_tilde'); rm(list = "m")

y_pred_cons <- y_pred$y_tilde %>% t() %>% as_tibble() %>%
  select(paste0("V",1:50)) %>%
  gather(iteration, y_tilde) 

path <- "stanout/LF/"
(files <- dir(path))
m <- readRDS(paste0(path, files[4]))
y_pred <- rstan::extract(m, 'y_tilde'); rm(list = "m")

y_pred_lf <- y_pred$y_tilde %>% t() %>% as_tibble() %>%
  select(paste0("V",1:50)) %>%
  gather(iteration, y_tilde) 

# Load posterior (LMM)
path <- "stanout/consonants/"
(files <- dir(path))
m <- readRDS(paste0(path, files[3]))
y_pred <- rstan::extract(m, 'y_tilde'); rm(list = "m")

y_pred_cons_lmm <- y_pred$y_tilde %>% t() %>% as_tibble() %>%
  select(paste0("V",1:50)) %>%
  gather(iteration, y_tilde) 

path <- "stanout/LF/"
(files <- dir(path))
m <- readRDS(paste0(path, files[3]))
y_pred <- rstan::extract(m, 'y_tilde'); rm(list = "m")

y_pred_lf_lmm <- y_pred$y_tilde %>% t() %>% as_tibble() %>%
  select(paste0("V",1:50)) %>%
  gather(iteration, y_tilde) 

# Save the predicted samples and data before plotting 
y_pred_cons %<>% mutate(Comp = "Consonants", Model = "MoG")
y_pred_lf %<>% mutate(Comp = "LF bigrams", Model = "MoG")
y_pred_lf_lmm %<>% mutate(Comp = "LF bigrams", Model = "LMM")
y_pred_cons_lmm %<>% mutate(Comp = "Consonants", Model = "LMM")

y_pred <- bind_rows(y_pred_cons, y_pred_lf,y_pred_cons_lmm, y_pred_lf_lmm) %>%
  mutate(type = paste(expression(tilde(y)))) %>%
  rename(y = y_tilde)

write_csv(y_pred, "results/y_tilde_LMM_MoG.csv")
