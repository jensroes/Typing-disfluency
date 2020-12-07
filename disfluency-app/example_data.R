# Packages
library(tidyverse)
source("functions/functions.R")
source("functions/get_data.R")

# Load df
path <- "data/"
#d <- get_data(path = path) %>% filter(component == "Consonants") %>% select(-component)
d <- get_data(path = path) %>% filter(component %in% c("Consonants", "LF"), rep == 1, target == 1) %>%
  select(subj, bg, IKI, component) 

write_csv(d, "data/data.csv")
