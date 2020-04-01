library(tidyverse)
library(rstan)
#dir("stanout/consonants")

m <- readRDS("stanout/consonants/MoGbyppt.rda")

m %>% as.data.frame(m)
names(m)
summary(m)                    
