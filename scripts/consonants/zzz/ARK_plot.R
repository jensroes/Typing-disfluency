library(tidyverse)
library(plyr)
library(magrittr)

# Sampling parameters
set.seed(125)

# Load df
d <- read_csv("data/ct.csv") %>% 
  filter(component == "Consonants", 
         !is.na(age), 
         age >= 50, 
         session == 1,
         !is.na(IKI), 
         IKI > 0, 
         target == 1) %>%
  select(bigram, IKI, subj) %>%
  mutate(subj = as.integer(factor(subj))); d

target <- "tjxggl pgkfkq dtdrtt npwdvf"
target


d %>% group_by(subj) %>%
  dplyr::mutate(id = 1:n()) -> d

d %>% group_by(bigram) %>%
  dplyr::count(bigram) %>%
  filter(n < 80) -> remove

d %>% 
  filter(!(bigram %in% remove$bigram)) %>%
  group_by(bigram) %>%
  dplyr::summarise(M = mean(IKI),
                   lower = M-1.96*sd(IKI)/sqrt(n()),
                   upper = M+1.96*sd(IKI)/sqrt(n())) %>%
  ggplot(aes(x = bigram, y= M, ymax = upper, ymin = lower)) +
  geom_point() +
  geom_errorbar(width = 0, size = .25) +
  geom_line() +
  theme_minimal()

