get_data <- function(path){
  library(tidyverse)
  d <- read_csv(paste0(path, "ct.csv")) %>% 
    filter(component %in% c("Consonants", "Tapping"), 
           !is.na(age), 
           age >= 18,
           age <= 25,
           session == 1) %>%
    mutate(subj = as.integer(factor(paste0(subj, logfile, session))),
           bg = bigram
    ) %>% 
    filter(!is.na(IKI), 
           IKI > 0 
           #target == 1
           ) %>%
    select(-logfile) %>%
    select(subj, bg, bigram, IKI, sex, age, component, target) %>%
    group_by(subj, component) %>%
    mutate(bigram = 1:n()) %>%
    ungroup() 

  d %>% count(subj, component) %>% 
    count(subj) %>% 
    arrange(n) %>% 
    filter(n == 1) %>% 
    pull(subj) -> remove_subj_1
  
  d %>% group_by(subj) %>%
    mutate(maxB = max(bigram)) %>%
    filter(maxB > 30) -> remove_subj_2
  
  set.seed(125)
  d <- d %>% filter(!(subj %in% c(remove_subj_1, remove_subj_2))) %>%
    filter(subj %in% sample(unique(subj), 100)) %>%
    mutate(subj = as.integer(factor(subj))) %>%
    arrange(subj, component, bigram) 
    
  return(d)
}

