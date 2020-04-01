get_missing_data <- function(path){
  library(tidyverse)
  d <- read_csv(paste0(path, "ct.csv")) %>% 
    filter(component %in% c("Consonants", "Tapping"), 
           !is.na(age), 
           age >= 18,
           age <= 25,
           session == 1) %>%
    mutate(logfile = paste0(subj, logfile),
           subj = as.integer(factor(logfile)),
           bg = bigram
    ) %>% 
    select(-logfile) %>%
    select(subj, bg, bigram, IKI, sex, age, component, target) 
  
  d %>% count(subj, component) %>% 
    count(subj) %>% 
    arrange(n) %>% 
    filter(n == 1) %>% 
    pull(subj) -> remove
  
  set.seed(125)
  d <- d %>% filter(!(subj %in% remove) ) %>%
    filter(subj %in% sample(unique(subj), 500)) %>%
    mutate(subj = as.integer(factor(subj))) %>%
    arrange(subj, component, bigram) %>%
    group_by(subj, component) %>%
    mutate(bigram = 1:n()) %>%
    ungroup() 
  
  d %>%
    mutate(missing_IKI = is.na(IKI), 
           IKI_0 = IKI <= 0, 
           non_target = target == 0) %>%
    group_by(subj,component) %>%
    summarise(missing_IKI = length(which(missing_IKI == TRUE)),
              IKI_0 = length(which(IKI_0 == TRUE)),
              non_target = length(which(non_target == TRUE)),
              N = n()) %>% ungroup() -> filtered_data
    
  return(filtered_data)
}

