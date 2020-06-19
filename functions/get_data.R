get_data <- function(path){
  library(tidyverse)
  d <- read_csv(paste0(path, "ct.csv")) %>% 
    filter(!is.na(age), 
           age >= 18,
           age <= 25,
           session == 1) %>%
    mutate(subj = as.integer(factor(paste0(subj, logfile, session))),
           bg = bigram
    ) %>% 
    filter(!is.na(IKI), 
           IKI > 0 
           ) %>%
    select(-logfile) %>%
    select(subj, bg, bigram, IKI, sex, age, component, target, rep, freq) %>%
    mutate(rep = ifelse(component %in% c("Tapping","Consonants"), 1, rep)) %>%
    group_by(subj, component, rep) %>%
    mutate(bigram = 1:n()) %>%
    ungroup() 

  maxC <- d %>% pull(component) %>% unique() %>% length()
  
  d %>% count(subj, component) %>% 
    count(subj) %>% 
    arrange(n) %>% 
    filter(n < maxC) %>% 
    pull(subj) -> remove_subj_1
  
  d %>% 
    filter(component != "Tapping") %>%
    group_by(subj, component) %>%
    summarise(maxB = max(bigram)) %>%
    ungroup() %>% 
    group_by(component) %>%
    mutate(medianB = mean(maxB)) %>%
    arrange(desc(maxB),component,subj) %>% 
    filter(maxB > 1.75*medianB) %>%
    pull(subj) -> remove_subj_2
  
  set.seed(125)
  d <- d %>% filter(!(subj %in% c(remove_subj_1, remove_subj_2))) %>%
    filter(subj %in% sample(unique(subj), 250)) %>%
    mutate(subj = as.integer(factor(subj))) %>%
    arrange(subj, component, bigram) 
    
  return(d)
}

