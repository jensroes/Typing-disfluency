lf <- read_csv("../results/LF_posterior_ARKMoG.csv")
cons <- read_csv("../results/consonants_posterior_ARKMoG.csv")

#lf <- read_csv("results/LF_posterior_ARKMoG.csv")
#cons <- read_csv("results/consonants_posterior_ARKMoG.csv")

lf$Comp <- "LF"
cons$Comp <- "Consonants"

ps <- bind_rows(lf, cons)

ps %>% select(starts_with("phi"), Comp) %>%
  gather(Param, value, -Comp) %>%
    group_by(Param, Comp) %>%
    summarise(M = mean(value),
            lo = quantile(value, probs = .025),
            up = quantile(value, probs = .975),
            p = mean(value < .5)) %>%
    ungroup() %>%
    separate(Param, into = c("Param", "id"), sep = "\\[") %>%
    mutate(id = gsub("\\]", "", id),
           id = gsub(",1", "", id),
         id = ifelse(is.na(id), 0, id)) -> d_phi

d_phi %>% filter(Param == "phi") %>% pivot_longer(M:up) %>% pull(value) -> phi_sum


ps %>% select(starts_with("prob"), Comp) %>%
  gather(Param, value, -Comp) %>%
  group_by(Param, Comp) %>%
  summarise(M = mean(value),
            lo = quantile(value, probs = .025),
            up = quantile(value, probs = .975),
            p = mean(value < .5)) %>%
  ungroup() %>%
  separate(Param, into = c("Param", "id"), sep = "\\[") %>%
  mutate(id = gsub("\\]", "", id),
         id = ifelse(is.na(id), 0, id)) -> d_theta

d_theta %>% filter(id == 0) %>% pivot_longer(M:p) %>% pull(value) -> theta_sum


ps %>% select(contains("beta"), contains("phi"), 
              Comp, -`phi[1]`, -beta, -tau_phi, -starts_with("beta2")) %>%
  gather(Param, value, -Comp) %>%
  separate(Param, into = c("Param", "id"), sep = "\\[") %>%
  mutate(id = gsub("\\]", "", id),
         id = gsub(",1", "", id)) %>%
  group_by(Param, id) %>%
  mutate(row_id = 1:n()) %>%
  ungroup() %>%
  pivot_wider(names_from = Param, values_from = value) %>%
  mutate(beta_s = ifelse(Comp == "LF", beta_s + (beta_s * phi_s), beta_s)) %>%
  mutate_at(vars(starts_with("beta")), exp) %>%
  rename(value = beta_s) %>%
  group_by(id, Comp) %>%
  summarise(M = mean(value),
            lo = quantile(value, probs = .025),
            up = quantile(value, probs = .975)) %>%
  ungroup() -> d_beta

ps %>% select(beta, Comp) %>%
  mutate_at(vars(starts_with("beta")), exp) %>%
  rename(value = beta) %>%
  group_by(Comp) %>%
  summarise(M = mean(value),
            lo = quantile(value, probs = .025),
            up = quantile(value, probs = .975)) %>%
  ungroup() -> d_beta2

d_beta2 %>% pivot_longer(M:up) %>% pull(value) %>% round(0) -> beta_sum

ps %>% select(starts_with("prob_s"), starts_with("beta_s"), Comp) %>%
  mutate_at(vars(starts_with("beta_s")), exp) %>%
  pivot_longer(-Comp) %>%
  separate(name, into = c("param", "subj"), sep = "_") %>%
  group_by(param, subj) %>%
  mutate(rep = 1:n()) %>%
  group_by(subj, param, Comp) %>%
  summarise(M = mean(value),
            lo = quantile(value, probs = .025),
            up = quantile(value, probs = .975)) %>%
  pivot_wider(names_from = param, values_from = c(M, lo, up)) -> pss

ps %>% select(beta, delta, Comp) %>%
  mutate(delta = exp(beta+delta) - exp(beta)) %>%
  group_by(Comp) %>%
  summarise(mean(delta), quantile(delta, probs = .025), quantile(delta, probs = .975)) %>%
  pivot_longer(-Comp) %>% pull(value) %>% round(0) -> delta_sum
