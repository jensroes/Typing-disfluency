ps <- read_csv("../results/LF_posterior_MoG.csv")

ps %>% select(starts_with("prob")) %>%
  gather(Param, value) %>%
  group_by(Param) %>%
  summarise(M = mean(value),
            lo = quantile(value, probs = .025),
            up = quantile(value, probs = .975),
            p = mean(value < .5)) %>%
  ungroup() %>%
  separate(Param, into = c("Param", "id"), sep = "\\[") %>%
  mutate(id = gsub("\\]", "", id),
         id = ifelse(is.na(id), 0, id)) %>%
  arrange(M) %>%
  mutate(id2 = ifelse(id == 0, 1, 0),
         id2 = 1:n(),
         id2 = factor(id2, levels = unique(id2)[id2], ordered = T)) -> d_theta

d_theta %>% filter(id == 0) %>% pivot_longer(M:p) %>% pull(value) -> theta_sum

ps %>% select(starts_with("beta")) %>%
  select(-starts_with("beta2")) %>%
  gather(Param, value) %>%
  mutate(value = exp(value)) %>%
  group_by(Param) %>%
  summarise(M = mean(value),
            lo = quantile(value, probs = .025),
            up = quantile(value, probs = .975)) %>%
  ungroup() %>%
  separate(Param, into = c("Param", "id"), sep = "\\[") %>%
  mutate(id = gsub("\\]", "", id),
         id = ifelse(is.na(id), 0, id)) %>%
  arrange(M) %>%
  mutate(id2 = ifelse(id == 0, 1, 0),
         id2 = 1:n(),
         id2 = factor(id2, levels = unique(id2)[id2], ordered = T)) -> d_beta

d_beta %>% filter(id == 0) %>% pivot_longer(M:up) %>% pull(value) -> beta_sum

ps %>% select(starts_with("prob_s"), starts_with("beta_s")) %>%
  mutate_at(vars(starts_with("beta_s")), exp) %>%
  pivot_longer(everything()) %>%
  separate(name, into = c("param", "subj"), sep = "_") %>%
  group_by(param, subj) %>%
  mutate(rep = 1:n()) %>%
  group_by(subj, param) %>%
  summarise(M = mean(value),
            lo = quantile(value, probs = .025),
            up = quantile(value, probs = .975)) %>%
  pivot_wider(names_from = param, values_from = c(M, lo, up)) -> pss


ps %>% select(beta, delta) %>%
  transmute(delta = exp(beta+delta) - exp(beta)) %>%
  summarise(mean(delta), quantile(delta, probs = .025), quantile(delta, probs = .975)) %>%
  pivot_longer(everything()) %>% pull(value) %>% round(0) -> delta_sum
