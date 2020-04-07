looc <- read_csv("../results/loo_results_consonants.csv") %>%
  mutate_if(is.numeric, round, 0) %>%
  unite(col = "elpd_diff", elpd_diff, se_diff, sep = " (") %>%
  unite(col = "elpd_loo", elpd_loo, se_elpd_loo, sep = " (") %>%
  mutate(elpd_diff = paste0(elpd_diff, ")")) %>%
  mutate(elpd_loo = paste0(elpd_loo, ")")) %>%
  mutate(Model = recode(Model, ARK1ppt = "M3",
                        LMMbigramintercepts = "M1",
                        LMMbigramslopes = "M2",
                        MoGpptsbigramintercepts = "M4",
                        ARKMoGppt = "M5")) %>%
  mutate(Type = recode(Model, M2 = "LMM",
                       M1 = "LMM",
                       M3 = "AR",
                       M4 = "MoG",
                       M5 = "AR + MoG")) %>%
#  mutate(Bigrams = recode(Model, M3 = "Autoregressor $\\phi$",
#                          M2 = "Random intercepts and slopes",
#                          M1 = "Random intercepts",
#                          M4 = "Random intercepts",
#                          M5 = "Autoregressor $\\phi$ or random intercepts")) %>%
  select(Model, Type, elpd_diff, elpd_loo)

names(looc)[c(3,4)] <- c("$\\Delta\\widehat{elpd}$", "$\\widehat{elpd}$")

mc <- read_csv("../results/loo_results_consonants_mog.csv") %>% 
  filter(Model == "MoGpptsbigramintercepts - MoGpptsbigraminterceptsdelta") %>% 
  mutate_if(is.numeric, round, 0)  %>% 
  mutate_if(is.numeric, abs)

mc2 <- read_csv("../results/loo_results_consonants_mog.csv") %>% 
  filter(Model == "MoGpptsbigraminterceptsdelta2 - MoGpptsbigramintercepts") %>% 
  mutate_if(is.numeric, round, 0)
