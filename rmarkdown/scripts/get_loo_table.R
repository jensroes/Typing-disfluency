loo_cons <- read_csv("../results/loo_results_consonants.csv") %>%
  mutate(ctc = "Consonants")

loo_lf <- read_csv("../results/loo_results_LF.csv") %>%
  mutate(ctc = "LF bigrams")

bind_rows(loo_cons, loo_lf) %>%
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
  select(ctc, Model, Type, elpd_diff, elpd_loo) -> looc

stringi::stri_sub(looc$elpd_loo, 4, 1) <- ","
stringi::stri_sub(looc[nchar(looc$elpd_diff) >= 10,]$elpd_diff, 3, 1) <- ","

names(looc)[c(4,5)] <- c("$\\Delta\\widehat{elpd}$", "$\\widehat{elpd}$")

mc_cons <- read_csv("../results/loo_results_consonants_mog.csv")

mc_lf <- read_csv("../results/loo_results_LF_mog.csv")


