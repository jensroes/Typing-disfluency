ps %>% select(-contains("_s")) %>%
  mutate(delta = beta + delta) %>%
  mutate_at(vars(beta, delta), exp) %>%
  mutate(delta = delta - beta) %>%
  gather(Param, value, -Comp) %>%
  group_by(Param, Comp) %>%
  summarise(M = mean(value),
            lo = quantile(value, probs = .025),
            up = quantile(value, probs = .975)) %>%
  ungroup() %>%
  mutate_at(vars(M, lo, up), round, 2) %>%
  mutate(M = ifelse(Param %in% c("beta", "delta"), round(M,0), M),
         lo = ifelse(Param %in% c("beta", "delta"), round(lo,0), lo),
         up = ifelse(Param %in% c("beta", "delta"), round(up,0), up)) %>%
  filter(Param != "sigma_diff") %>%
  mutate(var_comp = ifelse(grepl("sigma|tau",Param), TRUE, FALSE),
         Param = gsub("\\[1]", "", Param),
         description = "",
         description = ifelse(Param == "beta", "Fluent typing", description),
         description = ifelse(Param == "delta", "Disfluency slowdown", description),
         description = ifelse(Param == "phi", "Autoregressor", description),
         description = ifelse(Param == "prob", "Disfluency probability", description),
         description = ifelse(Param == "sigma", "Residual error", description),
#         description = ifelse(Param == "sigma_diff", "Difference between mixture components", description),
         description = ifelse(Param == "sigma_u", "Between-participants", description),
         description = ifelse(Param == "sigma_w", "Between-bigrams", description),
         description = ifelse(Param == "sigmap_e", "Disfluency", description),
         description = ifelse(Param == "sigma_e", "Fluent typing", description),
         description = ifelse(Param == "tau_phi", "Autoregressor", description),
         Param = recode(Param, "sigmap_e" = "sigma_{e'}", "tau_phi" = "eta", "prob" = "theta"), #"sigma_diff" = "sigma_{diff}"
         Param = ifelse(var_comp == TRUE, paste0(Param, "^2"), Param),
         Param = paste0("\\", Param)) -> param_table

#param_table %>% as.data.frame()
