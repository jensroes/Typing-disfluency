library(tidyverse)
library(stringi)
library(xtable)

path <- "slides/sig-27-2020/r_output/"
files <- dir(path, pattern = "loo")

map(files, ~read_csv(paste0(path,.x)) %>%
      mutate(Comp = gsub("loo_results_|.csv", "", .x))) %>%
  bind_rows() %>%
  mutate_if(is.numeric, round, 0) %>%
  unite(col = "elpd_diff", elpd_diff, se_diff, sep = " (") %>%
  unite(col = "elpd_loo", elpd_loo, se_elpd_loo, sep = " (") %>%
  mutate(elpd_diff = paste0(elpd_diff, ")"),
         elpd_loo = paste0(elpd_loo, ")"),
         Model = gsub("ppts|bigram|intercepts", "", Model),
         Description = recode(Model, MoG = "Log-normal + disfluency component",
                                     LMM = "Log-normal model")) %>%
  pivot_wider(names_from = Comp, values_from = c(elpd_diff, elpd_loo)) %>%
  select(Model, Description, ends_with("consonants"), ends_with("lf")) -> looc;looc

stri_sub(looc[nchar(looc$elpd_loo_lf) >= 10,]$elpd_loo_lf, 4, 1) <- ","
stri_sub(looc[nchar(looc$elpd_loo_consonants) >= 10,]$elpd_loo_consonants, 4, 1) <- ","

looc
names(looc)[c(3:6)] <- c("$\\Delta\\widehat{elpd}$", "$\\widehat{elpd}$")

print(xtable(looc), include.rownames=FALSE, )
