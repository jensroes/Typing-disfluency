library(loo)
library(tidyverse)
library(magrittr)

path <- "stanout/consonants/"
#path <- "stanout/"

(files <- dir(path, pattern = ".rda"))


## Check model fit
m <- readRDS(paste0(path, files[2] ))
traceplot(m, pars = "beta")

param <- c("beta",  "phi", "theta_corr", "theta_incorr", "delta_corr", "delta_incorr", 
           "sigma", "sigma_diff_corr", "sigma_diff_incorr", 
           "sigma_e", "sigmap_e_corr", "sigmap_e_incorr") 

plot(m, pars=param)
plot(m, show_density = TRUE, ci_level = 0.5, fill_color = "purple")
plot(m, plotfun = "hist", pars = "theta_corr", include = FALSE)
plot(m, plotfun = "trace", pars = param, inc_warmup = FALSE)
plot(m, plotfun = "rhat") + ggtitle(bquote(hat(italic("R"))))
pairs(m, pars = c("beta", "phi", "delta_corr", "delta_incorr"))
