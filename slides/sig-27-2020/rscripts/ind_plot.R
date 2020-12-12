library(tidyverse)
library(ggstance)
library(cowplot)

dir("results")
(files <- dir("results", pattern = "_MoG.csv"))
(files <- files[!grepl(pattern = "y_tilde", files)])

# Load posteriors
ps <- map(files, ~read_csv(paste0("results/", .x)) %>%
            mutate(Comp = gsub("_posterior_MoG.csv", "", .x))) %>%
  bind_rows() 

(files <- dir("results", pattern = "_LMM.csv"))

# Load posteriors
ps_lmm <- map(files, ~read_csv(paste0("results/", .x)) %>%
            mutate(Comp = gsub("_posterior_LMM.csv", "", .x))) %>%
  bind_rows() 

# Summarise posterior
ps_summary_lmm <- ps_lmm %>% 
  mutate_at(vars(contains("_s")), ~exp(beta + .x)) %>%
  select(contains("_s"), Comp) %>%
  pivot_longer(-Comp) %>%
  separate(name, into = c("Param", "id"), sep = "\\[") %>%
  mutate(id = gsub("\\]", "", id)) %>%
  group_by(Comp, id, Param) %>%
  summarise(M = mean(value),
            lo = quantile(value, probs = .025),
            up = quantile(value, probs = .975)) %>%
  pivot_wider(names_from = "Comp", values_from = c("M", "lo", "up")) %>%
  mutate(Model = "LMM")

# Summarise posterior
ps_summary_mog <- ps %>% select(contains("_s"), Comp) %>%
  mutate_at(vars(starts_with("beta")), exp) %>%
  pivot_longer(-Comp) %>%
  separate(name, into = c("Param", "id"), sep = "\\[") %>%
  mutate(id = gsub("\\]", "", id)) %>%
  group_by(Comp, id, Param) %>%
  summarise(M = mean(value),
            lo = quantile(value, probs = .025),
            up = quantile(value, probs = .975)) %>%
  pivot_wider(names_from = "Comp", values_from = c("M", "lo", "up")) %>%
  mutate(Model = "MoG")

ps_prob <- ps_summary_mog %>% filter(Param == "prob_s") %>% select(-Param)
ps_beta_mog <- ps_summary_mog %>% filter(Param == "beta_s") %>% select(-Param)
ps_beta <- ps_summary_lmm %>% filter(Param == "beta_s") %>% select(-Param) %>%
  bind_rows(ps_beta_mog)

ggplot(ps_beta, aes(x = M_consonants, y = M_LF, colour = Model, shape = Model)) +
  geom_point(size = 3, alpha = .35, show.legend = F) +
  geom_smooth(method = "lm",  se = F, fullrange = TRUE) +
  scale_colour_manual("",values = c("black", "firebrick3")) +
  scale_shape_manual("",values = c(21, 24)) +
  theme_bw(base_size = 18) +
  scale_x_continuous(limits = c(140, 600), breaks = seq(100, 600, 50)) +
  scale_y_continuous(limits = c(100, 375), breaks = seq(100, 600, 50)) +
  labs(y = "LF-bigrams task",
       x = "Consonants task") +
  theme(legend.position = "top",
        panel.grid = element_blank(),
        legend.justification = "right",
        legend.key.width = unit(1.75, "cm"),
        plot.margin = margin(0,.6,.1,.2, "cm"),
        axis.title = element_text(hjust = 0))

path <- "slides/sig-27-2020/gfx/"
ggsave(paste0(path,"beta_corr.pdf"), width = 8, height = 5, device = cairo_pdf)


ggplot(ps_prob, aes(x = M_consonants, y = M_LF)) +
  geom_point(size = 3, alpha = .35, shape = 24, colour = "firebrick3") +
  geom_smooth(method = "lm",  se = F, fullrange = TRUE, colour = "firebrick3") +
  theme_bw(base_size = 18) +
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, .2)) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, .2)) +
  labs(y = "LF-bigrams task",
       x = "Consonants task") +
  theme(legend.position = "top",
        panel.grid = element_blank(),
        legend.justification = "right",
        legend.key.width = unit(1.75, "cm"),
        plot.margin = margin(0,.6,.1,.2, "cm"),
        axis.title = element_text(hjust = 0))

path <- "slides/sig-27-2020/gfx/"
ggsave(paste0(path,"prob_corr.pdf"), width = 8, height = 4, device = cairo_pdf)

ps_prob %>% pivot_longer(cols = M_consonants:up_LF,
  names_to = c("sum", "Comp"),
  names_pattern = "(.*)_(.*)",
  values_to = "value") %>%
  pivot_wider(names_from = "sum", values_from = value) %>% select(-Model) %>%
  mutate(Comp = recode(Comp, LF = "LF bigrams", consonants = "Consonants")) %>%
  #filter(Comp == "LF bigrams") %>%
  ggplot(aes(y = M, ymin = lo, ymax = up, x = reorder(id, M), colour = Comp)) +
  geom_pointrange(alpha = .5, size = .5) +
  theme_bw(base_size = 18) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, .2)) +
  scale_colour_viridis_d(begin = .5, end = 0) +
  labs(y = bquote("Probability "*theta*" with 95% PIs"),
       x = "Participant id", colour = "Copy-task\ncomponent:") +
  theme(legend.position = "top",
        legend.justification = "right",
        legend.key.size = unit(40, "pt"),
        plot.margin = margin(0,.6,.1,.2, "cm"),
        panel.grid = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title = element_text(hjust = 0, size = 18))

path <- "slides/sig-27-2020/gfx/"
ggsave(paste0(path,"prob_ppts.pdf"), width = 8.5, height = 5.5, device = cairo_pdf)

ps_summary_mog %>% 
  pivot_longer(cols = M_consonants:up_LF,
               names_to = c("sum", "Comp"),
               names_pattern = "(.*)_(.*)",
               values_to = "value") %>%
  group_by(Comp, id, Param) %>%
  summarise(M = mean(value),
            lo = quantile(value, probs = .025),
            up = quantile(value, probs = .975)) %>%
  pivot_wider(names_from = "Param", values_from = c("M", "lo", "up")) %>%
  mutate(Comp = recode(Comp, LF = "LF bigrams", consonants = "Consonants")) %>%
  ggplot(aes(x = M_beta_s, y = M_prob_s, colour = Comp, shape = Comp)) +
  geom_point(size = 3, alpha = .35, show.legend = F) +
  geom_smooth(method = "lm",  se = F, fullrange = TRUE) +
  theme_bw(base_size = 18) +
  scale_colour_viridis_d(begin = .5, end = 0) +
  scale_shape_manual(values = c(21, 24)) +
  scale_x_continuous(limits = c(100, 600), breaks = seq(100, 600, 100)) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, .2)) +
  labs(y = bquote("Disfluency probability"~theta),
       x = bquote("Fluent keystroke transitions "*alpha),
       colour = "Copy-task\ncomponent:",
       shape = "Copy-task\ncomponent:") +
  theme(legend.position = "top",
        panel.grid = element_blank(),
        legend.justification = "right",
        legend.key.width = unit(1.75, "cm"),
        plot.margin = margin(0,.6,.1,.2, "cm"),
        axis.title = element_text(hjust = 0))

path <- "slides/sig-27-2020/gfx/"
ggsave(paste0(path,"betaprob_corr.pdf"), width = 8, height = 5, device = cairo_pdf)


ps_beta_mog %>% pivot_longer(cols = M_consonants:up_LF,
                         names_to = c("sum", "Comp"),
                         names_pattern = "(.*)_(.*)",
                         values_to = "value") %>%
  pivot_wider(names_from = "sum", values_from = value) %>% select(-Model) %>%
  mutate(Comp = recode(Comp, LF = "LF bigrams", consonants = "Consonants")) %>%
  ggplot(aes(y = M, ymin = lo, ymax = up, x = reorder(id, M), colour = Comp)) +
  geom_pointrange(alpha = .5, size = .5) +
  theme_bw(base_size = 18) +
  scale_y_continuous(breaks = seq(0, 800, 100)) +
  scale_colour_viridis_d(begin = .5, end = 0) +
  labs(y = bquote("Fluent IKIs "*beta*" with 95% PIs"),
       x = "Participant id", colour = "Copy-task\ncomponent:") +
  theme(legend.position = "top",
        legend.justification = "right",
        legend.key.size = unit(40, "pt"),
        plot.margin = margin(0,.6,.1,.2, "cm"),
        panel.grid = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title = element_text(hjust = 0, size = 18))

path <- "slides/sig-27-2020/gfx/"
ggsave(paste0(path,"betas_ppts.pdf"), width = 8.5, height = 5.5, device = cairo_pdf)


