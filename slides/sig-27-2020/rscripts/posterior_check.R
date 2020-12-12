library(loo)
library(tidyverse)
library(magrittr)
source("functions/functions.R")
source("functions/get_data.R")

# Load df
path <- "data/"
#d <- get_data(path = path) %>% filter(component == "Consonants") %>% select(-component)
d <- get_data(path = path) %>% filter(component %in% c("Consonants", "LF"), rep == 1) %>%
  select(subj, bg, IKI, component) %>%
  mutate(component = recode(component, LF = "LF bigrams task",
                                       Consonants = "Consonants task"),
         type = bquote("y[observed]")) %>%
  rename(Comp = component,
         y = IKI)

# Predicted data
y_pred <- read_csv("results/y_tilde_LMM_MoG.csv") %>%
  mutate(Comp = paste0(Comp, " task"))

type_labels <- c(bquote("LMM:"~tilde(italic("y"))), bquote("MoG:"~tilde(italic("y"))), bquote(italic("y")["obs"]))

# Posterior predictive checks
y_pred %<>% mutate(type = paste(Model, type, sep = "_"))

y_colours <- c("#69b3a2", "firebrick4", "black")

ggplot(y_pred, aes(x = y)) +
  geom_density(data = d, aes(x = y, group = NULL, colour = type), size = .75) +
  geom_density(alpha = .1, size = .2, aes(group = interaction(type,iteration), colour = type)) +
  facet_wrap( ~ Comp) +
  scale_x_continuous(limits = c(min(d$y), 1500)) +
  scale_colour_manual(values = y_colours, labels = parse(text = type_labels)) +
  theme_bw(base_size = 16) + labs(x = "IKIs [in msecs]", colour = "", caption = "Sample of 50 iterations") +
  theme(axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        axis.title.x = element_text(hjust = 0),
        strip.text = element_text(hjust = 0),
        strip.background = element_rect(colour = "white", fill = "white"),
        panel.grid = element_blank(),
        legend.position = "top",
        legend.justification = "right")

ggsave("slides/sig-27-2020/gfx/ppc.pdf", width = 8, height = 5)

ggplot(y_pred, aes(x = y)) +
  geom_density(data = d, aes(x = y, group = NULL, colour = type), size = .75) +
#  geom_density(alpha = .1, aes(group = interaction(type,iteration), colour = type)) +
  facet_wrap( ~ Comp) +
  scale_x_continuous(limits = c(min(d$y), 1500)) +
  scale_colour_manual(values = y_colours[3], labels = parse(text = type_labels[3])) +
  theme_bw(base_size = 16) + labs(x = "IKIs [in msecs]", colour = "", caption = "Sample of 50 iterations") +
  theme(axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        axis.title.x = element_text(hjust = 0),
        strip.text = element_text(hjust = 0),
        strip.background = element_rect(colour = "white", fill = "white"),
        panel.grid = element_blank(),
        legend.position = "top",
        legend.justification = "right")

ggsave("slides/sig-27-2020/gfx/ppc2.pdf", width = 8, height = 5)

ggplot(y_pred[y_pred$Model == "LMM",], aes(x = y)) +
  geom_density(data = d, aes(x = y, group = NULL, colour = type), size = .75) +
  geom_density(alpha = .1, size = .3, aes(group = interaction(type,iteration), colour = type)) +
  facet_wrap( ~ Comp) +
  scale_x_continuous(limits = c(min(d$y), 1500)) +
  scale_colour_manual(values = y_colours[c(1,3)], labels = parse(text = type_labels[c(1,3)])) +
  theme_bw(base_size = 16) + labs(x = "IKIs [in msecs]", colour = "", caption = "Sample of 50 iterations") +
  theme(axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        axis.title.x = element_text(hjust = 0),
        strip.text = element_text(hjust = 0),
        strip.background = element_rect(colour = "white", fill = "white"),
        panel.grid = element_blank(),
        legend.position = "top",
        legend.justification = "right")

ggsave("slides/sig-27-2020/gfx/ppc3.pdf", width = 8, height = 5)

d_cons_subjs <- d %>% filter(Comp == unique(Comp)[1]) %>% pull(subj)
d_lf_subjs <- d %>% filter(Comp == unique(Comp)[2]) %>% pull(subj)

y_pred_subj <- y_pred %>%
  group_by(Model, iteration) %>%
  mutate(subj = "",
         subj = ifelse(Comp == unique(Comp)[1], d_cons_subjs, subj),
         subj = ifelse(Comp == unique(Comp)[2], d_lf_subjs, subj)) %>%
  ungroup()

y_pred_sample <- y_pred_subj %>% filter(iteration %in% sample(paste0("V", 1:50), 50), subj %in% c(1,sample(1:250, 3))) %>%
  mutate(subj = paste0("Ppt id: ", subj))

d_plot <- d %>% mutate(subj = paste0("Ppt id: ", subj)) %>%
  filter(subj %in% y_pred_sample$subj)

ggplot(y_pred_sample, aes(x = y)) +
  geom_density(data = d_plot, 
               aes(x = y, group = NULL, colour = type), size = .75) +
  geom_density(alpha = .1, size = .3, aes(group = interaction(type,iteration), colour = type)) +
  facet_grid(Comp ~ subj, labeller = label_value) +
  scale_x_continuous(limits = c(min(d$y), 1500)) +
  scale_colour_manual(values = y_colours, labels = parse(text = type_labels)) +
  theme_bw(base_size = 16) + labs(x = "IKIs [in msecs]", colour = "", caption = "Sample of 50 iterations") +
  theme(axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        axis.title.x = element_text(hjust = 0),
        strip.text = element_text(hjust = 0),
        strip.background = element_rect(colour = "white", fill = "white"),
        panel.grid = element_blank(),
        legend.position = "top",
        legend.justification = "right")

ggsave("slides/sig-27-2020/gfx/ppc_ppt.pdf", width = 8, height = 6)

ggplot(y_pred_sample, aes(x = y)) +
  geom_density(data = d_plot, 
               aes(x = y, group = NULL, colour = type), size = .75) +
#  geom_density(alpha = .1, size = .3, aes(group = interaction(type,iteration), colour = type)) +
  facet_grid(Comp ~ subj, labeller = label_value) +
  scale_x_continuous(limits = c(min(d$y), 1500)) +
  scale_colour_manual(values = y_colours[3], labels = parse(text = type_labels[3])) +
  theme_bw(base_size = 16) + labs(x = "IKIs [in msecs]", colour = "", caption = "Sample of 50 iterations") +
  theme(axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        axis.title.x = element_text(hjust = 0),
        strip.text = element_text(hjust = 0),
        strip.background = element_rect(colour = "white", fill = "white"),
        panel.grid = element_blank(),
        legend.position = "top",
        legend.justification = "right")

ggsave("slides/sig-27-2020/gfx/ppc_ppt1.pdf", width = 8, height = 6)

ggplot(y_pred_sample[y_pred_sample$Model == "LMM",], aes(x = y)) +
  geom_density(data = d_plot, 
               aes(x = y, group = NULL, colour = type), size = .75) +
  geom_density(alpha = .1, size = .1, aes(group = interaction(type,iteration), colour = type)) +
  facet_grid(Comp ~ subj, labeller = label_value) +
  scale_x_continuous(limits = c(min(d$y), 1500)) +
  scale_colour_manual(values = y_colours[c(1,3)], labels = parse(text = type_labels[c(1,3)])) +
  theme_bw(base_size = 16) + labs(x = "IKIs [in msecs]", colour = "", caption = "Sample of 50 iterations") +
  theme(axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        axis.title.x = element_text(hjust = 0),
        strip.text = element_text(hjust = 0),
        strip.background = element_rect(colour = "white", fill = "white"),
        panel.grid = element_blank(),
        legend.position = "top",
        legend.justification = "right")

ggsave("slides/sig-27-2020/gfx/ppc_ppt2.pdf", width = 8, height = 6)

ggplot(y_pred_sample[y_pred_sample$Model == "MoG",], aes(x = y)) +
  geom_density(data = d_plot, 
               aes(x = y, group = NULL, colour = type), size = .75) +
  geom_density(alpha = .1, size = .1, aes(group = interaction(type,iteration), colour = type)) +
  facet_grid(Comp ~ subj, labeller = label_value) +
  scale_x_continuous(limits = c(min(d$y), 1500)) +
  scale_colour_manual(values = y_colours[c(2,3)], labels = parse(text = type_labels[c(2,3)])) +
  theme_bw(base_size = 16) + labs(x = "IKIs [in msecs]", colour = "", caption = "Sample of 50 iterations") +
  theme(axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        axis.title.x = element_text(hjust = 0),
        strip.text = element_text(hjust = 0),
        strip.background = element_rect(colour = "white", fill = "white"),
        panel.grid = element_blank(),
        legend.position = "top",
        legend.justification = "right")

ggsave("slides/sig-27-2020/gfx/ppc_ppt3.pdf", width = 8, height = 6)


#names(m)
#pars <- c("beta",  "theta", "delta", "sigma") #"sigmap_e", "sigma_e", "sigma_diff"
#print(m, pars, probs = c(.025, .975))
#traceplot(m, pars)
#ll <- extract_log_lik(m, merge_chains = F) 
#r_eff <- relative_eff(exp(ll)) 
#loo(ll, r_eff = r_eff, cores =2)

#Theta <- extract(m, 'theta')
#Theta <- unlist(Theta, use.names=FALSE)

