library(tidyverse)
library(ggstance)

dir("results")
(files <- dir("results", pattern = "_LMM|_MoG.csv"))
(files <- files[!grepl("y_tilde", files)])

# Load posteriors
ps <- map(files, ~read_csv(paste0("results/", .x)) %>%
        mutate(Comp = gsub("_posterior_MoG.csv|_posterior_LMM.csv", "", .x),
               Model = gsub("LF_posterior_|consonants_posterior_|.csv", "", .x)) %>%
        gather(Param, value, -Comp, -Model)) %>%
  bind_rows() %>%
  filter(Param %in% c("beta", "delta", "prob"))

ps2 <- ps %>%
  group_by(Param, Model) %>%
  mutate(id = 1:n()) %>%
  pivot_wider(names_from = Param, values_from = value) %>%
  mutate(delta = ifelse(Model == "MoG", beta + delta, delta),
         Comp = recode(Comp, consonants = "Consonants",
                             LF = "LF bigrams")) %>%
  mutate_at(vars(delta, beta), exp) %>% select(-id) %>%
  pivot_longer(beta:prob, names_to = "Param", values_to = "value") %>%
  unite("Param", Model:Param) %>%
  drop_na()

ps_summary <- ps2 %>%
  group_by(Param, Comp) %>%
  summarise(M = mean(value),
            lo = quantile(value, probs = .025),
            up = quantile(value, probs = .975)) %>%
  separate(Param, into = c("Model", "Param"))

  #legend_title <- "Copy-task\ncomponent:"
legend_title <- ""

legend_mapping <- c("MoG_delta" = bquote("MoG: Disfluency ("*alpha~+~delta*")"),
                    "MoG_beta" = bquote("MoG: Fluent typing ("*alpha*")"),
                    "LMM_beta" = bquote("LMM: Typing ("*alpha*")"))
p_ikis <- ps_summary %>% filter(Param %in% c("delta", "beta")) %>%
  unite("Param", Model:Param) %>%
  ggplot(aes(x = M, xmin = lo, xmax = up, y = Comp, 
             colour = Param, 
             linetype = Param,
             shape = Param)) +
  geom_pointrangeh(position = position_dodge(-.75), fatten = 7) +
  scale_colour_manual(values = c( "black", "firebrick3", "firebrick4"), labels = legend_mapping) +
  scale_linetype_manual(values = c("dashed", "solid", "solid"), labels = legend_mapping) +
  scale_shape_manual(values = c(21,23, 24), labels = legend_mapping) +
  theme_bw(base_size = 16) +
  scale_x_continuous(limits = c(140, 600), breaks = seq(100, 600, 50)) +
  labs(x = "IKIs [in msecs] with 95% PIs", y = "",
       colour = legend_title, shape = legend_title, linetype = legend_title) +
  theme(legend.position = "top",
        legend.direction = "vertical",
        legend.justification = "right",
        legend.key.width = unit(1.75, "cm"),
        legend.key.height = unit(.75, "cm"),
        plot.margin = margin(-.6,.6,.1,.2, "cm"),
        panel.grid = element_blank(),
        axis.title.x = element_text(hjust = 0),
        axis.title.y = element_text(angle = 360),
        axis.ticks.y = element_blank(),
        axis.text.y = element_text(size = 17)); p_ikis

p_prob <- ps_summary %>% filter(Param %in% c("prob")) %>%
  ggplot(aes(x = M, xmin = lo, xmax = up, y = Comp)) +
  geom_pointrangeh(fatten = 7, color  = "firebrick4", shape = 24) + 
  scale_y_discrete(labels = c("prob" = bquote(atop("Disfluency", "probability ("*theta*")")))) +
  theme_bw(base_size = 16) +
  scale_x_continuous(limits = c(.25, .85), breaks = seq(0, 1, .1)) +
  labs(x = bquote(atop("Disfluency probability "*theta*" with 95% PIs")), y = "",
       colour = legend_title, shape = legend_title, linetype = legend_title) +
  theme(plot.margin = margin(0.1,.6,.1,.2, "cm"),
        panel.grid = element_blank(),
        legend.position = "none",
        axis.title.x = element_text(hjust = 0),
        axis.title.y = element_text(angle = 360),
        axis.ticks.y = element_blank(),
        axis.text.y = element_text(size = 17));p_prob

# extract the legend from one of the plots
#legend <- cowplot::get_legend(
  # create some space to the left of the legend
#  p_ikis + theme(legend.box.margin = margin(0, 0, 0, -20))
#)

path <- "slides/sig-27-2020/gfx/"
ggsave(paste0(path,"pop_beta.pdf"), plot = p_ikis, width = 8, height = 5.5, device = cairo_pdf)
ggsave(paste0(path,"pop_prob.pdf"), plot = p_prob, width = 8, height = 4, device = cairo_pdf)

  
  