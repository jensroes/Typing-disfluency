set.seed(123)
sample_ppts <- sample(1:250,n_ppts)

d_beta %<>% mutate(id = as.character(id))
d_beta %>%
  filter(id != 0, id %in% sample_ppts) %>%
  mutate(Comp = recode(Comp, LF = "LF bigrams")) %>%
  ggplot(aes(x = M, xmin = lo, xmax =up, y = reorder(id, M), colour = Comp)) +
  geom_vline(xintercept = beta_sum[4], linetype = "dotted") +
  geom_vline(xintercept = beta_sum[1], linetype = "dotted") +
  geom_pointrangeh(position = posn.j, fatten = 1, alpha = .5) +
  scale_colour_manual(label, values = mycolours[c(7,4)]) +
  scale_x_continuous(breaks = seq(100, 1400, 100)) +
  labs(y = "Participant ID", x= "By-participant fluent-typing estimates [in msecs]") +
  theme(strip.text = element_text(hjust = 0),
        axis.text.y = element_text(size = 2),
        axis.title.y = element_text(size = 9)) +
  guides(colour = guide_legend(reverse = TRUE),
         fill = guide_legend(reverse = TRUE)) -> p_beta;p_beta

d_theta %>%
  group_by(Comp) %>%
  filter(id != 0, id %in% sample_ppts) %>%
  mutate(Comp = recode(Comp, LF = "LF bigrams")) %>%
  ggplot(aes(x = M, xmin = lo, xmax =up, y = reorder(id, M), colour = Comp)) +
  geom_vline(xintercept = theta_sum[5], linetype = "dotted") +
  geom_vline(xintercept = theta_sum[1], linetype = "dotted") +
  geom_pointrangeh(position = posn.j, fatten = 1, alpha = .5) +
  scale_x_continuous(limits = c(0,1), breaks = seq(0,1,.25)) +
  scale_colour_manual(label, values = mycolours[c(7,4)]) +
  labs(y = "Participant ID", x= "By-participant typing-disfluency probability") +
  theme(strip.text = element_text(hjust = 0),
        axis.text.y = element_text(size = 2)) -> p_theta;p_theta

legend <- get_legend(p_theta + theme(legend.position = "top", 
                                   legend.justification = "right",
                                   plot.margin = unit(c(0,0,0,0), "cm")))

p_ppts <- cowplot::plot_grid(NULL, legend,
                             p_beta + theme(legend.position = "none"), 
                             p_theta + theme(legend.position = "none"), 
                             nrow = 2, labels = c("","","A", "B"), 
                             rel_heights = c(1.5,10),
                             axis = "l", align = "b");p_ppts

left_join(d_beta,d_theta, by = c("Comp", "id")) %>%
  left_join(d_phi, by = c("Comp", "id")) %>%
  filter(Param.y != "phi") %>%
  select(-starts_with("Param"),-starts_with("p.")) %>%
  rename(prob = M.y,
         beta = M.x,
         phi = M) %>% select(id, Comp, beta, prob, phi) %>%
  pivot_longer(beta:phi, names_to = "Param", values_to = "value") %>%
  pivot_wider(names_from = Comp, values_from = value) -> d_cors

d_cors %>% filter(Param == "beta") %>%
  ggplot(aes(x = Consonants, y = LF)) +
  geom_point(alpha = .35, size = .35) +
  geom_smooth(method = "lm", se = F, color = "darkred") +
  scale_x_continuous(breaks = seq(100, 600, 100)) +
  scale_y_continuous(breaks = seq(100, 600, 100)) +
  labs(y = "LF bigrams", x = "", 
       subtitle = bquote("By-participant typing estimates [in msecs]")) +
  theme(plot.subtitle = element_text(size = 8),
        plot.margin = margin(r = 0, t = 5, b = 0, l = 5))-> p_cor_beta

d_cors %>% filter(Param == "phi") %>%
  ggplot(aes(x = Consonants, y = LF)) +
  geom_point(alpha = .35, size = .35) +
  geom_smooth(method = "lm", se = F, color = "darkred") +
  #  scale_x_continuous(breaks = seq(-.1, .1, .05)) +
  scale_y_continuous(breaks = seq(-.1, .1, .05)) +
  labs(y = "", x = "", 
       subtitle = bquote("By-participant autoregressor estimates")) +
  theme(plot.subtitle = element_text(size = 8),
        plot.margin = margin(r = 0, t = 5, b = 0, l = -5)) -> p_cor_phi;p_cor_phi

d_cors %>% filter(Param == "prob") %>%
  ggplot(aes(x = Consonants, y = LF)) +
  geom_point(alpha = .35, size = .35) +
  geom_smooth(method = "lm", se = F, color = "darkred") +
  scale_x_continuous(breaks = seq(0, 1, .25)) +
  scale_y_continuous(breaks = seq(0, .5, .1)) +
  labs(y = "", x = "", 
       subtitle = bquote("By-participant disfluency probability estimates")) +
  theme(plot.subtitle = element_text(size = 8),
        plot.margin = margin(r = 0, t = 5, b = 0, l = -5)) -> p_cor_theta;p_cor_theta


p_cors <- cowplot::plot_grid(p_cor_beta, p_cor_phi, p_cor_theta, 
                             labels = c("C", "D", "E"),
                             ncol = 3, axis = "l", align = "b")

p_cors <- ggdraw(add_sub(p_cors, "Consonants", 
                         vpadding=grid::unit(0,"lines"),
                         y=6, x=0.5, vjust=4.75, size = 9))


plots_post <- cowplot::plot_grid(p_ppts, p_cors, nrow = 2, 
                                 rel_heights = c(3,2.25), 
                                 axis = "l", align = "b");plots_post

