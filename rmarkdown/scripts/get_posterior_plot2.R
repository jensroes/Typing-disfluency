set.seed(123)
sample_ppts <- sample(1:250,n_ppts)
d_beta %>%
    filter(Comp == "LF", id %in% sample_ppts) %>%
    mutate(Comp = recode(Comp, LF = "LF bigrams")) %>%
    ggplot(aes(x = M, xmin = lo, xmax =up, y = reorder(id, M))) +
    geom_vline(xintercept = beta_sum[4], linetype = "dotted") +
    geom_pointrangeh(position = posn.j, fatten = 1, alpha = .5) +
    scale_x_continuous(breaks = seq(100, 1400, 50)) +
    facet_grid(~Comp) +
    labs(y = "Participant ID", x= "") +
    theme(strip.text = element_text(hjust = 0),
          axis.text.y = element_text(size = 3),
          axis.title.y = element_text(size = 9),
          plot.margin = margin(r = 0, t = 0, b = 0, l = 5)) -> p_lf;p_lf

d_beta %>%
  filter(Comp == "Consonants", id %in% sample_ppts) %>%
  ggplot(aes(x = M, xmin = lo, xmax =up, y = reorder(id, M))) +
  geom_vline(xintercept = beta_sum[1], linetype = "dotted") +
  geom_pointrangeh(position = posn.j, fatten = 1, alpha = .5) +
  scale_x_continuous(breaks = seq(100, 1400, 200), limits = c(100, 1000) ) +
  facet_grid(~Comp) +
  labs(y = "", x= "") +
  theme(strip.text = element_text(hjust = 0),
        axis.text.y = element_text(size = 3),
        plot.margin = margin(r = 0, t = 0, b = 0, l = -5)) -> p_cons;p_cons

p_beta <- cowplot::plot_grid(p_lf, p_cons)
p_beta <- ggdraw(add_sub(p_beta, "By-participant fluent-typing estimates [in msecs]", 
               vpadding=grid::unit(0,"lines"),y=6, x=0.5, vjust=4.75, size = 9)); p_beta

d_theta %>%
  group_by(Comp) %>%
  filter(Comp == "LF", id != 0, id %in% sample_ppts) %>%
  mutate(Comp = recode(Comp, LF = "LF bigrams")) %>%
  ggplot(aes(x = M, xmin = lo, xmax =up, y = reorder(id, M))) +
  geom_vline(xintercept = theta_sum[5], linetype = "dotted") +
  geom_pointrangeh(position = posn.j, fatten = 1, alpha = .5) +
  scale_x_continuous(limits = c(0,1), breaks = seq(0,1,.25)) +
  facet_grid(~Comp) +
  labs(y = "", x= "") +
  theme(strip.text = element_text(hjust = 0),
        axis.text.y = element_text(size = 3),
        plot.margin = margin(r = 0, t = 0, b = 0, l = 0.25)) -> p_lf;p_lf

d_theta %>%
  group_by(Comp) %>%
  filter(Comp == "Consonants", id != 0, id %in% sample_ppts) %>%
  ggplot(aes(x = M, xmin = lo, xmax =up, y = reorder(id, M))) +
  geom_vline(xintercept = theta_sum[1], linetype = "dotted") +
  geom_pointrangeh(position = posn.j, fatten = 1, alpha = .5) +
  scale_x_continuous(limits = c(0,1), breaks = seq(0,1,.25)) +
  facet_grid(~Comp) +
  labs(y = "", x= "") +
  theme(strip.text = element_text(hjust = 0),
        axis.text.y = element_text(size = 3),
        plot.margin = margin(r = 5, t = 0, b = 0, l = -5)) -> p_cons;p_cons

p_theta <- cowplot::plot_grid(p_lf, p_cons)
p_theta <- ggdraw(add_sub(p_theta, "By-participant typing-disfluency probability", 
                         vpadding=grid::unit(0,"lines"),y=6, x=0.5, vjust=4.75, size = 9))


  


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

p_ppts <- cowplot::plot_grid(p_beta, p_theta, 
                   nrow = 1, labels = c("A", "B"), 
                   axis = "l", align = "b");p_ppts

p_cors <- cowplot::plot_grid(p_cor_beta, p_cor_phi, p_cor_theta, 
                             labels = c("C", "D", "E"),
                             ncol = 3, axis = "l", align = "b")

p_cors <- ggdraw(add_sub(p_cors, "Consonants", 
                          vpadding=grid::unit(0,"lines"),
                         y=6, x=0.5, vjust=4.75, size = 9))


plots_post <- cowplot::plot_grid(p_ppts, p_cors, nrow = 2, 
                                 rel_heights = c(3,2.25), 
                                 axis = "l", align = "b");plots_post


ps %>% select(beta, delta, Comp) %>%
  group_by(Comp) %>%
  transmute(delta = exp(beta+delta) - exp(beta)) %>%
  ungroup() %>%
  mutate(Comp = recode(Comp, LF = "LF bigrams")) %>%
  ggplot(aes(x = delta, fill = Comp, linetype = Comp)) +
  geom_histogram(aes(y = ..density..), alpha = .75, binwidth = 5, position = position_dodge()) +
  scale_color_manual(label, values = mycolours[c(7,4)]) +
  scale_fill_manual(label, values = mycolours[c(7,4)]) +
  scale_linetype_manual(label, values = c("dashed", "dotted")) +
  scale_x_continuous(breaks = seq(50, 650, 50), limits = c(50, 600)) +
  labs(x = "Slowdown magnitude for\ndisfluent typing [in msecs]", 
       y = "Posterior-probability\ndensity") -> p_param
