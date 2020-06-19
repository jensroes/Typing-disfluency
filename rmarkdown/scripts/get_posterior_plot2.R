d_beta %>%
  filter(id != 0) %>%
  group_by(Comp) %>%
  arrange(Comp, M) %>%
  mutate(id = 1:n()) %>%
  ungroup() %>%
  mutate(Comp = recode(Comp, LF = "LF bigrams")) %>%
  ggplot(aes(x = M, xmin = lo, xmax =up, y = reorder(id, M), 
         color = Comp, 
         linetype = Comp,
         shape = Comp)) +
  geom_pointrangeh(position = posn.j, fatten = .25, alpha = .5) +
  scale_color_manual(label, values = mycolours[c(7,4)]) +
  scale_shape_manual(label, values = c(1, 2)) +
  scale_linetype_manual(label, values = c("dashed", "dotted")) +
  scale_x_continuous(limits = c(50, 800), breaks = seq(100, 800, 100)) +
  labs(y = "", 
       x = bquote(atop("Fluent typing IKIs"~hat(beta)~"[in msecs]"))) +
  theme(axis.text.y = element_blank()) -> p_beta

d_theta %>%
  filter(id != 0) %>%
  group_by(Comp) %>%
  mutate(id = as.numeric(factor(id))) %>%
  arrange(Comp, M) %>%
  mutate(id = 1:n()) %>%
  ungroup() %>%
  mutate(Comp = recode(Comp, LF = "LF bigrams")) %>%
  ggplot(aes(x = M, xmin = lo, xmax =up, y = reorder(id, M),
             color = Comp, 
             linetype = Comp,
             shape = Comp)) +
  geom_pointrangeh(position = posn.j, fatten = .25, alpha = .5) +
  scale_color_manual(label, values = mycolours[c(7,4)]) +
  scale_shape_manual(label, values = c(1, 2)) +
  scale_linetype_manual(label, values = c("dashed", "dotted")) +
  scale_x_continuous(limits = c(0,1), breaks = seq(0,1,.25)) +
  labs(y = "", 
       x = bquote(atop("Disfluency probability"~hat(theta)))) +
  theme(axis.text.y = element_blank()) -> p_theta


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
  scale_x_continuous(breaks = seq(50, 400, 50)) +
  labs(x = bquote("Disfluency slowdown"~hat(delta)~"[in msecs]"), 
       y = "Posterior probability density") -> p_param


#pss %>%
#  mutate(Comp = recode(Comp, LF = "LF bigrams")) %>%
#  pivot_wider(-ends_with("prob"), names_from = Comp, values_from = c(M_beta, lo_beta, up_beta)) %>%
#  ggplot(aes(x = M_beta_Consonants, y = `M_beta_LF bigrams`, 
#             ymin = lo_prob, ymax = up_prob,
#             xmin = lo_beta, xmax = up_beta
#)) +
#  geom_point(size = .2, alpha = .5) +
#  geom_smooth(method = "lm")

pss %>%
  mutate(Comp = recode(Comp, LF = "LF bigrams")) %>%
  ggplot(aes(x = M_beta, y = M_prob, 
             ymin = lo_prob, ymax = up_prob,
             xmin = lo_beta, xmax = up_beta,
             group = Comp)) +
  geom_point(aes(shape = Comp, colour = Comp), size = .2, alpha = .5) +
  stat_density_2d(geom = "polygon", aes(alpha = ..level.., fill = Comp), 
                  bins = 10, show.legend = F) +
  scale_y_continuous(limits = c(0,1)) +
  scale_x_continuous(breaks = seq(100, 700, 100), limits = c(75, 550)) +
  scale_color_manual(label, values = mycolours[c(7,4)]) +
  scale_fill_manual(label, values = mycolours[c(7,4)]) +
  scale_shape_manual(label, values = c(1, 2)) +
  geom_point(inherit.aes =F, aes(x = beta_sum[1], y = theta_sum[1]), shape=16, size = 1.5, colour = "black") +
  geom_point(inherit.aes =F, aes(x = beta_sum[4], y = theta_sum[5]), shape=17, size = 1.5, colour = "black") +
  geom_errorbar(inherit.aes = F, aes(ymin = theta_sum[2], ymax = theta_sum[3], x = beta_sum[1]), 
                colour = "black", size =.25, width = 0) +
  geom_errorbarh(inherit.aes = F, aes(xmin = beta_sum[2], xmax = beta_sum[3], y = theta_sum[1]), 
                 colour = "black", size =.25, height = 0) +
  geom_errorbar(inherit.aes = F, aes(ymin = theta_sum[6], ymax = theta_sum[7], x = beta_sum[4]), 
                colour = "black", size =.25, width = 0) +
  geom_errorbarh(inherit.aes = F, aes(xmin = beta_sum[5], xmax = beta_sum[6], y = theta_sum[5]), 
                 colour = "black", size =.25, height = 0) +
  labs(y = bquote("Disfluency probability"~hat(theta)),
       x = bquote("Fluent typing IKIs"~hat(beta)~"[in msecs]")) -> p_scatter;p_scatter

p_theta <- p_theta + guides(colour = guide_legend(reverse = TRUE),
                 shape = guide_legend(reverse = TRUE),
                 linetype =  guide_legend(reverse = TRUE))

legend <- get_legend(p_theta)

plot_post <- cowplot::plot_grid(NULL, legend,
                       p_beta + theme(legend.position = "none"), 
                       p_theta + theme(legend.position = "none"),
                       p_scatter + theme(legend.position = "none"), 
                       p_param + theme(legend.position = "none"),
                       rel_heights = c(.5,4,2.85),
                       axis = 'l',
                       align = 'b',
                       labels = c("","", "A", "B", "C", "D"),
                       ncol=2, nrow = 3)#;plot_post
