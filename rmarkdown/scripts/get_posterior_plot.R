ggplot(d_beta, aes(y = M, ymin = lo, ymax =up, x =id2, 
                   color = Param, 
                   shape = Param,
                   size = Param)) +
  geom_hline(yintercept = beta_sum[1], linetype = "dotted") +
  geom_errorbar(alpha = .5, width = 0, size = .4, show.legend = F) +
  geom_point(show.legend = F) +
  scale_color_manual(values = c("darkred", "black")) +
  scale_shape_manual(values = c(17, 19)) +
  scale_size_manual(values = c(1.5, .1)) +
  coord_flip() + theme_few(base_size = 10) +
  labs(x = "Participants", 
       y = bquote(atop("IKIs"~hat(beta)~"[in msecs]"))) +
  theme(axis.text.y = element_blank(),
        axis.ticks = element_blank()) -> p_beta

ggplot(d_theta, 
       aes(y = M, ymin = lo, ymax =up, x =id2, 
           color = Param, shape = Param, size = Param)) +
  geom_hline(yintercept = theta_sum[1], linetype = "dotted") +
  geom_errorbar(alpha = .5, width = 0, size = .4, show.legend = F) +
  geom_point(show.legend = F) +
  scale_color_manual(values = c("darkred", "black")) +
  scale_shape_manual(values = c(17, 19)) +
  scale_size_manual(values = c(1.5, .1)) +
  coord_flip() + theme_few(base_size = 10) +
  labs(x = "Participants", 
       y = bquote(atop("Disfluency probability"~hat(theta)))) +
  theme(axis.text.y = element_blank(),
        axis.ticks = element_blank()) -> p_theta


ggplot(pss, aes(x = M_beta, y = M_prob)) +
  geom_point(size = .5, alpha = .8) +
  geom_errorbar(aes(ymin = lo_prob, ymax = up_prob), show.legend = F, size =.05, alpha =.15) +
  geom_errorbarh(aes(xmin = lo_beta, xmax = up_beta), show.legend = F, size = .05, alpha =.15) +
  geom_vline(xintercept=beta_sum[1], linetype="dashed", alpha = 0.4, colour = "darkred") +
  geom_hline(yintercept=theta_sum[1], linetype="dashed", alpha = 0.4, colour = "darkred") +
  geom_errorbar(inherit.aes = F, aes(ymin = theta_sum[2], ymax = theta_sum[3], x = beta_sum[1]), 
                colour = "darkred", size =.1, width = 0) +
  geom_errorbarh(inherit.aes = F, aes(xmin = beta_sum[2], xmax = beta_sum[3], y = theta_sum[1]), 
                 colour = "darkred", size =.1, height = 0) +
  geom_point(inherit.aes =F, aes(x = beta_sum[1], y = theta_sum[1]), shape=17, colour = "darkred")+
  theme_few(base_size = 10) +
  labs(y = bquote("Disfluency probability"~hat(theta)),
       x = bquote("IKIs"~hat(beta)~"[in msecs]")) +
  theme(axis.ticks = element_blank()) -> p_scatter

ps %>% select(beta, delta) %>%
  transmute(delta = exp(beta+delta) - exp(beta)) %>%
  ggplot(aes(x = delta)) +
  geom_density(fill = "darkred", alpha = .45, colour = "darkred") +
  geom_segment(aes(x = quantile(delta, probs = .025), 
                   xend = quantile(delta, probs = .975), y=0, yend=0), size=1) +
  geom_point(aes(x = mean(delta), y =0), size = 2, color = "darkred", shape = 17) +
  labs(x = bquote("Disfluency slowdown"~hat(delta)~"[in msecs]"), 
       y = "Posterior density") +
  theme_few(base_size = 10) +
  theme(axis.ticks = element_blank()) -> p_param

plot_post <- plot_grid(p_beta, p_theta, p_scatter, p_param, align = "vh", labels = c("A", "B", "C", "D"), ncol=2, nrow = 2)