Ms <- c("Median" = "dotted", "Mean" = "longdash", "Mode" = "solid")
Mcol <- c("Median" = "darkolivegreen4", "Mean" = "turquoise4", "Mode" = "darkred")

d %<>% select(subj, bg, bigram, IKI)

d %>% group_by(subj) %>%
  mutate(bigram = 1:n()) %>%
  arrange(subj, bigram) %>%
  ungroup() %>%
  select(subj, bigram, bg) %>%
  pivot_wider(names_from = bigram, values_from = bg) %>% 
  unite("sequence", `1`:`20`, na.rm = T, remove = F, sep = " ") %>%  
  ungroup() %>%
  mutate(sequence = trimws(sequence)) %>%
  group_by(sequence) %>%
  mutate(n = n()) %>%
  mutate(p = n/500*100) -> d_wide

d %>% left_join(d_wide, by = "subj") %>% select(-n) %>% rename(p_seq = p) -> d2

d2 %>% ungroup() %>% group_by(subj) %>%
  mutate(bigram = 1:n()) %>%
  ungroup() %>%
  group_by(bigram) %>%
  mutate(N = n()) %>%
  filter(bigram < 21) %>% 
  summarise(Mode = dmode(IKI),
            Mean = mean(IKI),
            Median = median(IKI),
            lo = sd(IKI)/sqrt(n()),
            up = sd(IKI)/sqrt(n())) %>% 
  ungroup() -> d_summary


cons <- c("tjxgfl", "pgkfkq", "dtdrgt", "npwdvf")
d_summary$bgs <- c("tj", "jx", "xg", "gf", "fl", 
                   "pg","gk", "kf", "fk",  "kq",
                   "dt", "td", "dr", "rg", "gt", 
                   "np", "pw", "wd", "dv", "vf")

d_summary %<>%
  mutate(block = rep(1:4, each = 5)) %>%
  mutate(bgs = factor(bgs, levels = unique(bgs)[bigram], ordered = T))

d_summary %>%
  pivot_longer(Mode:Median, names_to = "CT", values_to = "values") %>%
  mutate(lo = values - lo,
         up = values + up) %>% 
  mutate(lo = ifelse(lo < 0, 0, lo)) %>% 
  mutate(consB = plyr::mapvalues(block, from = unique(block), to = cons)) %>%
  mutate(consB = factor(consB, levels = cons, ordered = T)) -> d_summary_long


d2 %>% ungroup() %>% group_by(subj) %>% mutate(bigram = 1:n()) %>% 
  filter(bigram < 21) %>% ungroup() %>%
  left_join(d_summary[,c("bigram", "bgs", "block")], by = "bigram") %>%
  select(-bg,-sequence,-`1`:-`p_seq`) %>%
  mutate(bgs = factor(bgs, levels = unique(bgs)[unique(bigram)], ordered = T)) %>% 
  mutate(consB = plyr::mapvalues(block, from = unique(block), to = cons)) %>%
  mutate(consB = factor(consB, levels = cons, ordered = T)) -> plot_df


ggplot(data = plot_df, aes(x = IKI, y = bgs, group = factor(subj))) +
  coord_flip() +
  facet_wrap(~consB, scales = "free_x", nrow = 1 ) +
  theme_few(base_size = 10) + 
  scale_x_log10(limits = c(50, 20000)) +
  geom_point(size = .1, show.legend = F, position = position_dodge(.1)) +
  geom_line(linetype = "solid", show.legend = F, size = .1, alpha = .15, orientation = "y", position = position_dodge(.1)) +
  labs(y = "", x = "") +
  scale_linetype_manual(name = "", values = Ms) +
  scale_color_colorblind("") +
  theme(strip.text = element_text(size = 12, hjust=0),
        axis.ticks = element_blank(),
        axis.text.x = element_blank(),
        legend.key.width = unit(1, "cm")) -> p.ppt;p.ppt

ggplot(data = plot_df) +
  coord_flip() +
  facet_wrap(~consB, scales = "free_x", nrow = 1 ) +
  theme_few(base_size = 10) + 
  scale_x_log10(limits = c(130, 1200), breaks = c(250, 500, 1000)) +
#  scale_x_log10(limits = c(50, 20000)) +
  geom_pointintervalh(data = d_summary_long, aes(x = values, y = bgs, xmin = lo, xmax = up, color = CT), 
                                position = position_dodgev(height = .75), fatten_point = 2.5, show.legend = F) +
  geom_line(data = d_summary_long, aes(x = values, y = bgs, group = CT, color = CT, linetype = CT), 
            position = position_dodge(.75), show.legend = T, orientation = "y") +
  labs(y = "Bigram (in order)", x = "") +
  scale_linetype_manual(name = "", values = Ms) +
  scale_colour_wsj(name = "") +
  theme(strip.text = element_blank(),
        axis.ticks = element_blank(),
        legend.position = "top",
        legend.justification = "right",
        legend.key.width = unit(1, "cm")) -> p.big


d.raw <- d2 %>% ungroup() %>% 
  group_by(subj) %>% mutate(bigram = 1:n()) %>% 
  filter(bigram < 21) %>% ungroup() %>% summarise(mean = mean(IKI), 
                                                  median = median(IKI), 
                                                  mode = dmode(IKI))


d2 %>% ungroup() %>% group_by(subj) %>% 
  mutate(bigram = 1:n()) %>% filter(bigram < 21) %>% 
  ungroup() %>%
  ggplot(aes( x = IKI)) +
  geom_segment(data=d.raw, inherit.aes = F, 
               aes(x=mean, y=-Inf, yend = Inf, xend = mean, linetype="Mean", color = "Mean"), size = .5) +
  geom_segment(data=d.raw, inherit.aes = F, 
               aes(x=median, y=-Inf, yend = Inf, xend = median, linetype="Median", color = "Median"), size = .5) +
  geom_segment(data=d.raw, inherit.aes = F, 
               aes(x=mode, y=-Inf, yend = Inf, xend = mode, linetype="Mode", color = "Mode"), size = .5) +
  stat_density(geom="line", size = .5, 
               colour = "black", 
               show.legend = F,
               position = "identity") +
  scale_x_log10() +
  theme_few(base_size = 10) +
  labs(x = "log-scaled IKIs [in msecs]", y = "Density") +
  facet_grid(~1) +
  theme(axis.ticks = element_blank(),
        legend.key.width = unit(1, "cm"),
        legend.key = element_blank(),
        axis.title.x = element_text(size = 10),
        axis.title.y = element_text(size = 10, margin=margin(0,28,0,0)),
        strip.text = element_text(colour = "transparent"),
        legend.justification = "top") +
  scale_linetype_manual(name = "", values = Ms) +
  scale_colour_wsj(name = "") -> p.log

plotA <- plot_grid(
  p.ppt + theme(legend.position="none", plot.margin = unit(c(0,.1,-.1,.1), "cm")),
  p.big + theme(legend.position="none", plot.margin = unit(c(-.2,.1,.1,.1), "cm")), 
  nrow = 2, axis = c("l"), align = 'v')


y.grob <- textGrob("log-scaled IKIs [in msecs]", gp=gpar(fontsize=10), rot=90)

plotAgrob <- arrangeGrob(plotA, left = y.grob)

pcol <- plot_grid(
  plotAgrob, 
  p.log + theme(legend.position="none", plot.margin = unit(c(-.1,.1,.1,.1), "cm")),
  nrow = 2,  rel_heights = c(2,1.15), 
  labels = c("A", "B")
)

legend <- get_legend(p.log + theme(legend.position = "bottom", 
                                   legend.justification = "right",
                                   plot.margin = unit(c(0,0,0,0), "cm")))

plot_all <- plot_grid(pcol,legend, ncol = 1, rel_heights = c(1,.1))

