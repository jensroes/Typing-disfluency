Ms <- c("Median" = "dotted", "Mean" = "dashed", "Mode" = "longdash")
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
            SE = 2*sd(IKI)/sqrt(n())) %>% 
  ungroup() -> d_summary


cons <- c("tjxgfl", "pgkfkq", "dtdrgt", "npwdvf")
d_summary$bgs <- c("tj", "jx", "xg", "gf", "fl",
                   "pg","gk", "kf", "fk",  "kq", 
                   "dt", "td", "dr", "rg", "gt", 
                   "np", "pw", "wd", "dv", "vf")

d_summary %<>%
  mutate(bgs = factor(bgs, levels = unique(bgs)[bigram], ordered = T))

d2 %>% ungroup() %>% group_by(subj) %>% mutate(bigram = 1:n()) %>% 
  filter(bigram < 21) %>% ungroup() %>%
  left_join(d_summary[,c("bigram", "bgs")], by = "bigram") %>%
  mutate(bgs = factor(bgs, levels = unique(bgs)[unique(bigram)], ordered = T)) %>% 
  ggplot() +
  geom_point(aes(y = IKI, x = bgs, group = factor(subj)), size =.5, show.legend = F, alpha = .1) +
  geom_line(aes(y = IKI, x = bgs, group = factor(subj)), linetype = "dashed", show.legend = F, alpha = .05) +
  geom_pointrange(data = d_summary, aes(y = Mode, x = bigram, ymin = Mode - SE, ymax= Mode + SE, color = "Mode"),fatten = 3) +
  geom_pointrange(data = d_summary, aes(y = Mean, x = bigram, ymin = Mean - SE, ymax= Mean + SE, color = "Mean"), fatten = 3) +
  geom_pointrange(data = d_summary, aes(y = Median, x = bigram, ymin = Median - SE, ymax= Median + SE, color = "Median"), fatten = 3) +
  geom_line(data = d_summary, aes(y = Mode, x = bgs, group = 1, color = "Mode", linetype = "Mode")) +
  geom_line(data = d_summary, aes(y = Median, x = bgs, group = 1, color = "Median", linetype = "Median")) +
  geom_line(data = d_summary, aes(y = Mean, x = bgs, group = 1, color = "Mean", linetype = "Mean")) +
  geom_vline(data = filter(d_summary), aes(xintercept = c(5.5)), linetype = "dotted", size = .25) +
  geom_vline(data = filter(d_summary), aes(xintercept = c(10.5)), linetype = "dotted", size = .25) +
  geom_vline(data = filter(d_summary), aes(xintercept = c(15.5)), linetype = "dotted", size = .25) +
  theme_few(base_size = 10) + scale_y_log10() +
  labs(x = "Bigram (in order)", y = "log-scaled IKIs [in msecs]") +
  scale_linetype_manual(name = "", values = Ms) +
  scale_color_manual(name = "", values = Mcol) +
  theme(strip.text = element_text(size = 12, hjust=0),
        axis.ticks = element_blank()) + 
  annotate(geom="text", x=c(1.75,7.75,12.75,17.75), y=17500, label=cons, size = 4) -> p.big

d.raw <- d2 %>% ungroup() %>% 
  group_by(subj) %>% mutate(bigram = 1:n()) %>% 
  filter(bigram < 21) %>% ungroup() %>% summarise(mean = mean(IKI), 
                                                  median = median(IKI), 
                                                  mode = dmode(IKI))


d2 %>% ungroup() %>% group_by(subj) %>% 
  mutate(bigram = 1:n()) %>% filter(bigram < 21) %>% 
  ungroup() %>%
  ggplot(aes( x = IKI)) +
  stat_density(geom="line", size = .2, colour = "black", show.legend = F, position = "identity") +
  scale_x_log10() +
  theme_void(base_size = 10) +
  labs(x = "", y = "") +
  geom_segment(data=d.raw, inherit.aes = F, 
               aes(x=mean, y=0, yend = .75, xend = mean, linetype="Mean", color = "Mean")) +
  geom_segment(data=d.raw, inherit.aes = F, 
               aes(x=median, y=0, yend = 1.1, xend = median, linetype="Median", color = "Median")) +
  geom_segment(data=d.raw, inherit.aes = F, 
               aes(x=mode, y=0, yend = 1.25, xend = mode, linetype="Mode", color = "Mode")) +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        legend.key.width = unit(1, "cm"),
        legend.justification = "top") +
  scale_linetype_manual(name = "", values = Ms) +
  scale_color_manual(name = "", values = Mcol) +
  coord_flip() -> p.log;p.log

pcol <- plot_grid(
  p.big + theme(legend.position="none", plot.margin = unit(c(.1, 0, .1, .1), "cm")),
  p.log + theme(legend.position="none", plot.margin = unit(c(.1, 0, .1, 0), "cm")),
  align = 'h', nrow = 1,
  rel_widths = c(4,.5)
)

legend <- get_legend(p.log + theme(legend.position = "top", legend.justification = "right"))

plot_all <- plot_grid(legend, pcol, ncol = 1, rel_heights = c(.1,1))
