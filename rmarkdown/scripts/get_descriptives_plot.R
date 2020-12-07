d %<>% select(subj, bg, bigram, IKI, component)

d %>% group_by(subj, component) %>%
  mutate(bigram = 1:n()) %>%
  arrange(subj, component, bigram) %>%
  ungroup() %>%
  select(subj, bigram, bg, component) %>%
  pivot_wider(names_from = bigram, values_from = bg) %>% 
  unite("sequence", `1`:`20`, na.rm = T, remove = F, sep = " ") %>%  
  ungroup() %>%
  mutate(sequence = trimws(sequence)) %>%
  group_by(sequence) %>%
  mutate(n = n()) %>%
  mutate(p = n/500*100) -> d_wide

d %>% left_join(d_wide, by = c("subj", "component")) %>%
  select(-n) %>% 
  rename(p_seq = p) -> d2


lf <- c("een", "chaotische", "cowboy")
lf_bgs <- c("ee", "en", 
            "ch", "ha", "ao", "ot", "ti", "is", "sc", "ch", "he", 
            "co", "ow", "wb", "bo", "oy")

cons <- c("tjxgfl", "pgkfkq", "dtdrgt", "npwdvf")

cons_bgs <- c("tj", "jx", "xg", "gf", "fl", 
              "pg","gk", "kf", "fk",  "kq",
              "dt", "td", "dr", "rg", "gt", 
              "np", "pw", "wd", "dv", "vf")

length(cons_bgs)
length(lf_bgs)

d2 %>% ungroup() %>% 
  group_by(subj, component) %>%
  mutate(bigram = 1:n()) %>%
  mutate(max = max(bigram)) %>%
  ungroup() %>%
  select(subj, bigram, IKI, component, max) %>%
  filter((component == "Consonants" & bigram <= length(cons_bgs)) | 
           (component == "LF" & bigram <= length(lf_bgs) )) %>% 
  group_by(bigram, component) %>%
  summarise(Mode = dmode(IKI),
            Mean = mean(IKI),
            Median = median(IKI),
            lo = sd(IKI)/sqrt(n()),
            up = sd(IKI)/sqrt(n())) %>% 
  ungroup() %>%
  arrange(component, bigram) -> d_summary

d_summary$bgs <- c(cons_bgs, lf_bgs)
d_summary$block <- c(rep(1:4, each = 5), c(rep(1,2), rep(2, 9), rep(3,5)))

d_summary %>% 
  pivot_longer(Mode:Median, names_to = "CT", values_to = "values") %>%
  mutate(lo = values - lo,
         up = values + up) %>% 
  mutate(lo = ifelse(lo < 0, 0, lo)) %>% 
  mutate(chunk = ifelse(component == "Consonants", cons[block],
                        ifelse(component == "LF", lf[block], "test"))) %>%
  mutate(chunk = factor(chunk, levels = c(cons,lf), ordered = T)) -> d_summary_long


d2 %>% ungroup() %>% 
  group_by(subj, component) %>%
  mutate(bigram = 1:n()) %>%
  mutate(max = max(bigram)) %>%
  ungroup() %>%
  select(subj, bigram, IKI, component, max) %>%
  filter((component == "Consonants" & bigram <= length(cons_bgs)) | 
           (component == "LF" & bigram <= length(lf_bgs) )) %>% 
  left_join(d_summary[,c("bigram", "bgs", "block", "component")], by = c("bigram", "component")) %>%
  mutate(chunk = ifelse(component == "Consonants", cons[block],
                        ifelse(component == "LF", lf[block], "test"))) %>%
  mutate(chunk = factor(chunk, levels = c(cons,lf), ordered = T)) -> plot_df


plot_df %>%
  filter(component == "Consonants") %>%
  mutate(bgs = factor(bgs, levels = unique(bgs)[unique(bigram)], ordered = T)) %>% 
  ggplot(aes(y = IKI, x = bgs, group = factor(subj))) +
  facet_grid(~chunk, scales = "free_x", space = "free_x") +
  scale_y_log10(limits = c(50, 10000)) +
  geom_point(size = .001, show.legend = F,
             alpha = .25,
             colour = mycolours[7],
             position = position_dodge(.1)) +
  geom_line(linetype = "solid", show.legend = F, 
            colour = mycolours[7],
            size = .1, alpha = .15,  
            position = position_dodge(.1)) +
  labs(y = "", x = "") +
  theme(strip.text = element_text(hjust=0)) -> p_cons;p_cons

d_summary_long %>%
  filter(component == "Consonants") %>%
  mutate(bgs = factor(bgs, levels = unique(bgs)[unique(bigram)], ordered = T)) %>% 
  ggplot(aes(y = values, x = bgs, ymin = lo, ymax = up, linetype = CT, group = CT)) +
  facet_grid(~chunk, scales = "free_x", space = "free_x") +
  geom_pointinterval(position = position_dodge(width = .75), fatten_point = 1, 
                      colour = mycolours[7],
                      show.legend = F) +
  geom_line(position = position_dodge(.75), 
            colour = mycolours[7],
            size = .25,
            show.legend = T) +
  labs(x = "Bigram (in order)", y = "") +
  scale_y_log10(breaks = c(200, 300, 500, 800)) +
  scale_linetype_manual(name = "", values = Ms) +
  theme(strip.text = element_blank()) -> p_cons_sum;p_cons_sum


plot_df %>%
  filter(component == "LF") %>%
  ggplot(aes(y = IKI, x = factor(bigram), group = factor(subj))) +
  facet_grid(~chunk, scales = "free_x", space = "free_x") +
  scale_y_log10(limits = c(20,2000)) +
  scale_x_discrete(breaks = 1:16, labels = lf_bgs) +
  geom_point(size = .001, show.legend = F, 
             alpha = .25,
             colour = mycolours[4],
             position = position_dodge(.1)) +
  geom_line(linetype = "solid", 
            colour = mycolours[4],
            show.legend = F, 
            size = .1, alpha = .15,  
            position = position_dodge(.1)) +
  labs(y = "", x = "") +
  theme(strip.text = element_text(hjust=0)) -> p_lf;p_lf


d_summary_long %>%
  filter(component == "LF") %>%
  ggplot(aes(y = values, x = factor(bigram), ymin = lo, ymax = up, linetype = CT, group = CT)) +
  facet_grid(~chunk, scales = "free_x", space = "free_x") +
  scale_y_log10(breaks = c(100, 175, 300)) +
  scale_x_discrete(breaks = 1:16, labels = lf_bgs) +
  geom_pointinterval(position = position_dodge(width = .75), 
                      colour = mycolours[4],
                      fatten_point = 1, show.legend = F) +
  geom_line(position = position_dodge(.75), 
            colour = mycolours[4],
            size = .25,
            show.legend = T) +
  labs(x = "Bigram (in order)", y = "") +
  scale_linetype_manual(name = "", values = Ms) +
  theme(strip.text = element_blank(),
        legend.key.width = unit(.75, "cm")) -> p_lf_sum;p_lf_sum


d.raw <- d2 %>%
  group_by(subj, component) %>%
  mutate(bigram = 1:n()) %>%
  mutate(max = max(bigram)) %>%
  ungroup() %>%
  select(subj, bigram, IKI, component, max) %>%
  filter((component == "Consonants" & bigram <= length(cons_bgs)) | 
           (component == "LF" & bigram <= length(lf_bgs) )) %>% 
  ungroup() %>% 
  mutate(component = recode(component, LF = "LF bigrams")) %>%
  group_by(component) %>%
  summarise(mean = mean(IKI),
            median = median(IKI), 
            mode = dmode(IKI))

d2 %>% 
  group_by(subj, component) %>%
  mutate(bigram = 1:n()) %>%
  mutate(max = max(bigram)) %>%
  ungroup() %>%
  select(subj, bigram, IKI, component, max) %>%
  filter((component == "Consonants" & bigram <= length(cons_bgs)) | 
           (component == "LF" & bigram <= length(lf_bgs) )) %>% 
  mutate(component = recode(component, LF = "LF bigrams")) %>%
  ggplot(aes( x = IKI, group = component, fill = component, colour = component)) +
  geom_segment(data=d.raw,  
               aes(x=mean, y=-Inf, yend = Inf, xend = mean, linetype="Mean"), size = .15) +
  geom_segment(data=d.raw,  
               aes(x=median, y=-Inf, yend = Inf, xend = median, linetype="Median"), size = .15) +
  geom_segment(data=d.raw, 
               aes(x=mode, y=-Inf, yend = Inf, xend = mode, linetype="Mode"), size = .15) +
  scale_fill_manual(label, values = mycolours[c(7,4)]) +
  scale_colour_manual(label, values = mycolours[c(7,4)]) +
  scale_x_log10(limits = c(20, 10000)) +
  geom_density(size = .25, alpha = .35) +
  labs(x = "log-scaled IKIs [in msecs]", y = "Density") +
  theme(legend.key.width = unit(.75, "cm"),
        legend.key = element_blank()) +
  scale_linetype_manual(name = "", values = Ms) +
  guides(colour = guide_legend(reverse = TRUE),
         fill = guide_legend(reverse = TRUE)) -> p.log;p.log


plotA <- plot_grid(p_cons + theme(legend.position="none",
                 axis.text.x = element_blank() ,
                 plot.margin = unit(c(0,.1,-.2,-.1), "cm")),
                   p_cons_sum + theme(legend.position="none",
                   axis.title.x = element_blank(),
                   plot.margin = unit(c(-.1,.1,.2,-.1), "cm")),
  nrow = 2, axis = "l", align = 'v');plotA

y.grob <- textGrob("Consonants", gp=gpar(fontsize=9), rot=270)

plotA <- arrangeGrob(plotA, right = y.grob)

plotB <- plot_grid(p_lf + theme(legend.position="none",
                                  axis.text.x = element_blank() ,
                                  plot.margin = unit(c(0,.1,-.2,.05), "cm")),
                   p_lf_sum + theme(legend.position="none",
                                      axis.title.x = element_blank(),
                                      plot.margin = unit(c(-.1,.1,.2,.05), "cm")),
                   nrow = 2, axis = "l", align = 'v');plotB

y.grob <- textGrob("LF bigrams", gp=gpar(fontsize=9), rot=270)

plotB <- arrangeGrob(plotB, right = y.grob)

plot_top <- plot_grid(plotB, plotA, nrow = 2, axis = "l")

y.grob <- textGrob("log-scaled IKIs [in msecs]", gp=gpar(fontsize=10), rot=90)

plot_top <- arrangeGrob(plot_top, left = y.grob)

legend <- get_legend(p.log + theme(legend.position = "bottom", 
                                   legend.justification = "right",
                                   plot.margin = unit(c(0,0,0,0), "cm")))
plots <- plot_grid(
  plot_top, 
  p.log + theme(legend.position="none", 
                plot.margin = unit(c(.3,.6,.1,.7), "cm")),
  ncol = 2, rel_widths = c(2,1.25), 
  labels = c("A", "B")
)

plots_all <- plot_grid(plots, legend, nrow = 2, rel_heights = c(3,.5))
