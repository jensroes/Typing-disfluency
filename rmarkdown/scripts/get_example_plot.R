

# Create the plot for a few ppts (first to all and then select)
d_example <- d %>% filter(component == "LF",
                          subj %in% c(74:75)) %>%
  select(-component) %>%
  mutate(logIKI = log(IKI)) %>%
  group_by(subj) %>%
  mutate(mean = round(mean(IKI),0),
         logmean = mean(logIKI),
         median = median(IKI),
         mode = dmode(IKI),
         sd = round(sd(IKI),0),
         logsd = sd(logIKI),
         logsd = sd(log(IKI)),
         n = n(),
         lognorm = dnorm(logIKI, logmean, logsd),
         norm = dnorm(IKI, mean, sd)) %>%
  ungroup() %>%
  mutate(maxnorm = max(norm),
         maxlognorm = max(lognorm),
         maxIKI = max(IKI),
         maxlogIKI = max(logIKI),
         subj = paste0("Participant id: ", subj))

d_summary <- d_example %>% group_by(subj) %>% select(mean:n, maxnorm:maxlogIKI) %>% unique()

grid <- with(d_example, seq(min(IKI), max(IKI), length = 100))
normaldens <- plyr::ddply(d_example, "subj", function(df) {
  data.frame( 
    IKI = grid,
    logIKI = log(grid),
    density = dnorm(grid, mean(df$IKI), sd(df$IKI)),
    logdensity = dnorm(log(grid), mean(log(df$IKI)), sd(log(df$IKI)))
  )
})

normalscale <- ggplot(d_example, aes(x = IKI))  + 
  geom_line(aes(y = ..density.., linetype = 'Empirical'), stat = 'density') +
  geom_line(aes(y = density, linetype = 'Normal'), data = normaldens) +
  scale_linetype_manual("Density", values = c("dashed", "solid")) +
  geom_vline(data = filter(d_summary, subj == unique(subj)[1]), aes(xintercept = mean), 
             linetype = "dotted", colour = "grey30") +
  geom_vline(data = filter(d_summary, subj == unique(subj)[2]), aes(xintercept = mean), 
             linetype = "dotted", colour = "grey30") +
  labs(y = "Density", x = "IKIs [in msecs]") +
  facet_wrap(~subj, nrow = 2) +
  theme(strip.text = element_text(hjust = 0)); normalscale


logscale <- ggplot(d_example, aes(x = logIKI)) +
  geom_line(aes(y = ..density.., linetype = 'Empirical'), stat = 'density') +
  geom_line(aes(y = logdensity, linetype = 'Normal'), data = normaldens) +
  scale_linetype_manual("Density", values = c("dashed", "solid")) +
  geom_vline(data = filter(d_summary, subj == unique(subj)[1]), aes(xintercept = logmean), 
             linetype = "dotted", colour = "grey30") +
  geom_vline(data = filter(d_summary, subj == unique(subj)[2]), aes(xintercept = logmean), 
             linetype = "dotted", colour = "grey30") +
  facet_wrap(~subj, nrow = 2) +
  labs(y = "", x = "log-scaled IKIs [in msecs]") +
  theme(strip.text = element_text(hjust = 0)) 

legend <- get_legend(logscale)

plot <- cowplot::plot_grid(NULL, legend,
                           normalscale + theme(legend.position = "none"),
                           logscale + theme(legend.position = "none"),
                           rel_heights = c(.5,4,4),
                           axis = 'l',
                           align = 'b',
                           labels = c("","", "A", "B"),
                           ncol=2, nrow = 2)
