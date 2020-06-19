ytilde <- read_csv("../results/consonants_posterior_MoG_ytilde.csv")

y <- d %>% filter(bigram != 1)

ytilde %<>%
  mutate(samples = 1:n()) %>%
  filter(samples %in% sample(samples, 500)) %>%
  pivot_longer(cols = starts_with("y_tilde"), names_to = "id", values_to = "y_tilde") 

ggplot() +
  geom_density(data = ytilde, aes(x = y_tilde, group = samples, linetype = "y_tilde"), 
               alpha = .25, colour = "darkred") + 
  geom_density(data= y, aes(x = IKI, linetype = "y"), fill = "darkred", alpha = .25) +
  scale_x_log10() +
  scale_linetype_manual(values = c("y_tilde" = "dotted", "y" = "solid"))
