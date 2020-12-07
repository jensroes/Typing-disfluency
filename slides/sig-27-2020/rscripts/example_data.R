library(magrittr)
library(tidyverse)
library(ggstance)
library(ggforce)
library(ggthemes)
library(ggExtra)
library(cowplot)

source("functions/functions.R")
source("functions/get_data.R")

theme_set(theme_few(base_size = 14) + theme(axis.ticks = element_blank(),
                                           legend.position = "top",
                                           legend.justification = "right"))

mycolours = c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
Ms <- c("Median" = "dotted", "Mean" = "longdash", "Mode" = "solid")
Mcol <- c("Median" = "darkolivegreen4", "Mean" = "turquoise4", "Mode" = "darkred")

posn.j <- position_dodge(width = 1)
label = "Copy-task\ncomponent: "

# Load df
path <- "data/"
#d <- get_data(path = path) %>% filter(component == "Consonants") %>% select(-component)
d <- get_data(path = path) %>% filter(component %in% c("Consonants", "LF"), rep == 1) %>%
  select(subj, bg, bigram, IKI, sex, age, component)

lf_bgs <- c("ee", "en", 
            "ch", "ha", "ao", "ot", "ti", "is", "sc", "ch", "he", 
            "co", "ow", "wb", "bo", "oy")#;length(lf_bgs)

d %<>% filter(component == "LF") %>% 
  mutate(subj = as.character(subj)) %>%
  group_by(bg, subj) %>%
  filter(bg %in% lf_bgs) %>%
  group_by(subj) %>%
  mutate(bigram = 1:n(),
         n2 = n()) %>%
  filter(n2 == 16) %>%
  ungroup()

set.seed(123)
# Create the plot for a few ppts (first to all and then select)
d_example <- d %>% filter(subj %in% sample(x = unique(subj), size = 6)) %>%
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
         subj = paste0("Participant id: ", subj, " (M=", round(mean,0), ", SD=", round(sd,0), ")"))

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



# Empirical and normal
normalscale <- ggplot(d_example, aes(x = IKI))  + 
  geom_line(aes(y = ..density.., linetype = 'Empirical', colour = 'Empirical'), stat = 'density', size = .75) +
  geom_line(aes(y = density, linetype = 'Normal', colour = 'Normal'), data = normaldens, size = .75) +
  scale_linetype_manual("Density", values = c("solid", "dashed")) +
  scale_colour_manual("Density", values = c("#69b3a2", "purple")) + #
  geom_vline(data = filter(d_summary, subj == unique(subj)[1]), aes(xintercept = mean), 
             linetype = "dotted", colour = "grey30") +
  labs(y = "Density", x = "IKIs [in msecs]") +
  facet_wrap(~subj, nrow = 2) +
  theme(strip.text = element_text(hjust = 0),
        axis.text.y = element_blank(),
        legend.key.width = unit(1,"cm")); normalscale

ggsave("slides/sig-27-2020/gfx/normalemp.pdf", width = 8.25, height = 5)


# Empirical
normalscale <- ggplot(d_example, aes(x = IKI))  + 
  geom_line(aes(y = ..density.., linetype = 'Empirical', colour = 'Empirical'), stat = 'density', size = .75) +
  scale_linetype_manual("Density", values = c("solid")) +
  scale_colour_manual("Density", values = c("#69b3a2")) + #
  geom_vline(data = filter(d_summary, subj == unique(subj)[1]), aes(xintercept = mean), 
             linetype = "dotted", colour = "grey30") +
  labs(y = "Density", x = "IKIs [in msecs]") +
  facet_wrap(~subj, nrow = 2) +
  theme(strip.text = element_text(hjust = 0),
        axis.text.y = element_blank(),
        legend.key.width = unit(1,"cm")); normalscale

ggsave("slides/sig-27-2020/gfx/emp.pdf", width = 8.25, height = 5)

# Normal empirical and normal on log scale

logscale <- ggplot(d_example, aes(x = logIKI)) +
  geom_line(aes(y = ..density.., linetype = 'Empirical', colour = 'Empirical'), stat = 'density', size = .75) +
  geom_line(aes(y = logdensity, linetype = 'Normal', colour = 'Normal'), data = normaldens, size = .75) +
  scale_linetype_manual("Density", values = c("solid", "dashed")) +
  scale_colour_manual("Density", values = c("#69b3a2", "purple")) + #
  geom_vline(data = filter(d_summary, subj == unique(subj)[1]), aes(xintercept = logmean), 
             linetype = "dotted", colour = "grey30") +
  labs(y = "Density", x = "log-scaled IKIs [in log(msecs)]") +
  facet_wrap(~subj, nrow = 2) +
  theme(strip.text = element_text(hjust = 0),
        axis.text.y = element_blank(),
        legend.key.width = unit(1,"cm"));logscale 

ggsave("slides/sig-27-2020/gfx/lognormal.pdf", width = 8.25, height = 5)
  

logscale <- ggplot(d_example, aes(x = logIKI)) +
  geom_line(aes(y = ..density.., linetype = 'Empirical', colour = 'Empirical'), stat = 'density', size = .75) +
  scale_linetype_manual("Density", values = c("solid")) +
  scale_colour_manual("Density", values = c("#69b3a2")) + #
  geom_vline(data = filter(d_summary, subj == unique(subj)[1]), aes(xintercept = logmean), 
             linetype = "dotted", colour = "grey30") +
  labs(y = "Density", x = "IKIs [in log(msecs)]") +
  facet_wrap(~subj, nrow = 2) +
  theme(strip.text = element_text(hjust = 0),
        axis.text.y = element_blank(),
        legend.key.width = unit(1,"cm"));logscale 

ggsave("slides/sig-27-2020/gfx/lognormal2.pdf", width = 8.25, height = 5)
