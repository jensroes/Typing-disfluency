library(magrittr)
library(tidyverse)
library(ggstance)
library(ggforce)
library(ggthemes)
library(ggExtra)
library(cowplot)
library(plotmm)
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

d_example <- d %>% mutate(logIKI = log(IKI)) 

ppt1 <- d_example %>% filter(subj == 1)
out <- mixtools::normalmixEM(ppt1$logIKI, k = 2)

grid <- with(ppt1, seq(min(IKI), max(IKI), length = 100))
normaldens <- plyr::ddply(ppt1, "subj", function(df) {
  data.frame( 
    IKI = grid,
    logIKI = log(grid),
    density = dnorm(grid, mean(ppt1$IKI), sd(ppt1$IKI)),
    logdensity = dnorm(log(grid), mean(log(ppt1$IKI)), sd(log(ppt1$IKI)))
  )
})

x <- out$x
x <- data.frame(x)

ggplot(ppt1, aes(x = logIKI)) +
  geom_line(aes(y = ..density.., linetype = 'observed', colour = 'observed'), stat = 'density', size = .75) +
  geom_line(aes(y = logdensity, linetype = 'normal', colour = 'normal'), data = normaldens, size = .75) +
  stat_function(geom = "line", fun = plotmm::plot_mix_comps_normal,
                         args = list(mu = out$mu[1], sigma = out$sigma[1], lam = out$lambda[1]),
                         aes(colour =  "mixture", linetype = "mixture"), size = .75) +
  stat_function(geom = "line", fun = plotmm::plot_mix_comps_normal,
                         args = list(mu = out$mu[2], sigma = out$sigma[2], lam = out$lambda[2]),
                aes(colour =  "mixture", linetype = "mixture"), size = .75) +
  scale_linetype_manual("Density", values = c("dashed", "dashed", "solid")) +
  scale_colour_manual("Density", values = c("#69b3a2", "purple", "black")) + 
  labs(y = "", x = "IKIs [in log(msecs)]", caption = "Participant id: 1") +
  theme(strip.text = element_text(hjust = 0),
        axis.text.y = element_blank(),
        axis.title = element_text(hjust = 0),
        legend.key.width = unit(.75,"cm")) +
  guides(color = guide_legend(reverse = T),
         linetype = guide_legend(reverse = T))

ggsave("slides/sig-27-2020/gfx/mixture.pdf", width = 4, height = 5)

ggplot(ppt1, aes(x = logIKI)) +
  geom_line(aes(y = ..density.., linetype = 'observed', colour = 'observed'), stat = 'density', size = .75) +
  geom_line(aes(y = logdensity, linetype = 'normal', colour = 'normal'), data = normaldens, size = .75) +
  scale_linetype_manual("Density", values = c("dashed", "solid")) +
  scale_colour_manual("Density", values = c( "purple", "black")) + 
  labs(y = "", x = "IKIs [in log(msecs)]", caption = "Participant id: 1") +
  theme(strip.text = element_text(hjust = 0),
        axis.text.y = element_blank(),
        axis.title = element_text(hjust = 0),
        legend.key.width = unit(.75,"cm")) +
  guides(color = guide_legend(reverse = T),
         linetype = guide_legend(reverse = T))

ggsave("slides/sig-27-2020/gfx/lognormal2.pdf", width = 4, height = 5)

ggplot(ppt1, aes(x = logIKI)) +
  geom_line(aes(y = ..density.., linetype = 'observed', colour = 'observed'), stat = 'density', size = .75) +
  scale_linetype_manual("Density", values = c("solid")) +
  scale_colour_manual("Density", values = c("black")) + 
  labs(y = "", x = "IKIs [in log(msecs)]", caption = "Participant id: 1") +
  theme(strip.text = element_text(hjust = 0),
        axis.text.y = element_blank(),
        axis.title = element_text(hjust = 0),
        legend.key.width = unit(.75,"cm")) +
  guides(color = guide_legend(reverse = T),
         linetype = guide_legend(reverse = T))

ggsave("slides/sig-27-2020/gfx/observed.pdf", width = 4, height = 5)

ggplot(ppt1, aes(x = logIKI)) +
  geom_line(aes(y = ..density.., linetype = 'observed', colour = 'observed'), stat = 'density', size = .75) +
#  geom_line(aes(y = logdensity, linetype = 'normal', colour = 'normal'), data = normaldens, size = .75) +
  stat_function(geom = "line", fun = plotmm::plot_mix_comps_normal,
                args = list(mu = out$mu[1], sigma = out$sigma[1], lam = out$lambda[1]),
                aes(colour =  "mixture", linetype = "mixture"), size = .75) +
  stat_function(geom = "line", fun = plotmm::plot_mix_comps_normal,
                args = list(mu = out$mu[2], sigma = out$sigma[2], lam = out$lambda[2]),
                aes(colour =  "mixture", linetype = "mixture"), size = .75) +
  scale_linetype_manual("Density", values = c("dashed", "solid")) +
  scale_colour_manual("Density", values = c("#69b3a2", "black")) + 
  labs(y = "", x = "IKIs [in log(msecs)]", caption = "Participant id: 1") +
  theme(strip.text = element_text(hjust = 0),
        axis.text.y = element_blank(),
        axis.title = element_text(hjust = 0),
        legend.key.width = unit(.75,"cm")) +
  guides(color = guide_legend(reverse = T),
         linetype = guide_legend(reverse = T))

ggsave("slides/sig-27-2020/gfx/mixture2.pdf", width = 4, height = 5)

