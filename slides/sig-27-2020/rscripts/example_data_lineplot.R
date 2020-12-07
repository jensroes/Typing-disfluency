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

lf <- c("een", "chaotische", "cowboy")
lf_bgs <- c("ee", "en", 
            "ch", "ha", "ao", "ot", "ti", "is", "sc", "ch", "he", 
            "co", "ow", "wb", "bo", "oy")#;length(lf_bgs)

d_example2 <- d %>% filter(component == "LF") %>% 
  mutate(subj = as.character(subj)) %>%
  group_by(bg, subj) %>%
  filter(bg %in% lf_bgs) %>%
  group_by(subj) %>%
  mutate(bigram = 1:n(),
         n2 = n()) %>%
  filter(n2 == 16) %>%
  ungroup()

#d_example2 %>% pull(subj) %>% unique()
#d_example2 %>% pull(bg) # %>% unique()

#d_example3 <- d_example2 %>% filter(subj %in% c(105, 240, 238))
set.seed(123)
d_example3 <- d_example2 %>% filter(subj %in% as.numeric(sample(x = unique(subj), size = 6)))

ggplot(d_example3, aes(y = IKI, x = bigram, 
                       linetype = subj, 
                       colour = subj,
                       group = subj, shape = subj)) +
  geom_point(position = position_dodge(.5), size = 2) +
  geom_line(position = position_dodge(.5)) +
  labs(shape = "Participant id", 
       linetype = "Participant id",
       colour = "Participant id",
       y = "IKIs [in msecs]", x = "Bigrams") +
  scale_shape_manual(values = 20:26) +
  scale_color_viridis_d() +
  scale_x_continuous(breaks = sort(unique(d_example3[d_example3$ subj == unique(d_example3$subj)[1],]$bigram)), 
                     labels = d_example3[d_example3$subj == unique(d_example3$subj)[1],]$bg) +
  theme(legend.key.width =  unit(1, "cm"))

ggsave("slides/sig-27-2020/gfx/trial.pdf", width = 8.25, height = 5)
