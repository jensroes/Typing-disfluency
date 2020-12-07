library(tidyverse)
library(grid)

casc <- tibble(stage = c(rep("Conceptual\nencoding",3), 
                         rep("Linguistic\nencoding",2), 
                         rep("Orthographic\nprocessing",2), 
                         rep("Motor\nplanning",2)),
               activity = c(rep("Planning", 2), "Buffering", "Planning",
                          rep("Buffering", 1), 
                          "Spelling difficulty",
                          rep("Planning", 3)), 
               start = c(0, 6,  10, 6,  10, 10, 18, 18, 20),
               end =   c(6, 10, 18, 10, 18, 18, 20, 20, 22)) %>% ungroup() %>%
  mutate(item = c(rep(4,3),rep(3,2),rep(2,2), rep(1,2)),
         end = end - .05) %>%
  group_by(stage, activity) %>%
  mutate(chunk = 1:(n()),
         chunk = as.character(chunk - 1)) %>% ungroup() %>%
  mutate(chunk = ifelse(chunk == 0, "italic('n')", paste0("italic('n')~+~", chunk)),
         chunk = ifelse(start == 0, paste0("chunk~", chunk), chunk),
         chunk = ifelse(item == 2 & activity == "Planning", "italic('n')~+~1", chunk),
         chunk = ifelse(item == 3 & activity == "Buffering", "italic('n')~+~1~(buffering)", chunk),
         chunk = ifelse(item == 4 & activity == "Buffering", "italic('n')~+~2~(buffering)", chunk),
         chunk = ifelse(activity == "Spelling difficulty", "italic('n')~(Spelling~difficulty)", chunk )) %>%
  ungroup()

# Remove arrow before buffering
casc %<>%
  mutate(removearrow = rep("solid", n()), 
         chunk_lead = lead(chunk),
         removearrow = ifelse(chunk_lead == "Buffering", "blank", removearrow),
         removearrow = ifelse(is.na(removearrow), "solid", removearrow)) %>%
  select(-chunk_lead)

p_buffer <- ggplot(casc, aes(x = start, xend = end, 
                 y = reorder(stage,item),
                 yend = reorder(stage,item),
                 colour = activity)) + 
  geom_segment(size = 12, show.legend = F) +
  geom_segment(aes(x = end, xend = end, 
                   y = item -.25, 
                   yend = item - .65, linetype = removearrow),
               show.legend = F, colour = "grey60", size = .25,
               arrow = arrow(length = unit(0.2, "cm"),
                             type = "closed")) +
  geom_text(aes(x = (start + end) / 2 , label = chunk), 
            size = 4, colour = "white", parse = TRUE) +
  #ggtitle("b. Interrupted planning cascade") +
  scale_colour_viridis_d("", begin = .1, end = .9) +
  scale_x_continuous(limits = c(0,22), breaks = 20,labels = "writing\nonset") +
  scale_linetype_manual(values = rev(unique(casc$removearrow))) +
  labs(x = "Time", y = "Planning\nstage") + theme_bw() +
  theme(axis.text.x = element_text(face = "italic", size = 14),
        axis.text.y = element_text(size = 14),
        axis.title.x = element_text(hjust = 0, size = 14),
        axis.title.y = element_text(size = 14, hjust = 1, angle = 360),
        axis.ticks = element_blank(),
        panel.grid = element_blank(),
        legend.position = "top",
        legend.justification = "right"); p_buffer

ggsave("slides/sig-27-2020/gfx/cascade_buffer.pdf", width = 8.25, height = 4.5)


