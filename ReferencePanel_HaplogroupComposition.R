##################################################
##### Reference panel haplogroup composition #####
##################################################

library(tidyverse)
library(stringr)
library(ggplot2)

Haplogroups <- read.delim("ReferencePanelHaplogroups.txt", header=TRUE)

major_haplogroups <- Haplogroups %>%
  mutate(Major_HG = str_extract(Haplogroup, "^[A-Z]")) %>%
  count(Major_HG, sort = TRUE)

major_haplogroups

haplo_colors_major <- c(
  "A" = "#5E81AC",
  "B" = "#81A1C1",
  "C" = "#4C72B0",
  "D" = "#3B5B92",
  "E" = "#6B8E23",
  "F" = "#4F81BD",
  "G" = "#2C7FB8",
  "H" = "#9E0142",
  "I" = "#B2182B",
  "J" = "#D6604D",
  "K" = "#F4A582",
  "L" = "#8C510A",
  "M" = "#3288BD",
  "N" = "#5E4FA2",
  "O" = "#BDBDBD",
  "P" = "#7B3294",
  "Q" = "#8E63B0",
  "R" = "#C51B7D",
  "S" = "#998EC3",
  "T" = "#F46D43",
  "U" = "#D53E4F",
  "V" = "#FDAE61",
  "W" = "#FEE08B",
  "X" = "#8073AC",
  "Y" = "#B2ABD2",
  "Z" = "#4575B4"
)

plot <- ggplot(major_haplogroups,
       aes(x = Major_HG,
           y = n,
           fill = Major_HG)) +
  
  geom_bar(stat = "identity") +
  
  scale_fill_manual(values = haplo_colors_major) +
  
  theme_minimal(base_size = 14) +
  
  labs(
    x = "Major Haplogroup",
    y = "Number of Individuals"
  ) +
  
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(face = "bold"),
    legend.position = "none"
  )

ggsave("Reference panel haplogroup composition.png", plot = plot, width = 12, height = 6, dpi = 300)
