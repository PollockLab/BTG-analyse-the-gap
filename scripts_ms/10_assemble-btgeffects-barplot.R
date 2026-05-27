# Script to assemble barplots for the Blitz the Gap effect figure

library(patchwork)
library(ggplot2)

A = readRDS("outputs/btg-effect-figs/barplot-btgeffect.rds")
B = readRDS("outputs/btg-effect-figs/barplot-bioblitzeffect.rds")
C = readRDS("outputs/btg-effect-figs/barplot-fundingeffect.rds")

A = A + theme(plot.margin = margin(r = 0),
              legend.text = element_text(size = 16),
              axis.text.x = element_text(size = 12),
              axis.title.y = element_text(size = 20),
              axis.title.x = element_text(size = 20))+
  coord_cartesian(ylim = c(0, 600000)) 
B = B + theme(axis.title.y = element_blank(),
              plot.margin = margin(l = 0),
              axis.text.x = element_text(size = 12),
              legend.text = element_text(size = 16),
              axis.title.x = element_text(size = 20)) +
  coord_cartesian(ylim = c(0, 600000)) 
C = C + theme(axis.title.y = element_blank(),
              axis.text.x = element_text(size = 12),
              plot.margin = margin(l = 0),
              legend.text = element_text(size = 16),
              axis.title.x = element_text(size = 20)) +
  coord_cartesian(ylim = c(0, 600000))

A + B + C + plot_annotation(tag_levels = "a")
ggsave("figures/btg-effect/barplot.png", width = 14.5, height = 6.8)
