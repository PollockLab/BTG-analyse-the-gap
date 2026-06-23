# Script to assemble barplots for the Blitz the Gap effect figure

library(patchwork)
library(ggplot2)

A = readRDS("outputs/btg-effect-figs/barplot-btgeffect-nolabels.rds")
B = readRDS("outputs/btg-effect-figs/barplot-bioblitzeffect-nolabels.rds")
C = readRDS("outputs/btg-effect-figs/barplot-fundingeffect-nolabels.rds")

A = A + theme(plot.margin = margin(r = 0.1, t = 0, b = 0),
              legend.text = element_text(size = 16),
              axis.text.x = element_blank(),
              axis.text.y = element_text(size = 18),
              axis.title.y = element_text(size = 20,hjust = .5),
              axis.title.x = element_blank(),
              legend.position = "none")+
  coord_cartesian(ylim = c(0, 600000)) 
B = B + theme(axis.title.y = element_blank(),
              plot.margin = margin(l = 0.1, t = 0, b = 0),
              axis.text.x = element_blank(),
              axis.text.y = element_blank(),
              legend.text = element_text(size = 16),
              axis.title.x = element_blank(),
              legend.position = "none") +
  coord_cartesian(ylim = c(0, 600000)) 
C = C + theme(axis.title.y = element_blank(),
              axis.text.x = element_blank(),
              axis.text.y = element_blank(),
              plot.margin = margin(l = 0.1, t = 0, b = 0),
              legend.text = element_text(size = 16),
              axis.title.x = element_blank(),
              legend.position = "none") +
  coord_cartesian(ylim = c(0, 600000)) 

A + B + C + plot_annotation(tag_levels = "a")
ggsave("figures/btg-effect/barplot.png", width = 15, height = 6)


A2 = readRDS("outputs/btg-effect-figs/barplot-btgeffect-ratio.rds")
B2 = readRDS("outputs/btg-effect-figs/barplot-bioblitzeffect-ratio.rds")
C2 = readRDS("outputs/btg-effect-figs/barplot-fundingeffect-ratio.rds")

A2 = A2 + theme(plot.margin = margin(r = 0.1, t = 1, b = 0),
              legend.text = element_text(size = 16),
              axis.text.y = element_text(size = 18),
              axis.text.x = element_text(size = 14, angle = 45, hjust = 1, vjust = 1.1),
              axis.title.y = element_text(size = 20,hjust = .5),
              axis.title.x = element_blank(),
              legend.position = "none") +
  coord_cartesian(ylim = c(0,12))+
  scale_y_sqrt(labels = scales::label_number(scale =1, suffix = "x"))
(B2 = B2 + theme(axis.title.y = element_blank(),
              plot.margin = margin(l = 0.1, t = 1, b = 0),
              axis.text.x = element_text(size = 14, angle = 45, hjust = 1, vjust = 1.1),
              axis.text.y = element_blank(),
              legend.text = element_text(size = 16),
              axis.title.x = element_blank(),
              legend.position = "none") +
  coord_cartesian(ylim = c(0,12))+
  scale_y_sqrt())
C2 = C2 + theme(axis.title.y = element_blank(),
                axis.text.x = element_text(size = 14, angle = 45, hjust = 1, vjust = 1.1),
                axis.text.y = element_blank(),
              plot.margin = margin(l = 0.1, t = 1, b = 0),
              legend.text = element_text(size = 16),
              axis.title.x = element_blank(),
              legend.position = "none") +
  coord_cartesian(ylim = c(0,12)) +
  # wrap long text in the x axis
  scale_x_discrete(labels = function(x) str_wrap(x, width = 10))+
  scale_y_sqrt()

(A2 + B2 + C2) + plot_annotation(tag_levels = "a")
ggsave("figures/btg-effect/barplot-ratio.png", width = 13.5, height = 8)


(A + B + C)/(A2 + B2 + C2) + plot_annotation(tag_levels = "a")
ggsave("figures/btg-effect/barplot-allpanels.png", width = 8.04, height = 8.05)
