library(tidyverse)
library(patchwork)

## read in plot objects
p_self_report <- read_rds("manuscript/figures/R_plot_objects/SPARK_language_diagnosis_ES-PGS_forest.rds")
p_reversion <- read_rds('manuscript/figures/R_plot_objects/SPARK_rare_reversion_diagnosis_forest.rds')
p_rev_hist <- read_rds('manuscript/figures/R_plot_objects/SPARK_rare_reversion_histograms_horizontal.rds')

## merge plots
plot_top <- (p_self_report / p_reversion) +
  plot_layout(heights = c(.5,.5, 1)) & 
  theme(plot.tag = element_text(face = 'bold', size = 24),
        axis.text = element_text(size = 18),
        # axis.text.x = element_text(angle = 18, hjust = 1),
        axis.title = element_text(size = 20),
        legend.text = element_text(size = 18),
        legend.title = element_text(size = 20),
        strip.text = element_text(size = 20))
plot_bottom <- (p_rev_hist) &
  theme(axis.text = element_text(size = 18),
        axis.title = element_text(size = 20),
        legend.text = element_text(size = 18),
        legend.title = element_text(size = 20),
        strip.text = element_text(size = 20))

## save individual plot components to put together in biorender
ggsave(plot_top, 
       filename = 'manuscript/figures/paper_figures/fig3_top.pdf', 
       device = 'pdf', 
       units = 'in', width = 14, height = 9,
       dpi = 300)
ggsave(plot_bottom, 
       filename = 'manuscript/figures/paper_figures/fig3_bottom.pdf', 
       device = 'pdf', 
       units = 'in', width = 14, height = 5,
       dpi = 300)       

ggsave(plot_top, 
       filename = 'manuscript/figures/paper_figures/fig3_top.png', 
       device = 'png', 
       units = 'in', width = 14, height = 9,
       dpi = 300)
ggsave(plot_bottom, 
       filename = 'manuscript/figures/paper_figures/fig3_bottom.png', 
       device = 'png', 
       units = 'in', width = 14, height = 5,
       dpi = 300)    