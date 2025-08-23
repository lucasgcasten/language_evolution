library(tidyverse)
library(patchwork)

## read in plot objects
p_self_report <- read_rds("manuscript/figures/R_plot_objects/SPARK_language_diagnosis_ES-PGS_forest.rds")
p_reversion <- read_rds('manuscript/figures/R_plot_objects/SPARK_rare_reversion_diagnosis_forest.rds')
p_cpd <- read_rds('manuscript/figures/R_plot_objects/ABCD_cephalopelvic_disproportion_HAQER_ES-PGS.rds')
p_rev_hist <- read_rds('manuscript/figures/R_plot_objects/SPARK_rare_reversion_histograms_horizontal.rds')

## merge plots
plots <-(p_self_report | p_reversion) / (p_cpd | p_rev_hist) + 
  plot_annotation(tag_levels = "A") & 
  theme(plot.tag = element_text(face = 'bold', size = 24),
        axis.text = element_text(size = 18),
        # axis.text.x = element_text(angle = 18, hjust = 1),
        axis.title = element_text(size = 20),
        legend.text = element_text(size = 18),
        legend.title = element_text(size = 20))

ggsave(plots, 
       filename = 'manuscript/figures/paper_figures/fig3.pdf', 
       device = 'pdf', 
       units = 'in', width = 20, height = 14,
       dpi = 300)

ggsave(plots, 
       filename = 'manuscript/figures/paper_figures/fig3.png', 
       device = 'png', 
       units = 'in', width = 20, height = 14,
       dpi = 300)