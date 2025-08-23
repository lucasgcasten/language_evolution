library(tidyverse)
library(patchwork)

## read in plot objects
p_es_pgs_forest <- read_rds('manuscript/figures/R_plot_objects/EpiSLI_factor_ES-PGS_forest.rds') + xlab(NULL)
p_es_pgs_forest_all <- read_rds('manuscript/figures/R_plot_objects/EpiSLI_factor_comparison_ES-PGS_forest.rds')
p_f1 <- read_rds('manuscript/figures/R_plot_objects/EpiSLI_factor1_HAQER-CP-PGS.rds')
p_f1_bg <- read_rds('manuscript/figures/R_plot_objects/EpiSLI_factor1_background_HAQER_ES-PGS.rds')
p_f1_matched <- read_rds('manuscript/figures/R_plot_objects/EpiSLI_factor1_matched_HAQER_ES-PGS.rds')

## merge plots (A-G)
layout = 
"
AAABBBBBB
DDDEEEFFF
"

plots <- p_es_pgs_forest  + p_es_pgs_forest_all +
  p_f1 + p_f1_bg + p_f1_matched +
  plot_layout(design = layout) + plot_annotation(tag_levels = "A") & 
  theme(plot.tag = element_text(face = 'bold', size = 24),
        axis.text = element_text(size = 18),
        # axis.text.x = element_text(angle = 18, hjust = 1),
        axis.title = element_text(size = 20),
        legend.text = element_text(size = 18),
        legend.title = element_text(size = 20))

ggsave(plots, 
       filename = 'manuscript/figures/paper_figures/fig2.pdf', 
       device = 'pdf', 
       units = 'in', width = 20, height = 14,
       dpi = 300)

ggsave(plots, 
       filename = 'manuscript/figures/paper_figures/fig2.png', 
       device = 'png', 
       units = 'in', width = 20, height = 14,
       dpi = 300)

## manually merge above figures in BioRender (figure alignment is wonky in R)

# ## merge plots (A-G)
# plots <- (p_es_pgs_forest | p_f1 | p_f3) / (p_sp_rev_dist | p_sp_rev_forest) + 
#   plot_annotation(tag_levels = "a") & theme(plot.tag = element_text(face = 'bold'))

# ggsave(plots, 
#        filename = 'manuscript/figures/paper_figures/fig2.pdf', 
#        device = 'pdf', 
#        units = 'in', width = 24, height = 14,
#        dpi = 300)
