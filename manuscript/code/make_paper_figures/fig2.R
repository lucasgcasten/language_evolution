library(tidyverse)
library(patchwork)

## read in plot objects
p_es_pgs_forest <- read_rds('manuscript/figures/R_plot_objects/EpiSLI_factor_ES-PGS_forest.rds')
p_f1 <- read_rds('manuscript/figures/R_plot_objects/EpiSLI_factor1_HAQER-CP-PGS.rds')
p_f2 <- read_rds('manuscript/figures/R_plot_objects/EpiSLI_factor2_HAQER-CP-PGS.rds')
p_f3 <- read_rds('manuscript/figures/R_plot_objects/EpiSLI_factor3_HAQER-CP-PGS.rds')
# p_gwas <- read_rds('manuscript/figures/R_plot_objects/HAQER_cogPerf_GWAS_zcsore_reversion_distribution.rds')
p_sp_f1 <- read_rds('manuscript/figures/R_plot_objects/SPARK_HAQER_F1_replication.rds')
p_sp_iq <- read_rds('manuscript/figures/R_plot_objects/SPARK_HAQER_IQ_replication.rds')
p_sp_rev_dist <- read_rds('manuscript/figures/R_plot_objects/SPARK_rare_reversion_histograms.rds')
p_sp_rev_forest <- read_rds('manuscript/figures/R_plot_objects/SPARK_rare_reversion_forest.rds')

## manually merge above figures in BioRender (figure alignment is wonky in R)

# ## merge plots (A-G)
# plots <- (p_es_pgs_forest | p_f1 | p_f3) / (p_sp_rev_dist | p_sp_rev_forest) + 
#   plot_annotation(tag_levels = "a") & theme(plot.tag = element_text(face = 'bold'))

# ggsave(plots, 
#        filename = 'manuscript/figures/paper_figures/fig2.pdf', 
#        device = 'pdf', 
#        units = 'in', width = 24, height = 14,
#        dpi = 300)
