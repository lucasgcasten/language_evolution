library(tidyverse)
library(patchwork)

cl <- c("#762776", "#e04468", "#dcc699") ## color palette (HAQER, HAR, RAND)

## read in plot objects
# ABCD stuff
p_abcd_icv <- read_rds('manuscript/figures/R_plot_objects/ABCD_HAQER-CP-PGS_ICV.rds')
p_abcd_icv_growth <- read_rds('manuscript/figures/R_plot_objects/ABCD_HAQER-CP-PGS_ICV_growth.rds')  
p_sib_bw <- read_rds('manuscript/figures/R_plot_objects/ABCD_HAQER-CP-PGS_sibling_birth_weight.rds')
p_sib_bw_bg <- read_rds('manuscript/figures/R_plot_objects/ABCD_background-CP-PGS_sibling_birth_weight.rds')
p_bw_lines <- read_rds('manuscript/figures/R_plot_objects/ABCD_HAQER-CP-PGS_sibling_birth_weight_lines.rds')
p_paired_bg <- read_rds('manuscript/figures/R_plot_objects/ABCD_HAQER-CP-PGS_sibling_csection.rds')
p_paired_haq <- read_rds('manuscript/figures/R_plot_objects/ABCD_background-CP-PGS_sibling_csection.rds')

# HAQER scQTL
p_scqtl_enr <- read_rds('manuscript/figures/R_plot_objects/HAQER_scQTL_enrichment.rds')
# EpiSLI validation
p_scqtl_valid <- read_rds('manuscript/figures/R_plot_objects/HAQER_scQTL_EpiSLI_validation_F1_F3_scatterplot.rds')

## merge plots
layout <- "
AABBCD
AABBEF
AABBGH
"

fig4 <- p_scqtl_enr +
       p_paired_bg +
       p_paired_haq +
       p_sib_bw_bg +
       p_sib_bw +
        p_abcd_icv +
        p_abcd_icv_growth +
        plot_layout(design = layout) + plot_annotation(tag_levels = 'a') & 
        theme(plot.tag = element_text(size = 28, face = 'bold'),
               axis.text = element_text(size = 14),
               axis.title = element_text(size = 16),
               legend.text = element_text(size = 20),
               legend.title = element_text(size = 20),
               strip.text = element_text(size = 20),
               plot.title = element_text(size = 20, face = 'bold'))

## save figure
ggsave(fig4,
       filename = 'manuscript/figures/paper_figures/fig4.pdf', 
       device = 'pdf', 
       units = 'in', width = 30, height = 16,
       dpi = 300)
