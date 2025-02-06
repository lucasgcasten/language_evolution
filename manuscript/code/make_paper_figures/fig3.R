library(tidyverse)
library(patchwork)

cl <- c("#762776", "#e04468", "#dcc699") ## color palette (HAQER, HAR, RAND)

## read in plot objects
# ABCD stuff
p_abcd_icv <- read_rds('manuscript/figures/R_plot_objects/ABCD_HAQER-CP-PGS_ICV.rds')
p_abcd_icv_growth <- read_rds('manuscript/figures/R_plot_objects/ABCD_HAQER-CP-PGS_ICV_growth.rds')  
p_sib_bw <- read_rds('manuscript/figures/R_plot_objects/ABCD_HAQER-CP-PGS_sibling_birth_weight.rds')
p_bw_lines <- read_rds('manuscript/figures/R_plot_objects/ABCD_HAQER-CP-PGS_sibling_birth_weight_lines.rds')
p_paired_bg <- read_rds('manuscript/figures/R_plot_objects/ABCD_HAQER-CP-PGS_sibling_csection.rds')
p_paired_haq <- read_rds('manuscript/figures/R_plot_objects/ABCD_background-CP-PGS_sibling_csection.rds')
# AADR stuff
tmp <- read_rds('manuscript/figures/R_plot_objects/AADR_background-CP-PGS_nean_data.rds')
p_nean_bg <- read_rds('manuscript/figures/R_plot_objects/AADR_background-CP-PGS_nean_dist.rds')
p_nean_haq <- read_rds('manuscript/figures/R_plot_objects/AADR_HAQER-CP-PGS_nean_dist.rds')
p_sel <- read_rds('manuscript/figures/R_plot_objects/AADR_HAQER-CP-PGS_selection.rds')
p_gen_corr <- read_rds('manuscript/figures/R_plot_objects/AADR_HAQER-CP-PGS_polygenic_correlation.rds') # + guides(color = guide_legend(reverse = TRUE))
## AADR stats for fig
es_pgs_selection <- read_csv('manuscript/supplemental_materials/stats/AADR_HAQER_ES-PGS_polygenic_selection_results.csv')
sel_pgs_haqer <- es_pgs_selection %>% 
    filter(x == 'cp_pgs.HAQER') %>% 
    mutate(lab = str_c('Selection coef. = ', round(selection_coefficient, digits = 2), ', p-val = ', formatC(p.value, digits = 2)))
sel_pgs_bg <- es_pgs_selection %>% 
    filter(x == 'cp_pgs.background_HAQER') %>% 
    mutate(lab = str_c('Selection coef. = ', round(selection_coefficient, digits = 2), ', p-val = ', formatC(p.value, digits = 2)))

## merge figures
layout = "
AD
BD
CD"

fig3 <- p_nean_haq +
        p_nean_bg +
        p_sel +
        p_gen_corr +
         plot_layout(design = layout) + plot_annotation(tag_levels = 'A') & 
         theme(plot.tag = element_text(size = 24, face = 'bold'),
               axis.text = element_text(size = 12),
               axis.title = element_text(size = 14),
               legend.text = element_text(size = 12),
               legend.title = element_text(size = 14),
               strip.text = element_text(size = 14))

## save figure
ggsave(fig3, 
       filename = 'manuscript/figures/paper_figures/fig3.pdf', 
       device = 'pdf', 
       units = 'in', width = 12.85, height = 14,
       dpi = 300)
