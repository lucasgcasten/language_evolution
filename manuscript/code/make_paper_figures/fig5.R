library(tidyverse)
library(patchwork)

## read in plot objects
# neanderthal + denisovan WGS PGS figures
p_nean_bg <- read_rds('manuscript/figures/R_plot_objects/AADR_background-CP-PGS_nean_dist.rds')
p_nean_haq <- read_rds('manuscript/figures/R_plot_objects/AADR_HAQER-CP-PGS_nean_dist.rds')
p_nean_matched <- read_rds('manuscript/figures/R_plot_objects/AADR_matched-CP-PGS_nean_dist.rds')
p_nean_merged <- read_rds("manuscript/figures/R_plot_objects/AADR_HAQER-CP-PGS_nean_dist_all_merged.rds")
# AADR stuff
p_sel <- read_rds('manuscript/figures/R_plot_objects/AADR_HAQER-CP-PGS_selection.rds')
## EpiSLI balancing selection
p_sfs <- read_rds("manuscript/figures/R_plot_objects/HAQER_SFS_bin_comparisons.rds")
p_f_comp <- read_rds("manuscript/figures/R_plot_objects/EpiSLI_fstat_comparisons.rds")
## ABCD cephalopelvic disprortion
p_cpd <- read_rds('manuscript/figures/R_plot_objects/ABCD_cephalopelvic_disproportion_HAQER_ES-PGS.rds') + scale_fill_manual(values = c('grey70', 'seagreen'))
## cross species validation figures
p_vl <- read_rds('manuscript/figures/R_plot_objects/vocal_learners_HAQER-like_no_primates.rds')
p_bs <- read_rds('manuscript/figures/R_plot_objects/brain_size_HAQER-like_no_primates.rds')
p_wt <- read_rds('manuscript/figures/R_plot_objects/birth_weight_ratio_HAQER-like_no_primates.rds')    

## merge figures
# layout = "AAAABBCCDD
# AAAAEEEEFF"
layout = "AAAABBBBBB
AAAACCCCDD
EEEFFFGGII"

fig5 <- (p_sel +  
        p_nean_merged +
        p_sfs +
        p_f_comp +
        p_cpd + 
        p_vl + 
        p_bs + 
        p_wt) +
        plot_layout(design = layout) + plot_annotation(tag_levels = 'A') & 
        theme(plot.tag = element_text(size = 24, face = 'bold'),
              axis.text = element_text(size = 18),
              axis.title = element_text(size = 20),
              legend.text = element_text(size = 18),
              legend.title = element_text(size = 20),
              strip.text = element_text(size = 20))

## save figure
ggsave(fig5, 
       filename = 'manuscript/figures/paper_figures/fig5.pdf', 
       device = 'pdf', 
       units = 'in', width = 28, height = 21,
       dpi = 300)
ggsave(fig5, 
       filename = 'manuscript/figures/paper_figures/fig5.png', 
       device = 'png', 
       units = 'in', width = 28, height = 21,
       dpi = 300)
