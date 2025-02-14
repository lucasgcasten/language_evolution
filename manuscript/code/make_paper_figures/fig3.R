library(tidyverse)
library(patchwork)

cl <- c("#762776", "#e04468", "#dcc699") ## color palette (HAQER, HAR, RAND)

## read in plot objects
# neanderthal + denisovan WGS PGS figures
p_nean_bg <- read_rds('manuscript/figures/R_plot_objects/AADR_background-CP-PGS_nean_dist.rds')
p_nean_haq <- read_rds('manuscript/figures/R_plot_objects/AADR_HAQER-CP-PGS_nean_dist.rds')
# AADR stuff
p_sel <- read_rds('manuscript/figures/R_plot_objects/AADR_HAQER-CP-PGS_selection.rds')

## merge figures
layout = "
AAABB
AAACC"

fig3 <- p_sel + 
        p_nean_haq +
        p_nean_bg +
        plot_layout(design = layout) + plot_annotation(tag_levels = 'a') & 
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
       units = 'in', width = 14, height = 8,
       dpi = 300)

## old figure (had some additional results that didn't contribute a ton to the overall story)
# layout = "
# AD
# BD
# CD"
# fig3 <- p_nean_haq +
#         p_nean_bg +
#         p_sel_greml +
#         p_gen_corr +
#          plot_layout(design = layout) + plot_annotation(tag_levels = 'a') & 
#          theme(plot.tag = element_text(size = 24, face = 'bold'),
#                axis.text = element_text(size = 12),
#                axis.title = element_text(size = 14),
#                legend.text = element_text(size = 12),
#                legend.title = element_text(size = 14),
#                strip.text = element_text(size = 14))