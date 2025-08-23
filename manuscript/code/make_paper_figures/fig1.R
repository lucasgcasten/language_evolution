library(tidyverse)
library(GGally)
library(gridExtra)
library(patchwork)

## read in plot objects
plot_a <- read_rds('manuscript/figures/R_plot_objects/EpiSLI_factor_loadings.rds') + xlab(NULL) + ylab(NULL)
plot_b <- read_rds('manuscript/figures/R_plot_objects/EpiSLI_factor_correlations.rds')
table_a <- read_rds('manuscript/figures/R_plot_objects/EpiSLI_factor_names.rds')
p_cbcl <- read_rds('manuscript/figures/R_plot_objects/EpiSLI_factor_CBCL_correlations.rds')
p_pgs <- read_rds('manuscript/figures/R_plot_objects/EpiSLI_factor_PGS_correlations.rds')
p_es_pgs_forest <- read_rds('manuscript/figures/R_plot_objects/EpiSLI_factor_ES-PGS_forest.rds')
p_f1 <- read_rds('manuscript/figures/R_plot_objects/EpiSLI_factor1_HAQER-CP-PGS.rds')
p_f2 <- read_rds('manuscript/figures/R_plot_objects/EpiSLI_factor2_HAQER-CP-PGS.rds')
p_f3 <- read_rds('manuscript/figures/R_plot_objects/EpiSLI_factor3_HAQER-CP-PGS.rds')

###############################
## merge figures to make fig 1
tmp <- tribble(~Factor, ~Description, 
  'F1', 'Core language',
  'F2', 'Receptive language', 
  'F3', 'Nonverbal IQ', 
  'F4', 'Early language', 
  'F5', 'Talkativeness', 
  'F6', 'Instruction comprehension', 
  'F7', 'Vocabulary')


layout <- "
AAAAABBBBB
CCCCDDEEEE
"

mytheme <- gridExtra::ttheme_default(
    core = list(fg_params=list(cex = 1.75)),
    colhead = list(fg_params=list(cex = 1.75)),
    rowhead = list(fg_params=list(cex = 1.75)))
tab2 = tableGrob(tmp,  rows = NULL, theme = mytheme)

## merge together with patchwork
fig1 <- plot_a +  
        ggmatrix_gtable(plot_b, gg = theme(axis.text = element_text(size = 12),
                                            axis.title = element_text(size = 12),
                                            legend.text = element_text(size = 10),
                                            legend.title = element_text(size = 12))) + 
         wrap_elements(tab2) +
         plot_spacer() +
         p_pgs +  
         plot_layout(design = layout) + plot_annotation(tag_levels = 'A') & 
         theme(plot.tag = element_text(size = 24, face = 'bold'),
               axis.text = element_text(size = 12),
               axis.title = element_text(size = 14),
               legend.text = element_text(size = 12),
               legend.title = element_text(size = 14),
               strip.text = element_text(size = 14))

## save figure
ggsave(fig1, 
       filename = 'manuscript/figures/paper_figures/fig1.pdf', 
       device = 'pdf', 
       units = 'in', width = 14.25, height = 14.25,
       dpi = 300)

ggsave(fig1, 
       filename = 'manuscript/figures/paper_figures/fig1.png', 
       device = 'png', 
       units = 'in', width = 14.25, height = 14.25,
       dpi = 300)
