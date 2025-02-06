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

##########################
## make GWAS plot
gwas_dat <- read_csv('manuscript/supplemental_materials/cogPerf_and_NVIQ_GWAS_HAQER_sumstats_reversions_data.csv')
## initialize empty list to save KS-test sumtats
p_list <- list()
## loop over each GWAS and do analysis / make figure
for(i in c('cogPerf')) { 
    message('Working on: ', i)

    ## subset to sumstats of interest
    wd <- gwas_dat %>% 
        filter(gwas == i) %>% 
        mutate(z = beta / se)
    hca_haqer <- wd %>% 
        filter(type == 'HCA reversion in HAQER') 
    other_haqer <- wd %>% 
        filter(type == 'Conserved base change in HAQER')

    ## do KS test comparing GWAS Z-score of reversions (homo allele > HCA allele) with other variants (conserved > new variant)
    ks_res <- broom::tidy(ks.test(hca_haqer$z, other_haqer$z)) %>% 
        mutate(gwas = i,
               n_ind_HAQER_SNP = nrow(wd),
               n_HCA_reversion_SNP = nrow(hca_haqer),
               n_conserved_variant_SNP = nrow(other_haqer),
               mean_HCA_reversion_zscore = mean(hca_haqer$z),
               mean_conserved_variant_zscore = mean(other_haqer$z)) %>% 
        relocate(gwas)

    ## gather data for ECDF figure
    sample1 <- hca_haqer$z ## SNP z-scores for HCA reversions
    sample2 <- other_haqer$z ## SNP z-scores for other variants types (conserved -> other allele)
    cdf1 <- ecdf(hca_haqer$z) 
    cdf2 <- ecdf(other_haqer$z) 

    # find min and max statistics to draw line between points of greatest distance
    minMax <- seq(min(sample1, sample2), max(sample1, sample2), length.out=length(sample1)) 
    x0 <- minMax[which( abs(cdf1(minMax) - cdf2(minMax)) == max(abs(cdf1(minMax) - cdf2(minMax))) )] 
    y0 <- cdf1(x0) 
    y1 <- cdf2(x0) 

    ## ECDF figure
    trait = case_when(i == 'cogPerf' ~ 'Cognitive performance',
                      i == 'NVIQ' ~ 'Nonverbal IQ')
    p_list[[i]] <- wd %>% 
        mutate(type = case_when(type == 'HCA reversion in HAQER' ~ 'Human-chimp ancestral allele reversion',
                                type == 'Conserved base change in HAQER' ~ 'Human-chimp conserved base change')) %>%
        group_by(type) %>% 
        mutate(n = n()) %>% 
        mutate(type = str_c(type, '\n(N=', n, ' SNPs)')) %>%
        ggplot(aes(x = z, group = type, color = type))+
        stat_ecdf(size=1.2) +
        xlab(str_c(trait, " GWAS zscore")) +
        ylab("ECDF") +
        # geom_segment(aes(x = x0[1], y = y0[1], xend = x0[1], yend = y1[1]),
                    # linetype = "dashed", color = "red", size = 1.2) +
        # geom_point(aes(x = x0[1] , y= y0[1]), color="red", size=3.5) +
        # geom_point(aes(x = x0[1] , y= y1[1]), color="red", size=3.5) +
        theme_classic() +
        theme(axis.text = element_text(size = 12),
              axis.title = element_text(size = 14),
              strip.text = element_text(size = 14),
              legend.text = element_text(size = 12),
              legend.title = element_text(size = 14),
              legend.position = c(0.7, 0.2),
              legend.box.background = element_rect(color = "black")) +
        geom_text(aes(x = 1.75, y = 0.4, label = str_c('D = ', round(ks_res$statistic, digits =2), ', p = ', formatC(ks_res$p.value, digits = 2))), check_overlap = TRUE, color = 'black', size = 5) +
        coord_cartesian(xlim = c(-2.5,3.25)) +
        labs(color = 'HAQER variant type:') +
        geom_vline(xintercept = 0, color = 'red', linetype = 'dashed', size = 1.075)
} 
p_gwas <- p_list[[i]]

## add empty space for Jake's figures
p_jake <- plot_spacer() + theme(plot.background = element_rect(fill = 'grey'))


## merge plots (A-G)
# plots <- (p_es_pgs_forest | p_f1 | p_f3 | p_gwas) / (p_jake | p_sp_rev_dist | p_sp_rev_forest) + 
#   plot_annotation(tag_levels = "A")
plots <- (p_es_pgs_forest | p_f1 | p_f3) / (p_sp_rev_dist | p_sp_rev_forest) + 
  plot_annotation(tag_levels = "A") & theme(plot.tag = element_text(face = 'bold'))

ggsave(plots, 
       filename = 'manuscript/figures/paper_figures/fig2.pdf', 
       device = 'pdf', 
       units = 'in', width = 24, height = 14,
       dpi = 300)
    
plots_top <- (p_es_pgs_forest | p_f1 | p_f3)
plots_bottom <-  (p_sp_rev_dist | p_sp_rev_forest) 

ggsave(plots, 
       filename = 'manuscript/figures/paper_figures/fig2.pdf', 
       device = 'pdf', 
       units = 'in', width = 24, height = 14,
       dpi = 300)

ggsave(plots_top, 
       filename = 'manuscript/figures/paper_figures/fig2_row1.pdf', 
       device = 'pdf', 
       units = 'in', width = 24, height = 7,
       dpi = 300)
    
ggsave(plots_bottom, 
       filename = 'manuscript/figures/paper_figures/fig2_row2.pdf', 
       device = 'pdf', 
       units = 'in', width = 20, height = 7,
       dpi = 300)


## patchwork is wonky (tries too hard to align everything)
# layout <- 
# "AAAAABCD
# EEEEEFGG"

# fig2 <- p_es_pgs_forest +  
#         p_f1 +
#         p_f3 +
#         p_gwas +
#         p_jake +
#         p_sp_rev_dist +
#         p_sp_rev_forest +  
#          plot_layout(design = layout) + plot_annotation(tag_levels = 'A') & 
#          theme(plot.tag = element_text(size = 12),
#                axis.text = element_text(size = 10),
#                axis.title = element_text(size = 12),
#                legend.text = element_text(size = 10),
#                legend.title = element_text(size = 12),
#                strip.text = element_text(size = 12))

# ## save figure
# ggsave(fig2, 
#        filename = 'manuscript/figures/paper_figures/fig2.pdf', 
#        device = 'pdf', 
#        units = 'in', width = 16, height = 10,
#        dpi = 300)
