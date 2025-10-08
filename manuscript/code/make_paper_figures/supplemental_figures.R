library(tidyverse)
library(patchwork)

cl <- c("#762776", "#e04468", "#dcc699") ## color palette (HAQER, HAR, RAND)

## read in enrichment data
head_enr <- read_csv('manuscript/supplemental_materials/stats/HAQER_birth_head_circ_gwas_Vogelezang2022_enrichment_stats.csv') %>% 
       mutate(type = 'Birth head circumference loci')
vl_enr <- read_csv('manuscript/supplemental_materials/stats/HAQER_vocal_learning_Wirthlin2024_enrichment_stats.csv') %>% 
       mutate(type = 'Mammalian vocal learning enhancers')
ca_enr <- read_csv("manuscript/supplemental_materials/stats/HAQER_chromatin_accessibility_evo_Li2023_enrichment_stats.csv")


## make figure panels
### birth head circumference GWAS enrichment
p_bhc <- head_enr %>% 
    mutate(evo_annot = factor(evo_annot, levels = c('RAND', 'HAR', 'HAQER'))) %>%
    arrange(evo_annot, desc(enrichment_p)) %>%
    mutate(set = str_replace_all(set, pattern = '__', replacement = ' ')) %>%
    mutate(set = factor(set, levels = unique(set))) %>% 
    drop_na() %>%
    mutate(enrichment_p = ifelse(enrichment_p > 0.85, .85, enrichment_p)) %>% ## add a small constant value so they show up on the plot
    ggplot(aes(x = -log10(enrichment_p), y = type)) +
    geom_bar(stat = 'identity', aes(fill = evo_annot), position = 'dodge') +
    geom_vline(xintercept = -log10(0.05), color = 'red', linetype = 'dashed', size = 1.075) +
    labs(fill = NULL) +
    theme_classic() +
    theme(axis.text = element_text(size = 12),
          axis.title = element_text(size = 14),
          strip.text = element_text(size = 14),
          legend.text = element_text(size = 12),
          plot.title = element_text(size = 14, hjust = .5),
          legend.position = 'bottom') +
    scale_fill_manual(values = rev(cl)) +
    xlab('-log10(enrichment p-value)') +
    ylab(NULL)

## make figure panels
p_vl <- vl_enr %>% 
    mutate(evo_annot = factor(evo_annot, levels = c('RAND', 'HAR', 'HAQER'))) %>%
    arrange(evo_annot, desc(enrichment_p)) %>%
    mutate(set = str_replace_all(set, pattern = '__', replacement = ' ')) %>%
    mutate(set = factor(set, levels = unique(set))) %>% 
    drop_na() %>%
    mutate(enrichment_p = ifelse(enrichment_p > 0.85, .85, enrichment_p)) %>% ## add a small constant value so they show up on the plot
    ggplot(aes(x = -log10(enrichment_p), y = type)) +
    geom_bar(stat = 'identity', aes(fill = evo_annot), position = 'dodge') +
    geom_vline(xintercept = -log10(0.05), color = 'red', linetype = 'dashed', size = 1.075) +
    labs(fill = NULL) +
    theme_classic() +
    theme(axis.text = element_text(size = 12),
          axis.title = element_text(size = 14),
          strip.text = element_text(size = 14),
          legend.text = element_text(size = 12),
          plot.title = element_text(size = 14, hjust = .5),
          legend.position = 'bottom') +
    scale_fill_manual(values = rev(cl)) +
    xlab('-log10(enrichment p-value)') +
    ylab(NULL)

## make chromatin accessibility evo figs
hs_df <- ca_enr %>% 
    filter(type == 'human_specific_chromatin_accessibility') %>%
    mutate(evo_annot = factor(evo_annot, levels = c('RAND', 'HAR', 'HAQER'))) %>%
    arrange(desc(evo_annot), desc(enrichment_p)) %>%
    mutate(chromatin_accessibility_set = factor(chromatin_accessibility_set, levels = unique(chromatin_accessibility_set))) %>% 
    drop_na() %>%
    mutate(enrichment_p = ifelse(enrichment_p > 0.65, .65, enrichment_p)) %>% ## add a small constant value so they show up on the plot
    mutate(class = factor(class, levels = c('GABA+', 'vGlut+', 'Non-neuronal')))

p_hs <- hs_df %>% 
    ggplot(aes(x = -log10(enrichment_p), y = chromatin_accessibility_set)) +
    geom_bar(stat = 'identity', aes(fill = evo_annot), position = 'dodge') +
    geom_vline(xintercept = -log10(0.05), color = 'red', linetype = 'dashed', size = 1.075) +
    ggforce::facet_col(facets = vars(class), 
                       scales = "free_y", 
                       space = "free") +
    labs(fill = NULL) +
    theme_classic() +
    theme(axis.text = element_text(size = 12),
          axis.title = element_text(size = 14),
          strip.text = element_text(size = 14),
          legend.text = element_text(size = 12),
          plot.title = element_text(size = 14, hjust = .5),
          legend.position = 'bottom') +
    scale_fill_manual(values = rev(cl)) +
    xlab('-log10(enrichment p-value)') +
    ylab('Cell-type') +
    ggtitle('Human-specific cCREs')

p_cons <- ca_enr %>%
    filter(type == 'conserved_chromatin_accessibility') %>%
    mutate(evo_annot = factor(evo_annot, levels = c('RAND', 'HAR', 'HAQER'))) %>%
    mutate(chromatin_accessibility_set = factor(chromatin_accessibility_set, levels = levels(hs_df$chromatin_accessibility_set))) %>% 
    drop_na() %>%
    mutate(enrichment_p = ifelse(enrichment_p > 0.65, .65, enrichment_p)) %>% ## add a small constant value so they show up on the plot
    mutate(class = case_when(str_detect(chromatin_accessibility_set, 'FOXP2|LAMP|MSN|PVALB|SST|VIP') ~ 'GABA+',
                             str_detect(chromatin_accessibility_set, '^CT|^IT|L6|NP') ~ 'vGlut+',
                             str_detect(chromatin_accessibility_set, 'ASCT|MGC|OGC|OPC') ~ 'Non-neuronal'),
           class = factor(class, levels = c('GABA+', 'vGlut+', 'Non-neuronal'))) %>% 
    ggplot(aes(x = -log10(enrichment_p), y = chromatin_accessibility_set)) +
    geom_bar(stat = 'identity', aes(fill = evo_annot), position = 'dodge') +
    geom_vline(xintercept = -log10(0.05), color = 'red', linetype = 'dashed', size = 1.075) +
    ggforce::facet_col(facets = vars(class), 
                       scales = "free_y", 
                       space = "free") +
    labs(fill = NULL) +
    theme_classic() +
    theme(axis.text = element_text(size = 12),
          axis.title = element_text(size = 14),
          strip.text = element_text(size = 14),
          legend.text = element_text(size = 12),
          plot.title = element_text(size = 14, hjust = .5),
          legend.position = 'bottom') +
    scale_fill_manual(values = rev(cl)) +
    xlab('-log10(enrichment p-value)') +
    ylab('Cell-type') +
    ggtitle('Conserved cCREs')


## read in scQTL plot
# HAQER scQTL
p_scqtl_enr <- read_rds('manuscript/figures/R_plot_objects/HAQER_scQTL_enrichment.rds')

## merge panels 
p_merge <- (p_hs + p_cons & theme(legend.position = 'none')) / (p_scqtl_enr & theme(legend.position = 'none')) / ((p_bhc + p_vl) + plot_layout(guides = 'collect') & theme(legend.position = 'bottom'))

## adjust formatting
fig_supp <- p_merge +
        plot_layout(heights = c(1,.7,.1)) +
        plot_annotation(tag_levels = 'A') & 
        theme(plot.tag = element_text(size = 24, face = 'bold'),
               axis.text = element_text(size = 14),
               axis.title = element_text(size = 16),
               legend.text = element_text(size = 20),
               legend.title = element_text(size = 20),
               strip.text = element_text(size = 20),
               plot.title = element_text(size = 20, face = 'bold'))

## save figure
ggsave(fig_supp,
       filename = 'manuscript/figures/paper_figures/supp_fig_HAQER_enrichment.pdf', 
       device = 'pdf', 
       units = 'in', width = 14, height = 21,
       dpi = 300)
ggsave(fig_supp,
       filename = 'manuscript/figures/paper_figures/supp_fig_HAQER_enrichment.png', 
       device = 'png', 
       units = 'in', width = 14, height = 21,
       dpi = 300)

#####################################
## neanderthal, ancient, modern PGSs
#####################################
p_nean_merged <- read_rds("manuscript/figures/R_plot_objects/AADR_HAQER-CP-PGS_nean_dist_all_merged.rds") & 
    scale_y_continuous(limits = c(-5,5)) &
    plot_annotation(tag_levels = 'A') & 
    theme(plot.tag = element_text(size = 24, face = 'bold'))
ggsave(p_nean_merged,
       filename = 'manuscript/figures/paper_figures/supp_fig_AADR_ES-PGS_boxplots.png', 
       device = 'png', 
       units = 'in', width = 14, height = 6,
       dpi = 300)