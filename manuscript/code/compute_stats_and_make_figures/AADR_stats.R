library(tidyverse)

cl <- c("#762776", "#e04468", "#dcc699") ## color palette (HAQER, HAR, RAND)

## read in data
kg_df <- read_csv('manuscript/supplemental_materials/1000_genomes_eur_ES-PGS_with_AADR_SNPs_data.csv') ## 1000 genomes Europeans
aadr_df <- read_csv('manuscript/supplemental_materials/AADR_data.csv') ## ancient humans + neanderthals

## make density plots with neanderthals vs 1000 Genomes Europeans
nean_id <- c('AltaiNea', 'Chagyrskaya-Phalanx', 'DenisovaPinky', 'Vindija33.19')
tmp <- aadr_df %>% 
    filter(IID %in% nean_id) %>% 
    mutate(type = 'Neanderthals and Denisovan')

p_nean_haq <- kg_df %>% 
    mutate(type = 'Modern Europeans') %>%
    ggplot(aes(x = cp_pgs.HAQER)) +
    geom_density(aes(fill = type), alpha = 0.7) +
    xlab('HAQER CP-PGS') +
    ylab('Density') +
    theme_classic() +
    theme(axis.text = element_text(size = 10),
          axis.title = element_text(size = 10),
          strip.text = element_text(size = 12),
          legend.text = element_text(size = 10),
          legend.position = 'bottom') +
    geom_vline(data = tmp, aes(xintercept = cp_pgs.HAQER, color = type), size = 1, linetype = 'dashed') +
    scale_color_manual(name = NULL, values = c("Neanderthals and Denisovan" = "black")) +
    scale_fill_manual(name = NULL, values = c('Modern Europeans' = cl[1])) +
    guides(fill = guide_legend(order = 1, override.aes = list(alpha = 0.7)),
           color = guide_legend(order = 2, override.aes = list(size = 1.1, linetype = 'solid')))

p_nean_bg <- kg_df %>% 
    mutate(type = 'Modern Europeans') %>%
    ggplot(aes(x = cp_pgs.background_HAQER)) +
    geom_density(aes(fill = type), alpha = 0.7) +
    xlab('Background CP-PGS') +
    ylab('Density') +
    theme_classic() +
    theme(axis.text = element_text(size = 10),
          axis.title = element_text(size = 10),
          strip.text = element_text(size = 12),
          legend.text = element_text(size = 10),
          legend.position = 'bottom') +
    geom_vline(data = tmp, aes(xintercept = cp_pgs.background_HAQER, color = type), size = 1, linetype = 'dashed') +
    scale_color_manual(name = NULL, values = c("Neanderthals and Denisovan" = "black")) +
    scale_fill_manual(name = NULL, values = c('Modern Europeans' = 'grey75')) +
    guides(fill = guide_legend(order = 1, override.aes = list(alpha = 0.7)),
           color = guide_legend(order = 2, override.aes = list(size = 1.1, linetype = 'solid')))

## save density plots
p_nean_bg %>% 
    ggsave(filename = 'manuscript/figures/AADR_neanderthal_background_ES-PGS.png',
           device = 'png', dpi = 300, bg = 'white',
           units = 'in', width = 5, height = 5)
p_nean_haq %>% 
    ggsave(filename = 'manuscript/figures/AADR_neanderthal_HAQER_ES-PGS.png',
           device = 'png', dpi = 300, bg = 'white',
           units = 'in', width = 5, height = 5)

#######################################################
## polygenic selection analysis (correlate age with PGS)
#######################################################
## log transform sample age for selection analysis
aadr_df$log_age <- log10(aadr_df$sample_age_years_before_1950)

## polygenic selection
aadr_polygenic_selection_stats <- aadr_df %>% 
    select(IID, log_age, matches('pgs')) %>% 
    pivot_longer(cols = matches('pgs')) %>%
    group_by(name) %>% 
    do(res = broom::tidy(cor.test(.$log_age, .$value))) %>% 
    unnest(res) %>% 
    arrange(p.value)

## flip correlation coefficient to get selection direction since we just correlated sample age with PGS
es_pgs_selection <- aadr_polygenic_selection_stats %>% 
    filter(str_detect(name, 'cp_pgs.')) %>% 
    mutate(selection_coefficient = -1 * estimate, statistic = -1 * statistic, conf.low = -1 * conf.low, conf.high = -1 * conf.high) %>% 
    mutate(pheno = 'sample_age_log',
           n = parameter + 2) %>%
    select(pheno, x = name, selection_coefficient, statistic, p.value, conf.low, conf.high, method, n)
es_pgs_selection %>%
    write_csv('manuscript/supplemental_materials/stats/AADR_HAQER_ES-PGS_polygenic_selection_results.csv')


sel_pgs_haqer <- es_pgs_selection %>% 
    filter(x == 'cp_pgs.HAQER') %>% 
    mutate(lab = str_c('Selection coef. = ', round(selection_coefficient, digits = 2), ', p-val = ', formatC(p.value, digits = 2)))
sel_pgs_bg <- es_pgs_selection %>% 
    filter(x == 'cp_pgs.background_HAQER') %>% 
    mutate(lab = str_c('Selection coef. = ', round(selection_coefficient, digits = 2), ', p-val = ', formatC(p.value, digits = 2)))

## make figure
p_sel <- aadr_df %>%
    select(IID, log_age, matches('cp_pgs')) %>%
    pivot_longer(cols = matches('cp_pgs'), names_to = 'x') %>%
    mutate(x = case_when(x == 'cp_pgs.HAQER' ~ 'HAQERs',
                         x == 'cp_pgs.background_HAQER' ~ 'Background')) %>%
    ggplot(aes(x = -1 * log_age, y = value, color = x)) +
    geom_smooth(size = 1.5) +
    geom_text(aes(x = -3.5, y = -1.75, label = sel_pgs_bg$lab), check_overlap = TRUE, size = 4, color = 'grey50') +
    geom_text(aes(x = -3.5, y = 1.5, label = sel_pgs_haqer$lab), check_overlap = TRUE, size = 4, color = '#762776') +
    xlab('Years ago') +
    ylab('CP-PGS') +
    scale_color_manual(values = c('grey75', "#762776")) +
    theme_classic() +
    labs(color = NULL) +
    theme(axis.text = element_text(size = 12),
          axis.title = element_text(size = 14),
          strip.text = element_text(size = 14),
          legend.text = element_text(size = 12),
          plot.title = element_text(size = 14, hjust = .5),
          legend.position = 'bottom') +
    scale_x_continuous(breaks = c(-5, -4, -3, -2), labels = c('100,000', '10,000', '1,000', '100'))
p_sel %>% 
    ggsave(filename = 'manuscript/figures/AADR_HAQER_polygenic_selection.png',
           device = 'png', dpi = 300, bg = 'white',
           units = 'in', width = 5, height = 5)

## genetic correlations
genetic_corr_res <- aadr_df %>% 
    select(IID, matches('pgs')) %>% 
    pivot_longer(cols = matches('pgs.genome_wide')) %>% 
    group_by(name) %>% 
    do(res = broom::tidy(lm(value ~ cp_pgs.HAQER + cp_pgs.background_HAQER, data = .))) %>% 
    unnest(res) %>% 
    filter(term != '(Intercept)') %>% 
    filter(str_detect(name, 'intra|birth|nrep|nviq|read|spell|wr|phoneme|word')) %>% 
    rename(beta = estimate, x = term, pheno = name) %>% 
    mutate(n = nrow(aadr_df))
genetic_corr_res %>% 
    write_csv('manuscript/supplemental_materials/stats/AADR_ES-PGS_polygenic_correlation_results.csv')

## make forest plot of HAQER vs background PGS correlations
p_gen_corr <- genetic_corr_res %>% 
    arrange(p.value) %>% 
    mutate(x = case_when(x == 'cp_pgs.HAQER' ~ 'HAQERs',
                         x == 'cp_pgs.background_HAQER' ~ 'Background'),
           pheno = str_remove_all(pheno, '_pgs.genome_wide'),
           type = case_when(str_detect(pheno, 'intra|birth') ~ 'Developmental',
                            str_detect(pheno, 'nrep|nviq|read|spell|phon') ~ 'Cognitive'),
            clean_name = case_when(pheno == 'birth_head_circ' ~ 'Birth head circumference',
                                    pheno == 'birth_weight' ~ 'Birth weight',
                                    pheno == 'intracranial_vol' ~ 'Adult intracranial volume',
                                    pheno == 'nread' ~ 'Nonword reading',
                                    pheno == 'nrep' ~ 'Nonword repetition',
                                    pheno == 'nviq' ~ 'Nonverbal IQ',
                                    pheno == 'spelling' ~ 'Spelling',
                                    pheno == 'word_reading' ~ 'Word reading',
                                    pheno == 'phoneme_awareness' ~ 'Phoneme awareness')
            ) %>%
    mutate(x = factor(x, levels = rev(c('HAQERs', 'Background'))),
           type = factor(type, levels = c('Developmental', 'Cognitive'))) %>%
    arrange(desc(x), beta) %>% 
    mutate(clean_name = factor(clean_name, levels = unique(clean_name))) %>%
    mutate(sig = ifelse(p.value < 0.05, TRUE, FALSE)) %>%
    arrange(beta) %>%
    ggplot(aes(x = beta, y = clean_name, color = x, group = x)) +
    geom_linerange(aes(xmin = beta - 1.96 * std.error, xmax = beta + 1.96 * std.error), size = 1.1, position = position_dodge2(.45)) +
    geom_point(size = 3.5, position = position_dodge2(.45), aes(shape = sig)) +
    scale_shape_manual(values = c(1, 16)) +
    geom_vline(xintercept = 0, color = 'red', linetype = 'dashed', size = 1.075) +
    ggforce::facet_col(facets = vars(type), 
                       scales = "free", 
                       space = "free") +
    labs(color = "CP-PGS:") +
    theme_classic() +
    theme(axis.text = element_text(size = 10),
          axis.title = element_text(size = 10),
          strip.text = element_text(size = 12),
          legend.text = element_text(size = 10),
          legend.position = 'bottom') +
    scale_color_manual(values = rev(c(cl[1], 'grey75'))) +
    ylab('Polygenic score') +
    xlab('Polygenic correlation beta (95% CI)') +
    guides(shape = 'none')

p_gen_corr %>%
    ggsave(filename = 'manuscript/figures/AADR_HAQER_polygenic_correlation.png',
           device = 'png', dpi = 300, bg = 'white',
           units = 'in', width = 5, height = 7)

#########################
## imputed genotype analysis of ancient Europeans
aadr_imp_df <- read_csv('manuscript/supplemental_materials/AADR_imputed_data_no_neanderthals.csv')
## log transform sample age for selection analysis
aadr_imp_df$log_age <- log10(aadr_imp_df$sample_age_years_before_1950)

## polygenic selection
aadr_imp_polygenic_selection_stats <- aadr_imp_df %>% 
    select(IID, log_age, matches('pgs')) %>% 
    pivot_longer(cols = matches('pgs')) %>%
    group_by(name) %>% 
    do(res = broom::tidy(cor.test(.$log_age, .$value))) %>% 
    unnest(res) %>% 
    arrange(p.value)
aadr_imp_polygenic_selection_stats %>% select(1:4)
plot(aadr_imp_df$log_age, aadr_imp_df$cp_pgs.imputed.HAQER)
## flip correlation coefficient to get selection direction since we just correlated sample age with PGS
es_pgs_imp_selection <- aadr_imp_polygenic_selection_stats %>% 
    filter(str_detect(name, 'cp_pgs.imputed.')) %>% 
    mutate(selection_coefficient = -1 * estimate, statistic = -1 * statistic, conf.low = -1 * conf.low, conf.high = -1 * conf.high) %>% 
    mutate(pheno = 'sample_age_log',
           n = parameter + 2) %>%
    select(pheno, x = name, selection_coefficient, statistic, p.value, conf.low, conf.high, method, n)
es_pgs_imp_selection %>%
    write_csv('manuscript/supplemental_materials/stats/AADR_imputed_no_neanderthals_HAQER_ES-PGS_polygenic_selection_results.csv')
sel_pgs_imp_haqer <- es_pgs_imp_selection %>% 
    filter(x == 'cp_pgs.imputed.HAQER') %>% 
    mutate(lab = str_c('Selection coef. = ', round(selection_coefficient, digits = 2), ', p-val = ', formatC(p.value, digits = 2)))
sel_pgs_imp_bg <- es_pgs_imp_selection %>% 
    filter(x == 'cp_pgs.imputed.background_HAQER') %>% 
    mutate(lab = str_c('Selection coef. = ', round(selection_coefficient, digits = 2), ', p-val = ', formatC(p.value, digits = 2)))

## make figure
p_sel_imp <- aadr_imp_df %>%
    select(IID, log_age, matches('cp_pgs')) %>%
    pivot_longer(cols = matches('cp_pgs'), names_to = 'x') %>%
    mutate(x = case_when(x == 'cp_pgs.imputed.HAQER' ~ 'HAQERs',
                         x == 'cp_pgs.imputed.background_HAQER' ~ 'Background')) %>%
    mutate(sel_pgs_imp_bg = sel_pgs_imp_bg$lab,
           sel_pgs_imp_haqer = sel_pgs_imp_haqer$lab) %>%
    ggplot(aes(x = -1 * log_age, y = value, color = x)) +
    geom_smooth(size = 1.5, method = 'lm') +
    geom_text(aes(x = -3.1, y = 0, label = sel_pgs_imp_bg), check_overlap = TRUE, size = 4, color = 'grey50') +
    geom_text(aes(x = -3.1, y = 4, label = sel_pgs_imp_haqer), check_overlap = TRUE, size = 4, color = '#762776') +
    xlab('Years ago') +
    ylab('CP-PGS') +
    scale_color_manual(values = c('grey75', "#762776")) +
    theme_classic() +
    labs(color = NULL) +
    theme(axis.text = element_text(size = 12),
          axis.title = element_text(size = 14),
          strip.text = element_text(size = 14),
          legend.text = element_text(size = 12),
          plot.title = element_text(size = 14, hjust = .5),
          legend.position = 'bottom') +
    scale_x_continuous(breaks = c(-5, -4, -3, -2), labels = c('100,000', '10,000', '1,000', '100'))
p_sel_imp %>% 
    ggsave(filename = 'manuscript/figures/AADR_HAQER_polygenic_selection_imputed.png',
           device = 'png', dpi = 300, bg = 'white',
           units = 'in', width = 5, height = 5)



#########################
## save plot objects
tmp %>% 
    write_rds('manuscript/figures/R_plot_objects/AADR_background-CP-PGS_nean_data.rds')
p_nean_bg %>% 
    write_rds('manuscript/figures/R_plot_objects/AADR_background-CP-PGS_nean_dist.rds')
p_nean_haq %>% 
    write_rds('manuscript/figures/R_plot_objects/AADR_HAQER-CP-PGS_nean_dist.rds')
p_sel %>% 
    write_rds('manuscript/figures/R_plot_objects/AADR_HAQER-CP-PGS_selection.rds')
p_gen_corr %>% 
    write_rds('manuscript/figures/R_plot_objects/AADR_HAQER-CP-PGS_polygenic_correlation.rds')
p_sel_imp %>% 
    write_rds('manuscript/figures/R_plot_objects/AADR_HAQER-CP-PGS_selection_imputed.rds')