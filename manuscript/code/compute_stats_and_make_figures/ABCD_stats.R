library(tidyverse)

## read in ABCD imaging data and ES-PGS
abcd_df <- read_csv('manuscript/supplemental_materials/ABCD_brain_imaging_data.csv')

##############################
## ABCD ICV analysis
##############################
## sample size
icv_n <- abcd_df %>% 
    select(IID, intracranial_vol_resid, cp_pgs.HAQER, cp_pgs.background_HAQER) %>% 
    drop_na() %>% 
    nrow(.)
## run ES-PGS
abcd_icv_es_pgs_stats <- broom::tidy(lm(intracranial_vol_resid ~ cp_pgs.HAQER + cp_pgs.background_HAQER, data = abcd_df)) %>% 
    filter(term != '(Intercept)') %>% 
    rename(x = term, beta = estimate) %>%
    mutate(pheno = 'intracranial_volume',
           n = icv_n) %>% 
    relocate(pheno)

## make labels for figure
icv_haq_lab = abcd_icv_es_pgs_stats %>% 
    filter(x == 'cp_pgs.HAQER') %>% 
    mutate(lab = str_c('ES-PGS beta = ', round(beta, digits = 2), ', p-val = ', formatC(p.value, digits = 2)))
icv_bg_lab = abcd_icv_es_pgs_stats %>% 
    filter(x == 'cp_pgs.background_HAQER') %>% 
    mutate(lab = str_c('ES-PGS beta = ', round(beta, digits = 2), ', p-val = ', formatC(p.value, digits = 2)))

## make figure
p_abcd_icv <- abcd_df %>% 
    select(IID, intracranial_vol_resid, cp_pgs.HAQER, cp_pgs.background_HAQER) %>% 
    drop_na() %>% 
    mutate(icv_haq_lab = icv_haq_lab$lab,
           icv_bg_lab = icv_bg_lab$lab) %>%
    pivot_longer(cols = matches('cp_pgs'), names_to = 'x') %>% 
    mutate(x = case_when(x == 'cp_pgs.HAQER' ~ 'HAQERs',
                         x == 'cp_pgs.background_HAQER' ~ 'Background')) %>%
    ggplot(aes(x = intracranial_vol_resid, y = value, color = x)) +
    geom_smooth(method = 'lm', size = 1.5) +
    geom_text(aes(x = -.1, y = .3, label = icv_bg_lab), check_overlap = TRUE, size = 4, color = 'grey50') +
    geom_text(aes(x = 1, y = -.2, label = icv_haq_lab), check_overlap = TRUE, size = 4, color = '#762776') +
    xlab('Intracranial volume') +
    ylab('CP-PGS') +
    scale_color_manual(values = c('grey75', "#762776")) +
    theme_classic() +
    labs(color = NULL) +
    theme(axis.text = element_text(size = 12),
          axis.title = element_text(size = 14),
          strip.text = element_text(size = 14),
          legend.text = element_text(size = 12),
          plot.title = element_text(size = 14, hjust = .5),
          legend.position = 'bottom')
## save figure
p_abcd_icv %>% 
    ggsave(filename = 'manuscript/figures/ABCD_ICV_HAQER_ES-PGS.png',
           device = 'png', dpi = 300, bg = 'white', 
           units = 'in', width = 5, height = 5)

##############################
## ABCD brain growth analysis
##############################
## sample size
icv_n <- abcd_df %>% 
    select(IID, intracranial_growth_resid, cp_pgs.HAQER, cp_pgs.background_HAQER) %>% 
    drop_na() %>% 
    nrow(.)
## run ES-PGS
abcd_icv_growth_es_pgs_stats <- broom::tidy(lm(intracranial_growth_resid ~ cp_pgs.HAQER + cp_pgs.background_HAQER, data = abcd_df)) %>% 
    filter(term != '(Intercept)') %>% 
    rename(x = term, beta = estimate) %>%
    mutate(pheno = 'intracranial_growth_adolescence',
           n = icv_n) %>% 
    relocate(pheno)

## make labels for figure
icv_haq_lab = abcd_icv_growth_es_pgs_stats %>% 
    filter(x == 'cp_pgs.HAQER') %>% 
    mutate(lab = str_c('ES-PGS beta = ', round(beta, digits = 2), ', p-val = ', formatC(p.value, digits = 2)))
icv_bg_lab = abcd_icv_growth_es_pgs_stats %>% 
    filter(x == 'cp_pgs.background_HAQER') %>% 
    mutate(lab = str_c('ES-PGS beta = ', round(beta, digits = 2), ', p-val = ', formatC(p.value, digits = 2)))

## make figure
p_abcd_icv_growth <- abcd_df %>% 
    select(IID, intracranial_growth_resid, cp_pgs.HAQER, cp_pgs.background_HAQER) %>% 
    drop_na() %>% 
    mutate(icv_haq_lab = icv_haq_lab$lab,
           icv_bg_lab = icv_bg_lab$lab) %>%
    pivot_longer(cols = matches('cp_pgs'), names_to = 'x') %>% 
    mutate(x = case_when(x == 'cp_pgs.HAQER' ~ 'HAQERs',
                         x == 'cp_pgs.background_HAQER' ~ 'Background')) %>%
    ggplot(aes(x = intracranial_growth_resid, y = value, color = x)) +
    geom_smooth(method = 'lm', size = 1.5) +
    geom_text(aes(x = -.1, y = .25, label = icv_bg_lab), check_overlap = TRUE, size = 4, color = 'grey50') +
    geom_text(aes(x = .5, y = -.15, label = icv_haq_lab), check_overlap = TRUE, size = 4, color = '#762776') +
    xlab('Postnatal intracranial growth') +
    ylab('CP-PGS') +
    scale_color_manual(values = c('grey75', "#762776")) +
    theme_classic() +
    labs(color = NULL) +
    theme(axis.text = element_text(size = 12),
          axis.title = element_text(size = 14),
          strip.text = element_text(size = 14),
          legend.text = element_text(size = 12),
          plot.title = element_text(size = 14, hjust = .5),
          legend.position = 'bottom')
## save figure
p_abcd_icv_growth %>% 
    ggsave(filename = 'manuscript/figures/ABCD_ICV_postnatal_growth_HAQER_ES-PGS.png',
           device = 'png', dpi = 300, bg = 'white', 
           units = 'in', width = 5, height = 5)

#################################################################
## sibling birth weight differences vs HAQER CP-PGS difference
#################################################################
abcd_df <- read_csv('manuscript/supplemental_materials/ABCD_sibling_birth_phenotypes_long.csv') %>% 
    filter(pheno == 'bw')
n_sib_bw <- nrow(abcd_df)
abcd_sib_bw_es_pgs_res <- broom::tidy(lm(sib_diff ~ sib_diff_cp_pgs.HAQER + sib_diff_cp_pgs.background_HAQER + age_diff + premature_diff + as.factor(sex_female_diff), data = abcd_df)) %>% 
    filter(term != '(Intercept)') %>% 
    mutate(pheno = 'sibling_difference_in_birth_weight',
           n = n_sib_bw) %>%
    rename(x = term) %>% 
    relocate(pheno) %>% 
    rename(beta = estimate)
## scatterplot with just HAQER CP-PGS
abcd_sib_bw_lab <- abcd_sib_bw_es_pgs_res %>%
    filter(x == 'sib_diff_cp_pgs.HAQER') %>% 
    mutate(lab = str_c('ES-PGS beta = ', round(beta, digits = 2), ', p-val = ', formatC(p.value, digits = 2), '\nN = ', n, ' sibling pairs'))
p_sib_bw <- abcd_df %>% 
    mutate(sib_diff_resid = scale(resid(lm(sib_diff ~  sib_diff_cp_pgs.background_HAQER + age_diff + premature_diff + as.factor(sex_female_diff))))[,1]) %>% 
    mutate(abcd_sib_bw_lab = abcd_sib_bw_lab$lab) %>%
    ggplot(aes(x = sib_diff_resid, y = sib_diff_cp_pgs.HAQER)) +
    geom_point(size = 0.5, color = "#762776") +
    geom_smooth(method = 'lm', size = 1.5, color = "#762776") +
    geom_text(aes(x = 0, y = 3, label = abcd_sib_bw_lab), check_overlap = TRUE, size = 4, color = "#762776") +
    xlab('Sib. difference in birth weight') +
    ylab('Sib. difference in HAQER CP-PGS') +
    theme_classic() +
    labs(color = NULL) +
    theme(axis.text = element_text(size = 12),
          axis.title = element_text(size = 14),
          strip.text = element_text(size = 14),
          legend.text = element_text(size = 12),
          plot.title = element_text(size = 14, hjust = .5),
          legend.position = 'bottom')
p_sib_bw %>% 
    ggsave(filename = 'manuscript/figures/ABCD_sibling_birth_weight_HAQER_ES-PGS.png',
           device = 'png', dpi = 300, bg = 'white', 
           units = 'in', width = 5, height = 4.5)

## scatterplot with just background CP-PGS
abcd_sib_bw_lab <- abcd_sib_bw_es_pgs_res %>%
    filter(x == 'sib_diff_cp_pgs.background_HAQER') %>% 
    mutate(lab = str_c('ES-PGS beta = ', round(beta, digits = 2), ', p-val = ', formatC(p.value, digits = 2), '\nN = ', n, ' sibling pairs'))
p_sib_bw_bg <- abcd_df %>% 
    mutate(sib_diff_resid = scale(resid(lm(sib_diff ~  sib_diff_cp_pgs.HAQER + age_diff + premature_diff + as.factor(sex_female_diff))))[,1]) %>% 
    mutate(abcd_sib_bw_lab = abcd_sib_bw_lab$lab) %>%
    ggplot(aes(x = sib_diff_resid, y = sib_diff_cp_pgs.background_HAQER)) +
    geom_point(size = 0.5, color = "grey75") +
    geom_smooth(method = 'lm', size = 1.5, color = "grey75") +
    geom_text(aes(x = 0, y = 3, label = abcd_sib_bw_lab), check_overlap = TRUE, size = 4, color = "grey50") +
    xlab('Sib. difference in birth weight') +
    ylab('Sib. difference in background CP-PGS') +
    theme_classic() +
    labs(color = NULL) +
    theme(axis.text = element_text(size = 12),
          axis.title = element_text(size = 14),
          strip.text = element_text(size = 14),
          legend.text = element_text(size = 12),
          plot.title = element_text(size = 14, hjust = .5),
          legend.position = 'bottom')
p_sib_bw_bg %>% 
    ggsave(filename = 'manuscript/figures/ABCD_sibling_birth_weight_background_ES-PGS.png',
           device = 'png', dpi = 300, bg = 'white', 
           units = 'in', width = 5, height = 4.5)

## show best fit lines for both HAQER + background on same figure
abcd_sib_bw_bg_lab <- abcd_sib_bw_es_pgs_res %>%
    filter(x == 'sib_diff_cp_pgs.background_HAQER') %>% 
    mutate(lab = str_c('ES-PGS beta = ', round(beta, digits = 2), ', p-val = ', formatC(p.value, digits = 2)))
abcd_sib_bw_haq_lab <- abcd_sib_bw_es_pgs_res %>%
    filter(x == 'sib_diff_cp_pgs.HAQER') %>% 
    mutate(lab = str_c('ES-PGS beta = ', round(beta, digits = 2), ', p-val = ', formatC(p.value, digits = 2)))

p_bw_lines <- abcd_df %>%
    mutate(sib_diff_resid = scale(resid(lm(sib_diff ~  sib_diff_cp_pgs.background_HAQER + age_diff + premature_diff + as.factor(sex_female_diff))))[,1]) %>% 
    mutate(abcd_sib_bw_bg_lab = abcd_sib_bw_bg_lab$lab,
           abcd_sib_bw_haq_lab = abcd_sib_bw_haq_lab$lab) %>%
    pivot_longer(cols = matches('cp_pgs'), names_to = 'x') %>% 
    mutate(x = case_when(x == 'sib_diff_cp_pgs.HAQER' ~ 'HAQERs',
                         x == 'sib_diff_cp_pgs.background_HAQER' ~ 'Background')) %>%
    ggplot(aes(x = sib_diff_resid, y = value, color = x)) +
    geom_smooth(method = 'lm', size = 1.5, aes(color = x)) +
    geom_text(aes(x = 0, y = .45, label = abcd_sib_bw_haq_lab), check_overlap = TRUE, size = 4, color = '#762776') +
    geom_text(aes(x = 0, y = -.5, label = abcd_sib_bw_bg_lab), check_overlap = TRUE, size = 4, color = 'grey50') +
    xlab('Sib. difference in birth weight') +
    ylab('Sib. difference in CP-PGS') +
    scale_color_manual(values = c('grey75', "#762776")) +
    theme_classic() +
    labs(color = NULL) +
    theme(axis.text = element_text(size = 12),
          axis.title = element_text(size = 14),
          strip.text = element_text(size = 14),
          legend.text = element_text(size = 12),
          plot.title = element_text(size = 14, hjust = .5),
          legend.position = 'bottom')
p_bw_lines %>%
    ggsave(filename = 'manuscript/figures/ABCD_sibling_birth_weight_HAQER_background_ES-PGS.png',
           device = 'png', dpi = 300, bg = 'white', 
           units = 'in', width = 5, height = 4.5)


#######################################
## ABCD c-section discordance analysis
#######################################
## read in data
sib_csec_wd <- read_csv('manuscript/supplemental_materials/ABCD_sibling_csection_data.csv')

## paired t-test of HAQER CP-PGS (sibling born via c-section compared to sibling not born via c-section)
paired_csec <- broom::tidy(t.test(sib_csec_wd$non_csec_sib_cp_pgs.HAQER, sib_csec_wd$csec_sib_cp_pgs.HAQER, paired = TRUE)) %>% 
    mutate(n = nrow(sib_csec_wd)) %>% 
    mutate(pheno = 'comparing_sibling_pairs_where_only_one_born_via_csection') %>%
    mutate(x = 'cp_pgs.HAQER') %>%
    mutate(lab = str_c('t-statistic = ', round(abs(statistic), 2), ', p-val = ', formatC(p.value, digits = 2), '\nN=', parameter + 1, ' sibling pairs')) %>% 
    relocate(pheno, x)

p_paired_haq <- sib_csec_wd %>% 
    mutate(lab = paired_csec$lab) %>%
    pivot_longer(cols = c(non_csec_sib_cp_pgs.HAQER, csec_sib_cp_pgs.HAQER)) %>% 
    mutate(type = ifelse(name == 'non_csec_sib_cp_pgs.HAQER', 'Normal delivery', 'C-section')) %>% 
    mutate(type = factor(type, levels = c( 'Normal delivery', 'C-section'))) %>%
    ggplot(aes(x = type, y = value)) +
    geom_violin(size = 1.1) +
    geom_boxplot(width = 0.4, aes(fill = type), size = 1.1, alpha = 0.7) +
    geom_point() +
    geom_line(aes(group = FID), color = 'grey60', alpha = .4) +
    geom_hline(yintercept = 0, color = 'red', linetype = 'dashed', size = 1.1) +
    ylab('HAQER CP-PGS') +
    xlab('Birth delivery method') +
    geom_text(aes(x = 1.5, y = 2.7, label = lab), size = 4, check_overlap = TRUE) +
    scale_fill_manual(values = c('grey70', 'chocolate1')) +
    theme_classic() +
    theme(axis.text = element_text(size = 12),
          axis.title = element_text(size = 14),
          strip.text = element_text(size = 14),
          legend.text = element_text(size = 12),
          plot.title = element_text(size = 14, hjust = .5),
          legend.position = 'none')

## background HAQER CP-PGS
paired_csec_bg <- broom::tidy(t.test(sib_csec_wd$non_csec_sib_cp_pgs.background_HAQER, sib_csec_wd$csec_sib_cp_pgs.background_HAQER, paired = TRUE)) %>% 
    mutate(n = nrow(sib_csec_wd)) %>% 
    mutate(pheno = 'comparing_sibling_pairs_where_only_one_born_via_csection') %>%
    mutate(x = 'cp_pgs.background_HAQER') %>%
    mutate(lab = str_c('t-statistic = ', round(abs(statistic), 2), ', p-val = ', formatC(p.value, digits = 2), '\nN=', parameter + 1, ' sibling pairs')) %>% 
    relocate(pheno, x)

p_paired_bg <- sib_csec_wd %>% 
    mutate(lab = paired_csec_bg$lab) %>%
    pivot_longer(cols = c(non_csec_sib_cp_pgs.background_HAQER, csec_sib_cp_pgs.background_HAQER)) %>% 
    mutate(type = ifelse(name == 'non_csec_sib_cp_pgs.background_HAQER', 'Normal delivery', 'C-section')) %>% 
    mutate(type = factor(type, levels = c( 'Normal delivery', 'C-section'))) %>%
    ggplot(aes(x = type, y = scale(value)[,1])) +
    geom_violin(size = 1.1) +
    geom_boxplot(width = 0.4, aes(fill = type), size = 1.1, alpha = 0.7) +
    geom_point() +
    geom_line(aes(group = FID), color = 'grey60', alpha = .4) +
    geom_hline(yintercept = 0, color = 'red', linetype = 'dashed', size = 1.1) +
    xlab('Birth delivery method') +
    ylab('Background CP-PGS') +
    geom_text(aes(x = 1.5, y = 3, label = lab), size = 4, check_overlap = TRUE) +
    scale_fill_manual(values = c('grey70', 'chocolate1')) +
    theme_classic() +
    theme(axis.text = element_text(size = 12),
          axis.title = element_text(size = 14),
          strip.text = element_text(size = 14),
          legend.text = element_text(size = 12),
          plot.title = element_text(size = 14, hjust = .5),
          legend.position = 'none')

## save c-section figures
p_paired_haq %>% 
    ggsave(filename = 'manuscript/figures/ABCD_sibling_csection_HAQER_ES-PGS.png',
           device = 'png', dpi = 300, bg = 'white', 
           units = 'in', width = 5, height = 5)
p_paired_bg %>% 
    ggsave(filename = 'manuscript/figures/ABCD_sibling_csection_HAQER_background_ES-PGS.png',
           device = 'png', dpi = 300, bg = 'white', 
           units = 'in', width = 5, height = 5)

## merge stats results and save to file
bind_rows(abcd_icv_es_pgs_stats, abcd_icv_growth_es_pgs_stats, abcd_sib_bw_es_pgs_res) %>% 
    write_csv('manuscript/supplemental_materials/stats/ABCD_HAQER_ES-PGS_results.csv')

bind_rows(paired_csec, paired_csec_bg) %>% 
    select(-lab) %>%
    write_csv('manuscript/supplemental_materials/stats/ABCD_c-section_HAQER_ES-PGS_results.csv')


################
## save figure objects
p_abcd_icv %>% 
    write_rds('manuscript/figures/R_plot_objects/ABCD_HAQER-CP-PGS_ICV.rds')
p_abcd_icv_growth  %>% 
    write_rds('manuscript/figures/R_plot_objects/ABCD_HAQER-CP-PGS_ICV_growth.rds')  
p_sib_bw %>%
    write_rds('manuscript/figures/R_plot_objects/ABCD_HAQER-CP-PGS_sibling_birth_weight.rds')
p_sib_bw_bg %>%
    write_rds('manuscript/figures/R_plot_objects/ABCD_background-CP-PGS_sibling_birth_weight.rds')
p_bw_lines %>% 
    write_rds('manuscript/figures/R_plot_objects/ABCD_HAQER-CP-PGS_sibling_birth_weight_lines.rds')
p_paired_bg %>% 
    write_rds('manuscript/figures/R_plot_objects/ABCD_HAQER-CP-PGS_sibling_csection.rds')
p_paired_haq %>% 
    write_rds('manuscript/figures/R_plot_objects/ABCD_background-CP-PGS_sibling_csection.rds')
