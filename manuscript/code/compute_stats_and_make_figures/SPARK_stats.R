library(tidyverse)

######################
## read in SPARK data
######################
spark_df <- read_csv('manuscript/supplemental_materials/SPARK_data.csv')

####################
## ES-PGS analysis
####################
## sent rep factor stats
n_mod <- spark_df %>% 
    select(IID, factor_sent_rep, age_years, sex, cp_pgs.HAQER, cp_pgs.background_HAQER) %>% 
    drop_na() %>% 
    nrow(.)
rep_fac_stat_es_pgs <- broom::tidy(lm(factor_sent_rep ~ cp_pgs.HAQER + cp_pgs.background_HAQER + age_years + as.factor(sex), data = spark_df)) %>% 
    filter(term == 'cp_pgs.HAQER') %>% 
    mutate(pheno = 'SPARK_core_language_factor',
           x = 'cp_pgs.HAQER',
           n = n_mod) %>% 
    relocate(pheno, x) %>% 
    select(-term)

## make scatterplot
tmp <- spark_df %>% 
    select(IID, factor_sent_rep, age_years, sex, cp_pgs.HAQER, cp_pgs.background_HAQER) %>% 
    drop_na() %>% 
    mutate(factor_sent_rep_resid = resid(lm(factor_sent_rep ~ age_years + as.factor(sex) + cp_pgs.background_HAQER)),
           factor_sent_rep_resid = scale(factor_sent_rep_resid)[,1]) %>% 
    relocate(factor_sent_rep_resid, .after = factor_sent_rep)
rep_fac_stat_cor <- broom::tidy(cor.test(tmp$factor_sent_rep_resid, tmp$cp_pgs.HAQER)) %>% 
    mutate(pheno = 'SPARK_core_language_factor',
           x = 'cp_pgs.HAQER',
           n = parameter + 2) %>% 
    mutate(lab = str_c('r = ', round(estimate, digits = 2), ', p-val = ', formatC(p.value, digits = 2), '\nN = ', n))

p_f1_rep <- tmp %>% 
    ggplot(aes(x = factor_sent_rep_resid, y = cp_pgs.HAQER)) +
    geom_point(size = 0.7) +
    geom_smooth(method = 'lm', size = 1.5, color = 'black') +
    xlab('SPARK core language score\nadjusted for background CP-PGS') +
    ylab('HAQER CP-PGS') +
    geom_text(aes(x = -1, y = 3.2, label = rep_fac_stat_cor$lab), size = 4, check_overlap = TRUE) +
    theme_classic() +  
    theme(axis.text = element_text(size = 12),
          axis.title = element_text(size = 14),
          legend.text = element_text(size = 12),
          legend.title = element_text(size = 14))
p_f1_rep %>%
    ggsave(filename = 'manuscript/figures/SPARK_replication_core_language_HAQER_ES-PGS.png', 
           device = 'png', dpi = 300, bg = 'white', 
           units = 'in', width = 5, height = 5)

###################
## ES-PGS dx stats
dx_cnt <- spark_df %>%
    filter(asd == TRUE) %>%
    drop_na(cp_pgs.HAQER) %>%
    pivot_longer(cols = matches('dx_')) %>% 
    group_by(name, value) %>% 
    count() %>% 
    drop_na() %>%
    mutate(type = ifelse(value == 0, 'n_control', 'n_case')) %>% 
    pivot_wider(id_cols = name, names_from = type, values_from = n) %>% 
    rename(pheno = name) %>% 
    mutate(n = n_control + n_case) %>% 
    relocate(n, .before = n_control) %>% 
    ungroup()
dx_mean = spark_df %>%
    filter(asd == TRUE) %>%
    drop_na(cp_pgs.HAQER) %>%
    pivot_longer(cols = matches('dx_')) %>% 
    group_by(name, value) %>% 
    summarise(mean_val = mean(cp_pgs.HAQER)) %>% 
    drop_na() %>%
    mutate(type = ifelse(value == 0, 'control_mean_cp_pgs.HAQER', 'case_mean_cp_pgs.HAQER')) %>% 
    pivot_wider(id_cols = name, names_from = type, values_from = mean_val) %>% 
    rename(pheno = name) %>% 
    ungroup()
dx_sd = spark_df %>%
    filter(asd == TRUE) %>%
    drop_na(cp_pgs.HAQER) %>%
    pivot_longer(cols = matches('dx_')) %>% 
    group_by(name, value) %>% 
    summarise(sd_val = sd(cp_pgs.HAQER)) %>% 
    drop_na() %>%
    mutate(type = ifelse(value == 0, 'control_std_dev_cp_pgs.HAQER', 'case_std_dev_cp_pgs.HAQER')) %>% 
    pivot_wider(id_cols = name, names_from = type, values_from = sd_val) %>% 
    rename(pheno = name) %>% 
    ungroup()
spark_dx_es_pgs_res <- spark_df %>%
    filter(asd == TRUE) %>%
    drop_na(cp_pgs.HAQER) %>%
    pivot_longer(cols = matches('dx_')) %>% 
    group_by(name) %>% 
    do(res = broom::tidy(glm(value ~ cp_pgs.HAQER + cp_pgs.background_HAQER + age_at_eval_years + as.factor(sex), family = 'binomial', data = .))) %>% 
    unnest(res) %>% 
    filter(term == 'cp_pgs.HAQER') %>% 
    arrange(p.value) %>% 
    rename(pheno = name) %>%
    mutate(x = 'cp_pgs.HAQER') %>% 
    relocate(pheno, x) %>% 
    inner_join(dx_cnt) %>%
    inner_join(dx_mean) %>%
    inner_join(dx_sd) %>%
    select(-term) %>% 
    filter(pheno %in% c('dx_dev_lang', 'dx_dev_id')) %>% 
    mutate(pheno = case_when(pheno == 'dx_dev_lang' ~ 'dx_dev_lang_disorder',
                             pheno == 'dx_dev_id' ~ 'dx_intellectual_disability',
                             TRUE ~ pheno))

## IQ ES-PGS stats
iq_cnt <- spark_df %>%
    filter(asd == TRUE) %>%
    filter(age_years >= 4) %>%
    pivot_longer(cols = matches('iq_score')) %>% 
    drop_na(name, value, cp_pgs.HAQER, cp_pgs.background_HAQER, sex, age_years) %>% 
    group_by(name) %>% 
    count() %>% 
    ungroup()
spark_iq_es_pgs_res <- spark_df %>%
    filter(asd == TRUE) %>%
    filter(age_years >= 4) %>%
    drop_na(cp_pgs.HAQER, cp_pgs.background_HAQER, sex, age_years) %>% 
    pivot_longer(cols = matches('iq_score')) %>% 
    group_by(name) %>% 
    do(res = broom::tidy(lm(value ~ cp_pgs.HAQER + cp_pgs.background_HAQER + as.factor(sex) + age_years, data = .))) %>% 
    unnest(res) %>% 
    filter(term == 'cp_pgs.HAQER') %>% 
    arrange(p.value) %>% 
    mutate(x = 'cp_pgs.HAQER') %>% 
    inner_join(iq_cnt) %>%
    rename(pheno = name) %>%
    relocate(pheno, x) %>% 
    select(-term)

spark_iq_es_pgs_res_lab <- spark_iq_es_pgs_res %>%  
    mutate(lab = str_c('ES-PGS beta = ', round(estimate, digits = 2), ', p-val = ', formatC(p.value, digits = 2), '\nN = ', prettyNum(n, big.mark = ',')))

spark_p_iq <- spark_df %>%
    filter(asd == TRUE) %>%
    filter(age_years >= 4) %>%
    pivot_longer(cols = matches('iq_score')) %>% 
    rename(pheno = name) %>%
    inner_join(spark_iq_es_pgs_res_lab) %>% 
    drop_na(cp_pgs.HAQER) %>%
    mutate(pheno_clean = case_when(pheno == 'viq_score' ~ 'Verbal IQ',
                                   pheno == 'nviq_score' ~ 'Nonverbal IQ',
                                   pheno == 'fsiq_score' ~ 'Full-scale IQ'),
           pheno_clean = factor(pheno_clean, levels = c('Verbal IQ', 'Nonverbal IQ', 'Full-scale IQ'))) %>%
    mutate(sig = ifelse(p.value < 0.05, TRUE, FALSE)) %>%
    ggplot(aes(x = value, y = cp_pgs.HAQER)) +
    geom_jitter(size = 0.3, alpha = 0.7, aes(color = sig)) +
    geom_smooth(method = 'lm', size = 1.5, aes(color = sig)) +
    facet_wrap(~ pheno_clean) +
    geom_text(aes(x = 85, y = 3.75, label = lab), check_overlap = TRUE, size = 3) +
    theme_classic() +  
    theme(axis.text = element_text(size = 12),
          axis.title = element_text(size = 14),
          legend.text = element_text(size = 12),
          legend.title = element_text(size = 14),
          strip.text = element_text(size = 14),
          legend.position = 'none') +
    scale_color_manual(values = c('grey75', 'black')) +
    xlab('IQ score') +
    ylab('HAQER CP-PGS')
spark_p_iq %>%
    ggsave(filename = 'manuscript/figures/SPARK_replication_IQ_HAQER_ES-PGS.png', 
           device = 'png', dpi = 300, bg = 'white', 
           units = 'in', width = 8, height = 3.5)

## merge ES-PGS sumstats to a table and save
bind_rows(rep_fac_stat_es_pgs, spark_dx_es_pgs_res, spark_iq_es_pgs_res) %>% 
    rename(beta = estimate) %>% 
    write_csv('manuscript/supplemental_materials/stats/SPARK_ES-PGS_HAQER_validations.csv')

############################
## HAQER reversion analysis
############################
## diagnosis results
dx_cnt_rev <- spark_df %>%
    filter(asd == TRUE) %>%
    drop_na(haqer_rare_variant_reversion_count, matches('dx')) %>%
    pivot_longer(cols = matches('dx_')) %>% 
    group_by(name, value) %>% 
    count() %>% 
    drop_na() %>%
    mutate(type = ifelse(value == 0, 'n_control', 'n_case')) %>% 
    pivot_wider(id_cols = name, names_from = type, values_from = n) %>% 
    rename(pheno = name) %>% 
    mutate(n = n_control + n_case) %>% 
    relocate(n, .before = n_control) %>% 
    ungroup()

reversions_spark_dx_res <- spark_df %>%
    filter(asd == TRUE) %>%
    drop_na(haqer_rare_variant_reversion_count, sex, age_years) %>%
    mutate(sex_female = ifelse(sex == 'Female', 1, 0)) %>%
    pivot_longer(cols = matches('dx_')) %>% 
    group_by(name) %>% 
    do(res = broom::tidy(glm(value ~ scale(haqer_rare_variant_reversion_count)[,1] + scale(har_rare_variant_reversion_count)[,1] + scale(rand_rare_variant_reversion_count)[,1] + scale(age_years)[,1] + scale(sex_female)[,1], 
                             family = 'binomial', data = .))) %>% 
    unnest(res) %>% 
    filter(str_detect(term, pattern = 'rare_variant_reversion')) %>% 
    rename(pheno = name, beta = estimate, x = term) %>% 
    mutate(x = str_split(x, pattern = '[(]', simplify = TRUE)[,2],
           x = str_split(x, pattern = '[)]', simplify = TRUE)[,1]) %>%
    arrange(p.value) %>% 
    relocate(pheno, x) %>% 
    inner_join(dx_cnt_rev) %>%
    filter(pheno %in% c('dx_dev_lang', 'dx_dev_id')) %>% 
    mutate(pheno = case_when(pheno == 'dx_dev_lang' ~ 'dx_dev_lang_disorder',
                             pheno == 'dx_dev_id' ~ 'dx_intellectual_disability',
                             TRUE ~ pheno))


## developmental milestones results
n_dev_rev <- spark_df %>% 
    filter(asd == TRUE) %>%
    filter(age_years >= 5) %>%
    select(IID, used_words_age_mos, combined_words_age_mos, combined_phrases_age_mos, walked_age_mos, haqer_rare_variant_reversion_count, har_rare_variant_reversion_count, rand_rare_variant_reversion_count, age_years, sex) %>%
    drop_na() %>% 
    nrow(.) 

reversions_dev_res <- spark_df %>% 
    filter(asd == TRUE) %>%
    filter(age_years >= 5) %>%
    select(IID, used_words_age_mos, combined_words_age_mos, combined_phrases_age_mos, walked_age_mos, haqer_rare_variant_reversion_count, har_rare_variant_reversion_count, rand_rare_variant_reversion_count, age_years, sex) %>%
    drop_na() %>%
    pivot_longer(cols = matches('_mos')) %>% 
    group_by(name) %>% 
    mutate(normed_val = qnorm((rank(value,na.last="keep")-0.5)/sum(!is.na(value)))) %>% ## normalize trait distribution for linear regression
    mutate(sex_female = ifelse(sex == 'Female', 1, 0)) %>%
    do(res = broom::tidy(lm(normed_val ~ scale(haqer_rare_variant_reversion_count)[,1] + scale(har_rare_variant_reversion_count)[,1] + scale(rand_rare_variant_reversion_count)[,1] + scale(age_years)[,1] + scale(sex_female)[,1], data = .))) %>% 
    unnest(res) %>%
    filter(str_detect(term, pattern = 'rare_variant_reversion')) %>% 
    rename(pheno = name, beta = estimate, x = term) %>% 
    mutate(n = n_dev_rev,
           x = str_split(x, pattern = '[(]', simplify = TRUE)[,2],
           x = str_split(x, pattern = '[)]', simplify = TRUE)[,1])

## save rare reversion sumstats to file
bind_rows(reversions_spark_dx_res, reversions_dev_res) %>% 
    write_csv('manuscript/supplemental_materials/stats/SPARK_rare_reversion_results.csv')

#########################
## make reversion figures
#########################
## histograms of reversion counts
hist_dat <- spark_df %>% 
    drop_na(haqer_rare_variant_reversion_count) %>%
    pivot_longer(cols = matches('reversion')) %>%
    group_by(name, value) %>% 
    count() 

p_dist_haq <- hist_dat[hist_dat$name == 'haqer_rare_variant_reversion_count',] %>%
    ggplot(aes(x =value, y = n, fill = name)) +
    geom_bar(stat = 'identity', position = 'dodge') +
    coord_cartesian(xlim = c(0, max(hist_dat$value)), ylim = c(0, max(hist_dat$n))) +
    xlab('HAQER rare reversions') +
    ylab('Count') +
    scale_fill_manual(values = c('#762776')) +
    theme_classic() +
    theme(axis.text = element_text(size = 14),
          axis.title = element_text(size = 16),
          legend.position = 'none')

p_dist_har <- hist_dat[hist_dat$name == 'har_rare_variant_reversion_count',] %>%
    ggplot(aes(x =value, y = n, fill = name)) +
    geom_bar(stat = 'identity', position = 'dodge') +
    coord_cartesian(xlim = c(0, max(hist_dat$value)), ylim = c(0, max(hist_dat$n))) +
    xlab('HAR rare reversions') +
    ylab('Count') +
    scale_fill_manual(values = c('#e04468')) +
    theme_classic() +
    theme(axis.text = element_text(size = 14),
          axis.title = element_text(size = 16),
          legend.position = 'none')

p_dist_rand <- hist_dat[hist_dat$name == 'rand_rare_variant_reversion_count',] %>%
    ggplot(aes(x =value, y = n, fill = name)) +
    geom_bar(stat = 'identity', position = 'dodge') +
    coord_cartesian(xlim = c(0, max(hist_dat$value)), ylim = c(0, max(hist_dat$n))) +
    xlab('RAND rare reversions') +
    ylab('Count') +
    scale_fill_manual(values = c('#dcc699')) +
    theme_classic() +
    theme(axis.text = element_text(size = 14),
          axis.title = element_text(size = 16),
          legend.position = 'none')

library(patchwork)
p_rev_dist <- p_dist_haq / p_dist_har / p_dist_rand # + plot_annotation(tag_levels = 'A') & theme(plot.tag = element_text(face = 'bold', size = 18))
ggsave(p_rev_dist, 
       filename = 'manuscript/figures/SPARK_rare_reversion_distributions.png', 
       dpi = 300, device = 'png', 
       units = 'in', width = 4.5, height = 8)

## forest plot of regression betas for language phenos
cl <- c("#762776", "#e04468", "#dcc699")

p_rev_forest <- bind_rows(reversions_spark_dx_res, reversions_dev_res) %>% 
    mutate(lab = str_c('N ASD cases = ', n)) %>% 
    mutate(pheno_clean = case_when(pheno == 'dx_dev_lang_disorder' ~ str_c('Developmental language disorder\n(', prettyNum(n_case, big.mark = ','), ' cases)'),
                                   pheno == 'dx_intellectual_disability' ~ str_c('Intellectual disability\n(', prettyNum(n_case, big.mark = ','), ' cases)'),
                                   pheno == 'walked_age_mos' ~ 'Age started walking',
                                   pheno == 'used_words_age_mos' ~ 'Age of first word',
                                   pheno == 'combined_words_age_mos' ~ 'Age combined words',
                                   pheno == 'combined_phrases_age_mos' ~ 'Age combined phrases'),
           x_clean = case_when(str_detect(x, 'haqer_') ~ 'HAQER rare reversions',
                               str_detect(x, 'har_') ~ 'HAR rare reversions',
                               str_detect(x, 'rand_') ~ 'RAND rare reversions')) %>% 
    mutate(type = case_when(str_detect(pheno, 'dx_') ~ str_c('Diagnosis (N = ', n, ')'),
                            str_detect(pheno, 'dx_', negate = TRUE) ~ str_c('Developmental milestones (N = ', prettyNum(n, big.mark = ','), ')'))
           ) %>%
    mutate(pheno = factor(pheno, levels = rev(c('dx_dev_lang_disorder', 'dx_intellectual_disability', 'used_words_age_mos', 'combined_words_age_mos', 'combined_phrases_age_mos', 'walked_age_mos')))) %>%
    arrange(pheno) %>% 
    mutate(pheno_clean = factor(pheno_clean, levels = unique(pheno_clean))) %>%
    mutate(x_clean = factor(x_clean, levels = rev(c('HAQER rare reversions', 'HAR rare reversions', 'RAND rare reversions')))) %>%
    mutate(sig = ifelse(p.value < 0.05, TRUE, FALSE)) %>%
    ggplot(aes(x = beta, y = pheno_clean, color = x_clean)) +
    geom_linerange(aes(xmin = beta - 1.96 * std.error, xmax = beta + 1.96 * std.error), size = 1.1, position = position_dodge(.45)) +
    geom_point(size = 3.5, position = position_dodge(.45), aes(shape = sig)) +
    scale_shape_manual(values = c(1,16)) +
    geom_vline(xintercept = 0, color = 'red', linetype = 'dashed', size = 1.075) +
    ggforce::facet_col(facets = vars(type), 
                       scales = "free", 
                       space = "free") +
    labs(color = NULL) +
    theme_classic() +
    theme(axis.text = element_text(size = 10),
          axis.title = element_text(size = 10),
          strip.text = element_text(size = 12),
          legend.text = element_text(size = 10),
          legend.position = 'bottom') +
    scale_color_manual(values = rev(cl)) +
    ylab(NULL) +
    xlab('Regression beta (95% CI)') +
    guides(shape = 'none')

p_rev_forest %>% 
    ggsave(filename = 'manuscript/figures/SPARK_rare_reversion_associations.png', 
           dpi = 300, device = 'png', 
           units = 'in', width = 8.2, height = 6)


########################
## save plot objects
p_f1_rep %>% 
    write_rds('manuscript/figures/R_plot_objects/SPARK_HAQER_F1_replication.rds')
spark_p_iq %>% 
    write_rds('manuscript/figures/R_plot_objects/SPARK_HAQER_IQ_replication.rds')
p_rev_dist  %>% 
    write_rds('manuscript/figures/R_plot_objects/SPARK_rare_reversion_histograms.rds')
p_rev_forest %>% 
    write_rds('manuscript/figures/R_plot_objects/SPARK_rare_reversion_forest.rds')
