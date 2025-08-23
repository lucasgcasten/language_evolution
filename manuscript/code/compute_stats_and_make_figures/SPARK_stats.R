library(tidyverse)

######################
## read in SPARK data
######################
spark_df <- read_csv('manuscript/supplemental_materials/SPARK_data.csv')

####################
## ES-PGS analysis
####################
## stats for self-reported language diagnosis count (language impairment + dyslexia + stutter + ...)
### get number of language related diagoses for everyone
spark_df$lang_dx_count <- spark_df$dx_lang_impairment + spark_df$dx_hearing_impairment + spark_df$dx_speech_impediment + spark_df$dx_speech_stutter + spark_df$dx_dyslexia
spark_df$psych_dx_count <- spark_df$dx_mdd + spark_df$dx_bipolar + spark_df$dx_anxiety + spark_df$dx_schizophrenia + spark_df$dx_ocd + spark_df$dx_tourettes + spark_df$dx_substance_abuse

n_self_report_lang_dx <- nrow(spark_df[spark_df$asd == FALSE & is.na(spark_df$lang_dx_count) == FALSE & is.na(spark_df$cp_pgs.background_HAQER_v2) == FALSE,])

### fit zero inflated regression models to predict number of language diagnoses
m_dx_null <- pscl::zeroinfl(lang_dx_count ~ cp_pgs.background_HAQER_v2 + cp_pgs.matched_control_HAQER_v2 + age_self_report_dx + as.factor(sex) + pc1 + pc2 + pc3 + pc4 + pc5, 
                            data = spark_df[spark_df$asd == FALSE,], 
                            dist = 'poisson')
m_dx_alt <- pscl::zeroinfl(lang_dx_count ~ cp_pgs.HAQER_v2 + cp_pgs.background_HAQER_v2 + cp_pgs.matched_control_HAQER_v2 + age_self_report_dx + as.factor(sex) + pc1 + pc2 + pc3 + pc4 + pc5, 
                           data = spark_df[spark_df$asd == FALSE,], 
                           dist = 'poisson')

m_dx_psych_null <- pscl::zeroinfl(psych_dx_count ~ cp_pgs.background_HAQER_v2 + cp_pgs.matched_control_HAQER_v2 + age_self_report_dx + as.factor(sex) + pc1 + pc2 + pc3 + pc4 + pc5, 
                            data = spark_df[spark_df$asd == FALSE,], 
                            dist = 'poisson')
m_dx_psych_alt <- pscl::zeroinfl(psych_dx_count ~ cp_pgs.HAQER_v2 + cp_pgs.background_HAQER_v2 + cp_pgs.matched_control_HAQER_v2 + age_self_report_dx + as.factor(sex) + pc1 + pc2 + pc3 + pc4 + pc5, 
                           data = spark_df[spark_df$asd == FALSE,], 
                           dist = 'poisson')

### gather model stats
mod_sum <- summary(m_dx_alt)
mod_sum_psych <- summary(m_dx_psych_alt)

count_results <- data.frame(
  component = "count",
  variable = rownames(mod_sum$coefficients$count),
  estimate = mod_sum$coefficients$count[,1],
  std_error = mod_sum$coefficients$count[,2],
  z_value = mod_sum$coefficients$count[,3],
  p_value = mod_sum$coefficients$count[,4],
  stringsAsFactors = FALSE
  )
zero_results <- data.frame(
  component = "zero",
  variable = rownames(mod_sum$coefficients$zero),
  estimate = mod_sum$coefficients$zero[,1],
  std_error = mod_sum$coefficients$zero[,2],
  z_value = mod_sum$coefficients$zero[,3],
  p_value = mod_sum$coefficients$zero[,4],
  stringsAsFactors = FALSE
  )


count_results_psych <- data.frame(
  component = "count",
  variable = rownames(mod_sum_psych$coefficients$count),
  estimate = mod_sum_psych$coefficients$count[,1],
  std_error = mod_sum_psych$coefficients$count[,2],
  z_value = mod_sum_psych$coefficients$count[,3],
  p_value = mod_sum_psych$coefficients$count[,4],
  stringsAsFactors = FALSE
  )
zero_results_psych <- data.frame(
  component = "zero",
  variable = rownames(mod_sum_psych$coefficients$zero),
  estimate = mod_sum_psych$coefficients$zero[,1],
  std_error = mod_sum_psych$coefficients$zero[,2],
  z_value = mod_sum_psych$coefficients$zero[,3],
  p_value = mod_sum_psych$coefficients$zero[,4],
  stringsAsFactors = FALSE
  )

### save stats
lang_res <- broom::tidy(lmtest::lrtest(m_dx_null, m_dx_alt)) %>% 
    mutate(group = 'psychiatric')
psych_res <-  broom::tidy(lmtest::lrtest(m_dx_psych_null, m_dx_psych_alt)) %>% 
    mutate(group = 'psychiatric')

bind_rows(lang_res, psych_res) %>%
    drop_na() %>%
    select(-df) %>%
    rename(df = X.Df) %>% 
    rename(model = term) %>% 
    mutate(added_x = 'cp_pgs.HAQER_v2') %>% 
    write_csv("manuscript/supplemental_materials/stats/SPARK_ES-PGS_HAQER_validations_self_reported_language_diagnosis_LRT.csv")

lang_coef_res <- bind_rows(zero_results, count_results) %>% 
    mutate(y = 'language_diagnosis_count_self_report')
psych_coef_res <- bind_rows(zero_results_psych, count_results_psych) %>% 
    mutate(y = 'psychiatric_diagnosis_count_self_report')

bind_rows(lang_coef_res, psych_coef_res) %>% 
    as_tibble() %>% 
    filter(variable == 'cp_pgs.HAQER_v2') %>% 
    rename(x = variable) %>%
    relocate(y, x, component) %>% 
    mutate(n = n_self_report_lang_dx) %>% 
    write_csv('manuscript/supplemental_materials/stats/SPARK_ES-PGS_HAQER_validations_self_reported_language_diagnosis.csv')

## SCQ stats for HAQERs
n_mod <- spark_df %>% 
    select(IID, SCQ_final_score, SCQ_age_at_eval_years, SCQ_sex, matches("^SCQ_q"), cp_pgs.HAQER_v2, cp_pgs.background_HAQER, cp_pgs.matched_control_HAQER_v2, str_c('pc', 1:5)) %>% 
    drop_na(SCQ_final_score, cp_pgs.HAQER_v2) %>% 
    pivot_longer(cols = matches('SCQ_q')) %>% 
    drop_na(value) %>%
    group_by(name, value) %>% 
    count() %>% 
    ungroup() %>%
    mutate(value = str_c(as.logical(value), '_n')) %>% 
    pivot_wider(names_from = 'value', values_from = 'n')
rep_scq_stat_es_pgs <- spark_df %>% 
    select(IID, SCQ_final_score, SCQ_age_at_eval_years, SCQ_sex, SCQ_q01_phrases, matches("HAQER_v2"), matches('pc')) %>% 
    drop_na(SCQ_final_score, cp_pgs.HAQER_v2) %>% 
    pivot_longer(cols = matches('^SCQ_q')) %>%
    group_by(name) %>% 
    do(res = broom::tidy(glm(value ~ cp_pgs.HAQER_v2 + cp_pgs.matched_control_HAQER_v2 + cp_pgs.background_HAQER_v2 + SCQ_age_at_eval_years + as.factor(SCQ_sex) + pc1 + pc2 + pc3 + pc4 + pc5, data = ., family = 'binomial'))) %>% 
    unnest(res) %>% 
    filter(term == 'cp_pgs.HAQER_v2') %>% 
    mutate(pheno = str_c(name),
           x = 'cp_pgs.HAQER_v2'
           ) %>% 
    inner_join(n_mod) %>%
    relocate(pheno, x) %>%  
    select(-c(term, name)) %>% 
    mutate(n = FALSE_n + TRUE_n)

## IQ ES-PGS stats
iq_cnt <- spark_df %>%
    filter(asd == TRUE) %>%
    # filter(age_years >= 4) %>%
    pivot_longer(cols = matches('iq_score')) %>% 
    drop_na(name, value, cp_pgs.HAQER_v2, cp_pgs.background_HAQER, sex, age_years) %>% 
    group_by(name) %>% 
    count() %>% 
    ungroup()
spark_iq_es_pgs_res <- spark_df %>%
    drop_na(cp_pgs.HAQER_v2, cp_pgs.background_HAQER, sex, age_years) %>% 
    pivot_longer(cols = matches('iq_score')) %>% 
    group_by(name) %>% 
    do(res = broom::tidy(lm(value ~ cp_pgs.HAQER_v2 + cp_pgs.matched_control_HAQER_v2 + cp_pgs.background_HAQER + as.factor(sex) + age_years + pc1 + pc2 + pc3 + pc4 + pc5, data = .))) %>% 
    unnest(res) %>% 
    filter(term == 'cp_pgs.HAQER_v2') %>% 
    arrange(p.value) %>% 
    mutate(x = 'cp_pgs.HAQER_v2') %>% 
    inner_join(iq_cnt) %>%
    rename(pheno = name) %>%
    relocate(pheno, x) %>% 
    select(-term)

spark_iq_es_pgs_res_lab <- spark_iq_es_pgs_res %>%  
    mutate(lab = str_c('ES-PGS beta = ', round(estimate, digits = 2), ', p-val = ', formatC(p.value, digits = 2), '\nN = ', prettyNum(n, big.mark = ',')))

spark_p_iq <- spark_df %>%
    filter(asd == TRUE) %>%
    # filter(age_years >= 4) %>%
    pivot_longer(cols = matches('iq_score')) %>% 
    rename(pheno = name) %>%
    inner_join(spark_iq_es_pgs_res_lab) %>% 
    drop_na(cp_pgs.HAQER_v2) %>%
    mutate(pheno_clean = case_when(pheno == 'viq_score' ~ 'Verbal IQ',
                                   pheno == 'nviq_score' ~ 'Nonverbal IQ',
                                   pheno == 'fsiq_score' ~ 'Full-scale IQ'),
           pheno_clean = factor(pheno_clean, levels = c('Verbal IQ', 'Nonverbal IQ', 'Full-scale IQ'))) %>%
    mutate(sig = ifelse(p.value < 0.05, TRUE, FALSE)) %>%
    ggplot(aes(x = value, y = cp_pgs.HAQER_v2)) +
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
bind_rows(rep_scq_stat_es_pgs, spark_iq_es_pgs_res) %>% 
    rename(beta = estimate) %>% 
    write_csv('manuscript/supplemental_materials/stats/SPARK_ES-PGS_HAQER_validations.csv')

############################
## HAQER reversion analysis
############################
## diagnosis results
dx_cnt_rev <- spark_df %>%
    filter(asd == TRUE) %>%
    filter(age_years >= 3) %>%
    select(-c(199:ncol(spark_df))) %>% ## drop self reported language diagnosis
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
    filter(age_years >= 3) %>%
    drop_na(haqer_rare_variant_reversion_count, sex, age_years) %>%
    mutate(sex_female = ifelse(sex == 'Female', 1, 0)) %>%
    pivot_longer(cols = matches('dx_dev_lang|dx_dev_id')) %>% 
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

p_rev_dist_h <- p_dist_haq | p_dist_har | p_dist_rand # + plot_annotation(tag_levels = 'A') & theme(plot.tag = element_text(face = 'bold', size = 18))

## forest plot of regression betas for language phenos
cl <- c("#762776", "#e04468", "#dcc699")

p_rev_forest <- bind_rows(reversions_spark_dx_res, reversions_dev_res) %>% 
    mutate(lab = str_c('N ASD cases = ', n)) %>% 
    mutate(pheno_clean = case_when(pheno == 'dx_dev_lang_disorder' ~ str_c('Dev. language disorder\n(', prettyNum(n_case, big.mark = ','), ' cases)'),
                                   pheno == 'dx_intellectual_disability' ~ str_c('Intellectual disability\n(', prettyNum(n_case, big.mark = ','), ' cases)'),
                                   pheno == 'walked_age_mos' ~ 'Age started walking',
                                   pheno == 'used_words_age_mos' ~ 'Age of first word',
                                   pheno == 'combined_words_age_mos' ~ 'Age combined words',
                                   pheno == 'combined_phrases_age_mos' ~ 'Age combined phrases'),
           x_clean = case_when(str_detect(x, 'haqer_') ~ 'HAQER rare reversions',
                               str_detect(x, 'har_') ~ 'HAR rare reversions',
                               str_detect(x, 'rand_') ~ 'RAND rare reversions')) %>% 
    mutate(type = case_when(str_detect(pheno, 'dx_') ~ str_c('Diagnosis (N = ', prettyNum(n, big.mark = ','), ')'),
                            str_detect(pheno, 'dx_', negate = TRUE) ~ str_c('Developmental milestones (N = ', prettyNum(n, big.mark = ','), ')'))
           ) %>%
    mutate(pheno = factor(pheno, levels = rev(c('dx_dev_lang_disorder', 'dx_intellectual_disability', 'used_words_age_mos', 'combined_words_age_mos', 'combined_phrases_age_mos', 'walked_age_mos')))) %>%
    arrange(pheno) %>% 
    mutate(pheno_clean = factor(pheno_clean, levels = unique(pheno_clean))) %>%
    mutate(x_clean = factor(x_clean, levels = rev(c('HAQER rare reversions', 'HAR rare reversions', 'RAND rare reversions')))) %>%
    arrange(type) %>% 
    mutate(type = factor(type, levels = rev(unique(type)))) %>%
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


p_rev_forest_dx <- reversions_spark_dx_res %>% 
    mutate(lab = str_c('N ASD cases = ', n)) %>% 
    mutate(pheno_clean = case_when(pheno == 'dx_dev_lang_disorder' ~ str_c('Dev. language disorder\n(', prettyNum(n_case, big.mark = ','), ' cases)'),
                                   pheno == 'dx_intellectual_disability' ~ str_c('Intellectual disability\n(', prettyNum(n_case, big.mark = ','), ' cases)'),
                                   pheno == 'walked_age_mos' ~ 'Age started walking',
                                   pheno == 'used_words_age_mos' ~ 'Age of first word',
                                   pheno == 'combined_words_age_mos' ~ 'Age combined words',
                                   pheno == 'combined_phrases_age_mos' ~ 'Age combined phrases'),
           x_clean = case_when(str_detect(x, 'haqer_') ~ 'HAQER rare reversions',
                               str_detect(x, 'har_') ~ 'HAR rare reversions',
                               str_detect(x, 'rand_') ~ 'RAND rare reversions')) %>% 
    mutate(type = case_when(str_detect(pheno, 'dx_') ~ str_c('Diagnosis (N = ', prettyNum(n, big.mark = ','), ')'),
                            str_detect(pheno, 'dx_', negate = TRUE) ~ str_c('Developmental milestones (N = ', prettyNum(n, big.mark = ','), ')'))
           ) %>%
    mutate(pheno = factor(pheno, levels = rev(c('dx_dev_lang_disorder', 'dx_intellectual_disability', 'used_words_age_mos', 'combined_words_age_mos', 'combined_phrases_age_mos', 'walked_age_mos')))) %>%
    arrange(pheno) %>% 
    mutate(pheno_clean = factor(pheno_clean, levels = unique(pheno_clean))) %>%
    mutate(x_clean = factor(x_clean, levels = rev(c('HAQER rare reversions', 'HAR rare reversions', 'RAND rare reversions')))) %>%
    arrange(type) %>% 
    mutate(type = factor(type, levels = rev(unique(type)))) %>%
    mutate(sig = ifelse(p.value < 0.05, TRUE, FALSE)) %>%
    mutate(type = 'SPARK parent reported diagnosis') %>%
    ggplot(aes(x = beta, y = pheno_clean, color = x_clean)) +
    geom_linerange(aes(xmin = beta - 1.96 * std.error, xmax = beta + 1.96 * std.error), size = 1.5, position = position_dodge(.45)) +
    geom_point(size = 5, position = position_dodge(.45), aes(shape = sig)) +
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
    xlab('Rare reversion beta (95% CI)') +
    guides(shape = 'none')


########################
## save plot objects
spark_p_iq %>% 
    write_rds('manuscript/figures/R_plot_objects/SPARK_HAQER_IQ_replication.rds')
p_rev_dist %>% 
    write_rds('manuscript/figures/R_plot_objects/SPARK_rare_reversion_histograms.rds')
p_rev_dist_h %>% 
    write_rds('manuscript/figures/R_plot_objects/SPARK_rare_reversion_histograms_horizontal.rds')
p_rev_forest %>% 
    write_rds('manuscript/figures/R_plot_objects/SPARK_rare_reversion_forest.rds')
p_rev_forest_dx %>% 
    write_rds('manuscript/figures/R_plot_objects/SPARK_rare_reversion_diagnosis_forest.rds')