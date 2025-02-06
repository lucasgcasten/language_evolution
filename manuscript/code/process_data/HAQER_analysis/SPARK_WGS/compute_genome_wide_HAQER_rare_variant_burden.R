library(tidyverse)

#####################################
## read in data
#####################################
## phenos
ph <- read_csv('/wdata/lcasten/spark/WGS/data/WGS_gathered_phenotypes.csv')

## read in VCF like DF to get samples and score each individual
tmp_list <- list()
for(i in 1:22) {
    ####
    cat('\n\n\n\n')
    message('==========================================')
    message('Chromosome: ', i)
    message('==========================================')
    ####
    tmp <- read_rds(str_c('/wdata/lcasten/sli_wgs/HAQER_variant_associations/data/SPARK_WGS/chr', i, '_HAQER_10Kb_flank.hg38.rds')) %>% 
        mutate(max_gnomAD_af = ifelse(is.na(max_gnomAD_af), 0, max_gnomAD_af),
               thousand_genomes = ifelse(is.na(thousand_genomes), 0, thousand_genomes))

    tmp_brd <- data.frame(sample = colnames(tmp)[str_detect(colnames(tmp), 'SP')]) %>% 
        as_tibble()

    tmp_rare <- tmp %>% 
        filter(af < 0.01 & max_gnomAD_af < 0.01 & thousand_genomes < 0.01) %>% 
        select(matches('SP'))
    
    tmp_brd$haqer_rare_variant_count = colSums(tmp_rare)
    tmp_rare <- tmp  %>% 
        filter(af < 0.01 & max_gnomAD_af < 0.01 & thousand_genomes < 0.01 & alt == HCA_allele & ref == neanderthal_allele) %>% 
        select(matches('SP'))
    tmp_brd$haqer_rare_variant_reversion_count = colSums(tmp_rare)
    tmp_list[[i]] <- tmp_brd
}

brd <- bind_rows(tmp_list) %>% 
    rename(spid = sample) %>%
    group_by(spid) %>% 
    summarise(haqer_rare_variant_count = sum(haqer_rare_variant_count),
              haqer_rare_variant_reversion_count = sum(haqer_rare_variant_reversion_count))
hist(brd$haqer_rare_variant_reversion_count)
brd

brd %>% 
    write_csv('/wdata/lcasten/sli_wgs/HAQER_variant_associations/data/SPARK_WGS/HAQER_10Kb_flank_rare_variant_burden_counts.csv')

brd2 <- brd  %>% 
    inner_join(pc) %>% 
    mutate(resid_haqer_rare_variant_reversion_count = scale(resid(lm(haqer_rare_variant_reversion_count ~ pc1 + pc2 + pc3 + pc4 + pc5)))[,1],
           resid_haqer_rare_variant_count = scale(resid(lm(haqer_rare_variant_count ~ pc1 + pc2 + pc3 + pc4 + pc5)))[,1])

#################
pc <- read_csv('/wdata/lcasten/spark/prs/HapMap3_plus/PCA/raw_KING_pca_results.csv') %>% 
    rename(spid = sample.ID) %>% 
    select(1:6)

##
tmp <- ph %>% 
    # inner_join(pc) %>%
    inner_join(brd) %>% 
    mutate(asd = as.numeric(asd))

tmp2 <- ph %>% 
    inner_join(brd2) %>% 
    mutate(asd = as.numeric(asd))

############################################################################################
## remove outliers from analysis (just to be safe, saw the same associations even when they're included)
############################################################################################
md <- mad(brd$haqer_rare_variant_reversion_count)
me <- median(brd$haqer_rare_variant_reversion_count)

kp <- brd %>% 
    filter(haqer_rare_variant_reversion_count <= me + (3 * md) & haqer_rare_variant_reversion_count >= me - (3 * md))
names(ph)


#############################################################
## language development milestones
#############################################################
hist(ph$used_words_age_mos)
hist(ph$combined_words_age_mos)
hist(ph$combined_phrases_age_mos)

word_age_res <- ph %>% 
    select(spid, matches('words|phrase')) %>% 
    filter(used_words_age_mos < 888 & combined_words_age_mos < 888 & combined_phrases_age_mos < 888) %>%
    drop_na() %>%
    pivot_longer(cols = where(is.numeric)) %>% 
    filter(str_detect(name, 'words|phrase')) %>% 
    drop_na(value) %>%
    filter(value < 888) %>%
    drop_na(value) %>%
    inner_join(kp) %>% 
    group_by(name) %>% 
    do(res = broom::tidy(cor.test(.$haqer_rare_variant_reversion_count, .$value))) %>% 
    unnest(res) %>% 
    arrange(p.value) %>% 
    filter(p.value < 0.05) %>% 
    mutate(lab = str_c('r = ', round(estimate, digits = 2), '\np-val = ', formatC(p.value, digits = 2)))

p_lang_milestone = ph %>% 
    select(spid, matches('words|phrase')) %>% 
    filter(used_words_age_mos < 888 & combined_words_age_mos < 888 & combined_phrases_age_mos < 888) %>%
    drop_na() %>%
    pivot_longer(cols = where(is.numeric)) %>% 
    filter(str_detect(name, 'words|phrase')) %>% 
    drop_na(value) %>%
    filter(value < 888) %>%
    drop_na(value) %>%
    inner_join(kp) %>% 
    inner_join(word_age_res) %>%
    mutate(name = case_when(name == 'used_words_age_mos' ~ 'First word',
                            name == 'combined_words_age_mos' ~ 'First time combining words',
                            name == 'combined_phrases_age_mos' ~ 'First time combining phrases'),
           name = factor(name, levels = c('First word', 'First time combining words', 'First time combining phrases'))) %>%
    ggplot(aes(y = haqer_rare_variant_reversion_count, x = value)) +
    geom_jitter(alpha = 0.7) +
    geom_smooth(method = 'lm') +
    facet_wrap(~name) +
    geom_text(aes(y = 23, x = 70, label = lab), size = 4, check_overlap = TRUE) +
    ylab('HAQER rare reversion count') +
    xlab('Age of language milestone (months)') +
    theme_classic() +
    theme(axis.text = element_text(size = 14),
          axis.title = element_text(size = 16),
          strip.text = element_text(size = 16))

quantile(tmp$combined_words_age_mos)
p_dat <- ph %>% 
    select(spid, matches('words|phrase')) %>% 
    filter(used_words_age_mos < 888 & combined_words_age_mos < 888 & combined_phrases_age_mos < 888) %>%
    drop_na() %>%
    pivot_longer(cols = where(is.numeric)) %>% 
    filter(str_detect(name, 'words|phrase')) %>% 
    drop_na(value) %>%
    filter(value < 888) %>%
    drop_na(value) %>%
    inner_join(kp) %>% 
    inner_join(word_age_res) %>%
    group_by(name) %>% 
    mutate(age_tile = ntile(value, n = 6),
           age_tile = case_when(value < 12 ~ '< 1 year old',
                                value >= 12 & value <= 36 ~ '1-2 years old',
                                value > 36 ~ '3 years and older'
                                ),
            age_tile = factor(age_tile, levels = c('< 1 year old', '1-2 years old', '3 years and older'))
            ) %>% 
    ungroup() %>%
    mutate(name = case_when(name == 'used_words_age_mos' ~ 'First word',
                            name == 'combined_words_age_mos' ~ 'First time combining words',
                            name == 'combined_phrases_age_mos' ~ 'First time combining phrases'),
           name = factor(name, levels = c('First word', 'First time combining words', 'First time combining phrases')),
            xlb = case_when(age_tile == 1 ~ 'Youngest 25%', 
                            age_tile == 2 ~ '25-50%', 
                            age_tile == 3 ~ '50-75%', 
                            age_tile == 4 ~ 'Oldest 25%'
                            ),
            xlb = factor(xlb, levels = c('Youngest 25%', '25-50%', '50-75%', 'Oldest 25%'))) %>% 
    group_by(name, ) %>% 
    mutate(mean_rev = mean(haqer_rare_variant_reversion_count)) %>% 
    ungroup()
p_wa <- p_dat %>%
    filter(name == 'First word') %>%
    ggplot(aes(y = haqer_rare_variant_reversion_count, x = as.factor(age_tile))) +
    geom_violin() +
    geom_boxplot(aes(fill = age_tile), width = 0.3, alpha = .9) +
    geom_smooth(method = 'lm') +
    # facet_wrap(~name) +
    # geom_text(aes(y = 23, x = 70, label = lab), size = 4, check_overlap = TRUE) +
    ylab('HAQER rare reversion count') +
    xlab('Age of first word') +
    theme_classic() +
    theme(axis.text = element_text(size = 14),
          axis.title = element_text(size = 16),
          strip.text = element_text(size = 16),
          legend.position = 'none') +
    geom_hline(yintercept = me, color = 'red', linetype = 'dashed', size = 1.075) +
    scale_fill_manual(values = c('#ffe6ad', '#e04468', '#762776'))
names(tmp)

quantile(tmp$combined_words_age_mos[tmp$combined_words_age_mos < 888], na.rm = TRUE)

p_wc <- p_dat %>%
    filter(name == 'First time combining words') %>%
    mutate(age_tile = ntile(value, n = 6),
           age_tile = case_when(value < 18 ~ '< 1.5 years old',
                                value >= 18 & value <= 48 ~ '1.5-3 years old',
                                value > 48 ~ '4 years and older'
                                ),
            age_tile = factor(age_tile, levels = c('< 1.5 years old', '1.5-3 years old', '4 years and older'))
            ) %>% 
    ggplot(aes(y = haqer_rare_variant_reversion_count, x = as.factor(age_tile))) +
    geom_violin() +
    geom_boxplot(aes(fill = age_tile), width = 0.3, alpha = .9) +
    geom_smooth(method = 'lm') +
    # facet_wrap(~name) +
    # geom_text(aes(y = 23, x = 70, label = lab), size = 4, check_overlap = TRUE) +
    ylab(NULL) +
    xlab('Age when first combined words') +
    theme_classic() +
    theme(axis.text = element_text(size = 14),
          axis.title = element_text(size = 16),
          strip.text = element_text(size = 16),
          legend.position = 'none') +
    geom_hline(yintercept = me, color = 'red', linetype = 'dashed', size = 1.075) +
    scale_fill_manual(values = c('#ffe6ad', '#e04468', '#762776'))

quantile(tmp$combined_phrases_age_mos[tmp$combined_phrases_age_mos < 888], na.rm = TRUE)
p_cp <- p_dat %>%
    filter(name == 'First time combining phrases') %>%
    mutate(age_tile = ntile(value, n = 6),
           age_tile = case_when(value < 24 ~ '< 2 years old',
                                value >= 24 & value <= 60 ~ '2-4 years old',
                                value > 60 ~ '5 years and older'
                                ),
            age_tile = factor(age_tile, levels = c('< 2 years old', '2-4 years old', '5 years and older'))
            ) %>%
    ggplot(aes(y = haqer_rare_variant_reversion_count, x = as.factor(age_tile))) +
    geom_violin() +
    geom_boxplot(aes(fill = age_tile), width = 0.3, alpha = .9) +
    geom_smooth(method = 'lm') +
    # facet_wrap(~name) +
    # geom_text(aes(y = 23, x = 70, label = lab), size = 4, check_overlap = TRUE) +
    ylab(NULL) +
    xlab('Age when first combined phrases') +
    theme_classic() +
    theme(axis.text = element_text(size = 14),
          axis.title = element_text(size = 16),
          strip.text = element_text(size = 16),
          legend.position = 'none') +
    geom_hline(yintercept = me, color = 'red', linetype = 'dashed', size = 1.075) +
    scale_fill_manual(values = c('#ffe6ad', '#e04468', '#762776'))
library(patchwork)
p <- p_wa + p_wc + p_cp + plot_annotation(tag_levels = 'A') & theme(plot.tag = element_text(size = 20, face = 'bold'), axis.title = element_text(face = 'bold'))
ggsave(p, filename = '/wdata/lcasten/sli_wgs/HCA_reversion/figures/SPARK_cog_lang_milestones_binned_rare_HAQER_reversions.png', dpi = 300, units = 'in', device = 'png', width = 20, height = 5.5)

#####################################################3
## diagnosis analysis
#####################################################
ph %>% 
    select(spid, sex, asd, age_at_registration_years, matches('dev_|adhd|mood'))  %>% 
    drop_na() %>%
    mutate(asd = as.numeric(asd)) %>%
    pivot_longer(cols = matches('asd|dev_|adhd|mood')) %>% 
    inner_join(brd) %>% 
    filter(haqer_rare_variant_reversion_count <= me + (3 * md) & haqer_rare_variant_reversion_count >= me - (3 * md)) %>%
    group_by(name) %>% 
    do(res = broom::tidy(glm(haqer_rare_variant_reversion_count ~ as.factor(value), data = .))) %>% 
    unnest(res) %>% 
    filter(str_detect(term, 'haqer_rare_variant_reversion_count|value')) %>%
    arrange(p.value) %>% 
    mutate(lab = str_c('Beta = ', round(estimate, digits = 2), '\np-val = ', formatC(p.value, digits = 2)))


res <- ph %>% 
    select(spid, sex, asd, age_at_registration_years, matches('dev_|adhd|mood'))  %>% 
    drop_na() %>%
    mutate(asd = as.numeric(asd)) %>%
    pivot_longer(cols = matches('asd|dev_|adhd|mood')) %>% 
    inner_join(brd) %>% 
    filter(haqer_rare_variant_reversion_count <= me + (3 * md) & haqer_rare_variant_reversion_count >= me - (3 * md)) %>%
    group_by(name) %>% 
    do(res = broom::tidy(glm(value ~ scale(haqer_rare_variant_reversion_count)[,1], family = 'binomial', data = .))) %>% 
    unnest(res) %>% 
    filter(str_detect(term, 'haqer_rare_variant_reversion_count')) %>%
    arrange(p.value) %>% 
    mutate(lab = str_c('Beta = ', round(estimate, digits = 2), '\np-val = ', formatC(p.value, digits = 2)))

p_res <- res %>% 
    mutate(xmin = estimate - 1.96 * std.error, xmax = estimate + 1.96 * std.error) %>%
    filter(p.value < 0.05) %>% 
    arrange(estimate) %>% 
    mutate(name = factor(name, levels = unique(name))) %>%
    ggplot(aes(x = estimate, y = name)) +
    geom_point(size = 4) +
    geom_linerange(aes(xmin = xmin, xmax = xmax), size = 1.1) +
    geom_vline(xintercept = 0, color = 'red', linetype = 'dashed', size = 1.075) +
    xlab('Regression beta for rare HAQER reversion count (95% CI)') +
    ylab('Diagnosis') +
    theme_classic() +
    theme(axis.text = element_text(size = 14),
          axis.title = element_text(size = 16))


## cognitive impariment analysis
wilcox.test(tmp$haqer_rare_variant_reversion_count ~ tmp$asd)
wilcox.test(tmp$haqer_rare_variant_reversion_count ~ tmp$dev_lang_dis)
wilcox.test(tmp$haqer_rare_variant_reversion_count ~ tmp$asd)

res_ci <- tmp %>% 
    mutate(cog_imp = case_when(dev_id == 1 | cognitive_impairment_at_enrollment == TRUE ~ 1,
                               TRUE ~ 0)) %>% 
    filter(asd == TRUE) %>% 
    inner_join(kp) %>% 
    pivot_longer(cols = c(dev_id, cog_imp, dev_lang, regress_lang_y_n))%>% 
    drop_na(value) %>%
    group_by(name) %>% 
    do(res = broom::tidy(t.test(.$haqer_rare_variant_reversion_count ~ .$value))) %>% 
    unnest(res) %>% 
    mutate(lab = str_c('t = ', round(-1 * statistic, digits =2), '\np-val = ', formatC(p.value, digits = 2))) %>% 
    filter(name != 'dev_id')

lower_ci <- function(mean, se, n, conf_level = 0.95){
  lower_ci <- mean - qt(1 - ((1 - conf_level) / 2), n - 1) * se
}
upper_ci <- function(mean, se, n, conf_level = 0.95){
  upper_ci <- mean + qt(1 - ((1 - conf_level) / 2), n - 1) * se
}
p_imp <- tmp %>% 
    mutate(cog_imp = case_when(dev_id == 1 | cognitive_impairment_at_enrollment == TRUE ~ 1,
                               TRUE ~ 0)) %>% 
    filter(asd == TRUE) %>% 
    inner_join(kp) %>% 
    pivot_longer(cols = c(dev_id, cog_imp, dev_lang))%>%  ## regress_lang_y_n
    drop_na(value) %>% 
    group_by(name, value) %>% 
    mutate(n = n()) %>% 
    mutate(value = str_c(as.logical(value), '\nN = ', n)) %>% 
    group_by(name, value) %>%
    summarise(smean = mean(haqer_rare_variant_reversion_count, na.rm = TRUE),
              ssd = sd(haqer_rare_variant_reversion_count, na.rm = TRUE),
              count = n()) %>%
    mutate(se = ssd / sqrt(count),
           lower_ci = lower_ci(smean, se, count),
           upper_ci = upper_ci(smean, se, count)) %>% 
    inner_join(res_ci) %>% 
    mutate(name = case_when(name == 'cog_imp' ~ 'Intellectual disability/cog. impairment',
                            name == 'dev_lang' ~ 'Language disorder',
                            name == 'regress_lang_y_n' ~ 'Language regression')) %>%
    ggplot(aes(x = value, y = smean)) +
    geom_point(size = 4) +
    geom_linerange(aes(ymin = lower_ci, ymax = upper_ci), size = 1.2) +
    facet_wrap(~ name, scales = 'free_x') +
    geom_text(aes(x = 1.5, y = 10.8, label = lab), check_overlap = TRUE, size = 5) +
    xlab('Diagnosis') +
    ylab('HAQER rare reversion count (95% CI)') +
    theme_classic() +
    theme(axis.text = element_text(size = 14),
          axis.title = element_text(size = 16),
          strip.text = element_text(size = 16))

cog_imp = ifelse(tmp$dev_id[tmp$asd == TRUE] == 1 | tmp$cognitive_impairment_at_enrollment[tmp$asd == TRUE] == TRUE, 1, 0)
t.test(tmp$haqer_rare_variant_reversion_count[tmp$asd == TRUE] ~ tmp$cognitive_impairment_at_enrollment[tmp$asd == TRUE])
t.test(tmp$haqer_rare_variant_reversion_count[tmp$asd == TRUE] ~ tmp$dev_lang[tmp$asd == TRUE])
t.test(tmp$haqer_rare_variant_reversion_count[tmp$asd == TRUE] ~ cog_imp)

wilcox.test(tmp2$resid_haqer_rare_variant_reversion_count ~ tmp2$dev_lang_dis)
wilcox.test(tmp2$resid_haqer_rare_variant_reversion_count ~ tmp2$dev_lang)


############################################
## relative cog / language level analysis
############################################
tdat <- ph %>% 
    select(spid, sex, asd, age_at_registration_years, matches('level'))  %>% 
    drop_na() %>%
    inner_join(kp) %>% 
    filter(age_at_registration_years >= 2)
lang_level <- aov(tdat$haqer_rare_variant_reversion_count ~ tdat$language_level_at_enrollment)
summary(lang_level)
tukey.test <- TukeyHSD(lang_level)
tukey.test

unique(tdat$language_level_at_enrollment)
tdat %>% 
    mutate(language_level_at_enrollment = factor(language_level_at_enrollment, levels = c("No words/does not speak", "Uses single words meaningfully (for example, to request)", "Combines 3 words together into short sentences", "Uses longer sentences of his/her own and is able to tell you something that happened"))) %>% 
    group_by(language_level_at_enrollment) %>%
    summarise(smean = mean(haqer_rare_variant_reversion_count, na.rm = TRUE),
              ssd = sd(haqer_rare_variant_reversion_count, na.rm = TRUE),
              count = n()) %>%
    mutate(se = ssd / sqrt(count),
           lower_ci = lower_ci(smean, se, count),
           upper_ci = upper_ci(smean, se, count)) %>% 
    ggplot(aes(x = language_level_at_enrollment, y = smean)) +
    geom_point(size = 4) +
    geom_linerange(aes(ymin = lower_ci, ymax = upper_ci), size = 1.2) +
    # geom_text(aes(x = 1.5, y = 10.8, label = lab), check_overlap = TRUE, size = 5) +
    xlab('Diagnosis') +
    ylab('HAQER rare reversion count (95% CI)') +
    theme_classic() +
    theme(axis.text = element_text(size = 14),
          axis.title = element_text(size = 16),
          strip.text = element_text(size = 16))

summary(aov(tdat$haqer_rare_variant_reversion_count ~ tdat$cog_age_level))
summary(aov(tdat$haqer_rare_variant_reversion_count ~ tdat$language_age_level))


##########################
## IQ analysis
##########################
broom::tidy(cor.test(tmp$haqer_rare_variant_reversion_count, tmp$fsiq_score))
broom::tidy(cor.test(tmp$haqer_rare_variant_reversion_count, tmp$viq_score))
broom::tidy(cor.test(tmp$haqer_rare_variant_reversion_count, tmp$nviq_score))

iq_res <- tmp %>% 
    select(spid, haqer_rare_variant_reversion_count, matches('iq_score')) %>% 
    filter(haqer_rare_variant_reversion_count <= me + (3 * md) & haqer_rare_variant_reversion_count >= me - (3 * md)) %>%
    drop_na() %>% 
    pivot_longer(cols = matches('iq_score')) %>% 
    group_by(name) %>% 
    do(res = broom::tidy(cor.test(.$value, .$haqer_rare_variant_reversion_count))) %>% 
    unnest(res) %>% 
    mutate(lab = str_c('r = ', round(estimate, digits = 2), '\np-val = ', formatC(p.value, digits = 2)))

p_iq <- tmp %>% 
    select(spid, haqer_rare_variant_reversion_count, matches('iq_score')) %>% 
    filter(haqer_rare_variant_reversion_count <= me + (3 * md) & haqer_rare_variant_reversion_count >= me - (3 * md)) %>%
    drop_na() %>% 
    pivot_longer(cols = matches('iq_score')) %>% 
    inner_join(iq_res) %>% 
    mutate(name = case_when(name == 'fsiq_score' ~ 'Full scale IQ',
                            name == 'nviq_score' ~ 'Nonverbal IQ',
                            name == 'viq_score' ~ 'Verbal IQ'),
           name = factor(name, levels = c('Full scale IQ', 'Nonverbal IQ', 'Verbal IQ'))) %>%
    ggplot(aes(x = haqer_rare_variant_reversion_count, y = value)) +
    geom_jitter(size = 2) +
    geom_smooth(method = 'lm') +
    facet_wrap(~name) +
    geom_text(aes(x = 12.5, y = 150, label = lab), size = 4, check_overlap = TRUE) +
    xlab('HAQER rare reversion count') +
    ylab('IQ score') +
    theme_classic() +
    theme(axis.text = element_text(size = 14),
          axis.title = element_text(size = 16),
          strip.text = element_text(size = 16))

##############################3
## make figures
###############################
p_iq %>% 
    ggsave(filename = '/wdata/lcasten/sli_wgs/HCA_reversion/figures/SPARK_IQ_rare_HAQER_reversions.png', dpi = 300, units = 'in', device = 'png', width = 15, height = 6)
p_lang_milestone %>% 
    ggsave(filename = '/wdata/lcasten/sli_wgs/HCA_reversion/figures/SPARK_language_milestones_rare_HAQER_reversions.png', dpi = 300, units = 'in', device = 'png', width = 15, height = 6)
p_imp%>% 
    ggsave(filename = '/wdata/lcasten/sli_wgs/HCA_reversion/figures/SPARK_cog_lang_impairments_rare_HAQER_reversions.png', dpi = 300, units = 'in', device = 'png', width = 10, height = 7)

brd %>% 
    filter(haqer_rare_variant_reversion_count <= me + (3 * md) & haqer_rare_variant_reversion_count >= me - (3 * md)) %>%
    ggplot(aes(x = haqer_rare_variant_reversion_count)) +
    geom_histogram()







###########################################################
##############################################################
###########################################################

ph %>% 
    select(spid, sex, asd, age_at_registration_years, matches('dev_|adhd|mood'))  %>% 
    drop_na() %>%
    mutate(asd = as.numeric(asd)) %>%
    pivot_longer(cols = matches('asd|dev_|adhd|mood')) %>% 
    inner_join(brd) %>% 
    inner_join(res) %>% 
    filter(str_detect(name, 'asd|lang|soc_prag|soc_anx|mood_or_anx|mood_anx|mood_dep|motor|ocd')) %>% 
    arrange(desc(estimate)) %>% 
    group_by(name, value) %>% 
    mutate(n = n()) %>% 
    mutate(value = str_c(as.logical(value), '\nN = ', n)) %>%
    ungroup() %>%
    filter(haqer_rare_variant_reversion_count <= me + (3 * md) & haqer_rare_variant_reversion_count >= me - (3 * md)) %>%
    mutate(name = factor(name, levels = unique(name))) %>% 
    ggplot(aes(x = value, y = haqer_rare_variant_reversion_count)) +
    geom_violin() +
    geom_boxplot(width = 0.3, aes(fill = value)) +
    geom_hline(yintercept = median(brd$haqer_rare_variant_reversion_count), color = 'red', linetype = 'dashed', size = 1.075) +
    facet_wrap(~ name, scales = 'free') +
    geom_text(aes(x = 1.5, y = 20, label = lab), check_overlap = TRUE, size = 5) +
    theme_classic() +
    theme(legend.position = 'none') +
    coord_cartesian(ylim = c(5,15))


ph %>% 
    select(spid, sex, asd, age_at_registration_years, matches('dev_|adhd|mood'))  %>% 
    drop_na() %>%
    mutate(asd = as.numeric(asd)) %>%
    pivot_longer(cols = matches('asd|dev_|adhd|mood')) %>% 
    inner_join(brd) %>% 
    inner_join(res) %>% 
    filter(str_detect(name, 'asd|lang|soc_prag|soc_anx|mood_or_anx|mood_anx|mood_dep|motor|ocd')) %>% 
    arrange(desc(estimate)) %>% 
    group_by(name, value) %>% 
    mutate(n = n()) %>% 
    mutate(value = str_c(as.logical(value), '\nN = ', n)) %>%
    ungroup() %>%
    mutate(name = factor(name, levels = unique(name))) %>% 
    filter(haqer_rare_variant_reversion_count <= me + 3 * md & haqer_rare_variant_reversion_count >= me - 3 * md) %>%
    ggplot(aes(x = haqer_rare_variant_reversion_count, fill = value)) +
    geom_density() +
    geom_vline(xintercept = median(brd$haqer_rare_variant_reversion_count), color = 'red', linetype = 'dashed', size = 1.075) +
    facet_wrap(~ name, scales = 'free') +
    # geom_text(aes(x = 1.5, y = 20, label = lab), check_overlap = TRUE, size = 5) +
    theme_classic() +
    theme(legend.position = 'none')



#######################
ph %>% 
    select(spid, sex, asd, age_at_registration_years, matches('dev_|adhd|mood'))  %>% 
    drop_na() %>%
    mutate(asd = as.numeric(asd)) %>%
    pivot_longer(cols = matches('asd|dev_|adhd|mood')) %>% 
    inner_join(brd2) %>% 
    group_by(name) %>% 
    do(res = broom::tidy(glm(value ~ haqer_rare_variant_reversion_count + pc1 + pc2 + pc3 + pc4 + pc5, family = 'binomial', data = .))) %>% 
    unnest(res) %>% 
    filter(str_detect(term, 'haqer_rare_variant_reversion_count')) %>%
    arrange(p.value)