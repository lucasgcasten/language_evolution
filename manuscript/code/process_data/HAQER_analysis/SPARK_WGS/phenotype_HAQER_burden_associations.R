library(tidyverse)

#####################################
## read in data
#####################################
## phenos
ph <- read_csv('/wdata/lcasten/spark/WGS/data/WGS_gathered_phenotypes.csv')

brd <- read_csv('/wdata/lcasten/sli_wgs/HAQER_variant_associations/data/SPARK_WGS/HAQER_HAR_RAND_10Kb_flank_rare_variant_burden_counts.csv')

#################
pc <- read_csv('/wdata/lcasten/spark/prs/HapMap3_plus/PCA/raw_KING_pca_results.csv') %>% 
    rename(spid = sample.ID) %>% 
    select(1:6)

##
tmp <- ph %>% 
    # inner_join(pc) %>%
    inner_join(brd) %>% 
    mutate(asd = as.numeric(asd))

############################################################################################
## remove outliers from analysis (just to be safe, saw the same associations even when they're included)
############################################################################################
md <- mad(brd$haqer_rare_variant_reversion_count)
me <- median(brd$haqer_rare_variant_reversion_count)

kp <- brd %>% 
    pivot_longer(cols = matches('reversion')) %>% 
    group_by(name) %>% 
    mutate(medn = median(value),
           madv = mad(value)) %>%
    ungroup() %>%
    filter(value <= medn + (2.5 * madv) & value >= medn - (2.5 * madv)) %>% 
    pivot_wider(id_cols = spid) %>% 
    drop_na()

# kp <- brd %>% 
#     filter(haqer_rare_variant_reversion_count <= me + (3 * md) & haqer_rare_variant_reversion_count >= me - (3 * md))
names(ph)


#############################################################
## language development milestones
#############################################################
hist(ph$used_words_age_mos)
hist(ph$combined_words_age_mos)
hist(ph$combined_phrases_age_mos)

p_forest_age_dat <- ph %>% 
    select(spid, matches('words|phrase')) %>% 
    filter(used_words_age_mos < 888 & combined_words_age_mos < 888 & combined_phrases_age_mos < 888) %>%
    drop_na() %>%
    distinct(spid, .keep_all = TRUE) %>%
    pivot_longer(cols = where(is.numeric)) %>% 
    filter(str_detect(name, 'words|phrase')) %>% 
    drop_na(value) %>%
    filter(value < 888) %>%
    drop_na(value) %>%
    inner_join(kp) %>% 
    group_by(name) %>% 
    mutate(n = n()) %>% 
    group_by(name, n) %>% 
    # do(res = broom::tidy(cor.test(.$haqer_rare_variant_reversion_count, .$value))) %>% 
    do(res = broom::tidy(lm(scale(value)[,1] ~ scale(haqer_rare_variant_reversion_count)[,1] + scale(har_rare_variant_reversion_count)[,1] + scale(rand_rare_variant_reversion_count)[,1], data = .))) %>% 
    unnest(res) 

word_age_res <- ph %>% 
    select(spid, matches('words|phrase')) %>% 
    filter(used_words_age_mos < 888 & combined_words_age_mos < 888 & combined_phrases_age_mos < 888) %>%
    drop_na() %>%
    distinct(spid, .keep_all = TRUE) %>%
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
    distinct(spid, .keep_all = TRUE) %>%
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
    distinct(spid, .keep_all = TRUE) %>%
    pivot_longer(cols = where(is.numeric)) %>% 
    filter(str_detect(name, 'words|phrase')) %>% 
    drop_na(value) %>%
    filter(value < 888) %>%
    drop_na(value) %>%
    inner_join(kp) %>% 
    inner_join(word_age_res) %>%
    group_by(name) %>% 
    mutate(# age_tile = ntile(value, n = 6),
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
av <- aov(haqer_rare_variant_reversion_count ~ age_tile, data = p_dat[p_dat$name == 'First word',])
summary(av)
TukeyHSD(av)

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

td <- p_dat %>%
    filter(name == 'First time combining words') %>%
    mutate(age_tile = ntile(value, n = 6),
           age_tile = case_when(value < 18 ~ '< 1.5 years old',
                                value >= 18 & value <= 48 ~ '1.5-3 years old',
                                value > 48 ~ '4 years and older'
                                ),
            age_tile = factor(age_tile, levels = c('< 1.5 years old', '1.5-3 years old', '4 years and older'))
            )
av <- aov(haqer_rare_variant_reversion_count ~ age_tile, data = td)
summary(av)
TukeyHSD(av)

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

td <- p_dat %>%
    filter(name == 'First time combining phrases') %>%
    mutate(age_tile = ntile(value, n = 6),
           age_tile = case_when(value < 24 ~ '< 2 years old',
                                value >= 24 & value <= 60 ~ '2-4 years old',
                                value > 60 ~ '5 years and older'
                                ),
            age_tile = factor(age_tile, levels = c('< 2 years old', '2-4 years old', '5 years and older'))
            )
av <- aov(haqer_rare_variant_reversion_count ~ age_tile, data = td)
summary(av)
TukeyHSD(av)

library(patchwork)
p <- p_wa + p_wc + p_cp + plot_annotation(tag_levels = 'A') & theme(plot.tag = element_text(size = 20, face = 'bold'), axis.title = element_text(face = 'bold'))
ggsave(p, filename = '/wdata/lcasten/sli_wgs/HCA_reversion/figures/SPARK_cog_lang_milestones_binned_rare_HAQER_reversions.png', dpi = 300, units = 'in', device = 'png', width = 20, height = 5.5)

##############################################
## handedness analysis
##############################################
td <- ph %>% 
    select(spid, sex, asd, age_at_registration_years, hand)  %>% 
    drop_na() %>%
    distinct(spid, .keep_all = TRUE) %>%
    mutate(asd = as.numeric(asd)) %>%
    inner_join(kp) %>% 
    mutate(rhand = ifelse(hand == 'right', 'right handed', 'not right handed'))
t.test(td$haqer_rare_variant_reversion_count ~ td$rhand) ##
table(td$rhand)
table(td$hand)
av <- aov(haqer_rare_variant_reversion_count ~ hand, data = td)
summary(av)
TukeyHSD(av)

##############################################
## asd diagnosis analysis
##############################################
td <- ph %>% 
    select(spid, sex, asd, age_at_registration_years, diagnosis, diagnosis_age)  %>% 
    drop_na() %>%
    distinct(spid, .keep_all = TRUE) %>%
    filter(asd == TRUE) %>%
    inner_join(kp) %>% 
    filter(diagnosis != 'Autism Spectrum Disorder')
table(td$diagnosis)
range(td$diagnosis_age)
sum(td$diagnosis_age > 18 * 12)
av <- aov(haqer_rare_variant_reversion_count ~ diagnosis, data = td)
summary(av)
TukeyHSD(av)


##############################################
## histogram of rare reversion dist
##############################################
hist_dat <- kp %>% 
    distinct(spid, .keep_all = TRUE) %>%
    pivot_longer(cols = matches('reversion')) %>%
    group_by(name, value) %>% 
    count() 
unique(hist_dat$name)
range(hist_dat$value)

p_dist_haq <- hist_dat[hist_dat$name == 'haqer_rare_variant_reversion_count',] %>%
    ggplot(aes(x =value, y = n, fill = name)) +
    geom_bar(stat = 'identity', position = 'dodge') +
    coord_cartesian(xlim = c(0, max(hist_dat$value)), ylim = c(0, max(hist_dat$n))) +
    xlab('HAQER rare reversion count') +
    ylab(NULL) +
    scale_fill_manual(values = c('#762776')) +
    theme_classic() +
    theme(axis.text = element_text(size = 14),
          axis.title = element_text(size = 16),
          legend.position = 'none')

p_dist_har <- hist_dat[hist_dat$name == 'har_rare_variant_reversion_count',] %>%
    ggplot(aes(x =value, y = n, fill = name)) +
    geom_bar(stat = 'identity', position = 'dodge') +
    coord_cartesian(xlim = c(0, max(hist_dat$value)), ylim = c(0, max(hist_dat$n))) +
    xlab('HAR rare reversion count') +
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
    xlab('RAND rare reversion count') +
    ylab(NULL) +
    scale_fill_manual(values = c('#ffe6ad')) +
    theme_classic() +
    theme(axis.text = element_text(size = 14),
          axis.title = element_text(size = 16),
          legend.position = 'none')
library(patchwork)
p <- p_dist_haq / p_dist_har / p_dist_rand + plot_annotation(tag_levels = 'A') & theme(plot.tag = element_text(face = 'bold', size = 18))
ggsave(p, filename = '/wdata/lcasten/sli_wgs/HCA_reversion/figures/SPARK_rare_reversions_distributions.png', dpi = 300, units = 'in', device = 'png', width = 4.5, height = 8)

td <- ph %>% 
    select(spid, sex, asd, age_at_registration_years, diagnosis_age)  %>% 
    drop_na() %>%
    filter(asd == TRUE) %>%
    mutate(dx_before = ifelse(diagnosis_age < 12 * 5, 1, 0)) %>%
    inner_join(kp)
sum(td$diagnosis_age < 12 * 3)
mod <- glm(dx_before ~ scale(haqer_rare_variant_reversion_count)[,1] + scale(har_rare_variant_reversion_count)[,1] + scale(rand_rare_variant_reversion_count)[,1], data = td, family = 'binomial')
summary(mod)
t.test(td$haqer_rare_variant_reversion_count ~ td$dx_before)

boxplot(td$haqer_rare_variant_reversion_count ~ td$dx_before)
t.test(td$har_rare_variant_reversion_count ~ td$dx_before)
t.test(td$rand_rare_variant_reversion_count ~ td$dx_before)

mod <- lm(scale(diagnosis_age)[,1] ~ scale(haqer_rare_variant_reversion_count)[,1] + scale(har_rare_variant_reversion_count)[,1] + scale(rand_rare_variant_reversion_count)[,1], data = td)
summary(mod)

#####################################################3
## diagnosis analysis
#####################################################
p_forest_dat_dx <- ph %>% 
    select(spid, sex, asd, age_at_registration_years, matches('dev_|adhd|mood'))  %>% 
    drop_na() %>%
    mutate(asd = as.numeric(asd)) %>%
    pivot_longer(cols = matches('asd|dev_|adhd|mood')) %>% 
    inner_join(kp) %>% 
    distinct(spid, name, .keep_all = TRUE) %>%
    # filter(haqer_rare_variant_reversion_count <= me + (3 * md) & haqer_rare_variant_reversion_count >= me - (3 * md)) %>%
    group_by(name) %>% 
    mutate(n = n()) %>% 
    group_by(name, n) %>% 
    do(res = broom::tidy(glm(value ~ scale(haqer_rare_variant_reversion_count)[,1] + scale(har_rare_variant_reversion_count)[,1] + scale(rand_rare_variant_reversion_count)[,1], family = 'binomial', data = .))) %>% 
    unnest(res) %>% 
    # filter(str_detect(term, 'rand_rare_variant_reversion_count|value')) %>%
    arrange(p.value) %>% 
    mutate(lab = str_c('Beta = ', round(estimate, digits = 2), '\np-val = ', formatC(p.value, digits = 2)))
p_forest_dat_dx %>% filter(term != '(Intercept)') %>% filter(str_detect(term, 'haq'))
p_forest_dat_dx %>% filter(term != '(Intercept)') %>% filter(str_detect(term, 'har'))
p_forest_dat_dx %>% filter(term != '(Intercept)') %>% filter(str_detect(term, 'rand'))

res <- ph %>% 
    select(spid, sex, asd, age_at_registration_years, matches('dev_|adhd|mood'))  %>% 
    drop_na() %>%
    distinct(spid, .keep_all = TRUE) %>%
    mutate(asd = as.numeric(asd)) %>%
    pivot_longer(cols = matches('asd|dev_|adhd|mood')) %>% 
    inner_join(kp) %>% 
    # filter(haqer_rare_variant_reversion_count <= me + (3 * md) & haqer_rare_variant_reversion_count >= me - (3 * md)) %>%
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
kp_samp <- ph %>% 
    select(spid, sex, asd, age_at_registration_years, matches('dev_|adhd|mood'))  %>% 
    drop_na() %>%
    mutate(asd = as.numeric(asd)) %>%
    pivot_longer(cols = matches('asd|dev_|adhd|mood')) %>% 
    inner_join(kp) %>% 
    distinct(spid)
res_ci <- ph %>% 
    mutate(cog_imp = case_when(dev_id == 1 | cognitive_impairment_at_enrollment == TRUE ~ 1,
                               TRUE ~ 0)) %>% 
    # filter(asd == TRUE) %>% 
    inner_join(kp) %>% 
    filter(spid %in% kp_samp$spid) %>%
    distinct(spid, .keep_all = TRUE) %>%
    pivot_longer(cols = c(dev_id, cog_imp, matches('dev_lang'), regress_lang_y_n))%>%  #
    drop_na(value) %>%
    group_by(name) %>%
    mutate(n = n()) %>% 
    group_by(name, n) %>% 
    # do(res = broom::tidy(t.test(.$haqer_rare_variant_reversion_count ~ .$value))) %>% 
    do(res = broom::tidy(glm(value ~ scale(haqer_rare_variant_reversion_count)[,1] + scale(har_rare_variant_reversion_count)[,1] + scale(rand_rare_variant_reversion_count)[,1], family = 'binomial', data = .))) %>% 
    unnest(res) %>% 
    mutate(lab = str_c('Beta = ', round(estimate, digits =2), '\np-val = ', formatC(p.value, digits = 2))) %>% 
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

#####################################################
## more diagnosis analysis
#####################################################
dx2 <- read_csv('/sdata/Simons/SPARK/DATA/phenotypes/SPARK_Collection_v13_DataRelease_2024-09-24/basic_medical_screening-2024-09-24.csv')
names(dx2)
addt_dx_res <- dx2 %>% 
    select(spid = subject_sp_id, sex, asd, age_at_eval_years, matches('croceph|eating_disorder|neuro_|pers_dis|schiz|sleep_dx|tics|visaud|behav_|birth_')) %>% 
    pivot_longer(cols = -c(spid, sex, asd, age_at_eval_years)) %>% 
    # filter(asd == TRUE) %>%
    mutate(value = ifelse(is.na(value), 0, value))  %>%
    inner_join(kp) %>% 
    group_by(name) %>%
    # do(res = broom::tidy(glm(value ~ as.factor(asd) + scale(haqer_rare_variant_reversion_count)[,1] + scale(har_rare_variant_reversion_count)[,1] + scale(rand_rare_variant_reversion_count)[,1] + as.factor(sex) + age_at_eval_years, data = ., family = 'binomial'))) %>% 
    do(res = broom::tidy(glm(value ~ as.factor(asd) + scale(haqer_rare_variant_reversion_count)[,1] + as.factor(sex) + age_at_eval_years, data = ., family = 'binomial'))) %>% 
    unnest(res) %>% 
    filter(term != '(Intercept)') %>% 
    arrange(p.value)
addt_dx_res %>% 
    filter(str_detect(term, 'haqer'))
addt_dx_res %>% 
    filter(str_detect(term, 'har'))
addt_dx_res %>% 
    filter(str_detect(term, 'rand'))

#####################################################3
## vineland analysis
#####################################################
p_forest_dat_vineland <- ph %>% 
    select(spid, sex, asd, age_at_registration_years, matches('vineland_'))  %>% 
    drop_na() %>%
    mutate(asd = as.numeric(asd)) %>%
    pivot_longer(cols = matches('vineland_')) %>% 
    inner_join(kp) %>% 
    distinct(spid, name, .keep_all = TRUE) %>%
    # filter(haqer_rare_variant_reversion_count <= me + (3 * md) & haqer_rare_variant_reversion_count >= me - (3 * md)) %>%
    group_by(name) %>% 
    mutate(n = n()) %>% 
    group_by(name, n) %>% 
    do(res = broom::tidy(glm(value ~ scale(haqer_rare_variant_reversion_count)[,1] + scale(har_rare_variant_reversion_count)[,1] + scale(rand_rare_variant_reversion_count)[,1], data = .))) %>% 
    unnest(res) %>% 
    # filter(str_detect(term, 'rand_rare_variant_reversion_count|value')) %>%
    arrange(p.value) %>% 
    mutate(lab = str_c('Beta = ', round(estimate, digits = 2), '\np-val = ', formatC(p.value, digits = 2)))


############################################
## relative cog / language level analysis
############################################
tdat <- ph %>% 
    select(spid, sex, asd, age_at_registration_years, matches('level'))  %>% 
    drop_na() %>%
    distinct(spid, .keep_all = TRUE) %>%
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
table(ph$iq_invalid)
iq13 <- read_csv('/sdata/Simons/SPARK/DATA/phenotypes/SPARK_Collection_v13_DataRelease_2024-09-24/iq-2024-09-24.csv')
iq_res <- ph %>%
    # iq13 %>% 
    # select(spid = subject_sp_id, iq_test_age_months = age_test_date_months, iq_invalid = invalid, matches('iq_score')) %>%
    inner_join(kp) %>%
    # filter(asd == TRUE) %>%
    # filter(iq_invalid == 0) %>% 
    select(spid, iq_test_age_months, matches('reversion_count'), matches('iq_score')) %>% #asd, sex
    drop_na(fsiq_score, viq_score, nviq_score) %>% 
    filter(iq_test_age_months >= 5 * 12) %>%
    # filter(haqer_rare_variant_reversion_count <= me + (3 * md) & haqer_rare_variant_reversion_count >= me - (3 * md)) %>%
    distinct(spid, .keep_all = TRUE) %>%
    pivot_longer(cols = matches('iq_score')) %>% 
    group_by(name) %>% 
    do(res = broom::tidy(cor.test(.$value, .$haqer_rare_variant_reversion_count, method = 'p'))) %>% 
    # do(res = broom::tidy(glm(value ~  scale(haqer_rare_variant_reversion_count)[,1] + scale(har_rare_variant_reversion_count)[,1] + scale(rand_rare_variant_reversion_count)[,1], data = .))) %>% 
    unnest(res) %>% 
    mutate(lab = str_c('r = ', round(estimate, digits = 2), ', p-val = ', formatC(p.value, digits = 2), ', N = ', parameter + 2))

p_iq <- ph %>%
    inner_join(kp) %>%
    # filter(asd == TRUE) %>%
    # filter(iq_invalid == 0) %>% 
    select(spid, asd, iq_test_age_months, sex, matches('reversion_count'), matches('iq_score')) %>%
    drop_na(fsiq_score, viq_score, nviq_score) %>% 
    filter(iq_test_age_months >= 5 * 12) %>%
    # filter(haqer_rare_variant_reversion_count <= me + (3 * md) & haqer_rare_variant_reversion_count >= me - (3 * md)) %>%
    distinct(spid, .keep_all = TRUE) %>%
    pivot_longer(cols = matches('iq_score')) %>% 
    inner_join(iq_res) %>% 
    mutate(name = case_when(name == 'fsiq_score' ~ 'Full scale IQ',
                            name == 'nviq_score' ~ 'Nonverbal IQ',
                            name == 'viq_score' ~ 'Verbal IQ'),
           name = factor(name, levels = c('Full scale IQ', 'Nonverbal IQ', 'Verbal IQ'))) %>%
    mutate(sig = ifelse(p.value < 0.05, TRUE, FALSE)) %>%
    ggplot(aes(x = haqer_rare_variant_reversion_count, y = value)) +
    geom_jitter(size = 2, aes(color = sig)) +
    geom_smooth(method = 'lm', size = 1.5, alpha = 0.25) +
    facet_wrap(~name) +
    geom_text(aes(x = 12.5, y = 150, label = lab), size = 6, check_overlap = TRUE) +
    xlab('HAQER rare reversion count') +
    ylab('Clinical IQ score') +
    theme_classic() +
    theme(axis.text = element_text(size = 16),
          axis.title = element_text(size = 18),
          strip.text = element_text(size = 18),
          legend.position = 'none') +
    scale_color_manual(values = c('grey85', 'black'))

###################################################
## SCQ, DCDQ, RSR-R scores
###################################################
bh <- read_csv('/sdata/Simons/SPARK/DATA/phenotypes/SPARK_Collection_v13_DataRelease_2024-09-24/core_descriptive_variables-2024-09-24.csv')
names(bh)
bh_res <- bh %>%
    # iq13 %>% 
    select(spid = subject_sp_id, age_at_registration_months, sex, matches('total_final_score|dcdq_|vineland_abc_ss_latest|diagnosis_age')) %>%
    inner_join(kp) %>%
    # filter(asd == TRUE) %>%
    # filter(iq_invalid == 0) %>% 
    # select(spid, iq_test_age_months, matches('reversion_count'), matches('iq_score')) %>% #asd, sex
    # drop_na(fsiq_score, viq_score, nviq_score) %>% 
    # filter(iq_test_age_months >= 5 * 12) %>%
    # filter(haqer_rare_variant_reversion_count <= me + (3 * md) & haqer_rare_variant_reversion_count >= me - (3 * md)) %>%
    distinct(spid, .keep_all = TRUE) %>%
    pivot_longer(cols = matches('total_final_score|dcdq_|vineland_abc_ss_latest|diagnosis_age')) %>% 
    drop_na() %>%
    group_by(name) %>% 
    # do(res = broom::tidy(cor.test(.$value, .$haqer_rare_variant_reversion_count, method = 'p'))) %>% 
    do(res = broom::tidy(glm(value ~  scale(haqer_rare_variant_reversion_count)[,1] + scale(har_rare_variant_reversion_count)[,1] + scale(rand_rare_variant_reversion_count)[,1] + as.factor(sex) + scale(age_at_registration_months)[,1], data = .))) %>% 
    unnest(res) %>% 
    filter(str_detect(term, 'rare_'))

p_iq <- ph %>%
    inner_join(kp) %>%
    # filter(asd == TRUE) %>%
    # filter(iq_invalid == 0) %>% 
    select(spid, asd, iq_test_age_months, sex, matches('reversion_count'), matches('iq_score')) %>%
    drop_na(fsiq_score, viq_score, nviq_score) %>% 
    filter(iq_test_age_months >= 5 * 12) %>%
    # filter(haqer_rare_variant_reversion_count <= me + (3 * md) & haqer_rare_variant_reversion_count >= me - (3 * md)) %>%
    distinct(spid, .keep_all = TRUE) %>%
    pivot_longer(cols = matches('iq_score')) %>% 
    inner_join(iq_res) %>% 
    mutate(name = case_when(name == 'fsiq_score' ~ 'Full scale IQ',
                            name == 'nviq_score' ~ 'Nonverbal IQ',
                            name == 'viq_score' ~ 'Verbal IQ'),
           name = factor(name, levels = c('Full scale IQ', 'Nonverbal IQ', 'Verbal IQ'))) %>%
    mutate(sig = ifelse(p.value < 0.05, TRUE, FALSE)) %>%
    ggplot(aes(x = haqer_rare_variant_reversion_count, y = value)) +
    geom_jitter(size = 2, aes(color = sig)) +
    geom_smooth(method = 'lm', size = 1.5, alpha = 0.25) +
    facet_wrap(~name) +
    geom_text(aes(x = 12.5, y = 150, label = lab), size = 6, check_overlap = TRUE) +
    xlab('HAQER rare reversion count') +
    ylab('Clinical IQ score') +
    theme_classic() +
    theme(axis.text = element_text(size = 16),
          axis.title = element_text(size = 18),
          strip.text = element_text(size = 18),
          legend.position = 'none') +
    scale_color_manual(values = c('grey85', 'black'))


###################################################
## ASR scores
###################################################
bh <- read_csv('/sdata/Simons/SPARK/DATA/phenotypes/SPARK_Collection_v13_DataRelease_2024-09-24/asr-2024-09-24.csv')
names(bh)
bh_res <- bh %>%
    # iq13 %>% 
    select(spid = subject_sp_id, asd, age_at_eval_months, sex, matches('_t_|q079')) %>%
    inner_join(kp) %>%
    # filter(asd == TRUE) %>%
    # filter(iq_invalid == 0) %>% 
    # select(spid, iq_test_age_months, matches('reversion_count'), matches('iq_score')) %>% #asd, sex
    # drop_na(fsiq_score, viq_score, nviq_score) %>% 
    # filter(iq_test_age_months >= 5 * 12) %>%
    # filter(haqer_rare_variant_reversion_count <= me + (3 * md) & haqer_rare_variant_reversion_count >= me - (3 * md)) %>%
    distinct(spid, .keep_all = TRUE) %>%
    pivot_longer(cols = matches('_t_|q079')) %>% 
    drop_na() %>%
    filter(str_detect(name, '_t_', negate = TRUE)) %>%
    group_by(name) %>% 
    # do(res = broom::tidy(cor.test(.$value, .$haqer_rare_variant_reversion_count, method = 'p'))) %>% 
    do(res = broom::tidy(glm(value ~  scale(haqer_rare_variant_reversion_count)[,1] + scale(har_rare_variant_reversion_count)[,1] + scale(rand_rare_variant_reversion_count)[,1] + as.factor(asd) + as.factor(sex) + scale(age_at_eval_months)[,1], 
                             data = ., family = 'quasipoisson'))) %>% 
    unnest(res) %>% 
    filter(str_detect(term, 'rare_')) %>% 
    arrange(p.value)

bh_res %>% 
    filter(str_detect(term, 'haqer')) %>% 
    head(n = 20)
bh_res %>% 
    filter(str_detect(name, 'speech'))

p_iq <- ph %>%
    inner_join(kp) %>%
    # filter(asd == TRUE) %>%
    # filter(iq_invalid == 0) %>% 
    select(spid, asd, iq_test_age_months, sex, matches('reversion_count'), matches('iq_score')) %>%
    drop_na(fsiq_score, viq_score, nviq_score) %>% 
    filter(iq_test_age_months >= 5 * 12) %>%
    # filter(haqer_rare_variant_reversion_count <= me + (3 * md) & haqer_rare_variant_reversion_count >= me - (3 * md)) %>%
    distinct(spid, .keep_all = TRUE) %>%
    pivot_longer(cols = matches('iq_score')) %>% 
    inner_join(iq_res) %>% 
    mutate(name = case_when(name == 'fsiq_score' ~ 'Full scale IQ',
                            name == 'nviq_score' ~ 'Nonverbal IQ',
                            name == 'viq_score' ~ 'Verbal IQ'),
           name = factor(name, levels = c('Full scale IQ', 'Nonverbal IQ', 'Verbal IQ'))) %>%
    mutate(sig = ifelse(p.value < 0.05, TRUE, FALSE)) %>%
    ggplot(aes(x = haqer_rare_variant_reversion_count, y = value)) +
    geom_jitter(size = 2, aes(color = sig)) +
    geom_smooth(method = 'lm', size = 1.5, alpha = 0.25) +
    facet_wrap(~name) +
    geom_text(aes(x = 12.5, y = 150, label = lab), size = 6, check_overlap = TRUE) +
    xlab('HAQER rare reversion count') +
    ylab('Clinical IQ score') +
    theme_classic() +
    theme(axis.text = element_text(size = 16),
          axis.title = element_text(size = 18),
          strip.text = element_text(size = 18),
          legend.position = 'none') +
    scale_color_manual(values = c('grey85', 'black'))

#########################################
## LSR items 
#########################################
dm = read_csv('/wdata/lcasten/spark/research_match/language_merged_results/data/ID_mapping_and_demographics.csv')
bh <- readxl::read_xlsx('/sdata/Simons/SPARK/ResearchMatch/language_2023/Speech_and_Communication_Questionnaire_Responses_v4.xlsx')
tmp <- bh %>%
    # iq13 %>% 
    select(spid = Subject_ParticipantId, matches('^q[0-9]')) %>%
    select(-matches('_text|_deceased_age|neurologic'))
msng <- rowSums(is.na(tmp[,-1]))

msng_cutoff <- round(.05 * (ncol(tmp) - 1), digits = 0)
msng_cutoff 
msng_cutoff = 1
to_drop <- msng > msng_cutoff
to_kp <- msng <= msng_cutoff
tmp <- tmp[to_kp, ]
sum(msng > msng_cutoff)

for(i in 2:ncol(tmp)){
  tmp[is.na(tmp[,i]), i] <- round(mean(unlist(tmp[,i]), na.rm = TRUE), digits = 0)
}

names(tmp)
hist(rowSums(tmp[,-1]))
tmp$total_score <- rowSums(tmp[,-1]) ## to filter out people who didn't put in effort (selected almost all 0's or 2's which should be impossible with all the reverse scorig)
median(tmp$total_score)
mad(tmp$total_score)

lsr_kp <- tmp %>% 
    select(1, total_score) %>% 
    arrange(total_score) %>%
    # filter(total_score > 0) 
    filter(total_score >= median(tmp$total_score, na.rm = TRUE) - 2.5 * mad(tmp$total_score, na.rm = TRUE) & total_score <= median(tmp$total_score, na.rm = TRUE) + 2.5 * mad(tmp$total_score, na.rm = TRUE))
names(bh)

bh_res <- tmp %>%
    # bh %>%
    # iq13 %>% 
    # select(spid = Subject_ParticipantId, matches('^q[0-9]')) %>%
    select(-matches('_text')) %>%
    # drop_na() %>%
    inner_join(dm) %>% 
    inner_join(kp) %>%
    distinct(spid, .keep_all = TRUE) %>%
    # filter(diagnosis == 'ASD') %>%
    filter(spid %in% lsr_kp$spid) %>%
    pivot_longer(cols = matches('^q[0-9]')) %>% 
    drop_na() %>%
    group_by(name) %>% 
    # do(res = broom::tidy(cor.test(.$value, .$haqer_rare_variant_reversion_count, method = 'p'))) %>% 
    do(res = broom::tidy(glm(value ~  scale(haqer_rare_variant_reversion_count)[,1] + scale(har_rare_variant_reversion_count)[,1] + scale(rand_rare_variant_reversion_count)[,1] + as.factor(sex) + scale(age)[,1] + as.factor(diagnosis), 
                             data = ., family = 'quasipoisson'))) %>% 
    # do(res = broom::tidy(glm(value ~  scale(haqer_rare_variant_reversion_count)[,1] + as.factor(sex) + scale(age)[,1] + as.factor(diagnosis), 
    #                          data = ., family = 'quasipoisson'))) %>% 
    unnest(res) %>% 
    filter(str_detect(term, 'rare_')) %>% 
    arrange(p.value)

bh_res %>% 
    filter(str_detect(term, 'haqer')) %>% 
    head(n = 20)

bh_res %>% 
    filter(str_detect(term, 'har')) %>% 
    head(n = 20)

p_iq <- ph %>%
    inner_join(kp) %>%
    # filter(asd == TRUE) %>%
    # filter(iq_invalid == 0) %>% 
    select(spid, asd, iq_test_age_months, sex, matches('reversion_count'), matches('iq_score')) %>%
    drop_na(fsiq_score, viq_score, nviq_score) %>% 
    filter(iq_test_age_months >= 5 * 12) %>%
    # filter(haqer_rare_variant_reversion_count <= me + (3 * md) & haqer_rare_variant_reversion_count >= me - (3 * md)) %>%
    distinct(spid, .keep_all = TRUE) %>%
    pivot_longer(cols = matches('iq_score')) %>% 
    inner_join(iq_res) %>% 
    mutate(name = case_when(name == 'fsiq_score' ~ 'Full scale IQ',
                            name == 'nviq_score' ~ 'Nonverbal IQ',
                            name == 'viq_score' ~ 'Verbal IQ'),
           name = factor(name, levels = c('Full scale IQ', 'Nonverbal IQ', 'Verbal IQ'))) %>%
    mutate(sig = ifelse(p.value < 0.05, TRUE, FALSE)) %>%
    ggplot(aes(x = haqer_rare_variant_reversion_count, y = value)) +
    geom_jitter(size = 2, aes(color = sig)) +
    geom_smooth(method = 'lm', size = 1.5, alpha = 0.25) +
    facet_wrap(~name) +
    geom_text(aes(x = 12.5, y = 150, label = lab), size = 6, check_overlap = TRUE) +
    xlab('HAQER rare reversion count') +
    ylab('Clinical IQ score') +
    theme_classic() +
    theme(axis.text = element_text(size = 16),
          axis.title = element_text(size = 18),
          strip.text = element_text(size = 18),
          legend.position = 'none') +
    scale_color_manual(values = c('grey85', 'black'))

#############################
## LSR factors
#############################
lsr_fac <- read_csv("/wdata/gsnyder2/language/code/adjusted_factor_scores")
lsr_res <- lsr_fac %>% 
    relocate(spid = Subject_ParticipantId) %>% 
    inner_join(dm) %>% 
    inner_join(kp) %>%
    distinct(spid, .keep_all = TRUE) %>%
    # filter(diagnosis == 'ASD') %>%
    filter(spid %in% lsr_kp$spid) %>%
    pivot_longer(cols = matches('^MR')) %>% 
    drop_na() %>%
    group_by(name) %>% 
    # do(res = broom::tidy(cor.test(.$value, .$haqer_rare_variant_reversion_count, method = 'p'))) %>% 
    do(res = broom::tidy(glm(value ~  scale(haqer_rare_variant_reversion_count)[,1] + scale(har_rare_variant_reversion_count)[,1] + scale(rand_rare_variant_reversion_count)[,1] + as.factor(sex) + scale(age)[,1] + as.factor(diagnosis), 
                             data = .))) %>% 
    # do(res = broom::tidy(glm(value ~  scale(haqer_rare_variant_reversion_count)[,1] + as.factor(sex) + scale(age)[,1] + as.factor(diagnosis), 
    #                          data = ., family = 'quasipoisson'))) %>% 
    unnest(res) %>% 
    filter(str_detect(term, 'rare_')) %>% 
    arrange(p.value)

bh_res %>% 
    filter(str_detect(term, 'haqer')) %>% 
    head(n = 20)

##############################3
## make figures
###############################
p_iq %>% 
    ggsave(filename = '/wdata/lcasten/sli_wgs/HCA_reversion/figures/SPARK_IQ_rare_HAQER_reversions.png', dpi = 300, units = 'in', device = 'png', width = 15, height = 6)
p_lang_milestone %>% 
    ggsave(filename = '/wdata/lcasten/sli_wgs/HCA_reversion/figures/SPARK_language_milestones_rare_HAQER_reversions.png', dpi = 300, units = 'in', device = 'png', width = 15, height = 6)
p_imp %>% 
    ggsave(filename = '/wdata/lcasten/sli_wgs/HCA_reversion/figures/SPARK_cog_lang_impairments_rare_HAQER_reversions.png', dpi = 300, units = 'in', device = 'png', width = 10, height = 7)


## ==============================================================================
## make forest plot with all phenos and other rare reversion scores as covars
## ==============================================================================
kp_dx <- p_forest_dat_dx %>% 
    filter(term != '(Intercept)') %>% 
    filter(str_detect(name, '_dep|mood_anx|dev_lang_dis|soc_prag|asd|motor|bipol|ocd')) %>% 
    distinct(name)

p_forest <- bind_rows(res_ci[res_ci$name == 'cog_imp',], p_forest_dat_dx[p_forest_dat_dx$name %in% kp_dx$name,], p_forest_age_dat) %>% 
    filter(term != '(Intercept)') %>% 
    arrange(p.value) %>% 
    mutate(pred = case_when(str_detect(term, 'haqer') ~ 'HAQERs',
                            str_detect(term, 'har') ~ 'HARs',
                            str_detect(term, 'rand') ~ 'RAND')) %>% 
    mutate(xmn = estimate - (1.96 * std.error),
           xmx = estimate + (1.96 * std.error)) %>% 
    mutate(type = ifelse(str_detect(name, '_mos'), 'Language milestones', 'Diagnosis')) %>% 
    group_by(type) %>% 
    mutate(maxn = max(n)) %>%
    ungroup() %>%
    mutate(type = str_c(type, ' (N = ', maxn, ')')) %>%
    mutate(pred = factor(pred, levels = c('HAQERs', 'HARs', 'RAND'))) %>%
    arrange(desc(type), pred, estimate) %>%
    mutate(name = factor(name, levels = unique(name))) %>%
    mutate(sig = ifelse(p.value < 0.05, TRUE, FALSE)) %>% 
    distinct(name, pred, .keep_all = TRUE) %>%
    relocate(name, pred, type, estimate, xmn, xmx, n) %>% 
    arrange(desc(type), name, pred, estimate) %>%
    mutate(pred = factor(pred, levels = rev(c('HAQERs', 'HARs', 'RAND')))) %>%
    ggplot(aes(x = estimate, y = name, color = pred, group = pred)) +
    geom_point(aes(x = estimate, y = name, shape = sig), position = position_dodge(width = 0.4), size = 4) +
    geom_linerange(aes(xmin = xmn, xmax = xmx, y = name, color = pred), position = position_dodge2(width = 0.4), size = 1.2) +
    # facet_grid(rows =  vars(type), scales = 'free_y', space = 'free') +
    ggforce::facet_col( ~ type, scales = 'free_y', space = 'free') +
    geom_vline(xintercept = 0, color = 'black', linetype = 'dashed', size = 1.05) +
    scale_color_manual(values = c('#dcc695', '#e04468', '#762776')) +
    scale_shape_manual(values = c(1, 19)) +
    theme_classic() +
    theme(axis.text = element_text(size = 14),
          axis.title = element_text(size = 16),
          strip.text = element_text(size = 16),
          legend.position = 'bottom',
          legend.text = element_text(size = 14),
          legend.title = element_text(size = 16)) +
    xlab('Regression beta (95% CI)') +
    ylab(NULL) +
    labs(shape = 'p < 0.05:', color = 'Rare reversion count:') +
    guides(shape = "none")
des <- "
14
24
34
"
p <- (p_dist_haq / p_dist_har / p_dist_rand) + (p_forest) + 
    plot_layout(design = des, widths = c(0.65, 1)) + 
    plot_annotation(tag_levels = 'A') & theme(plot.tag = element_text(face = 'bold', size = 18), axis.text = element_text(size = 16), axis.title = element_text(size = 18), legend.text = element_text(size = 16), legend.title = element_text(size = 18), strip.text = element_text(size = 18))
p
ggsave(p, filename = '/wdata/lcasten/sli_wgs/paper_figures/SPARK_rare_HAQER_reversions_comparison.png', dpi = 300, units = 'in', device = 'png', width = 16, height = 12)

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