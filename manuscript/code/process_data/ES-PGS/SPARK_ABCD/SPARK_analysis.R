##
library(tidyverse)

# unrel_samples <- read_table('/sdata-new/public-data/HCP/data/derivatives/phg000989.v1.MappingHumanConnectome_Marchini.genotype-imputed-data.c1.GRU-IRB-PUB/autosomes.qc.prs.unrelated.rel.id', col_names = FALSE)

fam <- read_table('/wdata/lcasten/spark_abcd_array_merge/topmed_hg19/autosomes.qc.fam', col_names = FALSE)
# kin <- read_table('/sdata-new/public-data/HCP/data/derivatives/phg000989.v1.MappingHumanConnectome_Marchini.genotype-imputed-data.c1.GRU-IRB-PUB/autosomes.qc.prs.rel.gz', col_names = FALSE)
# kin <- as.matrix(kin)
# colnames(kin) <- fam$X2
# rownames(kin) <- fam$X2
# dim(kin)
# kin[1:5,1:5]

##
ph <- read_csv('/wdata/lcasten/spark/research_match/language_merged_results/data/merged_raw.qc.csv')
names(ph)

# ph <- ph %>% 
#     select(Subject, Gender, fMRI_Lang_PctCompl, CogFluidComp_AgeAdj, CogEarlyComp_AgeAdj, CogCrystalComp_AgeAdj, ReadEng_AgeAdj, PicVocab_AgeAdj, ProcSpeed_AgeAdj, Flanker_AgeAdj, CardSort_AgeAdj, PicSeq_AgeAdj, matches('Language_Task'), matches('_Vol$|_Thck$|Area$'))
# ph %>% 
#     names()


cor.test(ph$PIC.unique_word_count, ph$COWAT.unique_word_count)


###################################
## ================================
## analysis
## ================================
###################################
##
dat <- read_csv('/wdata/lcasten/sli_wgs/prs/pathway_prs/SPARK_ABCD.gathered_pgs_pc_corrected_long_full_data.cogPerf.complement.csv')
unique(dat$pgs_name)

pgs_base_wide <- dat %>% 
#   filter(IID %in% unrel_samples$X2) %>%
#   filter(str_detect(pgs_name, pattern = 'human_evo')) %>% 
  filter(str_detect(pgs_name, pattern = '1$')) %>%
  mutate(pgs = str_split(pgs_name, pattern = '[.]', simplify = TRUE)[,1]) %>%
  select(IID, pgs, pgs_genome_wide_baseline) %>% 
  distinct()
pgs_long_both <- dat %>% 
#   filter(str_detect(IID, pattern = 'sample')) %>% 
#   filter(IID %in% unrel_samples$X2) %>%
#   filter(str_detect(pgs_name, pattern = 'human_evo|chimp')) %>% 
  filter(str_detect(pgs_name, pattern = '1$')) %>%
  select(-pgs_raw) %>%
  select(-pgs_genome_wide_baseline) %>%
  select(-cohort) %>%
  # mutate(pgs = str_split(pgs_name, pattern = '[.]', simplify = TRUE)[,1]) %>%
  mutate(pgs = str_split(pgs_name, pattern = '[.]', simplify = TRUE)[,1],
         gs = str_split(pgs_name, pattern = '[.]', simplify = TRUE)[,2],
        #  anno_tmp = str_split(pgs_name, pattern = '[.]', simplify = TRUE)[,3],
         anno = str_remove_all(gs, pattern = '_5e_08|_0.0005|_0.05|_0.2|_1|_0'),
         thr = str_remove_all(gs, pattern = anno)) %>%
  # select(IID, pgs, anno, pgs_pathway_corrected_for_genome_wide_burden) %>%
  # pivot_wider(id_cols = -matches('anno|pgs_pathway'), names_from = anno, values_from = 'pgs_pathway_corrected_for_genome_wide_burden')
  select(IID, pgs, anno, matches('pgs_pc_corrected'))
pgs_evo_wide <- dat %>% 
#   filter(str_detect(IID, pattern = 'sample')) %>% 
#   filter(IID %in% unrel_samples$X2) %>%
#   filter(str_detect(pgs_name, pattern = 'human_evo|chimp')) %>% 
  filter(str_detect(pgs_name, pattern = '1$')) %>%
  select(-pgs_raw) %>%
  select(-pgs_genome_wide_baseline) %>%
  select(-cohort) %>%
  # mutate(pgs = str_split(pgs_name, pattern = '[.]', simplify = TRUE)[,1]) %>%
  mutate(pgs = str_split(pgs_name, pattern = '[.]', simplify = TRUE)[,1],
         gs = str_split(pgs_name, pattern = '[.]', simplify = TRUE)[,2],
        #  anno_tmp = str_split(pgs_name, pattern = '[.]', simplify = TRUE)[,3],
         anno = str_remove_all(gs, pattern = '_5e_08|_0.0005|_0.05|_0.2|_1|_0'),
         thr = str_remove_all(gs, pattern = anno)) %>%
  # select(IID, pgs, anno, pgs_pathway_corrected_for_genome_wide_burden) %>%
  # pivot_wider(id_cols = -matches('anno|pgs_pathway'), names_from = anno, values_from = 'pgs_pathway_corrected_for_genome_wide_burden')
  select(IID, pgs, anno, matches('pgs_pc_corrected')) %>%
  pivot_wider(id_cols = -matches('anno|pgs_pathway'), names_from = anno, values_from = matches('pgs_pc_corrected'))
pgs_wide <- pgs_base_wide %>% 
  inner_join(pgs_evo_wide)
pgs_long <- pgs_wide %>% 
  pivot_longer(cols = -c(1:2), names_to = 'pgs_nm', values_to = 'pgs_val')
pgs_long
unique(pgs_long$pgs_nm)

##########################
pgs_long %>% 
    group_by(pgs_nm) %>% 
    summarise(min(pgs_val),
              max(pgs_val))
unique(pgs_long$pgs_nm)
pgs_long_qc <- pgs_long %>% 
    filter(str_detect(pgs_nm, pattern = 'genome|UCE|_div_|HAQER|singleton|NeanderthalSelectiveSweep')) %>%
    filter(IID %in% ph$subject_sp_id) %>% # distinct(IID)
    group_by(pgs_nm) %>% 
    mutate(med = median(pgs_val),
           mad = mad(pgs_val)) %>% 
    # filter(abs(pgs_val) <= 3) %>%
    # filter(pgs_val >= med - 2.5 * mad & pgs_val <= med + 2.5 * mad) %>%
    group_by(IID) %>% 
    mutate(n_valid_pgs = n()) %>% 
    ungroup() %>% 
    filter(n_valid_pgs == max(n_valid_pgs))
length(unique(pgs_long_qc$IID))

pgs_long_qc %>%
    group_by(pgs_nm) %>% 
    mutate(pgs_val = scale(pgs_val)[,1]) %>% 
    summarise(min(pgs_val),
              max(pgs_val))


##################################
##

## factor analysis of raw scores
fc <- read_csv('/wdata/lcasten/screener_paper/data/factors/v2/factor_scores.csv')
eur <- read_tsv('/wdata/trthomas/array/merged_2022_ABCD_iWES1_WGS_2-4/master.tsv') %>%
  filter(cluster_allPCs == 1)
pc <- read_csv('/wdata/lcasten/spark/prs/HapMap3_plus/PCA/raw_KING_pca_results.csv')
names(pc)[1] <- 'IID'
ph2 <- read_csv('/wdata/lcasten/screener_paper/data/datasets/v2/final_data/sentence_repetition_bigram_accuracy_summary_scores_merged.csv')
##
res <- ph %>% 
    # filter(subject_sp_id %in% eur$IID) %>%
    left_join(fc) %>%
    pivot_longer(cols = matches('^[A-Z]', ignore.case = FALSE)) %>%
   rename(IID = subject_sp_id) %>% 
    group_by(name) %>% 
    mutate(norm_value = qnorm((rank(value,na.last="keep")-0.5)/sum(!is.na(value)))) %>%
  inner_join(pgs_long_qc) %>%
  inner_join(pc) %>%
  group_by(name, pgs, pgs_nm) %>% 
  do(res = broom::tidy(cor.test(.$norm_value, .$pgs_val, method = 'p'))) %>%
  # do(res = broom::tidy(lm(norm_value ~ pgs_val + age*as.factor(sex) + pc1 + pc2 + pc3 + pc4 + pc5, data = .))) %>%
  unnest(res) %>% 
  # filter(term == 'pgs_val') %>%
  arrange(p.value) %>% 
  mutate(lab = str_c('Pearson r = ', round(estimate, digits = 2), ', p-val = ', formatC(p.value, digits = 2), '\nN = ', parameter + 2))

kp_res <- res %>% 
    filter(str_detect(pgs_nm, pattern = 'genome|complement', negate = TRUE)) %>% 
    filter(str_detect(pgs_nm, 'HAQER|chimp|singleton')) %>% 
    filter(str_detect(name, 'bigram|Factor')) %>% 
    filter(name == 'Factor2' & str_detect(pgs_nm, 'HAQER'))
    # select(name, pgs_nm, statistic, p.value) %>% 

kp_res <- res %>% 
    # filter(str_detect(pgs_nm, pattern = 'genome|complement', negate = TRUE)) %>% 
    filter(name == 'Factor2') 

tmp_gw = pgs_wide %>% 
  select(IID, genome_wide_cog_pgs = pgs_genome_wide_baseline)

kp_res_resid = ph %>% 
    # filter(subject_sp_id %in% eur$IID) %>%
    left_join(fc) %>%
    pivot_longer(cols = matches('^[A-Z]', ignore.case = FALSE)) %>%
   rename(IID = subject_sp_id) %>% 
    group_by(name) %>% 
    mutate(norm_value = qnorm((rank(value,na.last="keep")-0.5)/sum(!is.na(value)))) %>%
    # inner_join(tmp_gw) %>%
    inner_join(pgs_long_both) %>%
    drop_na(name, norm_value, pgs_pc_corrected_complement) %>%
    rename(pgs_nm = anno) %>%
    group_by(name, pgs_nm) %>% 
    mutate(resid_norm_value = resid(lm(norm_value ~ pgs_pc_corrected_complement)),
           resid_norm_value = scale(resid_norm_value)[,1]) %>% 
  inner_join(pc) %>%
  drop_na(name, pgs, pgs_nm) %>%
  group_by(name, pgs, pgs_nm) %>% 
  do(res = broom::tidy(cor.test(.$resid_norm_value, .$pgs_pc_corrected, method = 'p'))) %>%
  # do(res = broom::tidy(lm(norm_value ~ pgs_val + age*as.factor(sex) + pc1 + pc2 + pc3 + pc4 + pc5, data = .))) %>%
  unnest(res) %>% 
  # filter(term == 'pgs_val') %>%
  arrange(p.value) %>% 
  mutate(lab = str_c('Pearson r = ', round(estimate, digits = 2), ', p-val = ', formatC(p.value, digits = 2), '\nN = ', parameter + 2))
kp_res_resid %>% 
  filter(pgs_nm == 'HAQER')
kp_res_resid %>% 
  filter(pgs_nm == 'human_singleton_density_score_top5pct')

##
ph %>% 
    # filter(subject_sp_id %in% eur$IID) %>%
    left_join(fc) %>%
    pivot_longer(cols = matches('^[A-Z]', ignore.case = FALSE)) %>%
   rename(IID = subject_sp_id) %>% 
    group_by(name) %>% 
    mutate(norm_value = qnorm((rank(value,na.last="keep")-0.5)/sum(!is.na(value)))) %>%
  inner_join(pgs_long_qc) %>%
  inner_join(pc) %>% 
  inner_join(kp_res) %>% 
  write_csv('/wdata/lcasten/sli_wgs/prs/replication/SPARK_Lingo_verbal_memory_ES-PGS_data.csv')

ph %>% 
    # filter(subject_sp_id %in% eur$IID) %>%
    left_join(fc) %>%
    pivot_longer(cols = matches('^[A-Z]', ignore.case = FALSE)) %>%
   rename(IID = subject_sp_id) %>% 
    group_by(name) %>% 
    mutate(norm_value = qnorm((rank(value,na.last="keep")-0.5)/sum(!is.na(value)))) %>%
    # inner_join(tmp_gw) %>%
    inner_join(pgs_long_both) %>%
    drop_na(name, norm_value, pgs_pc_corrected_complement) %>%
    rename(pgs_nm = anno) %>%
    group_by(name, pgs_nm) %>% 
    mutate(resid_norm_value = resid(lm(norm_value ~ pgs_pc_corrected_complement)),
           resid_norm_value = scale(resid_norm_value)[,1]) %>% 
  inner_join(pc) %>%
  drop_na(name, pgs, pgs_nm) %>% 
  inner_join(kp_res_resid) %>% 
  ungroup() %>% 
  # filter(name == 'Factor2' & pgs_nm == 'HAQER')
  write_csv('/wdata/lcasten/sli_wgs/prs/replication/SPARK_Lingo_verbal_memory_ES-PGS_genome_wide_corrected_data.csv')

mod_dat <- fc %>% 
  rename(IID = subject_sp_id) %>%
  inner_join(pgs_wide)
m1 <- lm(Factor2 ~ pgs_genome_wide_baseline, data = mod_dat)
m2 <- lm(Factor2 ~ pgs_genome_wide_baseline + pgs_pc_corrected_HAQER, data = mod_dat)
summary(m1)
summary(m2)
anova(m1, m2)

##
pgs_gw <- pgs_long %>% 
  filter(pgs_nm == 'pgs_genome_wide_baseline') %>% 
  select(IID, pgs_cog_gw = pgs_val)
ph_dx <- read_csv('/wdata/lcasten/spark/SparkDataRelease-2024-04-10/SparkDataRelease-2024-04-10/basic_medical_screening-2024-04-10.csv')
ph_dx <- ph_dx %>% 
  select(subject_sp_id, age_at_eval_years, sex, asd, behav_adhd, matches('ld|dev_id|lang|speech')) %>% 
  drop_na(age_at_eval_years, sex) %>%
  filter(asd == TRUE)
table(ph_dx$asd)
ph_dx[is.na(ph_dx)] <- 0
res_ld <- ph_dx %>% 
  mutate(speech_language_dis = case_when(dev_lang == 1 | dev_lang_dis == 1 | dev_speech == 1 ~ 1,
                                        TRUE ~ 0),
         language_dis = case_when(dev_lang == 1 | dev_lang_dis == 1 ~ 1,
                                        TRUE ~ 0),
         cog_or_lang_dis_factor = case_when(dev_id == 1 ~ 'ID',
                                            dev_id == 0 & language_dis == 1 ~ 'Language impairment without ID',
                                            TRUE ~ 'No impairment'),
         cog_or_lang_dis_factor = factor(cog_or_lang_dis_factor, levels = c('No impairment', 'Language impairment without ID', 'ID'))) %>%
  select(-cog_or_lang_dis_factor) %>%
  pivot_longer(cols = -c(1:4)) %>% 
  rename(IID = subject_sp_id) %>%
  inner_join(pgs_long_both) %>% 
  rename(pgs_nm = anno) %>%
  filter(pgs_nm != 'pgs_genome_wide_baseline') %>% 
#   filter(abs(pgs_val) <= 3) %>%
  group_by(name, pgs, pgs_nm) %>% 
  do(res = broom::tidy(glm(value ~ pgs_pc_corrected + pgs_pc_corrected_complement + as.factor(sex)*age_at_eval_years, data = .))) %>% 
  unnest(res) %>%
  filter(term == 'pgs_pc_corrected') %>% 
  arrange(p.value)
res_ld %>% 
  filter(str_detect(pgs_nm, pattern = 'HAQER|chimp|singleton')) %>%# select(1,3,estimate,p.value) %>% as.data.frame()
  filter(p.value < 0.05) %>% 
  select(-c(term, pgs))

wd <- ph_dx %>% 
  # filter(age_at_eval_years >= 18) %>%
  mutate(speech_language_dis = case_when(dev_lang == 1 | dev_lang_dis == 1 | dev_speech == 1 ~ 1,
                                        TRUE ~ 0),
         language_dis = case_when(dev_lang == 1 | dev_lang_dis == 1 ~ 1,
                                        TRUE ~ 0),
          language_dis_no_ci = case_when(language_dis == 1 & dev_id == 0 ~ 1,
                                        TRUE ~ 0),                          
         cog_or_lang_dis_factor = case_when(dev_id == 1 ~ 'ID',
                                            dev_id == 0 & language_dis == 1 ~ 'Language impairment without ID',
                                            TRUE ~ 'No impairment'),
         cog_or_lang_dis_factor = factor(cog_or_lang_dis_factor, levels = c('No impairment', 'Language impairment without ID', 'ID'))) %>% 
  rename(IID = subject_sp_id) %>%
  inner_join(pgs_wide)

names(wd)

##
m1 <- glm(dev_lang ~ pgs_pc_corrected_complement_HAQER + as.factor(sex)*age_at_eval_years, data = wd, family = 'binomial')
summary(m1)
m2 <-  glm(dev_lang ~ pgs_pc_corrected_complement_HAQER + as.factor(sex)*age_at_eval_years + pgs_pc_corrected_HAQER, data = wd, family = 'binomial')
summary(m2)
haqer_li_res <- broom::tidy(m2) %>% 
  filter(term == 'pgs_pc_corrected_HAQER')
anova(m1, m2, test = "Chisq")
bp_lab = str_c('Beta = ', round(haqer_li_res$estimate, digits = 3), ', p-val = ', formatC(haqer_li_res$p.value, digits = 2))

##
wilcox.test(wd$HAQER ~ wd$language_dis_no_ci)
t.test(wd$HAQER ~ wd$language_dis)
# bp_lab <- broom::tidy(t.test(wd$HAQER ~ wd$dev_lang)) %>% 
#   mutate(lab = str_c('t = ', round(statistic, 2), ', p-val = ', formatC(p.value, digits = 2)))
p <- wd %>% 
  group_by(dev_lang) %>% 
  mutate(n = n()) %>% 
  ungroup() %>%
  mutate(ld = str_c(as.logical(dev_lang),'\nN=', n)) %>% 
  ggplot(aes(x = ld, y = HAQER)) +
  geom_violin() +
  geom_boxplot(width = 0.5, aes(fill = ld)) +
  geom_hline(yintercept = 0, color = 'red2', linetype = 'dashed', size = 1.025) +
  geom_text(aes(x = 1.5, y = 3.5, label = bp_lab), size = 3.5, check_overlap = TRUE, color = 'grey30') +
  xlab('Developmental language disorder in ASD cases') +
  ylab('Cognitive performance HAQER ES-PGS') +
  theme_classic() +
  theme(legend.position = 'none') +
  coord_cartesian(ylim = c(-3.8,3.8))

p %>% 
  ggsave(filename = '/wdata/lcasten/sli_wgs/paper_figures/SPARK_language_impairment_replication.png', device = 'png', dpi = 300, units  ='in', width = 7, height = 7)


#####################
## try selection analysis by looking at birth year
byr_sp <- read_csv('/sdata/Simons/SPARK/DATA/phenotypes/SPARK_Collection_v13_DataRelease_2024-09-24/individuals_registration-2024-09-24.csv') %>% 
  select(IID = subject_sp_id, year_reg = registration_year, age_yr = age_at_registration_years)
names(byr_sp)
byr_sp <- byr_sp %>% 
  mutate(birth_year = year_reg - age_yr)
range(byr_sp$birth_year, na.rm = TRUE)

byr_sp %>% 
  inner_join(pgs_long_both) %>% 
  rename(pgs_nm = anno) %>%
  filter(pgs_nm != 'pgs_genome_wide_baseline') %>% 
  #   filter(abs(pgs_val) <= 3) %>%
  mutate(name = 'birth_year') %>%
  group_by(name, pgs, pgs_nm) %>% 
  do(res = broom::tidy(glm(birth_year ~ pgs_pc_corrected + pgs_pc_corrected_complement, data = .))) %>% 
  unnest(res) %>%
  filter(str_detect(term, 'pgs_pc_corrected')) %>% 
  arrange(p.value) %>% 
  filter(str_detect(pgs_nm, 'HAQER')) %>% 
  select(-pgs_nm) %>% 
    mutate(lab = str_c('Beta = ', round(estimate, digits = 2), ', p = ', formatC(p.value, digits = 2))) %>% 
    mutate(term_clean = ifelse(term == 'pgs_pc_corrected', 'HAQER CP-PGS', 'Background CP-PGS')) %>% 
    ggplot(aes(x = estimate, y = term_clean)) +
    geom_point(size = 4) +
    geom_linerange(aes(xmin = estimate - 1.96 * std.error, xmax = estimate + 1.96 * std.error), size = 1.2) +
    xlab('Regression beta for birth year') +
    geom_vline(xintercept = 0, color = 'red', linetype = 'dashed', size = 1.075) +
    ylab(NULL)

#####################
## analysis of IQ
iq <- read_csv('/sdata/Simons/SPARK/DATA/phenotypes/SPARK_Collection_v13_DataRelease_2024-09-24/iq-2024-09-24.csv') %>% 
  select(IID = subject_sp_id, age = age_test_date_months, invalid, matches('_score')) %>%
 filter(age >= 5 * 12) # %>%
#  filter(invalid == 0)

res_iq <- iq %>%
  # drop_na() %>%
  pivot_longer(cols = -c(1,2, 3)) %>% 
  inner_join(pgs_long_both) %>% 
  rename(pgs_nm = anno) %>%
  filter(pgs_nm != 'pgs_genome_wide_baseline') %>% 
  #   filter(abs(pgs_val) <= 3) %>%
  group_by(name, pgs, pgs_nm) %>% 
  do(res = broom::tidy(glm(value ~ pgs_pc_corrected + pgs_pc_corrected_complement, data = .))) %>% 
  unnest(res) %>%
  filter(term == 'pgs_pc_corrected') %>% 
  arrange(p.value)

res_iq %>% 
  filter(str_detect(pgs_nm, pattern = 'HAQER|chimp|singleton')) # %>%# select(1,3,estimate,p.value) %>% as.data.frame()
  filter(p.value < 0.05) %>% 
  select(-c(term, pgs))


##########################
## ES-PGS modeling
## do anova analysis for regressions w/ each of the annos
#############################
pgs_wide

coi <- names(pgs_wide)[str_detect(names(pgs_wide), pattern = 'pgs_pc_corrected') & str_detect(names(pgs_wide), 'complement', negate = TRUE)]

tmp <- fc %>% 
  select(IID = subject_sp_id, Factor2) %>%
  pivot_longer(cols = -c(1)) %>% 
  # rename(IID = sample) %>%
  inner_join(pgs_wide)

iter = 0
res_list = list()
for(c in coi){
  cat(sprintf('\n\n\n\n'))
  message(c)
  for(i in 2:2){
      message('Factor ', i)
      iter = iter + 1
      tmp2 <- tmp %>% 
        filter(name == str_c('Factor', i)) %>% 
        select(IID, value, matches(str_remove_all(c, pattern = 'pgs_pc_corrected_')))
      bdat <- tmp2 %>% 
        select(IID, value, matches('complement'))

      baseline <- lm(value ~ . - IID, data = bdat)
      baseline_plus_anno <- lm(value ~ . - IID, data = tmp2)
      baseline_rsq = summary(baseline)$r.squared
      baseline_plus_anno_rsq = summary(baseline_plus_anno)$r.squared
      baseline_plus_anno_coef <- broom::tidy(baseline_plus_anno) %>% 
        filter(term != '(Intercept)' & str_detect(term, 'complement', negate = TRUE))

      res <- broom::tidy(anova(baseline, baseline_plus_anno)) %>% 
        drop_na() %>% 
        mutate(factor = str_c('Factor', i)) %>% 
        mutate(c = str_remove_all(c, 'pgs_pc_corrected_')) %>% 
        relocate(factor, c) %>% 
        mutate(baseline_rsq = baseline_rsq,
               baseline_plus_anno_rsq = baseline_plus_anno_rsq,
               annotation_beta = baseline_plus_anno_coef$estimate,
               annotation_std_err = baseline_plus_anno_coef$std.error,
               annotation_pval = baseline_plus_anno_coef$p.value)
      res_list[[iter]] <- res
  }
}

## get SNP info
snp_counts <- read_csv('/wdata/lcasten/sli_wgs/prs/pathway_prs/SPARK_ABCD.prsice_info.cogPerf.complement.csv') %>% 
  filter(Threshold == 1) %>% 
  select(anno = Set, Num_SNP)  %>% 
  mutate(model = str_remove_all(anno, pattern = 'complement_')) %>% 
  mutate(complement = ifelse(str_detect(anno, pattern = 'complement_'), 'complement_PGS_n_snp', 'annotation_PGS_n_snp')) %>%
  pivot_wider(id_cols = model, names_from = complement, values_from = Num_SNP)
snp_counts %>% 
  as.data.frame()

##
unique(bind_rows(res_list)$c)
bind_rows(res_list) %>% 
  arrange(p.value) %>% 
  select(-term) %>% 
  rename(model = c, p.value_model_comparison = p.value) %>%
  inner_join(snp_counts) %>% 
  mutate(model = factor(model, levels = c('pgs_genome_wide_baseline', 'consPrimates_65mya', 'consPrimates_UCE', 'LinAR_Simiformes', 'LinAR_Catarrhini', 'LinAR_Hominoidea', 'LinAR_Hominidae', 'LinAR_Homininae', 'LinAR_Human', 'human_chimp_div_DMG', 'HAQER', 'HAR', 'NeanderthalSelectiveSweep', 'human_singleton_density_score_top5pct'))) %>%
  drop_na(model) %>%
  arrange(factor, model) %>%
  relocate(complement_PGS_n_snp, .before = annotation_PGS_n_snp) %>% 
  filter(p.value_model_comparison < 0.05) %>% 
  select(-c(3:6)) %>% 
  relocate(matches('annotation_')) %>% 
  filter(str_detect(factor, '1|2')) %>% 
  filter(annotation_beta > 0) %>% 
  relocate(model, factor)

f2_res <- bind_rows(res_list) %>%
  arrange(p.value) %>% 
  select(-term) %>% 
  rename(model = c, p.value_model_comparison = p.value) %>%
  inner_join(snp_counts) %>% 
  mutate(model = factor(model, levels = c('pgs_genome_wide_baseline', 'consPrimates_65mya', 'consPrimates_UCE', 'LinAR_Simiformes', 'LinAR_Catarrhini', 'LinAR_Hominoidea', 'LinAR_Hominidae', 'LinAR_Homininae', 'LinAR_Human', 'human_chimp_div_DMG', 'HAQER', 'HAR', 'NeanderthalSelectiveSweep', 'human_singleton_density_score_top5pct'))) %>%
  drop_na(model) %>%
  arrange(factor, model) %>%
  relocate(complement_PGS_n_snp, .before = annotation_PGS_n_snp) %>%
  mutate(rsq_gain_per_1000indSNP = (baseline_plus_anno_rsq - baseline_rsq) / (annotation_PGS_n_snp / 1000))

## repeat for binary phenos
tmp <- wd %>% 
  select(IID, behav_adhd, dev_id, dev_lang, dev_lang_dis, dev_ld, dev_speech) %>%
  pivot_longer(cols = -c(1)) %>% 
  # rename(IID = sample) %>%
  inner_join(pgs_wide)

iter = 0
res_list = list()
for(c in coi){
  cat(sprintf('\n\n\n\n'))
  message(c)
  for(p in unique(tmp$name)){
      message(p)
      iter = iter + 1
      tmp2 <- tmp %>% 
        filter(name == p) %>% 
        select(IID, value, matches(str_remove_all(c, pattern = 'pgs_pc_corrected_')))
      bdat <- tmp2 %>% 
        select(IID, value, matches('complement'))

      baseline <- glm(value ~ . - IID, data = bdat, family = 'binomial')
      baseline_plus_anno <- glm(value ~ . - IID, data = tmp2, family = 'binomial')
      baseline_aic = summary(baseline)$aic
      baseline_plus_anno_aic = summary(baseline_plus_anno)$aic
      baseline_plus_anno_coef <- broom::tidy(baseline_plus_anno) %>% 
        filter(term != '(Intercept)' & str_detect(term, 'complement', negate = TRUE))

      res <- broom::tidy(anova(baseline, baseline_plus_anno, test = "Chisq")) %>% 
        drop_na() %>% 
        mutate(factor = p) %>% 
        mutate(c = str_remove_all(c, 'pgs_pc_corrected_')) %>% 
        relocate(factor, c) %>% 
        mutate(baseline_aic = baseline_aic,
               baseline_plus_anno_aic = baseline_plus_anno_aic,
               annotation_beta = baseline_plus_anno_coef$estimate,
               annotation_std_err = baseline_plus_anno_coef$std.error,
               annotation_pval = baseline_plus_anno_coef$p.value)
      res_list[[iter]] <- res
  }
}

li_res <- bind_rows(res_list) %>%
  arrange(p.value) %>% 
  select(-term) %>% 
  rename(model = c, p.value_model_comparison = p.value) %>%
  inner_join(snp_counts) %>% 
  mutate(model = factor(model, levels = c('pgs_genome_wide_baseline', 'consPrimates_65mya', 'consPrimates_UCE', 'LinAR_Simiformes', 'LinAR_Catarrhini', 'LinAR_Hominoidea', 'LinAR_Hominidae', 'LinAR_Homininae', 'LinAR_Human', 'human_chimp_div_DMG', 'HAQER', 'HAR', 'NeanderthalSelectiveSweep', 'human_singleton_density_score_top5pct'))) %>%
  drop_na(model) %>%
  arrange(factor, model) %>%
  relocate(complement_PGS_n_snp, .before = annotation_PGS_n_snp)


##
f2_res %>% 
  # filter(model == 'human_singleton_density_score_top5pct') %>% relocate(matches('annotation'))
  write_csv('/wdata/lcasten/sli_wgs/prs/pathway_prs/results/SPARK_evo_prs_lm_comparison_results.csv')

li_res %>% 
  # filter(model == 'human_singleton_density_score_top5pct') %>% relocate(matches('annotation'))
  write_csv('/wdata/lcasten/sli_wgs/prs/pathway_prs/results/SPARK_evo_prs_glm_comparison_results.csv')






############################
##
library(ggpubr)
summary(aov(HAQER ~ cog_or_lang_dis_factor + as.factor(sex) + age_at_eval_years, data = wd))
table(wd$cog_or_lang_dis_factor)
wd %>% 
  mutate(lfactor = case_when(cog_or_lang_dis_factor == 'No impairment' ~ 'No impairment\nN=10,344',
                             cog_or_lang_dis_factor == 'Language impairment without ID' ~ 'Language impairment without ID\nN=13,797',
                             cog_or_lang_dis_factor == 'ID' ~ 'ID\nN=6092'),
         lfactor = factor(lfactor, levels = c('No impairment\nN=10,344', 'Language impairment without ID\nN=13,797', 'ID\nN=6092'))) %>%
    # pivot_longer(cols = c('pgs_pc_corrected_human_singleton_density_score_top5pct', 'pgs_pc_corrected_human_chimp_div_DMG', 'pgs_pc_corrected_HAQER')) %>%
    # mutate(name2 = case_when(name == 'pgs_pc_corrected_human_singleton_density_score_top5pct' ~ 'PGS in recent selection loci (2-3 Kya)',
    #                         name == 'pgs_pc_corrected_human_chimp_div_DMG' ~ 'PGS in human-chimp divergent genes (12 Mya)', 
    #                         name == 'pgs_pc_corrected_HAQER' ~ 'PGS in HAQERs (600 Kya)'),
    #        name2 = factor(name2, levels = c('PGS in human-chimp divergent genes (12 Mya)', 'PGS in HAQERs (600 Kya)', 'PGS in recent selection loci (2-3 Kya)'))) %>%
    ggboxplot(x = "lfactor", y = "HAQER",
          fill = "lfactor", palette = "OrRd")+
  stat_compare_means(method = "t.test", label = 'p.signif', comparisons = list( c('No impairment\nN=10,344', 'Language impairment without ID\nN=13,797'), c('No impairment\nN=10,344', 'ID\nN=6092'), c('Language impairment without ID\nN=13,797', 'ID\nN=6092')))+ # Add pairwise comparisons p-value
  stat_compare_means(method = "anova") + # method = "anova"
  geom_hline(yintercept = 0, color = 'black', linetype = 'dashed', size = 1.075) +
  theme(legend.position = 'none') +
  xlab(NULL) +
  ylab('Cognitive performance HAQER ES-PGS') +
  theme(axis.title = element_text(size = 15, face = 'bold'),
        strip.text = element_text(size = 15, face = 'bold'))


summary(aov(human_singleton_density_score_top5pct ~ cog_or_lang_dis_factor, data = wd))
summary(aov(human_chimp_div_DMG ~ cog_or_lang_dis_factor, data = wd))

