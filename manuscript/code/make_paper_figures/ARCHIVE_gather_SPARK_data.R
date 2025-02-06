library(tidyverse)

###################################
## get unrel Europeans in SPARK + ABCD for LD calculations in ES-PGS with all 160k samples
eur1 <- read_tsv('/wdata/trthomas/array/merged_2022_ABCD_iWES1_WGS_2-4/master.tsv') %>% 
    filter(cluster_allPCs == 1) %>% 
    select(IID) %>% 
    unlist() %>% 
    unname()
eur2 <- read_tsv('/genome/SPARK_WES/iWES3/metadata/SPARK.iWES_v3.2024_08.ancestry.tsv') %>% 
    filter() %>% 
    select() %>% 
    unlist() %>% 
    unname()
eur2
eur <- unique(c(eur1, eur2))
unrel <- read_table('/Dedicated/jmichaelson-wdata/lcasten/SPARK_ABCD_merge_v2/PCA/unrelated_samples.txt', col_names = FALSE)

unrel %>% 
    filter(X2 %in% eur) %>% 
    write_delim(, col_names = FALSE, delim = ' ')

#######################################################
## gather data for 88k samples in v1 of ES-PGS
#######################################################
espgs <- read_csv('/wdata/lcasten/sli_wgs/prs/pathway_prs/SPARK_ABCD.gathered_pgs_pc_corrected_long_full_data.cogPerf.complement.csv')
dat <- espgs
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
names(pgs_wide) <- str_replace_all(names(pgs_wide), pattern = 'pgs_pc_corrected_', 'cp_pgs.')
names(pgs_wide) <- str_replace_all(names(pgs_wide), pattern = 'complement', 'background')
names(pgs_wide)[3] <- 'cp_pgs.genome_wide'
pgs_wide <- pgs_wide[,-2]

###########################
## rare reversions
###########################
reversions <- read_csv('/wdata/lcasten/sli_wgs/HAQER_variant_associations/data/SPARK_WGS/HAQER_HAR_RAND_10Kb_flank_rare_variant_burden_counts.csv') %>% 
    select(IID = spid, matches('reversion')) %>% 
    filter(str_detect(IID, pattern = '[:]', negate = TRUE))
reversions
reversions <- reversions %>% 
    pivot_longer(cols = -1) %>% 
    group_by(name) %>% 
    mutate(med = median(value),
           md = mad(value)) %>% 
    filter(value >= med - 2.5 * md & value <= med + 2.5 * md) %>% 
    pivot_wider(id_cols = IID) %>% 
    drop_na()

#################
## get phenos
#################
list.files('/sdata/Simons/SPARK/DATA/phenotypes/SPARK_Collection_v13_DataRelease_2024-09-24/')

## basic covariates
core <- read_csv('/sdata/Simons/SPARK/DATA/phenotypes/SPARK_Collection_v13_DataRelease_2024-09-24/individuals_registration-2024-09-24.csv')
core_cov <- core %>% 
    select(IID = subject_sp_id, sex, age_years = age_at_registration_years, asd) %>% 
    drop_na()
asd <- core$subject_sp_id[core$asd == TRUE]

##
cimp <- read_csv('/sdata/Simons/SPARK/DATA/phenotypes/SPARK_Collection_v13_DataRelease_2024-09-24/core_descriptive_variables-2024-09-24.csv')
cimp <- cimp$subject_sp_id[cimp$cognitive_impairment_latest == TRUE]

## milestones
ms <- read_csv('/sdata/Simons/SPARK/DATA/phenotypes/SPARK_Collection_v13_DataRelease_2024-09-24/background_history_child-2024-09-24.csv') %>% 
    select(IID = subject_sp_id, matches('_mos')) %>% 
    pivot_longer(cols = matches('_mos')) %>% 
    mutate(value = ifelse(value == 888, NA_real_, value)) %>% 
    pivot_wider(id_cols = IID)

## sent rep factor from Lingo
fc <- read_csv('/wdata/lcasten/screener_paper/data/factors/v2/factor_scores.csv')
fac_rep <- fc %>% 
    select(IID = subject_sp_id, factor_sent_rep = Factor2)

## diagnoses
dx2 <- read_csv('/sdata/Simons/SPARK/DATA/phenotypes/SPARK_Collection_v13_DataRelease_2024-09-24/basic_medical_screening-2024-09-24.csv')
names(dx2)
addt_dx_res <- dx2 %>% 
    select(spid = subject_sp_id, sex, asd, age_at_eval_years, matches('mood_|dev_|lang|ld|croceph|eating_disorder|neuro_|pers_dis|schiz|sleep_dx|tics|visaud|behav_|birth_')) %>% 
    pivot_longer(cols = -c(spid, sex, asd, age_at_eval_years)) %>% 
    mutate(value = ifelse(is.na(value), 0, value)) %>%
    mutate(name = str_c('dx_', name)) %>% 
    pivot_wider(id_cols = -matches('name|value')) %>% 
    rename(IID = spid)

## IQ
iq <- read_csv('/sdata/Simons/SPARK/DATA/phenotypes/SPARK_Collection_v13_DataRelease_2024-09-24/iq-2024-09-24.csv') %>% 
    select(IID = subject_sp_id, age_test_date_months, matches('_score')) %>% 
    filter(age_test_date_months >= 12 * 5)

##


##############
## MERGE TO ONE MASSIVE SPARK TABLE
spark_df <- core_cov %>% 
    left_join(pgs_wide) %>% 
    left_join(reversions) %>%
    left_join(fac_rep) %>%
    left_join(ms) %>%
    left_join(addt_dx_res) %>%
    left_join(iq) %>% 
    rename(iq_test_age_months = age_test_date_months)
names(spark_df)
colSums(is.na(spark_df))

spark_df %>% 
    filter(is.na(cp_pgs.genome_wide) == FALSE | is.na(rand_rare_variant_reversion_count) == FALSE) %>% ## subset to people with either PGS or WGS reversion data
    write_csv('manuscript/supplemental_materials/SPARK_data.csv')
