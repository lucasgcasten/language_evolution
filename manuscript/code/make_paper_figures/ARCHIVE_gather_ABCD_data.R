library(tidyverse)
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

############################################
############################################
## ABCD data
############################################
############################################
## Function to read in abcd data text files, returns a two-element list with the data and description of column names (`dict`)
read_abcd <- function(filename, guess_max = 2000){
  require(readr)
  require(dplyr)
  require(tidyr)
  
  if(!file.exists(filename)){stop('file does not exist.')}
  
  data <- readr::read_tsv(filename, 
                          col_names = readr::read_tsv(filename, n_max = 1, col_names = FALSE, col_types = readr::cols()) %>% as.character(),
                          skip = 2, 
                          guess_max = guess_max, 
                          locale = locale(date_format = "%m/%d/%Y"),
                          na = c('NA', '999' ,'', '.'), 
  ) %>%
    dplyr::mutate_at(dplyr::vars(dplyr::ends_with("_id")), as.character) %>%
    as.data.frame()
  
  dict <- readr::read_tsv(filename, n_max = 1, col_types = readr::cols()) %>%
    tidyr::gather(colname, description) %>%
    as.data.frame()
  
  return(list(data = data, dict = dict))
}

abcd_height <- read_abcd('/Dedicated/jmichaelson-sdata/ABCD/abcd_release_4_0/abcd_ant01.txt')
abcd_height$dict
ht <- abcd_height$data %>% 
  select(subjectkey, interview_age, anthroheightcalc,,anthroweight1lb) %>%
  drop_na() %>%
  as_tibble() %>% 
  arrange(desc(interview_age)) %>%
  # distinct(subjectkey, .keep_all = TRUE) %>%
  rename(IID = subjectkey) %>% 
  filter(anthroheightcalc >= 36)
abcd_icv <- read_abcd("/Dedicated/jmichaelson-sdata/ABCD/abcd_release_4_0/abcd_smrip10201.txt")
# abcd_icv <- read_abcd("/Dedicated/jmichaelson-sdata/ABCD/abcd_release_5_0/abcd_smrip10201.txt")

icv <- abcd_icv$data %>%
  as_tibble() %>%
  select(subjectkey, interview_age, sex, smri_vol_scs_intracranialv, smri_vol_scs_wholeb) %>% 
  arrange(desc(interview_age)) %>% 
  drop_na() %>%
  distinct(subjectkey, .keep_all = TRUE)

vol_bg <- abcd_icv$data %>%
  as_tibble() %>%
  select(subjectkey, interview_age, sex, smri_vol_scs_intracranialv, smri_vol_scs_wholeb, matches('parahpal|caudate|putamen|pallidum|hpus|amygd|_aal|_aar|smri_vol_scs_subcorticalgv')) %>% 
  arrange(desc(interview_age)) %>% 
  drop_na() %>%
  distinct(subjectkey, .keep_all = TRUE) %>% 
  mutate(basal_ganglia_vol = smri_vol_scs_caudatelh + smri_vol_scs_putamenlh + smri_vol_scs_pallidumlh + smri_vol_scs_aal + smri_vol_scs_caudaterh + smri_vol_scs_putamenrh + smri_vol_scs_pallidumrh + smri_vol_scs_aar)
hist(vol_bg$basal_ganglia_vol)


icv_counts <- abcd_icv$data %>%
  as_tibble() %>% 
  group_by(subjectkey) %>% 
  count()

icv_youngest <- abcd_icv$data %>%
  as_tibble() %>%
  select(subjectkey, interview_age, sex, smri_vol_scs_intracranialv, smri_vol_scs_wholeb) %>% 
  arrange(interview_age) %>% 
  drop_na() %>%
  distinct(subjectkey, .keep_all = TRUE)


icv_oldest <- abcd_icv$data %>%
  as_tibble() %>%
  select(subjectkey, interview_age, sex, smri_vol_scs_intracranialv, smri_vol_scs_wholeb) %>% 
  arrange(desc(interview_age)) %>% 
  drop_na() %>%
  distinct(subjectkey, .keep_all = TRUE) 
names(icv_oldest)[-1] <- str_c(names(icv_oldest)[-1], '_latest')

ht_oldest <- ht
names(ht_oldest)[-1] <- str_c(names(ht_oldest)[-1], '_latest')

ht_youngest <- ht

brain_growth_ph <- icv_youngest %>% 
  inner_join(icv_oldest) %>% 
  rename(IID = subjectkey) %>%
  inner_join(ht_oldest) %>%
  inner_join(ht_youngest) %>%
  filter(interview_age < interview_age_latest) %>% 
  filter(smri_vol_scs_intracranialv_latest - smri_vol_scs_intracranialv >= 0) %>% 
  mutate(icv_growth_raw = scale(smri_vol_scs_intracranialv_latest - smri_vol_scs_intracranialv)[,1],
         wb_growth_raw = scale(smri_vol_scs_wholeb_latest - smri_vol_scs_wholeb)[,1],
         age_diff = interview_age_latest - interview_age,
         ht_diff = anthroheightcalc_latest - anthroheightcalc,
         wt_diff = anthroweight1lb_latest - anthroweight1lb) %>% 
  mutate(icv_growth_resid = scale(resid(lm(icv_growth_raw ~ smri_vol_scs_intracranialv_latest + smri_vol_scs_intracranialv + interview_age + as.factor(sex) + age_diff)))[,1],
         wb_growth_resid = scale(resid(lm(wb_growth_raw ~ smri_vol_scs_wholeb_latest + smri_vol_scs_wholeb + interview_age + as.factor(sex) + age_diff)))[,1]) %>%  ##  + ht_diff + wt_diff
  select(IID, matches('growth')) %>% 
  mutate(icv_growth_resid = qnorm((rank(icv_growth_resid,na.last="keep")-0.5)/sum(!is.na(icv_growth_resid))),
         wb_growth_resid = qnorm((rank(wb_growth_resid,na.last="keep")-0.5)/sum(!is.na(wb_growth_resid))))
cor(brain_growth_ph[,-1])

brain_growth_ph %>% 
  write_csv('/wdata/lcasten/sli_wgs/prs/pathway_prs/ABCD_brain_growth_phenos.csv')

pgs_abcd <- read_csv('/wdata/lcasten/sli_wgs/prs/pathway_prs/SPARK_ABCD.gathered_pgs_pc_corrected_long_full_data.cogPerf.complement.csv')
# unrel <- read_table('/genome/SPARK_ABCD_imputed_genotype_merge/autosomes_all.GRM.unrelated_cutoff_0.05.singleton.txt', col_names = FALSE)
# unrel <- read_table('/genome/SPARK_ABCD_imputed_genotype_merge/autosomes_all.GRM.unrelated_cutoff_0.05.singleton.txt', col_names = FALSE)
unrel <- read_table('/wdata/lcasten/spark_abcd_array_merge/autosomes_all.GRM.unrelated_cutoff_0.05.singleton.txt', col_names = FALSE)

length(unique(pgs_abcd$IID))

## association with brain vol and PGS at most recent time point
vol_dat <- icv %>%
  rename(IID = subjectkey) %>%
  mutate(vol = qnorm((rank(smri_vol_scs_intracranialv,na.last="keep")-0.5)/sum(!is.na(smri_vol_scs_intracranialv)))) %>%
  mutate(med = median(smri_vol_scs_intracranialv),
         md = mad(smri_vol_scs_intracranialv)) %>%
  inner_join(ht) %>%
  mutate(resid_vol = resid(lm(smri_vol_scs_intracranialv ~ interview_age + anthroheightcalc + anthroweight1lb)),
         resid_vol = scale(resid_vol)[,1]) %>%
  mutate(vol = qnorm((rank(resid_vol,na.last="keep")-0.5)/sum(!is.na(resid_vol)))) %>%
  inner_join(pgs_abcd) 

vol_dat %>% 
  filter(IID %in% unrel$X2) %>%
#   distinct(IID)
  write_csv('/wdata/lcasten/sli_wgs/ancient_DNA/data/ABCD_unrelated_brain_volumes.csv')

icv_res <- vol_dat %>% 
  filter(IID %in% unrel$X2) %>%
  group_by(pgs_name) %>% # count()
  mutate(age_sq = (interview_age / 12) ^2) %>%
  do(res = broom::tidy(glm(scale(smri_vol_scs_intracranialv)[,1] ~ as.factor(sex) + interview_age + anthroheightcalc + anthroweight1lb + pgs_pc_corrected + pgs_pc_corrected_complement, data = .))) %>%
  # do(res = broom::tidy(cor.test(.$resid_vol, .$pgs_pc_corrected, method = 'p'))) %>%
  unnest(res) %>% 
  # filter(term == 'pgs_pc_corrected') %>% 
  filter(str_detect(term, 'pgs_pc_corrected')) %>% 
  arrange(p.value) %>% 
  mutate(lab = str_c('Beta = ', round(estimate, digits = 2), ', p-val = ', formatC(p.value, digits = 2)))
icv_res %>% 
  filter(str_detect(pgs_name, 'HAQER_1'))

## assoc w/ subcortical regions
vol_dat <- vol_bg %>%
  rename(IID = subjectkey) %>%
  inner_join(ht) %>%
  pivot_longer(cols = matches('parahpal|caudate|putamen|pallidum|hpus|amygd|_aal|_aar|smri_vol_scs_subcorticalgv|basal_gang')) %>% 
  group_by(name) %>%
  mutate(resid_vol = resid(lm(value ~ smri_vol_scs_intracranialv + interview_age + anthroheightcalc + anthroweight1lb)),
         resid_vol = scale(resid_vol)[,1]) %>%
  mutate(vol = qnorm((rank(resid_vol,na.last="keep")-0.5)/sum(!is.na(resid_vol)))) %>%
  inner_join(pgs_abcd) 
unique(vol_dat$name)

vol_dat %>% 
  filter(IID %in% unrel$X2) %>%
#   distinct(IID)
  write_csv('/wdata/lcasten/sli_wgs/ancient_DNA/data/ABCD_unrelated_subcortical_brain_volumes.csv')

subcort_res <- vol_dat %>% 
  filter(IID %in% unrel$X2) %>%
  group_by(name, pgs_name) %>% # count()
  mutate(age_sq = (interview_age / 12) ^2) %>%
  do(res = broom::tidy(glm(scale(value)[,1] ~ smri_vol_scs_intracranialv + anthroheightcalc + anthroweight1lb + as.factor(sex) + interview_age + pgs_pc_corrected + pgs_pc_corrected_complement, data = .))) %>%
  # do(res = broom::tidy(cor.test(.$resid_vol, .$pgs_pc_corrected, method = 'p'))) %>%
  unnest(res) %>% 
  # filter(term == 'pgs_pc_corrected') %>% 
  filter(str_detect(term, 'pgs_pc_corrected')) %>% 
  arrange(p.value) %>% 
  mutate(lab = str_c('Beta = ', round(estimate, digits = 2), ', p-val = ', formatC(p.value, digits = 2)))
subcort_res %>% 
  filter(str_detect(pgs_name, 'HAQER_1')) %>% 
  select(name, term, estimate, p.value) %>% 
  filter(str_detect(term, 'complement', negate = TRUE))


## association with brain vol and PGS at most recent time point
vol_dat <- icv %>%
  rename(IID = subjectkey) %>%
  mutate(vol = qnorm((rank(smri_vol_scs_intracranialv,na.last="keep")-0.5)/sum(!is.na(smri_vol_scs_intracranialv)))) %>%
  mutate(med = median(smri_vol_scs_intracranialv),
         md = mad(smri_vol_scs_intracranialv)) %>%
  inner_join(ht) %>%
  mutate(resid_vol = resid(lm(smri_vol_scs_intracranialv ~ interview_age + anthroheightcalc + anthroweight1lb + as.factor(sex))),
         resid_vol = scale(resid_vol)[,1]) %>%
  mutate(vol = qnorm((rank(resid_vol,na.last="keep")-0.5)/sum(!is.na(resid_vol)))) %>%
  inner_join(pgs_abcd) 

vol_dat %>% 
  filter(IID %in% unrel$X2) %>%
  # distinct(IID)
  write_csv('/wdata/lcasten/sli_wgs/ancient_DNA/data/ABCD_unrelated_brain_volumes.csv')

icv_res <- vol_dat %>% 
  filter(IID %in% unrel$X2) %>%
  group_by(pgs_name) %>% # count()
  mutate(age_sq = (interview_age / 12) ^2) %>%
  do(res = broom::tidy(glm(scale(smri_vol_scs_intracranialv)[,1] ~ as.factor(sex) + interview_age + anthroheightcalc + anthroweight1lb + pgs_pc_corrected + pgs_pc_corrected_complement, data = .))) %>%
  # do(res = broom::tidy(cor.test(.$resid_vol, .$pgs_pc_corrected, method = 'p'))) %>%
  unnest(res) %>% 
  # filter(term == 'pgs_pc_corrected') %>% 
  filter(str_detect(term, 'pgs_pc_corrected')) %>% 
  arrange(p.value) %>% 
  mutate(lab = str_c('Beta = ', round(estimate, digits = 2), ', p-val = ', formatC(p.value, digits = 2)))
icv_res %>% 
  filter(str_detect(pgs_name, 'HAQER_1'))

## assoc with brain growth phenotype over ABCD
vol_dat <- brain_growth_ph %>%
  # inner_join(ht) %>%
  inner_join(pgs_abcd) 

hist(brain_growth_ph$icv_growth_resid)
icv_growth_res <- vol_dat %>% 
  filter(IID %in% unrel$X2) %>% # distinct(IID)
  group_by(pgs_name) %>%
  # mutate(age_sq = (interview_age / 12) ^2) %>%
  do(res = broom::tidy(glm(icv_growth_resid ~  pgs_pc_corrected + pgs_pc_corrected_complement, data = .))) %>%
  # do(res = broom::tidy(cor.test(.$resid_vol, .$pgs_pc_corrected, method = 'p'))) %>%
  unnest(res) %>% 
  filter(str_detect(term, 'pgs_pc_corrected')) %>% 
  arrange(p.value) %>% 
  mutate(lab = str_c('Beta = ', round(estimate, digits = 2), ', p-val = ', formatC(p.value, digits = 2)))
icv_growth_res %>% 
  filter(str_detect(pgs_name, 'HAQER_1'))

vol_dat %>%
  filter(IID %in% unrel$X2) %>%
  write_csv('/wdata/lcasten/sli_wgs/ancient_DNA/data/ABCD_unrelated_brain_growth.csv')


icv_growth = read_csv('/wdata/lcasten/sli_wgs/ancient_DNA/data/ABCD_unrelated_brain_growth.csv') %>%
    select(IID, intracranial_growth_resid = icv_growth_resid) %>% 
    distinct()
icv_df = read_csv('/wdata/lcasten/sli_wgs/ancient_DNA/data/ABCD_unrelated_brain_volumes.csv') %>%
    select(1:4, matches('anthro'), resid_vol) %>% 
    distinct()

tmp_vol <- icv_df %>% 
    left_join(icv_growth) %>% 
    rename(mri_age = interview_age)
tmp_vol

pgs_haq <- pgs_abcd %>% 
    filter(str_detect(pgs_name, 'HAQER_1$')) %>% 
    select(IID, cp_pgs.HAQER = pgs_pc_corrected, cp_pgs.background_HAQER = pgs_pc_corrected_complement)

abcd_dat <- tmp_vol %>% 
    select(IID, mri_age, sex, intracranial_vol_resid = resid_vol, intracranial_growth_resid) %>%
    inner_join(pgs_haq)

broom::tidy(lm(intracranial_vol_resid ~ cp_pgs.HAQER + cp_pgs.background_HAQER, data = abcd_dat))
broom::tidy(lm(intracranial_growth_resid ~ cp_pgs.HAQER + cp_pgs.background_HAQER, data = abcd_dat))

abcd_dat %>% 
    write_csv('manuscript/supplemental_materials/ABCD_brain_imaging_data.csv')

## figures
p_vol <- icv %>%
  rename(IID = subjectkey) %>%
  filter(IID %in% unrel$X2) %>% 
  inner_join(ht) %>%
  mutate(resid_vol = resid(lm(smri_vol_scs_intracranialv ~ interview_age + anthroheightcalc + anthroweight1lb + as.factor(sex))),
         resid_vol = scale(resid_vol)[,1]) %>%
  mutate(vol = qnorm((rank(resid_vol,na.last="keep")-0.5)/sum(!is.na(resid_vol)))) %>%
  inner_join(pgs_abcd) %>%
  filter(str_detect(pgs_name, 'HAQER_1')) %>%
  pivot_longer(cols = matches('pgs_pc_corrected'), names_to = 'term') %>% 
  inner_join(icv_res) %>%
  select(IID, resid_vol, term, value, lab) %>%
  mutate(term_clean = ifelse(term == 'pgs_pc_corrected', 'HAQER CP-PGS', 'Background CP-PGS')) %>%
  ggplot(aes(x = resid_vol, y = value)) +
  # geom_point(alpha = 0.2, size = .5, color = '#762776') +
  geom_smooth(method = 'lm', aes(color = term_clean), size = 1., alpha = 0.155) +
  geom_text(aes(x = -.25, y = .25, label = icv_res$lab[icv_res$pgs_name == 'cogPerf.HAQER_1' & icv_res$term == 'pgs_pc_corrected_complement']), color = 'hotpink3', check_overlap = TRUE, size  = 6) +
  geom_text(aes(x = 2, y = -.25, label = icv_res$lab[icv_res$pgs_name == 'cogPerf.HAQER_1' & icv_res$term == 'pgs_pc_corrected']), color = 'dodgerblue', check_overlap = TRUE, size  = 6) +
  theme_classic() +
  theme(axis.text = element_text(size = 16),
        axis.title = element_text(size = 18),
        legend.position = 'bottom') +
  xlab('Intracranial volume') +
  ylab('CP-PGS') +
  labs(color = NULL)


p_growth <- vol_dat %>%
  filter(IID %in% unrel$X2) %>% 
  # inner_join(ht) %>%
  # inner_join(pgs_abcd) %>%
  filter(str_detect(pgs_name, 'HAQER_1')) %>%
  pivot_longer(cols = matches('pgs_pc_corrected'), names_to = 'term') %>% 
  inner_join(icv_growth_res) %>%
  select(IID, icv_growth_resid, term, value, lab) %>%
  mutate(term_clean = ifelse(term == 'pgs_pc_corrected', 'HAQER CP-PGS', 'Background CP-PGS')) %>%
  ggplot(aes(x = icv_growth_resid, y = value)) +
  geom_smooth(method = 'lm', aes(color = term_clean), size = 1.5, alpha = 0.15) +
  geom_text(aes(x = -.25, y = .15, label = icv_growth_res$lab[icv_growth_res$pgs_name == 'cogPerf.HAQER_1' & icv_growth_res$term == 'pgs_pc_corrected_complement']), color = 'hotpink3', check_overlap = TRUE, size  = 6) +
  geom_text(aes(x = 2, y = -.15, label = icv_growth_res$lab[icv_growth_res$pgs_name == 'cogPerf.HAQER_1' & icv_growth_res$term == 'pgs_pc_corrected']), color = 'dodgerblue', check_overlap = TRUE, size  = 6) +
  theme_classic() +
  theme(axis.text = element_text(size = 16),
        axis.title = element_text(size = 18),
        legend.position = 'bottom') +
  xlab('Adolescent intracranial growth') +
  ylab('CP-PGS') +
  labs(color = NULL)


######################################################################################
#########################
## within family (between full siblings) analysis of birth weight and complications
######################################################################################
sp_espgs <- pgs_abcd

kp <- unrel %>% 
    filter(str_detect(X2, pattern = '^SP', negate = TRUE))
cbcl <- read_abcd('/sdata/ABCD/abcd_release_4_0/abcd_cbcls01.txt')

dhx <- read_abcd('/sdata/ABCD/abcd_release_4_0/dhx01.txt')

cbcl_abcd <- cbcl$data %>% 
    as_tibble() %>%
    select(IID = subjectkey, age = interview_age, sex, matches('_t$')) %>%  
    drop_na() %>%
    arrange(desc(age)) %>% 
    distinct(IID, .keep_all = TRUE) %>% 
    rename(total_prob = cbcl_scr_syn_totprob_t) %>% 
    filter(IID %in% kp$X2)

cbcl_abcd %>% 
    pivot_longer(cols = matches('_t$')) %>% 
    inner_join(sp_espgs)


names( dhx$data)
dhx$dict %>% 
    as_tibble() %>% 
    filter(str_detect(description, 'clamp|oxemia|lang|gest|weeks')) %>% 
    as.data.frame()
dhx$dict %>% 
    as_tibble() %>% 
    filter(str_detect(colname, 'devhx_21_p')) %>% 
    as.data.frame()

dhx <- dhx$data %>% 
    as_tibble() %>%
    select(IID = subjectkey, age = interview_age, sex, birth_weight_lbs, birth_weight_oz, rel_lang = devhx_21_p, csec = devhx_13_3_p, preeclampsia_toxemia = devhx_10c3_p, premature = devhx_12_p) %>%  
    mutate(premature = ifelse(is.na(premature), 0, premature),
           #birth_weight_oz = ifelse(is.na(birth_weight_oz), 0, birth_weight_oz)
           ) %>% # relocate(premature)    %>% filter(premature != 0)
    drop_na() %>%
    arrange(desc(age)) %>% 
    distinct(IID, .keep_all = TRUE)

pre_res <- dhx %>%
    # filter(IID %in% kp$X2) %>%
    mutate(bw = birth_weight_lbs * 16 + birth_weight_oz) %>%
  pivot_longer(cols = -c(1:3, 9) ,names_to = 'birth', values_to = 'birth_val') %>% 
  filter(birth_val < 999) %>%
  group_by(birth) %>% 
#   mutate(resid_val = scale(resid(lm(cbcl_val ~ t)))[,1]) %>%
  ungroup() %>%
  inner_join(sp_espgs) %>% 
  group_by(birth, pgs_name) %>%
    do(res = broom::tidy(glm(premature ~ pgs_pc_corrected + pgs_pc_corrected_complement + age + as.factor(sex), data = .))) %>% 
    unnest(res) %>% 
    # filter(term == 'pgs_pc_corrected') %>% 
    filter(str_detect(term, 'pgs_')) %>%
    arrange(p.value) %>% 
    filter(pgs_name == 'cogPerf.HAQER_1') %>%
    filter(birth == 'bw') %>% 
    mutate(lab = str_c('Beta = ', round(estimate, 2), ', p = ', formatC(p.value, digits = 2))) %>% 
    mutate(term = ifelse(term == 'pgs_pc_corrected', 'HAQER CP-PGS', 'Background CP-PGS'))


bw_res <- dhx %>%
    # filter(IID %in% kp$X2) %>%
    mutate(bw = birth_weight_lbs * 16 + birth_weight_oz) %>%
  filter(bw >= 64) %>%
  pivot_longer(cols = -c(1:3, 9) ,names_to = 'birth', values_to = 'birth_val') %>% 
  filter(birth_val < 999) %>%
  group_by(birth) %>% 
#   mutate(resid_val = scale(resid(lm(cbcl_val ~ t)))[,1]) %>%
  ungroup() %>%
  inner_join(sp_espgs) %>% 
  group_by(birth, pgs_name) %>%
    do(res = broom::tidy(glm(birth_val ~ pgs_pc_corrected + pgs_pc_corrected_complement + age*as.factor(sex) + premature, data = .))) %>% 
    unnest(res) %>% 
    # filter(term == 'pgs_pc_corrected') %>% 
    filter(str_detect(term, 'pgs_')) %>%
    arrange(p.value) %>% 
    filter(pgs_name == 'cogPerf.HAQER_1') %>%
    filter(birth == 'bw') %>% 
    mutate(lab = str_c('Beta = ', round(estimate, 2), ', p = ', formatC(p.value, digits = 2))) %>% 
    mutate(term = ifelse(term == 'pgs_pc_corrected', 'HAQER CP-PGS', 'Background CP-PGS'))

p <- dhx %>% 
    inner_join(sp_espgs) %>% 
    filter(str_detect(pgs_name, 'HAQER')) %>% 
    mutate(bw = birth_weight_lbs * 16 + birth_weight_oz) %>%
    mutate(bw_resid = scale(resid(lm(bw ~ age + as.factor(sex))))) %>%
    pivot_longer(cols = c('pgs_pc_corrected', 'pgs_pc_corrected_complement')) %>%
    mutate(term = ifelse(name == 'pgs_pc_corrected', 'HAQER CP-PGS', 'Background CP-PGS')) %>%
    inner_join(bw_res) %>%
    mutate(term = str_remove_all(term, pattern = ' CP-PGS')) %>%
    ggplot(aes(x = value, y = bw_resid, color = term)) +
    geom_smooth(method = 'lm', size = 1.5, alpha = 0.15) +
    ylab('Birth weight zscore') +
    xlab('CP-PGS') +
    labs(color = NULL) +
    geom_text(data = bw_res[bw_res$term == 'HAQER CP-PGS', ], aes(x = 0, y = .15, label = lab), color = 'dodgerblue', size  = 4) +
    geom_text(data = bw_res[bw_res$term != 'HAQER CP-PGS', ], aes(x = 0, y = -0.15, label = lab), color = 'red', size  = 4)
# p %>% 
#     ggsave(filename = '/wdata/lcasten/sli_wgs/paper_figures/subplots/ABCD_birthweight_espgs.png', device = 'png', units = 'in', dpi = 300, width = 6, height = 6)

#####################################
## within family sibling analysis of bw
## does the sibling with higher HAQER CP-PGS have higher birthweight
eur <- read_tsv('/wdata/trthomas/array/merged_2022_ABCD_iWES1_WGS_2-4/master.tsv') %>% 
    filter(cluster_allPCs == 1)
# sib <- read_table('/genome/SPARK_ABCD_imputed_genotype_merge/autosomes_all.GRM.unrelated_cutoff_0.05.family.txt', col_names = FALSE)
sib <- read_table('/wdata/lcasten/spark_abcd_array_merge/autosomes_all.GRM.unrelated_cutoff_0.05.family.txt', col_names = FALSE)
sib <- sib %>% 
    filter(str_detect(X2, 'SP', negate = TRUE) & str_detect(X4, 'SP', negate = TRUE)) %>% 
    filter(X5 >= 0.35 & X5 <= .65)
sib$FID <- 1:nrow(sib)
sib_long <- sib %>% 
    select(FID, X2, X4) %>% 
    pivot_longer(cols = -FID) %>% 
    mutate(sibn = ifelse(name == 'X2', 'sib1', 'sib2')) %>% 
    rename(IID = value) %>% 
    select(FID, IID, sibn)
sib_long %>% 
    group_by(IID) %>% 
    count() %>% 
    ungroup() %>% 
    arrange(desc(n))

sib_long <- sib_long %>% 
    distinct(IID, .keep_all = TRUE)

length(unique(sib_long$IID))

haq_sp_espgs <- sp_espgs %>% 
    filter(str_detect(pgs_name, 'HAQER_1$')) %>%
    select(IID, haqer_cp_pgs = pgs_pc_corrected, background_cp_pgs = pgs_pc_corrected_complement)

sib_dhx <- dhx %>% 
    # filter(IID %in% eur$IID) %>%
    inner_join(haq_sp_espgs) %>%
    mutate(bw = birth_weight_lbs * 16 + birth_weight_oz) %>% 
    mutate(sex_female = ifelse(sex == 'F', 1, 0)) %>% 
    select(-sex) %>%
    pivot_longer(cols = -IID) %>%
    inner_join(sib_long) %>% 
    mutate(name = str_c(sibn, '_', name)) %>%
    pivot_wider(id_cols = FID) # %>% 
    # filter(sib1_premature == 0 & sib2_premature == 0)

cor.test(sib_dhx$sib1_haqer_cp_pgs, sib_dhx$sib2_haqer_cp_pgs)
cor.test(sib_dhx$sib1_background_cp_pgs, sib_dhx$sib2_background_cp_pgs)

p_sib_haq <- sib_dhx %>% 
  ggplot(aes(x = sib1_haqer_cp_pgs, sib2_haqer_cp_pgs)) +
  geom_point() +
  geom_smooth(method = 'lm', alpha = 0.15) +
  xlab('Sibling 1 HAQER CP-PGS') +
  ylab('Sibling 2 HAQER CP-PGS')

p_sib_bg <- sib_dhx %>% 
  ggplot(aes(x = sib1_background_cp_pgs, sib2_background_cp_pgs)) +
  geom_point() +
  geom_smooth(method = 'lm', alpha = 0.15) +
  xlab('Sibling 1 Background CP-PGS') +
  ylab('Sibling 2 Background CP-PGS')  

# p_sib_haq + p_sib_bg

sib_dhx %>% 
    filter(abs(sib1_age - sib2_age) > 0) %>% 
    select(matches('_pgs')) %>% 
    cor()
sib_dhx %>% 
    filter(abs(sib1_age - sib2_age) > 0) %>% 
    select(matches('_pgs')) %>% 
    mutate(haq_diff = sib2_haqer_cp_pgs - sib1_haqer_cp_pgs,
           gw_diff = sib2_background_cp_pgs - sib1_background_cp_pgs) %>% 
    select(matches('diff')) %>% 
    ggplot(aes(x = gw_diff)) +
    geom_histogram()

no_csec <- sib_dhx %>% 
    filter(sib1_csec == 0 & sib2_csec == 0)
both_csec <- sib_dhx %>% 
    filter(sib1_csec == 1 & sib2_csec == 1)

sib_diff_dhx <- sib_dhx %>% 
    pivot_longer(cols = matches('sib1_'), names_to = 'sib1_name', values_to = 'sib1_val') %>% 
    pivot_longer(cols = matches('sib2_'), names_to = 'sib2_name', values_to = 'sib2_val') %>% 
    mutate(name1 = str_remove_all(sib1_name, 'sib1_'),
           name2 = str_remove_all(sib2_name, 'sib2_')) %>% 
    filter(name1 == name2) %>% 
    mutate(sib_diff = sib2_val - sib1_val) %>% 
    select(FID, name1, sib_diff)
length(unique(sib_diff_dhx$FID))


sib_pgs <- sib_diff_dhx %>% 
    filter(str_detect(name1, 'cp_pgs')) %>% 
    pivot_wider(id_cols = FID, names_from = name1, values_from = sib_diff)   
hist(sib_pgs$background_cp_pgs)

sib_cov <-  sib_diff_dhx %>% 
    filter(str_detect(name1, 'age|sex|premat')) %>% 
    pivot_wider(id_cols = FID, names_from = name1, values_from = sib_diff)   

##
bw_n <- sib_diff_dhx %>% 
    filter(str_detect(name1, 'cp_pgs', negate = TRUE)) %>% 
    inner_join(sib_pgs) %>% 
    inner_join(sib_cov) %>%
    # filter(age == 0 & sex_female == 0) %>% # distinct(FID) ## same sex fraternal twins
    filter(age != 0) %>% 
    filter(name1 == 'bw') %>% 
    select(sib_diff, haqer_cp_pgs, background_cp_pgs, premature, sex_female) %>% 
    drop_na() %>% 
    nrow()
    
bw_res <- sib_diff_dhx %>% 
    filter(str_detect(name1, 'cp_pgs', negate = TRUE)) %>% 
    inner_join(sib_pgs) %>% 
    inner_join(sib_cov) %>%
    # filter(age == 0 & sex_female == 0) %>% # distinct(FID) ## same sex fraternal twins
    filter(age != 0) %>% # distinct(FID) ## only siblings, no twins
    group_by(name1) %>% 
    do(res = broom::tidy(lm(sib_diff ~ haqer_cp_pgs + background_cp_pgs + premature + as.factor(sex_female) + age, data = .))) %>% 
    # do(res = broom::tidy(cor.test(.$sib_diff, .$haqer_cp_pgs , method = 's'))) %>% 
    unnest(res) %>% 
    filter(term != '(Intercept)') %>%
    arrange(p.value) %>% 
    filter(str_detect(term, 'pgs')) %>% 
    mutate(lab = str_c('Beta = ', round(estimate, digits = 2), ', p = ', formatC(p.value, digits = 2), '\nN sibling pairs = ', bw_n)) %>% 
    filter(name1 == 'bw') %>% 
    filter(term == 'haqer_cp_pgs')

p_bw <- sib_diff_dhx %>% 
    filter(str_detect(name1, 'cp_pgs', negate = TRUE)) %>% 
    inner_join(sib_pgs) %>% 
    inner_join(sib_cov) %>%
    # filter(age == 0 & sex_female == 0) %>% # distinct(FID) ## same sex fraternal twins
    filter(age != 0) %>% 
    inner_join(bw_res) %>%
    distinct() %>%
    ggplot(aes(x = scale(sib_diff)[,1], y = haqer_cp_pgs)) +
    geom_point(size = 0.7) +
    geom_smooth(method = 'lm', size = 1.5, color = 'dodgerblue') +
    geom_text(aes(x = 0, y = 3.4, label = lab), size = 6, check_overlap = TRUE) +
    xlab('Sibling difference in birth weight') +
    ylab('Sibling difference in\nHAQER CP-PGS') +
    theme(legend.position = 'none')


## save sibling difference data
both_csec %>% 
  write_csv('/wdata/lcasten/sli_wgs/ancient_DNA/data/ABCD_sibling_pairs_both_born_by_csection.csv')

sib_long %>% 
  pivot_wider(id_cols = 1, names_from = 'sibn', values_from = IID) %>% 
  inner_join(sib_diff_dhx) %>% 
# sib_diff_dhx %>% 
    filter(str_detect(name1, 'cp_pgs', negate = TRUE)) %>% 
    inner_join(sib_pgs) %>% 
    inner_join(sib_cov)  %>% # distinct(FID)
    rename(age_diff = age, premature_diff = premature, sex_female_diff = sex_female, pheno = name1, sib_diff_cp_pgs.HAQER = haqer_cp_pgs, sib_diff_cp_pgs.background_HAQER = background_cp_pgs) %>%
    filter(age_diff != 0) %>% # group_by(pheno) %>% do(res = broom::tidy(cor.test(.$sib_diff, .$sib_diff_cp_pgs.HAQER))) %>% unnest(res) %>% arrange(p.value)
    write_csv('/wdata/lcasten/sli_wgs/ancient_DNA/data/ABCD_sibling_birth_phenotypes.csv')

sib_long %>% 
  pivot_wider(id_cols = 1, names_from = 'sibn', values_from = IID) %>% 
  inner_join(sib_diff_dhx) %>% 
# sib_diff_dhx %>% 
    filter(str_detect(name1, 'cp_pgs', negate = TRUE)) %>% 
    inner_join(sib_pgs) %>% 
    inner_join(sib_cov)  %>% # distinct(FID)
    rename(age_diff = age, premature_diff = premature, sex_female_diff = sex_female, pheno = name1, sib_diff_cp_pgs.HAQER = haqer_cp_pgs, sib_diff_cp_pgs.background_HAQER = background_cp_pgs) %>%
    filter(age_diff != 0) %>%
    write_csv('manuscript/supplemental_materials/ABCD_sibling_birth_phenotypes_long.csv')


read_csv('manuscript/supplemental_materials/ABCD_sibling_birth_phenotypes_long.csv') %>% 
    filter(pheno == 'bw') %>% 
    group_by(pheno) %>%
    do(res = broom::tidy(lm(sib_diff ~ sib_diff_cp_pgs.HAQER + sib_diff_cp_pgs.background_HAQER + age_diff + as.factor(sex_female_diff) + premature_diff, data = .))) %>% 
    unnest(res) %>% 
    filter(term == 'sib_diff_cp_pgs.HAQER')



## do logit for csec and complications
csec_res <- sib_diff_dhx %>% 
    filter(str_detect(name1, 'cp_pgs', negate = TRUE)) %>% 
    inner_join(sib_pgs) %>% 
    inner_join(sib_cov) %>%
    filter(str_detect(name1, 'csec|preec')) %>% 
        mutate(sn = sign(sib_diff),
           sn = ifelse(sn == 0, 1, sn),
           sib_diff = abs(sib_diff),
           haqer_cp_pgs = sn * haqer_cp_pgs,
           background_cp_pgs = sn * background_cp_pgs) %>% 
    # filter(age == 0 & sex_female == 0) %>% # distinct(FID) ## same sex fraternal twins
    filter(age != 0) %>% # distinct(FID) ## only siblings, no twins
    # group_by(name1, sib_diff) %>% count() 
    filter(! FID %in% both_csec$FID) %>% ## remove sibling pairs where both were born w/o c-section
    filter(premature == 0) %>% ## remove siblings born at different gestational ages
    # filter(! FID %in% no_csec$FID) %>%
    group_by(name1) %>%
    do(res = broom::tidy(glm(sib_diff ~ haqer_cp_pgs + background_cp_pgs + premature + as.factor(abs(sex_female)), data = ., family = 'binomial'))) %>% 
    # do(res = broom::tidy(cor.test(.$sib_diff, .$haqer_cp_pgs , method = 's'))) %>% 
    unnest(res) %>% 
    filter(term != '(Intercept)') %>%
    arrange(p.value) %>% 
    filter(str_detect(term, 'pgs')) %>% 
    mutate(lab = str_c('Beta = ', round(estimate, digits = 2), ', p = ', formatC(p.value, digits = 2))) %>% 
    filter(name1 == 'csec')

p_csec <- sib_diff_dhx %>% 
    filter(str_detect(name1, 'cp_pgs', negate = TRUE)) %>% 
    inner_join(sib_pgs) %>% 
    inner_join(sib_cov) %>%
    filter(str_detect(name1, 'csec|preec')) %>% 
    mutate(sn = sign(sib_diff),
           sn = ifelse(sn == 0, 1, sn),
           sib_diff = abs(sib_diff),
           haqer_cp_pgs = sn * haqer_cp_pgs,
           background_cp_pgs = sn * background_cp_pgs) %>% 
    # select(name1, sib_diff, haqer_cp_pgs, background_cp_pgs, sn) %>%
    # filter(age == 0 & sex_female == 0) %>% # distinct(FID) ## same sex fraternal twins
    filter(age != 0) %>% 
    inner_join(csec_res) %>%
    filter(term == 'haqer_cp_pgs') %>%
    filter(! FID %in% both_csec$FID) %>% ## remove sibling pairs where both were born w/o c-section
    group_by(name1, sib_diff) %>% 
    mutate(n = n()) %>% 
    # mutate(val = str_c(as.logical(sib_diff), '\nN=', n)) %>% 
    mutate(sib_val = ifelse(sib_diff == 1, 'Only 1 sibling', 'Neither')) %>% 
    # mutate(val = str_c(as.logical(sib_diff), '\nN=', n)) %>% 
    mutate(val = str_c(sib_val, '\nN=', n)) %>% 
    ggplot(aes(x = val, y = haqer_cp_pgs)) +
    geom_violin(size = 1.5) +
    geom_boxplot(size = 1.5, width = .4, aes(fill = val)) +
    geom_hline(yintercept = 0, color = 'red', linetype = 'dashed', size = 1.1) + 
    geom_text(aes(x = 1.5, y = 3, label = lab), size =5, check_overlap = TRUE) +
    xlab('Born via c-section') +
    ylab('Sibling difference in HAQER CP-PGS') +
    theme(legend.position = 'none')

# p <- p_bw + p_csec
# ggsave(p, filename = '/wdata/lcasten/sli_wgs/paper_figures/subplots/ABCD_between_sibling_pgs_birth_complications.png', device = 'png', units = 'in', dpi = 300, width = 14, height = 7, bg = 'white')

sib_dhx %>% 
  write_csv('/wdata/lcasten/sli_wgs/ancient_DNA/data/ABCD_sibling_birth_pheno_raw.csv')

## compare csec sib to non-csec sibling
sib_csec_wd <- sib_dhx %>% 
    mutate(csec_sum = sib1_csec + sib2_csec) %>% 
    filter(csec_sum == 1) %>% 
    filter(abs(sib1_age - sib2_age) > 0) %>% 
    # filter(sib2_sex_female - sib1_sex_female == 0) %>%
    # filter(abs(sib1_premature - sib2_premature) <= 5) %>% ## remove siblings born at different gestational ages
    inner_join(sib_cov) %>%
    mutate(csec_haqer_pgs = ifelse(sib1_csec == 1, sib1_haqer_cp_pgs, sib2_haqer_cp_pgs),
           non_csec_haqer_pgs = ifelse(sib1_csec == 0, sib1_haqer_cp_pgs, sib2_haqer_cp_pgs),
           csec_bg_pgs = ifelse(sib1_csec == 1, sib1_background_cp_pgs, sib2_background_cp_pgs),
           non_csec_bg_pgs = ifelse(sib1_csec == 0, sib1_background_cp_pgs, sib2_background_cp_pgs),
           diff_sex = abs(sex_female),
           diff_premature = sib1_premature - sib2_premature) %>% 
    select(FID, csec_haqer_pgs, non_csec_haqer_pgs, csec_bg_pgs, non_csec_bg_pgs, diff_sex, diff_premature)


## paired analysis of siblings and HAQER CP-PGS

paired_csec <- broom::tidy(t.test(sib_csec_wd$non_csec_haqer_pgs, sib_csec_wd$csec_haqer_pgs, paired = TRUE)) %>% 
    mutate(lab = str_c('t-statistic = ', round(abs(statistic), 2), ', p-val = ', formatC(p.value, digits = 2), '\nN=', parameter + 1, ' sibling pairs'))

p_paired_haq <- sib_csec_wd %>% 
    pivot_longer(cols = c(non_csec_haqer_pgs, csec_haqer_pgs)) %>% 
    mutate(type = ifelse(name == 'non_csec_haqer_pgs', 'Normal delivery', 'C-section')) %>% 
    mutate(type = factor(type, levels = c( 'Normal delivery', 'C-section'))) %>%
    ggplot(aes(x = type, y = value)) +
    geom_violin(size = 1.1) +
    geom_boxplot(width = 0.4, aes(fill = type), size = 1.1, alpha = 0.7) +
    geom_point() +
    geom_line(aes(group = FID), color = 'grey60', alpha = 0.8) +
    geom_hline(yintercept = 0, color = 'red', linetype = 'dashed', size = 1.1) +
    ylab('HAQER CP-PGS') +
    xlab('Birth delivery method') +
    geom_text(aes(x = 1.5, y = 2.7, label = paired_csec$lab), size  = 6, check_overlap = TRUE) +
    theme(legend.position = 'none') +
    scale_fill_manual(values = c('grey70', 'chocolate1'))

## paired analysis of siblings and HAQER CP-PGS
csec_wd_paired <- read_csv('manuscript/supplemental_materials/ABCD_sibling_birth_phenotypes_long.csv') %>% 
    filter(pheno == 'csec') %>%
    filter(abs(sib_diff) == 1) %>% 
    select(FID, sib1, sib2) %>% 
    inner_join(sib_csec_wd) %>% 
    rename(csec_sib_cp_pgs.HAQER = csec_haqer_pgs, non_csec_sib_cp_pgs.HAQER = non_csec_haqer_pgs, csec_sib_cp_pgs.background_HAQER = csec_bg_pgs, non_csec_sib_cp_pgs.background_HAQER = non_csec_bg_pgs) %>% 
    select(1:7)
t.test(csec_wd_paired$csec_sib_cp_pgs.HAQER, csec_wd_paired$non_csec_sib_cp_pgs.HAQER, paired = TRUE)
t.test(csec_wd_paired$csec_sib_cp_pgs.background_HAQER, csec_wd_paired$non_csec_sib_cp_pgs.background_HAQER, paired = TRUE)
csec_wd_paired %>% 
    write_csv('manuscript/supplemental_materials/ABCD_sibling_csection_data.csv')

paired_csec_bg <- broom::tidy(t.test(sib_csec_wd$non_csec_bg_pgs, sib_csec_wd$csec_bg_pgs, paired = TRUE)) %>% 
    mutate(lab = str_c('t-statistic = ', round(abs(statistic), 2), ', p-val = ', formatC(p.value, digits = 2), '\nN=', parameter + 1, ' sibling pairs'))

p_paired_bg <- sib_csec_wd %>% 
    pivot_longer(cols = c(non_csec_bg_pgs, csec_bg_pgs)) %>% 
    mutate(type = ifelse(name == 'non_csec_bg_pgs', 'Normal delivery', 'C-section')) %>% 
    mutate(type = factor(type, levels = c( 'Normal delivery', 'C-section'))) %>%
    ggplot(aes(x = type, y = scale(value)[,1])) +
    geom_violin(size = 1.1) +
    geom_boxplot(width = 0.4, aes(fill = type), size = 1.1, alpha = 0.7) +
    geom_point() +
    geom_line(aes(group = FID), color = 'grey60', alpha = 0.8) +
    geom_hline(yintercept = 0, color = 'red', linetype = 'dashed', size = 1.1) +
    ylab('Background CP-PGS') +
    xlab('Birth delivery method') +
    geom_text(aes(x = 1.5, y = 2.7, label = paired_csec_bg$lab), size  = 6, check_overlap = TRUE) +
    theme(legend.position = 'none') +
    scale_fill_manual(values = c('grey70', 'chocolate1'))


