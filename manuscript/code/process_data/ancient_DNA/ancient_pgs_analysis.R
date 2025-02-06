library(tidyverse)
library(patchwork)

##
source('/wdata/lcasten/functions/simple_ggplot_theme_presentation.R')

##########################################
## sample metadata
smeta <- read_tsv('/Dedicated/jmichaelson-wdata/lcasten/tools/ref_data/allen_ancient_dna_resource/v54.1.p1_1240K_public/v54.1.p1_1240K_public.anno') %>%
  janitor::clean_names() %>%
  rename(IID = genetic_id) %>% 
  select(IID, group_id, locality, political_entity, lat, long, molecular_sex, sample_age_years_before_1950 = date_mean_in_bp_in_years_before_1950_ce_ox_cal_mu_for_a_direct_radiocarbon_date_and_average_of_range_for_a_contextual_date, age_at_death = age_at_death_from_physical_anthropology)
# smeta <- read_tsv('/Dedicated/jmichaelson-wdata/lcasten/tools/ref_data/allen_ancient_dna_resource/v62/v62.0_1240k_public.anno') %>%
#     janitor::clean_names() %>%
#     rename(IID = 1) %>% 
#     select(IID, group_id, locality, political_entity, lat, long, molecular_sex, sample_age_years_before_1950 = date_mean_in_bp_in_years_before_1950_ce_ox_cal_mu_for_a_direct_radiocarbon_date_and_average_of_range_for_a_contextual_date, age_at_death = age_at_death_morphological_sex_from_physical_anthropology)

names(smeta)
range(smeta$sample_age_years_before_1950)
smeta %>% 
  filter(sample_age_years_before_1950 > 0) %>% 
  distinct(group_id) %>% as.data.frame()
smeta %>% 
  filter(sample_age_years_before_1950 > 50000) %>% 
  select(IID)

## genetic PCs
pc_new <- read_delim('/Dedicated/jmichaelson-wdata/lcasten/sli_wgs/ancient_DNA/data/merged.pruned_snp.1000GenomesEur_PCA_projection.eigenvec', col_names = FALSE, delim = ' ')
names(pc_new) <- c('FID', 'IID', str_c('new_pc', 1:20))
pca <- pc_new %>% 
  select(IID, str_c('new_pc', 1:5))
colnames(pca) = str_remove_all(colnames(pca), 'new_')

##
kg_pc <- read_csv('/Dedicated/jmichaelson-wdata/lcasten/sli_wgs/PCA/PCA_results.merged.1000_genomes_EUR.SLI_WGS.qc.demo.csv') %>% 
  filter(population == '1000 Genomes') %>% 
  rename(IID = sample.ID)

## compute distance in 1000 genomes Eur
pc_kg <- pca %>% 
  filter(IID %in% kg_pc$IID) %>% 
  as.data.frame()
rownames(pc_kg) <- pc_kg$IID
pc_kg$IID <- NULL
pc_kg <- pc_kg[,1:5]
pc_kg$mdist <- mahalanobis(pc_kg, center = colMeans(pc_kg), cov = cov(pc_kg))

## compute distance in other samples based on 1000 genomes eur
pc_non_kg <- pca %>% 
  filter(! IID %in% kg_pc$IID) %>% 
  as.data.frame()
rownames(pc_non_kg) <- pc_kg$IID
samp_id <- pc_non_kg$IID
pc_non_kg$IID <- NULL
pc_non_kg <- pc_non_kg[,1:5]
pc_non_kg$mdist <- mahalanobis(pc_non_kg, center = colMeans(pc_kg[,1:5]), cov = cov(pc_kg[,1:5]))

## samples of interest
non_outlier_samples <- pc_non_kg %>% 
  as.data.frame() %>% 
  as_tibble() %>% 
  mutate(IID = samp_id) %>% 
  relocate(IID) %>%
  filter(mdist < median(pc_kg$mdist) + 2.5 * mad(pc_kg$mdist) & mdist > median(pc_kg$mdist) - 2.5 * mad(pc_kg$mdist)) %>%
  select(IID, mdist) %>% 
  inner_join(smeta) %>% 
  filter(sample_age_years_before_1950 > 0)

addt_samples <- pc_non_kg %>% 
  as.data.frame() %>% 
  as_tibble() %>% 
  mutate(IID = samp_id) %>% 
  relocate(IID) %>%
  # filter(str_detect(IID, 'REF|indij|enisov|ltai|eanderth|primat|hagyrsk')) %>%
  filter(str_detect(IID, 'REF|indija33|enisovaPi|ltaiNea$|primat|hagyrskaya-Ph')) %>%
  select(IID, mdist) %>% 
  left_join(smeta) # %>% 
# filter(sample_age_years_before_1950 > 0 | is.na(sample_age_years_before_1950))
addt_samples <- addt_samples %>% 
  filter(str_detect(group_id, 'Altaian', negate = TRUE) | is.na(group_id))

soi <- c(non_outlier_samples$IID, addt_samples$IID)


##############################
## basic analysis 
##############################
unrel <- read_table('/Dedicated/jmichaelson-wdata/lcasten/sli_wgs/ancient_DNA/data/unrelated_merged_ancient_european_samples.fam.rel.id', col_names= FALSE)
soi_qc <- c(unrel$X2, addt_samples$IID)


kg_samp <- kg_pc$IID
nean_samp <- addt_samples %>% 
  filter(str_detect(IID, 'REF', negate = TRUE)) %>% 
  select(IID) %>% 
  unlist() %>% 
  unname()
nean_samp = c(nean_samp, 'Chagyrskaya')

ref_samp <- addt_samples %>% 
  filter(str_detect(IID, 'REF')) %>% 
  select(IID) %>%
  filter(IID != 'Href.REF') %>% 
  unlist() %>% 
  unname()
length(soi_qc)

smeta_nean <- smeta %>% 
  filter(str_detect(IID, 'Vind|Altai|Deniso|Chagy')) %>% # select(IID)
  filter(sample_age_years_before_1950 > 0) %>% # select(IID, sample_age_years_before_1950)
  distinct(sample_age_years_before_1950, .keep_all = TRUE) %>%
  filter(str_detect(IID, 'Denisova3', negate = TRUE)) %>% # select(IID) ## wrong denisova, so drop
  mutate(IID = case_when(IID == 'AltaiNeanderthal.DG' ~ 'AltaiNea',
                         IID == 'Denisova11_noUDG.SG' ~ 'DenisovaPinky',
                         IID == 'Vindija_snpAD.DG' ~ 'Vindija33.19',
                         IID == 'Chagyrskaya_noUDG.SG' ~ 'Chagyrskaya')) %>% # select(IID, sample_age_years_before_1950)
  drop_na()

smeta2 <- bind_rows(smeta, smeta_nean)

drop_samp <- nean_samp[str_detect(nean_samp, 'DG|SG') & str_detect(nean_samp, 'Altai|Vindija|Chagyr|Denisova3')]

old_samp <- smeta2 %>% 
  filter(sample_age_years_before_1950 > 0)

##################################
## ES-PGS
##################################
pgs_corrected <- read_csv('/Dedicated/jmichaelson-wdata/lcasten/sli_wgs/ancient_DNA/data/cogPerf_ES-PGS_1000GenomesEur_corrected.csv')

list.files('/Dedicated/jmichaelson-wdata/lcasten/sli_wgs/ancient_DNA/data', pattern = 'ES-PGS_1000GenomesEur_corrected.csv')
pgs_corrected_piw <- read_csv('/Dedicated/jmichaelson-wdata/lcasten/sli_wgs/ancient_DNA/data/female_pelvic_inlet_width_ES-PGS_1000GenomesEur_corrected.csv') %>% 
  rename(pgs_pc_corrected_fpiw = pgs_pc_corrected) %>% 
  select(-PGS)%>% 
  filter(pgs_name == 'Base_1') %>% 
  select(-pgs_name)
pgs_corrected_pw <- read_csv('/Dedicated/jmichaelson-wdata/lcasten/sli_wgs/ancient_DNA/data/female_waist_circumference_bmi_adj_ES-PGS_1000GenomesEur_corrected.csv') %>% 
  rename(pgs_pc_corrected_fpw = pgs_pc_corrected) %>% 
  select(-PGS)%>% 
  filter(pgs_name == 'Base_1') %>% 
  select(-pgs_name)

pgs_corrected_pwgt <- read_csv('/Dedicated/jmichaelson-wdata/lcasten/sli_wgs/ancient_DNA/data/placental_weight_ES-PGS_1000GenomesEur_corrected.csv') %>% 
  rename(pgs_pc_corrected_pwgt = pgs_pc_corrected) %>% 
  select(-PGS)%>% 
  filter(pgs_name == 'Base_1') %>% 
  select(-pgs_name)
pgs_corrected_bw <- read_csv('/Dedicated/jmichaelson-wdata/lcasten/sli_wgs/ancient_DNA/data/fetal_birth_weight_ES-PGS_1000GenomesEur_corrected.csv') %>% 
  rename(pgs_pc_corrected_bw = pgs_pc_corrected) %>% 
  select(-PGS)%>% 
  filter(pgs_name == 'Base_1') %>% 
  select(-pgs_name)
pgs_corrected_bhc <- read_csv('/Dedicated/jmichaelson-wdata/lcasten/sli_wgs/ancient_DNA/data/birth_head_circumference_ES-PGS_1000GenomesEur_corrected.csv') %>% 
  rename(pgs_pc_corrected_bhc = pgs_pc_corrected) %>% 
  select(-PGS)%>% 
  filter(pgs_name == 'Base_1') %>% 
  select(-pgs_name)
pgs_corrected_ihc <- read_csv('/Dedicated/jmichaelson-wdata/lcasten/sli_wgs/ancient_DNA/data/infant_head_circumference_ES-PGS_1000GenomesEur_corrected.csv') %>% 
  rename(pgs_pc_corrected_ihc = pgs_pc_corrected) %>% 
  select(-PGS) %>% 
  filter(pgs_name == 'Base_1') %>% 
  select(-pgs_name)
pgs_corrected_gd <- read_csv('/Dedicated/jmichaelson-wdata/lcasten/sli_wgs/ancient_DNA/data/gestational_duration_ES-PGS_1000GenomesEur_corrected.csv') %>% 
  rename(pgs_pc_corrected_gd = pgs_pc_corrected) %>% 
  select(-PGS) %>% 
  filter(pgs_name == 'Base_1') %>% 
  select(-pgs_name)
pgs_corrected_icv <- read_csv('/Dedicated/jmichaelson-wdata/lcasten/sli_wgs/ancient_DNA/data/intracranial_volume_adult_meta_ES-PGS_1000GenomesEur_corrected.csv') %>% 
  rename(pgs_pc_corrected_icv = pgs_pc_corrected) %>% 
  select(-PGS) %>% 
  filter(pgs_name == 'Base_1') %>% 
  select(-pgs_name)
pgs_corrected_nread <- read_csv('/Dedicated/jmichaelson-wdata/lcasten/sli_wgs/ancient_DNA/data/genlang_NREAD_ES-PGS_1000GenomesEur_corrected.csv') %>% 
  rename(pgs_pc_corrected_nread = pgs_pc_corrected) %>% 
  select(-PGS) %>% 
  filter(pgs_name == 'Base_1') %>% 
  select(-pgs_name)
pgs_corrected_nrep <- read_csv('/Dedicated/jmichaelson-wdata/lcasten/sli_wgs/ancient_DNA/data/genlang_NREP_ES-PGS_1000GenomesEur_corrected.csv') %>% 
  rename(pgs_pc_corrected_nrep = pgs_pc_corrected) %>% 
  select(-PGS) %>% 
  filter(pgs_name == 'Base_1') %>% 
  select(-pgs_name)
pgs_corrected_nviq <- read_csv('/Dedicated/jmichaelson-wdata/lcasten/sli_wgs/ancient_DNA/data/genlang_NVIQ_ES-PGS_1000GenomesEur_corrected.csv') %>% 
  rename(pgs_pc_corrected_nviq = pgs_pc_corrected) %>% 
  select(-PGS) %>% 
  filter(pgs_name == 'Base_1') %>% 
  select(-pgs_name)
pgs_corrected_pa <- read_csv('/Dedicated/jmichaelson-wdata/lcasten/sli_wgs/ancient_DNA/data/genlang_PA_ES-PGS_1000GenomesEur_corrected.csv') %>% 
  rename(pgs_pc_corrected_pa = pgs_pc_corrected) %>% 
  select(-PGS) %>% 
  filter(pgs_name == 'Base_1') %>% 
  select(-pgs_name)
pgs_corrected_sp <- read_csv('/Dedicated/jmichaelson-wdata/lcasten/sli_wgs/ancient_DNA/data/genlang_SP_ES-PGS_1000GenomesEur_corrected.csv') %>% 
  rename(pgs_pc_corrected_sp = pgs_pc_corrected) %>% 
  select(-PGS) %>% 
  filter(pgs_name == 'Base_1') %>% 
  select(-pgs_name)
pgs_corrected_wr <- read_csv('/Dedicated/jmichaelson-wdata/lcasten/sli_wgs/ancient_DNA/data/genlang_WR_ES-PGS_1000GenomesEur_corrected.csv') %>% 
  rename(pgs_pc_corrected_wr = pgs_pc_corrected) %>% 
  select(-PGS) %>% 
  filter(pgs_name == 'Base_1') %>% 
  select(-pgs_name)

## loop over personality and psych PGS
new_files <- list.files('/Dedicated/jmichaelson-wdata/lcasten/sli_wgs/ancient_DNA/data', pattern = 'schiz|depres|PTSD|bipol|autism|GPC|Neur|addict|ADH|anx|brain_|reproduct|circumference_GC|vocal', full.names = TRUE)
new_files <- new_files[str_detect(new_files, 'imputed|unrelated', negate = TRUE)]
new_files <- new_files[str_detect(new_files, 'csv$')]
new_files

pgs_list <- list()
for ( f in new_files) {
  tmp <- read_csv(f, show_col_types = FALSE) %>% 
    select(-PGS) %>% 
    filter(pgs_name == 'Base_1') %>% 
    select(-pgs_name)
  ph <- str_remove_all(basename(f), '_ES-PGS_1000GenomesEur_corrected.csv|_MVP_GPC1|_MVP_UKB')
  # names(tmp)[2] <- str_c('pgs_pc_corrected_', ph)
  pgs_list[[basename(f)]] <- tmp %>% 
    mutate(trait =  str_c('pgs_pc_corrected_', ph))
}
pgs_other <- bind_rows(pgs_list) %>% 
  pivot_wider(id_cols = IID, names_from = 'trait', values_from = 'pgs_pc_corrected')
names(pgs_other)[-1]

pgs_corrected %>% 
  filter(IID %in% soi_qc) %>% 
  filter(pgs_name == 'HAQER_1') %>% 
  inner_join(pgs_corrected_piw) %>% 
  inner_join(pgs_corrected_pw) %>% 
  inner_join(pgs_corrected_pwgt) %>% 
  inner_join(pgs_corrected_bw) %>% 
  inner_join(pgs_corrected_bhc) %>% 
  inner_join(pgs_corrected_ihc) %>% 
  inner_join(pgs_corrected_gd) %>%
  inner_join(pgs_corrected_icv) %>%
  inner_join(pgs_other) %>%
  select(matches('pgs_pc_c')) %>% 
  cor()

##
wd_tmp <- pgs_corrected %>% 
  filter(str_detect(pgs_name, 'HAQER')) %>% 
  select(IID, pgs_name, pgs_pc_corrected) %>% 
  pivot_wider(id_cols = IID, names_from = pgs_name, values_from = pgs_pc_corrected) %>%
  inner_join(pgs_corrected_piw) %>% 
  inner_join(pgs_corrected_pw) %>% 
  inner_join(pgs_corrected_pwgt) %>% 
  inner_join(pgs_corrected_bw) %>% 
  inner_join(pgs_corrected_bhc) %>% 
  inner_join(pgs_corrected_ihc) %>% 
  inner_join(pgs_corrected_gd) %>% 
  inner_join(pgs_other) %>%
  filter(IID %in% soi_qc) %>%
  # filter(IID %in% kg_samp) %>%
  filter(! IID %in% drop_samp) %>%
  filter(str_detect(IID, pattern = 'Href|Gorilla|Chimp|Ancest', negate = TRUE)) %>%
  select(IID, matches('HAQER'), pgs_pc_corrected_fpiw, pgs_pc_corrected_fpw, pgs_pc_corrected_pwgt, pgs_pc_corrected_bw, pgs_pc_corrected_bhc, pgs_pc_corrected_ihc, pgs_pc_corrected_gd, pgs_pc_corrected_brain_surface_area, matches('brain_|reproductive|circumference_GC|vocal')) 

cor.test(wd_tmp$HAQER_1, wd_tmp$pgs_pc_corrected_fpiw)
cor.test(wd_tmp$HAQER_1, wd_tmp$pgs_pc_corrected_fpw)

cor.test(wd_tmp$complement_HAQER_1, wd_tmp$pgs_pc_corrected_fpiw)
cor.test(wd_tmp$complement_HAQER_1, wd_tmp$pgs_pc_corrected_fpw)

broom::tidy(lm(pgs_pc_corrected_ihc ~ complement_HAQER_1 + HAQER_1, data = wd_tmp))
broom::tidy(lm(pgs_pc_corrected_bhc ~ complement_HAQER_1 + HAQER_1, data = wd_tmp))
broom::tidy(lm(pgs_pc_corrected_fpiw ~ complement_HAQER_1 + HAQER_1, data = wd_tmp))
broom::tidy(lm(pgs_pc_corrected_fpw ~ complement_HAQER_1 + HAQER_1, data = wd_tmp))
broom::tidy(lm(pgs_pc_corrected_brain_surface_area ~ complement_HAQER_1 + HAQER_1, data = wd_tmp))
broom::tidy(lm(pgs_pc_corrected_waist_circumference_GCST90302888 ~ complement_HAQER_1 + HAQER_1, data = wd_tmp))
broom::tidy(lm(pgs_pc_corrected_reproductive_CL ~ complement_HAQER_1 + HAQER_1, data = wd_tmp))
broom::tidy(lm(pgs_pc_corrected_reproductive_NEB ~ complement_HAQER_1 + HAQER_1, data = wd_tmp))
broom::tidy(lm(pgs_pc_corrected_reproductive_NEBmen ~ complement_HAQER_1 + HAQER_1, data = wd_tmp))
broom::tidy(lm(pgs_pc_corrected_reproductive_NEBwomen ~ complement_HAQER_1 + HAQER_1, data = wd_tmp))
broom::tidy(lm(pgs_pc_corrected_vocal_pitch ~ complement_HAQER_1 + HAQER_1, data = wd_tmp))

pgs_corrected %>% 
  filter(IID %in% soi_qc) %>% # distinct(IID)
  filter(str_detect(pgs_name, 'HAQER')) %>% 
  select(IID, pgs_name, pgs_pc_corrected) %>% 
  pivot_wider(id_cols = IID, names_from = pgs_name, values_from = pgs_pc_corrected) %>%
  inner_join(pgs_corrected_piw) %>% 
  inner_join(pgs_corrected_pw) %>% 
  inner_join(pgs_corrected_pwgt) %>% 
  inner_join(pgs_corrected_bw) %>% 
  inner_join(pgs_corrected_bhc) %>% 
  inner_join(pgs_corrected_ihc) %>% 
  inner_join(pgs_corrected_gd) %>% 
  inner_join(pgs_corrected_icv) %>%
  inner_join(pgs_corrected_nviq) %>%
  inner_join(pgs_corrected_nread) %>%
  inner_join(pgs_corrected_nrep) %>%
  inner_join(pgs_corrected_sp) %>%
  inner_join(pgs_corrected_pa) %>%
  inner_join(pgs_corrected_wr) %>%
  inner_join(pgs_other) %>%
  filter(! IID %in% drop_samp) %>%
  filter(str_detect(IID, pattern = 'Href|Gorilla|Chimp|Ancest', negate = TRUE)) %>%
  filter(IID %in% old_samp$IID | IID == 'Chagyrskaya-Phalanx') %>% # select(IID)
  select(IID, matches('HAQER'), pgs_pc_corrected_fpiw, pgs_pc_corrected_fpw, pgs_pc_corrected_pwgt, pgs_pc_corrected_bw, pgs_pc_corrected_bhc, pgs_pc_corrected_ihc, pgs_pc_corrected_gd, pgs_pc_corrected_icv, pgs_pc_corrected_nviq, pgs_pc_corrected_nread, pgs_pc_corrected_nrep, pgs_pc_corrected_sp, pgs_pc_corrected_pa, pgs_pc_corrected_wr, matches('brain_|reproduc|waist')) %>% 
  pivot_longer(cols = matches('pgs_pc_corr')) %>%
  group_by(name) %>%
  do(res = broom::tidy(lm(value ~ HAQER_1 + complement_HAQER_1, data = .))) %>% 
  unnest(res) %>% 
  filter(term != '(Intercept)') %>%
  arrange(p.value) %>% 
  filter(term == 'HAQER_1')

## genlang reading polygenic scores
p_pgs_gl <- pgs_corrected %>% 
  filter(str_detect(pgs_name, 'HAQER')) %>% 
  select(IID, pgs_name, pgs_pc_corrected) %>% 
  pivot_wider(id_cols = IID, names_from = pgs_name, values_from = pgs_pc_corrected) %>%
  inner_join(pgs_corrected_piw) %>% 
  inner_join(pgs_corrected_pw) %>% 
  inner_join(pgs_corrected_pwgt) %>% 
  inner_join(pgs_corrected_bw) %>% 
  inner_join(pgs_corrected_bhc) %>% 
  inner_join(pgs_corrected_ihc) %>% 
  inner_join(pgs_corrected_gd) %>% 
  inner_join(pgs_corrected_icv) %>%
  inner_join(pgs_corrected_nviq) %>%
  inner_join(pgs_corrected_nread) %>%
  inner_join(pgs_corrected_nrep) %>%
  inner_join(pgs_corrected_sp) %>%
  inner_join(pgs_corrected_pa) %>%
  inner_join(pgs_corrected_wr) %>%
  filter(IID %in% soi_qc) %>% # distinct(IID)
  filter(! IID %in% drop_samp) %>%
  filter(str_detect(IID, pattern = 'Href|Gorilla|Chimp|Ancest', negate = TRUE)) %>%
  filter(IID %in% old_samp$IID | IID == 'Chagyrskaya-Phalanx') %>% # select(IID)
  select(IID, matches('HAQER'), pgs_pc_corrected_fpiw, pgs_pc_corrected_fpw, pgs_pc_corrected_pwgt, pgs_pc_corrected_bw, pgs_pc_corrected_bhc, pgs_pc_corrected_ihc, pgs_pc_corrected_gd, pgs_pc_corrected_icv, pgs_pc_corrected_nviq, pgs_pc_corrected_nread, pgs_pc_corrected_nrep, pgs_pc_corrected_sp, pgs_pc_corrected_pa, pgs_pc_corrected_wr) %>% 
  pivot_longer(cols = matches('pgs_pc_corr')) %>%
  group_by(name) %>%
  do(res = broom::tidy(lm(value ~ HAQER_1 + complement_HAQER_1, data = .))) %>% 
  unnest(res) %>% 
  filter(term != '(Intercept)') %>%
  arrange(p.value) %>% 
  mutate(clean_name = case_when(str_detect(name, '_bw$') ~ 'Birth weight PGS',
                                str_detect(name, '_ihc$') ~ 'Infant head circumference PGS',
                                str_detect(name, '_wr$') ~ 'Word reading PGS',
                                str_detect(name, '_nread$') ~ 'Non-word reading PGS',
                                str_detect(name, '_icv$') ~ 'Adult intracranial volume PGS',
                                str_detect(name, '_nrep$') ~ 'Non-word repetition PGS',
                                str_detect(name, '_pa$') ~ 'Phoneme awareness PGS',
                                str_detect(name, '_nviq$') ~ 'Nonverbal IQ PGS',
                                str_detect(name, '_sp$') ~ 'Spelling PGS',
                                str_detect(name, '_gd$') ~ 'Gestational duration PGS',
                                str_detect(name, '_bhc$') ~ 'Birth head circumference PGS',
                                TRUE ~ NA_character_)) %>% 
  drop_na() %>%
  mutate(term = ifelse(term == 'HAQER_1', 'HAQER CP-PGS', 'Background CP-PGS')) %>% 
  mutate(term = factor(term, levels = rev(c('HAQER CP-PGS', 'Background CP-PGS')))) %>% 
  arrange(desc(term), estimate) %>% 
  mutate(clean_name = factor(clean_name, levels = unique(clean_name))) %>% 
  arrange(term) %>%
  filter(str_detect(clean_name, 'Non|Phen|verbal|Spelling|read|honeme')) %>%
  ggplot(aes(x = estimate, y = clean_name, color = term)) +
  geom_point(size  = 6, position = position_dodge2(width = 0.3)) +
  geom_linerange(aes(xmin = estimate - 1.96 * std.error, xmax = estimate + 1.96 * std.error), position = position_dodge2(width = 0.3), size = 1.2) +
  geom_vline(xintercept = 0, color = 'red', linetype = 'dashed', size = 1.1) +
  xlab('Regression beta (95% CI)') +
  ylab(NULL) +
  labs(color = NULL)

## pgs brain
p_pgs_br <- pgs_corrected %>% 
  filter(str_detect(pgs_name, 'HAQER')) %>% 
  select(IID, pgs_name, pgs_pc_corrected) %>% 
  pivot_wider(id_cols = IID, names_from = pgs_name, values_from = pgs_pc_corrected) %>%
  inner_join(pgs_corrected_piw) %>% 
  inner_join(pgs_corrected_pw) %>% 
  inner_join(pgs_corrected_pwgt) %>% 
  inner_join(pgs_corrected_bw) %>% 
  inner_join(pgs_corrected_bhc) %>% 
  inner_join(pgs_corrected_ihc) %>% 
  inner_join(pgs_corrected_gd) %>% 
  inner_join(pgs_corrected_icv) %>%
  inner_join(pgs_corrected_nviq) %>%
  inner_join(pgs_corrected_nread) %>%
  inner_join(pgs_corrected_nrep) %>%
  inner_join(pgs_corrected_sp) %>%
  inner_join(pgs_corrected_pa) %>%
  inner_join(pgs_corrected_wr) %>%
  inner_join(pgs_other) %>%
  filter(IID %in% soi_qc) %>% # distinct(IID)
  filter(! IID %in% drop_samp) %>%
  filter(str_detect(IID, pattern = 'Href|Gorilla|Chimp|Ancest', negate = TRUE)) %>%
  filter(IID %in% old_samp$IID | IID == 'Chagyrskaya-Phalanx') %>% # select(IID)
  select(IID, matches('HAQER'), pgs_pc_corrected_pwgt, pgs_pc_corrected_bw, pgs_pc_corrected_bhc, pgs_pc_corrected_icv, matches('brain_')) %>% 
  pivot_longer(cols = matches('pgs_pc_corr')) %>%
  filter(str_detect(name, 'caud|amyg|thalam|brainstem|surface|pallid|icv|putam|Folding')) %>%
  group_by(name) %>%
  do(res = broom::tidy(lm(value ~ HAQER_1 + complement_HAQER_1, data = .))) %>% 
  unnest(res) %>% 
  filter(term != '(Intercept)') %>%
  arrange(p.value) %>% 
  mutate(clean_name = str_remove_all(name, pattern = 'brain_|pgs_pc_corrected_')) %>% 
  mutate(clean_name = ifelse(clean_name == 'icv', 'intracranial_vol', clean_name)) %>%
  drop_na() %>%
  mutate(term = ifelse(term == 'HAQER_1', 'HAQER CP-PGS', 'Background CP-PGS')) %>% 
  mutate(term = factor(term, levels = rev(c('HAQER CP-PGS', 'Background CP-PGS')))) %>% 
  arrange(desc(term), estimate) %>% 
  mutate(clean_name = factor(clean_name, levels = unique(clean_name))) %>% 
  arrange(term) %>%
  # filter(str_detect(clean_name, 'Non|Phen|verbal|Spelling|read|honeme')) %>%
  ggplot(aes(x = estimate, y = clean_name, color = term)) +
  geom_point(size  = 6, position = position_dodge2(width = 0.3)) +
  geom_linerange(aes(xmin = estimate - 1.96 * std.error, xmax = estimate + 1.96 * std.error), position = position_dodge2(width = 0.3), size = 1.2) +
  geom_vline(xintercept = 0, color = 'red', linetype = 'dashed', size = 1.1) +
  xlab('Regression beta (95% CI)') +
  ylab(NULL) +
  labs(color = NULL)

## PGC neuropsych and personality traits
##
p_pgs_psych <- pgs_corrected %>% 
  filter(str_detect(pgs_name, 'HAQER')) %>% 
  select(IID, pgs_name, pgs_pc_corrected) %>% 
  pivot_wider(id_cols = IID, names_from = pgs_name, values_from = pgs_pc_corrected) %>%
  inner_join(pgs_other) %>%
  filter(IID %in% soi_qc) %>%
  filter(! IID %in% drop_samp) %>%
  filter(str_detect(IID, pattern = 'Href|Gorilla|Chimp|Ancest', negate = TRUE)) %>%
  filter(IID %in% old_samp$IID | IID == 'Chagyrskaya-Phalanx') %>% # select(IID)
  select(IID, matches('HAQER'), matches('pgs_pc_corrected_')) %>% 
  pivot_longer(cols = matches('pgs_pc_corr')) %>%
  group_by(name) %>%
  do(res = broom::tidy(lm(value ~ HAQER_1 + complement_HAQER_1, data = .))) %>% 
  unnest(res) %>% 
  filter(term != '(Intercept)') %>%
  arrange(p.value) %>% 
  mutate(clean_name = case_when(str_detect(name, '_ADHD$') ~ 'ADHD PGS',
                                str_detect(name, '_addiction$') ~ 'Addiction PGS',
                                str_detect(name, '_anxiety$') ~ 'Anxiety PGS',
                                str_detect(name, '_PTSD$') ~ 'PTSD PGS',
                                str_detect(name, '_autism$') ~ 'Autism PGS',
                                str_detect(name, '_open$') ~ 'Openness PGS',
                                str_detect(name, '_schizophrenia$') ~ 'Schizophrenia PGS',
                                str_detect(name, '_bipolar$') ~ 'Bipolar PGS',
                                str_detect(name, '_consc$') ~ 'Conscientiousness PGS',
                                str_detect(name, '_Neuro$') ~ 'Neuroticism PGS',
                                str_detect(name, '_depression$') ~ 'Depression PGS',
                                str_detect(name, '_agree$') ~ 'Agreeableness PGS',
                                str_detect(name, '_extra$') ~ 'Extraversion PGS',
                                TRUE ~ NA_character_)) %>% 
  drop_na() %>%
  mutate(term = ifelse(term == 'HAQER_1', 'HAQER CP-PGS', 'Background CP-PGS')) %>% 
  mutate(term = factor(term, levels = rev(c('HAQER CP-PGS', 'Background CP-PGS')))) %>% 
  arrange(desc(term), estimate) %>% 
  mutate(clean_name = factor(clean_name, levels = unique(clean_name))) %>% 
  arrange(term) %>%
  ggplot(aes(x = estimate, y = clean_name, color = term)) +
  geom_point(size  = 6, position = position_dodge2(width = 0.3)) +
  geom_linerange(aes(xmin = estimate - 1.96 * std.error, xmax = estimate + 1.96 * std.error), position = position_dodge2(width = 0.3), size = 1.2) +
  geom_vline(xintercept = 0, color = 'red', linetype = 'dashed', size = 1.1) +
  xlab('Regression beta (95% CI)') +
  ylab(NULL) +
  labs(color = NULL)

##
p_pgs <- pgs_corrected %>% 
  filter(IID %in% soi_qc) %>% # distinct(IID)
  filter(str_detect(pgs_name, 'HAQER')) %>% 
  select(IID, pgs_name, pgs_pc_corrected) %>% 
  pivot_wider(id_cols = IID, names_from = pgs_name, values_from = pgs_pc_corrected) %>%
  inner_join(pgs_corrected_piw) %>% 
  inner_join(pgs_corrected_pw) %>% 
  inner_join(pgs_corrected_pwgt) %>% 
  inner_join(pgs_corrected_bw) %>% 
  inner_join(pgs_corrected_bhc) %>% 
  inner_join(pgs_corrected_ihc) %>% 
  inner_join(pgs_corrected_gd) %>% 
  inner_join(pgs_corrected_icv) %>%
  inner_join(pgs_corrected_nviq) %>%
  inner_join(pgs_corrected_nread) %>%
  inner_join(pgs_corrected_nrep) %>%
  inner_join(pgs_corrected_sp) %>%
  inner_join(pgs_corrected_pa) %>%
  inner_join(pgs_corrected_wr) %>%
  inner_join(select(pgs_other, IID, pgs_pc_corrected_brain_surface_area)) %>%
  filter(! IID %in% drop_samp) %>%
  filter(str_detect(IID, pattern = 'Href|Gorilla|Chimp|Ancest', negate = TRUE)) %>%
  filter(IID %in% old_samp$IID | IID == 'Chagyrskaya-Phalanx') %>% # select(IID)
  select(IID, matches('HAQER'), pgs_pc_corrected_fpiw, pgs_pc_corrected_fpw, pgs_pc_corrected_pwgt, pgs_pc_corrected_bw, pgs_pc_corrected_bhc, pgs_pc_corrected_ihc, pgs_pc_corrected_gd, pgs_pc_corrected_icv, pgs_pc_corrected_nviq, pgs_pc_corrected_nread, pgs_pc_corrected_nrep, pgs_pc_corrected_sp, pgs_pc_corrected_pa, pgs_pc_corrected_wr, pgs_pc_corrected_brain_surface_area, pgs_pc_corrected_brain_surface_area) %>% 
  pivot_longer(cols = matches('pgs_pc_corr')) %>%
  group_by(name) %>%
  do(res = broom::tidy(lm(value ~ HAQER_1 + complement_HAQER_1, data = .))) %>% 
  unnest(res) %>% 
  filter(term != '(Intercept)') %>%
  arrange(p.value) %>% 
  mutate(clean_name = case_when(str_detect(name, '_bw$') ~ 'Birth weight PGS',
                                # str_detect(name, '_ihc$') ~ 'Infant head circumference PGS',
                                str_detect(name, '_wr$') ~ 'Word reading PGS',
                                str_detect(name, '_nread$') ~ 'Non-word reading PGS',
                                str_detect(name, '_icv$') ~ 'Adult intracranial volume PGS',
                                str_detect(name, '_nrep$') ~ 'Non-word repetition PGS',
                                str_detect(name, '_pa$') ~ 'Phoneme awareness PGS',
                                str_detect(name, '_nviq$') ~ 'Nonverbal IQ PGS',
                                str_detect(name, '_sp$') ~ 'Spelling PGS',
                                str_detect(name, '_gd$') ~ 'Gestational duration PGS',
                                str_detect(name, '_bhc$') ~ 'Birth head circumference PGS',
                                str_detect(name, '_surface_area') ~ 'Brain surface area PGS',
                                TRUE ~ NA_character_)) %>% 
  drop_na() %>%
  mutate(term = ifelse(term == 'HAQER_1', 'HAQER CP-PGS', 'Background CP-PGS')) %>% 
  mutate(term = factor(term, levels = rev(c('HAQER CP-PGS', 'Background CP-PGS')))) %>% 
  arrange(desc(term), estimate) %>% 
  mutate(clean_name = factor(clean_name, levels = unique(clean_name))) %>% 
  arrange(term) %>%
  filter(str_detect(clean_name, 'Non|Phen|verbal|Spelling|read|honeme', negate = TRUE)) %>%
  ggplot(aes(x = estimate, y = clean_name, color = term)) +
  geom_point(size  = 6, position = position_dodge2(width = 0.3)) +
  geom_linerange(aes(xmin = estimate - 1.96 * std.error, xmax = estimate + 1.96 * std.error), position = position_dodge2(width = 0.3), size = 1.2) +
  geom_vline(xintercept = 0, color = 'red', linetype = 'dashed', size = 1.1) +
  xlab('Regression beta (95% CI)') +
  ylab(NULL) +
  labs(color = NULL)
p_pgs + p_pgs_psych + p_pgs_gl

p_pgs %>% 
  ggsave(filename = '/wdata/lcasten/sli_wgs/paper_figures/subplots/ancient_dna_ES-PGS_associations_with_birth_pgs.png', device = 'png', units = 'in', dpi = 300, width = 12, height = 6)
length(unique(pgs_corrected$IID))


##############################
## basic analysis 
##############################
# unrel <- read_table('/Dedicated/jmichaelson-wdata/lcasten/sli_wgs/ancient_DNA/data/unrelated_merged_ancient_european_samples.fam.rel.id', col_names= FALSE)
# addt_samples <- addt_samples %>% 
#   filter(str_detect(group_id, 'Altaian', negate = TRUE) | is.na(group_id))

# kg_samp <- kg_pc$IID
# soi_qc <- c(unrel$X2, addt_samples$IID)
# nean_samp <- addt_samples %>% 
#   filter(str_detect(IID, 'REF', negate = TRUE)) %>% 
#   select(IID) %>% 
#   unlist() %>% 
#   unname()
# ref_samp <- addt_samples %>% 
#   filter(str_detect(IID, 'REF')) %>% 
#   select(IID) %>%
#   filter(IID != 'Href.REF') %>% 
#   unlist() %>% 
#   unname()
# length(soi_qc)

# smeta_nean <- smeta %>% 
#   filter(str_detect(IID, 'Vind|Altai|Deniso|Chagy')) %>% # select(IID)
#   filter(sample_age_years_before_1950 > 0) %>% 
#   distinct(sample_age_years_before_1950, .keep_all = TRUE) %>%
#   filter(str_detect(IID, 'Denisova3', negate = TRUE)) %>% select(IID) ## wrong denisova, so drop
#   mutate(IID = case_when(IID == 'AltaiNeanderthal.DG' ~ 'AltaiNea',
#                          IID == 'Denisova11_noUDG.SG' ~ 'DenisovaPinky',
#                          IID == 'Vindija_snpAD.DG' ~ 'Vindija33.19',
#                          IID == 'Chagyrskaya_noUDG.SG' ~ 'Chagyrskaya'))

# smeta2 <- bind_rows(smeta, smeta_nean)

pgs_wd <- pgs_corrected %>% 
  filter(IID %in% kg_samp | IID %in% soi_qc) %>% 
  mutate(type = case_when(IID %in% kg_samp ~ '1000 Genomes Europeans',
                          IID %in% nean_samp ~ 'Neanderthals',
                          IID %in% unrel$X2 ~ 'Ancient Europeans',
                          IID %in% ref_samp ~ 'Primate reference genomes')) %>% 
  relocate(type, .after = IID)

pgs_wd %>% 
  filter(type == 'Neanderthals')

pgs_wd %>% 
  filter(pgs_name == 'Base_1') %>%
  ggplot(aes(x = pgs_pc_corrected, fill = type)) +
  geom_density(alpha = 0.7)

pgs_wd %>% 
  filter(pgs_name == 'HAQER_1') %>%
  ggplot(aes(x = pgs_pc_corrected, fill = type)) +
  geom_density()

##
pgs_wd2 <- pgs_corrected %>% 
  inner_join(pgs_corrected_piw) %>% 
  inner_join(pgs_corrected_pw) %>%
  inner_join(pgs_corrected_pwgt) %>%
  inner_join(pgs_corrected_bw) %>%
  inner_join(pgs_corrected_ihc) %>%
  inner_join(pgs_corrected_bhc) %>%
  inner_join(pgs_corrected_gd) %>% 
  inner_join(pgs_corrected_icv) %>%
  inner_join(pgs_corrected_nviq) %>%
  inner_join(pgs_corrected_nread) %>%
  inner_join(pgs_corrected_nrep) %>%
  inner_join(pgs_corrected_sp) %>%
  inner_join(pgs_corrected_pa) %>%
  inner_join(pgs_corrected_wr) %>%
  inner_join(pgs_other) %>%
  filter(IID %in% kg_samp | IID %in% soi_qc) %>% 
  mutate(type = case_when(IID %in% kg_samp ~ '1000 Genomes Europeans',
                          IID %in% nean_samp ~ 'Neanderthals',
                          IID %in% unrel$X2 ~ 'Ancient Europeans',
                          IID %in% ref_samp ~ 'Primate reference genomes')) %>% 
  mutate(IID = case_when(IID == 'Chagyrskaya-Phalanx' ~ 'Chagyrskaya',
                         TRUE ~ IID)) %>%
  inner_join(smeta2)

drop_samp <- nean_samp[str_detect(nean_samp, 'DG|SG') & str_detect(nean_samp, 'Altai|Vindija|Chagyr|Denisova3')]

#################################################
## density plot of neanderthals vs modern humans
#################################################
tmp_dat <- pgs_corrected %>% 
  inner_join(pgs_corrected_piw) %>% 
  inner_join(pgs_corrected_pw) %>%
  inner_join(pgs_corrected_pwgt) %>%
  inner_join(pgs_corrected_bw) %>%
  inner_join(pgs_corrected_ihc) %>%
  inner_join(pgs_corrected_bhc) %>%
  inner_join(pgs_corrected_gd) %>% 
  inner_join(pgs_corrected_icv) %>%
  inner_join(pgs_corrected_nviq) %>%
  inner_join(pgs_corrected_nread) %>%
  inner_join(pgs_corrected_nrep) %>%
  inner_join(pgs_corrected_sp) %>%
  inner_join(pgs_corrected_pa) %>%
  inner_join(pgs_corrected_wr) %>%
  inner_join(pgs_other) %>%
  filter(IID %in% kg_samp | IID %in% soi_qc) %>% 
  filter(! IID %in% drop_samp) %>%
  filter(str_detect(IID, pattern = 'Href|Gorilla|Chimp|Ancest', negate = TRUE)) %>%
  # filter(IID %in% old_samp$IID | IID == 'Chagyrskaya-Phalanx') %>% # select(IID)
  mutate(type = case_when(IID %in% kg_samp ~ '1000 Genomes Europeans',
                          IID %in% nean_samp ~ 'Neanderthals',
                          IID %in% unrel$X2 ~ 'Ancient Europeans',
                          IID %in% ref_samp ~ 'Primate reference genomes'))

tmp_nean <- tmp_dat %>% 
  # filter(sample_age_years_before_1950 > 100) %>%
  filter(! IID %in% drop_samp) %>%
  filter(IID != 'Denisova11_noUDG.SG') %>%
  filter(type %in% c('Neanderthals')) %>%
  filter(pgs_name %in% c('HAQER_1', 'complement_HAQER_1'))
tmp_nean %>% 
  select(IID, matches('pgs_pc_corrected_')) %>% 
  distinct() %>% 
  pivot_longer(cols = -1) %>% 
  group_by(name) %>% 
  summarise(mean(value),
            min(value),
            max(value)) %>% 
  arrange(desc(`mean(value)`))

tmp_kg <- tmp_dat %>% 
  filter(! IID %in% drop_samp) %>%
  filter(type %in% c('1000 Genomes Europeans')) %>%
  filter(pgs_name %in% c('HAQER_1', 'complement_HAQER_1'))

p_haq_neand <- tmp_kg %>% 
  filter(pgs_name == 'HAQER_1') %>% 
  ggplot(aes(x = pgs_pc_corrected)) +
  geom_density(alpha = 0.85, fill = '#00bfc4') +
  geom_vline(data = tmp_nean[tmp_nean$pgs_name == 'HAQER_1',], aes(xintercept = pgs_pc_corrected), color = 'black', size = 1.25) +
  xlab('HAQER CP-PGS') +
  ylab('Density') +
  geom_text(aes(x = .75, y = 0.425, label = 'Neanderthals'), check_overlap = TRUE, size  = 6)

p_bg_neand <- tmp_kg %>% 
  filter(pgs_name == 'complement_HAQER_1') %>% 
  ggplot(aes(x = pgs_pc_corrected)) +
  geom_density(alpha = 0.65, aes(fill = 'na')) +
  geom_vline(data = tmp_nean[tmp_nean$pgs_name == 'complement_HAQER_1',], aes(xintercept = pgs_pc_corrected), color = 'black', size = 1.25) +
  xlab('Background CP-PGS') +
  ylab('Density') +
  geom_text(aes(x = -1, y = 0.425, label = 'Neanderthals'), check_overlap = TRUE, size  = 6) +
  theme(legend.position = 'none')

ord <- tmp_nean %>% 
  select(IID, matches('pgs_pc_corrected_')) %>% 
  distinct() %>% 
  pivot_longer(cols = -1) %>% 
  group_by(name) %>% 
  summarise(mean(value),
            min(value),
            max(value)) %>% 
  arrange(desc(`mean(value)`))
ord %>% 
  filter(str_detect(name, 'waist|fpiw|fp'))
tmp_dat %>% 
  filter(! IID %in% drop_samp) %>%
  # filter(IID != 'Denisova11_noUDG.SG') %>%
  filter(str_detect(type, 'Neand|1000')) %>% 
  select(IID, type, matches('pgs_pc_corrected_')) %>% 
  distinct() %>% 
  pivot_longer(cols = -c(1:2)) %>% 
  mutate(clean_name = case_when(str_detect(name, '_bw$') ~ 'Birth weight PGS',
                                # str_detect(name, '_ihc$') ~ 'Infant head circumference PGS',
                                str_detect(name, '_wr$') ~ 'Word reading PGS',
                                str_detect(name, '_nread$') ~ 'Non-word reading PGS',
                                str_detect(name, '_icv$') ~ 'Adult intracranial volume PGS',
                                str_detect(name, '_nrep$') ~ 'Non-word repetition PGS',
                                str_detect(name, '_pa$') ~ 'Phoneme awareness PGS',
                                str_detect(name, '_nviq$') ~ 'Nonverbal IQ PGS',
                                str_detect(name, '_sp$') ~ 'Spelling PGS',
                                str_detect(name, '_gd$') ~ 'Gestational duration PGS',
                                str_detect(name, '_bhc$') ~ 'Birth head circumference PGS',
                                TRUE ~ NA_character_)) %>% 
  mutate(name = factor(name, levels = ord$name)) %>% 
  arrange(name) %>% 
  mutate(clean_name = factor(clean_name, levels = unique(clean_name))) %>%
  drop_na() %>%
  ggplot(aes(x = value, fill = type)) +
  geom_density(alpha = 0.7) +
  facet_wrap(~ clean_name, scales = 'free') +
  xlab('Polygenic score')

tmp_dat %>% 
  filter(! IID %in% drop_samp) %>%
  # filter(IID != 'Denisova11_noUDG.SG') %>%
  filter(str_detect(type, 'Neand|1000')) %>% 
  select(IID, type, matches('pgs_pc_corrected_')) %>% 
  distinct() %>% 
  pivot_longer(cols = -c(1:2)) %>% 
  mutate(clean_name = case_when(str_detect(name, '_ADHD$') ~ 'ADHD PGS',
                                str_detect(name, '_addiction$') ~ 'Addiction PGS',
                                str_detect(name, '_anxiety$') ~ 'Anxiety PGS',
                                str_detect(name, '_PTSD$') ~ 'PTSD PGS',
                                str_detect(name, '_autism$') ~ 'Autism PGS',
                                str_detect(name, '_open$') ~ 'Openness PGS',
                                str_detect(name, '_schizophrenia$') ~ 'Schizophrenia PGS',
                                str_detect(name, '_bipolar$') ~ 'Bipolar PGS',
                                str_detect(name, '_consc$') ~ 'Conscientiousness PGS',
                                str_detect(name, '_Neuro$') ~ 'Neuroticism PGS',
                                str_detect(name, '_depression$') ~ 'Depression PGS',
                                str_detect(name, '_agree$') ~ 'Agreeableness PGS',
                                str_detect(name, '_extra$') ~ 'Extraversion PGS',
                                TRUE ~ NA_character_)) %>% 
  mutate(name = factor(name, levels = ord$name)) %>% 
  arrange(name) %>% 
  mutate(clean_name = factor(clean_name, levels = unique(clean_name))) %>%
  drop_na() %>%
  ggplot(aes(x = value, fill = type)) +
  geom_density(alpha = 0.7) +
  facet_wrap(~ clean_name, scales = 'free') +
  xlab('Polygenic score')


tmp_dat %>% 
  filter(! IID %in% drop_samp) %>%
  # filter(IID != 'Denisova11_noUDG.SG') %>%
  filter(str_detect(type, 'Neand|1000')) %>% 
  select(IID, type, matches('pgs_pc_corrected_')) %>% 
  distinct() %>% 
  select(IID, type, matches('brain_|icv')) %>% 
  pivot_longer(cols = -c(1:2)) %>% 
  mutate(clean_name = str_remove_all(name, pattern = 'pgs_pc_corrected_brain_|pgs_pc_corrected_')) %>% 
  mutate(name = factor(name, levels = ord$name)) %>% 
  arrange(name) %>% 
  filter(str_detect(name, 'caud|amyg|thalam|brainstem|surface|pallid|icv|putam')) %>%
  mutate(clean_name = factor(clean_name, levels = unique(clean_name))) %>%
  drop_na() %>%
  ggplot(aes(x = value, fill = type)) +
  geom_density(alpha = 0.7) +
  facet_wrap(~ clean_name, scales = 'free', nrow = 2) +
  xlab('Polygenic score')

tmp_dat %>% 
  filter(! IID %in% drop_samp) %>%
  # filter(IID != 'Denisova11_noUDG.SG') %>%
  filter(str_detect(type, 'Neand|1000')) %>% 
  select(IID, type, matches('pgs_pc_corrected_')) %>% 
  distinct() %>% 
  select(IID, type, matches('brain_|icv')) %>% 
  pivot_longer(cols = -c(1:2)) %>% 
  mutate(clean_name = str_remove_all(name, pattern = 'pgs_pc_corrected_brain_|pgs_pc_corrected_')) %>% 
  mutate(name = factor(name, levels = ord$name)) %>% 
  arrange(name) %>% 
  filter(str_detect(name, 'caud|amyg|thalam|brainstem|pallid|putam|Volume', negate = TRUE)) %>%
  mutate(clean_name = factor(clean_name, levels = unique(clean_name))) %>%
  drop_na() %>%
  ggplot(aes(x = value, fill = type)) +
  geom_density(alpha = 0.7) +
  facet_wrap(~ clean_name, scales = 'free', nrow = 3) +
  xlab('Polygenic score')


#########################
## polygenic selection
#########################
unique(pgs_corrected$pgs_name)

pgs_res <- pgs_corrected %>% 
  filter(IID %in% kg_samp | IID %in% soi_qc) %>% 
  mutate(IID = case_when(IID == 'Chagyrskaya-Phalanx' ~ 'Chagyrskaya',
                         TRUE ~ IID)) %>%
  mutate(type = case_when(IID %in% kg_samp ~ '1000 Genomes Europeans',
                          IID %in% nean_samp ~ 'Neanderthals',
                          IID %in% unrel$X2 ~ 'Ancient Europeans',
                          IID %in% ref_samp ~ 'Primate reference genomes')) %>% 
  inner_join(smeta2) %>% 
  select(-PGS) %>%
  pivot_wider(id_cols = -matches('pgs_|PGS'), names_from = pgs_name, values_from = pgs_pc_corrected) %>%#  names()
  mutate(pgs_name = 'HAQER') %>%
  group_by(pgs_name) %>% 
  filter(sample_age_years_before_1950 > 100) %>%
  filter(! IID %in% drop_samp) %>%
  filter(str_detect(IID, pattern = 'Href|REF|Denisova11_noUDG.SG', negate = TRUE)) %>%
  # filter(type == 'Ancient Europeans') %>%
  # do(res = broom::tidy(cor.test(.$pgs_pc_corrected, log10((.$sample_age_years_before_1950 + 74) / 1000)), method = 'p')) %>% 
  # do(res = broom::tidy(cor.test(.$pgs_pc_corrected, .$sample_age_years_before_1950, method = 's'))) %>% 
  do(res = broom::tidy(glm(scale(log10((.$sample_age_years_before_1950 + 74) / 1000))[,1] ~ HAQER_1 + complement_HAQER_1, data = .))) %>%
  unnest(res) %>% 
  arrange(p.value) %>% 
  filter(str_detect(pgs_name, 'Base|HAQER|^HAR')) %>% 
  mutate(lab = str_c('r = ', round(-1 * estimate, digits = 2), ', p-val = ', formatC(p.value, digits = 2)))


## compute polygenic selective pressures across birth traits
pgs_res_birth <- pgs_wd2 %>% 
  filter(pgs_name == 'Base_1') %>%
  select(-c(pgs_name, PGS, pgs_pc_corrected)) %>% 
  distinct() %>% 
  pivot_longer(cols = matches('pgs_pc_corrected'), names_to = 'pgs_name', values_to = 'pgs_pc_corrected') %>%
  filter(sample_age_years_before_1950 > 100) %>%
  filter(! IID %in% drop_samp) %>%
  filter(str_detect(IID, pattern = 'Href|Gorilla|Chimp|Ancest|Denisova11_noUDG.SG', negate = TRUE)) %>%
  group_by(pgs_name) %>% 
  # filter(type == 'Ancient Europeans') %>%
  # filter(sample_age_years_before_1950 < 50000) %>%
  do(res = broom::tidy(cor.test(.$pgs_pc_corrected, .$sample_age_years_before_1950), method = 's')) %>% 
  # do(res = broom::tidy(cor.test(.$pgs_pc_corrected, log10((.$sample_age_years_before_1950 + 74) / 1000)), method = 'p')) %>% 
  unnest(res) %>% 
  arrange(p.value) %>%
  rename(name = pgs_name) %>% 
  mutate(clean_name = case_when(str_detect(name, '_bw$') ~ 'Birth weight PGS',
                                str_detect(name, '_ihc$') ~ 'Infant head circumference PGS',
                                str_detect(name, '_wr$') ~ 'Word reading PGS',
                                str_detect(name, '_nread$') ~ 'Non-word reading PGS',
                                str_detect(name, '_icv$') ~ 'Adult intracranial volume PGS',
                                str_detect(name, '_nrep$') ~ 'Non-word repetition PGS',
                                str_detect(name, '_pa$') ~ 'Phoneme awareness PGS',
                                str_detect(name, '_nviq$') ~ 'Nonverbal IQ PGS',
                                str_detect(name, '_sp$') ~ 'Spelling PGS',
                                str_detect(name, '_gd$') ~ 'Gestational duration PGS',
                                str_detect(name, '_bhc$') ~ 'Birth head circumference PGS',
                                str_detect(name, '_ADHD$') ~ 'ADHD PGS',
                                str_detect(name, '_addiction$') ~ 'Addiction PGS',
                                str_detect(name, '_anxiety$') ~ 'Anxiety PGS',
                                str_detect(name, '_PTSD$') ~ 'PTSD PGS',
                                str_detect(name, '_autism$') ~ 'Autism PGS',
                                str_detect(name, '_open$') ~ 'Openness PGS',
                                str_detect(name, '_schizophrenia$') ~ 'Schizophrenia PGS',
                                str_detect(name, '_bipolar$') ~ 'Bipolar PGS',
                                str_detect(name, '_consc$') ~ 'Conscientiousness PGS',
                                str_detect(name, '_Neuro$') ~ 'Neuroticism PGS',
                                str_detect(name, '_depression$') ~ 'Depression PGS',
                                str_detect(name, '_agree$') ~ 'Agreeableness PGS',
                                str_detect(name, '_extra$') ~ 'Extraversion PGS',
                                TRUE ~ NA_character_)) %>% 
  mutate(lab = str_c('r = ', round(-1 * estimate, digits = 2), ', p-val = ', formatC(p.value, digits = 2))) %>% ## flip direction to match time
  relocate(clean_name) %>% 
  rename(pgs_name = name) %>% 
  select(-pgs_name) %>% 
  drop_na()

pgs_res_tmp <- pgs_wd2 %>% 
  filter(pgs_name == 'Base_1') %>%
  select(-c(pgs_name, PGS, pgs_pc_corrected)) %>% 
  distinct() %>% 
  pivot_longer(cols = matches('pgs_pc_corrected'), names_to = 'pgs_name', values_to = 'pgs_pc_corrected') %>%
  filter(sample_age_years_before_1950 > 100) %>%
  filter(! IID %in% drop_samp) %>%
  filter(str_detect(IID, pattern = 'Href|Gorilla|Chimp|Ancest|Denisova11_noUDG.SG', negate = TRUE)) %>%
  group_by(pgs_name) %>% 
  # filter(type == 'Ancient Europeans') %>%
  # filter(sample_age_years_before_1950 < 50000) %>%
  do(res = broom::tidy(cor.test(.$pgs_pc_corrected, .$sample_age_years_before_1950), method = 's')) %>% 
  # do(res = broom::tidy(cor.test(.$pgs_pc_corrected, log10((.$sample_age_years_before_1950 + 74) / 1000)), method = 'p')) %>% 
  unnest(res) %>% 
  arrange(p.value) %>%
  rename(name = pgs_name) %>% 
  mutate(clean_name = str_remove_all(name, pattern = 'pgs_pc_corrected_')) %>% 
  mutate(lab = str_c('r = ', round(-1 * estimate, digits = 2), ', p-val = ', formatC(p.value, digits = 2))) %>% ## flip direction to match time
  relocate(clean_name) %>% 
  rename(pgs_name = name) %>% 
  select(-pgs_name) %>% 
  drop_na()

res_sel_pelvis <- pgs_res_tmp %>%
  filter(str_detect(clean_name, 'waist|pw|piw')) %>%
  filter(str_detect(clean_name, 'pwgt', negate = TRUE))

pgs_res


pgs_res_birth %>% 
  filter(p.value < 0.05)
pgs_res_cog <- pgs_res_birth %>% 
  filter(str_detect(clean_name, 'Phoneme|Non|read|Spell'))

pgs_res_psych <- pgs_res_birth %>% 
  filter(str_detect(clean_name, 'ADHD|Addict|Anxiety|PTSD|Autism|Openn|Schiz|Bipolar|Conscie|Neurotic|Depress|Agree|Extra'))

pgs_res_birth <- pgs_res_birth %>% 
  filter(str_detect(clean_name, 'Phoneme|Non|read|Spell|ADHD|Addict|Anxiety|PTSD|Autism|Openn|Schiz|Bipolar|Conscie|Neurotic|Depress|Agree|Extra', negate = TRUE))


##
nean <- -1 * log10(100) ## neanderthal sample age
stone_age <- -1 * log10(3) ## stone age ended ~ 5,000 years ago
ren <- -1 * log10(.7) ## renaissance started ~ 700 years ago

pgs_wd2 <- pgs_wd2 %>% 
  group_by(pgs_name) %>% 
  filter(sample_age_years_before_1950 > 100) %>% 
  filter(str_detect(pgs_name, 'Base|^HAQER|^HAR')) %>% 
  filter(! IID %in% drop_samp) %>% 
  filter(str_detect(IID, pattern = 'Href|Gorilla|Chimp|Ancest|Denisova11_noUDG.SG', negate = TRUE))

##############################
## plot selection associations with pelvis measures
p_pelvis_evo_dat <- pgs_wd2 %>% 
  # group_by(pgs_name) %>% 
  filter(pgs_name == 'Base_1') %>% 
  filter(sample_age_years_before_1950 > 100) %>% 
  ungroup() %>%
  select(-c(pgs_name, PGS, pgs_pc_corrected)) %>% 
  pivot_longer(cols = matches('pgs_pc_corr')) %>%
  mutate(clean_name = str_remove_all(name, 'pgs_pc_corrected_')) %>%
  drop_na() %>% 
  inner_join(res_sel_pelvis) %>% 
  mutate(clean_name = case_when(str_detect(clean_name, 'fpw') ~ 'Waist circumference in females (BMI adj.) PGS',
                                str_detect(clean_name, 'waist_cir') ~ 'Waist circumference PGS',
                                str_detect(clean_name, 'fpiw') ~ 'Pelvic inlet width in females PGS',
                                TRUE ~ NA_character_)) %>% 
  mutate(pgs_name_clean = clean_name)
p_pelvis_evo_dat %>%
  rename(pgs_pc_corrected = value) %>%
  mutate(pgs_name_clean = str_remove_all(pgs_name_clean, ' PGS')) %>%
  mutate(pgs_name_clean = str_c(pgs_name_clean, ' (', lab, ')')) %>%
  group_by(pgs_name_clean) %>% 
  mutate(pgs_pc_corrected = scale(pgs_pc_corrected)[,1]) %>% 
  ungroup() %>%
  ggplot(aes(x = -1 * log10((sample_age_years_before_1950 + 74) / 1000), y = pgs_pc_corrected)) +
  geom_smooth(size = 1.5, aes(color = pgs_name_clean), alpha = 0.15) +
  scale_x_continuous(name = 'Years ago', breaks = c(seq(-2, 1, by = 1)), labels = c('100,000', '10,000', '1,000', '100')) + ## as.character(rev(10^seq(-1,2, by = 1)) * 1000)
  ylab('PGS') +
  geom_text(aes(x = nean + 0.25, y = 5, label = 'Neanderthals'), size  = 6, check_overlap = TRUE, color = 'grey50') +
  geom_text(aes(x = stone_age, y = 5, label = 'Stone age ends'), size  = 6, check_overlap = TRUE, color = 'grey50') +
  # geom_text(aes(x = -.5, y = 2, label = lab_bw$lab), check_overlap = TRUE, size  = 6, color = 'turquoise4') +
  # geom_text(aes(x = -.5, y = -2.5, label = lab_pwgt$lab), check_overlap = TRUE, size  = 6, color = 'cornflowerblue') +
  # geom_text(aes(x = -.5, y = 2, label = lab_bhc$lab), check_overlap = TRUE, size  = 6, color = 'midnightblue') +
  # geom_text(aes(x = -.5, y = 2, label = lab_ihc$lab), check_overlap = TRUE, size  = 6, color = 'darkgoldenrod') +
  theme_classic() +
  theme(axis.text = element_text(size = 16),
        axis.title = element_text(size = 18),
        legend.position = c(0.65, 0.6),
        legend.text = element_text(size = 10)) +
  labs(color = NULL)



##
p_evo_dat <- pgs_corrected %>% 
  filter(IID %in% kg_samp | IID %in% soi_qc) %>% 
  mutate(IID = case_when(IID == 'Chagyrskaya-Phalanx' ~ 'Chagyrskaya',
                         TRUE ~ IID)) %>%
  inner_join(smeta2) %>% 
  mutate(type = case_when(IID %in% kg_samp ~ '1000 Genomes Europeans',
                          IID %in% nean_samp ~ 'Neanderthals',
                          IID %in% unrel$X2 ~ 'Ancient Europeans',
                          IID %in% ref_samp ~ 'Primate reference genomes')) %>% 
  inner_join(smeta2) %>% 
  group_by(pgs_name) %>% 
  filter(sample_age_years_before_1950 > 100) %>%
  filter(! IID %in% drop_samp) %>%
  filter(str_detect(IID, pattern = 'Href|Gorilla|Chimp|Ancest|Denisova11_noUDG.SG', negate = TRUE))%>% 
  inner_join(pgs_res) %>%
  mutate(pgs_name_clean = case_when(pgs_name == 'Base_1' ~ 'Genome wide CP-PGS',
                                    pgs_name == 'HAQER_1' ~ 'HAQER CP-PGS',
                                    pgs_name == 'HAR_1' ~ 'HAR CP-PGS',
                                    pgs_name == 'complement_HAQER_1' ~ 'Background CP-PGS')) 

p_evo_dat2 <- pgs_wd2 %>% 
  # group_by(pgs_name) %>% 
  filter(pgs_name == 'Base_1') %>% 
  filter(sample_age_years_before_1950 > 100) %>% 
  ungroup() %>%
  select(-c(pgs_name, PGS, pgs_pc_corrected)) %>% 
  pivot_longer(cols = matches('pgs_pc_corr')) %>%
  mutate(clean_name = case_when(str_detect(name, '_bw$') ~ 'Birth weight PGS',
                                str_detect(name, '_ihc$') ~ 'Infant head circumference PGS',
                                str_detect(name, '_wr$') ~ 'Word reading PGS',
                                str_detect(name, '_nread$') ~ 'Non-word reading PGS',
                                str_detect(name, '_icv$') ~ 'Adult intracranial volume PGS',
                                str_detect(name, '_nrep$') ~ 'Non-word repetition PGS',
                                str_detect(name, '_pa$') ~ 'Phoneme awareness PGS',
                                str_detect(name, '_nviq$') ~ 'Nonverbal IQ PGS',
                                str_detect(name, '_sp$') ~ 'Spelling PGS',
                                str_detect(name, '_gd$') ~ 'Gestational duration PGS',
                                str_detect(name, '_bhc$') ~ 'Birth head circumference PGS',
                                TRUE ~ NA_character_)) %>% 
  drop_na() %>% 
  inner_join(pgs_res_birth) %>%
  rename(pgs_name = clean_name) %>% 
  filter(! IID %in% drop_samp) %>%
  mutate(pgs_name_clean = pgs_name)

p_evo_dat_cog <- pgs_wd2 %>% 
  # group_by(pgs_name) %>% 
  filter(pgs_name == 'Base_1') %>% 
  filter(sample_age_years_before_1950 > 100) %>% 
  ungroup() %>%
  select(-c(pgs_name, PGS, pgs_pc_corrected)) %>% 
  pivot_longer(cols = matches('pgs_pc_corr')) %>%
  mutate(clean_name = case_when(str_detect(name, '_bw$') ~ 'Birth weight PGS',
                                str_detect(name, '_ihc$') ~ 'Infant head circumference PGS',
                                str_detect(name, '_wr$') ~ 'Word reading PGS',
                                str_detect(name, '_nread$') ~ 'Non-word reading PGS',
                                str_detect(name, '_icv$') ~ 'Adult intracranial volume PGS',
                                str_detect(name, '_nrep$') ~ 'Non-word repetition PGS',
                                str_detect(name, '_pa$') ~ 'Phoneme awareness PGS',
                                str_detect(name, '_nviq$') ~ 'Nonverbal IQ PGS',
                                str_detect(name, '_sp$') ~ 'Spelling PGS',
                                str_detect(name, '_gd$') ~ 'Gestational duration PGS',
                                str_detect(name, '_bhc$') ~ 'Birth head circumference PGS',
                                TRUE ~ NA_character_)) %>% 
  drop_na() %>% 
  inner_join(pgs_res_cog) %>%
  rename(pgs_name = clean_name) %>% 
  filter(! IID %in% drop_samp) %>%
  mutate(pgs_name_clean = pgs_name)


p_evo_dat_psych <- pgs_wd2 %>% 
  # group_by(pgs_name) %>% 
  filter(pgs_name == 'Base_1') %>% 
  filter(sample_age_years_before_1950 > 100) %>% 
  ungroup() %>%
  select(-c(pgs_name, PGS, pgs_pc_corrected)) %>% 
  pivot_longer(cols = matches('pgs_pc_corr')) %>%
  mutate(clean_name = case_when(str_detect(name, '_ADHD$') ~ 'ADHD PGS',
                                str_detect(name, '_addiction$') ~ 'Addiction PGS',
                                str_detect(name, '_anxiety$') ~ 'Anxiety PGS',
                                str_detect(name, '_PTSD$') ~ 'PTSD PGS',
                                str_detect(name, '_autism$') ~ 'Autism PGS',
                                str_detect(name, '_open$') ~ 'Openness PGS',
                                str_detect(name, '_schizophrenia$') ~ 'Schizophrenia PGS',
                                str_detect(name, '_bipolar$') ~ 'Bipolar PGS',
                                str_detect(name, '_consc$') ~ 'Conscientiousness PGS',
                                str_detect(name, '_Neuro$') ~ 'Neuroticism PGS',
                                str_detect(name, '_depression$') ~ 'Depression PGS',
                                str_detect(name, '_agree$') ~ 'Agreeableness PGS',
                                str_detect(name, '_extra$') ~ 'Extraversion PGS',
                                TRUE ~ NA_character_)) %>% 
  drop_na() %>% 
  inner_join(pgs_res_psych) %>%
  rename(pgs_name = clean_name) %>% 
  filter(! IID %in% drop_samp) %>%
  mutate(pgs_name_clean = pgs_name)


range(p_evo_dat$pgs_pc_corrected)
# scale_color_manual(values = c('#dcc695', '#e04468', '#762776')) +

####################
## make plot with only lines and no points
lab_haq <- pgs_res %>% 
  filter(pgs_name == 'HAQER_1')
lab_gw <- pgs_res %>% 
  filter(pgs_name == 'complement_HAQER_1')

p_sel_summary <- p_evo_dat %>%
  filter(str_detect(pgs_name_clean, 'HAR|Base|Genome', negate = TRUE)) %>%
  mutate(pgs_name_clean = str_remove_all(pgs_name_clean, ' CP-PGS')) %>%
  ggplot(aes(x = -1 * log10((sample_age_years_before_1950 + 74) / 1000), y = pgs_pc_corrected)) +
  geom_smooth(size = 1.5, aes(color = pgs_name_clean), alpha = 0.15) +
  # facet_wrap(~ pgs_name) +
  scale_x_continuous(name = 'Years ago', breaks = c(seq(-2, 1, by = 1)), labels = c('100,000', '10,000', '1,000', '100')) + ## as.character(rev(10^seq(-1,2, by = 1)) * 1000)
  ylab('CP-PGS') +
  geom_text(aes(x = nean + 0.25, y = 3.6, label = 'Neanderthals'), size  = 6, check_overlap = TRUE, color = 'grey50') +
  geom_text(aes(x = stone_age, y = 3.6, label = 'Stone age ends'), size  = 6, check_overlap = TRUE, color = 'grey50') +
  # geom_text(aes(x = -0.75, y = 4.2, label = lab), check_overlap = TRUE, size  = 6) +
  geom_text(aes(x = -.5, y = 2, label = lab_haq$lab), check_overlap = TRUE, size  = 6, color = '#00bfc4') +
  geom_text(aes(x = -.5, y = -2.5, label = lab_gw$lab), check_overlap = TRUE, size  = 6, color = '#f8766d') +
  theme_classic() +
  # theme(axis.text = element_text(size = 16),
  #       axis.title = element_text(size = 18),
  #       legend.position = "inside",
  #       legend.position.inside = c(0.75, 0.1)) +
  theme(axis.text = element_text(size = 16),
        axis.title = element_text(size = 18),
        legend.position = c(0.75, 0.1)) +
  labs(color = NULL) +
  scale_color_manual(values = c('#f8766d', '#00bfc4'))

####
## repeat for birth stuff\
unique(pgs_res_birth$clean_name)
lab_bw <- pgs_res_birth %>% 
  filter(clean_name == 'Birth weight PGS')
lab_pwgt <- pgs_res_birth %>% 
  filter(clean_name == 'Placental weight PGS')
lab_ihc <- pgs_res_birth %>% 
  filter(clean_name == 'Infant head circumference PGS')
lab_bhc <- pgs_res_birth %>% 
  filter(clean_name == 'Birth head circumference PGS')


p_sel_summary_birth <- p_evo_dat2 %>%
  rename(pgs_pc_corrected = value) %>%
  mutate(pgs_name_clean = str_remove_all(pgs_name_clean, ' PGS')) %>%
  mutate(pgs_name_clean = str_c(pgs_name_clean, ' (', lab, ')')) %>%
  group_by(pgs_name) %>% 
  mutate(pgs_pc_corrected = scale(pgs_pc_corrected)[,1]) %>% 
  ungroup() %>%
  filter(str_detect(pgs_name_clean, 'verbal|reading|Non|honeme|pelling|Infant', negate = TRUE)) %>% 
  ggplot(aes(x = -1 * log10((sample_age_years_before_1950 + 74) / 1000), y = pgs_pc_corrected)) +
  geom_smooth(size = 1.5, aes(color = pgs_name_clean), alpha = 0.15) +
  scale_x_continuous(name = 'Years ago', breaks = c(seq(-2, 1, by = 1)), labels = c('100,000', '10,000', '1,000', '100')) + ## as.character(rev(10^seq(-1,2, by = 1)) * 1000)
  ylab('PGS') +
  geom_text(aes(x = nean + 0.25, y = 7, label = 'Neanderthals'), size  = 6, check_overlap = TRUE, color = 'grey50') +
  geom_text(aes(x = stone_age, y = 7, label = 'Stone age ends'), size  = 6, check_overlap = TRUE, color = 'grey50') +
  # geom_text(aes(x = -.5, y = 2, label = lab_bw$lab), check_overlap = TRUE, size  = 6, color = 'turquoise4') +
  # geom_text(aes(x = -.5, y = -2.5, label = lab_pwgt$lab), check_overlap = TRUE, size  = 6, color = 'cornflowerblue') +
  # geom_text(aes(x = -.5, y = 2, label = lab_bhc$lab), check_overlap = TRUE, size  = 6, color = 'midnightblue') +
  # geom_text(aes(x = -.5, y = 2, label = lab_ihc$lab), check_overlap = TRUE, size  = 6, color = 'darkgoldenrod') +
  theme_classic() +
  theme(axis.text = element_text(size = 16),
        axis.title = element_text(size = 18),
        legend.position = c(0.65, 0.6),
        legend.text = element_text(size = 10)) +
  labs(color = NULL)


p_sel_summary_lang <- p_evo_dat_cog %>%
  rename(pgs_pc_corrected = value) %>%
  mutate(pgs_name_clean = str_remove_all(pgs_name_clean, ' PGS')) %>%
  mutate(pgs_name_clean = str_c(pgs_name_clean, ' (', lab, ')')) %>%
  group_by(pgs_name) %>% 
  mutate(pgs_pc_corrected = scale(pgs_pc_corrected)[,1]) %>% 
  ungroup() %>%
  filter(str_detect(pgs_name_clean, 'verbal|reading|Non|honeme|pelling')) %>% 
  ggplot(aes(x = -1 * log10((sample_age_years_before_1950 + 74) / 1000), y = pgs_pc_corrected)) +
  geom_smooth(size = 1.5, aes(color = pgs_name_clean), alpha = 0.15) +
  scale_x_continuous(name = 'Years ago', breaks = c(seq(-2, 1, by = 1)), labels = c('100,000', '10,000', '1,000', '100')) + ## as.character(rev(10^seq(-1,2, by = 1)) * 1000)
  ylab('PGS') +
  geom_text(aes(x = nean + 0.25, y = 7, label = 'Neanderthals'), size  = 6, check_overlap = TRUE, color = 'grey50') +
  geom_text(aes(x = stone_age, y = 7, label = 'Stone age ends'), size  = 6, check_overlap = TRUE, color = 'grey50') +
  # geom_text(aes(x = -.5, y = 2, label = lab_bw$lab), check_overlap = TRUE, size  = 6, color = 'turquoise4') +
  # geom_text(aes(x = -.5, y = -2.5, label = lab_pwgt$lab), check_overlap = TRUE, size  = 6, color = 'cornflowerblue') +
  # geom_text(aes(x = -.5, y = 2, label = lab_bhc$lab), check_overlap = TRUE, size  = 6, color = 'midnightblue') +
  # geom_text(aes(x = -.5, y = 2, label = lab_ihc$lab), check_overlap = TRUE, size  = 6, color = 'darkgoldenrod') +
  theme_classic() +
  theme(axis.text = element_text(size = 16),
        axis.title = element_text(size = 18),
        legend.position = c(0.65, 0.12),
        legend.text = element_text(size = 10)) +
  labs(color = NULL)


p_sel_summary_psych <- p_evo_dat_psych %>%
  rename(pgs_pc_corrected = value) %>%
  mutate(pgs_name_clean = str_remove_all(pgs_name_clean, ' PGS')) %>%
  mutate(pgs_name_clean = str_c(pgs_name_clean, ' (', lab, ')')) %>%
  group_by(pgs_name) %>% 
  mutate(pgs_pc_corrected = scale(pgs_pc_corrected)[,1]) %>% 
  ungroup() %>%
  filter(str_detect(pgs_name_clean, 'Open|Agree|Neuro|Consc|Extra', negate = TRUE)) %>% 
  ggplot(aes(x = -1 * log10((sample_age_years_before_1950 + 74) / 1000), y = pgs_pc_corrected)) +
  # ggplot(aes(x = -1 * sample_age_years_before_1950, y = pgs_pc_corrected)) +
  geom_smooth(size = 1.5, aes(color = pgs_name_clean), alpha= 0.15) +
  scale_x_continuous(name = 'Years ago', breaks = c(seq(-2, 1, by = 1)), labels = c('100,000', '10,000', '1,000', '100')) + ## as.character(rev(10^seq(-1,2, by = 1)) * 1000)
  ylab('PGS') +
  geom_text(aes(x = nean + 0.25, y = 7, label = 'Neanderthals'), size  = 6, check_overlap = TRUE, color = 'grey50') +
  geom_text(aes(x = stone_age, y = 7, label = 'Stone age ends'), size  = 6, check_overlap = TRUE, color = 'grey50') +
  # geom_text(aes(x = -.5, y = 2, label = lab_bw$lab), check_overlap = TRUE, size  = 6, color = 'turquoise4') +
  # geom_text(aes(x = -.5, y = -2.5, label = lab_pwgt$lab), check_overlap = TRUE, size  = 6, color = 'cornflowerblue') +
  # geom_text(aes(x = -.5, y = 2, label = lab_bhc$lab), check_overlap = TRUE, size  = 6, color = 'midnightblue') +
  # geom_text(aes(x = -.5, y = 2, label = lab_ihc$lab), check_overlap = TRUE, size  = 6, color = 'darkgoldenrod') +
  theme_classic() +
  theme(axis.text = element_text(size = 16),
        axis.title = element_text(size = 18),
        legend.position = c(0.65, 0.12),
        legend.text = element_text(size = 10)) +
  labs(color = NULL)

p_sel_summary_pers <- p_evo_dat_psych %>%
  rename(pgs_pc_corrected = value) %>%
  mutate(pgs_name_clean = str_remove_all(pgs_name_clean, ' PGS')) %>%
  mutate(pgs_name_clean = str_c(pgs_name_clean, ' (', lab, ')')) %>%
  group_by(pgs_name) %>% 
  mutate(pgs_pc_corrected = scale(pgs_pc_corrected)[,1]) %>% 
  ungroup() %>%
  filter(str_detect(pgs_name_clean, 'Open|Agree|Neuro|Consc|Extra')) %>% 
  ggplot(aes(x = -1 * log10((sample_age_years_before_1950 + 74) / 1000), y = pgs_pc_corrected)) +
  # ggplot(aes(x = -1 * sample_age_years_before_1950, y = pgs_pc_corrected)) +
  geom_smooth(size = 1.5, aes(color = pgs_name_clean), alpha= 0.15) +
  scale_x_continuous(name = 'Years ago', breaks = c(seq(-2, 1, by = 1)), labels = c('100,000', '10,000', '1,000', '100')) + ## as.character(rev(10^seq(-1,2, by = 1)) * 1000)
  ylab('PGS') +
  geom_text(aes(x = nean + 0.25, y = 7, label = 'Neanderthals'), size  = 6, check_overlap = TRUE, color = 'grey50') +
  geom_text(aes(x = stone_age, y = 7, label = 'Stone age ends'), size  = 6, check_overlap = TRUE, color = 'grey50') +
  # geom_text(aes(x = -.5, y = 2, label = lab_bw$lab), check_overlap = TRUE, size  = 6, color = 'turquoise4') +
  # geom_text(aes(x = -.5, y = -2.5, label = lab_pwgt$lab), check_overlap = TRUE, size  = 6, color = 'cornflowerblue') +
  # geom_text(aes(x = -.5, y = 2, label = lab_bhc$lab), check_overlap = TRUE, size  = 6, color = 'midnightblue') +
  # geom_text(aes(x = -.5, y = 2, label = lab_ihc$lab), check_overlap = TRUE, size  = 6, color = 'darkgoldenrod') +
  theme_classic() +
  theme(axis.text = element_text(size = 16),
        axis.title = element_text(size = 18),
        legend.position = c(0.65, 0.12),
        legend.text = element_text(size = 10)) +
  labs(color = NULL)

#####################
## try association with age of death
pgs_wd2 %>% 
  group_by(pgs_name) %>% 
  filter(sample_age_years_before_1950 > 100) %>%
  ungroup() %>%
  mutate(death_age = str_split(age_at_death, '[;] ', simplify = TRUE)[,2],
         death_age = str_split(death_age, ' ', simplify = TRUE)[,1],
         death_age = str_split(death_age, '[-]', simplify = TRUE)[,2],
         death_age = str_remove_all(death_age, '[A-Z]|[a-z]'), #) %>% distinct(IID, death_age) %>% select(death_age) %>% table()#,
         death_age = as.numeric(death_age)) %>%
  drop_na(death_age) %>%
  # filter(death_age >= 15) %>%
  # filter(type == 'Ancient Europeans') %>%
  group_by(pgs_name) %>%
  do(res = broom::tidy(cor.test(.$pgs_pc_corrected, .$death_age, method = 'p'))) %>% 
  unnest(res) %>% 
  arrange(p.value) %>% 
  filter(str_detect(pgs_name, 'Base|^HAQER|^HAR')) %>% 
  mutate(lab = str_c('r = ', round(-1 * estimate, digits = 2), ', p-val = ', formatC(p.value, digits = 2)))

pgs_wd2 %>% 
  ungroup() %>%
  filter(pgs_name == 'Base_1') %>% 
  select(-c(pgs_name, PGS, pgs_pc_corrected)) %>% 
  pivot_longer(cols = matches('pgs_pc'), names_to = 'pgs_name', values_to = 'pgs_pc_corrected') %>% 
  mutate(pgs_name = str_remove_all(pgs_name, 'pgs_pc_corrected_')) %>%
  group_by(pgs_name) %>% 
  filter(sample_age_years_before_1950 > 100) %>%
  ungroup() %>%
  mutate(death_age = str_split(age_at_death, '[;] ', simplify = TRUE)[,2],
         death_age = str_split(death_age, ' ', simplify = TRUE)[,1],
         death_age = str_split(death_age, '[-]', simplify = TRUE)[,2],
         death_age = str_remove_all(death_age, '[A-Z]|[a-z]'), #) %>% distinct(IID, death_age) %>% select(death_age) %>% table()#,
         death_age = as.numeric(death_age)) %>%
  drop_na(death_age) %>%
  # filter(death_age >= 15) %>%
  # filter(type == 'Ancient Europeans') %>%
  group_by(pgs_name) %>%
  do(res = broom::tidy(cor.test(.$pgs_pc_corrected, .$death_age, method = 's'))) %>% 
  unnest(res) %>% 
  arrange(p.value) %>% 
  # filter(str_detect(pgs_name, 'Base|^HAQER|^HAR')) %>% 
  mutate(lab = str_c('r = ', round(-1 * estimate, digits = 2), ', p-val = ', formatC(p.value, digits = 2)))


###############################################
## look for associations with archaic ancestry
###############################################
pc_anc <- pca %>% 
  filter(! IID %in% kg_pc$IID)
aa <- read_csv('/Dedicated/jmichaelson-wdata/lcasten/sli_wgs/ancient_DNA/data/archaic_ancestry_scores_RinkerAlleles.csv') # %>% 
#  filter(NDA_ancestry < 0.125)
aa$nda_resid_ancestry <- scale(resid(lm(NDA_ancestry ~ RHA_ancestry + RAA_ancestry, data = aa)))[,1]
aa %>% 
  filter(NDA_ancestry > .25)
hist(aa$nda_resid_ancestry)

aa_long <- aa %>% 
  pivot_longer(cols = matches('ancestry'))

# aa_long <- aa_long %>% 
#   inner_join(smeta2) %>% 
#   filter(IID %in% unrel$X2) %>% 
#   mutate(lat = as.numeric(lat),
#          long = as.numeric(long)) %>% 
#   drop_na(lat, long, value) %>% 
#   group_by(name) %>% 
#   mutate(value_loc_adj = scale(resid(lm(value ~ lat + long)))[,1]) %>% 
#   relocate(value_loc_adj, .after = value)

aa_long %>% 
  # select(-c(lat, long)) %>%
  inner_join(pgs_wd2) %>% 
  filter(sample_age_years_before_1950 > 100) %>%
  filter(value < 0.25) %>%
  # filter(sample_age_years_before_1950 < 10000) %>%
  inner_join(pc_anc) %>%
  group_by(name) %>%
  mutate(resid_anc = scale(resid(lm(value ~ pc1 + pc2 + pc3 + pc4 + pc5)))[,1]) %>%
  group_by(pgs_name, name) %>%
  do(res = broom::tidy(cor.test(.$pgs_pc_corrected, .$value, method = 'p'))) %>% 
  unnest(res) %>% 
  arrange(p.value) %>% 
  filter(str_detect(pgs_name, 'Base|^HAQER|^HAR')) %>% 
  mutate(lab = str_c('r = ', round(-1 * estimate, digits = 2), ', p-val = ', formatC(p.value, digits = 2)))

aa_long %>% 
  # select(-c(lat, long)) %>%
  inner_join(pgs_wd2) %>% 
  filter(sample_age_years_before_1950 > 100) %>%
  # filter(sample_age_years_before_1950 < 10000) %>%
  inner_join(pc_anc) %>% 
  filter(name == 'NDA_ancestry') %>% 
  filter(pgs_name == 'HAQER_1') %>% 
  ggplot(aes(x = pgs_pc_corrected, y = value )) +
  geom_point() +
  geom_smooth()

aa_age <- aa_long %>% 
  inner_join(smeta2) %>% 
  filter(IID %in% pgs_wd2$IID) %>%
  filter(sample_age_years_before_1950 > 100) 

aa_age %>%
  group_by(name) %>%
  do(res = broom::tidy(cor.test(log10(.$sample_age_years_before_1950), .$value, method = 'p'))) %>% 
  unnest(res) %>% 
  arrange(p.value) %>% 
  mutate(lab = str_c('r = ', round(-1 * estimate, digits = 2), ', p-val = ', formatC(p.value, digits = 2)))

range(aa_age$sample_age_years_before_1950)
aa_age %>%
  inner_join(pc_anc) %>%
  filter(sample_age_years_before_1950 < 14000) %>%
  group_by(name) %>%
  mutate(resid_anc = scale(resid(lm(value ~ pc1 + pc2 + pc3 + pc4 + pc5)))[,1]) %>%
  group_by(name) %>%
  do(res = broom::tidy(cor.test(log10(.$sample_age_years_before_1950), .$value, method = 'p'))) %>% 
  unnest(res) %>% 
  arrange(p.value) %>% 
  mutate(lab = str_c('r = ', round(-1 * estimate, digits = 2), ', p-val = ', formatC(p.value, digits = 2)))


aa_age %>%
  filter(name == 'NDA_ancestry') %>%
  arrange(desc(log10(sample_age_years_before_1950))) %>% 
  # filter(sample_age_years_before_1950 < 10000) %>%
  ggplot(aes(x = value, y = log10(sample_age_years_before_1950))) +
  geom_point() +
  geom_smooth()

##############################
## look at brain size changes over time
#############################
bs <- readxl::read_xlsx('/Dedicated/jmichaelson-wdata/lcasten/sli_wgs/ancient_DNA/data/DeSilva2023_ancient_hominin_brain_size.xlsx') %>% 
  janitor::clean_names()
bs2 <- bs %>% 
  mutate(age_ma = as.numeric(age_ma)) %>%
  mutate(years_ago = age_ma * 1000000) %>% 
  select(species, specimen, cc, years_ago) %>% 
  drop_na() %>% 
  filter(years_ago < 100000 & years_ago > 100 & species != 'H. floresiensis')

table(bs2$species)
brain_cor <- broom::tidy(cor.test(-1 * log10(bs2$years_ago), bs2$cc)) %>% 
  mutate(lab = str_c('r = ', round(estimate, digits = 2), ', p-val = ', formatC(p.value, digits = 2)))

p_bs <- bs2 %>% 
  ggplot(aes(x = (-1 * log10(years_ago / 1000)), y = cc)) +
  geom_jitter(size = 0.7) +
  geom_smooth(method = 'lm', color = 'black', size = 1.5) +
  scale_x_continuous(name = 'Years ago', breaks = c(seq(-2, 1, by = 1)), labels = c('100,000', '10,000', '1,000', '100')) + ## as.character(rev(10^seq(-1,2, by = 1)) * 1000)
  ylab('Brain size (CC)') +
  geom_text(aes(x = -0.75, y = 1900, label = brain_cor$lab), check_overlap = TRUE, size  = 6) +
  theme_classic() +
  theme(axis.text = element_text(size = 16),
        axis.title = element_text(size = 18))

p_sel_summary + p_sel_summary_birth + p_sel_summary_lang

##################################################
## try and show HAQER CP-PGS predicts ICV in ABCD and metabolism
###########
### Functions for a) processing and formatting ABCD raw text files and b) optionally saving the data as an Rdata file.

parseABCD = function(dir, shortname, viewcol=TRUE){
  
  # get the file path using input 'shortname'.
  fn = list.files(path=dir, pattern=shortname, recursive=TRUE, include.dirs=T, full.names=TRUE)
  
  # in the case of multiple matching files...
  if(length(fn)!=1){
    print(fn)
    user.inp = readline("multiple matches to shortname. Which file should I use?\nselection #: ")
    user.inp = as.numeric(user.inp)
    fn = fn[user.inp]
    print(paste("parsing: ", fn, "\n", sep=""))
  }
  
  #read in table.
  dat = read.table(fn, stringsAsFactors=F, header=T, na.strings=c("", 999, NA), encoding="UTF-8", sep="\t", fill=T)
  coldat = structure(dat[1,], names=colnames(dat))
  coldat = as.character(coldat[1,])
  dat = dat[-1,]
  rownames(dat) = NULL
  
  if(viewcol){
    print(head(coldat, 10))
    #error handling.
    user.inp = readline("Does this look correct? Shall I finish processing/saving the RDA file?\n(y/n): ")
    if(user.inp != 'y'){
      stop("aborted")
    }
  }
  
  #save rdata file.
  save(dat, coldat, file=gsub("\\.txt", ".rda" , fn))
}


## Function to read in abcd data text files
## Requires filepath as input
## Returns a two-element list with the data and description of column names (`dict`)
## Uses readr for (mostly) smart column typing
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
                          na = c('NA', '999' ,''), 
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
abcd_icv <- read_abcd("/Dedicated/jmichaelson-sdata/ABCD/abcd_release_5_0/abcd_smrip10201.txt")
abcd_icv$dict
abcd_icv$dict %>% 
  as_tibble() %>%
  filter(str_detect(description, 'audate|asal|allid|utamen|ippocamp|mygda|acc|ubcort')) %>% 
  as.data.frame()

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
  filter(interview_age <= interview_age_latest) %>% 
  filter(smri_vol_scs_intracranialv_latest - smri_vol_scs_intracranialv > 0) %>% 
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

# ##
# pgs_abcd <- read_csv('~/LSS/jmichaelson-wdata/lcasten/sli_wgs/prs/pathway_prs/ABCD_gathered_pgs_pc_corrected_long_full_data.cogPerf.complement.csv') %>%
#     filter(str_detect(pgs_name, '_1$'))
pgs_abcd <- read_csv('/wdata/lcasten/sli_wgs/prs/pathway_prs/SPARK_ABCD.gathered_pgs_pc_corrected_long_full_data.cogPerf.complement.csv')
# unrel <- read_table('/genome/SPARK_ABCD_imputed_genotype_merge/autosomes_all.GRM.unrelated_cutoff_0.05.singleton.txt', col_names = FALSE)
unrel <- read_table('/genome/SPARK_ABCD_imputed_genotype_merge/autosomes_all.GRM.unrelated_cutoff_0.05.singleton.txt', col_names = FALSE)

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
  # distinct(IID)
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
  mutate(resid_vol = resid(lm(smri_vol_scs_intracranialv ~ interview_age + anthroheightcalc + anthroweight1lb)),
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


icv_growth_res %>% 
  filter(str_detect(pgs_name, 'HAQER_1'))


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

# abcd_pc <- read_delim('/Dedicated/jmichaelson-wdata/lcasten/spark_abcd_array_merge/topmed_hg19/ABCD_unrelated_europeans.qc.prs.eigenvec', col_names = FALSE, delim = ' ')
kp <- read_table('/genome/SPARK_ABCD_imputed_genotype_merge/autosomes_all.GRM.unrelated_cutoff_0.05.singleton.txt', col_names = FALSE) %>% 
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
sib <- read_table('/genome/SPARK_ABCD_imputed_genotype_merge/autosomes_all.GRM.unrelated_cutoff_0.05.family.txt', col_names = FALSE)
sib <- sib %>% 
    filter(str_detect(X2, 'SP', negate = TRUE), str_detect(X4, 'SP', negate = TRUE)) %>% 
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

p_sib_haq + p_sib_bg

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
    rename(age_diff = age, premature_diff = premature, sex_female_diff = sex_female, pheno = name1) %>%
    write_csv('/wdata/lcasten/sli_wgs/ancient_DNA/data/ABCD_sibling_birth_phenotypes.csv')

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
    filter(sib1_premature - sib2_premature <= 3) %>% ## remove siblings born at different gestational ages
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

p_paired_haq + p_paired_bg



########################################################
## merge plots linking HAQERs to evolution + brain size
########################################################
layout <- "
AABBCCDD
EEFFGGHH
IIJJKKLL
"

p_merged <- p_bg_neand + p_haq_neand + p_sel_summary + p_sel_summary_birth + ## row 1
  p_sel_summary_lang + p_vol + p_pgs + p_pgs_gl + ## row 2
  p_bw + p_growth +p_paired_bg + p_paired_haq + ## row 3
  plot_layout(design = layout) + 
  plot_annotation(tag_levels = 'A') & theme(plot.tag = element_text(face = 'bold', size = 18), axis.text = element_text(size = 16), axis.title = element_text(size = 18, face = 'bold'), legend.text = element_text(size = 14), legend.title = element_text(size = 16), strip.text = element_text(size = 18))

ggsave(p_merged, filename = '/Dedicated/jmichaelson-wdata/lcasten/sli_wgs/paper_figures/HAQER_selective_pressures.png', device = 'png', units = 'in', dpi = 300, width = 38, height = 26)

p_bg_neand + p_haq_neand + p_sel_summary  & theme(plot.tag = element_text(face = 'bold', size = 18), axis.text = element_text(size = 16), axis.title = element_text(size = 18, face = 'bold'), 
  legend.text = element_text(size = 16), legend.title = element_text(size = 16), strip.text = element_text(size = 18))


p_sel_summary_birth + ## row 1
  p_sel_summary_lang & theme(plot.tag = element_text(face = 'bold', size = 18), axis.text = element_text(size = 16), axis.title = element_text(size = 18, face = 'bold'), 
  legend.text = element_text(size = 16), legend.title = element_text(size = 16), strip.text = element_text(size = 18))

p_pgs + p_pgs_gl

p_vol + p_growth& theme(plot.tag = element_text(face = 'bold', size = 18), axis.text = element_text(size = 16), axis.title = element_text(size = 18, face = 'bold'), 
  legend.text = element_text(size = 16), legend.title = element_text(size = 16), strip.text = element_text(size = 18))

p_bw + p_paired_bg + p_paired_haq

p_sel_summary_pers + p_sel_summary_psych & theme(plot.tag = element_text(face = 'bold', size = 18), axis.text = element_text(size = 16), axis.title = element_text(size = 18, face = 'bold'), 
  legend.text = element_text(size = 16), legend.title = element_text(size = 16), strip.text = element_text(size = 18))


# ###########################
# ##
# ###################
# ## metabolism assoc?
# ## male eq: 10 × weight (in kilograms) + 6.25 × height (in centimeters) – 5 × age (in years) + 5
# ## female eq: 10 × weight (in kilograms) + 6.25 × height (in centimeters) – 5 × age (in years) -161
# met <- abcd_height$data %>% 
#   select(subjectkey, interview_age, sex, anthroheightcalc, anthroweight1lb) %>%
#   drop_na() %>%
#   as_tibble() %>% 
#   arrange(desc(interview_age)) %>%
#   distinct(subjectkey, .keep_all = TRUE) %>%
#   filter(anthroheightcalc > 36) %>%
#   rename(IID = subjectkey) %>% 
#   mutate(wt_kg = anthroweight1lb / 2.205,
#          ht_cm = anthroheightcalc * 2.54,
#          age_yrs = interview_age / 12,
#          offset = case_when(sex == 'F' ~ -161,
#                             sex == 'M' ~ 5,
#                             TRUE ~ NA)) %>%
#   mutate(bmr = 10 * wt_kg + 6.25 * ht_cm - 5 * age_yrs + offset)

# met %>% filter(bmr < 500)
# hist(met$bmr)

# met %>% 
#   inner_join(pgs_abcd) %>% 
#   group_by(pgs_name) %>%
#   do(res = broom::tidy(cor.test(.$bmr, .$pgs_pc_corrected))) %>%
#   unnest(res) %>%
#   arrange(p.value)


# ##
# bw <- read_abcd("/Dedicated/jmichaelson-sdata/ABCD/abcd_release_4_0/dhx01.txt")$data %>% 
#   select(subjectkey, sex, interview_age, 12,14) %>% 
#   drop_na() %>%
#   as_tibble() %>%
#   mutate(oz = 16 * birth_weight_lbs,
#          birth_wt = oz + birth_weight_oz) %>%
#   arrange(desc(interview_age)) %>%
#   distinct(subjectkey, .keep_all = TRUE) %>%
#   rename(IID = subjectkey)

# ####################################
# ## look at neanderthal  ancestry
# ####################################
# ##
# nda <-read_table('/Dedicated/jmichaelson-wdata/lcasten/sli_wgs/evolution/NeanderScore_common.Rinker_class-NDA.profile') %>% 
#   mutate(SCORESUM = SCORESUM / CNT) %>%
#   select(IID, SCORESUM_NDA = SCORESUM)
# rha <- read_table('/Dedicated/jmichaelson-wdata/lcasten/sli_wgs/evolution/NeanderScore_common.Rinker_class-RHA.profile') %>% 
#   mutate(SCORESUM = SCORESUM / CNT) %>%
#   select(IID, SCORESUM_RHA = SCORESUM)
# raa <- read_table('/Dedicated/jmichaelson-wdata/lcasten/sli_wgs/evolution/NeanderScore_common.Rinker_class-RAA.profile') %>% 
#   mutate(SCORESUM = SCORESUM / CNT) %>%
#   select(IID, SCORESUM_RAA = SCORESUM)
# rns <- read_csv('/Dedicated/jmichaelson-wdata/lcasten/sli_wgs/evolution/NeanderScore_rare_missense.csv')
# pgs 

# nda %>%
#   inner_join(rha) %>%
#   inner_join(raa) %>%
#   inner_join(pgs_corrected) %>%
#   group_by(pgs_name) %>%
#   do(res = broom::tidy(cor.test(.$SCORESUM_NDA, .$pgs_pc_corrected, method = 's'))) %>%
#   unnest(res) %>%
#   arrange(p.value) %>%
#   filter(str_detect(pgs_name, 'HAQER|HAR|Base'))







############################
# p_evo_gw <- p_evo_dat %>%
#   filter(pgs_name_clean == 'Genome wide CP-PGS') %>%
#   ggplot(aes(x = -1 * log10((sample_age_years_before_1950 + 74) / 1000), y = pgs_pc_corrected)) +
#   geom_jitter(size = 0.7) +
#   geom_smooth(method = 'lm', size = 1.5, color = 'black') +
#   # facet_wrap(~ pgs_name) +
#   scale_x_continuous(name = 'Years ago', breaks = c(seq(-2, 1, by = 1)), labels = c('100,000', '10,000', '1,000', '100')) + ## as.character(rev(10^seq(-1,2, by = 1)) * 1000)
#   ylab('Genome wide CP-PGS') +
#   geom_text(aes(x = nean + 0.25, y = 3.7, label = 'Neanderthals'), size = 3, check_overlap = TRUE, color = 'grey50') +
#   geom_text(aes(x = stone_age, y = 3.7, label = 'Stone age ends'), size = 2.5, check_overlap = TRUE, color = 'grey50') +
#   # geom_text(aes(x = ren, y = 3.7, label = 'Renaissance'), size = 2.5, check_overlap = TRUE, color = 'grey50') +
#   geom_text(aes(x = -0.75, y = 4.2, label = lab), check_overlap = TRUE, size  = 6) +
#   theme_classic() +
#   theme(axis.text = element_text(size = 16),
#         axis.title = element_text(size = 18)) +
#   coord_cartesian(ylim = c(-4.179814, 4.553474))

# p_evo_haq <- p_evo_dat %>%
#   filter(pgs_name_clean == 'HAQER CP-PGS') %>%
#   ggplot(aes(x = -1 * log10((sample_age_years_before_1950 + 74) / 1000), y = pgs_pc_corrected)) +
#   geom_jitter(color = '#762776', size = 0.7) +
#   geom_smooth(method = 'lm', size = 1.5, color = 'black') +
#   # facet_wrap(~ pgs_name) +
#   scale_x_continuous(name = 'Years ago', breaks = c(seq(-2, 1, by = 1)), labels = c('100,000', '10,000', '1,000', '100')) + ## as.character(rev(10^seq(-1,2, by = 1)) * 1000)
#   ylab('HAQER CP-PGS') +
#   geom_text(aes(x = nean + 0.25, y = 3.7, label = 'Neanderthals'), size = 3, check_overlap = TRUE, color = 'grey50') +
#   geom_text(aes(x = stone_age, y = 3.7, label = 'Stone age ends'), size = 2.5, check_overlap = TRUE, color = 'grey50') +
#   # geom_text(aes(x = ren, y = 3.7, label = 'Renaissance'), size = 2.5, check_overlap = TRUE, color = 'grey50') +
#   geom_text(aes(x = -0.75, y = 4.2, label = lab), check_overlap = TRUE, size  = 6) +
#   theme_classic() +
#   theme(axis.text = element_text(size = 16),
#         axis.title = element_text(size = 18))+
#   coord_cartesian(ylim = c(-4.179814, 4.553474))

# p_evo_har <- p_evo_dat %>%
#   filter(pgs_name_clean == 'HAR CP-PGS') %>%
#   ggplot(aes(x = -1 * log10((sample_age_years_before_1950 + 74) / 1000), y = pgs_pc_corrected)) +
#   geom_jitter(color = '#e04468', size = 0.7) +
#   geom_smooth(method = 'lm', size = 1.5, color = 'black') +
#   # facet_wrap(~ pgs_name) +
#   scale_x_continuous(name = 'Years ago', breaks = c(seq(-2, 1, by = 1)), labels = c('100,000', '10,000', '1,000', '100')) + ## as.character(rev(10^seq(-1,2, by = 1)) * 1000)
#   ylab('HAR CP-PGS') +
#   geom_text(aes(x = nean + 0.25, y = 3.7, label = 'Neanderthals'), size = 3, check_overlap = TRUE, color = 'grey50') +
#   geom_text(aes(x = stone_age, y = 3.7, label = 'Stone age ends'), size = 3, check_overlap = TRUE, color = 'grey50') +
#   # geom_text(aes(x = ren, y = 3.7, label = 'Renaissance'), size = 3, check_overlap = TRUE, color = 'grey50') +
#   geom_text(aes(x = -0.75, y = 4.2, label = lab), check_overlap = TRUE, size  = 6) +
#   theme_classic() +
#   theme(axis.text = element_text(size = 16),
#         axis.title = element_text(size = 18))+
#   coord_cartesian(ylim = c(-4.179814, 4.553474))

# p <- p_evo_gw + p_evo_haq + p_evo_har + 
#   plot_annotation(tag_levels = 'A') & theme(plot.tag = element_text(face = 'bold', size = 18), axis.text = element_text(size = 16), axis.title = element_text(size = 18), legend.text = element_text(size = 16), legend.title = element_text(size = 18), strip.text = element_text(size = 18))
# p
# ggsave(p, filename = '/Dedicated/jmichaelson-wdata/lcasten/sli_wgs/paper_figures/AADR_neanderthals_CP-ES-PGS_evo.png', device = 'png', units = 'in', dpi = 300, width = 12, height = 4)
