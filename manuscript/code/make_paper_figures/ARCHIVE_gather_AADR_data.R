library(tidyverse)

##
smeta <- read_tsv('/Dedicated/jmichaelson-wdata/lcasten/tools/ref_data/allen_ancient_dna_resource/v54.1.p1_1240K_public/v54.1.p1_1240K_public.anno') %>%
  janitor::clean_names() %>%
  rename(IID = genetic_id) %>% 
  select(IID, group_id, locality, political_entity, lat, long, molecular_sex, sample_age_years_before_1950 = date_mean_in_bp_in_years_before_1950_ce_ox_cal_mu_for_a_direct_radiocarbon_date_and_average_of_range_for_a_contextual_date, age_at_death = age_at_death_from_physical_anthropology)

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
pgs_corrected %>% 
  filter(pgs_name %in% c('HAQER_1', 'complement_HAQER_1')) %>% 
  filter(IID %in% kg_samp) %>% 
  mutate(pgs_name = case_when(pgs_name == 'HAQER_1' ~ 'cp_pgs.HAQER',
                              pgs_name == 'complement_HAQER_1' ~ 'cp_pgs.background_HAQER')) %>%
  select(IID, pgs_name, pgs_pc_corrected) %>% 
  pivot_wider(id_cols = IID, names_from = pgs_name, values_from = pgs_pc_corrected) %>% 
  write_csv('manuscript/supplemental_materials/1000_genomes_eur_ES-PGS_with_AADR_SNPs_data.csv')

list.files('/Dedicated/jmichaelson-wdata/lcasten/sli_wgs/ancient_DNA/data', pattern = 'ES-PGS_1000GenomesEur_corrected.csv')
pgs_corrected_piw <- read_csv('/Dedicated/jmichaelson-wdata/lcasten/sli_wgs/ancient_DNA/data/female_pelvic_inlet_width_ES-PGS_1000GenomesEur_corrected.csv') %>% 
  rename(pgs_pc_corrected_female_pelvic_inlet_width = pgs_pc_corrected) %>% 
  select(-PGS)%>% 
  filter(pgs_name == 'Base_1') %>% 
  select(-pgs_name)
pgs_corrected_pw <- read_csv('/Dedicated/jmichaelson-wdata/lcasten/sli_wgs/ancient_DNA/data/female_waist_circumference_bmi_adj_ES-PGS_1000GenomesEur_corrected.csv') %>% 
  rename(pgs_pc_corrected_female_waist_circ_BMI_adj = pgs_pc_corrected) %>% 
  select(-PGS)%>% 
  filter(pgs_name == 'Base_1') %>% 
  select(-pgs_name)

pgs_corrected_pwgt <- read_csv('/Dedicated/jmichaelson-wdata/lcasten/sli_wgs/ancient_DNA/data/placental_weight_ES-PGS_1000GenomesEur_corrected.csv') %>% 
  rename(pgs_pc_corrected_placental_weight = pgs_pc_corrected) %>% 
  select(-PGS)%>% 
  filter(pgs_name == 'Base_1') %>% 
  select(-pgs_name)
pgs_corrected_bw <- read_csv('/Dedicated/jmichaelson-wdata/lcasten/sli_wgs/ancient_DNA/data/fetal_birth_weight_ES-PGS_1000GenomesEur_corrected.csv') %>% 
  rename(pgs_pc_corrected_birth_weight = pgs_pc_corrected) %>% 
  select(-PGS)%>% 
  filter(pgs_name == 'Base_1') %>% 
  select(-pgs_name)
pgs_corrected_bhc <- read_csv('/Dedicated/jmichaelson-wdata/lcasten/sli_wgs/ancient_DNA/data/birth_head_circumference_ES-PGS_1000GenomesEur_corrected.csv') %>% 
  rename(pgs_pc_corrected_birth_head_circ = pgs_pc_corrected) %>% 
  select(-PGS)%>% 
  filter(pgs_name == 'Base_1') %>% 
  select(-pgs_name)
pgs_corrected_ihc <- read_csv('/Dedicated/jmichaelson-wdata/lcasten/sli_wgs/ancient_DNA/data/infant_head_circumference_ES-PGS_1000GenomesEur_corrected.csv') %>% 
  rename(pgs_pc_corrected_infant_head_circ = pgs_pc_corrected) %>% 
  select(-PGS) %>% 
  filter(pgs_name == 'Base_1') %>% 
  select(-pgs_name)
pgs_corrected_gd <- read_csv('/Dedicated/jmichaelson-wdata/lcasten/sli_wgs/ancient_DNA/data/gestational_duration_ES-PGS_1000GenomesEur_corrected.csv') %>% 
  rename(pgs_pc_corrected_gestation_duration = pgs_pc_corrected) %>% 
  select(-PGS) %>% 
  filter(pgs_name == 'Base_1') %>% 
  select(-pgs_name)
pgs_corrected_icv <- read_csv('/Dedicated/jmichaelson-wdata/lcasten/sli_wgs/ancient_DNA/data/intracranial_volume_adult_meta_ES-PGS_1000GenomesEur_corrected.csv') %>% 
  rename(pgs_pc_corrected_intracranial_vol = pgs_pc_corrected) %>% 
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
  rename(pgs_pc_corrected_phoneme_awareness = pgs_pc_corrected) %>% 
  select(-PGS) %>% 
  filter(pgs_name == 'Base_1') %>% 
  select(-pgs_name)
pgs_corrected_sp <- read_csv('/Dedicated/jmichaelson-wdata/lcasten/sli_wgs/ancient_DNA/data/genlang_SP_ES-PGS_1000GenomesEur_corrected.csv') %>% 
  rename(pgs_pc_corrected_spelling = pgs_pc_corrected) %>% 
  select(-PGS) %>% 
  filter(pgs_name == 'Base_1') %>% 
  select(-pgs_name)
pgs_corrected_wr <- read_csv('/Dedicated/jmichaelson-wdata/lcasten/sli_wgs/ancient_DNA/data/genlang_WR_ES-PGS_1000GenomesEur_corrected.csv') %>% 
  rename(pgs_pc_corrected_word_reading = pgs_pc_corrected) %>% 
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

pgs_all <- pgs_corrected %>% 
  filter(IID %in% soi_qc) %>% 
  filter(pgs_name %in% c('HAQER_1', 'complement_HAQER_1')) %>% 
  select(IID, pgs_name, pgs_pc_corrected) %>% 
  mutate(pgs_name = str_remove_all(pgs_name, pattern = '_1$')) %>% 
  mutate(pgs_name = str_c('cp_pgs.', pgs_name)) %>%
  pivot_wider(id_cols = IID, names_from = pgs_name, values_from = pgs_pc_corrected) %>% 
  rename(cp_pgs.background_HAQER = cp_pgs.complement_HAQER) %>%
  inner_join(pgs_corrected_piw) %>% 
  inner_join(pgs_corrected_pw) %>% 
  inner_join(pgs_corrected_pwgt) %>% 
  inner_join(pgs_corrected_bw) %>% 
  inner_join(pgs_corrected_bhc) %>% 
  inner_join(pgs_corrected_ihc) %>% 
  inner_join(pgs_corrected_gd) %>%
  inner_join(pgs_corrected_icv) %>%
  inner_join(pgs_corrected_nread) %>%
  inner_join(pgs_corrected_nrep) %>%
  inner_join(pgs_corrected_nviq) %>%
  inner_join(pgs_corrected_pa) %>%
  inner_join(pgs_corrected_wr) %>%
  inner_join(pgs_corrected_sp) %>%
  inner_join(pgs_other)

names(pgs_all)[-c(1:3)] <- str_remove_all(names(pgs_all)[-c(1:3)], pattern = 'pgs_pc_corrected_')
names(pgs_all)[-c(1:3)] <-str_c(names(pgs_all)[-c(1:3)] , '_pgs.genome_wide')

smeta2 %>% 
  inner_join(pca) %>%
  inner_join(pgs_all) %>% 
  filter(str_detect(IID, 'REF')) %>% 
  select(IID, matches('pgs')) %>% 
  pivot_longer(cols = -c(IID)) %>% 
  arrange(desc(abs(value))) %>% 
  filter(IID == 'Href.REF') %>% 
  select(-1) %>%
  as.data.frame()

tmp <- smeta2 %>% 
  mutate(IID = ifelse(IID == 'Chagyrskaya', 'Chagyrskaya-Phalanx', IID)) %>%
  inner_join(pca) %>%
  inner_join(pgs_all) %>% 
  filter(str_detect(IID, 'REF', negate = TRUE))

tmp %>% 
    arrange(desc(sample_age_years_before_1950))
tmp %>% 
    write_csv('/wdata/lcasten/sli_wgs/ancient_DNA/AADR_neanderthal_PGS_gathered.csv')

##########################
## read in imputed ES-PGS
imp <- read_csv('/wdata/lcasten/sli_wgs/ancient_DNA/data/cogPerf_ES-PGS_1000GenomesEur_corrected.imputed.csv') %>% 
  filter(! IID %in% c(smeta_nean$IID, 'Chagyrskaya-Phalanx')) %>%
  filter(pgs_name %in% c('HAQER_1', 'complement_HAQER_1')) %>% 
  mutate(pgs_name = str_remove_all(pgs_name, '_1$'),
         pgs_name = str_replace_all(pgs_name, pattern = 'complement_', replacement = 'background_'),
         pgs_name = str_c('cp_pgs.imputed.', pgs_name)) %>%
  pivot_wider(id_cols = IID, names_from = 'pgs_name', values_from = 'pgs_pc_corrected')
tmp2 <- imp %>% 
  inner_join(tmp)
hist(tmp2$cp_pgs.imputed.HAQER)
hist(tmp2$cp_pgs.imputed.background_HAQER)

cor.test(tmp2$cp_pgs.imputed.HAQER, tmp2$cp_pgs.HAQER)
cor.test(tmp2$cp_pgs.imputed.HAQER, log10(tmp2$sample_age_years_before_1950))

tmp2 %>% 
  select(-matches('pgs.genome_wide|pgs.HAQ|pgs.back|pc')) %>% 
  relocate(cp_pgs.imputed.HAQER, cp_pgs.imputed.background_HAQER, .after = age_at_death) %>% 
  write_csv('manuscript/supplemental_materials/AADR_imputed_data_no_neanderthals.csv')
  