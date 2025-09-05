library(tidyverse)

message('Running ES-PGS analysis now...')

####################################
## read in ES-PGS and phenotype data
####################################
dat <- read_csv('ES-PGS_example/example_data/gathered_es-pgs_corrected.cogPerf.csv')

fc = read_csv('ES-PGS_example/example_data/pheno_factors_resid.csv')

####################
## reformat data
####################
pgs_base_wide <- dat %>% 
  filter(str_detect(IID, pattern = 'sample')) %>%
  mutate(pgs = str_split(pgs_name, pattern = '[.]', simplify = TRUE)[,1]) %>%
  select(IID, pgs, pgs_genome_wide) %>% 
  distinct()

pgs_evo_wide <- dat %>% 
  filter(str_detect(IID, pattern = 'sample')) %>%
  select(-pgs_genome_wide) %>%
  select(-cohort) %>%
  mutate(pgs = str_split(pgs_name, pattern = '[.]', simplify = TRUE)[,1],
         gs = str_split(pgs_name, pattern = '[.]', simplify = TRUE)[,2],
         anno = str_remove_all(gs, pattern = '_1$')) %>%
  select(IID, pgs, anno, matches('pgs_pc_corrected')) %>%
  distinct() %>%
  pivot_wider(id_cols = -matches('anno|pgs_pathway'), 
              names_from = anno, 
              values_from = matches('pgs_pc_corrected'))
names(pgs_evo_wide) <- str_remove_all(names(pgs_evo_wide), pattern = '_pc_corrected')

pgs_wide <- pgs_base_wide %>% 
  inner_join(pgs_evo_wide)

## merge phenotype data with ES-PGS
tmp <- fc %>% 
  pivot_longer(cols = -c(1)) %>% 
  rename(IID = sample) %>%
  inner_join(pgs_wide)

##########################################
## ES-PGS analysis (reduced vs full model)
##########################################
## this code chunk will loop over each phenotype and ES-PGS annotation
## cols of interest 
coi <- names(pgs_wide)[str_detect(names(pgs_wide), pattern = 'pgs_') & str_detect(names(pgs_wide), 'complement|pgs_genome_wide|matched', negate = TRUE)]
## initilialize variables
iter = 0
res_list = list()

## adjust phenotype info to match your data!!
for(c in coi){
  cat(sprintf('\n\n'))
  message(c)
  for(i in 1:7){
      message('ES-PGS for phenotype: Factor ', i)
      iter = iter + 1
      tmp2 <- tmp %>% 
        filter(name == str_c('Factor', i)) %>% ## you'll need to change to match your phenotype names
        select(IID, value, matches(str_c(str_remove_all(c, pattern = 'pgs_'), '$')))
      bdat <- tmp2 %>% 
        select(IID, value, matches('complement|matched'))
      ## run stats
      baseline <- lm(value ~ ., data = bdat[,-1])
      baseline_plus_anno <- lm(value ~ ., data = tmp2[,-1])
      baseline_rsq = summary(baseline)$r.squared
      baseline_plus_anno_rsq = summary(baseline_plus_anno)$r.squared
      baseline_complement_coef <- broom::tidy(baseline_plus_anno) %>% 
        filter(term != '(Intercept)' & str_detect(term, 'complement'))
      baseline_matched_coef <- broom::tidy(baseline_plus_anno) %>% 
        filter(term != '(Intercept)' & str_detect(term, 'matched'))
      baseline_plus_anno_coef <- broom::tidy(baseline_plus_anno) %>% 
        filter(term != '(Intercept)' & str_detect(term, 'complement|matched', negate = TRUE))
      ## compare reduced vs full models and gather stats
      res <- broom::tidy(anova(baseline, baseline_plus_anno)) %>% 
        drop_na() %>% 
        mutate(factor = str_c('F', i)) %>% 
        mutate(evo = str_remove_all(c, 'pgs_')) %>% 
        relocate(factor, evo) %>% 
        mutate(baseline_rsq = baseline_rsq,
               baseline_plus_anno_rsq = baseline_plus_anno_rsq,
               annotation_beta = baseline_plus_anno_coef$estimate,
               annotation_std_err = baseline_plus_anno_coef$std.error,
               annotation_pval = baseline_plus_anno_coef$p.value,
               matched_beta = baseline_matched_coef$estimate,
               matched_std_err = baseline_matched_coef$std.error,
               matched_pval = baseline_matched_coef$p.value,
               background_beta = baseline_complement_coef$estimate,
               background_std_err = baseline_complement_coef$std.error,
               background_pval = baseline_complement_coef$p.value,
               )
      ## save results to list
      res_list[[iter]] <- res
      message('model comparison p = ', formatC(res$p.value, digits = 2), ', ES-PGS term beta = ', round(res$annotation_beta, digits = 3))
  }
}

## get SNP info for each ES-PGS
snp_counts <- read_table('ES-PGS_example/example_data/ES-PGS_raw.prsice') %>% 
  filter(Threshold == 1) %>% 
  select(anno = Set, Num_SNP)  %>% 
  mutate(model = str_remove_all(anno, pattern = 'complement_|matched_')) %>% 
  mutate(complement = case_when(str_detect(anno, pattern = 'complement_') ~ 'complement_PGS_n_snp', 
                                str_detect(anno, pattern = 'matched_') ~ 'matched_PGS_n_snp',
                                TRUE ~ 'annotation_PGS_n_snp')) %>%
  pivot_wider(id_cols = model, names_from = complement, values_from = Num_SNP) %>% 
  mutate(model = str_replace_all(model, pattern = '-', replacement = '_'))

##
espgs_res <- bind_rows(res_list) %>%
  arrange(p.value) %>% 
  select(-term) %>% 
  rename(model = evo, p.value_model_comparison = p.value) %>%
  inner_join(snp_counts) %>% 
  mutate(model = factor(model, levels = c('HAQERs', 'HARs'))) %>%
  drop_na(model) %>%
  arrange(factor, model) %>%
  relocate(complement_PGS_n_snp, .before = annotation_PGS_n_snp) %>%
  mutate(rsq_gain_per_1000indSNP = (baseline_plus_anno_rsq - baseline_rsq) / (annotation_PGS_n_snp / 1000)) %>% 
  rename(phenotype = factor, ES_PGS_model = model, reduced_model_rsq = baseline_rsq, full_model_rsq = baseline_plus_anno_rsq, background_PGS_n_snp = complement_PGS_n_snp) %>%
  mutate(relative_rsq_gain_per_1000indSNP = rsq_gain_per_1000indSNP / reduced_model_rsq)

names(espgs_res) <- str_replace_all(names(espgs_res), pattern = 'annotation_', replacement = 'ES_PGS_')

espgs_res %>%
  write_csv('ES-PGS_example/example_data/ES-PGS_lm_comparison_results.csv')
 
message('Done, saving ES-PGS model results to: "ES-PGS_example/example_data/ES-PGS_lm_comparison_results.csv')