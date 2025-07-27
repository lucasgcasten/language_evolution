library(tidyverse)


##############################################
##
dat <- read_csv('/Dedicated/jmichaelson-wdata/lcasten/sli_wgs/prs/pathway_prs/gathered_pgs_pc_corrected_long_full_data.cogPerf.complement_rand_controls.csv')

fc = read_csv('/wdata/common/SLI_WGS/public/phenotype/factors/pheno_factors_resid.csv')

unique(dat$pgs_name)
pgs_base_wide <- dat %>% 
  filter(str_detect(IID, pattern = 'sample')) %>%
  filter(str_detect(pgs_name, pattern = 'human_evo')) %>% 
  filter(str_detect(pgs_name, pattern = '_1$')) %>%
  mutate(pgs = str_split(pgs_name, pattern = '[.]', simplify = TRUE)[,1]) %>%
  select(IID, pgs, pgs_genome_wide_baseline) %>% 
  distinct()

unique(dat$pgs_name)

pgs_evo_wide <- dat %>% 
  filter(str_detect(IID, pattern = 'sample')) %>%
  filter(str_detect(pgs_name, pattern = 'human_evo|chimp')) %>% 
  filter(str_detect(pgs_name, pattern = '_1$')) %>%
  select(-matches('pgs_raw')) %>%
  select(-pgs_genome_wide_baseline) %>%
  select(-cohort) %>%
  mutate(pgs_name = str_replace_all(pgs_name, pattern = 'human_specific_evolution_brain_expression.', replacement = 'human_specific_evolution_brain_expression_')) %>%
  mutate(pgs = str_split(pgs_name, pattern = '[.]', simplify = TRUE)[,1],
         gs = str_split(pgs_name, pattern = '[.]', simplify = TRUE)[,2],
         anno_tmp = str_split(pgs_name, pattern = '[.]', simplify = TRUE)[,3],
         anno = str_remove_all(anno_tmp, pattern = '_5e_08$|_0.0005$|_0.05$|_0.2$|_1$|_0$'),
         thr = str_remove_all(anno_tmp, pattern = anno)) %>% # relocate(anno, thr, gs)
  select(IID, pgs, anno, matches('pgs_pc_corrected')) %>%
  distinct() %>%
  pivot_wider(id_cols = -matches('anno|pgs_pathway'), names_from = anno, values_from = matches('pgs_pc_corrected')) # %>% names()
names(pgs_evo_wide)
pgs_wide <- pgs_base_wide %>% 
  inner_join(pgs_evo_wide)

pgs_wide2 <- pgs_wide %>% 
  select(-c(pgs, pgs_genome_wide_baseline))
names(pgs_wide2) <- str_replace_all(names(pgs_wide2), pattern = 'pgs_pc_corrected_', replacement = 'cp_pgs_')
pgs_wide2 %>% 
  select(-matches('consVert|consMammal|consPrimates_65mya|exons|LinAR_Human|diverged_ExN|diverged_InN|diverged_NonN|introgressed')) %>% 
  relocate(cp_pgs_NeanderthalSelectiveSweep, .after = cp_pgs_HAR) %>% 
  relocate(cp_pgs_complement_NeanderthalSelectiveSweep, .after = cp_pgs_complement_HAR) %>% 
  write_csv('/wdata/lcasten/sli_wgs/prs/pathway_prs/results/es_pgs_data.csv')
names(pgs_wide2)

pgs_long <- pgs_wide %>% 
  pivot_longer(cols = -c(1:2), names_to = 'pgs_nm', values_to = 'pgs_val') %>% 
  mutate(pgs_nm = str_remove_all(pgs_nm, pattern = 'pgs_pc_corrected_')) %>% 
  distinct()

pgs_complement_long <- pgs_long %>% 
  filter(str_detect(pgs_nm, pattern = 'complement')) %>% 
  rename(complement_pgs_val = pgs_val) %>% 
  mutate(pgs_nm = str_remove_all(pgs_nm, pattern = 'complement_')) %>% 
  select(-pgs)

pgs_matched_long <- pgs_long %>% 
  filter(str_detect(pgs_nm, pattern = 'random_matched_control')) %>% 
  rename(matched_pgs_val = pgs_val) %>% 
  mutate(pgs_nm = str_remove_all(pgs_nm, pattern = 'random_matched_control_regions_')) %>% 
  select(-pgs)

unique(pgs_long$pgs_nm)

###########################
fc %>% 
  pivot_longer(cols = -c(1)) %>%
  rename(IID = sample) %>%
  inner_join(pgs_base_wide) %>%
  group_by(name) %>% 
  do(res = broom::tidy(cor.test(.$value, .$pgs_genome_wide_baseline, method = 'p'))) %>%
  unnest(res) %>% 
  arrange(p.value) %>% 
  mutate(lab = str_c('Pearson r = ', round(estimate, digits = 2), ', p-val = ', formatC(p.value, digits = 2)))


res <- fc %>% 
  pivot_longer(cols = -c(1)) %>%
  rename(IID = sample) %>%
  inner_join(pgs_long) %>%
  group_by(name, pgs, pgs_nm) %>% 
  do(res = broom::tidy(cor.test(.$value, .$pgs_val, method = 'p'))) %>%
  unnest(res) %>% 
  arrange(p.value) %>% 
  mutate(lab = str_c('Pearson r = ', round(estimate, digits = 2), ', p-val = ', formatC(p.value, digits = 2)))

unique(res$pgs_nm)

tmpp <- fc %>% 
  pivot_longer(cols = -c(1)) %>%
  rename(IID = sample) %>%
  # inner_join(select(pgs_wide, IID, pgs_genome_wide_baseline)) %>% 
  inner_join(pgs_complement_long) %>%
  inner_join(pgs_matched_long) %>%
  filter(pgs_nm == 'HAQER_v2') %>% 
  inner_join(pgs_long)

summary(lm(value ~ complement_pgs_val + matched_pgs_val + pgs_val, data = tmpp[tmpp$name == 'Factor1',]))

res2 <- fc %>% 
  pivot_longer(cols = -c(1)) %>%
  rename(IID = sample) %>%
  # inner_join(select(pgs_wide, IID, pgs_genome_wide_baseline)) %>% 
  inner_join(pgs_complement_long) %>%
  inner_join(pgs_matched_long) %>%
  group_by(name, pgs_nm) %>%
  mutate(resid_value = resid(lm(value ~ complement_pgs_val + matched_pgs_val)),
         resid_value = scale(resid_value)[,1]) %>% 
  ungroup() %>%
  inner_join(pgs_long) %>%
  filter(pgs_nm != 'pgs_genome_wide_baseline') %>%
  group_by(name, pgs, pgs_nm) %>% 
  do(res = broom::tidy(cor.test(.$resid_value, .$pgs_val, method = 'p'))) %>%
  unnest(res) %>% 
  arrange(p.value) %>% 
  mutate(lab = str_c('Pearson r = ', round(estimate, digits = 2), ', p-val = ', formatC(p.value, digits = 2)))

res2 %>%
  # filter(p.value < 0.05) %>% 
  filter(str_detect(pgs_nm, 'intro|Sweep')) %>%
  select(1,3,estimate, statistic) %>%
  filter(str_detect(name, '1|2|3'))



foxp2_targets = c('L4_6_RORB_2', 'L5_6_THEMIS_1')
ct_res <- res2 %>% 
  select(1, 3:6) %>% 
  select(-statistic) %>% 
  filter(str_detect(pgs_nm, pattern = 'human_specific')) %>% 
  mutate(pgs_nm = str_remove_all(pgs_nm, pattern = 'human_specific_evolution_brain_expression_')) %>% 
  filter(pgs_nm == 'L3_5_RORB_1') %>% 
  filter(str_detect(name, '1|2|3'))

ct_atac_res <- res2 %>% 
  select(1, 3:6) %>% 
  select(-statistic) %>% 
  filter(str_detect(pgs_nm, pattern = 'human_specific_brain_CRE')) %>% 
  mutate(pgs_nm = str_remove_all(pgs_nm, pattern = 'human_specific_brain_CRE_')) %>% 
  filter(pgs_nm == 'L3_5_RORB_1') %>% 
  filter(str_detect(name, '1|2|3'))

res %>% 
  select(name, pgs_nm, estimate, p.value) %>% 
  filter(str_detect(pgs_nm, 'HAQER'))
res %>% 
  select(name, pgs_nm, estimate, p.value) %>% 
  filter(str_detect(pgs_nm, 'chimp'))
res %>% 
  select(name, pgs_nm, estimate, p.value) %>% 
  filter(str_detect(pgs_nm, 'LinAR')) %>%
  filter(str_detect(pgs_nm, 'complement', negate = TRUE)) %>% 
  filter(str_detect(name, '1|2'))


res %>% 
  select(name, pgs_nm, estimate, p.value) %>% 
  filter(p.value < 0.05) %>% 
  filter(str_detect(pgs_nm, pattern = 'complement', negate = TRUE)) %>% 
  filter(str_detect(name, '1|2'))

res %>% 
  select(name, pgs_nm, estimate, p.value) %>% 
  # filter(p.value < 0.05) %>% 
  filter(str_detect(pgs_nm, pattern = 'single|UCE', negate = FALSE)) %>% filter(str_detect(pgs_nm, pattern = 'complement_', negate = TRUE))
unique(res$pgs_nm)


p_dat <- fc %>% 
  rename(IID = sample) %>% 
  pivot_longer(cols = -c(1)) %>%
  inner_join(res) %>%
  mutate(name = str_remove_all(name, pattern = 'actor')) %>%
  inner_join(pgs_long) %>% 
  filter(name %in% c('F1', 'F2')) %>%
  filter(pgs_nm %in% c('pgs_genome_wide_baseline', 'consPrimates_65mya', 'consPrimates_UCE', 'LinAR_Simiformes', 'LinAR_Catarrhini', 'LinAR_Hominoidea', 'LinAR_Hominidae', 'LinAR_Homininae', 'LinAR_Human', 'human_chimp_div_DMG', 'HAQER', 'HAR', 'NeanderthalSelectiveSweep', 'human_singleton_density_score_top5pct')) %>% 
  mutate(pgs_nm = str_replace_all(pgs_nm, pattern = 'consPrimates_65mya', replacement = 'primate conserved loci (65 Mya)'),
         pgs_nm = str_replace_all(pgs_nm, pattern = 'consPrimates_UCE', replacement = 'primate UCEs (65 Mya)'),
         pgs_nm = str_replace_all(pgs_nm, pattern = 'LinAR_Simiformes', replacement = 'Simian accelerated regions  (> 30 Mya)'),
         pgs_nm = str_replace_all(pgs_nm, pattern = 'LinAR_Catarrhini', replacement = 'Old world monkey accelerated regions (30 Mya)'),
         pgs_nm = str_replace_all(pgs_nm, pattern = 'LinAR_Hominoidea', replacement = 'Ape accelerated regions (20 Mya)'),
         pgs_nm = str_replace_all(pgs_nm, pattern = 'LinAR_Hominidae', replacement = 'Great ape accelerated regions (15 Mya)'),
         pgs_nm = str_replace_all(pgs_nm, pattern = 'LinAR_Homininae', replacement = 'Hominid accelerated regions (9 Mya)'),
         pgs_nm = str_replace_all(pgs_nm, pattern = 'LinAR_Human', replacement = 'Human lineage accelerated regions (6 Mya)'),
         pgs_nm = str_replace_all(pgs_nm, pattern = 'human_chimp_div_DMG', replacement = 'human-chimp divergent genes (6 Mya)'),
         pgs_nm = str_replace_all(pgs_nm, pattern = 'NeanderthalSelectiveSweep', replacement = 'Neanderthal selective sweep loci (50 Kya)'),
         pgs_nm = str_replace_all(pgs_nm, pattern = 'HAR', replacement = 'HARs (600 Kya)'),
         pgs_nm = str_replace_all(pgs_nm, pattern = 'HAQER', replacement = 'HAQERs (600 Kya)'),
         pgs_nm = str_replace_all(pgs_nm, pattern = 'human_singleton_density_score_top5pct', replacement = 'recent selection loci (2-3 Kya)')
        ) %>%
  mutate(type = case_when(pgs_nm == 'pgs_genome_wide_baseline' ~ str_c('Genome-wide PGS'),
                          pgs_nm != 'pgs_genome_wide_baseline' ~ str_c('PGS in ', pgs_nm)),
         type = factor(type, levels = c('Genome-wide PGS', 'PGS in primate conserved loci (65 Mya)', 'PGS in primate UCEs (65 Mya)', 'PGS in Simian accelerated regions  (> 30 Mya)', 'PGS in Old world monkey accelerated regions (30 Mya)', 'PGS in Ape accelerated regions (20 Mya)', 'PGS in Great ape accelerated regions (15 Mya)', 'PGS in Hominid accelerated regions (9 Mya)', 'PGS in human-chimp divergent genes (6 Mya)', 'PGS in Human lineage accelerated regions (6 Mya)', 'PGS in HAQERs (600 Kya)', 'PGS in HARs (600 Kya)', 'PGS in Neanderthal selective sweep loci (50 Kya)', 'PGS in recent selection loci (2-3 Kya)')))
unique(p_dat$type)
p_dat %>% 
  filter(str_detect(type, 'Simian'))
p_dat %>% 
  filter(is.na(type))

p_dat %>% 
    write_csv('/wdata/lcasten/sli_wgs/prs/pathway_prs/results/fig_2_pathway_pgs_data.complement.csv')

###############
p_dat2 <- fc %>% 
  rename(IID = sample) %>% 
  pivot_longer(cols = -c(1)) %>%
  inner_join(pgs_complement_long) %>%
  group_by(name, pgs_nm) %>%
  mutate(resid_value = resid(lm(value ~ complement_pgs_val)),
         resid_value = scale(resid_value)[,1]) %>% 
  ungroup() %>%
  inner_join(pgs_long) %>%
  filter(pgs_nm != 'pgs_genome_wide_baseline') %>%
  inner_join(res2) %>%
  mutate(name = str_remove_all(name, pattern = 'actor')) %>%
  inner_join(pgs_long) %>% 
  filter(name %in% c('F1', 'F2')) %>%
  filter(pgs_nm %in% c('pgs_genome_wide_baseline', 'consPrimates_65mya', 'consPrimates_UCE', 'LinAR_Simiformes', 'LinAR_Catarrhini', 'LinAR_Hominoidea', 'LinAR_Hominidae', 'LinAR_Homininae', 'LinAR_Human', 'human_chimp_div_DMG', 'HAQER', 'HAR', 'NeanderthalSelectiveSweep', 'human_singleton_density_score_top5pct')) %>% 
  mutate(pgs_nm = str_replace_all(pgs_nm, pattern = 'consPrimates_65mya', replacement = 'primate conserved loci (65 Mya)'),
         pgs_nm = str_replace_all(pgs_nm, pattern = 'consPrimates_UCE', replacement = 'primate UCEs (65 Mya)'),
         pgs_nm = str_replace_all(pgs_nm, pattern = 'LinAR_Simiformes', replacement = 'Simian accelerated regions  (> 30 Mya)'),
         pgs_nm = str_replace_all(pgs_nm, pattern = 'LinAR_Catarrhini', replacement = 'Old world monkey accelerated regions (30 Mya)'),
         pgs_nm = str_replace_all(pgs_nm, pattern = 'LinAR_Hominoidea', replacement = 'Ape accelerated regions (20 Mya)'),
         pgs_nm = str_replace_all(pgs_nm, pattern = 'LinAR_Hominidae', replacement = 'Great ape accelerated regions (15 Mya)'),
         pgs_nm = str_replace_all(pgs_nm, pattern = 'LinAR_Homininae', replacement = 'Hominid accelerated regions (9 Mya)'),
         pgs_nm = str_replace_all(pgs_nm, pattern = 'LinAR_Human', replacement = 'Human lineage accelerated regions (6 Mya)'),
         pgs_nm = str_replace_all(pgs_nm, pattern = 'human_chimp_div_DMG', replacement = 'human-chimp divergent genes (6 Mya)'),
         pgs_nm = str_replace_all(pgs_nm, pattern = 'NeanderthalSelectiveSweep', replacement = 'Neanderthal selective sweep loci (50 Kya)'),
         pgs_nm = str_replace_all(pgs_nm, pattern = 'HAR', replacement = 'HARs (600 Kya)'),
         pgs_nm = str_replace_all(pgs_nm, pattern = 'HAQER', replacement = 'HAQERs (600 Kya)'),
         pgs_nm = str_replace_all(pgs_nm, pattern = 'human_singleton_density_score_top5pct', replacement = 'recent selection loci (2-3 Kya)')
        ) %>%
  mutate(type = case_when(pgs_nm == 'pgs_genome_wide_baseline' ~ str_c('Genome-wide PGS'),
                          pgs_nm != 'pgs_genome_wide_baseline' ~ str_c('PGS in ', pgs_nm)),
         type = factor(type, levels = c('Genome-wide PGS', 'PGS in primate conserved loci (65 Mya)', 'PGS in primate UCEs (65 Mya)', 'PGS in Simian accelerated regions  (> 30 Mya)', 'PGS in Old world monkey accelerated regions (30 Mya)', 'PGS in Ape accelerated regions (20 Mya)', 'PGS in Great ape accelerated regions (15 Mya)', 'PGS in Hominid accelerated regions (9 Mya)', 'PGS in human-chimp divergent genes (6 Mya)', 'PGS in Human lineage accelerated regions (6 Mya)', 'PGS in HAQERs (600 Kya)', 'PGS in HARs (600 Kya)', 'PGS in Neanderthal selective sweep loci (50 Kya)', 'PGS in recent selection loci (2-3 Kya)')))
unique(p_dat2$type)

unique(p_dat2$pgs_nm)
p_dat2 %>% 
    write_csv('/wdata/lcasten/sli_wgs/prs/pathway_prs/results/fig_2_pathway_pgs_data_corrected.complement.csv')



p <- fc %>% 
  rename(IID = sample) %>% 
  pivot_longer(cols = -c(1)) %>%
  inner_join(pgs_complement_long) %>%
  group_by(name, pgs_nm) %>%
  mutate(resid_value = resid(lm(value ~ complement_pgs_val)),
         resid_value = scale(resid_value)[,1]) %>% 
  ungroup() %>%
  inner_join(pgs_long) %>%
  filter(pgs_nm != 'pgs_genome_wide_baseline') %>%
  inner_join(res2) %>%
  mutate(name = str_remove_all(name, pattern = 'actor')) %>%
  inner_join(pgs_long) %>% 
  filter(name %in% c('F1')) %>%
  filter(str_detect(pgs_nm, 'HAQER')) %>% 
  inner_join(rename(snp_counts, pgs_nm = model)) %>% 
  mutate(pgs_nm = str_c(pgs_nm, '\n(N snp = ', annotation_PGS_n_snp,')')) %>%
  ggplot(aes(x = resid_value, y = pgs_val)) +
  geom_point(size = 1.5) +
  geom_smooth(method = 'lm') +
  xlab('Cognitive performance ES-PGS') +
  ylab('Factor score adjusted\nfor background PGS') +
  facet_wrap(~ pgs_nm) +
  theme_classic() +
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 16),
        strip.text = element_text(size = 16)) +
  geom_text(aes(x = 1.8, y = -3.3, label = lab), size = 3.5, check_overlap = TRUE)
p %>% 
  ggsave(filename = '/wdata/lcasten/sli_wgs/paper_figures/subplots/HAQER_ES-PGS_size_cutoff_comparison.png', device = 'png', width = 15, height = 6, dpi = 300)

############
unique(p_dat3$pgs_nm)
p_dat3 <- fc %>% 
  rename(IID = sample) %>% 
  pivot_longer(cols = -c(1)) %>%
  inner_join(pgs_complement_long) %>%
  group_by(name, pgs_nm) %>%
  mutate(resid_value = resid(lm(value ~ complement_pgs_val)),
         resid_value = scale(resid_value)[,1]) %>% 
  ungroup() %>%
  inner_join(pgs_long) %>%
  filter(pgs_nm != 'pgs_genome_wide_baseline') %>%
  inner_join(res2) %>%
  mutate(name = str_remove_all(name, pattern = 'actor')) %>%
  filter(str_detect(pgs_nm, 'human_specific_')) %>%
  inner_join(pgs_long) %>% 
  filter(name %in% c('F1', 'F2'))
unique(p_dat3$pgs_nm)

p_dat3 %>% 
    write_csv('/wdata/lcasten/sli_wgs/prs/pathway_prs/results/fig_2_pathway_pgs_data_corrected_cell_type_specific.complement.csv')




####
p_dat2 %>% 
  filter(name == 'F1') %>%
  ggplot(aes(x = pgs_val, y = resid_value)) +
  geom_point(size = 2) +
  geom_smooth(method = 'lm') +
  facet_wrap(~ type, nrow = 1) +
  xlab('Cognitive performance polygenic score') +
  ylab('F1 score') +
  geom_text(aes(x = -2, y = 2.7, label = lab), size = 5, check_overlap = TRUE)

p_f1 <- p_dat %>% 
  filter(name == 'F1') %>%
  ggplot(aes(x = pgs_val, y = value)) +
  geom_point(size = 2) +
  geom_smooth(method = 'lm') +
  facet_wrap(~ type, nrow = 1) +
  xlab('Cognitive performance polygenic score') +
  ylab('F1 score') +
  geom_text(aes(x = -2, y = 2.7, label = lab), size = 5, check_overlap = TRUE)
p_f2 <- p_dat %>% 
  filter(name == 'F2') %>%
  ggplot(aes(x = pgs_val, y = value)) +
  geom_point(size = 2) +
  geom_smooth(method = 'lm') +
  facet_wrap(~ type, nrow = 1) +
  xlab('Cognitive performance polygenic score') +
  ylab('F2 score') +
  geom_text(aes(x = -2.5, y = 1.9, label = lab), size = 5, check_overlap = TRUE)

###################################


#############################
## do anova analysis for regressions w/ each of the annos
pgs_wide
coi <- names(pgs_wide)[str_detect(names(pgs_wide), pattern = 'pgs_pc_corrected') & str_detect(names(pgs_wide), 'complement', negate = TRUE)]

tmp <- fc %>% 
  pivot_longer(cols = -c(1)) %>% 
  rename(IID = sample) %>%
  inner_join(pgs_wide)

iter = 0
res_list = list()
for(c in coi){
  cat(sprintf('\n\n\n\n'))
  message(c)
  for(i in 1:7){
      message('Factor ', i)
      iter = iter + 1
      tmp2 <- tmp %>% 
        filter(name == str_c('Factor', i)) %>% 
        select(IID, value, matches(str_c(str_remove_all(c, pattern = 'pgs_pc_corrected_'), '$')))
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
snp_counts <- read_csv('/Dedicated/jmichaelson-wdata/lcasten/sli_wgs/prs/pathway_prs/prsice_info.complement_rand_controls.csv') %>% 
  filter(Threshold == 1) %>% 
  select(anno = Set, Num_SNP)  %>% 
  mutate(model = str_remove_all(anno, pattern = 'complement_')) %>% 
  mutate(complement = ifelse(str_detect(anno, pattern = 'complement_'), 'complement_PGS_n_snp', 'annotation_PGS_n_snp')) %>%
  pivot_wider(id_cols = model, names_from = complement, values_from = Num_SNP) %>% 
  mutate(model = str_replace_all(model, pattern = 'human_specific_evolution_brain_expression.', replacement = 'human_specific_evolution_brain_expression_'),
         model = str_replace_all(model, pattern = '-', replacement = '_'))
unique(snp_counts$model)
snp_counts %>% 
  as.data.frame()

##
unique(bind_rows(res_list)$c)
bind_rows(res_list) %>% filter(str_detect(factor, '1|2|3')) %>% filter(str_detect(c, 'human_specific')) %>% arrange(p.value) %>%  relocate(c, factor, annotation_beta, p.value) %>% mutate(c = str_remove_all(c, 'human_specific_evolution_brain_expression_|human_specific_brain_')) %>% mutate(factor = str_remove_all(factor, 'Factor'))
bind_rows(res_list) %>% filter(str_detect(factor, '1|2|3')) %>% filter(str_detect(c, 'human_specific')) %>% arrange(p.value) %>%  relocate(c, factor, annotation_beta, p.value) %>% mutate(c = str_remove_all(c, 'human_specific_evolution_brain_expression_|human_specific_brain_')) %>% mutate(factor = str_remove_all(factor, 'Factor')) %>% filter(str_detect(c, 'FOXP2'))
bind_rows(res_list) %>% filter(str_detect(factor, '1|2|3')) %>% filter(str_detect(c, 'HAQER')) %>% arrange(p.value) %>%  relocate(c, factor, annotation_beta, p.value) %>% mutate(c = str_remove_all(c, 'human_specific_evolution_brain_expression_|human_specific_brain_')) %>% mutate(factor = str_remove_all(factor, 'Factor')) %>% select(1:4, matches('rsq')) %>% arrange(factor, p.value) %>% mutate(rsq_gain_pct = (baseline_plus_anno_rsq - baseline_rsq) / baseline_rsq) %>% select(-matches('baseline'))# %>% filter(str_detect(c, 'FOXP2'))


bind_rows(res_list) %>% 
  arrange(p.value) %>% 
  select(-term) %>% 
  rename(model = c, p.value_model_comparison = p.value) %>%
  inner_join(snp_counts) %>% 
  # filter(str_detect(model, 'human_')) %>% as.data.frame
  # mutate(model = factor(model, levels = c('pgs_genome_wide_baseline', 'consPrimates_65mya', 'consPrimates_UCE', 'LinAR_Simiformes', 'LinAR_Catarrhini', 'LinAR_Hominoidea', 'LinAR_Hominidae', 'LinAR_Homininae', 'LinAR_Human', 'human_chimp_div_DMG', 'HAQER', 'HAR', 'NeanderthalSelectiveSweep', 'human_singleton_density_score_top5pct'))) %>%
  drop_na(model) %>%
  arrange(factor, model) %>%
  relocate(complement_PGS_n_snp, .before = annotation_PGS_n_snp) %>% 
  # filter(p.value_model_comparison < 0.05) %>% 
  select(-c(3:6)) %>% 
  relocate(matches('annotation_')) %>% 
  filter(str_detect(factor, '1|2')) %>% 
  # filter(annotation_beta > 0) %>% 
  relocate(model, factor) %>% 
  arrange(p.value_model_comparison)

bind_rows(res_list) %>%
  arrange(p.value) %>% 
  select(-term) %>% 
  rename(model = c, p.value_model_comparison = p.value) %>%
  inner_join(snp_counts) %>% 
  mutate(model = factor(model, levels = c('pgs_genome_wide_baseline', 'consPrimates_65mya', 'consPrimates_UCE', 'LinAR_Simiformes', 'LinAR_Catarrhini', 'LinAR_Hominoidea', 'LinAR_Hominidae', 'LinAR_Homininae', 'LinAR_Human', 'human_chimp_div_DMG', 'HAQER', 'HAR', 'NeanderthalSelectiveSweep', 'human_singleton_density_score_top5pct'))) %>%
  drop_na(model) %>%
  arrange(factor, model) %>%
  relocate(complement_PGS_n_snp, .before = annotation_PGS_n_snp) %>%
  mutate(rsq_gain_per_1000indSNP = (baseline_plus_anno_rsq - baseline_rsq) / (annotation_PGS_n_snp / 1000)) %>% 
  write_csv('/wdata/lcasten/sli_wgs/prs/pathway_prs/results/evo_prs_lm_comparison_results.csv')

##
bind_rows(res_list) %>%
  arrange(p.value) %>% 
  select(-term) %>% 
  rename(model = c, p.value_model_comparison = p.value) %>%
  inner_join(snp_counts) %>% 
  # mutate(model = factor(model, levels = c('pgs_genome_wide_baseline', 'consPrimates_65mya', 'consPrimates_UCE', 'LinAR_Simiformes', 'LinAR_Catarrhini', 'LinAR_Hominoidea', 'LinAR_Hominidae', 'LinAR_Homininae', 'LinAR_Human', 'human_chimp_div_DMG', 'HAQER', 'HAR', 'NeanderthalSelectiveSweep', 'human_singleton_density_score_top5pct'))) %>%
  filter(str_detect(model, 'human_specific_evolution_brain_expression_')) %>% 
  mutate(model = str_remove_all(model, 'human_specific_evolution_brain_expression_')) %>%
  drop_na(model) %>%
  arrange(factor, model) %>%
  relocate(complement_PGS_n_snp, .before = annotation_PGS_n_snp) %>%
  mutate(rsq_gain_per_1000indSNP = (baseline_plus_anno_rsq - baseline_rsq) / (annotation_PGS_n_snp / 1000)) %>% 
  write_csv('/wdata/lcasten/sli_wgs/prs/pathway_prs/results/evo_prs_cell_type_lm_comparison_results.csv')

bind_rows(res_list) %>%
  arrange(p.value) %>% 
  select(-term) %>% 
  rename(model = c, p.value_model_comparison = p.value) %>%
  inner_join(snp_counts) %>% 
  # mutate(model = factor(model, levels = c('pgs_genome_wide_baseline', 'consPrimates_65mya', 'consPrimates_UCE', 'LinAR_Simiformes', 'LinAR_Catarrhini', 'LinAR_Hominoidea', 'LinAR_Hominidae', 'LinAR_Homininae', 'LinAR_Human', 'human_chimp_div_DMG', 'HAQER', 'HAR', 'NeanderthalSelectiveSweep', 'human_singleton_density_score_top5pct'))) %>%
  filter(str_detect(model, 'human_specific_brain_CRE')) %>% 
  mutate(model = str_remove_all(model, 'human_specific_brain_CRE_')) %>%
  drop_na(model) %>%
  arrange(factor, model) %>%
  relocate(complement_PGS_n_snp, .before = annotation_PGS_n_snp) %>%
  mutate(rsq_gain_per_1000indSNP = (baseline_plus_anno_rsq - baseline_rsq) / (annotation_PGS_n_snp / 1000)) %>% 
  write_csv('/wdata/lcasten/sli_wgs/prs/pathway_prs/results/evo_prs_cell_type_scATAC_lm_comparison_results.csv')

dat_p_ct <- fc %>% 
  pivot_longer(cols = -c(1)) %>%
  rename(IID = sample) %>%
  # inner_join(select(pgs_wide, IID, pgs_genome_wide_baseline)) %>% 
  inner_join(pgs_complement_long) %>%
  group_by(name, pgs_nm) %>%
  mutate(resid_value = resid(lm(value ~ complement_pgs_val)),
         resid_value = scale(resid_value)[,1]) %>% 
  ungroup() %>%
  inner_join(pgs_long) %>%
  filter(pgs_nm != 'pgs_genome_wide_baseline') %>% 
  filter(str_detect(pgs_nm, 'L3_5_RORB_1')) %>% 
  mutate(pgs_nm = str_remove_all(pgs_nm, 'human_specific_evolution_brain_expression_')) %>%
  inner_join(ct_res)
dat_p_ct %>% 
  mutate(lab = str_c('r = ', round(estimate, digits = 2), ', p-val = ', formatC(p.value, digits = 2))) %>%
  mutate(sig = ifelse(p.value < 0.05, TRUE, FALSE)) %>%
  ggplot(aes(x = resid_value, y = pgs_val)) +
  geom_point(aes(color = sig)) +
  geom_smooth(method = 'lm') +
  geom_text(aes(x = 0, y = 2.5, label = lab), check_overlap = TRUE, size = 5) +
  facet_wrap(~name) +
  xlab('Factor score (adjusted for complement PGS)') +
  ylab('L3-5 RORB-1 cog. ES-PGS') +
  theme_classic() +
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 16, face = 'bold'),
        strip.text = element_text(size = 16, face = 'bold'),
        legend.position = 'none') +
  scale_color_manual(values = c('grey85', 'black'))

##
bind_rows(res_list) %>% filter(str_detect(c, 'UCE'))

bind_rows(res_list) %>%
  arrange(p.value) %>% 
  select(-term) %>% 
  rename(model = c, p.value_model_comparison = p.value) %>%
  inner_join(snp_counts) %>% 
  mutate(model = factor(model, levels = c('pgs_genome_wide_baseline', 'consPrimates_65mya', 'consPrimates_UCE', 'LinAR_Simiformes', 'LinAR_Catarrhini', 'LinAR_Hominoidea', 'LinAR_Hominidae', 'LinAR_Homininae', 'LinAR_Human', 'human_chimp_div_DMG', 'HAQER', 'HAR', 'NeanderthalSelectiveSweep', 'human_singleton_density_score_top5pct'))) %>%
  drop_na(model) %>%
  arrange(factor, model) %>%
  relocate(complement_PGS_n_snp, .before = annotation_PGS_n_snp) %>%
  mutate(rsq_gain_per_1000indSNP = (baseline_plus_anno_rsq - baseline_rsq) / (annotation_PGS_n_snp / 1000))  %>% 
  arrange(desc(rsq_gain_per_1000indSNP)) %>% 
  select(-c(3:6)) %>% 
  filter(str_detect(factor, '1|2'))

## draft plot of Rsq gain
bind_rows(res_list) %>% 
  arrange(p.value) %>% 
  select(-term) %>% 
  rename(model = c, p.value_model_comparison = p.value) %>%
  inner_join(snp_counts) %>% 
  mutate(model = factor(model, levels = c('pgs_genome_wide_baseline', 'consPrimates_65mya', 'consPrimates_UCE', 'LinAR_Simiformes', 'LinAR_Catarrhini', 'LinAR_Hominoidea', 'LinAR_Hominidae', 'LinAR_Homininae', 'LinAR_Human', 'human_chimp_div_DMG', 'HAQER', 'HAR', 'NeanderthalSelectiveSweep', 'human_singleton_density_score_top5pct'))) %>%
  drop_na(model) %>%
  mutate(rsq_gain_per_1000indSNP = (baseline_plus_anno_rsq - baseline_rsq) / (annotation_PGS_n_snp / 1000)) %>% 
  mutate(model = str_replace_all(model, pattern = 'consPrimates_65mya', replacement = 'primate conserved loci (65 Mya)'),
         model = str_replace_all(model, pattern = 'consPrimates_UCE', replacement = 'primate UCEs (65 Mya)'),
         model = str_replace_all(model, pattern = 'LinAR_Simiformes', replacement = 'Simian accelerated regions  (> 30 Mya)'),
         model = str_replace_all(model, pattern = 'LinAR_Catarrhini', replacement = 'Old world monkey accelerated regions (30 Mya)'),
         model = str_replace_all(model, pattern = 'LinAR_Hominoidea', replacement = 'Ape accelerated regions (20 Mya)'),
         model = str_replace_all(model, pattern = 'LinAR_Hominidae', replacement = 'Great ape accelerated regions (15 Mya)'),
         model = str_replace_all(model, pattern = 'LinAR_Homininae', replacement = 'Hominid accelerated regions (9 Mya)'),
         model = str_replace_all(model, pattern = 'LinAR_Human', replacement = 'Human lineage accelerated regions (6 Mya)'),
         model = str_replace_all(model, pattern = 'human_chimp_div_DMG', replacement = 'human-chimp divergent genes (6 Mya)'),
         model = str_replace_all(model, pattern = 'NeanderthalSelectiveSweep', replacement = 'Neanderthal selective sweep loci (50 Kya)'),
         model = str_replace_all(model, pattern = 'HAR', replacement = 'HARs (600 Kya)'),
         model = str_replace_all(model, pattern = 'HAQER', replacement = 'HAQERs (600 Kya)'),
         model = str_replace_all(model, pattern = 'human_singleton_density_score_top5pct', replacement = 'recent selection loci (2-3 Kya)')
        ) %>%
  mutate(rsq_gain_per_1000indSNP = ifelse(annotation_beta < 0 | p.value_model_comparison > 0.05, 0, rsq_gain_per_1000indSNP)) %>%
  mutate(model = factor(model, levels = c('Genome-wide PGS', 'primate conserved loci (65 Mya)', 'primate UCEs (65 Mya)', 'Simian accelerated regions  (> 30 Mya)', 'Old world monkey accelerated regions (30 Mya)', 'Ape accelerated regions (20 Mya)', 'Great ape accelerated regions (15 Mya)', 'Hominid accelerated regions (9 Mya)', 'human-chimp divergent genes (6 Mya)', 'Human lineage accelerated regions (6 Mya)', 'HAQERs (600 Kya)', 'HARs (600 Kya)', 'Neanderthal selective sweep loci (50 Kya)', 'recent selection loci (2-3 Kya)'))) %>% 
  filter(model %in% c('primate UCEs (65 Mya)', 'Simian accelerated regions  (> 30 Mya)', 'Old world monkey accelerated regions (30 Mya)', 'Ape accelerated regions (20 Mya)', 'Great ape accelerated regions (15 Mya)', 'Hominid accelerated regions (9 Mya)', 'human-chimp divergent genes (6 Mya)', 'Human lineage accelerated regions (6 Mya)', 'HAQERs (600 Kya)', 'HARs (600 Kya)', 'Neanderthal selective sweep loci (50 Kya)', 'recent selection loci (2-3 Kya)')) %>%
  filter(factor %in% c('Factor1', 'Factor2')) %>% 
  mutate(factor = str_remove_all(factor, 'actor')) %>%
  drop_na() %>%
  ggplot(aes(x = model, y = rsq_gain_per_1000indSNP / baseline_rsq, color = factor, group = factor)) +
  geom_point(size = 5) +
  geom_line(size = 1.1) +
  xlab('ES-PGS annotation') +
  ylab('Relative R-squared improvement per\n1000 independent SNPs') +
  scale_y_continuous(labels = scales::percent) +
  geom_hline(yintercept = 0, color = 'red2', linetype = 'dashed', size = 1.075) +
  labs(color = 'Language factor') +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 15, face = 'bold'),
        legend.position = 'bottom')
# tab1 <- bind_rows(res_list) %>% 
#   arrange(p.value) %>% 
#   select(-term) %>% 
#   rename(model = c) %>% 
#   inner_join(snp_counts)
#   select(factor, model, statistic, p.value, baseline_rsq, baseline_plus_anno_rsq) %>% 
#   gridExtra::tableGrob(
#     rows = NULL,
#     cols = c('Factor', 'Annotation', 'Statistic', 'P', 'Baseline R sq.', 'Baseline + Annotation R sq.'), 
#     theme = gridExtra::ttheme_minimal(
#       core = list(
#         bg_params = list(fill = c("grey99", "grey96"), lwd = 1.5, col = "white"),
#         fg_params=list(hjust=0, x=0)),
#       colhead = list(fg_params=list(hjust=0, x=0))
#       )
#     )
# plot(tab1)
