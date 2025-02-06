library(tidyverse)


p_dat1 <- read_csv('/wdata/lcasten/sli_wgs/prs/pgs_correlation_heatmap_data.csv') %>% 
    filter(str_detect(pgs_name, pattern = 'sex', negate = TRUE)) %>% 
    mutate(fdr = p.adjust(p.value, method = 'fdr')) %>% 
    mutate(type = case_when(str_detect(pgs_name, pattern = 'empath') ~ 'Behavioral/SES',
                            TRUE ~ type))

# p_dat <- read_csv('/wdata/lcasten/sli_wgs/prs/pathway_prs/results/fig_2_pathway_pgs_data.csv')
raw_pgs <- read_csv('/wdata/lcasten/sli_wgs/prs/LDPred2-inf-v2/gathered_normed_pgs.long.csv')

# ph %>% 
#   pivot_longer(cols = matches('Factor')) %>% 
#   inner_join(raw_pgs) %>% 
#   filter(str_detect(IID, 'sample')) %>%
#   inner_join(p_dat1) %>% 
#   filter(str_detect(clean_name, 'Cognitive performance$|Schizo')) %>% 
#   filter(str_detect(clean_name, 'Cognitive performance$')) %>% 
#   filter(str_detect(name, '1|2|3')) %>% 
#   mutate(clean_ph = case_when(str_detect(name, '1') ~ 'Core language (F1)',
#                               str_detect(name, '2') ~ 'Receptive language (F2)',
#                               str_detect(name, '3') ~ 'Nonverbal IQ (F3)'),
#          clean_ph = factor(clean_ph, levels = c('Core language (F1)', 'Receptive language (F2)', 'Nonverbal IQ (F3)')),
#          lab = str_c('r = ', round(estimate, digits = 2), ', p-val = ', formatC(p.value, digits = 2))) %>% 
#   ggplot(aes(x = pgs_pc_corrected, y = value)) +
#   geom_point(size = 1.5) +
#   geom_smooth(method = 'lm') +
#   facet_wrap(~ clean_ph, nrow = 1) +
#   geom_text(aes(x = 0, y = 3, label = lab), check_overlap = TRUE, size = 3) +
#   xlab('Cognitive performance polygenic score') +
#   ylab('Factor score') +
#   theme_classic()

# ph %>% 
#   pivot_longer(cols = matches('Factor')) %>% 
#   inner_join(raw_pgs) %>% 
#   filter(str_detect(IID, 'sample')) %>%
#   inner_join(p_dat1) %>% 
#   filter(str_detect(clean_name, 'Cognitive performance$|Schizo')) %>% 
#   filter(str_detect(clean_name, 'Schizo')) %>% 
#   filter(str_detect(name, '1|2|3')) %>% 
#   mutate(clean_ph = case_when(str_detect(name, '1') ~ 'Core language (F1)',
#                               str_detect(name, '2') ~ 'Receptive language (F2)',
#                               str_detect(name, '3') ~ 'Nonverbal IQ (F3)'),
#          clean_ph = factor(clean_ph, levels = c('Core language (F1)', 'Receptive language (F2)', 'Nonverbal IQ (F3)')),
#          lab = str_c('r = ', round(estimate, digits = 2), ', p-val = ', formatC(p.value, digits = 2))) %>% 
#   ggplot(aes(x = pgs_pc_corrected, y = value)) +
#   geom_point(size = 1.5) +
#   geom_smooth(method = 'lm') +
#   facet_wrap(~ clean_ph, nrow = 1) +
#   geom_text(aes(x = 0, y = 3, label = lab), check_overlap = TRUE, size = 3) +
#   xlab('Schizophrenia polygenic score') +
#   ylab('Factor score') +
#   theme_classic()

p_dat <- read_csv('/wdata/lcasten/sli_wgs/prs/pathway_prs/results/fig_2_pathway_pgs_data_corrected.complement.csv') # %>% distinct(name)
p_dat2 <- read_csv('/wdata/lcasten/sli_wgs/prs/replication/SPARK_Lingo_verbal_memory_ES-PGS_genome_wide_corrected_data.csv')
p_dat3 <- read_csv('/wdata/lcasten/sli_wgs/prs/pathway_prs/results/fig_2_pathway_pgs_data_corrected_cell_type_specific.complement.csv')

unique(p_dat$pgs_nm)

p_dat1 %>% 
  filter(str_detect(clean_name, 'Cognitive performa')) # %>% 


## heatmap
p_ht <- p_dat1 %>%
  # distinct(clean_name) %>% unlist() %>% unname()
  mutate(type = ifelse(clean_name == 'Height', 'Misc.', type)) %>%
  filter(str_detect(clean_name, 'birth|Alzheim|Cannabis|Brainstem|Caudate|Amygdala|Pallidum|Putamen|accumbe|Thalamus|education skills|BMI', negate = TRUE)) %>% # distinct(clean_name) %>% unlist() %>% unname()
  group_by(name) %>%
  mutate(fdr = p.adjust(p.value, method = 'fdr')) %>%
  ungroup() %>%
  mutate(sig = case_when(fdr < 0.05 ~ '**',
                         fdr > 0.05 & p.value < 0.05 ~ '*',
                         TRUE ~ '')) %>%
  arrange(desc(clean_name)) %>% 
  mutate(clean_name = factor(clean_name, levels = unique(clean_name))) %>%
  mutate(clean_name = fct_relevel(.f = clean_name, "Cerebellar volume", after = 0)) %>%
  mutate(type = ifelse(type == 'Behavioral/SES', 'Personality', type),
         type = ifelse(clean_name %in% c('Townsend deprivation index', 'Income'), 'SES', type)) %>%
  mutate(type = factor(type, levels = c('Cognitive', 'Neuro/Psychiatric', 'Personality', 'SES', 'Brain MRI', 'Misc.'))) %>%
  ggplot(aes(x = str_remove_all(name, pattern = 'actor'), y = clean_name, fill = estimate)) +
  geom_tile() +
  geom_text(aes(label = sig), hjust = 0.5, size = 6.5) +
  facet_grid(rows = vars(type), scales = 'free_y', space = 'free_y') +
  # scale_fill_gradient2(low = 'darkviolet', mid = 'white', high = 'green3', midpoint = 0) +
  scale_fill_gradient2(low = 'dodgerblue', mid = 'white', high = 'chocolate1', midpoint = 0) +
  xlab(NULL) +
  ylab('Polygenic score') +
  labs(fill = 'Correlation:') +
  scale_x_discrete(expand = c(0,0), position = 'top') + 
  theme_classic() +  
  theme(legend.key.height = unit(1, "cm"),
        panel.border = element_rect(color = "black", fill = NA, size = 1),
        strip.background = element_rect(colour = "black", fill="white", size = 1),
        strip.placement = 'outside', 
        legend.position = 'right')

## scatterplots
p_dat <- p_dat %>% 
  mutate(pgs_nm = str_c('PGS in ', pgs_nm)) %>%
  filter(pgs_nm != 'PGS in primate conserved loci (65 Mya)') %>%
    mutate(type = case_when(pgs_nm == 'PGS in primate conserved loci (65 Mya)' ~ 'PGS in primate\nconserved loci',
                            pgs_nm == 'PGS in primate UCEs (65 Mya)' ~ 'PGS in primate UCEs',
                            pgs_nm == 'PGS in human-chimp divergent genes (6 Mya)' ~ 'PGS in human-chimp\ndivergent genes',
                            pgs_nm == 'PGS in HAQERs (600 Kya)' ~ 'PGS in HAQERs',
                            pgs_nm == 'PGS in Neanderthal selective sweep loci (50 Kya)' ~ 'PGS in Neanderthal\nselective sweep loci',
                            pgs_nm == 'PGS in recent selection loci (2-3 Kya)' ~ 'PGS in recent\nselection loci',
                            pgs_nm == 'PGS in pgs_genome_wide_baseline' ~ 'Genome-wide PGS',
                            TRUE ~ pgs_nm)) %>%
    mutate(type = factor(type, levels = c('PGS in primate\nconserved loci', 'PGS in primate UCEs', 'PGS in human-chimp\ndivergent genes', 'PGS in HAQERs', 'PGS in HARs', 'PGS in Neanderthal\nselective sweep loci', 'PGS in recent\nselection loci', 'Genome-wide PGS'))) %>% 
    mutate(lab = str_replace_all(lab, pattern = ', ', replacement = '\n')) %>%
    filter(str_detect(type, pattern = 'HAR', negate = TRUE))

pf1 <- p_dat %>% 
  filter(name %in% c('F1')) %>% 
  filter(str_detect(pgs_nm, 'chimp')) %>% 
  mutate(sig = case_when(p.value < 0.05 ~ TRUE,
                         TRUE ~ FALSE)) %>%
  ggplot(aes(x = pgs_val, y = value)) +
  geom_point(size = 2, aes(color = sig)) +
  geom_smooth(aes(color = sig), method = 'lm', size = 1.5) +
  # facet_wrap(~ type, nrow = 1) + # , scales  = 'free_x'
  xlab('Human-chimp divergent genes\nCognitive performance ES-PGS') +
  ylab('Core language (F1) score\nadjusted for background PGS') +
  theme_classic() +
  geom_text(aes(x = -2.325, y = 2.85, label = lab), size = 3.5, check_overlap = TRUE) +
  scale_color_manual(values = c('black')) +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 1),
        strip.background = element_rect(colour = "black", fill="white", size = 1),
        strip.placement = 'outside',
        strip.text = element_text(size = 12),
        legend.position = 'none')
pf1 %>% 
  ggsave(filename = '/wdata/lcasten/sli_wgs/paper_figures/subplots/human_chimp_divergence_ES-PGS_F1.png', device = 'png', units = 'in', dpi = 300, width = 5, height = 5)

pf3 <- p_dat %>% 
  filter(name %in% c('F3')) %>% 
  filter(str_detect(pgs_nm, 'chimp')) %>% 
  mutate(sig = case_when(p.value < 0.05 ~ TRUE,
                         TRUE ~ FALSE)) %>%
  ggplot(aes(x = pgs_val, y = value)) +
  geom_point(size = 2, aes(color = sig)) +
  geom_smooth(aes(color = sig), method = 'lm') +
  # facet_wrap(~ type, nrow = 1) + # , scales  = 'free_x'
  xlab('Cognitive performance ES-PGS') +
  ylab('NVIQ (F3) score\nadjusted for background PGS') +
  theme_classic() +
  geom_text(aes(x = -2.325, y = 2.85, label = lab), size = 3.5, check_overlap = TRUE) +
  scale_color_manual(values = c('grey90')) +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 1),
        strip.background = element_rect(colour = "black", fill="white", size = 1),
        strip.placement = 'outside',
        strip.text = element_text(size = 12),
        legend.position = 'none')

# p_dat2 <- p_dat %>% 
#   mutate(pgs_nm = str_c('PGS in ', pgs_nm))

p_f1 <- p_dat %>% 
  mutate(sig = case_when(p.value < 0.05 ~ TRUE,
                         TRUE ~ FALSE)) %>%
  filter(name == 'F1') %>%
  filter(type != 'Genome-wide PGS') %>%
  ggplot(aes(x = pgs_val, y = resid_value)) +
  geom_point(size = 2, aes(color = sig)) +
  geom_smooth(aes(color = sig), method = 'lm') +
  facet_wrap(~ type, nrow = 1) + # , scales  = 'free_x'
  xlab('Cognitive performance ES-PGS') +
  ylab('Core language (F1) score\nadjusted for background PGS') +
  theme_classic() +
  geom_text(aes(x = -2.325, y = 2.85, label = lab), size = 3.5, check_overlap = TRUE) +
  scale_color_manual(values = c('grey90', 'deeppink2')) +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 1),
        strip.background = element_rect(colour = "black", fill="white", size = 1),
        strip.placement = 'outside',
        strip.text = element_text(size = 12),
        legend.position = 'none') # +
#  coord_cartesian(ylim = c(-2.3, 3))
p_f2 <- p_dat %>% 
  mutate(sig = case_when(p.value < 0.05 ~ TRUE,
                         TRUE ~ FALSE)) %>%
  filter(name == 'F2') %>%
  filter(type != 'Genome-wide PGS') %>%
  ggplot(aes(x = pgs_val, y = resid_value)) +
  geom_point(size = 2, aes(color = sig)) +
  geom_smooth(aes(color = sig), method = 'lm') +
  facet_wrap(~ type, nrow = 1) + # , scales  = 'free_x'
  xlab('Cognitive performance ES-PGS') +
  ylab('Receptive language (F2) score\nadjusted for background PGS') +
  theme_classic() +  
  geom_text(aes(x = -2.3, y = 2.11, label = lab), size = 3.5, check_overlap = TRUE) +
  scale_color_manual(values = c('grey90', 'deeppink2')) +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 1),
        strip.background = element_rect(colour = "black", fill="white", size = 1),
        strip.placement = 'outside',
        strip.text = element_text(size = 12),
        legend.position = 'none')

p_rep <- p_dat2 %>% 
  filter(name == 'Factor2') %>%
  mutate(lab = str_c('Pearson r = ', round(estimate, digits = 2), '\np-val = ', formatC(p.value, digits = 2))) %>%
  mutate(sig = case_when(p.value < 0.05 ~ TRUE,
                         TRUE ~ FALSE)) %>%
  mutate(type = case_when(pgs_nm == 'consPrimates_UCE' ~ 'PGS in primate UCEs',
                          pgs_nm == 'human_chimp_div_DMG' ~ 'PGS in human-chimp\ndivergent genes',
                          pgs_nm == 'HAQER' ~ 'PGS in HAQERs',
                          pgs_nm == 'NeanderthalSelectiveSweep' ~ 'PGS in Neanderthal\nselective sweep loci',
                          pgs_nm == 'human_singleton_density_score_top5pct' ~ 'PGS in recent\nselection loci',
                          pgs_nm == 'pgs_genome_wide_baseline' ~ 'Genome-wide PGS',
                          TRUE ~ pgs_nm)) %>%
    mutate(type = factor(type, levels = c('PGS in primate\nconserved loci', 'PGS in primate UCEs', 'PGS in human-chimp\ndivergent genes', 'PGS in HAQERs', 'PGS in HARs', 'PGS in Neanderthal\nselective sweep loci', 'PGS in recent\nselection loci', 'Genome-wide PGS'))) %>% 
    mutate(lab = str_replace_all(lab, pattern = ', ', replacement = '\n')) %>%
    filter(str_detect(type, pattern = 'HAR', negate = TRUE)) %>%
    filter(pgs_nm != 'pgs_genome_wide_baseline') %>%
  # filter(name == 'F1') %>%
  ggplot(aes(x = pgs_pc_corrected, y = resid_norm_value)) +
  geom_point(size = 1.25, alpha = 0.7, aes(color = sig)) +
  geom_smooth(aes(color = sig), method = 'lm') +
  facet_wrap(~ type, nrow = 1) + # , scales  = 'free_x'
  xlab('Cognitive performance ES-PGS') +
  ylab('Replication cohort core language score\nadjusted for background PGS') +
  theme_classic() +
  geom_text(aes(x = -2.2, y = 2.5, label = lab), size = 3.5, check_overlap = TRUE) +
  scale_color_manual(values = c('grey90', 'darkgoldenrod3')) +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 1),
        strip.background = element_rect(colour = "black", fill="white", size = 1),
        strip.placement = 'outside',
        strip.text = element_text(size = 12),
        legend.position = 'none') # +
 # coord_cartesian(ylim = c(-2.6, 2.6))

mod_res <- read_csv('/wdata/lcasten/sli_wgs/prs/pathway_prs/results/evo_prs_lm_comparison_results.csv')
unique(mod_res$model)

mod_res %>% 
  filter(str_detect(factor, '1|2')) %>% 
  filter(p.value_model_comparison < 0.05 & annotation_beta > 0) %>% 
  select(factor, model, statistic, annotation_beta)

mod_res %>% 
  select(-c(3:6)) %>% 
  filter(model == 'consPrimates_UCE') %>% 
  as.data.frame

####
mod_res %>% 
  mutate(model = factor(model, levels = c('pgs_genome_wide_baseline', 'consPrimates_65mya', 'consPrimates_UCE', 'LinAR_Simiformes', 'LinAR_Catarrhini', 'LinAR_Hominoidea', 'LinAR_Hominidae', 'LinAR_Homininae', 'LinAR_Human', 'human_chimp_div_DMG', 'HAQER', 'HAR', 'NeanderthalSelectiveSweep', 'human_singleton_density_score_top5pct'))) %>%
  drop_na(model) %>%
  mutate(rsq_gain_per_1000indSNP = (baseline_plus_anno_rsq - baseline_rsq) / (annotation_PGS_n_snp / 1000)) %>% filter(factor == 'Factor1' & model == 'HAQER') %>% relocate(rsq_gain_per_1000indSNP) %>% 
  mutate(rsq_gain_per_1000indSNP = rsq_gain_per_1000indSNP / baseline_rsq)


####
p_dat_time <- mod_res %>% 
  mutate(model = factor(model, levels = c('pgs_genome_wide_baseline', 'consPrimates_65mya', 'consPrimates_UCE', 'LinAR_Simiformes', 'LinAR_Catarrhini', 'LinAR_Hominoidea', 'LinAR_Hominidae', 'LinAR_Homininae', 'LinAR_Human', 'human_chimp_div_DMG', 'HAQER', 'HAR', 'NeanderthalSelectiveSweep', 'human_singleton_density_score_top5pct'))) %>%
  drop_na(model) %>%
  mutate(rsq_gain_per_1000indSNP = (baseline_plus_anno_rsq - baseline_rsq) / (annotation_PGS_n_snp / 1000)) %>% # filter(factor == 'Factor1' & model == 'HAQER') %>% relocate(rsq_gain_per_1000indSNP)
  mutate(model = str_replace_all(model, pattern = 'consPrimates_65mya', replacement = 'primate conserved loci (65 Mya)'),
         model = str_replace_all(model, pattern = 'consPrimates_UCE', replacement = 'primate UCEs (65 Mya)'),
         model = str_replace_all(model, pattern = 'LinAR_Simiformes', replacement = 'Simian accelerated regions  (> 30 Mya)'),
         model = str_replace_all(model, pattern = 'LinAR_Catarrhini', replacement = 'Old world monkey accelerated regions (30 Mya)'),
         model = str_replace_all(model, pattern = 'LinAR_Hominoidea', replacement = 'Ape accelerated regions (20 Mya)'),
         model = str_replace_all(model, pattern = 'LinAR_Hominidae', replacement = 'Great ape accelerated regions (15 Mya)'),
         model = str_replace_all(model, pattern = 'LinAR_Homininae', replacement = 'Hominid accelerated regions (9 Mya)'),
         model = str_replace_all(model, pattern = 'LinAR_Human', replacement = 'Human lineage accelerated regions (8 Mya)'),
         model = str_replace_all(model, pattern = 'human_chimp_div_DMG', replacement = 'human-chimp divergent genes (8 Mya)'),
         model = str_replace_all(model, pattern = 'NeanderthalSelectiveSweep', replacement = 'Neanderthal selective sweep loci (50 Kya)'),
         model = str_replace_all(model, pattern = 'HAR', replacement = 'HARs (600 Kya)'),
         model = str_replace_all(model, pattern = 'HAQER', replacement = 'HAQERs (600 Kya)'),
         model = str_replace_all(model, pattern = 'human_singleton_density_score_top5pct', replacement = 'recent selection loci (2-3 Kya)')
        ) %>%
  mutate(years_ago = case_when(model == 'primate UCEs (65 Mya)' ~ -65000000,
                               model == 'Simian accelerated regions  (> 30 Mya)' ~ -35000000,
                               model == 'Old world monkey accelerated regions (30 Mya)' ~ -30000000,
                               model == 'Ape accelerated regions (20 Mya)' ~ -20000000,
                               model == 'Great ape accelerated regions (15 Mya)' ~ -14000000,
                               model == 'Hominid accelerated regions (9 Mya)' ~ -9000000,
                               model == 'human-chimp divergent genes (8 Mya)' ~ -8000000,
                               model == 'HAQERs (600 Kya)' ~ -7400000,
                               model == 'HARs (600 Kya)' ~ -7000000,
                               model == 'Neanderthal selective sweep loci (50 Kya)' ~ -50000,
                               model == 'recent selection loci (2-3 Kya)' ~ -2000),
          years_ago_upper = case_when(model == 'primate UCEs (65 Mya)' ~ -70000000,
                               model == 'Simian accelerated regions  (> 30 Mya)' ~ -42000000,
                               model == 'Old world monkey accelerated regions (30 Mya)' ~ -30000000,
                               model == 'Ape accelerated regions (20 Mya)' ~ -20000000,
                               model == 'Great ape accelerated regions (15 Mya)' ~ -14000000,
                               model == 'Hominid accelerated regions (9 Mya)' ~ -9000000,
                               model == 'human-chimp divergent genes (8 Mya)' ~ -8000000,
                               model == 'HAQERs (600 Kya)' ~ -8000000,
                               model == 'HARs (600 Kya)' ~ -8000000,
                               model == 'Neanderthal selective sweep loci (50 Kya)' ~ -50000,
                               model == 'recent selection loci (2-3 Kya)' ~ -3000),
          years_ago_lower = case_when(model == 'primate UCEs (65 Mya)' ~ -60000000,
                               model == 'Simian accelerated regions  (> 30 Mya)' ~ -35000000,
                               model == 'Old world monkey accelerated regions (30 Mya)' ~ -25000000,
                               model == 'Ape accelerated regions (20 Mya)' ~ -18000000,
                               model == 'Great ape accelerated regions (15 Mya)' ~ -12000000,
                               model == 'Hominid accelerated regions (9 Mya)' ~ -8500000,
                               model == 'human-chimp divergent genes (8 Mya)' ~ -6000000,
                               model == 'HAQERs (600 Kya)' ~ -6000000,
                               model == 'HARs (600 Kya)' ~ -6000000,
                               model == 'Neanderthal selective sweep loci (50 Kya)' ~ -40000,
                               model == 'recent selection loci (2-3 Kya)' ~ -2000)) %>%
  mutate(years_since_split = case_when(model == 'primate UCEs (65 Mya)' ~ 1,
                               model == 'Simian accelerated regions  (> 30 Mya)' ~ 65000000-35000000,
                               model == 'Old world monkey accelerated regions (30 Mya)' ~ 35000000-30000000,
                               model == 'Ape accelerated regions (20 Mya)' ~ 30000000-20000000,
                               model == 'Great ape accelerated regions (15 Mya)' ~ 20000000-14000000,
                               model == 'Hominid accelerated regions (9 Mya)' ~ 14000000-9000000,
                               model == 'human-chimp divergent genes (8 Mya)' ~ 9000000-8000000,
                               model == 'HAQERs (600 Kya)' ~ 8000000-7400000,
                               model == 'HARs (600 Kya)' ~ 7400000-7000000,
                               model == 'Neanderthal selective sweep loci (50 Kya)' ~ 7000000-50000,
                               model == 'recent selection loci (2-3 Kya)' ~ 7000000-2000),
          years_since_pca = case_when(model == 'primate UCEs (65 Mya)' ~ 1,
                               model == 'Simian accelerated regions  (> 30 Mya)' ~ 65000000-35000000,
                               model == 'Old world monkey accelerated regions (30 Mya)' ~ 65000000-30000000,
                               model == 'Ape accelerated regions (20 Mya)' ~ 65000000-20000000,
                               model == 'Great ape accelerated regions (15 Mya)' ~ 65000000-14000000,
                               model == 'Hominid accelerated regions (9 Mya)' ~ 65000000-9000000,
                               model == 'human-chimp divergent genes (8 Mya)' ~ 65000000-8000000,
                               model == 'HAQERs (600 Kya)' ~ 65000000-7400000,
                               model == 'HARs (600 Kya)' ~ 65000000-7000000,
                               model == 'Neanderthal selective sweep loci (50 Kya)' ~ 65000000-50000,
                               model == 'recent selection loci (2-3 Kya)' ~ 65000000-2000),
          paper = case_when(model == 'primate UCEs (65 Mya)' ~ 'Kuderna, 2024',
                               model == 'Simian accelerated regions  (> 30 Mya)' ~ 'Bi, 2023',
                               model == 'Old world monkey accelerated regions (30 Mya)' ~ 'Bi, 2023',
                               model == 'Ape accelerated regions (20 Mya)' ~ 'Bi, 2023',
                               model == 'Great ape accelerated regions (15 Mya)' ~ 'Bi, 2023',
                               model == 'Hominid accelerated regions (9 Mya)' ~ 'Bi, 2023',
                               model == 'human-chimp divergent genes (8 Mya)' ~ 'Human-chimp divergence - Gokhman, 2020',
                               model == 'HAQERs (600 Kya)' ~ 'HAQERs - Magnan, 2022',
                               model == 'HARs (600 Kya)' ~ 'Whalen, 2023',
                               model == 'Neanderthal selective sweep loci (50 Kya)' ~ 'Green, 2010',
                               model == 'recent selection loci (2-3 Kya)' ~ 'Field, 2016')) %>%
  filter(model %in% c('primate UCEs (65 Mya)', 'Simian accelerated regions  (> 30 Mya)', 'Old world monkey accelerated regions (30 Mya)', 'Ape accelerated regions (20 Mya)', 'Great ape accelerated regions (15 Mya)', 'Hominid accelerated regions (9 Mya)', 'human-chimp divergent genes (8 Mya)', 'Human lineage accelerated regions (8 Mya)', 'HAQERs (600 Kya)', 'HARs (600 Kya)', 'Neanderthal selective sweep loci (50 Kya)', 'recent selection loci (2-3 Kya)')) %>%
  mutate(rsq_gain_per_1000indSNP = ifelse(annotation_beta < 0 | p.value_model_comparison > 0.05, 0, rsq_gain_per_1000indSNP)) %>%
  mutate(model = factor(model, levels = c('Genome-wide PGS', 'primate conserved loci (65 Mya)', 'primate UCEs (65 Mya)', 'Simian accelerated regions  (> 30 Mya)', 'Old world monkey accelerated regions (30 Mya)', 'Ape accelerated regions (20 Mya)', 'Great ape accelerated regions (15 Mya)', 'Hominid accelerated regions (9 Mya)', 'human-chimp divergent genes (8 Mya)', 'Human lineage accelerated regions (8 Mya)', 'HAQERs (600 Kya)', 'HARs (600 Kya)', 'Neanderthal selective sweep loci (50 Kya)', 'recent selection loci (2-3 Kya)'))) %>% 
  filter(model %in% c('primate UCEs (65 Mya)', 'Simian accelerated regions  (> 30 Mya)', 'Old world monkey accelerated regions (30 Mya)', 'Ape accelerated regions (20 Mya)', 'Great ape accelerated regions (15 Mya)', 'Hominid accelerated regions (9 Mya)', 'human-chimp divergent genes (8 Mya)', 'Human lineage accelerated regions (8 Mya)', 'HAQERs (600 Kya)', 'HARs (600 Kya)', 'Neanderthal selective sweep loci (50 Kya)', 'recent selection loci (2-3 Kya)')) %>%
  filter(factor %in% c('Factor1', 'Factor2', 'Factor3')) %>% 
  mutate(factor = str_remove_all(factor, 'actor')) %>%
  mutate(factor = case_when(factor == 'F1' ~ 'Core language (F1)',
                            factor == 'F2' ~ 'Receptive lang. (F2)',
                            factor == 'F3' ~ 'NVIQ (F3)')) %>%
  filter(factor != 'NVIQ (F3)') %>%
  drop_na(model, years_since_pca) %>%
  mutate(sig = case_when(annotation_beta > 0 & p.value_model_comparison < 0.05 ~ TRUE,
                         TRUE ~ FALSE))
p_dat_time
unique(p_dat_time$model)
p_dat_time %>% 
  filter(str_detect(factor, '1')) %>% 
  relocate(annotation_beta, p.value_model_comparison, rsq_gain_per_1000indSNP)

p_dat_time2 <- p_dat_time %>%
  drop_na(model, years_since_pca, factor) %>% 
  # mutate(range_upper = log10(years_ago_upper + max(abs(years_ago_upper)) + 1),
        #  range_lower = log10(years_ago_lower + max(abs(years_ago_upper)) + 1)) %>%
  mutate(paper = str_remove_all(paper, pattern = 'HAQERs - |Human-chimp divergence - ')) %>%
  mutate(lab = ifelse(paper == 'Magnan, 2022', 'HAQERs', NA_character_),
         lab = case_when(str_detect(model, 'HAR') ~ 'HARs',
                         str_detect(model, 'divergent') ~ 'Human-chimp divergence',
                         TRUE ~ lab)) %>%
  select(model, factor, rsq_gain_per_1000indSNP, baseline_rsq, annotation_beta, p.value_model_comparison, years_ago, years_ago_lower, years_ago_upper, paper, lab) %>% 
  mutate(years_ago_upper = -1 * years_ago_upper / 1000000,
         years_ago_lower = -1 * years_ago_lower / 1000000,
         years_ago = ((years_ago_upper + years_ago_lower) / 2),
         improvement = rsq_gain_per_1000indSNP / baseline_rsq) %>% 
  mutate(improvement = ifelse(annotation_beta >= 0, improvement, 0)) %>%
  mutate(sig = case_when(improvement > 0 & p.value_model_comparison < 0.05 ~ TRUE,
                         TRUE ~ FALSE)) %>%
  filter(factor == 'Core language (F1)') %>%
  filter(str_detect(model, 'divergen|uman acceler|HAR|chimp'))
p_dat_time2 %>% 
  select(model, factor, paper, lab)

p_evo <- p_dat_time %>%
  drop_na(model, years_since_pca, factor) %>% 
  # mutate(range_upper = log10(years_ago_upper + max(abs(years_ago_upper)) + 1),
        #  range_lower = log10(years_ago_lower + max(abs(years_ago_upper)) + 1)) %>%
  mutate(paper = str_remove_all(paper, pattern = 'HAQERs - |Human-chimp divergence - ')) %>%
  mutate(lab = ifelse(paper == 'Magnan, 2022', 'HAQERs', NA_character_)) %>%
  select(model, factor, rsq_gain_per_1000indSNP, baseline_rsq, annotation_beta, p.value_model_comparison, years_ago, years_ago_lower, years_ago_upper, paper, lab) %>% 
  mutate(years_ago_upper = -1 * years_ago_upper / 1000000,
         years_ago_lower = -1 * years_ago_lower / 1000000,
         years_ago = ((years_ago_upper + years_ago_lower) / 2),
         improvement = rsq_gain_per_1000indSNP / baseline_rsq) %>% 
  mutate(improvement = ifelse(annotation_beta >= 0, improvement, 0)) %>%
  mutate(sig = case_when(improvement > 0 & p.value_model_comparison < 0.05 ~ TRUE,
                         TRUE ~ FALSE)) %>%
  filter(factor == 'Core language (F1)') %>%
  arrange(desc(improvement)) %>% 
  group_by(years_ago) %>% 
  slice_head(n = 1) %>%
  ungroup() %>%
  ggplot(aes(x = years_ago, y = improvement)) +
  geom_line(color = "grey65", size = 1, linetype = 'dashed') +      # Line for evolutionary events
  geom_linerange(aes(xmin = years_ago_lower, xmax = years_ago_upper), size = 1.2, color = 'black') +
  geom_point(size = 5, aes(shape = sig, color = paper)) +
  scale_shape_manual(values = c(21, 16)) +
  ggimage::geom_image(aes(x = 65, y = .3, image = '/wdata/lcasten/sli_wgs/primate.png'), size = .37, by="height") +
  ggimage::geom_image(aes(x = 10, y = .3, image = '/wdata/lcasten/sli_wgs/human_ancestor.png'), size = .375, by="height") +
  ggimage::geom_image(aes(x = 0, y = .3, image = '/wdata/lcasten/sli_wgs/human.png'), size = .4, by="height") +
  scale_x_reverse(limits = c(70, 0), breaks = seq(0, 70, 10)) +       # Reverse the x-axis for time scale
  labs(x = "Evolutionary timeline (millions of years ago)",   # X-axis label
       y = 'ES-PGS relative R-squared improvement\nper 1000 independent SNPs for F1',            # Y-axis label
       ) +
  theme_minimal() +
  scale_y_continuous(labels = scales::percent) +
  ggrepel::geom_text_repel(aes(label = lab), nudge_y = 0.0185, size = 5) +
  labs(shape = "ES-PGS model improvement p < 0.05:", color = 'Annotation source:') +
  theme(legend.position = 'bottom') +
  geom_point(data = p_dat_time2, size = 5, aes(shape = sig, color = paper)) +
  geom_linerange(data = p_dat_time2, aes(xmin = years_ago_lower, xmax = years_ago_upper), size = 1.2, color = 'black') +
  ggrepel::geom_text_repel(data = p_dat_time2[p_dat_time2$lab != 'HARs',], aes(label = lab), nudge_y = 0.02, size = 5) +
  ggrepel::geom_text_repel(data = p_dat_time2[p_dat_time2$lab == 'HARs',], aes(label = lab), nudge_y = 0.0025, nudge_x = 2.25, size = 5)

#####################
unique(p_dat3$pgs_nm)
ct_res <- read_csv('/wdata/lcasten/sli_wgs/prs/pathway_prs/results/evo_prs_cell_type_lm_comparison_results.csv')
ct_res %>% 
  filter(factor == 'Factor1') %>% 
  select(-factor) %>% 
  filter(str_detect(model, '^L[0-9]')) %>% 
  select(p.value_model_comparison, annotation_beta, model) %>% 
  arrange(p.value_model_comparison)
ct_dat <- p_dat3 %>% 
  filter(str_detect(pgs_nm, pattern = 'human_specific_evolution_brain_expression_')) %>% 
  mutate(pgs_nm = str_remove_all(pgs_nm, pattern = 'human_specific_evolution_brain_expression_')) %>% 
  filter(pgs_nm %in% c('L3_5_RORB_1', 'L3_5_RORB_2', 'L4_6_RORB_2', 'L5_6_THEMIS_1', 'L5_6_THEMIS_2', 'L3_5_RORB_3')) %>% 
  mutate(type = case_when(pgs_nm %in% c('L3_5_RORB_1', 'L3_5_RORB_2', 'L5_6_THEMIS_2') ~ ' (human-specific FOXP2 target DGE)',
                          pgs_nm %in% c('L4_6_RORB_2', 'L5_6_THEMIS_1') ~ ' (human-specific FOXP2 DGE)',
                          TRUE ~ ' (Control cell-type)'))

p_ct <- ct_dat %>% 
  filter(name == 'F1') %>%
  mutate(sig = case_when(p.value < 0.05 ~ TRUE,
                         TRUE ~ FALSE)) %>%
  mutate(pgs_nm = str_c(pgs_nm, type)) %>%
  mutate(pgs_nm = factor(pgs_nm, levels = c('L3_5_RORB_3 (Control cell-type)', 'L4_6_RORB_2 (human-specific FOXP2 DGE)', 'L5_6_THEMIS_1 (human-specific FOXP2 DGE)', 'L3_5_RORB_1 (human-specific FOXP2 target DGE)', 'L3_5_RORB_2 (human-specific FOXP2 target DGE)', 'L5_6_THEMIS_2 (human-specific FOXP2 target DGE)'))) %>%
  ggplot(aes(x = pgs_val, y = resid_value)) +
  geom_point(size = 2, aes(color = sig)) +
  geom_smooth(aes(color = sig), method = 'lm') +
  facet_wrap(~ pgs_nm, nrow = 2) + # , scales  = 'free_x'
  xlab('Cognitive performance ES-PGS') +
  ylab('Core language (F1) score adjusted for background PGS') +
  theme_classic() +  
  geom_text(aes(x = -2.3, y = 3, label = lab), size = 3.5, check_overlap = TRUE) +
  scale_color_manual(values = c('grey90', 'blue2')) +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 1),
        strip.background = element_rect(colour = "black", fill="white", size = 1),
        strip.placement = 'outside',
        legend.position = 'none',
        axis.title = element_text(size = 14), 
        axis.text = element_text(size = 12), 
        strip.text = element_text(size = 14), 
        legend.text=element_text(size=12), 
        legend.title = element_text(size = 14), 
        plot.title = element_text(size = 16))

unique(ct_dat$pgs_nm)




##################################
## merge plots to make figure
##################################
library(patchwork)
## +  plot_layout(widths = c(2, 2, 1), nrow = 1) 
# p <- p_ht + ((p_f1 / p_f2 / p_rep) / p_evo_time)
# p <- (p_ht + (p_f1 / p_f2 / p_rep) + plot_layout(widths = c(2,6))) / p_evo_time +  plot_layout(widths = c(1, .7), heights = c(2, 1))  + plot_annotation(tag_levels = 'A') & theme(plot.tag = element_text(size = 12))
# p <- (p_ht + ((p_f1 / p_f2 / p_rep) / p_evo_time) + plot_layout(widths = c(2,6))) + plot_annotation(tag_levels = 'A') & theme(plot.tag = element_text(size = 14), axis.title = element_text(size = 12), axis.text = element_text(size = 11), strip.text = element_text(size = 12), legend.text=element_text(size=12), legend.title = element_text(size = 14))
p <- (p_ht + ((p_f1 / p_f2 / p_rep) / p_evo) + plot_layout(widths = c(3,6))) + plot_annotation(tag_levels = 'A') & theme(plot.tag = element_text(size = 16), axis.title = element_text(size = 14), axis.text = element_text(size = 12), strip.text = element_text(size = 12), legend.text=element_text(size=12), legend.title = element_text(size = 14))

# p[[2]] <- p[[2]] + theme(axis.title.y = element_text(margin = margin(t = 0, r = -200, b = 0, l = 0)), legend.position = 'bottom')  #theme(axis.title.y = element_text(hjust = 0))

## plot w/ bar plot for relative perf increase and replication figure
# p2 <- (p_ht + (p_f1 / p_f2) + plot_layout(widths = c(1,6))) / p_perf +  plot_layout(widths = c(1, .7), heights = c(2, 1))  + plot_annotation(tag_levels = 'A') & theme(plot.tag = element_text(size = 12))
# p
# p3 <- (p_ht + (p_f1 / p_f2) + plot_layout(widths = c(1,6))) + plot_annotation(tag_levels = 'A') & theme(plot.tag = element_text(size = 12))

ggsave(p, filename = '/wdata/lcasten/sli_wgs/paper_figures/fig2_pgs_correlation.png', device = 'png', dpi = 300, units = 'in', width = 29, height = 19)
ggsave(p_ct, filename = '/wdata/lcasten/sli_wgs/paper_figures/fig_cell_type_DEG_pgs_correlation.png', device = 'png', dpi = 300, units = 'in', width = 17, height = 12)





############################################
## make plot for cell-type stuff
ct_res <- read_csv('/wdata/lcasten/sli_wgs/prs/pathway_prs/results/evo_prs_cell_type_lm_comparison_results.csv')
unique(ct_res$model)

ct_res %>% 
  filter(factor == 'Factor1') %>%
  filter(model == 'L3_5_RORB_1') %>% 
  relocate(matches('annotation'))

p_deg <- ct_res %>% 
  filter(factor == 'Factor1') %>%
  mutate(stat = sign(annotation_beta) * -log10(p.value_model_comparison)) %>% 
  arrange(stat) %>% 
  mutate(model = factor(model, levels = unique(model))) %>%
  mutate(sig = ifelse(p.value_model_comparison < 0.05, TRUE, FALSE)) %>%
  mutate(type = case_when(str_detect(model, pattern = '^L[0-9]|Upper') ~ 'Excitatory',
                          str_detect(model, pattern = 'VIP|SST|PVALB|LAMP') ~ 'Inhibitory',
                          str_detect(model, pattern = 'Oligo|OPC|Astro|Micro|OL') ~ 'Non-neuronal')) %>%
  drop_na(type) %>%
  ggplot(aes(x = stat, y = model, fill = sig)) +
  geom_bar(stat = 'identity') +
  geom_vline(xintercept = -log10(0.05), color = 'red', linetype = 'dashed', size = 1.075) +
  # geom_vline(xintercept = -1 * -log10(0.05), color = 'red', linetype = 'dashed', size = 1.075) +
  ggtitle('Human specific differentially expressed genes') +
  xlab('F1 ES-PGS signed -log10(p-value)') +
  ylab('Cell-type') +
  scale_fill_manual(values = c('grey85', 'black')) +
  theme_classic() +
  theme(legend.position = 'none') +
  facet_grid(rows = vars(type), space = 'free_y', scales = 'free_y')

atac_res <- read_csv('/wdata/lcasten/sli_wgs/prs/pathway_prs/results/evo_prs_cell_type_scATAC_lm_comparison_results.csv')
unique(atac_res$model)
p_dar <- atac_res %>% 
  filter(factor == 'Factor1') %>% 
  mutate(stat = sign(annotation_beta) * -log10(p.value_model_comparison)) %>% 
  arrange(stat) %>% 
  mutate(model = factor(model, levels = unique(model))) %>%
  mutate(sig = ifelse(p.value_model_comparison < 0.05, TRUE, FALSE)) %>%
  mutate(type = case_when(str_detect(model, pattern = '^L[0-9]|Upper') ~ 'Excitatory',
                          str_detect(model, pattern = 'VIP|SST|PVALB|LAMP') ~ 'Inhibitory',
                          str_detect(model, pattern = 'Oligo|OPC|Astro|Micro') ~ 'Non-neuronal')) %>%
  ggplot(aes(x = stat, y = model, fill = sig)) +
  geom_bar(stat = 'identity') +
  geom_vline(xintercept = -log10(0.05), color = 'red', linetype = 'dashed', size = 1.075) +
  # geom_vline(xintercept = -1 * -log10(0.05), color = 'red', linetype = 'dashed', size = 1.075) +
  ggtitle('Human specific differentially accessible regions') +
  xlab('ES-PGS signed -log10(p-value)') +
  ylab('Cell-type') +
  scale_fill_manual(values = c('grey85', 'black')) +
  theme_classic() +
  theme(legend.position = 'none') +
  facet_grid(rows = vars(type), space = 'free_y', scales = 'free_y')

library(patchwork)
p_ct <- p_deg + p_dar + plot_annotation(tag_levels = 'A') & theme(plot.tag = element_text(size = 16), axis.title = element_text(size = 14), axis.text = element_text(size = 12), strip.text = element_text(size = 12), legend.text=element_text(size=12), legend.title = element_text(size = 14), plot.title = element_text(size = 16))
ggsave(p_ct, filename = '/wdata/lcasten/sli_wgs/paper_figures/fig_es_pgs_cell_type.png', device = 'png', units = 'in', dpi = 300, width = 16, height = 10)



##################################




########################################
## old stuff
# p_x_val <- p_dat_time %>% 
#   select(model, matches('years')) %>% 
#   distinct() %>% 
#   drop_na(years_ago) %>% 
#   arrange(years_ago) %>% 
#   mutate(xvals = log10(years_since_pca))
# p_x_val[p_x_val$model == 'HAQERs (600 Kya)',] %>% 
#   as.data.frame()


# p_perf <- mod_res %>% 
#   mutate(pgs_nm = case_when(model == 'consPrimates_65mya' ~ 'primate conserved\nloci (65 Mya)',
#                             model == 'consPrimates_UCE' ~ 'primate UCEs (65 Mya)',
#                             model == 'human_chimp_div_DMG' ~ 'human-chimp\ndivergent genes (8 Mya)',
#                             model == 'HAQER' ~ 'HAQERs (600 Kya)',
#                             model == 'HAR' ~ 'HARs (600 Kya)',
#                             model == 'NeanderthalSelectiveSweep' ~ 'Neanderthal selective\nsweep loci (50 Kya)',
#                             model == 'human_singleton_density_score_top5pct' ~ 'recent selection\nloci (2-3 Kya)'),
#          pgs_nm = str_c('PGS in ', pgs_nm)) %>%
#   mutate(pgs_nm = factor(pgs_nm, levels = c('PGS in primate conserved\nloci (65 Mya)', 'PGS in primate UCEs (65 Mya)', 'PGS in human-chimp\ndivergent genes (8 Mya)', 'PGS in HAQERs (600 Kya)', 'PGS in HARs (600 Kya)', 'PGS in Neanderthal selective\nsweep loci (50 Kya)', 'PGS in recent selection\nloci (2-3 Kya)'))) %>% 
#   filter(factor %in% c('Factor1', 'Factor2')) %>% 
#   mutate(factor = str_remove_all(factor, 'actor')) %>%
#   drop_na() %>%
#   filter(str_detect(pgs_nm, pattern = 'HAR|primate conserved', negate = TRUE)) %>% 
#   ggplot(aes(x = pgs_nm, y = rsq_gain_per_1000indSNP, color = factor, group = factor)) +
#   geom_point(size = 5) +
#   geom_line(size = 1.1) +
#   xlab('ES-PGS annotation') +
#   ylab('R-squared improvement per\n1000 independent SNPs') +
#   scale_y_continuous(labels = scales::percent) +
#   geom_hline(yintercept = 0, color = 'red2', linetype = 'dashed', size = 1.075) +
#   labs(color = 'Language factor') +
#   theme_classic() +
#   theme(legend.position = 'bottom') 
# p_perf2 <- mod_res %>% 
#   mutate(pgs_nm = case_when(model == 'consPrimates_65mya' ~ 'primate conserved\nloci (65 Mya)',
#                             model == 'consPrimates_UCE' ~ 'primate UCEs (65 Mya)',
#                             model == 'human_chimp_div_DMG' ~ 'human-chimp\ndivergent genes (8 Mya)',
#                             model == 'HAQER' ~ 'HAQERs (600 Kya)',
#                             model == 'HAR' ~ 'HARs (600 Kya)',
#                             model == 'NeanderthalSelectiveSweep' ~ 'Neanderthal selective\nsweep loci (50 Kya)',
#                             model == 'human_singleton_density_score_top5pct' ~ 'recent selection\nloci (2-3 Kya)'),
#          pgs_nm = str_c('PGS in ', pgs_nm)) %>%
#   mutate(pgs_nm = factor(pgs_nm, levels = c('PGS in primate conserved\nloci (65 Mya)', 'PGS in primate UCEs (65 Mya)', 'PGS in human-chimp divergent\ngenes (8 Mya)', 'PGS in Neanderthal selective\nsweep loci (50 Kya)', 'PGS in HAQERs (600 Kya)', 'PGS in HARs (600 Kya)', 'PGS in recent selection\nloci (2-3 Kya)'))) %>% 
#   filter(factor %in% c('Factor1', 'Factor2')) %>% 
#   mutate(factor = str_remove_all(factor, 'actor')) %>%
#   filter(model == 'HAQER') %>%
#   drop_na() %>%
#   ggplot(aes(x = factor, y = rsq_gain_per_1000indSNP / baseline_rsq)) +
#   geom_bar(stat= 'identity', width = 0.7) +
#   xlab('Cognitive performance HAQER ES-PGS') +
#   ylab('Relative R-squared improvement per\n1000 independent SNPs') +
#   scale_y_continuous(labels = scales::percent) +
#   theme_classic() 
# unique(p_dat1$type)
# p_evo_time <- p_dat_time %>%
#   drop_na(model, years_since_pca, factor) %>% 
#   mutate(range_upper = log10(years_ago_upper + max(abs(years_ago_upper)) + 1),
#          range_lower = log10(years_ago_lower + max(abs(years_ago_upper)) + 1)) %>%
#   mutate(paper = str_remove_all(paper, pattern = 'HAQERs - |Human-chimp divergence - ')) %>%
#   mutate(lab = ifelse(paper == 'Magnan, 2022', 'HAQERs', NA_character_)) %>%
#   # ggplot(aes(x = mode l, y = rsq_gain_per_1000indSNP / baseline_rsq, color = factor, group = factor)) +
#   ggplot(aes(x = log10(years_since_pca), y = rsq_gain_per_1000indSNP / baseline_rsq)) +
#   geom_linerange(aes(xmin = -1 * range_lower, xmax = -1 * range_upper)) +
#   # geom_line(size = 1.1, alpha = 0.75) +
#   xlab('Human evolutionary timeline') +
#   ylab('ES-PGS relative R-squared improvement\nper 1000 independent SNPs') +
#   # scale_x_continuous(, labels = c('primate UCEs (65 Mya)', 'Simian accelerated regions  (> 30 Mya)', 'Old world monkey accelerated regions (30 Mya)', 'Ape accelerated regions (20 Mya)', 'Great ape accelerated regions (15 Mya)', 'Hominid accelerated regions (9 Mya)', 'human-chimp divergent genes (8 Mya)', 'Human lineage accelerated regions (8 Mya)', 'HAQERs (600 Kya)', 'HARs (600 Kya)', 'Neanderthal selective sweep loci (50 Kya)', 'recent selection loci (2-3 Kya)')) +
#   scale_y_continuous(labels = scales::percent) +
#   # geom_hline(yintercept = 0, color = 'red2', linetype = 'dashed', size = 1.075) +
#   theme_classic() +
#   theme(
#         axis.text.y.right = element_blank(),
#         axis.ticks.y.right = element_blank(),
#         axis.line.y.right = element_blank(),
#         legend.position = 'bottom',
#         strip.background.x = element_blank(),
#         strip.text.x = element_blank()) +
#   # ggrepel::geom_label_repel(aes(label = paper), nudge_y = .021, color = 'black') +
#   ggrepel::geom_text_repel(aes(label = lab), nudge_x = 0.01, size = 4.5) +
#   # geom_text(aes(y = 0.4, label = paper), check_overlap = TRUE) +
#   # geom_text(data = data.frame(years_since_pca = c(7.760422+ 0.01), rsq_gain_per_1000indSNP = c(.4)), aes(label = 'TEST')) +
#   geom_blank(data = data.frame(years_since_pca = c(1.5, 29950000), rsq_gain_per_1000indSNP = Inf, sig = FALSE)) +
#   # ggforce::facet_col(factor ~ years_since_pca > 2, , space = "free", scales = "free_x", strip.position = 'top') +
#   facet_grid(cols = vars(years_since_pca > 2), rows = vars(factor), space = "free_x", scales = "free_x", drop = TRUE) +
#   scale_x_continuous(breaks = p_x_val$xvals, labels = c('65 Mya', '35 Mya', '30 Mya', '20 Mya', '14 Mya', '', '', '', '7 Mya', '50 Kya', '')) +
#   # ggbreak::scale_x_break(c(log10(1.5), log10(29000000)), space = 0.4) +
#   labs(shape = "ES-PGS model improvement p-val < 0.05:", color = 'Annotation source:') +
#   geom_linerange(xmin = 0.0005, xmax = 6, size = 1.1, alpha = 0.75) +
#   geom_point(size = 5, aes(shape = sig, color = paper), alpha = 0.8) +
#   scale_shape_manual(values = c(1, 16)) # +
#   # geom_text(aes(x = .1, y = 0.4, label = 'Primate common ancestor'), check_overlap = TRUE) +
#   # geom_text(aes(x = 7.71, y = 0.4, label = 'Great ape split')) +
#   # geom_text(aes(x = 7.76, y = 0.4, label = 'Human-chimp divergence')) +
#   # geom_text(aes(x = 7.81, y = 0.4, label = 'Neanderthal selective sweep')) +
#   # geom_hline(data=filter(p_dat_time, years_since_pca < 2), aes(yintercept=0), colour="black") 
# ggsave(p_evo_time, filename = '/wdata/lcasten/sli_wgs/paper_figures/fig2_evo_time_improvement.png', device = 'png', dpi = 300, units = 'in', width = 14, height = 7)


# ###
# data <- data.frame(
#   time = c(70, 60, 50, 40, 30, 20, 10, 0),  # Time in millions of years ago
#   importance = c(1, 3, 2, 4, 5, 3, 4, 2)   # Importance of events
#   )
# ggplot(data, aes(x = time, y = importance)) +
#   geom_line(color = "blue", size = 1) +      # Line for evolutionary events
#   geom_point(color = "red", size = 3) +      # Points on the line
#   scale_x_reverse(limits = c(70, 0)) +       # Reverse the x-axis for time scale
#   labs(x = "Time (millions of years ago)",   # X-axis label
#        y = "Importance of Event",            # Y-axis label
#        title = "Evolutionary Events Importance Over Time") +
#   theme_minimal()
