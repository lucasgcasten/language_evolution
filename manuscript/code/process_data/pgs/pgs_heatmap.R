library(tidyverse)

#####
pc <- read_csv('/Dedicated/jmichaelson-wdata/lcasten/sli_wgs/PCA/PCA_results.merged.1000_genomes_EUR.SLI_WGS.qc.demo.csv')
names(pc)[1] <- 'IID'

#### get PGS files ####
files = list.files('/Dedicated/jmichaelson-wdata/lcasten/sli_wgs/prs/LDPred2-inf-v2', 
                   pattern = 'PGS.tsv', recursive = TRUE)
length(files)

pgs_list <- list()

for (f in files) {
  ##
  i = which(f == files)
  if (i %% 1 == 0) {
    print(str_c(i, '/', length(files)))
  }

  ##
  filename <- str_c('/Dedicated/jmichaelson-wdata/lcasten/sli_wgs/prs/LDPred2-inf-v2/', f)
  pgs <- read_table(file = filename,
                  show_col_types = FALSE)

  tmp <- pgs %>%
    inner_join(pc)
  ##
  kg <- tmp %>%
    filter(population == '1000 Genomes')
  dg <- tmp %>%
    filter(population != '1000 Genomes')
  ##
  mod <- lm(PGS ~ pc1 + pc2 + pc3 + pc4 + pc5, 
            data = kg)
  ##
  preds_kg <- predict(mod, newdata = kg)
  preds_all <- predict(mod, newdata = dg)
  ##
  resids_kg = kg$PGS - preds_kg 
  resids_dg <- dg$PGS - preds_all
  ##
  kg_mean = mean(resids_kg)
  kg_sd = sd(resids_kg)
  z_kg <- (resids_kg - kg_mean) / kg_sd
  z = (resids_dg - kg_mean) / kg_sd
  dg$pgs_pc_corrected = z
  kg$pgs_pc_corrected = z_kg
  
  tmp2 <- bind_rows(dg, kg) %>% 
    select(IID, pgs_name, PGS, pgs_pc_corrected)
  pgs_list[[i]] <- tmp2
}

##
pgs_all <- bind_rows(pgs_list)

unique(pgs_all$pgs_name)

##
pgs_all %>% 
  write_csv('/Dedicated/jmichaelson-wdata/lcasten/sli_wgs/prs/LDPred2-inf-v2/gathered_normed_pgs.long.csv')

##
unique(pgs_all$pgs_name)

##
fc <- read_csv('/wdata/common/SLI_WGS/public/phenotype/factors/pheno_factors_resid.csv')
fc <- fc %>%
  pivot_longer(cols = -1)
sample_map = read_csv('/wdata/common/SLI_WGS/public/phenotype/factors/pheno_factors.csv') %>%
  rename(IID = id) %>%
  select(1:2) %>%
  distinct()
fc <- sample_map %>%
  inner_join(fc) %>%
  rename(id_phenotype = IID, IID = sample)


##
res <- pgs_all %>% 
  inner_join(fc) %>% 
  group_by(pgs_name, name) %>% 
  do(res = broom::tidy(cor.test(.$pgs_pc_corrected, .$value, method = 'p'))) %>% 
  unnest(res) %>% 
  arrange(p.value)

res %>% 
  filter(p.value < 0.05) %>% 
  distinct(pgs_name) %>% 
  as.data.frame()

res %>% 
  group_by(name) %>% 
  slice_head(n = 3) %>% 
  select(1:2) %>% 
  as.data.frame

unique(res$pgs_name)

res %>% 
  filter(str_detect(pgs_name, pattern = 'height|BIG5|_volume$|_skills$')) %>% 
  filter(str_detect(pgs_name, pattern = 'cerebe', negate = TRUE)) %>% 
  filter(str_detect(pgs_name, pattern = 'height')) %>% 
  select(1:3, p.value) 

res %>% 
  filter(str_detect(pgs_name, pattern = 'addiction|neurodev|epilepsy$')) %>% 
  select(1:3, p.value) 

res %>% 
  filter(str_detect(pgs_name, pattern = 'PGC|SSGAC|cognitive|functional|income|GPC2|pitch|Hand|hand|^age_|pgc|2019_Hill_householdIncome|IEU_Neale_TownsendDeprivationIndex_ukb-a-44|_bmi|executive|ENIGMA|insomnia|antisocial|Edu|ssgac|pgc|cog_gFactor-UKB-2020|height')) %>% 
  filter(name == 'Factor1') %>% 
  select(1,3, p.value) 


res_subset <- res %>% 
  filter(str_detect(pgs_name, pattern = 'PGC|SSGAC|cognitive|functional|income|GPC2|pitch|Hand|hand|^age_|pgc|2019_Hill_householdIncome|IEU_Neale_TownsendDeprivationIndex_ukb-a-44|_bmi|executive|ENIGMA|insomnia|antisocial|Edu|ssgac|pgc|cog_gFactor-UKB-2020|BIG5|_volume|height_Yengo|Cogntive_skills|neurodev|addiction|epilepsy$'))
unique(res_subset$pgs_name)


source('/wdata/lcasten/functions/simple_ggplot_theme_presentation.R')

pgs_nm <- c('2018_ssgac_cognitive_performance', '2021_pgc_schizophrenia', '2017-Warrier-MolecularPsychiatry-cognitive_empathy', 'hatoum2023-executive_functioning', 'functionalConnectivity_Yeo7networkParcellation_Ventral.Attention', 'age_first_birth',
            '2018_giantUKBB_meta_bmi', 'ADHD-PGC-2019', 'IEU_Neale_TownsendDeprivationIndex_ukb-a-44', '2019_Hill_householdIncome', 'ENIGMA_globalSurfaceArea', 'broad_antisocial_behavior',
            'CTGlab_2019_insomnia', 'functionalConnectivity_Yeo7networkParcellation_Frontoparietal', 'alcohol_dependency-PGC-2018', 'depression-PGC-2019', 'functionalConnectivity_Yeo7networkParcellation_Default', 'PTSD-PGC-2019',
            'functionalConnectivity_Yeo7networkParcellation_Limbic', 'age_first_sexual_intercourse', 'GPC2_NEUROTICISM', '2021_PGC_bipolar_I_II', 'functionalConnectivity_Yeo7networkParcellation_Dorsal.Attention', 'fitzgerald2022_cognitive_resilience',
            'functionalConnectivity_Yeo7networkParcellation_Visual', 'vocal_pitch_median_f0', 'cannabis_use_disorder-PGC-2020', '2019_pgc_autism', 'functionalConnectivity_Yeo7networkParcellation_global', 'ENIGMA_globalThickness',
            '2019-PGC-anorexia', 'PGC2_exclUKBB_excl23andMe_alzheimers', 'functionalConnectivity_Yeo7networkParcellation_Somatomotor', 'GPC2_EXTRAVERSION', '2022_ssgac_educational_attainment', 'UKBB_fastGWA_leftHanded',
            'UKBB2019-non_heterosexual_behavior', 'cog_gFactor-UKB-2020', 'extra_BIG5', 'consc_BIG5', 'Neuro_BIG5', 'open_BIG5', 'agree_BIG5', 'accumbens_volume', 'putamen_volume', 'thalamus_volume', 'brainstem_volume', 'pallidum_volume', 'caudate_volume', 'amygdala_volume', 'subcortical_volume', 'cerebellar_volume', 'height_Yengo', 'Cogntive_skills', 'NCogntive_skills', 'rare_neurodevelopmental_condition', 'epilepsy', 'addiction_Hatoum2023'
            )
clean_nm <- c('Cognitive performance', 'Schizophrenia', 'Empathy', 'Executive functioning', 'FC ventral attention', 'Age at first birth',
              'BMI', 'ADHD', 'Townsend deprivation index', 'Income', 'Global cortical surface area', 'Antisocial behavior', 
              'Insomnia', 'FC frontoparietal', 'Alcohol dependency', 'Depression', 'FC default mode', 'PTSD',
              'FC limbic', 'Age at first sexual intercourse', 'Neuroticism', 'Bipolar', 'FC dorsal attention', 'Cognitive resilience',
              'FC visual', 'Vocal pitch', 'Cannabis use disorder', 'Autism', 'FC global', 'Global cortical thickness',
              'Anorexia', 'Alzheimers', 'FC somatomotor', 'Extraversion', 'Educational attainment', 'Left handed',
              'Non-heterosexual behavior', 'g Factor', 'Extraversion', 'Conscientiousness', 'Neuroticism', 'Openness', 'Agreeableness', 'Nucleus accumbens volume', 'Putamen volume', 'Thalamus volume', 'Brainstem volume', 'Pallidum volume', 'Caudate volume', 'Amygdala volume', 'Subcortical volume', 'Cerebellar volume', 'Height', 'Cognitive education skills', 'Non-cognitive education skills',
              'Neurodevelopmental condition', 'Epilepsy', 'Addiction'
              )
clean_names <- data.frame(pgs_name = pgs_nm, clean_name = clean_nm)

p_dat <- res_subset %>% 
  filter(str_detect(pgs_name, pattern = 'GPC2', negate = TRUE)) %>%
  mutate(type = case_when(str_detect(pgs_name, pattern = 'edu|cognitiv|exec|empath|Edu|gFactor|Cognitive') ~ 'Cognitive',
                          str_detect(pgs_name, pattern = 'Income|Townsend|GPC|antisoc|age|BIG5') ~ 'Behavioral/SES',
                          str_detect(pgs_name, pattern = 'autism|bipolar|anorex|alzheim|cannabis|PTSD|depression|alcohol|insomni|ADHD|schizo') ~ 'Neuro/Psychiatric',
                          str_detect(pgs_name, pattern = 'functional|ENIGMA|volume') ~ 'Brain MRI',
                          str_detect(pgs_name, pattern = 'bmi|vocal|Hand|Height') ~ 'Misc.')) %>% 
  inner_join(clean_names) %>%  
  filter(str_detect(pgs_name, pattern = 'resilien', negate = TRUE)) %>%
  group_by(name) %>%
  mutate(fdr = p.adjust(p.value, method = 'fdr')) %>% 
  ungroup() %>%
  mutate(sig = case_when(fdr <= 0.05 ~ '**',
                         fdr > 0.05 & p.value < 0.05 ~ '*',
                         TRUE ~ '')) %>% 
  arrange(desc(clean_name)) %>% 
  mutate(clean_name = factor(clean_name, levels = unique(clean_name))) %>% 
  mutate(type = factor(type, levels = c('Cognitive', 'Neuro/Psychiatric', 'Behavioral/SES', 'Brain MRI', 'Misc.')))
p <- p_dat %>%
  ggplot(aes(x = str_remove_all(name, pattern = 'actor'), y = clean_name, fill = estimate)) +
  geom_tile() +
  geom_text(aes(label = sig), size = 5) +
  facet_grid(rows = vars(type), scales = 'free_y', space = 'free_y') +
  scale_fill_gradient2(low = 'dodgerblue', mid = 'white', high = 'chocolate1', midpoint = 0) +
  xlab('Factor') +
  ylab('Polygenic score') +
  labs(fill = 'Correlation:') +
  theme(legend.key.size = unit(1, "cm"),
        panel.border = element_rect(color = "black", fill = NA, size = 1))

p %>%
  ggsave(filename = '/wdata/lcasten/sli_wgs/prs/pgs_correlation_heatmap.png', device = 'png', dpi = 300, units = 'in', width = 10, height = 17)

p_dat %>% 
    write_csv('/wdata/lcasten/sli_wgs/prs/pgs_correlation_heatmap_data.csv')