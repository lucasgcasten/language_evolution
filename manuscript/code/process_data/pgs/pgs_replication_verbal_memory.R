library(tidyverse)

##
setwd('/wdata/lcasten/screener_paper')

###############################################3
###############################################3
################################################
fs = read_rds(file = 'data/factors/merged_factors_whisper_raw_fa2.rds')
fs
#### use X factors ####
f_number = ncol(fs$scores)
# pheno_fact <- fs$fitted_imputed
fs$scores[1:5,]
pheno_fact <- fs$scores %>%
  as.data.frame() %>%
  rownames_to_column(var = 'screener_id') %>%
  as_tibble()
names(pheno_fact)[-1] <- str_replace_all(names(pheno_fact)[-1], pattern = '^MR', replacement = 'Factor')
pheno_fact
spark = read_csv('/wdata/lcasten/spark/research_match/language_merged_results/data/ID_mapping_and_demographics.csv')
spark_demo = spark %>%
  select(screener_id, subject_sp_id = spid, age, sex, diagnosis) %>%
  mutate(sex_male = case_when(sex == 'female' ~ 0,
                              sex == 'male' ~ 1,
                              TRUE ~ NA_real_)
         ) %>%
  mutate(age_months = age * 12) %>%
  select(screener_id, subject_sp_id, age_months, sex_male, diagnosis)

spark_demo_repeat = spark_demo %>%
  mutate(screener_id = str_c(screener_id, 'r'))

spark_demo <- bind_rows(spark_demo, spark_demo_repeat)

#### merge demographics ####
# demo = bind_rows(spark_demo)
demo = spark_demo
demo

#### merge demo with phenotype info ####
f7 = demo %>%
  inner_join(pheno_fact) %>%
  filter(str_detect(screener_id, pattern = 'r$', negate = T))
table(f7$sex_male)
table(f7$diagnosis)
f7
# f7 %>%
#   relocate(Factor1, Factor2, Factor3, Factor4, Factor5,.after = diagnosis) %>%
#   write_csv('/wdata/lcasten/screener_paper/data/factors/v2/factor_scores.csv')
### pgs 
eur <- read_tsv('/wdata/trthomas/array/merged_2022_ABCD_iWES1_WGS_2-4/master.tsv') %>%
  filter(cluster_allPCs == 1)
pgs <- read_csv('/wdata/lcasten/spark/prs/HapMap3_plus/LDPred2-inf-full/gathered_LDPred2-inf_pc_corrected_long.csv') %>%
  filter(IID %in% f7$subject_sp_id) %>%
  filter(IID %in% eur$IID)
pgs <- pgs %>%
  rename(subject_sp_id = IID)
pgs_names <- distinct(pgs, clean_name, pgs_name)
res_pgs <- f7  %>%
  pivot_longer(cols = -c(1:5)) %>%
  inner_join(pgs) %>%
  # filter(str_detect(pgs_name, pattern = 'PGC|pgc|ssgac|SSGAC')) %>%
  group_by(name, pgs_name) %>%
  do(res = broom::tidy(cor.test(.$value, .$pgs_pc_corrected))) %>%
  unnest(res) %>%
  arrange(p.value)


res_f2_pgs <- res_pgs %>%
  filter(name == 'Factor2') %>%
  mutate(padj = p.adjust(p.value, method = 'fdr')) %>%
    mutate(pval_lab = str_c('R = ', formatC(estimate, digits = 2), ', p-val = ', formatC(p.value, digits = 2), ', N = ', parameter + 2)) %>%
    ungroup() %>%
    mutate(pgs_name = factor(pgs_name, levels = rev(unique(pgs_name))))
spark_sig <- res_f2_pgs %>%
  filter(p.value < 0.05) %>%
  select(1:5) %>%
  as.data.frame()


###########################
## repeat for SLI
fc <- read_csv('/wdata/common/SLI_WGS/public/phenotype/factors/pheno_factors_resid.csv')
fc <- fc %>%
  pivot_longer(cols = -1)
sample_map = read_csv('/wdata/common/SLI_WGS/public/phenotype/factors/pheno_factors.csv') %>%
  rename(IID = id) %>%
  select(1:2) %>%
  distinct()
fc <- sample_map %>%
  inner_join(fc) %>%
  rename(id_phenotype = IID, id_genotype = sample)

pgs2 <-read_csv('/sdata/SLI_WGS/prs/LDPred2-inf/LDPred2_gathered_pgs_pc_corrected_long.csv')
# pgs = pgs %>%
#   mutate(IID = str_split(IID, pattern = '-', simplify = TRUE)[,1])
cor_res = fc %>%
  inner_join(pgs2)  %>%
  group_by(pgs_name, name) %>%
  do(res = broom::tidy(cor.test(.$pgs_pc_corrected, .$value))) %>%
  unnest(res) %>%
  arrange(p.value)
clean_names = pgs2 %>%
  select(pgs_name, pgs_clean_name) %>%
  distinct()

sli_sig <- cor_res %>%
  inner_join(clean_names) %>%
  relocate(pgs_clean_name) %>%
  filter(name == 'Factor1') %>%
  filter(p.value < 0.05) %>%
    mutate(pval_lab = str_c('R = ', formatC(estimate, digits = 2), ', p-val = ', formatC(p.value, digits = 2), ', N = ', parameter + 2))
kp_pgs <- intersect(spark_sig$pgs_name, sli_sig$pgs_name)

sli_order <- sli_sig %>%
    filter(pgs_name %in% kp_pgs) %>%
    filter(pgs_name %in% c('PGI-CP1_single_gwide_sumstats.txt', '2022_ssgac_educational_attainment', 'wang2022-physical_activity_during_leisure_time', 'hatoum2023-executive_functioning'))

##############################
## make figures
##### SPARK
source('/wdata/lcasten/functions/simple_ggplot_theme.R')
p_pgs <- f7  %>%
  pivot_longer(cols = -c(1:5)) %>%
  ungroup() %>%
  inner_join(pgs) %>%
  inner_join(filter(res_f2_pgs, pgs_name %in% c('PGI-CP1_single_gwide_sumstats.txt', '2022_ssgac_educational_attainment', 'wang2022-physical_activity_during_leisure_time', 'hatoum2023-executive_functioning'))) %>%
    arrange(p.value) %>%
    mutate(sample = 'replication (SPARK)') %>%
    inner_join(pgs_names) %>%
        mutate(clean_name = factor(clean_name, levels = sli_order$pgs_clean_name)) %>%
  ggplot(aes(x = pgs_pc_corrected, y = value)) +
  geom_point(alpha = 0.4) +
  geom_smooth(method = 'lm') +
  facet_wrap(~ clean_name, scales = 'free_x') +
  geom_text(aes(x = 0, y = 3, label = pval_lab), check_overlap = T, size = 4.5) +
  xlab('Polygenic score') +
  ylab('Verbal memory factor score') +
  ggtitle('Replication cohort (SPARK)')

p_pgs %>%
    write_rds(file = '/wdata/lcasten/sli_wgs/prs/SPARK_pgs_factors2.rds')

## repeat for SLI
p_pgs_sli <- fc %>%
  inner_join(pgs2) %>%
  inner_join(filter(sli_sig, pgs_name %in% c('PGI-CP1_single_gwide_sumstats.txt', '2022_ssgac_educational_attainment', 'wang2022-physical_activity_during_leisure_time', 'hatoum2023-executive_functioning'))) %>%
    arrange(p.value) %>%
    mutate(sample = 'epiSLI WGS') %>%
    inner_join(pgs_names) %>%
        mutate(clean_name = factor(clean_name, levels = sli_order$pgs_clean_name)) %>%
  ggplot(aes(x = pgs_pc_corrected, y = value)) +
  geom_point(alpha = 0.6) +
  geom_smooth(method = 'lm') +
  facet_wrap(~ clean_name, scales = 'free_x') +
  geom_text(aes(x = 0, y = 3, label = pval_lab), check_overlap = T, size = 4.5) +
  xlab('Polygenic score') +
  ylab('Verbal memory factor score') +
  ggtitle('epiSLI WGS')

p_pgs_sli %>%
    write_rds(file = '/wdata/lcasten/sli_wgs/prs/SLI_pgs_factor1.rds')


## merge plots
library(patchwork)
p_pgs_merged <- p_pgs_sli + p_pgs + plot_annotation(tag_levels = 'A') & theme(plot.tag = element_text(face = 'bold', size = 22))
ggsave(plot = p_pgs_merged, filename = '/wdata/lcasten/sli_wgs/prs/pgs_replication.png', 
       dpi = 300, units = 'in', device = 'png', width = 18, height = 10)
