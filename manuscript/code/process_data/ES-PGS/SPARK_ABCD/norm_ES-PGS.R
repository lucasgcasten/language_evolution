library(tidyverse)

#####
pc <- read_csv('/Dedicated/jmichaelson-wdata/lcasten/spark/prs/HapMap3_plus/PCA/raw_KING_pca_results.csv')
names(pc)[1] <- 'IID'

##
unrel <- read_table("/Dedicated/jmichaelson-wdata/lcasten/spark/prs/pathway/SPARK_ABCD_unrelated_europeans_for_LD.txt", col_names = FALSE)

#### get PGS files ####
files = list.files('/Dedicated/jmichaelson-wdata/lcasten/sli_wgs/prs/pathway_prs', 
                   pattern = 'all_score', recursive = TRUE)
files <- files[str_detect(files, pattern = 'SPARK_ABCD')]
files <- files[str_detect(files, pattern = 'human_evolution|custom_gene_sets')]
files <- files[str_detect(files, pattern = 'cogPerf')]
files <- files[str_detect(files, pattern = 'complement')]


pgs_list <- list()
prsice_list = list()
for (f in files) {
  ##
  i = which(f == files)
  if (i %% 1 == 0) {
    print(str_c(i, '/', length(files)))
  }
  trait = str_split(f, pattern = '[.]', simplify = TRUE)[,1]
  gs <- str_split(f, pattern = '[.]', simplify = TRUE)[,2]

  ##
  filename <- str_c('/Dedicated/jmichaelson-wdata/lcasten/sli_wgs/prs/pathway_prs/', f)
  pgs <- read_table(file = filename,
                  show_col_types = FALSE)
  prsice <- read_table(file = str_replace_all(filename, pattern = 'all_score', replacement = 'prsice'),
                  show_col_types = FALSE) %>% 
            mutate(trait = trait, gs = gs) %>% 
            relocate(trait, gs)
  prsice_list[[i]] <- prsice

  names(pgs)[-c(1:2)] <- str_c(trait, '.', gs, '.', names(pgs)[-c(1:2)])
  names(pgs) <- str_replace_all(names(pgs), pattern = '-', replacement = '_')
  tmp <- pgs %>%
    inner_join(pc)
  ##
  kg <- tmp %>%
    filter(IID %in% unrel$X2)
  dg <- tmp %>%
    filter(! IID %in% unrel$X2)
  pgs_list_tmp <- list()
  ##
  for(j in 3:ncol(pgs)){
    nm <- colnames(kg)[j]
    kg$PGS <- unname(unlist(kg[,colnames(kg)[j]]))
    dg$PGS <- unname(unlist(dg[,colnames(kg)[j]]))
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
    
    tmp2 <- bind_rows(dg, kg)
    tmp2 <- tmp2 %>%
      mutate(pgs_name = nm) %>%
      select(IID, pgs_name, PGS, pgs_pc_corrected)
    pgs_list_tmp[[j]] <- tmp2
  }
  pgs_list[[i]] <- bind_rows(pgs_list_tmp)
}

########################################################
##
##
pgs_all <- bind_rows(pgs_list)
prsice_all <- bind_rows(prsice_list)
prsice_sets <- prsice_all %>% 
  filter(Set != 'Base') %>% 
  mutate(thr = case_when(Threshold == 5e-08 ~ '5e_08'))
range(prsice_sets$Num_SNP)

## try and get excess burden (amount of burden in set not explained by base/genome-wide score)
base_all <- pgs_all %>% 
  filter(str_detect(pgs_name, pattern = 'Base')) %>% 
  mutate(trait = str_split(pgs_name, pattern = '[.]', simplify = TRUE)[,1],
         gs = str_split(pgs_name, pattern = '[.]', simplify = TRUE)[,2])

base_wide <- base_all %>% 
  pivot_wider(id_cols = 1, names_from = pgs_name, values_from = pgs_pc_corrected)
gs_wide <- pgs_all %>% 
  filter(str_detect(pgs_name, pattern = 'Base', negate = TRUE)) %>% 
  mutate(trait = str_split(pgs_name, pattern = '[.]', simplify = TRUE)[,1],
         gs = str_split(pgs_name, pattern = '[.]', simplify = TRUE)[,2]) %>% 
  pivot_wider(id_cols = 1, names_from = pgs_name, values_from = pgs_pc_corrected)


excess_burden_list = list()
for(i in 2:ncol(gs_wide)){
  cat(sprintf('\n\n\n\n'))
  message('====================================')
  message(i, '/', ncol(gs_wide))
  tmp <- gs_wide[,c(1,i)]
  trait <- str_split(names(tmp)[2], pattern = '[.]', simplify = TRUE)[,1] 
  gs <- str_split(names(tmp)[2], pattern = '[.]', simplify = TRUE)[,2] 

  thr <- case_when(str_detect(names(tmp)[2], pattern = '_5e_08$') ~ '_5e_08',
                   str_detect(names(tmp)[2], pattern = '_0.0005$') ~ '_0.0005',
                   str_detect(names(tmp)[2], pattern = '_0.05$') ~ '_0.05',
                   str_detect(names(tmp)[2], pattern = '_0.2$') ~ '_0.2',
                   str_detect(names(tmp)[2], pattern = '_1$') ~ '_1')
  base_name = str_c(trait, '.', gs, '.Base', thr)
  base_tmp <- base_wide[,c("IID", base_name)]

  nm = names(tmp)[2]
  df <- tmp %>% 
    inner_join(base_tmp)
  names(df)[2:3] <- c('y', 'x')
  kg <- df %>% 
    filter(IID %in% unrel$X2)
  dg <- df %>% 
    filter(! IID %in% unrel$X2)
  mod <- lm(y ~ x, 
              data = kg)
    ##
    preds_kg <- predict(mod, newdata = kg)
    preds_all <- predict(mod, newdata = dg)
    ##
    resids_kg = kg$y - preds_kg 
    resids_dg <- dg$y - preds_all
    ##
    kg_mean = mean(resids_kg)
    kg_sd = sd(resids_kg)
    z_kg <- (resids_kg - kg_mean) / kg_sd
    z = (resids_dg - kg_mean) / kg_sd
    dg$pgs_genome_wide_corrected = z
    kg$pgs_genome_wide_corrected = z_kg
    
    tmp2 <- bind_rows(dg, kg)
    tmp2 <- tmp2 %>%
      mutate(pgs_name = nm) %>%
      select(IID, pgs_name, pgs_genome_wide = x, pgs_pathway_corrected_for_genome_wide_burden = pgs_genome_wide_corrected)
    excess_burden_list[[i]] <- tmp2
}
excess_burden <- bind_rows(excess_burden_list)
unique(excess_burden$pgs_name)

length(unique(pgs_all$pgs_name))
unique(pgs_all$pgs_name)
length(unique(pgs_all$IID))

kp = sample(unique(pgs_all$pgs_name), size = 10)
kp
pgs_all  %>% 
  filter(pgs_name %in% kp) %>% 
  inner_join(pc) %>%
  ggplot(aes(x = pgs_pc_corrected)) +
  geom_histogram() +
  facet_wrap(~ pgs_name)

pgs_all  %>% 
  filter(str_detect(pgs_name, 'cogPerf') & str_detect(pgs_name, 'HAQER')) %>% 
  inner_join(pc) %>%
  ggplot(aes(x = pgs_pc_corrected)) +
  geom_histogram() +
  facet_wrap(~ pgs_name)

##
f <- '/Dedicated/jmichaelson-wdata/lcasten/spark_abcd_array_merge/topmed_hg19/autosomes.qc.fam'
fam <- read_table(file = f, col_names = FALSE)
fam = fam[,1:2]
names(fam) <- c('FID', 'IID')

dat <- fam %>%
  inner_join(pgs_all) %>%
  rename(pgs_raw = PGS) %>%
  inner_join(excess_burden) %>%
  relocate(pgs_genome_wide_baseline = pgs_genome_wide, .after = pgs_name) %>%
  inner_join(select(pc, IID)) %>%
  mutate(cohort = 'SPARK_ABCD') %>%
  relocate(cohort, .after = IID) %>% 
  mutate(pgs_name = str_remove_all(pgs_name, pattern = 'SPARK_ABCD.'))

## make df w/ cols for complement + anno scores
datc <- dat %>% 
  filter(str_detect(pgs_name, 'complement_')) %>% 
  mutate(pgs_name = str_remove_all(pgs_name, pattern = 'complement_')) 
names(datc)[6:8] <- str_c(names(datc)[6:8], '_complement')
names(datc)

data <- dat %>% 
  filter(str_detect(pgs_name, 'complement_', negate = TRUE)) %>% 
  inner_join(datc)

unique(datc$pgs_name)[1:10]
unique(data$pgs_name)[1:10]

data %>%
  write_csv('/Dedicated/jmichaelson-wdata/lcasten/sli_wgs/prs/pathway_prs/SPARK_ABCD.gathered_pgs_pc_corrected_long_full_data.cogPerf.complement.csv')
prsice_sets %>% 
  write_csv('/Dedicated/jmichaelson-wdata/lcasten/sli_wgs/prs/pathway_prs/SPARK_ABCD.prsice_info.cogPerf.complement.csv')
prsice_sets %>% 
    filter(Threshold == 1) %>% 
    head(n = 15)
# read_csv('/Dedicated/jmichaelson-wdata/lcasten/sli_wgs/prs/pathway_prs/gathered_pgs_pc_corrected_long_full_data.cogPerf.csv')

