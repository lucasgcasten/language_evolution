library(tidyverse)

## ===================================================
## subset modern PGS data to cog SNPs used by PRSet
## ===================================================
bim <- read_table('/Dedicated/jmichaelson-wdata/lcasten/sli_wgs/prs/1000_genomes_EUR_merge/merged.1000_genomes_EUR.SLI_WGS.qc.bim', col_names = FALSE)
names(bim) <- c('chromosome', 'snp', 'cm', 'start', 'minor', 'major')

##
cog <- read_table('/wdata/lcasten/sli_wgs/prs/pathway_prs/cogPerf.human_evolution_complement.snp')
cog

##
bim %>% 
    filter(snp %in% cog$SNP) %>% 
    select(snp) %>% 
    distinct() %>% 
    write_tsv('/wdata/lcasten/sli_wgs/prs/pathway_prs/cog_snp_PRSet.rsid', col_names = FALSE)

##
cmd <- "/wdata/lcasten/tools/plink --bfile /wdata/lcasten/sli_wgs/prs/1000_genomes_EUR_merge/merged.1000_genomes_EUR.SLI_WGS.qc --extract /wdata/lcasten/sli_wgs/prs/pathway_prs/cog_snp_PRSet.rsid --keep-allele-order --make-bed --out /wdata/lcasten/sli_wgs/neanderthal_analysis/data/EpiSLI_1000Genomes_Eur_cog_snp"
print(cmd)
system(cmd)


## ===================================================
## read in PLINK formatted cog SNPs
## ===================================================
ref <- BEDMatrix::BEDMatrix("/wdata/lcasten/sli_wgs/neanderthal_analysis/data/EpiSLI_1000Genomes_Eur_cog_snp.bed")
nean <- BEDMatrix::BEDMatrix("/wdata/lcasten/tools/ref_data/archaic_hominini/hg19_neanderthal_filtered_vcf/vcf/merged_all_chromosomes.cog_perf_independent_snp.bed")

ref_bim <- read_table('/wdata/lcasten/sli_wgs/neanderthal_analysis/data/EpiSLI_1000Genomes_Eur_cog_snp.bim', col_names = FALSE)
nean_bim <- read_table('/wdata/lcasten/tools/ref_data/archaic_hominini/hg19_neanderthal_filtered_vcf/vcf/merged_all_chromosomes.cog_perf_independent_snp.bim', col_names = FALSE)

## find overlapping SNPs with matching rsid, positions, and REF alleles
kp_snp <- ref_bim %>% 
    inner_join(select(nean_bim, -X5)) %>% 
    select(X2) %>% 
    unlist() %>% 
    unname()
length(kp_snp)

## convert 1000 genomes data to tibble
ref_df <- ref %>% 
    as.matrix() %>% 
    as.data.frame() %>% 
    rownames_to_column('IID') %>% 
    as_tibble()
## clean rsid colnames
colnames(ref_df) <- str_split(colnames(ref_df), pattern = '_', simplify = TRUE)[,1]


## ===================================================
## reformat genotype data and find overlapping SNPs
## ===================================================
nean_df <- nean %>% 
    as.matrix() %>% 
    as.data.frame()
## remove duplicate SNPs
nean_df <- nean_df[, !duplicated(colnames(nean_df))]
## conver to tibble
nean_df <- nean_df %>% 
    rownames_to_column('IID') %>% 
    as_tibble()
colnames(nean_df) <- str_split(colnames(nean_df), pattern = '_', simplify = TRUE)[,1]

## get overlapping SNPs
keep_snp <- intersect(colnames(nean_df), colnames(ref_df))
length(keep_snp)
nean_df <- nean_df %>% 
    select(any_of(keep_snp))
ref_df <- ref_df %>% 
    select(any_of(keep_snp))

## merge 1000 Genomes w/ Neanderthals
wd <- ref_df %>% 
    bind_rows(nean_df)


## ======================================================
## "impute" missing genotypes with mean value (i.e., MAF)
## ======================================================
wd_t <- data.table::transpose(wd[,-c(1)])

colnames(wd_t) <- str_split(wd$IID, pattern = '_', simplify = TRUE)[,1]
rownames(wd_t) <- colnames(wd)[-1]

wd_imp <- wd_t %>% 
    rownames_to_column(var = 'snp') %>% 
    as_tibble()
maf <- rowMeans(wd_imp[,-1], na.rm = TRUE)
wd_imp$maf <- maf # / 2

## double check nothing weird happened in MAF calculation
mean(wd_imp$maf)
min(wd_imp$maf)
max(wd_imp$maf)

##
samples <- colnames(wd_imp)[str_detect(colnames(wd_imp), pattern = 'HG|NA|sample|Vind|Alta|Denisov|Chagyr')]

## mean impute missing values
for(s in samples) {
    ind <- which(is.na(wd_imp[,s]))
    x <- wd_imp[,s] %>% 
        unlist() %>% 
        unname()
    tmp_geno <- ifelse(is.na(x), maf, x)
    wd_imp[,s] <- tmp_geno
}

## ===================================================
## compute background PGS
## ===================================================
## independent SNPs used in ES-PGS for EpiSLI
ind_snp <- read_table('/wdata/lcasten/sli_wgs/prs/pathway_prs/cogPerf.human_evolution_complement.snp')
ind_snp <- ind_snp %>% 
    select(snp = SNP, complement_HAQER, HAQER, HAR) %>% 
    filter(complement_HAQER == 1 | HAQER == 1 | HAR == 1)
## CP sumstats
cp <- read_table('/sdata/gwas_summary_stats/ssgac/2018/raw/GWAS_CP_all.txt.gz')
sumstats <- cp %>% 
    select(snp = MarkerName, CHR, POS, A1, A2, BETA = Beta, SE, p = Pval) %>% 
    filter(snp %in% kp_snp)
## merge sumstats w/ CP betas and flip sign where necessary
wd_df <- bim %>% 
    inner_join(ind_snp) %>%
    inner_join(sumstats) %>% 
    mutate(BETA = case_when(A1 == major & A2 == minor ~ -1 * BETA,
                            A1 == minor & A2 == major ~ BETA,
                            TRUE ~ NA)) %>% 
    inner_join(wd_imp)

## get betas and genotypes for background variants of interest
betas_gw <- wd_df$BETA[wd_df$complement_HAQER == 1]
genos_gw <- wd_df %>%
    filter(complement_HAQER == 1) %>%
    select(matches('^NA|^HG|^sample|^Deni|Vind|Altai|Chagyr')) %>% 
    as.data.frame() 
gw_pgs <- data.frame(IID = colnames(genos_gw)) %>% 
    mutate(pgs = NA) %>% 
    as_tibble()

## compute PGS
for(i in 1:nrow(gw_pgs)) {
    if(i %% 25 == 0) {
        message(i, '/', nrow(gw_pgs))
    }
    s = gw_pgs$IID[i]
    tpgs <- sum(genos_gw[,s] * betas_gw)
    gw_pgs$pgs[i] <- tpgs
}

## check that the manual PGS is the same as what we get from PRSet in 1000 Genomes + EpiSLI
prset <- read_table('/wdata/lcasten/sli_wgs/prs/pathway_prs/cogPerf.human_evolution_complement.all_score')

prset %>% 
    select(IID, complement_HAQER_1) %>% 
    inner_join(gw_pgs) %>% 
    select(-IID) %>% 
    cor()

gw_pgs2 <- gw_pgs %>% 
    # inner_join(pc_dat) %>% 
    mutate(pgs = scale(pgs)[,1])

## ==============================================
## compute ES-PGS for HAQERs
## ==============================================
## get betas
betas_haq <- wd_df$BETA[wd_df$HAQER == 1]
## get genotypes
genos_haq <- wd_df %>%
    filter(HAQER == 1) %>%
    select(matches('^NA|^HG|^sample|^Deni|Vind|Altai|Chagyr')) %>% 
    as.data.frame()

## PGS calculation
haq_pgs <- data.frame(IID = colnames(genos_haq)) %>% 
    mutate(pgs = NA) %>% 
    as_tibble()

for(i in 1:nrow(haq_pgs)) {
    if(i %% 25 == 0) {
        message(i, '/', nrow(haq_pgs))
    }
    s = haq_pgs$IID[i]
    tpgs <- sum(genos_haq[,s] * betas_haq)
    haq_pgs$pgs[i] <- tpgs
}

haq_pgs$pgs <- scale(haq_pgs$pgs)[,1]

## check that the manual PGS is the same as what we get from PRSet in 1000 Genomes + EpiSLI
prset <- read_table('/wdata/lcasten/sli_wgs/prs/pathway_prs/cogPerf.human_evolution_complement.all_score')

prset %>% 
    select(IID, HAQER_1) %>% 
    inner_join(haq_pgs) %>% 
    select(-IID) %>% 
    cor()


## ============================================================================================
## compute ES-PGS for HAQERs excluding alleles w/ any missingness in archaic humans
## ============================================================================================
## get betas
tail(names(wd_df))
wd_nean_mat <- as.data.frame(wd) %>% 
    filter(IID %in% c("AltaiNea_AltaiNea", "Chagyrskaya-Phalanx_Chagyrskaya-Phalanx", "DenisovaPinky_DenisovaPinky", "Vindija33.19_Vindija33.19"))
rownames(wd_nean_mat) <- wd_nean_mat$IID 
wd_nean_mat$IID <- NULL
wd_no_msg <- t(wd_nean_mat)
wd_no_msg <- wd_no_msg[rowSums(is.na(wd_no_msg)) == 0,]
haq_snp_no_msg_nean <- rownames(wd_no_msg)  
haq_rsid_no_msg <- wd_df$snp[wd_df$HAQER == 1 & wd_df$snp %in% haq_snp_no_msg_nean]
betas_haq_no_msg <- wd_df$BETA[wd_df$HAQER == 1 & wd_df$snp %in% haq_snp_no_msg_nean]

## get genotypes
genos_haq_no_msg <- wd_df %>%
    filter(snp %in% haq_rsid_no_msg) %>%
    select(matches('^NA|^HG|^sample|^Deni|Vind|Altai|Chagyr')) %>% 
    as.data.frame()

## PGS calculation
haq_pgs_no_msg <- data.frame(IID = colnames(genos_haq_no_msg)) %>% 
    mutate(pgs = NA) %>% 
    as_tibble()

for(i in 1:nrow(haq_pgs_no_msg)) {
    if(i %% 25 == 0) {
        message(i, '/', nrow(haq_pgs_no_msg))
    }
    s = haq_pgs_no_msg$IID[i]
    tpgs <- sum(genos_haq_no_msg[,s] * betas_haq_no_msg)
    haq_pgs_no_msg$pgs[i] <- tpgs
}

haq_pgs_no_msg$pgs <- scale(haq_pgs_no_msg$pgs)[,1]

## check that the manual PGS is the same as what we get from PRSet in 1000 Genomes + EpiSLI
prset <- read_table('/wdata/lcasten/sli_wgs/prs/pathway_prs/cogPerf.human_evolution_complement.all_score')

prset %>% 
    select(IID, HAQER_1) %>% 
    inner_join(haq_pgs_no_msg) %>% 
    select(-IID) %>% 
    cor()


## ================================
## gather PGS
## ================================
gw_pgs_clean <- gw_pgs %>%
    rename(cp_pgs.background = pgs) %>% 
    mutate(type = case_when(str_detect(IID, 'sample') ~ 'EpiSLI',
                            str_c(IID, '_', IID) %in% nean_df$IID ~ 'Neanderthals and Denisovan',
                            TRUE ~ '1000 Genomes')) 

haq_pgs_clean <- haq_pgs %>%
    rename(cp_pgs.HAQER = pgs) %>% 
    mutate(type = case_when(str_detect(IID, 'sample') ~ 'EpiSLI',
                            str_c(IID, '_', IID) %in% nean_df$IID ~ 'Neanderthals and Denisovan',
                            TRUE ~ '1000 Genomes'))

haq_pgs_no_msg_clean <- haq_pgs_no_msg %>%
    rename(cp_pgs.HAQER_no_missingness_archaic_subset = pgs) %>% 
    mutate(type = case_when(str_detect(IID, 'sample') ~ 'EpiSLI',
                            str_c(IID, '_', IID) %in% nean_df$IID ~ 'Neanderthals and Denisovan',
                            TRUE ~ '1000 Genomes'))  

pgs_wd <- gw_pgs_clean %>% 
    inner_join(haq_pgs_clean) %>%  
    inner_join(haq_pgs_no_msg_clean) %>%
    relocate(type, .after = IID)

## plot distributions across samples
pgs_wd %>% 
    pivot_longer(cols = matches('cp_pgs')) %>% 
    group_by(name) %>% 
    mutate(z = scale(value)[,1]) %>% 
    ggplot(aes(x = z, fill = type)) +
    geom_density(alpha = .6) +
    facet_wrap(~ name)

pgs_wd %>% 
    write_csv("/wdata/lcasten/sli_wgs/neanderthal_analysis/raw_ES-PGS_data.csv")
tail(pgs_wd)

#############################
## ===================================================
## compute genetic PCs - didn't use in figure since applying PCs 
## to 100,000 year old Neanderthals may not be valid (seems to make their PGS even more extreme)
## ===================================================
wd_tmp <- wd_t %>% 
    rownames_to_column(var = 'snp') %>% 
    as_tibble()

wd_tmp_haq <- wd_tmp %>% 
    filter(snp %in% wd_df$snp[wd_df$HAQER == 1])

mis_df <- data.frame(IID = colnames(wd_tmp), n_mis = colSums(is.na(wd_tmp))) %>% 
    as_tibble()
mis_df %>% 
    arrange(desc(n_mis))

mis_df_haq <- data.frame(IID = colnames(wd_tmp_haq), n_mis = colSums(is.na(wd_tmp_haq))) %>% 
    as_tibble()
mis_df_haq %>% 
    arrange(desc(n_mis))

pc_dat <- wd_imp %>% 
    select(matches('NA|HG')) %>% 
    as.data.frame() %>% 
    as.matrix()
pc_dat <- t(pc_dat)
colnames(pc_dat) <- wd_imp$snp
dim(pc_dat)

pc_dat_other <- wd_imp %>% 
    select(matches('sample|Vind|Denisov|Altai|Chagyr')) %>% 
    as.data.frame() %>% 
    as.matrix()
pc_dat_other <- t(pc_dat_other)
colnames(pc_dat_other) <- wd_imp$snp
pc <- prcomp(pc_dat, center = TRUE, .scale = TRUE)
summary(pc)
pc_proj = predict(pc, pc_dat_other)
str(pc_proj, 2)
pc_proj[1:5,1:5]

pc_df <- pc$x %>% 
    as.data.frame() %>% 
    rownames_to_column(var = 'IID') %>% 
    select(IID, str_c('PC', 1:5)) %>% 
    as_tibble() %>% 
    mutate(type = '1000 genomes')
pc_proj_df <- pc_proj %>% 
    as.data.frame() %>% 
    rownames_to_column(var = 'IID') %>% 
    select(IID, str_c('PC', 1:5)) %>% 
    as_tibble() %>% 
    mutate(type = ifelse(str_detect(IID, 'sample'), 'EpiSLI', 'Neanderthal'))
pc_dat <- bind_rows(pc_df,pc_proj_df)

pc_dat %>% 
    write_csv("/wdata/lcasten/sli_wgs/git/neanderthal_1000Genomes_PCs.csv")

## get resid based on 1000 genomes
td <- pgs_wd %>%
    inner_join(select(pc_dat, -type))
kg_pgs <- td %>% 
    filter(type == '1000 Genomes')

mod <- lm(cp_pgs.HAQER ~ PC1 + PC2 + PC3 + PC4 + PC5, data = kg_pgs)
prd <- predict(mod, td)
resid_tmp <- td$cp_pgs.HAQER - prd
mean_kg <- mean(resid_tmp[td$type == '1000 Genomes'])
sd_kg <- sd(resid_tmp[td$type == '1000 Genomes'])
resid_z <- (resid_tmp - mean_kg) / sd_kg
td$cp_pgs.HAQER_resid <- resid_z

### repeat for HAQER CP-PGS using only SNPs genotyped in all archaic samples (exclude variants missing in ANY neanderthal)
mod <- lm(cp_pgs.HAQER_no_missingness_archaic_subset ~ PC1 + PC2 + PC3 + PC4 + PC5, data = kg_pgs)
prd <- predict(mod, td)
resid_tmp <- td$cp_pgs.HAQER_no_missingness_archaic_subset - prd
mean_kg <- mean(resid_tmp[td$type == '1000 Genomes'])
sd_kg <- sd(resid_tmp[td$type == '1000 Genomes'])
resid_z <- (resid_tmp - mean_kg) / sd_kg
td$cp_pgs.HAQER_no_missingness_archaic_subset_resid <- resid_z
tail(td) %>% 
    select(-matches('PC')) %>% 
    select(matches('resid'))

### repeat for background
mod <- lm(cp_pgs.background ~ PC1 + PC2 + PC3 + PC4 + PC5, data = kg_pgs)
prd <- predict(mod, td)
resid_tmp <- td$cp_pgs.background - prd
mean_kg <- mean(resid_tmp[td$type == '1000 Genomes'])
sd_kg <- sd(resid_tmp[td$type == '1000 Genomes'])
resid_z <- (resid_tmp - mean_kg) / sd_kg
td$cp_pgs.background_resid <- resid_z

td %>% 
    tail() %>% 
    select(-matches('^PC|type'))

## save data
td %>% 
    rename(cp_pgs.HAQER_SNPs_with_no_archaic_missingness = cp_pgs.HAQER_no_missingness_archaic_subset, cp_pgs.HAQER_SNPs_with_no_archaic_missingness_resid = cp_pgs.HAQER_no_missingness_archaic_subset_resid) %>%
    write_csv("/wdata/lcasten/sli_wgs/neanderthal_analysis/ES-PGS_data.csv")


td %>% 
    rename(cp_pgs.HAQER_SNPs_with_no_archaic_missingness = cp_pgs.HAQER_no_missingness_archaic_subset, cp_pgs.HAQER_SNPs_with_no_archaic_missingness_resid = cp_pgs.HAQER_no_missingness_archaic_subset_resid) %>%
    write_csv("manuscript/supplemental_materials/archaic_human_data.csv")

##
table(pc_dat$type)

lab <- read_csv('/wdata/lcasten/sli_wgs/PCA/PCA_results.merged.1000_genomes_EUR.SLI_WGS.qc.demo.csv') %>% 
    rename(IID = sample.ID) %>% 
    select(1:8)
pc_dat <- pc_dat %>% 
    left_join(lab)
cor.test(pc_dat$PC1, pc_dat$pc1)
cor.test(pc_dat$PC4, pc_dat$pc4)
cor.test(pc_dat$PC5, pc_dat$pc5)
pc_dat <- pc_dat %>% 
    mutate(population_description = ifelse(is.na(population_description), 'Neanderthal', population_description))
table(pc_dat$population_description)

p <- pc_dat %>%
  filter(population_description != 'epiSLI' & population_description != 'Neanderthal') %>%
  ggplot(aes(x = PC1, y = PC2, color = population_description, size = 1)) +
  labs(color = NULL) +
  geom_point(size = 2) +
  geom_point(data = filter(pc_dat, population_description == 'epiSLI'), aes(color = 'epiSLI'), size = 3, shape = 'square', alpha = 0.7) +
  geom_point(data = filter(pc_dat, population_description == 'Neanderthal'), aes(color = 'Neanderthal'), size = 3, shape = 'triangle', alpha = 0.7) +
  scale_color_manual(
    values = c(
      "Iberian populations in Spain" = "coral4",
      'Toscani in Italy' = 'chocolate1',
      "Finnish in Finland" = "gold3",
      "British in England and Scotland" = "darkgreen",
      "Utah residents with Northern and Western European ancestry" = "deepskyblue",
      'epiSLI' = 'black',
      'Neanderthal' = 'red2')) +
  guides(color=guide_legend(nrow=4,byrow=TRUE)) +
  theme_classic() +
  theme(legend.position = 'bottom',
        legend.text = element_text(size = 12),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 16))

p %>% 
    ggsave(filename = "manuscript/figures/neanderthal_1000Genomes_EpiSLI_pca.png",
          device = 'png', bg = 'white', dpi = 300,
          width = 10, height = 10)