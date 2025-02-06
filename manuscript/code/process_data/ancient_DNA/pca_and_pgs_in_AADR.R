library(tidyverse)

## ===============================================================
## step 0 get overlapping SNPs (make sure ref / alt alleles match)
## ===============================================================
## check AADR alleles
aadr <- read_table('/Dedicated/jmichaelson-wdata/lcasten/tools/ref_data/allen_ancient_dna_resource/v54.1.p1_1240K_public/plink/v54.1.p1_1240K_public.bim', col_names = FALSE)
sli_1kg <- read_table('/Dedicated/jmichaelson-wdata/lcasten/sli_wgs/prs/1000_genomes_EUR_merge/merged.1000_genomes_EUR.SLI_WGS.qc.bim', col_names = FALSE)
names(sli_1kg)[-2] <- str_c('kg_', names(sli_1kg)[-2])

## add in RSID to neanderthal data where its missing
nean <- read_table('/Dedicated/jmichaelson-wdata/lcasten/tools/ref_data/archaic_hominini/hg19_neanderthal_filtered_vcf/vcf/merged_all_chromosomes.snp.bim.og', col_names = FALSE)
nean2 <- nean %>% 
    rename(kg_X1 = X1, kg_X4 = X4, nean_X2 = X2) %>% 
    left_join(sli_1kg) %>% 
    mutate(drop = case_when(nchar(X5) > 1 | nchar(X6) > 1 | X6 != kg_X6 ~ TRUE,
                            TRUE ~ FALSE)) %>%
    rename(X1 = kg_X1, X4 = kg_X4, kg_X2 = X2)

nean_clean <- nean %>% 
    mutate(idx = 1:nrow(.)) %>%
    left_join(rename(nean2, X2 = nean_X2)) %>% 
    mutate(X2 = ifelse(X2 == '.' & is.na(kg_X2) == FALSE, kg_X2, X2))

# file.copy(from = '/Dedicated/jmichaelson-wdata/lcasten/tools/ref_data/archaic_hominini/hg19_neanderthal_filtered_vcf/vcf/merged_all_chromosomes.snp.bim', to = '/Dedicated/jmichaelson-wdata/lcasten/tools/ref_data/archaic_hominini/hg19_neanderthal_filtered_vcf/vcf/merged_all_chromosomes.snp.bim.og')

tmp_pos <- nean_clean %>% 
    distinct(idx, .keep_all = TRUE)

nean_clean %>% 
    distinct(idx, .keep_all = TRUE) %>%
    select(str_c('X', 1:6)) %>% 
    write_delim('/Dedicated/jmichaelson-wdata/lcasten/tools/ref_data/archaic_hominini/hg19_neanderthal_filtered_vcf/vcf/merged_all_chromosomes.snp.bim', col_names = FALSE, delim = ' ')

aadr %>% 
    inner_join(sli_1kg) %>% 
    filter(kg_X5 == X6) %>% 
    nrow(.) ## all the AADR alleles need to be flipped (>99% of the ref/alt alleles are swapped)

# cmd <- "/Dedicated/jmichaelson-wdata/lcasten/tools/plink --bfile /Dedicated/jmichaelson-wdata/lcasten/tools/ref_data/allen_ancient_dna_resource/v54.1.p1_1240K_public/plink/v54.1.p1_1240K_public --flip /Dedicated/jmichaelson-wdata/lcasten/tools/ref_data/allen_ancient_dna_resource/v54.1.p1_1240K_public/plink/v54.1.p1_1240K_public.bim --make-bed --out /Dedicated/jmichaelson-wdata/lcasten/tools/ref_data/allen_ancient_dna_resource/v54.1.p1_1240K_public/plink/v54.1.p1_1240K_public.flipped_snp"
# system(cmd)

## ===================================================
## step 1 merge datasets and extract overlapping SNPs
## ===================================================
aadr_clean <- read_table('/Dedicated/jmichaelson-wdata/lcasten/tools/ref_data/allen_ancient_dna_resource/v54.1.p1_1240K_public/plink/v54.1.p1_1240K_public.bim', col_names = FALSE)
nean_clean <- read_table('/Dedicated/jmichaelson-wdata/lcasten/tools/ref_data/archaic_hominini/hg19_neanderthal_filtered_vcf/vcf/merged_all_chromosomes.snp.bim', col_names = FALSE)
sli_kg_clean <- sli_1kg
colnames(sli_kg_clean) = str_remove_all(colnames(sli_kg_clean), 'kg_')
# aadr_clean %>% 
#     inner_join(sli_kg_clean) %>% 
#     inner_join(nean_clean)
kp_snp <- Reduce(intersect, list(aadr_clean$X2, sli_kg_clean$X2, nean_clean$X2))
bad_snp <- read_table('/Dedicated/jmichaelson-wdata/lcasten/sli_wgs/ancient_DNA/data/merged-merge.missnp', col_names = FALSE)

data.frame(snp = kp_snp) %>% 
    as_tibble() %>%
    filter(! snp %in% bad_snp$X1) %>%
    filter(snp != '.') %>%
    write_tsv('/Dedicated/jmichaelson-wdata/lcasten/sli_wgs/ancient_DNA/data/overlapping_SNPs.snplist', col_names = FALSE)

cmd <- "/Dedicated/jmichaelson-wdata/lcasten/tools/plink --bfile /Dedicated/jmichaelson-wdata/lcasten/sli_wgs/prs/1000_genomes_EUR_merge/merged.1000_genomes_EUR.SLI_WGS.qc --extract /Dedicated/jmichaelson-wdata/lcasten/sli_wgs/ancient_DNA/data/overlapping_SNPs.snplist --make-bed --out /Dedicated/jmichaelson-wdata/lcasten/sli_wgs/ancient_DNA/data/1000Genomes_EpiSLI"
cmd1 <- "/Dedicated/jmichaelson-wdata/lcasten/tools/plink --bfile /Dedicated/jmichaelson-wdata/lcasten/tools/ref_data/allen_ancient_dna_resource/v54.1.p1_1240K_public/plink/v54.1.p1_1240K_public --extract /Dedicated/jmichaelson-wdata/lcasten/sli_wgs/ancient_DNA/data/overlapping_SNPs.snplist --make-bed --out /Dedicated/jmichaelson-wdata/lcasten/sli_wgs/ancient_DNA/data/AADR_v54"
cmd2 <- "/Dedicated/jmichaelson-wdata/lcasten/tools/plink --bfile /Dedicated/jmichaelson-wdata/lcasten/tools/ref_data/archaic_hominini/hg19_neanderthal_filtered_vcf/vcf/merged_all_chromosomes.snp --extract /Dedicated/jmichaelson-wdata/lcasten/sli_wgs/ancient_DNA/data/overlapping_SNPs.snplist --make-bed --out /Dedicated/jmichaelson-wdata/lcasten/sli_wgs/ancient_DNA/data/neanderthals"
system(cmd)
system(cmd1)
system(cmd2)

files <- c('/Dedicated/jmichaelson-wdata/lcasten/sli_wgs/ancient_DNA/data/1000Genomes_EpiSLI', 
           '/Dedicated/jmichaelson-wdata/lcasten/sli_wgs/ancient_DNA/data/AADR_v54', 
           '/Dedicated/jmichaelson-wdata/lcasten/sli_wgs/ancient_DNA/data/neanderthals')
write_lines(files, '/Dedicated/jmichaelson-wdata/lcasten/sli_wgs/ancient_DNA/data/filelist.txt')
cmd <- "/Dedicated/jmichaelson-wdata/lcasten/tools/plink --merge-list /Dedicated/jmichaelson-wdata/lcasten/sli_wgs/ancient_DNA/data/filelist.txt --extract /Dedicated/jmichaelson-wdata/lcasten/sli_wgs/ancient_DNA/data/overlapping_SNPs.snplist --allow-no-sex --make-bed --out /Dedicated/jmichaelson-wdata/lcasten/sli_wgs/ancient_DNA/data/merged"
system(cmd)

sli_1kg %>% filter(X2 == 'rs3094315')
aadr_clean %>% filter(X2 == 'rs3094315')
nean_clean %>% filter(X2 == 'rs3094315')

merged_samp <- read_table('/Dedicated/jmichaelson-wdata/lcasten/sli_wgs/ancient_DNA/data/merged.fam', col_names = FALSE) 
unique(merged_samp$X6)
table(merged_samp$X6)
merged_samp$X6 <- 1
merged_samp$X1 = ifelse(is.na(merged_samp$X1), merged_samp$X2, merged_samp$X1)

##
file.copy(from = '/Dedicated/jmichaelson-wdata/lcasten/sli_wgs/ancient_DNA/data/merged.fam', to = '/Dedicated/jmichaelson-wdata/lcasten/sli_wgs/ancient_DNA/data/merged.fam.og', overwrite = TRUE)
merged_samp %>% 
    write_delim('/Dedicated/jmichaelson-wdata/lcasten/sli_wgs/ancient_DNA/data/merged.fam', col_names = FALSE, delim = ' ')
kg_samp <- read_table('/Dedicated/jmichaelson-wdata/lcasten/sli_wgs/prs/1000_genomes_EUR_merge/merged.1000_genomes_EUR.SLI_WGS.qc.fam', col_names = FALSE) %>% 
    filter(str_detect(X2, 'sample', negate = TRUE))

merged_samp  %>% 
    filter(X2 %in% kg_samp$X2) %>% 
    select(X1, X2) %>% 
    write_delim('/Dedicated/jmichaelson-wdata/lcasten/sli_wgs/ancient_DNA/data/1000Genomes_eur.samples', col_names = FALSE, delim = ' ')

## ===================================================
## step 2: prune SNPs
## ===================================================
cmd <- "/Dedicated/jmichaelson-wdata/lcasten/tools/plink --bfile /Dedicated/jmichaelson-wdata/lcasten/sli_wgs/ancient_DNA/data/merged --keep /Dedicated/jmichaelson-wdata/lcasten/sli_wgs/ancient_DNA/data/1000Genomes_eur.samples --allow-no-sex --maf 0.05 --indep-pairwise 1000 100 0.8 --out /Dedicated/jmichaelson-wdata/lcasten/sli_wgs/ancient_DNA/data/merged.pruned_snp"
system(cmd)

## ===================================================
## step 3: PCA
## ===================================================
## make cluster file (1000 Genomes Eur vs everyone else)
merged_samp %>% 
    select(1:2) %>% 
    mutate(cluster = ifelse(X2 %in% kg_samp$X2, 1, 2)) %>% 
        write_delim('/Dedicated/jmichaelson-wdata/lcasten/sli_wgs/ancient_DNA/data/merged.cluster_file_for_pca.txt', col_names = FALSE, delim = ' ')

## PCA on 1000 genomes eur with pruned SNPs then project everyone else on that
# 
cmd <- "/Dedicated/jmichaelson-wdata/lcasten/tools/plink --bfile /Dedicated/jmichaelson-wdata/lcasten/sli_wgs/ancient_DNA/data/merged --allow-no-sex --extract /Dedicated/jmichaelson-wdata/lcasten/sli_wgs/ancient_DNA/data/merged.pruned_snp.prune.in --pca --pca-cluster-names 1 --within /Dedicated/jmichaelson-wdata/lcasten/sli_wgs/ancient_DNA/data/merged.cluster_file_for_pca.txt --out  /Dedicated/jmichaelson-wdata/lcasten/sli_wgs/ancient_DNA/data/merged.pruned_snp.1000GenomesEur_PCA_projection" # cluster = 1 means they will be used in PCA, cluster = 2 means they will be projected onto that
system(cmd)

## check PCs
kg_pc <- read_csv('/Dedicated/jmichaelson-wdata/lcasten/sli_wgs/PCA/PCA_results.merged.1000_genomes_EUR.SLI_WGS.qc.demo.csv') %>% 
    rename(IID = sample.ID) %>% 
    filter(population == '1000 Genomes') %>% 
    select(1, population_description, 4:8)
pc_new <- read_delim('/Dedicated/jmichaelson-wdata/lcasten/sli_wgs/ancient_DNA/data/merged.pruned_snp.1000GenomesEur_PCA_projection.eigenvec', col_names = FALSE, delim = ' ')
names(pc_new) <- c('FID', 'IID', str_c('new_pc', 1:20))
pc_new %>% 
    select(1:7) %>% 
    inner_join(kg_pc) %>% 
    select(matches('pc')) %>% 
    cor()

pc_new %>% 
    select(1:7) %>% 
    inner_join(kg_pc) %>% 
    ggplot(aes(x = pc1, y = pc2, color = population_description)) +
    geom_point()

pc_new %>% 
    select(1:7) %>% 
    inner_join(kg_pc) %>% 
    ggplot(aes(x = new_pc1, y = new_pc2, color = population_description)) +
    geom_point()


## ========================================================================
## ancestry outlier detection (mahalanobis distance of genetic PCs)
## ========================================================================
smeta <- read_tsv('/Dedicated/jmichaelson-wdata/lcasten/tools/ref_data/allen_ancient_dna_resource/v54.1.p1_1240K_public/v54.1.p1_1240K_public.anno') %>%
    janitor::clean_names() %>%
    rename(IID = genetic_id) %>% 
    select(IID, group_id, locality, political_entity, lat, long, molecular_sex, sample_age_years_before_1950 = date_mean_in_bp_in_years_before_1950_ce_ox_cal_mu_for_a_direct_radiocarbon_date_and_average_of_range_for_a_contextual_date, age_at_death = age_at_death_from_physical_anthropology)
smeta 
anc_samp = smeta %>% 
  filter(sample_age_years_before_1950 > 100)
anc_samp
# cmd <- "/Dedicated/jmichaelson-wdata/lcasten/tools/plink --bfile /Dedicated/jmichaelson-wdata/lcasten/sli_wgs/ancient_DNA/data/merged --allow-no-sex --extract /Dedicated/jmichaelson-wdata/lcasten/sli_wgs/ancient_DNA/data/merged.pruned_snp.prune.in  --genome --out /Dedicated/jmichaelson-wdata/lcasten/sli_wgs/ancient_DNA/data/merged_MDS"
# system(cmd)
# /Dedicated/jmichaelson-wdata/lcasten/tools/plink --bfile /Dedicated/jmichaelson-wdata/lcasten/sli_wgs/ancient_DNA/data/merged --read-genome MDS_merge2.genome --within /Dedicated/jmichaelson-wdata/lcasten/sli_wgs/ancient_DNA/data/merged.cluster_file_for_pca.txt --cluster --mds-plot 10 --out MDS_merge2
pca <- pc_new %>% 
    select(IID, str_c('new_pc', 1:5))
colnames(pca) = str_remove_all(colnames(pca), 'new_')

## compute distance in 1000 genomes Eur
pc_kg <- pca %>% 
    filter(IID %in% kg_pc$IID) %>% 
    as.data.frame()
rownames(pc_kg) <- pc_kg$IID
pc_kg$IID <- NULL
pc_kg <- pc_kg[,1:5]
pc_kg$mdist <- mahalanobis(pc_kg, center = colMeans(pc_kg), cov = cov(pc_kg))
pc_kg_p <- 1 - pchisq(pc_kg$mdist, df = 5)
hist(pc_kg_p)
sum(pc_kg_p < 0.001)
.05 / 503

## compute distance in other samples based on 1000 genomes eur
pc_non_kg <- pca %>% 
    filter(! IID %in% kg_pc$IID) %>% 
    as.data.frame()
rownames(pc_non_kg) <- pc_kg$IID
samp_id <- pc_non_kg$IID
pc_non_kg$IID <- NULL
pc_non_kg <- pc_non_kg[,1:5]
pc_non_kg$mdist <- mahalanobis(pc_non_kg, center = colMeans(pc_kg[,1:5]), cov = cov(pc_kg[,1:5]))
pc_non_kg_p <- 1 - pchisq(pc_non_kg$mdist, df = 5)
bonferroni_cutoff_anc = .05 / nrow(anc_samp)
sum(pc_non_kg_p < bonferroni_cutoff_anc)
pc_non_kg$pass_mahalanobis_qc_cutoff <- pc_non_kg_p < bonferroni_cutoff_anc

smeta %>% 
  select(political_entity, lat, long)
non_outlier_samples <- pc_non_kg %>% 
    as.data.frame() %>% 
    as_tibble() %>% 
    mutate(IID = samp_id) %>% 
    relocate(IID) %>%
    filter(mdist < median(pc_kg$mdist) + 2.5 * mad(pc_kg$mdist) & mdist > median(pc_kg$mdist) - 2.5 * mad(pc_kg$mdist)) %>%
    # filter(pass_mahalanobis_qc_cutoff == TRUE) %>%
    select(IID, mdist) %>% 
    inner_join(smeta) %>% 
    mutate(lat = as.numeric(lat),
           long = as.numeric(long)) %>% 
    # filter(long >= 25) %>% 
    # filter(lat >= 35) %>%
    filter(sample_age_years_before_1950 > 0)
hist(as.numeric(smeta$long))
hist(as.numeric(smeta$lat))

non_outlier_samples
hist(non_outlier_samples$long)
table(non_outlier_samples$political_entity)

addt_samples <- pc_non_kg %>% 
    as.data.frame() %>% 
    as_tibble() %>% 
    mutate(IID = samp_id) %>% 
    relocate(IID) %>%
    # filter(str_detect(IID, 'REF|indij|enisov|ltai|eanderth|primat')) %>%
    filter(str_detect(IID, 'REF|indija33|enisovaPi|ltaiNea$|primat|hagyrskaya-Ph')) %>%
    select(IID, mdist) %>% 
    left_join(smeta) # %>% 
    # filter(sample_age_years_before_1950 > 0 | is.na(sample_age_years_before_1950))

soi <- c(non_outlier_samples$IID, addt_samples$IID)
soi_fam <- read_delim('/Dedicated/jmichaelson-wdata/lcasten/sli_wgs/ancient_DNA/data/merged.fam', col_names = FALSE, delim = ' ') %>% 
    filter(X2 %in% soi)
soi_fam %>% 
    select(1:2) %>%
    write_delim('/Dedicated/jmichaelson-wdata/lcasten/sli_wgs/ancient_DNA/data/merged_ancient_european_samples.fam', col_names = FALSE, delim = ' ') 

cmd <- "/Dedicated/jmichaelson-wdata/lcasten/tools/plink --bfile /Dedicated/jmichaelson-wdata/lcasten/sli_wgs/ancient_DNA/data/merged --allow-no-sex --keep /Dedicated/jmichaelson-wdata/lcasten/sli_wgs/ancient_DNA/data/merged_ancient_european_samples.fam --geno 0.5 --extract /Dedicated/jmichaelson-wdata/lcasten/sli_wgs/ancient_DNA/data/merged.pruned_snp.prune.in --rel-cutoff 0.1875 --out /Dedicated/jmichaelson-wdata/lcasten/sli_wgs/ancient_DNA/data/unrelated_merged_ancient_european_samples.fam" #  --mind 0.5
system(cmd)

## ===================================================
## step 4: run ES-PGS with PRSet
## ===================================================
cmd <- 'bash /Dedicated/jmichaelson-wdata/lcasten/sli_wgs/ancient_DNA/code/ES-PGS_all_merged_samples.sh'
system(cmd)
pgs_og <- read_table('/Dedicated/jmichaelson-wdata/lcasten/sli_wgs/prs/pathway_prs/cogPerf.human_evolution.all_score')
pgs_raw <- read_delim('/Dedicated/jmichaelson-wdata/lcasten/sli_wgs/prs/pathway_prs/1000Genomes_EpiSLI_AADR_neanderthal_merged.human_evolution_complement.all_score', ' ')
pgs_og  %>% 
    select(IID, og_Base_1 = Base_1, og_haqer_1 = HAQER_1, og_har_1 = HAR_1) %>% 
    inner_join(select(pgs_raw, IID, Base_1, HAQER_1, HAR_1)) %>% 
    select(-IID) %>% 
    cor()

wd <- pgs_raw %>% 
    inner_join(pca)

##
pgs <- pgs_raw
trait = 'cogPerf'
prsice <- read_table(file = str_replace_all('/Dedicated/jmichaelson-wdata/lcasten/sli_wgs/prs/pathway_prs/1000Genomes_EpiSLI_AADR_neanderthal_merged.human_evolution_complement.all_score', pattern = 'all_score', replacement = 'prsice'),
                show_col_types = FALSE)
prsice %>% 
    filter(Set %in% c('Base', 'HAQER', 'HAR'))

names(pgs) <- str_replace_all(names(pgs), pattern = '-', replacement = '_')
tmp <- pgs %>%
    inner_join(pca)
##
kg <- tmp %>%
    filter(IID %in% kg_pc$IID)
dg <- tmp %>%
    filter(! IID %in% kg_pc$IID)
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
pgs_corrected <- bind_rows(pgs_list_tmp)

##
pgs_corrected %>% 
    write_csv('/Dedicated/jmichaelson-wdata/lcasten/sli_wgs/ancient_DNA/data/cogPerf_ES-PGS_1000GenomesEur_corrected.csv')


#####################
##
## ================================================================
## repeat with all the other pgs
## ================================================================
files = list.files('/Dedicated/jmichaelson-wdata/lcasten/sli_wgs/prs/pathway_prs', pattern = '1000Genomes_EpiSLI', full.names = TRUE)
files = files[str_detect(files, pattern = 'all_score$')]
files = files[str_detect(files, pattern = 'imputed', negate = TRUE)]
files = files[str_detect(files, pattern = 'reproductive_|circumference_GC|vocal')]
files

for(f in files) {
    ph <- basename(f)
    ph <- str_remove_all(ph, pattern = '1000Genomes_EpiSLI_AADR_neanderthal_merged.human_evolution_complement.|.all_score')
    pgs_raw <- read_delim(f, delim = ' ')

    wd <- pgs_raw %>% 
        inner_join(pca)

    ##
    pgs <- pgs_raw

    prsice <- read_table(file = str_replace_all(f, pattern = 'all_score', replacement = 'prsice'),
                        show_col_types = FALSE)
    prsice %>% 
        filter(Set %in% c('Base', 'HAQER', 'HAR'))

    names(pgs) <- str_replace_all(names(pgs), pattern = '-', replacement = '_')
    tmp <- pgs %>%
    inner_join(pca)
    ##
    kg <- tmp %>%
        filter(IID %in% kg_pc$IID)
    dg <- tmp %>%
       filter(! IID %in% kg_pc$IID)
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
    pgs_corrected <- bind_rows(pgs_list_tmp)

    ##
    pgs_corrected %>% 
        write_csv(str_c('/Dedicated/jmichaelson-wdata/lcasten/sli_wgs/ancient_DNA/data/', ph, '_ES-PGS_1000GenomesEur_corrected.csv'))
}
