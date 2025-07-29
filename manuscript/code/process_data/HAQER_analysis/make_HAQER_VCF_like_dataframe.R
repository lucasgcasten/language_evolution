##################################
library(tidyverse)

##################################
## fix bim variant IDs
# bim <- read_table('/Dedicated/jmichaelson-wdata/lcasten/sli_wgs/HAQER_variant_associations/data/HAQER.sli.seq.merge.norm.vep.vcfanno.hg19.10Kb_flank.bim', col_names = FALSE)
# file.copy(from = '/Dedicated/jmichaelson-wdata/lcasten/sli_wgs/HAQER_variant_associations/data/HAQER.sli.seq.merge.norm.vep.vcfanno.hg19.10Kb_flank.bim', to = '/Dedicated/jmichaelson-wdata/lcasten/sli_wgs/HAQER_variant_associations/data/HAQER.sli.seq.merge.norm.vep.vcfanno.hg19.10Kb_flank.bim.og')
# bim$new_id <- str_c(bim$X1, '_', bim$X4, '_', bim$X6, '_', bim$X5)
# bim %>% 
#     select(X1, new_id, X3, X4, X5, X6) %>% 
#     write_delim('/Dedicated/jmichaelson-wdata/lcasten/sli_wgs/HAQER_variant_associations/data/HAQER.sli.seq.merge.norm.vep.vcfanno.hg19.10Kb_flank.bim', col_names = FALSE, delim = ' ')

## read in PLINK genotypes
bim <- read_table("/Dedicated/jmichaelson-wdata/lcasten/sli_wgs/HAQER_variant_associations/data/HAQER.sli.seq.merge.norm.vep.vcfanno.hg19.10Kb_flank.bim", col_names = FALSE) %>% 
    select(chromosome = X1, pos = X4)
geno <- BEDMatrix::BEDMatrix("/Dedicated/jmichaelson-wdata/lcasten/sli_wgs/HAQER_variant_associations/data/HAQER.sli.seq.merge.norm.vep.vcfanno.hg19.10Kb_flank.bed") ## genotypes for common + rare variants in HAQERs (only filter of MAC >= 1 and missingness < 2%)
geno <- geno %>% 
    as.matrix() %>% 
    as.data.frame()
geno <- geno[, !duplicated(colnames(geno))]
geno <- geno %>%
    rownames_to_column(var = 'IID') %>% 
    as_tibble() %>% 
    mutate(IID = str_split(IID, pattern = '_', simplify = TRUE)[,1])
message('There are ', ncol(geno) - 1, ' variants')


## read in HAQER BED file with IDs
# haqer <- read_tsv('/wdata/lcasten/tools/ref_data/human_accelerated_regions/HAQERs/Mangan-Cell2022/HAQER.hg19.bed', col_names = FALSE) ## HAQER annotations
haqer <- read_tsv('/wdata/lcasten/sli_wgs/prs/gene_sets/HAQERs_v2.hg19.bed', col_names = FALSE) ## HAQER annotations
haqer$X4 <- str_c("HAQER", 1:nrow(haqer))
haqer$X2 = haqer$X2 - 10000
haqer$X3 = haqer$X3 + 10000
haqer2 <- haqer %>% 
    mutate(X1 = as.numeric(str_remove_all(X1, pattern = 'chr'))) %>% 
    rename(chromosome = X1, start = X2, end = X3, haqer_id = X4) %>% 
    drop_na()
haqer2

## compute MAF in EpiSLI
maf <- data.frame(variant = colnames(geno[,-1]), mac = colSums(geno[,-1], na.rm = TRUE)) %>% 
    as_tibble() %>%
    # arrange(desc(mac)) %>% 
    mutate(maf = mac / 700) %>% 
    mutate(flip = ifelse(maf >= 0.99, TRUE, FALSE))
range(maf$maf)

##################################
## make annotated variant table
fathmm_annot_all = read_csv('/Dedicated/jmichaelson-wdata/lcasten/sli_wgs/HAQER_variant_associations/data/FATHMM-XF_all_HAQER_variants.10Kb.csv')
vep_annot_all <- read_delim('/Dedicated/jmichaelson-wdata/lcasten/sli_wgs/HAQER_variant_associations/data/HAQER_v2_all_variants_10Kb_flank.annotations.vep.vcfanno.ID.annotations.txt', delim = ' ', na = c('NA', '', '.', 'NULL', 'null')) %>% 
    select(-matches('clin|pheno'))
cadd_annot_all <- read_tsv('/Dedicated/jmichaelson-wdata/lcasten/sli_wgs/HAQER_variant_associations/data/HAQER_v2_all_variants_10Kb_flank.annotations.vep.vcfanno.ID.anno.CADD.tsv', col_names = FALSE) %>% 
    rename(ID = X1, cadd = X2)


## core annotations
annot <- data.frame(name = colnames(geno)[-1]) %>% 
    as_tibble() %>%
    mutate(chromosome = as.numeric(str_split(name, pattern = '_', simplify = TRUE)[,1]),
           pos = as.numeric(str_split(name, pattern = '_', simplify = TRUE)[,2]),
           ref = str_split(name, pattern = '_', simplify = TRUE)[,3],
           alt = str_split(name, pattern = '_', simplify = TRUE)[,4],
           start = pos,
           end = pos) %>% 
    relocate(chromosome, start, end) %>%
    fuzzyjoin::genome_inner_join(haqer2) %>% 
    select(chromosome = chromosome.x, pos, ref, alt, ID = name, haqer_id, haqer_start = start.y, haqer_end = end.y) %>% 
    mutate(distance_from_variant_to_haqer_start = ifelse(pos > haqer_start + 10000 & pos < haqer_end - 10000, 0, pos - (haqer_start + 10000)),
           distance_from_variant_to_haqer_end = ifelse(pos > haqer_start + 10000 & pos < haqer_end - 10000, 0, haqer_end - 10000 - pos)) %>%
    # relocate(variant_type, distance_from_variant_to_haqer_start, distance_from_variant_to_haqer_end) %>% 
    mutate(distance_to_nearest_haqer = case_when(abs(distance_from_variant_to_haqer_start) < abs(distance_from_variant_to_haqer_end) ~ abs(distance_from_variant_to_haqer_start),
                                                 abs(distance_from_variant_to_haqer_start) > abs(distance_from_variant_to_haqer_end) ~ abs(distance_from_variant_to_haqer_end),
                                                 abs(distance_from_variant_to_haqer_start) == abs(distance_from_variant_to_haqer_end) ~ abs(distance_from_variant_to_haqer_start))) %>%
    mutate(haqer_variant_type = case_when(pos > haqer_start + 10000 & pos < haqer_end - 10000 ~ 'in_HAQER',
                                    abs(distance_from_variant_to_haqer_start) < abs(distance_from_variant_to_haqer_end) ~ 'upstream_of_HAQER',
                                    abs(distance_from_variant_to_haqer_start) > abs(distance_from_variant_to_haqer_end) ~ 'downstream_of_HAQER'
                                    )) %>%
    mutate(haqer_start = haqer_start + 10000,
           haqer_end = haqer_end - 10000) %>% 
    select(-c(distance_from_variant_to_haqer_start, distance_from_variant_to_haqer_end)) %>% 
    arrange(distance_to_nearest_haqer) %>% 
    distinct(ID, .keep_all = TRUE) %>% 
    arrange(chromosome, pos) %>% 
    relocate(ID) %>% 
    mutate(ID = str_remove_all(ID, pattern = str_c('_', alt, '$'))) # ,ID = str_c('chr', ID))
hist(annot$distance_to_nearest_haqer)

########################
## ------------------------------------
## GWAS sumstats
## ------------------------------------
## ====================
########################
cog <- data.table::fread('/sdata/gwas_summary_stats/ssgac/2018/raw/GWAS_CP_all.txt.gz') %>% 
    as_tibble() %>%
    select(chromosome = CHR, pos = POS, cog_A1 = A1, cog_A2 = A2, cog_beta = Beta, cog_se = SE, cog_p = Pval)
ea <- data.table::fread('/sdata/gwas_summary_stats/ssgac/2022/EA4_additive_excl_23andMe.txt.gz') %>% 
    as_tibble() %>%
    select(chromosome = Chr, pos = BP, edu_A1 = Effect_allele, edu_A2 = Other_allele, edu_beta = Beta, edu_se = SE, edu_p = P)


########################
## ------------------------------------
## ancestral allele annotations
## ------------------------------------
## ====================
########################
## inferred HCA allele from HAQER paper
# hca <- data.table::fread('/wdata/lcasten/tools/ref_data/human_accelerated_regions/HAQERs/Mangan-Cell2022/autosome.polar.HAQERs_10Kb_flank.hg19.ancestral_allele_hg19.tsv')
hca <- data.table::fread("/wdata/lcasten/tools/ref_data/human_accelerated_regions/HAQERs/Mangan-Cell2022/autosome.polar.recompressed.hg19.ancestral_allele_hg19.reformatted.csv")
dim(hca)
hca <- hca %>% 
    as_tibble() %>%
    # mutate(HCA_allele = str_remove_all(INFO, pattern = 'AA=')) %>% 
    # mutate(chromosome = as.numeric(str_remove_all(`#CHROM`, 'chr'))) %>% 
    select(chromosome, pos, HCA_allele)

## neanderthal ref allele
# nean_all <- read_delim('/wdata/lcasten/tools/ref_data/archaic_hominini/hg19_neanderthal_filtered_vcf/vcf/merged_all_chromosomes.HCA_pos.HCA_pos.alleles.txt', col_names = FALSE, na = c('NA', '', '.', 'NULL', 'null'), delim = ' ')

# nean_all %>% 
#     filter(X5 > 0) %>% 
#     mutate(nean_allele = case_when(is.na(X4) ~ X3,
#                                    is.na(X4) == FALSE & X6 <= 0.5 ~ X3,
#                                    is.na(X4) == FALSE & X6 > 0.5 ~ X4)) %>% 
#     select(chromosome = X1, pos = X2, nean_allele) %>% 
#     distinct(chromosome, pos, .keep_all = TRUE) %>% 
#     drop_na() %>% 
#     write_csv('/wdata/lcasten/tools/ref_data/archaic_hominini/hg19_neanderthal_filtered_vcf/vcf/neanderthal_ref_allele.HCA_overlap.csv')
# hca_all %>% 
#     left_join(nean_all) %>% 
#     drop_na(nean_allele) %>% filter(nean_allele != HCA_allele)
#     write_csv('/wdata/lcasten/tools/ref_data/human_accelerated_regions/HAQERs/Mangan-Cell2022/HCA_and_neanderthal_alleles.hg19.csv')

nean <- read_delim('/wdata/lcasten/tools/ref_data/archaic_hominini/hg19_neanderthal_filtered_vcf/vcf/merged_all_chromosomes.HAQERs_v2.alleles.txt', col_names = FALSE, na = c('NA', '', '.', 'NULL', 'null'), delim = ' ')
table(nean$X5)
nean <- nean %>% 
    filter(X5 > 0)
unique(nean$X4)
nean_allele <- nean %>% 
    tidytext::unnest_tokens(input = X4, output = 'alt', token = 'regex', pattern = '[,]', to_lower = FALSE)

nean_ref <- nean %>% 
    mutate(nean_allele = case_when(is.na(X4) ~ X3,
                                   is.na(X4) == FALSE & X6 <= 0.5 ~ X3,
                                   is.na(X4) == FALSE & X6 > 0.5 ~ X4)) %>% 
    select(chromosome = X1, pos = X2, nean_allele)
nean_ref <- nean_ref %>% 
    distinct(chromosome, pos, .keep_all = TRUE) %>% 
    drop_na()
nean_ref
rm(nean)
unique(nean_ref$nean_allele)

## primate alleles
primates <- read_csv('/wdata/lcasten/tools/ref_data/human_accelerated_regions/HAQERs/Mangan-Cell2022/5wayPrimate/HAQER_alleles/Mangan_MSA_alleles_v2.csv')
primates <- primates %>% 
    mutate(chromosome = as.numeric(str_remove_all(chr, 'chr'))) %>% 
    select(-chr) %>% 
    rename(pos = pos_hg19)

## merge all ancestral / species ref alleles
anc_alleles <- bim %>% 
    left_join(mutate(hca, pos = as.numeric(pos))) %>% 
    left_join(nean_ref) %>% 
    left_join(primates) %>%
    # left_join(chimp_ref) %>% 
    # left_join(gor_ref) %>% 
    # left_join(mac_ref) %>% 
    # left_join(ms_ref) %>% 
    # left_join(zf_ref)
    distinct(chromosome, pos, .keep_all = TRUE) %>% 
    rename(neanderthal_allele = nean_allele)
anc_alleles
colSums(is.na(anc_alleles))

anc_alleles %>% 
    filter(human_allele == neanderthal_allele)

#################################
## merge all annotations to one dataframe
annot_wd <- vep_annot_all %>% 
    inner_join(mutate(annot, ID = str_c('chr', ID))) %>% 
    left_join(cadd_annot_all) %>% 
    left_join(fathmm_annot_all) %>% 
    mutate(ID = str_remove_all(ID, 'chr')) %>%
    mutate(variant = str_c(ID, '_', alt)) %>% 
    inner_join(select(maf, variant, AF_EpiSLI = maf)) %>%
    relocate(ID, chromosome, pos, ref, alt, rsid, gene, worst_consequence, impact, cadd, fathmm = fathmm_score, AF_EpiSLI, max_AF, max_AF_pop) %>%
    mutate(max_AF = ifelse(is.na(max_AF), 0, max_AF))

annos <- annot_wd %>% 
    left_join(cog) %>% 
    left_join(ea) %>% 
    left_join(anc_alleles) %>% 
    relocate(hg38_ref_allele = human_allele, HCA_allele, neanderthal_allele, chimp_allele, bonobo_allele, gorilla_allele, orangutan_allele = orangatuan_allele, .after = alt) %>% 
    relocate(pos_hg38, .after = pos)
annos <- annos %>% 
    mutate(cog_beta = case_when(cog_A1 == ref & cog_A2 == alt ~ -1 * cog_beta,
                                TRUE ~ cog_beta),
           edu_beta = case_when(edu_A1 == ref & edu_A2 == alt ~ -1 * edu_beta,
                                TRUE ~ edu_beta)
           ) %>% 
    select(-c(cog_A1, cog_A2, edu_A1, edu_A2))
names(annos)

annos %>% 
    filter(ref != HCA_allele | is.na(HCA_allele))

## look at AF and discrepancies
annot_wd %>% 
    relocate(AF_EpiSLI, max_AF)

annot_wd %>% 
    filter(AF_EpiSLI > 0.1 & max_AF == 0)

obs_var <- annot_wd %>% 
    filter(max_AF > 0)
table(annot_wd$worst_consequence)

annot_wd %>% 
    filter(max_AF > 0 & max_AF < 0.01 & AF_EpiSLI < 0.01 & haqer_variant_type == 'in_HAQER')

## flip genotype DF to be variants x samples
geno_vcf <- data.table::transpose(geno[,-1])
colnames(geno_vcf) <- geno$IID 
geno_vcf <- geno_vcf %>% 
    as.data.frame() %>% 
    as_tibble() %>% 
    mutate(variant = colnames(geno)[-1]) %>% 
    relocate(variant)
geno_vcf

## merge annotations with genotypes (VCF style dataframe)
dat <- annos %>% 
    inner_join(geno_vcf)
names(dat)

## save as RDS
colSums(is.na(dat))
dat %>% 
    relocate(HCA_allele, .after = neanderthal_allele) %>%
    # as.data.frame() %>% 
    select(-variant) %>%
    # select(-hg38_ref_allele) %>%
    write_rds('/Dedicated/jmichaelson-wdata/lcasten/sli_wgs/HAQER_variant_associations/data/HAQER_v2_all_variants_10Kb_flank.vcf.rds')
dat %>% 
    filter(chimp_allele == '-' & is.na(HCA_allele) == TRUE)

## norm to hg19
tdat <- readRDS("/wdata/lcasten/sli_wgs/HAQER_variant_associations/data/HAQER_all_variants.hg19_normed.vcf.rds") 
tdat %>% filter(distance_to_nearest_haqer == 0) %>% nrow() 
names(dat)
dat2 <- dat %>% 
    filter(distance_to_nearest_haqer == 0)
hg19 <- data.table::fread('/wdata/lcasten/sli_wgs/hg19_ref_alt.tsv')
hg19 <- hg19 %>% 
    rename(chromosome = `#CHROM`, pos = POS, rsid = ID) %>% 
    mutate(chromosome = as.numeric(str_remove_all(chromosome, pattern = 'chr')))
dat3 <- dat2 %>% 
    rename(major = ref, minor = alt) %>% 
    inner_join(select(hg19, -rsid)) %>%
    filter(REF == major | REF == minor) %>% 
    filter(ALT == major | ALT == minor) %>% 
    pivot_longer(cols = matches('sample')) %>% 
    mutate(to_flip = ifelse(REF == minor & ALT == major, TRUE, FALSE),
            value = case_when(REF == major & ALT == minor ~ value,
                             REF == minor & ALT == major ~ -1 * value + 2,
                             TRUE ~ value)) %>% 
    pivot_wider(id_cols = -matches("^name$|^value$")) %>% 
    distinct(ID, .keep_all = TRUE) %>% 
    relocate(REF, ALT, .before = major) %>% 
    rename(ref = REF, alt = ALT)
dim(dat3)
table(dat3$to_flip)

colnames(tdat)
colnames(dat3)

setdiff(colnames(tdat), colnames(dat3))
setdiff( colnames(dat3), colnames(tdat))

id1 <- dat3 %>% drop_na(HCA_allele, neanderthal_allele, gorilla_allele) %>% select(ID)

dat3 %>% 
    as.data.frame() %>% 
    select(-c(hg38_ref_allele, variant, to_flip)) %>%
    write_rds("/wdata/lcasten/sli_wgs/HAQER_variant_associations/data/HAQER_v2_all_variants.hg19_normed.vcf.rds")

# tmp <- read_rds("/wdata/lcasten/sli_wgs/HAQER_variant_associations/data/HAQER_v2_all_variants.hg19_normed.vcf.rds")
# id2 <- tmp %>% drop_na(HCA_allele, neanderthal_allele, gorilla_allele) %>% as_tibble() %>% select(ID)
# tmp %>% 
#     filter(! ID %in% dat3$ID) %>% 
#     as_tibble()
# setdiff(id1$ID, id2$ID)

dat3 %>% 
    drop_na(neanderthal_allele, orangutan_allele)

dat %>% 
    filter(ref != hg38_ref_allele & hg38_ref_allele != alt)
dat %>% 
    filter(alt == hg38_ref_allele)

################################33

# ## chimp ref allele
# haqer_merge <- haqer %>% 
#     select(chromosome = X1, start = X2, end = X3)
# chimp <- read_csv('/sdata/SLI_WGS/callsets/ensemble/annotations/pt-alleles.csv') %>% 
#     select(1:2, ref)
# names(chimp) = c('chromosome', 'start', 'chimp_allele')
# chimp <- chimp %>% 
#     mutate(end = start) %>% 
#     relocate(end, .after = start)
# chimp_ref <- chimp %>% 
#     fuzzyjoin::genome_inner_join(haqer_merge) %>% 
#     select(chromosome = chromosome.x, pos = start.x, chimp_allele) %>% 
#     mutate(chromosome = as.numeric(str_remove_all(chromosome, 'chr'))) %>% 
#     distinct(chromosome, pos, .keep_all = TRUE)
# rm(chimp)

# ## gorilla ref allele
# gor <- read_csv('/wdata/lcasten/sli_wgs/evolution/gg5-gorilla-alleles.csv') %>% 
#     select(1:2, ref)
# names(gor) = c('chromosome', 'start', 'gorilla_allele')
# gor <- gor %>% 
#     mutate(end = start) %>% 
#     relocate(end, .after = start)
# gor_ref <- gor %>% 
#     fuzzyjoin::genome_inner_join(haqer_merge) %>% 
#     select(chromosome = chromosome.x, pos = start.x, gorilla_allele) %>% 
#     mutate(chromosome = as.numeric(str_remove_all(chromosome, 'chr')))
# rm(gor)

# ## macaque ref allele
# mac <- read_csv('/wdata/lcasten/sli_wgs/evolution/mm8-macaque-alleles.csv') %>% 
#     select(1:2, ref)
# names(mac) = c('chromosome', 'start', 'macaque_allele')
# mac <- mac %>% 
#     mutate(end = start) %>% 
#     relocate(end, .after = start)
# mac_ref <- mac %>% 
#     fuzzyjoin::genome_inner_join(haqer_merge) %>% 
#     select(chromosome = chromosome.x, pos = start.x, macaque_allele) %>% 
#     mutate(chromosome = as.numeric(str_remove_all(chromosome, 'chr')))
# rm(mac)

# ## mouse ref allele
# ms <- read_csv('/wdata/lcasten/sli_wgs/evolution/mm10-alleles.csv') %>% 
#     select(1:2, ref)
# names(ms) = c('chromosome', 'start', 'mouse_allele')
# ms <- ms %>% 
#     mutate(end = start) %>% 
#     relocate(end, .after = start)
# ms_ref <- ms %>% 
#     fuzzyjoin::genome_inner_join(haqer_merge) %>% 
#     select(chromosome = chromosome.x, pos = start.x, mouse_allele) %>% 
#     mutate(chromosome = as.numeric(str_remove_all(chromosome, 'chr')))
# rm(ms)

# ## zebrafish ref allele
# zf <- read_csv('/wdata/lcasten/sli_wgs/evolution/dr10-zebrafish-alleles.csv') %>% 
#     select(1:2, ref)
# names(zf) = c('chromosome', 'start', 'zebrafish_allele')
# zf <- zf %>% 
#     mutate(end = start) %>% 
#     relocate(end, .after = start)
# zf_ref <- zf %>% 
#     fuzzyjoin::genome_inner_join(haqer_merge) %>% 
#     select(chromosome = chromosome.x, pos = start.x, zebrafish_allele) %>% 
#     mutate(chromosome = as.numeric(str_remove_all(chromosome, 'chr')))
# rm(zf)
# zf_ref %>% drop_na()
# unique(zf_ref$zebrafish_allele)
