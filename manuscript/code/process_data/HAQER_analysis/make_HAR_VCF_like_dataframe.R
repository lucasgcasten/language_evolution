##################################
library(tidyverse)

##################################
## fix bim variant IDs
# bim <- read_table('/Dedicated/jmichaelson-wdata/lcasten/sli_wgs/HAQER_variant_associations/data/HAR.sli.seq.merge.norm.vep.vcfanno.hg19.10Kb_flank.bim', col_names = FALSE)
# file.copy(from = '/Dedicated/jmichaelson-wdata/lcasten/sli_wgs/HAQER_variant_associations/data/HAR.sli.seq.merge.norm.vep.vcfanno.hg19.10Kb_flank.bim', to = '/Dedicated/jmichaelson-wdata/lcasten/sli_wgs/HAQER_variant_associations/data/HAR.sli.seq.merge.norm.vep.vcfanno.hg19.10Kb_flank.bim.og')
# bim$new_id <- str_c(bim$X1, '_', bim$X4, '_', bim$X6, '_', bim$X5)
# bim %>% 
#     select(X1, new_id, X3, X4, X5, X6) %>% 
#     write_delim('/Dedicated/jmichaelson-wdata/lcasten/sli_wgs/HAQER_variant_associations/data/HAR.sli.seq.merge.norm.vep.vcfanno.hg19.10Kb_flank.bim', col_names = FALSE, delim = ' ')

## read in PLINK genotypes
geno <- BEDMatrix::BEDMatrix("/Dedicated/jmichaelson-wdata/lcasten/sli_wgs/HAQER_variant_associations/data/HAR.sli.seq.merge.norm.vep.vcfanno.hg19.10Kb_flank.bed") ## genotypes for common + rare variants in HARs (only filter of MAC >= 1 and missingness < 2%)
geno <- geno %>% 
    as.matrix() %>% 
    as.data.frame()
geno <- geno[, !duplicated(colnames(geno))]
geno <- geno %>%
    rownames_to_column(var = 'IID') %>% 
    as_tibble() %>% 
    mutate(IID = str_split(IID, pattern = '_', simplify = TRUE)[,1])
message('There are ', ncol(geno) - 1, ' variants')


## read in HAR BED file with IDs
har <- read_tsv('/wdata/lcasten/tools/ref_data/human_accelerated_regions/pollard_lab/nchaes_merged_hg19_sorted.bed', col_names = FALSE) ## HAR annotations
har$X2 = har$X2 - 10000
har$X3 = har$X3 + 10000
har2 <- har %>% 
    mutate(X1 = as.numeric(str_remove_all(X1, pattern = 'chr'))) %>% 
    rename(chromosome = X1, start = X2, end = X3, har_id = X4) %>% 
    drop_na()
har2

## compute MAF in EpiSLI
maf <- data.frame(variant = colnames(geno[,-1]), mac = colSums(geno[,-1], na.rm = TRUE)) %>% 
    as_tibble() %>%
    # arrange(desc(mac)) %>% 
    mutate(maf = mac / 700) %>% 
    mutate(flip = ifelse(maf >= 0.99, TRUE, FALSE))
range(maf$maf)

##################################
## make annotated variant table
read_csv('/Dedicated/jmichaelson-wdata/lcasten/sli_wgs/HAQER_variant_associations/data/FATHMM-XF_all_HAR_variants.10Kb.csv', n_max = 10)

fathmm_annot_coding <- read_tsv('/wdata/lcasten/sli_wgs/HAQER_variant_associations/data/HAR_all_variants_10Kb_flank.annotations.vep.vcfanno.ID.anno.FATHMM-XF-coding.tsv', col_names = FALSE) %>% 
    rename(ID = X1, fathmm = X2) %>% 
    mutate(type = 'coding')
fathmm_annot_noncoding <- read_tsv('/wdata/lcasten/sli_wgs/HAQER_variant_associations/data/HAR_all_variants_10Kb_flank.annotations.vep.vcfanno.ID.anno.FATHMM-XF-noncoding.tsv', col_names = FALSE) %>% 
    rename(ID = X1, fathmm = X2) %>% 
    mutate(type = 'noncoding')
fathmm_annot_all <- bind_rows(fathmm_annot_coding, fathmm_annot_noncoding) %>% 
    arrange(desc(fathmm)) %>% 
    group_by(ID) %>% 
    slice_head(n = 1) 
table(fathmm_annot_all$type)
fathmm_annot_all <- fathmm_annot_all %>% 
    select(-type) %>% 
    ungroup()

vep_annot_all <- read_delim('/Dedicated/jmichaelson-wdata/lcasten/sli_wgs/HAQER_variant_associations/data/HAR_all_variants_10Kb_flank.annotations.vep.vcfanno.ID.annotations.txt', delim = ' ', na = c('NA', '', '.', 'NULL', 'null')) %>% 
    select(-matches('clin|pheno'))
cadd_annot_all <- read_tsv('/Dedicated/jmichaelson-wdata/lcasten/sli_wgs/HAQER_variant_associations/data/HAR_all_variants_10Kb_flank.annotations.vep.vcfanno.ID.anno.CADD.tsv', col_names = FALSE) %>% 
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
    fuzzyjoin::genome_inner_join(har2) %>% 
    select(chromosome = chromosome.x, pos, ref, alt, ID = name, har_id, har_start = start.y, har_end = end.y) %>% 
    mutate(distance_from_variant_to_har_start = ifelse(pos > har_start + 10000 & pos < har_end - 10000, 0, pos - (har_start + 10000)),
           distance_from_variant_to_har_end = ifelse(pos > har_start + 10000 & pos < har_end - 10000, 0, har_end - 10000 - pos)) %>%
    # relocate(variant_type, distance_from_variant_to_har_start, distance_from_variant_to_har_end) %>% 
    mutate(distance_to_nearest_har = case_when(abs(distance_from_variant_to_har_start) < abs(distance_from_variant_to_har_end) ~ abs(distance_from_variant_to_har_start),
                                                 abs(distance_from_variant_to_har_start) > abs(distance_from_variant_to_har_end) ~ abs(distance_from_variant_to_har_end),
                                                 abs(distance_from_variant_to_har_start) == abs(distance_from_variant_to_har_end) ~ abs(distance_from_variant_to_har_start))) %>%
    mutate(har_variant_type = case_when(pos > har_start + 10000 & pos < har_end - 10000 ~ 'in_HAR',
                                    abs(distance_from_variant_to_har_start) < abs(distance_from_variant_to_har_end) ~ 'upstream_of_HAR',
                                    abs(distance_from_variant_to_har_start) > abs(distance_from_variant_to_har_end) ~ 'downstream_of_HAR'
                                    )) %>%
    mutate(har_start = har_start + 10000,
           har_end = har_end - 10000) %>% 
    select(-c(distance_from_variant_to_har_start, distance_from_variant_to_har_end)) %>% 
    arrange(distance_to_nearest_har) %>% 
    distinct(ID, .keep_all = TRUE) %>% 
    arrange(chromosome, pos) %>% 
    relocate(ID) %>% 
    mutate(ID = str_remove_all(ID, pattern = str_c('_', alt, '$')),
           ID = str_c('chr', ID))
hist(annot$distance_to_nearest_har)

## merge all annotations to one dataframe
annot_wd <- vep_annot_all %>% 
    inner_join(annot) %>% 
    left_join(cadd_annot_all) %>% 
    left_join(fathmm_annot_all) %>% 
    mutate(ID = str_remove_all(ID, 'chr')) %>%
    mutate(variant = str_c(ID, '_', alt)) %>% 
    inner_join(select(maf, variant, AF_EpiSLI = maf)) %>%
    relocate(ID, chromosome, pos, ref, alt, rsid, gene, worst_consequence, impact, cadd, fathmm, AF_EpiSLI, max_AF, max_AF_pop) %>%
    mutate(max_AF = ifelse(is.na(max_AF), 0, max_AF))

## look at AF and discrepancies
annot_wd %>% 
    relocate(AF_EpiSLI, max_AF)

annot_wd %>% 
    filter(AF_EpiSLI > 0.1 & max_AF == 0)

obs_var <- annot_wd %>% 
    filter(max_AF > 0)
table(annot_wd$worst_consequence)

annot_wd %>% 
    filter(max_AF > 0 & max_AF < 0.01 & AF_EpiSLI < 0.01 & har_variant_type == 'in_HAR')

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
dat <- annot_wd %>% 
    inner_join(geno_vcf)
names(dat)

## save as RDS
dat %>% 
    as.data.frame() %>% 
    select(-variant) %>%
    write_rds('/Dedicated/jmichaelson-wdata/lcasten/sli_wgs/HAQER_variant_associations/data/HAR_all_variants_10Kb_flank.vcf.rds')
