library(tidyverse)

##
anno_all <- read_delim('/Dedicated/jmichaelson-sdata/SLI_WGS/callsets/ensemble/sli.seq.merge.variants.annotations.vep.vcfanno.ID.annotations.txt', delim = ' ', na = c('NA', '.', '', 'null', 'na', 'NULL'))
# table(anno_all$max_AF_pop)

anno_all$max_AF = ifelse(is.na(anno_all$max_AF), 0, anno_all$max_AF)
anno_all

cadd_all <- read_tsv('/Dedicated/jmichaelson-sdata/SLI_WGS/callsets/ensemble/sli.seq.merge.variants.annotations.vep.vcfanno.ID.anno.CADD.tsv', col_names = FALSE) %>% 
    rename(variant = X1, cadd = X2)

bim <- read_delim('/Dedicated/jmichaelson-wdata/lcasten/sli_wgs/regenie/data/sli.seq.merge.norm.vep.vcfanno.filt.bim', col_names = FALSE)

sli_af <- read_table('/Dedicated/jmichaelson-wdata/lcasten/sli_wgs/regenie/data/sli.seq.merge.norm.vep.vcfanno.filt.frq.frq') %>% 
    filter(NCHROBS >= 0.98 * max(NCHROBS))

## rare vars
rare_vars <- anno_all %>% 
    filter(max_AF < 0.01) %>% 
    filter(ID %in% sli_af$SNP)

rare_vars_obs <- rare_vars %>% 
    filter(max_AF > 0)

## make SNPlist files for REGENIE (to restrict variants being analyzed)
rare_vars %>% 
    select(ID) %>% 
    distinct() %>%
    write_tsv('/Dedicated/jmichaelson-wdata/lcasten/sli_wgs/regenie/data/updated_rare_variant_id.snplist', col_names = FALSE)

rare_vars_obs %>% 
    select(ID) %>% 
    distinct() %>%
    write_tsv('/Dedicated/jmichaelson-wdata/lcasten/sli_wgs/regenie/data/updated_rare_variant_id_observed_in_ref_pops.snplist', col_names = FALSE)


## make PLINK score files
sli_af %>% 
    filter(SNP %in% rare_vars$ID) %>%
    select(SNP, A1, A2) %>% 
    mutate(beta_A1 = 1) %>% 
    write_delim('/Dedicated/jmichaelson-wdata/lcasten/sli_wgs/regenie/data/all_rare_variant_weights.txt', delim = ' ')

sli_af %>% 
    filter(SNP %in% rare_vars_obs$ID) %>%
    select(SNP, A1, A2) %>% 
    mutate(beta_A1 = 1) %>% 
    write_delim('/Dedicated/jmichaelson-wdata/lcasten/sli_wgs/regenie/data/previously_observed_rare_variant_weights.txt', delim = ' ')