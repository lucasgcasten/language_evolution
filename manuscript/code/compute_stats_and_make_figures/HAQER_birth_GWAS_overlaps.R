##
library(tidyverse)

###########################################################################################
## run enrichment analysis for previously associated birth weight / preeclampia SNPs
###########################################################################################
## hg19 chromosome info
chrom_limits = read_table('manuscript/supplemental_materials/hg19.chrom.sizes', col_names = FALSE)
names(chrom_limits) = c('#chrom', 'chrom_size')

## annotation coords
outbed <- 'manuscript/supplemental_materials/HAQER.hg19.sorted_autosomes_non_overlapping.bed'
outbed_har <- 'manuscript/supplemental_materials/HAR.hg19.sorted_autosomes_non_overlapping.bed'
outbed_rand <- 'manuscript/supplemental_materials/RAND.hg19.sorted_autosomes.bed'
outbed_uce <- 'manuscript/supplemental_materials/UCE.hg19.sorted_autosomes.bed'

## neurodevelopmental scQTL files (filtered to snpXgene paires with p < 5e-04)
files = list.files('manuscript/supplemental_materials/single_cell_QTL-Jerber-NatureGenetics-2021', pattern = '.txt', full.names = TRUE)

## GWAS hits 
# Fetal_BW_European_meta.NG2019.txt, head_circumference_at_birth_reformatted.txt, placental_weight_reformatted.txt
gw <- read_table("/sdata/gwas_summary_stats/anthropomorphic/birth_early_development/EGG_consortium/Fetal_BW_European_meta.NG2019.txt")
gw2 <- read_table("/sdata/gwas_summary_stats/anthropomorphic/birth_early_development/EGG_consortium/head_circumference_at_birth_reformatted.txt")
gw3 <- read_table("/sdata/gwas_summary_stats/anthropomorphic/birth_early_development/EGG_consortium/placental_weight_reformatted.txt")

gw4 <- read_tsv("/sdata/gwas_summary_stats/birth/preeclampsia_JAMAcardiology2023_GCST90269903/GCST90269903.tsv.gz")
gw4 <- gw4 %>% 
    filter(effect_allele_frequency >= .01 & effect_allele_frequency <= .99)

gene_coords <- read_tsv("/wdata/lcasten/tools/ref_data/hg19/gene_map.tsv")

goi <- gene_coords %>% 
    filter(symbol %in% vl$Gene) %>% 
    # filter(symbol %in% vl_genes) %>% 
    filter(`#chrom` %in% str_c("chr", 1:22))
outf = '/wdata/lcasten/tmp.bed'
outf2 = '/wdata/lcasten/tmp2.bed'

p_cutoff <- .05 / 9462395 ## number of tests
gw %>% 
    filter(eaf >= .01 & eaf <= .99) %>% 
    filter(p < 5.284074e-09) %>%
    select(chr, pos) %>%
    mutate(chromEnd = pos) %>% 
    arrange(chr, pos) %>%
    mutate(chr = str_c('chr', chr)) %>%
    rename(chromStart = pos) %>%
    inner_join(select(chrom_limits, chr = `#chrom`, chrom_size)) %>%
    mutate(chromStart = pmax(0, chromStart - 10000),
           chromEnd = pmin(chromEnd + 10000, chrom_size)) %>%
    write_tsv(outf, col_names = FALSE)

gw4 %>% 
    filter(p_value < .05 / nrow(gw4)) %>%
    select(chr = chromosome, pos = base_pair_location) %>%
    mutate(chromEnd = pos) %>% 
    arrange(chr, pos) %>%
    mutate(chr = str_c('chr', chr)) %>%
    rename(chromStart = pos) %>%
    inner_join(select(chrom_limits, chr = `#chrom`, chrom_size)) %>%
    mutate(chromStart = pmax(0, chromStart - 10000),
           chromEnd = pmin(chromEnd + 10000, chrom_size)) %>%
    write_tsv(outf, col_names = FALSE)

## 
hapmap3 <- read_rds("/wdata/lcasten/tools/LDPred2/map_hm3_plus.rds")

efn <- round(4 / (1 / 16743 + 1 / 280081), digits = 0) # 16743 European ancestry cases, 280081 controls
hapmap3 %>% 
    rename(chromosome = chr, base_pair_location = pos) %>% 
    select(rsid, chromosome, base_pair_location, a0, a1) %>%
    inner_join(gw4) %>% 
    select(SNP = rsid, chromosome, position = base_pair_location, A1 = effect_allele, A2 = other_allele, beta, se = standard_error, p = p_value, MAF = effect_allele_frequency) %>% 
    mutate(N = efn) %>%
    write_tsv("/Dedicated/jmichaelson-sdata/gwas_summary_stats/birth/preeclampsia_JAMAcardiology2023_GCST90269903/GCST90269903_hapmap3plus.tsv")

gw %>%  
    select(SNP = rsid, chromosome = chr, position = pos, A1 = ea, A2 = nea, beta, se, p, MAF = eaf, N = n) %>% 
    write_tsv("/Dedicated/jmichaelson-sdata/gwas_summary_stats/anthropomorphic/birth_early_development/EGG_consortium/Fetal_BW_European_meta_reformatted.tsv")


## merge overlapping regions for enrichment analysis
cmd <- str_c('bedtools merge -i ', outf, ' > ', outf2)
system(cmd)

## run enrichment: overlapEnrichments method elements1.lift elements2.lift noGap.lift out.txt
set_name = 'preeclampsia'
outres = str_c("/wdata/lcasten/tmp_res-", set_name, '.HAQER_enrichment.txt')
outres_har <- str_c("/wdata/lcasten/tmp_res-", set_name, '.HAR_enrichment.txt')
outres_rand <- str_c("/wdata/lcasten/tmp_res-", set_name, '.RAND_enrichment.txt')
outres_uce <- str_c("/wdata/lcasten/tmp_res-", set_name, '.UCE_enrichment.txt')
## HAQER enrichment
cmd <- str_c('~/go/bin/overlapEnrichments normalApproximate ', outf2, ' ', outbed, ' manuscript/supplemental_materials/hg19.chrom.sizes.bed ', outres)
system(cmd)
## HAR enrichment
cmd2 <- str_c('~/go/bin/overlapEnrichments normalApproximate ', outf2, ' ', outbed_har, ' manuscript/supplemental_materials/hg19.chrom.sizes.bed ', outres_har)
system(cmd2)
## RAND enrichment
cmd2 <- str_c('~/go/bin/overlapEnrichments normalApproximate ', outf2, ' ', outbed_rand, ' manuscript/supplemental_materials/hg19.chrom.sizes.bed ', outres_rand)
system(cmd2)
## UCE enrichment
cmd2 <- str_c('~/go/bin/overlapEnrichments normalApproximate ', outf2, ' ', outbed_uce, ' manuscript/supplemental_materials/hg19.chrom.sizes.bed ', outres_uce)
system(cmd2)
## print results
message('Results RAND:')
system(str_c('cat ', outres_rand, ' | cut -f 9,10,11'))
message('Results HAQERs:')
system(str_c('cat ', outres, ' | cut -f 9,10,11'))
message('Results HARs:')
system(str_c('cat ', outres_har, ' | cut -f 9,10,11'))
message('Results UCEs:')
system(str_c('cat ', outres_uce, ' | cut -f 9,10,11'))

## delete temp file
unlink(outf)
unlink(outf2)

###############
## make figure
###############
