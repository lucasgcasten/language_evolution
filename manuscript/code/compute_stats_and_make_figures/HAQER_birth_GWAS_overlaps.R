##
library(tidyverse)

###########################################################################################
## run enrichment analysis for previously associated birth head circumference SNPs
###########################################################################################
## hg19 chromosome info
chrom_limits = read_table('manuscript/supplemental_materials/hg19.chrom.sizes', col_names = FALSE)
names(chrom_limits) = c('#chrom', 'chrom_size')

## annotation coords
outbed <- 'manuscript/supplemental_materials/HAQER_conservative_set.hg19.bed'
outbed_har <- 'manuscript/supplemental_materials/HAR.hg19.sorted_autosomes_non_overlapping.bed'
outbed_rand <- 'manuscript/supplemental_materials/RAND.hg19.sorted_autosomes.bed'
outbed_uce <- 'manuscript/supplemental_materials/UCE.hg19.sorted_autosomes.bed'

## make tmp file paths for intermediate files
outf = '/wdata/lcasten/tmp.bed'
outf2 = '/wdata/lcasten/tmp2.bed'

## GWAS sumstats from EGG consortiums birth head circumference GWAS
gwas <- data.table::fread("/sdata/gwas_summary_stats/anthropomorphic/birth_early_development/EGG_consortium/head_circumference_at_birth_reformatted.txt")

## filter to common SNPs since this is a smaller GWAS (MAF >= 5%) and suggestive hits (p < 5e-05)
## add 100Kb flank to each hit
gwas %>% 
    filter(Freq1 >= .05 & Freq1 <= .95) %>% 
    filter(p < 5e-5) %>%
    select(chr, pos) %>%
    mutate(chromEnd = pos) %>% 
    arrange(chr, pos) %>%
    mutate(chr = str_c('chr', chr)) %>%
    rename(chromStart = pos) %>%
    inner_join(select(chrom_limits, chr = `#chrom`, chrom_size)) %>%
    mutate(chromStart = pmax(0, chromStart - 100000),
           chromEnd = pmin(chromEnd + 100000, chrom_size)) %>%
    write_tsv(outf, col_names = FALSE)

## merge overlapping hit regions for enrichment analysis
cmd <- str_c('bedtools merge -i ', outf, ' > ', outf2)
system(cmd)

## run enrichment: overlapEnrichments method elements1.lift elements2.lift noGap.lift out.txt
set_name = 'birth_head_circ'
outres = str_c("manuscript/supplemental_materials/stats/birth_head_circ_GWAS_region_enrichment/", set_name, '.HAQER_enrichment.txt')
outres_har <- str_c("manuscript/supplemental_materials/stats/birth_head_circ_GWAS_region_enrichment/", set_name, '.HAR_enrichment.txt')
outres_rand <- str_c("manuscript/supplemental_materials/stats/birth_head_circ_GWAS_region_enrichment/", set_name, '.RAND_enrichment.txt')
outres_uce <- str_c("manuscript/supplemental_materials/stats/birth_head_circ_GWAS_region_enrichment/", set_name, '.UCE_enrichment.txt')
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

## delete temp files
unlink(outf)
unlink(outf2)

#########################
## gather stats and save
#########################
files <- list.files('manuscript/supplemental_materials/stats/birth_head_circ_GWAS_region_enrichment', full.names = TRUE)

res_list = list()
for(f in files) {
    res_list[[basename(f)]] <- read_table(f) %>% 
        mutate(set = 'birth_head_circumference_EGG_Vogelezang2022',
               evo_annot = basename(Filename2)) %>% 
        relocate(evo_annot, set) %>% 
        select(-matches('Filename')) %>% 
        mutate(evo_annot = str_split(evo_annot, pattern = '[.]', simplify = TRUE)[,1]) %>% 
        select(evo_annot, set, enrichment_method = `#Method`, n_elements_evo_annot = LenElements2, n_elements_gwas_set = LenElements1, n_overlapping_elements = OverlapCount, n_expected_overlaps = ExpectedOverlap, enrichment = Enrichment, enrichment_p = EnrichPValue)
}

gwas_enr <- bind_rows(res_list) %>% 
    mutate(evo_annot = case_when(str_detect(evo_annot, 'HAQER') ~ 'HAQER',
                                 TRUE ~ evo_annot))

## save results
gwas_enr %>% 
    filter(evo_annot != 'UCE') %>%
    write_csv('manuscript/supplemental_materials/stats/HAQER_birth_head_circ_gwas_Vogelezang2022_enrichment_stats.csv')
