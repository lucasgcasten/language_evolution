###########################################################################################
## run enrichment analysis for previously associated vocal learning loci across species
###########################################################################################
##
library(tidyverse)

## hg19 chromosome info
chrom_limits = read_table('manuscript/supplemental_materials/hg19.chrom.sizes', col_names = FALSE)
names(chrom_limits) = c('#chrom', 'chrom_size')

## annotation coordinates in BED format
outbed <- 'manuscript/supplemental_materials/HAQER_conservative_set.hg19.bed'
outbed_har <- 'manuscript/supplemental_materials/HAR.hg19.sorted_autosomes_non_overlapping.bed'
outbed_rand <- 'manuscript/supplemental_materials/RAND.hg19.sorted_autosomes.bed'
outbed_uce <- 'manuscript/supplemental_materials/UCE.hg19.sorted_autosomes.bed'

## add 1Mb flanks to vocal learning associated OCRs identified in Wirthlin et al 2024 Science paper
outf = 'manuscript/supplemental_materials/cross_species_vocal_learning_OCRs_MCX_PV_Wirthlin2024_table_s9_hg19_liftover_1Mb_flanks.bed'
outf2 = '/wdata/lcasten/tmp2.bed'
read_tsv("manuscript/supplemental_materials/cross_species_vocal_learning_OCRs_MCX_PV_Wirthlin2024_table_s9_hg19_liftover.bed", col_names = FALSE) %>% 
    inner_join(rename(chrom_limits, X1 = `#chrom`)) %>%
    mutate(X2 = pmax(0, X2 - 1000000),
           X3 = pmin(chrom_size, X3 + 1000000)) %>% 
    write_tsv(outf, col_names = FALSE) 


## merge overlapping regions for enrichment analysis
cmd <- str_c('bedtools merge -i ', outf, ' > ', outf2)
system(cmd)

## run enrichment: overlapEnrichments method elements1.lift elements2.lift noGap.lift out.txt
set_name = 'vocal_learning_OCRs_MCX_PV_Wirthlin2024'
outres = str_c("manuscript/supplemental_materials/stats/HAQER_vocal_learning_region_enrichment/", set_name, '.HAQER_enrichment.txt')
outres_har <- str_c("manuscript/supplemental_materials/stats/HAQER_vocal_learning_region_enrichment/", set_name, '.HAR_enrichment.txt')
outres_rand <- str_c("manuscript/supplemental_materials/stats/HAQER_vocal_learning_region_enrichment/", set_name, '.RAND_enrichment.txt')
outres_uce <- str_c("manuscript/supplemental_materials/stats/HAQER_vocal_learning_region_enrichment/", set_name, '.UCE_enrichment.txt')
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

#########################
## gather stats and save
#########################
files <- list.files('manuscript/supplemental_materials/stats/HAQER_vocal_learning_region_enrichment', full.names = TRUE)

res_list = list()
for(f in files) {
    res_list[[basename(f)]] <- read_table(f) %>% 
        mutate(set = 'vocal_learning_OCRs_MCX_PV_Wirthlin2024',
               evo_annot = basename(Filename2)) %>% 
        relocate(evo_annot, set) %>% 
        select(-matches('Filename')) %>% 
        mutate(evo_annot = str_split(evo_annot, pattern = '[.]', simplify = TRUE)[,1]) %>% 
        select(evo_annot, set, enrichment_method = `#Method`, n_elements_evo_annot = LenElements2, n_elements_scQTL_set = LenElements1, n_overlapping_elements = OverlapCount, n_expected_overlaps = ExpectedOverlap, enrichment = Enrichment, enrichment_p = EnrichPValue)
}

haqer_vl <- bind_rows(res_list) %>% 
    mutate(evo_annot = case_when(str_detect(evo_annot, 'HAQER') ~ 'HAQER',
                                 TRUE ~ evo_annot))

## save results
haqer_vl %>% 
    filter(evo_annot != 'UCE') %>% 
    write_csv('manuscript/supplemental_materials/stats/HAQER_vocal_learning_Wirthlin2024_enrichment_stats.csv')
