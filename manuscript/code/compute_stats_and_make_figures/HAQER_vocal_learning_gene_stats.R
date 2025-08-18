##
library(tidyverse)

###########################################################################################
## run enrichment analysis for previously associated vocal learning genes across species
###########################################################################################
## hg19 chromosome info
chrom_limits = read_table('manuscript/supplemental_materials/hg19.chrom.sizes', col_names = FALSE)
names(chrom_limits) = c('#chrom', 'chrom_size')

## annotation coords
outbed <- 'manuscript/supplemental_materials/HAQER.hg19.sorted_autosomes_non_overlapping.bed'
outbed_har <- 'manuscript/supplemental_materials/HAR.hg19.sorted_autosomes_non_overlapping.bed'
outbed_rand <- 'manuscript/supplemental_materials/RAND.hg19.sorted_autosomes.bed'
outbed_uce <- 'manuscript/supplemental_materials/UCE.hg19.sorted_autosomes.bed'

## vl genes
vl <- read_csv("/wdata/lcasten/tools/ref_data/vocal_learning_cross_species-WirthlinScience2024/science.abn3263_data_s1_to_s10 (1)/science.abn3263_data_s2.csv") %>% 
    filter(P < (.05 / nrow(.)))

gene_coords <- read_tsv("/wdata/lcasten/tools/ref_data/hg19/gene_map.tsv")

goi <- gene_coords %>% 
    filter(symbol %in% vl$Gene) %>% 
    # filter(symbol %in% vl_genes) %>% 
    filter(`#chrom` %in% str_c("chr", 1:22))
outf = '/wdata/lcasten/tmp.bed'
outf2 = '/wdata/lcasten/tmp2.bed'
goi %>% 
    select(1:3) %>% 
    mutate(chromStart = pmax(0, chromStart - 1000000),
           chromEnd = chromEnd + 1000000) %>%
    write_tsv(outf)

## merge overlapping regions for enrichment analysis
cmd <- str_c('bedtools merge -i ', outf, ' > ', outf2)
system(cmd)

## run enrichment: overlapEnrichments method elements1.lift elements2.lift noGap.lift out.txt
set_name = 'Wirthlin_vocal_learning_genes'
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





##################
## try for tacit MM10 coords
read_tsv("/wdata/lcasten/tools/ref_data/vocal_learning_cross_species-WirthlinScience2024/table_s9_MCX_OCRs_hg19_liftover_coords.bed", col_names = FALSE) %>% 
    inner_join(rename(chrom_limits, X1 = `#chrom`)) %>%
    mutate(X2 = pmax(0, X2 - 1000000),
           X3 = pmin(chrom_size, X3 + 1000000)) %>% 
    write_tsv(outf, col_names = FALSE) 
outf = '/wdata/lcasten/tools/ref_data/vocal_learning_cross_species-WirthlinScience2024/table_s9_MCX_OCRs_hg19_liftover_coords_flanks.bed'

outf2 = '/wdata/lcasten/tmp2.bed'

## merge overlapping regions for enrichment analysis
cmd <- str_c('bedtools merge -i ', outf, ' > ', outf2)
system(cmd)

## run enrichment: overlapEnrichments method elements1.lift elements2.lift noGap.lift out.txt
set_name = 'Wirthlin_mm10_TACIT_MCX_liftover'
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

## --------------------------------------------------------------------------------------------------
## analyze both genes + OCRs associated w/ vocal learning in Wirthlin 2024 paper (with 1Mb flanks)
## --------------------------------------------------------------------------------------------------
outf = '/wdata/lcasten/tools/ref_data/vocal_learning_cross_species-WirthlinScience2024/RERconverge_genes_and_MCX_OCRs_hg19.bed'
outf2 = '/wdata/lcasten/tmp2.bed'

## merge vocal learning associated OCRs w/ genes associated with vocal learning
read_tsv("/wdata/lcasten/tools/ref_data/vocal_learning_cross_species-WirthlinScience2024/table_s9_MCX_OCRs_hg19_liftover_coords.bed", col_names = FALSE) %>% 
    bind_rows(select(goi, X1 = 1, X2 = 2, X3 = 3)) %>%
    inner_join(rename(chrom_limits, X1 = `#chrom`)) %>%
    mutate(X2 = pmax(0, X2 - 1000000),
           X3 = pmin(chrom_size, X3 + 1000000)) %>% 
    mutate(chr = as.numeric(str_remove_all(X1, 'chr'))) %>% 
    arrange(chr, X2, X3) %>%
    select(-c(chr, chrom_size)) %>%
    write_tsv(outf, col_names = FALSE) 

## merge overlapping regions for enrichment analysis
cmd <- str_c('bedtools merge -i ', outf, ' > ', outf2)
system(cmd)

## run enrichment: overlapEnrichments method elements1.lift elements2.lift noGap.lift out.txt
set_name = 'Wirthlin_mm10_TACIT_MCX_and_RERconverge_genes'
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