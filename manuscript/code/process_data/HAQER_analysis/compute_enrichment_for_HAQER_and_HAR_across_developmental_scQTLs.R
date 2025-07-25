##
library(tidyverse)

##############################################################
anno_table <- read_csv('/wdata/lcasten/sli_wgs/regenie/gene_sets/full_anno_table/gene_set_annotations.csv')
chrom_limits = read_table('/wdata/lcasten/tools/ref_data/hg19/hg19.chrom.sizes', col_names = FALSE)
names(chrom_limits) = c('#chrom', 'chrom_size')

###########################################################
## get list of constrained genes to remove from analysis
###########################################################
gene_coords <- read_tsv('/wdata/lcasten/tools/ref_data/hg19/gene_map.tsv') %>% 
    select(symbol, `#chrom`, chromStart, chromEnd) %>% 
    drop_na() %>% 
    filter(str_detect(`#chrom`, '[0-9]'))
gene_coords

pli <- read_table('/wdata/lcasten/tools/ref_data/gnomad/gnomad.v2.1.1.lof_metrics.by_gene.txt') %>% 
    select(gene, matches('pLI')) %>% 
    filter(pLI > 0.9)
loeuf <- read_tsv('/wdata/lcasten/tools/ref_data/gnomad/gnomad.v4.0.constraint_metrics.tsv') %>% 
    filter(lof.oe_ci.upper < 0.35 | lof.pLI > 0.9) %>% 
    distinct(gene)
constrained = gene_coords %>% 
    filter(symbol %in% loeuf$gene)

#################################
## run enrichment analysis
#################################

##
outbed = '/wdata/lcasten/tools/ref_data/human_accelerated_regions/HAQERs/Mangan-Cell2022/HAQER.hg19.sorted_autosomes_non_overlapping.bed'
outbed_har <- '/wdata/lcasten/tools/ref_data/human_accelerated_regions/pollard_lab/nchaes_merged_hg19_sorted.autosomes_non_overlapping.bed'

files = list.files('/wdata/lcasten/tools/ref_data/single_cell_QTL-Jerber-NatureGenetics-2021/eqtl_summary_stats_filtered', pattern = '.txt', full.names = TRUE)

##
for(f in files) {
    ph <- basename(f)
    ph = str_remove_all(ph, pattern = '.txt')
    set_name = str_c('scQTL.', ph)
    tmp <- read_tsv(f) 
    
    tmp <- tmp %>% 
        mutate(chromEnd = snp_position) %>%
        mutate(snp_position = snp_position - 1) %>%
        rename(`#chrom` = snp_chromosome, chromStart = snp_position) %>% 
        relocate(`#chrom`, chromStart, chromEnd)

    cat('\n\n\n\n')
    message('============================================')
    message('Set: ', set_name)
    ##
    outf <- str_c('/wdata/lcasten/sli_wgs/HAQER_enrichment/data/', ph, '.hg19.bed')

    tmp %>%
        mutate(`#chrom` = str_c('chr', `#chrom`)) %>%
        # filter(! symbol %in% constrained$symbol) %>% ## remove CRE influencing constrained genes
        inner_join(chrom_limits) %>%
        mutate(chr = as.numeric(str_remove_all(`#chrom`, 'chr'))) %>% 
        arrange(chr, chromStart, chromEnd) %>%
        select(`#chrom`, chromStart, chromEnd) %>% 
        write_tsv(outf, col_names = FALSE)

    ##   
    outf2 = str_c('/wdata/lcasten/sli_wgs/HAQER_enrichment/data/', set_name, '.non_overlapping.hg19.bed')
    cmd <- str_c('bedtools merge -i ', outf, ' > ', outf2)
    system(cmd)

    ## run enrichment: overlapEnrichments method         elements1.lift elements2.lift          noGap.lift                                    out.txt
    outres = str_c("/wdata/lcasten/sli_wgs/HAQER_enrichment/data/", set_name, '.txt')
    outres_har <- str_c("/wdata/lcasten/sli_wgs/HAQER_enrichment/data/", set_name, '.HAR_enrichment.txt')
    ## HAQER enrichment
    cmd <- str_c('~/go/bin/overlapEnrichments normalApproximate ', outf2, ' ', outbed, ' /wdata/lcasten/tools/ref_data/hg19/hg19.chrom_limits.bed ', outres)
    system(cmd)
    ## HAR enrichment
    cmd2 <- str_c('~/go/bin/overlapEnrichments normalApproximate ', outf2, ' ', outbed_har, ' /wdata/lcasten/tools/ref_data/hg19/hg19.chrom_limits.bed ', outres_har)
    system(cmd2)
    ## print results
    message('Results HAQERs:')
    system(str_c('cat ', outres, ' | cut -f 9,10,11'))
    message('Results HARs:')
    system(str_c('cat ', outres_har, ' | cut -f 9,10,11'))
}

#################################
## compute for bulk fetal eQTLs
#################################
qtl <- read_table('/Dedicated/jmichaelson-wdata/lcasten/tools/ref_data/fetal_human_brain_eQTL-OBrien-GenomeBiology-2018/fetal_brain_qtl.significant.txt')
var_loc <- read_tsv('/Dedicated/jmichaelson-wdata/lcasten/tools/ref_data/fetal_human_brain_eQTL-OBrien-GenomeBiology-2018/fetal_brain_qtl.significant.unique_rsid.loc.txt', col_names = FALSE)
fetal_brain_qtl <- var_loc %>% 
    select(variant_id = X5, `#chrom` = X2, chromStart = X3, chromEnd = X4) %>% 
    distinct() %>%
    inner_join(qtl) %>% 
    distinct(variant_id, .keep_all = TRUE)
fetal_brain_qtl %>% 
    write_tsv('/Dedicated/jmichaelson-wdata/lcasten/tools/ref_data/fetal_human_brain_eQTL-OBrien-GenomeBiology-2018/fetal_brain_qtl.significant.hg19.tsv')

ph <- 'fetal_brain_bulk_eQTL'
set_name = ph
tmp <- fetal_brain_qtl %>% 
    select(2:4)

cat('\n\n\n\n')
message('============================================')
message('Set: ', set_name)
##
outf <- str_c('/wdata/lcasten/sli_wgs/HAQER_enrichment/data/', ph, '.hg19.bed')

tmp %>%
    # filter(! symbol %in% constrained$symbol) %>% ## remove CRE influencing constrained genes
    inner_join(chrom_limits) %>%
    mutate(chr = as.numeric(str_remove_all(`#chrom`, 'chr'))) %>% 
    arrange(chr, chromStart, chromEnd) %>%
    select(`#chrom`, chromStart, chromEnd) %>% 
    write_tsv(outf, col_names = FALSE)

##   
outf2 = str_c('/wdata/lcasten/sli_wgs/HAQER_enrichment/data/', set_name, '.non_overlapping.hg19.bed')
cmd <- str_c('bedtools merge -i ', outf, ' > ', outf2)
system(cmd)

## run enrichment: overlapEnrichments method         elements1.lift elements2.lift          noGap.lift                                    out.txt
outres = str_c("/wdata/lcasten/sli_wgs/HAQER_enrichment/data/", set_name, '.txt')
outres_har <- str_c("/wdata/lcasten/sli_wgs/HAQER_enrichment/data/", set_name, '.HAR_enrichment.txt')
## HAQER enrichment
cmd <- str_c('~/go/bin/overlapEnrichments normalApproximate ', outf2, ' ', outbed, ' /wdata/lcasten/tools/ref_data/hg19/hg19.chrom_limits.bed ', outres)
system(cmd)
## HAR enrichment
cmd2 <- str_c('~/go/bin/overlapEnrichments normalApproximate ', outf2, ' ', outbed_har, ' /wdata/lcasten/tools/ref_data/hg19/hg19.chrom_limits.bed ', outres_har)
system(cmd2)
## print results
message('Results HAQERs:')
system(str_c('cat ', outres, ' | cut -f 9,10,11'))
message('Results HARs:')
system(str_c('cat ', outres_har, ' | cut -f 9,10,11'))



#################################
## compute for mouse brain bulk eQTLs
#################################
qtl <- read_tsv('/wdata/lcasten/tools/ref_data/mouse_brain_eQTL_Gonzales-NatureComm-2018/all_eQTL_converted_to_hg19.bed', col_names = FALSE)

ph <- 'mouse_brain_bulk_eQTL'
set_name = ph
tmp <- qtl

cat('\n\n\n\n')
message('============================================')
message('Set: ', set_name)
##
outf <- str_c('/wdata/lcasten/sli_wgs/HAQER_enrichment/data/', ph, '.hg19.bed')

tmp %>%
    dplyr::rename(`#chrom` = X1, chromStart = X2, chromEnd = X3) %>%
    # mutate(chromStart = chromStart - 100, 
    #        chromEnd = chromEnd + 100) %>%
    # filter(! symbol %in% constrained$symbol) %>% ## remove CRE influencing constrained genes
    inner_join(chrom_limits) %>%
    mutate(chr = as.numeric(str_remove_all(`#chrom`, 'chr'))) %>% 
    arrange(chr, chromStart, chromEnd) %>%
    select(`#chrom`, chromStart, chromEnd) %>% 
    write_tsv(outf, col_names = FALSE)

##   
outf2 = str_c('/wdata/lcasten/sli_wgs/HAQER_enrichment/data/', set_name, '.non_overlapping.hg19.bed')
cmd <- str_c('bedtools merge -i ', outf, ' > ', outf2)
system(cmd)

## run enrichment: overlapEnrichments method         elements1.lift elements2.lift          noGap.lift                                    out.txt
outres = str_c("/wdata/lcasten/sli_wgs/HAQER_enrichment/data/", set_name, '.txt')
outres_har <- str_c("/wdata/lcasten/sli_wgs/HAQER_enrichment/data/", set_name, '.HAR_enrichment.txt')
## HAQER enrichment
cmd <- str_c('~/go/bin/overlapEnrichments normalApproximate ', outf2, ' ', outbed, ' /wdata/lcasten/tools/ref_data/hg19/hg19.chrom_limits.bed ', outres)
system(cmd)
## HAR enrichment
cmd2 <- str_c('~/go/bin/overlapEnrichments normalApproximate ', outf2, ' ', outbed_har, ' /wdata/lcasten/tools/ref_data/hg19/hg19.chrom_limits.bed ', outres_har)
system(cmd2)
## print results
message('Results HAQERs:')
system(str_c('cat ', outres, ' | cut -f 9,10,11'))
message('Results HARs:')
system(str_c('cat ', outres_har, ' | cut -f 9,10,11'))



########################################
## psychENCODE trimester 1 eQTLs
########################################
qtl <- read_table('/Dedicated/jmichaelson-wdata/lcasten/tools/ref_data/developmental_eQTL-GandalLab/Gandal_lab_eQTL/T1-all_signifcant.txt', col_names = FALSE) %>% 
    rename(variant_id = X2)
var_loc <- read_tsv('/Dedicated/jmichaelson-wdata/lcasten/tools/ref_data/developmental_eQTL-GandalLab/Gandal_lab_eQTL/T1-all_signifcant.rsid.loc.txt', col_names = FALSE)
fetal_brain_qtl <- var_loc %>% 
    select(variant_id = X5, `#chrom` = X2, chromStart = X3, chromEnd = X4) %>% 
    distinct() %>%
    inner_join(qtl) %>% 
    distinct(variant_id, .keep_all = TRUE)
fetal_brain_qtl %>% 
    rename(qtl_gene = X1, distance_to_gene = X3, qtl_pval = X4, qtl_beta = X5) %>%
    write_tsv('/Dedicated/jmichaelson-wdata/lcasten/tools/ref_data/developmental_eQTL-GandalLab/Gandal_lab_eQTL/fetal_brain_qtl_trimester1.significant.hg19.tsv')

ph <- 'fetal_brain_eQTL_trimester1'
set_name = ph
tmp <- fetal_brain_qtl %>% 
    select(2:4)

cat('\n\n\n\n')
message('============================================')
message('Set: ', set_name)
##
outf <- str_c('/wdata/lcasten/sli_wgs/HAQER_enrichment/data/', ph, '.hg19.bed')

tmp %>%
    # filter(! symbol %in% constrained$symbol) %>% ## remove CRE influencing constrained genes
    inner_join(chrom_limits) %>%
    mutate(chr = as.numeric(str_remove_all(`#chrom`, 'chr'))) %>% 
    arrange(chr, chromStart, chromEnd) %>%
    select(`#chrom`, chromStart, chromEnd) %>% 
    write_tsv(outf, col_names = FALSE)

##   
outf2 = str_c('/wdata/lcasten/sli_wgs/HAQER_enrichment/data/', set_name, '.non_overlapping.hg19.bed')
cmd <- str_c('bedtools merge -i ', outf, ' > ', outf2)
system(cmd)

## run enrichment: overlapEnrichments method         elements1.lift elements2.lift          noGap.lift                                    out.txt
outres = str_c("/wdata/lcasten/sli_wgs/HAQER_enrichment/data/", set_name, '.txt')
outres_har <- str_c("/wdata/lcasten/sli_wgs/HAQER_enrichment/data/", set_name, '.HAR_enrichment.txt')
## HAQER enrichment
cmd <- str_c('~/go/bin/overlapEnrichments normalApproximate ', outf2, ' ', outbed, ' /wdata/lcasten/tools/ref_data/hg19/hg19.chrom_limits.bed ', outres)
system(cmd)
## HAR enrichment
cmd2 <- str_c('~/go/bin/overlapEnrichments normalApproximate ', outf2, ' ', outbed_har, ' /wdata/lcasten/tools/ref_data/hg19/hg19.chrom_limits.bed ', outres_har)
system(cmd2)
## print results
message('Results HAQERs:')
system(str_c('cat ', outres, ' | cut -f 9,10,11'))
message('Results HARs:')
system(str_c('cat ', outres_har, ' | cut -f 9,10,11'))


########################################
## psychENCODE trimester 2 eQTLs
########################################
qtl <- read_table('/Dedicated/jmichaelson-wdata/lcasten/tools/ref_data/developmental_eQTL-GandalLab/Gandal_lab_eQTL/T2-all_signifcant.txt', col_names = FALSE) %>% 
    rename(variant_id = X2)
var_loc <- read_tsv('/Dedicated/jmichaelson-wdata/lcasten/tools/ref_data/developmental_eQTL-GandalLab/Gandal_lab_eQTL/T2-all_signifcant.rsid.loc.txt', col_names = FALSE)
fetal_brain_qtl <- var_loc %>% 
    select(variant_id = X5, `#chrom` = X2, chromStart = X3, chromEnd = X4) %>% 
    distinct() %>%
    inner_join(qtl) %>% 
    distinct(variant_id, .keep_all = TRUE)
fetal_brain_qtl %>% 
    rename(qtl_gene = X1, distance_to_gene = X3, qtl_pval = X4, qtl_beta = X5) %>%
    write_tsv('/Dedicated/jmichaelson-wdata/lcasten/tools/ref_data/developmental_eQTL-GandalLab/Gandal_lab_eQTL/fetal_brain_qtl_trimester2.significant.hg19.tsv')

ph <- 'fetal_brain_eQTL_trimester2'
set_name = ph
tmp <- fetal_brain_qtl %>% 
    select(2:4)

cat('\n\n\n\n')
message('============================================')
message('Set: ', set_name)
##
outf <- str_c('/wdata/lcasten/sli_wgs/HAQER_enrichment/data/', ph, '.hg19.bed')

tmp %>%
    # filter(! symbol %in% constrained$symbol) %>% ## remove CRE influencing constrained genes
    inner_join(chrom_limits) %>%
    mutate(chr = as.numeric(str_remove_all(`#chrom`, 'chr'))) %>% 
    arrange(chr, chromStart, chromEnd) %>%
    select(`#chrom`, chromStart, chromEnd) %>% 
    write_tsv(outf, col_names = FALSE)

##   
outf2 = str_c('/wdata/lcasten/sli_wgs/HAQER_enrichment/data/', set_name, '.non_overlapping.hg19.bed')
cmd <- str_c('bedtools merge -i ', outf, ' > ', outf2)
system(cmd)

## run enrichment: overlapEnrichments method         elements1.lift elements2.lift          noGap.lift                                    out.txt
outres = str_c("/wdata/lcasten/sli_wgs/HAQER_enrichment/data/", set_name, '.txt')
outres_har <- str_c("/wdata/lcasten/sli_wgs/HAQER_enrichment/data/", set_name, '.HAR_enrichment.txt')
## HAQER enrichment
cmd <- str_c('~/go/bin/overlapEnrichments normalApproximate ', outf2, ' ', outbed, ' /wdata/lcasten/tools/ref_data/hg19/hg19.chrom_limits.bed ', outres)
system(cmd)
## HAR enrichment
cmd2 <- str_c('~/go/bin/overlapEnrichments normalApproximate ', outf2, ' ', outbed_har, ' /wdata/lcasten/tools/ref_data/hg19/hg19.chrom_limits.bed ', outres_har)
system(cmd2)
## print results
message('Results HAQERs:')
system(str_c('cat ', outres, ' | cut -f 9,10,11'))
message('Results HARs:')
system(str_c('cat ', outres_har, ' | cut -f 9,10,11'))


########################################
## psychENCODE trimester 1 eQTLs (permuted)
########################################
qtl <- read_table('/Dedicated/jmichaelson-wdata/lcasten/tools/ref_data/developmental_eQTL-GandalLab/Gandal_lab_eQTL/tri1_perm_eqtl_all_assoc.txt') %>% 
    rename(variant_id = sid) %>% 
    filter(npval < pval_nominal_threshold)
var_loc <- read_tsv('/Dedicated/jmichaelson-wdata/lcasten/tools/ref_data/developmental_eQTL-GandalLab/Gandal_lab_eQTL/tri1_perm_eqtl_all_assoc.rsid.loc.txt', col_names = FALSE)
fetal_brain_qtl <- var_loc %>% 
    select(variant_id = X5, `#chrom` = X2, chromStart = X3, chromEnd = X4) %>% 
    distinct() %>%
    inner_join(qtl) %>% 
    distinct(variant_id, .keep_all = TRUE)
fetal_brain_qtl %>% 
    write_tsv('/Dedicated/jmichaelson-wdata/lcasten/tools/ref_data/developmental_eQTL-GandalLab/Gandal_lab_eQTL/tri1_perm_eqtl_all_assoc.significant.hg19.tsv')

ph <- 'fetal_brain_eQTL_trimester1_permuted'
set_name = ph
tmp <- fetal_brain_qtl %>% 
    select(2:4)

cat('\n\n\n\n')
message('============================================')
message('Set: ', set_name)
##
outf <- str_c('/wdata/lcasten/sli_wgs/HAQER_enrichment/data/', ph, '.hg19.bed')

tmp %>%
    # filter(! symbol %in% constrained$symbol) %>% ## remove CRE influencing constrained genes
    inner_join(chrom_limits) %>%
    mutate(chr = as.numeric(str_remove_all(`#chrom`, 'chr'))) %>% 
    arrange(chr, chromStart, chromEnd) %>%
    select(`#chrom`, chromStart, chromEnd) %>% 
    write_tsv(outf, col_names = FALSE)

##   
outf2 = str_c('/wdata/lcasten/sli_wgs/HAQER_enrichment/data/', set_name, '.non_overlapping.hg19.bed')
cmd <- str_c('bedtools merge -i ', outf, ' > ', outf2)
system(cmd)

## run enrichment: overlapEnrichments method         elements1.lift elements2.lift          noGap.lift                                    out.txt
outres = str_c("/wdata/lcasten/sli_wgs/HAQER_enrichment/data/", set_name, '.txt')
outres_har <- str_c("/wdata/lcasten/sli_wgs/HAQER_enrichment/data/", set_name, '.HAR_enrichment.txt')
## HAQER enrichment
cmd <- str_c('~/go/bin/overlapEnrichments normalApproximate ', outf2, ' ', outbed, ' /wdata/lcasten/tools/ref_data/hg19/hg19.chrom_limits.bed ', outres)
system(cmd)
## HAR enrichment
cmd2 <- str_c('~/go/bin/overlapEnrichments normalApproximate ', outf2, ' ', outbed_har, ' /wdata/lcasten/tools/ref_data/hg19/hg19.chrom_limits.bed ', outres_har)
system(cmd2)
## print results
message('Results HAQERs:')
system(str_c('cat ', outres, ' | cut -f 9,10,11'))
message('Results HARs:')
system(str_c('cat ', outres_har, ' | cut -f 9,10,11'))


########################################
## psychENCODE trimester 2 eQTLs (permuted)
########################################
qtl <- read_table('/Dedicated/jmichaelson-wdata/lcasten/tools/ref_data/developmental_eQTL-GandalLab/Gandal_lab_eQTL/tri2_perm_eqtl_all_assoc.txt') %>% 
    rename(variant_id = sid) %>% 
    filter(npval < pval_nominal_threshold)
var_loc <- read_tsv('/Dedicated/jmichaelson-wdata/lcasten/tools/ref_data/developmental_eQTL-GandalLab/Gandal_lab_eQTL/tri2_perm_eqtl_all_assoc.rsid.loc.txt', col_names = FALSE)
fetal_brain_qtl <- var_loc %>% 
    select(variant_id = X5, `#chrom` = X2, chromStart = X3, chromEnd = X4) %>% 
    distinct() %>%
    inner_join(qtl) %>% 
    distinct(variant_id, .keep_all = TRUE)
fetal_brain_qtl %>% 
    write_tsv('/Dedicated/jmichaelson-wdata/lcasten/tools/ref_data/developmental_eQTL-GandalLab/Gandal_lab_eQTL/tri2_perm_eqtl_all_assoc.significant.hg19.tsv')

ph <- 'fetal_brain_eQTL_trimester2_permuted'
set_name = ph
tmp <- fetal_brain_qtl %>% 
    select(2:4)

cat('\n\n\n\n')
message('============================================')
message('Set: ', set_name)
##
outf <- str_c('/wdata/lcasten/sli_wgs/HAQER_enrichment/data/', ph, '.hg19.bed')

tmp %>%
    # filter(! symbol %in% constrained$symbol) %>% ## remove CRE influencing constrained genes
    inner_join(chrom_limits) %>%
    mutate(chr = as.numeric(str_remove_all(`#chrom`, 'chr'))) %>% 
    arrange(chr, chromStart, chromEnd) %>%
    select(`#chrom`, chromStart, chromEnd) %>% 
    write_tsv(outf, col_names = FALSE)

##   
outf2 = str_c('/wdata/lcasten/sli_wgs/HAQER_enrichment/data/', set_name, '.non_overlapping.hg19.bed')
cmd <- str_c('bedtools merge -i ', outf, ' > ', outf2)
system(cmd)

## run enrichment: overlapEnrichments method         elements1.lift elements2.lift          noGap.lift                                    out.txt
outres = str_c("/wdata/lcasten/sli_wgs/HAQER_enrichment/data/", set_name, '.txt')
outres_har <- str_c("/wdata/lcasten/sli_wgs/HAQER_enrichment/data/", set_name, '.HAR_enrichment.txt')
## HAQER enrichment
cmd <- str_c('~/go/bin/overlapEnrichments normalApproximate ', outf2, ' ', outbed, ' /wdata/lcasten/tools/ref_data/hg19/hg19.chrom_limits.bed ', outres)
system(cmd)
## HAR enrichment
cmd2 <- str_c('~/go/bin/overlapEnrichments normalApproximate ', outf2, ' ', outbed_har, ' /wdata/lcasten/tools/ref_data/hg19/hg19.chrom_limits.bed ', outres_har)
system(cmd2)
## print results
message('Results HAQERs:')
system(str_c('cat ', outres, ' | cut -f 9,10,11'))
message('Results HARs:')
system(str_c('cat ', outres_har, ' | cut -f 9,10,11'))



########################
## look at results
########################
files_dev = list.files('/wdata/lcasten/sli_wgs/HAQER_enrichment/data', pattern = 'scQTL.')
files_dev = files_dev[str_detect(files_dev, '.txt$')]

files_adult = list.files('/wdata/lcasten/sli_wgs/HAQER_enrichment/data', pattern = 'scQTL.')
files_adult = files_adult[str_detect(files_adult, '.txt$')]
files_adult = files_adult[str_detect(files_adult, '.D[0-9]', negate = TRUE)]
files_adult = files_adult[str_detect(files_adult, '.adult', negate = TRUE)]
files_adult


res_list <- list()

for(i in 1:length(files_dev)) {
    set <- basename(files_dev[i])
    set = str_remove_all(set, '.txt')
    type = ifelse(str_detect(set, 'HAR_enrichment'), 'HAR', 'HAQER')
    ct = str_split(set, pattern = '.qtl', simplify = TRUE)[,1]
    ct = str_remove_all(ct, pattern = ".qtl|scQTL.")

    res_list[[i]] <- read_table(str_c('/wdata/lcasten/sli_wgs/HAQER_enrichment/data/', files_dev[i]), show_col_types = FALSE) %>% 
        mutate(type = type) %>% 
        mutate(cell_type = ct) %>%
        relocate(cell_type, type, EnrichPValue, DepletePValue)
}

#################
dev_neur_enrich = bind_rows(res_list)
unique(dev_neur_enrich$cell_type)

##################
p_enrich_dev <- dev_neur_enrich %>% 
    # filter(type == 'HAQER') %>% 
    arrange(EnrichPValue) %>% 
    filter(str_detect(cell_type, '_treated', negate = TRUE)) %>%
    mutate(time_point = str_split(cell_type, '[.]', simplify = TRUE)[,1],
           time_point = factor(time_point, levels = c('D11', 'D30', 'D52')),
           time_point = str_replace_all(time_point, pattern = '^D', replacement = 'Day '),
           cell_type = str_split(cell_type, '[.]', simplify = TRUE)[,2]) %>%
    select(cell_type, time_point, type, matches('PValue')) %>% 
    mutate(cell_type_clean = case_when(cell_type == 'Epen1' ~ 'Ependymal',
                                       cell_type == 'DA' ~ 'Dopaminergic',
                                       cell_type == 'FPP' ~ 'Floor plate progenitors',
                                       cell_type == 'Sert' ~ 'Serotinergic',
                                       cell_type == 'Astro' ~ 'Astrocytes',
                                       cell_type == 'pseudobulk' ~ 'Any (pseudobulk)',
                                       cell_type == 'P_FPP' ~ 'Proliferating floor plate progenitors')) %>% 
    mutate(cell_type = cell_type_clean) %>%
    drop_na() %>%
    mutate(p = ifelse(EnrichPValue < DepletePValue, EnrichPValue, DepletePValue),
           stat = -log10(p),
           stat = ifelse(EnrichPValue < DepletePValue, stat, stat * -1)) %>%
    arrange(type, time_point, EnrichPValue) %>% 
    mutate(cell_type = factor(cell_type, levels = rev(unique(cell_type)))) %>%
    mutate(type = factor(type, levels = c('HAR', 'HAQER'))) %>%
    mutate(type_fill = case_when(p < 0.05 & type == 'HAQER' ~ 'HAQER_sig',
                                 p > 0.05 & type == 'HAQER' ~ 'HAQER_ns',
                                 p < 0.05 & type == 'HAR' ~ 'HAR_sig',
                                 p > 0.05 & type == 'HAR' ~ 'HAR_ns')) %>%
    mutate(sig = ifelse(p < 0.05, 1, 0.8)) %>%
    # mutate(EnrichPValue = ifelse(EnrichPValue >= 0.9 & str_detect(evo_type_clean, 'Human'), 0.9, EnrichPValue),
    #        EnrichPValue = ifelse(EnrichPValue >= 0.1 & str_detect(evo_type_clean, 'Human', negate = TRUE), 0.1, EnrichPValue)) %>%
    arrange(p) %>% # filter(str_detect(evo_type_clean, 'Conserv') & type == 'HAQER')
    ggplot(aes(x = -log10(EnrichPValue), y = cell_type, fill = type)) + # , alpha = sig
    geom_bar(stat = 'identity', position = 'dodge') +
    xlab('-log10(enrichment p-value)') +
    ylab(NULL) +
    facet_grid(rows = vars(time_point), scales = 'free_y', space = 'free') +
    geom_vline(xintercept = -log10(0.05), color = 'grey', linetype = 'dashed', size = 1.075) +
    # geom_vline(xintercept = -1 * -log10(0.05), color = 'grey', linetype = 'dashed', size = 1.075) +
    # facet_grid(cols = vars(evo_type_clean), scales = 'free') + # , space = 'free'
    scale_fill_manual(values = c('goldenrod', 'purple2')) +
    scale_alpha_continuous(guide=FALSE) +
    theme_classic() +
    labs(fill = NULL) +
    theme(legend.position = 'bottom',
          axis.text = element_text(size = 14),
          axis.title = element_text(size = 16),
          strip.text = element_text(size = 16),
          legend.text = element_text(size = 16),
          plot.title = element_text(size = 18, hjust = 0.5),
          panel.border = element_rect(colour = "black", fill=NA)) +
    ggtitle('Differentiating neuron scQTLs')

##################
## differentiate by day
##################
p_enrich_dev_day <- dev_neur_enrich %>% 
    filter(type == 'HAQER') %>% 
    arrange(EnrichPValue) %>% 
    filter(str_detect(cell_type, '_treated', negate = TRUE)) %>%
    mutate(time_point = str_split(cell_type, '[.]', simplify = TRUE)[,1],
           time_point = factor(time_point, levels = c('D11', 'D30', 'D52')),
           time_point = str_replace_all(time_point, pattern = '^D', replacement = 'Day '),
           cell_type = str_split(cell_type, '[.]', simplify = TRUE)[,2]) %>%
    mutate(time_point = factor(time_point, levels = rev(c('Day 11', 'Day 30', 'Day 52')))) %>%
    select(cell_type, time_point, type, matches('PValue')) %>% 
    mutate(cell_type_clean = case_when(cell_type == 'Epen1' ~ 'Ependymal',
                                       cell_type == 'DA' ~ 'Dopaminergic',
                                       cell_type == 'FPP' ~ 'Floor plate progenitors',
                                       cell_type == 'Sert' ~ 'Serotinergic',
                                       cell_type == 'Astro' ~ 'Astrocytes',
                                       cell_type == 'pseudobulk' ~ 'Any (pseudobulk)',
                                       cell_type == 'P_FPP' ~ 'Proliferating floor plate progenitors')) %>% 
    mutate(cell_type = cell_type_clean) %>%
    drop_na() %>%
    mutate(p = ifelse(EnrichPValue < DepletePValue, EnrichPValue, DepletePValue),
           stat = -log10(p),
           stat = ifelse(EnrichPValue < DepletePValue, stat, stat * -1)) %>%
    arrange(type, time_point, EnrichPValue) %>% 
    mutate(cell_type = factor(cell_type, levels = unique(cell_type))) %>%
    mutate(type = factor(type, levels = c('HAR', 'HAQER'))) %>%
    mutate(type_fill = case_when(p < 0.05 & type == 'HAQER' ~ 'HAQER_sig',
                                 p > 0.05 & type == 'HAQER' ~ 'HAQER_ns',
                                 p < 0.05 & type == 'HAR' ~ 'HAR_sig',
                                 p > 0.05 & type == 'HAR' ~ 'HAR_ns')) %>%
    mutate(sig = ifelse(p < 0.05, 1, 0.8)) %>%
    # mutate(EnrichPValue = ifelse(EnrichPValue >= 0.9 & str_detect(evo_type_clean, 'Human'), 0.9, EnrichPValue),
    #        EnrichPValue = ifelse(EnrichPValue >= 0.1 & str_detect(evo_type_clean, 'Human', negate = TRUE), 0.1, EnrichPValue)) %>%
    arrange(p) %>% # filter(str_detect(evo_type_clean, 'Conserv') & type == 'HAQER')
    ggplot(aes(x = -log10(EnrichPValue), y = cell_type, fill = time_point)) + # , alpha = sig
    geom_bar(stat = 'identity', position = position_dodge(preserve = "single")) +
    xlab('-log10(enrichment p-value)') +
    ylab(NULL) +
    geom_vline(xintercept = -log10(0.05), color = 'grey', linetype = 'dashed', size = 1.075) +
    # geom_vline(xintercept = -1 * -log10(0.05), color = 'grey', linetype = 'dashed', size = 1.075) +
    # facet_grid(cols = vars(evo_type_clean), scales = 'free') + # , space = 'free'
    scale_fill_manual(values = rev(c('yellow2', '#d36b04', 'darkviolet'))) +
    scale_alpha_continuous(guide=FALSE) +
    theme_classic() +
    labs(fill = NULL) +
    theme(legend.position = 'bottom',
          axis.text = element_text(size = 14),
          axis.title = element_text(size = 16),
          strip.text = element_text(size = 16),
          legend.text = element_text(size = 16),
          plot.title = element_text(size = 18, hjust = 0.5),
          panel.border = element_rect(colour = "black", fill=NA)) +
    ggtitle('Differentiating neuron scQTLs') +
    guides(fill = guide_legend(reverse=T))
p_enrich_dev_day

#################################
## repeat for adult cell-types
#################################
files_adult = list.files('/wdata/lcasten/sli_wgs/HAQER_enrichment/data', pattern = 'scQTL.')
files_adult = files_adult[str_detect(files_adult, '.txt$')]
files_adult = files_adult[str_detect(files_adult, '.D[0-9]', negate = TRUE)]
files_adult = files_adult[str_detect(files_adult, '.adult', negate = TRUE)]
files_adult


res_list <- list()

for(i in 1:length(files_adult)) {
    set <- basename(files_adult[i])
    set = str_remove_all(set, '.txt')
    type = ifelse(str_detect(set, 'HAR_enrichment'), 'HAR', 'HAQER')
    ct = str_split(set, pattern = '_sig_QTLs', simplify = TRUE)[,1]
    ct = str_remove_all(ct, pattern = "_sig_QTLs|scQTL.")

    res_list[[i]] <- read_table(str_c('/wdata/lcasten/sli_wgs/HAQER_enrichment/data/', files_adult[i]), show_col_types = FALSE) %>% 
        mutate(type = type) %>% 
        mutate(cell_type = ct) %>%
        relocate(cell_type, type, EnrichPValue, DepletePValue)
}

#################
adult_neur_enrich = bind_rows(res_list)
unique(adult_neur_enrich$cell_type)

##################
p_enrich_adult <- adult_neur_enrich %>% 
    # filter(type == 'HAQER') %>% 
    arrange(EnrichPValue) %>% 
    mutate(EnrichPValue = ifelse(EnrichPValue > 1, 1, EnrichPValue)) %>%
    mutate(cell_type = str_replace_all(cell_type, pattern = '__', replacement = '_')) %>%
    select(cell_type, type, matches('PValue')) %>% 
    # mutate(cell_type_clean = case_when(cell_type == 'Epen1' ~ 'Ependymal',
    #                                    cell_type == 'DA' ~ 'Dopaminergic',
    #                                    cell_type == 'FPP' ~ 'Floor plate progenitors',
    #                                    cell_type == 'Sert' ~ 'Serotinergic',
    #                                    cell_type == 'Astro' ~ 'Astrocytes',
    #                                    cell_type == 'pseudobulk' ~ 'Any (pseudobulk)',
    #                                    cell_type == 'P_FPP' ~ 'Proliferating floor plate progenitors')) %>% 
    # mutate(cell_type = cell_type_clean) %>%
    drop_na() %>%
    mutate(p = ifelse(EnrichPValue < DepletePValue, EnrichPValue, DepletePValue),
           stat = -log10(p),
           stat = ifelse(EnrichPValue < DepletePValue, stat, stat * -1)) %>%
    arrange(type, EnrichPValue) %>% 
    mutate(cell_type = factor(cell_type, levels = rev(unique(cell_type)))) %>%
    mutate(type = factor(type, levels = c('HAR', 'HAQER'))) %>%
    mutate(type_fill = case_when(p < 0.05 & type == 'HAQER' ~ 'HAQER_sig',
                                 p > 0.05 & type == 'HAQER' ~ 'HAQER_ns',
                                 p < 0.05 & type == 'HAR' ~ 'HAR_sig',
                                 p > 0.05 & type == 'HAR' ~ 'HAR_ns')) %>%
    mutate(sig = ifelse(p < 0.05, 1, 0.8)) %>%
    # mutate(EnrichPValue = ifelse(EnrichPValue >= 0.9 & str_detect(evo_type_clean, 'Human'), 0.9, EnrichPValue),
    #        EnrichPValue = ifelse(EnrichPValue >= 0.1 & str_detect(evo_type_clean, 'Human', negate = TRUE), 0.1, EnrichPValue)) %>%
    arrange(p) %>% # filter(str_detect(evo_type_clean, 'Conserv') & type == 'HAQER')
    ggplot(aes(x = -log10(EnrichPValue), y = cell_type, fill = type)) + # , alpha = sig
    geom_bar(stat = 'identity', position = 'dodge') +
    xlab('-log10(enrichment p-value)') +
    ylab(NULL) +
    geom_vline(xintercept = -log10(0.05), color = 'grey', linetype = 'dashed', size = 1.075) +
    # geom_vline(xintercept = -1 * -log10(0.05), color = 'grey', linetype = 'dashed', size = 1.075) +
    # facet_grid(cols = vars(evo_type_clean), scales = 'free') + # , space = 'free'
    scale_fill_manual(values = c('goldenrod', 'purple2')) +
    scale_alpha_continuous(guide=FALSE) +
    theme_classic() +
    labs(fill = NULL) +
    theme(legend.position = 'bottom',
          axis.text = element_text(size = 14),
          axis.title = element_text(size = 16),
          strip.text = element_text(size = 16),
          legend.text = element_text(size = 16),
          plot.title = element_text(size = 18, hjust = 0.5)) +
    ggtitle('Adult scQTLs')


p_enrich_dev + p_enrich_adult