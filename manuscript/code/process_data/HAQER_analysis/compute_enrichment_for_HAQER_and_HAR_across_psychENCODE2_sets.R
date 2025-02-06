##
library(tidyverse)

##############################################################
## make BED files of coords for each gene set of interest
haqer <- read_tsv("/wdata/lcasten/tools/ref_data/human_accelerated_regions/HAQERs/Mangan-Cell2022/HAQER.hg38.bed", col_names = FALSE)
haqer %>% 
    mutate(chr = str_remove_all(X1, 'chr'),
           chr = as.numeric(chr)) %>% 
    drop_na(chr) %>% 
    arrange(chr, X2, X3) %>% 
    select(1:3) %>% 
    write_tsv('/wdata/lcasten/tools/ref_data/human_accelerated_regions/HAQERs/Mangan-Cell2022/HAQER.hg38.sorted_autosomes.bed', col_names = FALSE)

gene_coords <- read_tsv('/wdata/lcasten/tools/ref_data/hg38/gene_map.tsv')
anno_table <- read_csv('/wdata/lcasten/sli_wgs/regenie/gene_sets/full_anno_table/gene_set_annotations.csv')
# wget https://hgdownload.soe.ucsc.edu/goldenpath/hg38/bigZips/hg38.chrom.sizes
chrom_limits <- read_table('/wdata/lcasten/tools/ref_data/hg38/hg38.chrom.sizes', col_names = FALSE) %>% 
    filter(X1 %in% str_c('chr', 1:22)) %>% 
    rename(`#chrom` = X1, chrom_size = X2)
chrom_limits %>% 
    mutate(start = 1) %>% 
    relocate(`#chrom`, start, chrom_size) %>% 
    write_tsv('/wdata/lcasten/tools/ref_data/hg38/hg38.chrom.sizes.bed', col_names = FALSE)

## remove self-overlapping regions in HAQERs
haqer_path = "/wdata/lcasten/tools/ref_data/human_accelerated_regions/HAQERs/Mangan-Cell2022/HAQER.hg38.sorted_autosomes.bed"
outbed = '/wdata/lcasten/tools/ref_data/human_accelerated_regions/HAQERs/Mangan-Cell2022/HAQER.hg38.sorted_autosomes_non_overlapping.bed'
cmd <- str_c('bedtools merge -i ', haqer_path, ' > ', outbed)
system(cmd)

##############################
##################################################
#############################
## repeat for HARs
#############################
## remove self-overlapping regions in HARs
har <- "/wdata/lcasten/tools/ref_data/human_accelerated_regions/HAQERs/Mangan-Cell2022/other_regions/har.hg38.bed"
read_tsv(har, col_names = FALSE) %>% 
    mutate(chr = as.numeric(str_remove_all(X1, 'chr'))) %>% 
    drop_na(chr) %>% 
    arrange(chr, X2, X3) %>%
    select(1:3) %>% 
    write_tsv('/wdata/lcasten/tools/ref_data/human_accelerated_regions/HAQERs/Mangan-Cell2022/other_regions/har.hg38_sorted.autosomes.bed', col_names = FALSE)
har_path = "/wdata/lcasten/tools/ref_data/human_accelerated_regions/HAQERs/Mangan-Cell2022/other_regions/har.hg38_sorted.autosomes.bed"
outbed = '/wdata/lcasten/tools/ref_data/human_accelerated_regions/HAQERs/Mangan-Cell2022/other_regions/har.hg38_sorted.autosomes_non_overlapping.bed'
cmd <- str_c('bedtools merge -i ', har_path, ' > ', outbed)
system(cmd)


###########################################################
## get list of constrained genes to remove from analysis
###########################################################
list.files('/wdata/lcasten/tools/ref_data/hg38')
gene_coords <- read_tsv('/wdata/lcasten/tools/ref_data/hg38/hg38_hgnc.tsv') %>% 
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
files = list.files('/wdata/lcasten/tools/ref_data/psychENCODE2/TFBS/bulkATAC_ChromBPNet', pattern = '.bed$', full.names = TRUE)
files

##
outbed = '/wdata/lcasten/tools/ref_data/human_accelerated_regions/HAQERs/Mangan-Cell2022/HAQER.hg38.sorted_autosomes_non_overlapping.bed'
outbed_har <- '/wdata/lcasten/tools/ref_data/human_accelerated_regions/HAQERs/Mangan-Cell2022/other_regions/har.hg38_sorted.autosomes_non_overlapping.bed'

##
for(f in files) {
    ph <- basename(f)
    ph = str_remove_all(ph, pattern = '.profile_scores.bw.|.bb.bed.annotated.bed')
    ph = str_c(ph, '.FOXP2')
    set_name = ph
    tmp <- read_tsv(f, col_names = FALSE) 
    
    if(str_detect(ph, 'VLPFC', negate = TRUE)) {
        tmp <- tmp %>% 
            # mutate(tf = ifelse(str_detect(ph, 'VLPFC_glia', negate = TRUE), str_split(X7, pattern = ' ', simplify = TRUE)[,2], X7)) %>%
            mutate(tf = str_split(X7, pattern = ' ', simplify = TRUE)[,2]) %>%
            filter(tf == 'FOXP2') %>% 
            filter(str_detect(X1, '[0-9]')) %>% 
            rename(`#chrom` = X1, chromStart = X2, chromEnd = X3)
    }
    if(str_detect(ph, 'VLPFC')) {
        tmp <- tmp %>% 
            # mutate(tf = ifelse(str_detect(ph, 'VLPFC_glia', negate = TRUE), str_split(X7, pattern = ' ', simplify = TRUE)[,2], X7)) %>%
            mutate(tf = X7) %>%
            filter(tf == 'FOXP2') %>% 
            filter(str_detect(X1, '[0-9]')) %>% 
            rename(`#chrom` = X1, chromStart = X2, chromEnd = X3)
    }

    cat('\n\n\n\n')
    message('============================================')
    message('Set: ', set_name)
    ##
    outf <- str_c('/wdata/lcasten/sli_wgs/HAQER_enrichment/data/', ph, '.hg38.bed')

    tmp %>%
        # filter(! symbol %in% constrained$symbol) %>% ## remove CRE influencing constrained genes
        inner_join(chrom_limits) %>%
        # mutate(chromStart = ifelse(chromStart - 10000 >= 0, chromStart - 10000, 1),
            #    chromEnd = ifelse(chromEnd + 10000 <= chrom_size, chromEnd + 10000, chrom_size)) %>% 
        mutate(chr = as.numeric(str_remove_all(`#chrom`, 'chr'))) %>% 
        arrange(chr, chromStart, chromEnd) %>%
        select(`#chrom`, chromStart, chromEnd) %>% 
        write_tsv(outf, col_names = FALSE)

    ##   
    outf2 = str_c('/wdata/lcasten/sli_wgs/HAQER_enrichment/data/', set_name, '.non_overlapping.hg38.bed')
    cmd <- str_c('bedtools merge -i ', outf, ' > ', outf2)
    system(cmd)

    ## run enrichment: overlapEnrichments method         elements1.lift elements2.lift          noGap.lift                                    out.txt
    outres = str_c("/wdata/lcasten/sli_wgs/HAQER_enrichment/data/", set_name, '.txt')
    outres_har <- str_c("/wdata/lcasten/sli_wgs/HAQER_enrichment/data/", set_name, '.HAR_enrichment.txt')
    ## HAQER enrichment
    cmd <- str_c('~/go/bin/overlapEnrichments exact ', outf2, ' ', outbed, ' /wdata/lcasten/tools/ref_data/hg38/hg38.chrom.sizes.bed ', outres)
    system(cmd)
    ## HAR enrichment
    cmd2 <- str_c('~/go/bin/overlapEnrichments exact ', outf2, ' ', outbed_har, ' /wdata/lcasten/tools/ref_data/hg38/hg38.chrom.sizes.bed ', outres_har)
    system(cmd2)
    ## print results
    message('Results HAQERs:')
    system(str_c('cat ', outres, ' | cut -f 9,10,11'))
    message('Results HARs:')
    system(str_c('cat ', outres_har, ' | cut -f 9,10,11'))
}


#################################
## run enrichment analysis for Bryois 2022 adult scQTLs
#################################
files = list.files('/wdata/lcasten/tools/ref_data/sc_eQTL_adult_brain-Bryois-NatureNeuroscience-2022', pattern = '.tsv$', full.names = TRUE)
files

##
outbed = '/wdata/lcasten/tools/ref_data/human_accelerated_regions/HAQERs/Mangan-Cell2022/HAQER.hg38.sorted_autosomes_non_overlapping.bed'
outbed_har <- '/wdata/lcasten/tools/ref_data/human_accelerated_regions/HAQERs/Mangan-Cell2022/other_regions/har.hg38_sorted.autosomes_non_overlapping.bed'

##
for(f in files) {
    ph <- basename(f)
    ph = str_split(ph, pattern = '.Bryois2022NN_', simplify = TRUE)[,1]
    set_name = str_c('scQTL.adult.', ph)
    tmp <- read_tsv(f, show_col_types = FALSE) %>% 
            filter(str_detect(chrom, '[0-9]')) %>% 
            select(2:3) %>%
            mutate(chromStart = position - 1) %>%
            rename(`#chrom` = chrom, chromEnd = position) %>% 
            select(`#chrom`, chromStart, chromEnd) %>% 
            mutate(`#chrom` = str_c('chr', `#chrom`))
   
    cat('\n\n\n\n')
    message('============================================')
    message('Set: ', set_name)
    message('There are ', nrow(tmp), ' eQTLs for this cell type')
    ##
    outf <- str_c('/wdata/lcasten/sli_wgs/HAQER_enrichment/data/', ph, '.hg38.bed')

    tmp %>%
        # filter(! symbol %in% constrained$symbol) %>% ## remove CRE influencing constrained genes
        inner_join(chrom_limits) %>%
        # mutate(chromStart = ifelse(chromStart - 10000 >= 0, chromStart - 10000, 1),
            #    chromEnd = ifelse(chromEnd + 10000 <= chrom_size, chromEnd + 10000, chrom_size)) %>% 
        mutate(chr = as.numeric(str_remove_all(`#chrom`, 'chr'))) %>% 
        arrange(chr, chromStart, chromEnd) %>%
        select(`#chrom`, chromStart, chromEnd) %>% 
        write_tsv(outf, col_names = FALSE)

    ##   
    outf2 = str_c('/wdata/lcasten/sli_wgs/HAQER_enrichment/data/', set_name, '.non_overlapping.hg38.bed')
    cmd <- str_c('bedtools merge -i ', outf, ' > ', outf2)
    system(cmd)

    ## run enrichment: overlapEnrichments method         elements1.lift elements2.lift          noGap.lift                                    out.txt
    outres = str_c("/wdata/lcasten/sli_wgs/HAQER_enrichment/data/", set_name, '.txt')
    outres_har <- str_c("/wdata/lcasten/sli_wgs/HAQER_enrichment/data/", set_name, '.HAR_enrichment.txt')
    ## HAQER enrichment
    cmd <- str_c('~/go/bin/overlapEnrichments normalApproximate ', outf2, ' ', outbed, ' /wdata/lcasten/tools/ref_data/hg38/hg38.chrom.sizes.bed ', outres)
    system(cmd)
    ## HAR enrichment
    cmd2 <- str_c('~/go/bin/overlapEnrichments normalApproximate ', outf2, ' ', outbed_har, ' /wdata/lcasten/tools/ref_data/hg38/hg38.chrom.sizes.bed ', outres_har)
    system(cmd2)
    ## print results
    message('Results HAQERs:')
    system(str_c('cat ', outres, ' | cut -f 9,10,11'))
    message('Results HARs:')
    system(str_c('cat ', outres_har, ' | cut -f 9,10,11'))
}




#################################
## run enrichment analysis for MetaBrain, Klein 2023 bulk eQTLs
#################################
files = list.files('/Dedicated/jmichaelson-wdata/lcasten/tools/ref_data/MetaBrain_eQTL-Klein-NatureGenetics-2023', pattern = '.significant.txt$', full.names = TRUE)
files

##
outbed = '/wdata/lcasten/tools/ref_data/human_accelerated_regions/HAQERs/Mangan-Cell2022/HAQER.hg38.sorted_autosomes_non_overlapping.bed'
outbed_har <- '/wdata/lcasten/tools/ref_data/human_accelerated_regions/HAQERs/Mangan-Cell2022/other_regions/har.hg38_sorted.autosomes_non_overlapping.bed'

##
for(f in files) {
    ph <- basename(f)
    ph = str_split(ph, pattern = '[-]', simplify = TRUE)[,4]
    set_name = str_c('MetaBrain.', ph)
    tmp <- read_tsv(f, show_col_types = FALSE) %>% 
            filter(str_detect(SNPChr, '[0-9]')) %>% 
            select(SNPChr, SNPPos) %>%
            mutate(chromStart = SNPPos - 1) %>%
            rename(`#chrom` = SNPChr, chromEnd = SNPPos) %>% 
            select(`#chrom`, chromStart, chromEnd) %>% 
            mutate(`#chrom` = str_c('chr', `#chrom`))
   
    cat('\n\n\n\n')
    message('============================================')
    message('Set: ', set_name)
    message('There are ', nrow(tmp), ' eQTLs for this brain region')
    ##
    outf <- str_c('/wdata/lcasten/sli_wgs/HAQER_enrichment/data/', ph, '.hg38.bed')

    tmp %>%
        # filter(! symbol %in% constrained$symbol) %>% ## remove CRE influencing constrained genes
        inner_join(chrom_limits) %>%
        # mutate(chromStart = ifelse(chromStart - 10000 >= 0, chromStart - 10000, 1),
            #    chromEnd = ifelse(chromEnd + 10000 <= chrom_size, chromEnd + 10000, chrom_size)) %>% 
        mutate(chr = as.numeric(str_remove_all(`#chrom`, 'chr'))) %>% 
        arrange(chr, chromStart, chromEnd) %>%
        select(`#chrom`, chromStart, chromEnd) %>% 
        write_tsv(outf, col_names = FALSE)

    ##   
    outf2 = str_c('/wdata/lcasten/sli_wgs/HAQER_enrichment/data/', set_name, '.non_overlapping.hg38.bed')
    cmd <- str_c('bedtools merge -i ', outf, ' > ', outf2)
    system(cmd)

    ## run enrichment: overlapEnrichments method         elements1.lift elements2.lift          noGap.lift                                    out.txt
    outres = str_c("/wdata/lcasten/sli_wgs/HAQER_enrichment/data/", set_name, '.txt')
    outres_har <- str_c("/wdata/lcasten/sli_wgs/HAQER_enrichment/data/", set_name, '.HAR_enrichment.txt')
    ## HAQER enrichment
    cmd <- str_c('~/go/bin/overlapEnrichments normalApproximate ', outf2, ' ', outbed, ' /wdata/lcasten/tools/ref_data/hg38/hg38.chrom.sizes.bed ', outres)
    system(cmd)
    ## HAR enrichment
    cmd2 <- str_c('~/go/bin/overlapEnrichments normalApproximate ', outf2, ' ', outbed_har, ' /wdata/lcasten/tools/ref_data/hg38/hg38.chrom.sizes.bed ', outres_har)
    system(cmd2)
    ## print results
    message('Results HAQERs:')
    system(str_c('cat ', outres, ' | cut -f 9,10,11'))
    message('Results HARs:')
    system(str_c('cat ', outres_har, ' | cut -f 9,10,11'))
}




#################################
## run enrichment analysis for scQTLs
#################################
files = list.files('/wdata/lcasten/tools/ref_data/psychENCODE2/cell_type_eQTL', pattern = '.dat$', full.names = TRUE)
files

##
outbed = '/wdata/lcasten/tools/ref_data/human_accelerated_regions/HAQERs/Mangan-Cell2022/HAQER.hg38.sorted_autosomes_non_overlapping.bed'
outbed_har <- '/wdata/lcasten/tools/ref_data/human_accelerated_regions/HAQERs/Mangan-Cell2022/other_regions/har.hg38_sorted.autosomes_non_overlapping.bed'

##
for(f in files) {
    ph <- basename(f)
    ph = str_remove_all(ph, pattern = '.dat')
    set_name = str_c('scQTL.', ph)
    tmp <- read_table(f, col_names = FALSE)  %>% 
            filter(str_detect(X9, '[0-9]')) %>% 
            select(9:11) %>%
            rename(`#chrom` = X9, chromStart = X10, chromEnd = X11) %>% 
            mutate(chromStart = ifelse(chromStart == chromEnd, chromStart - 1, chromStart))
   
    cat('\n\n\n\n')
    message('============================================')
    message('Set: ', set_name)
    ##
    outf <- str_c('/wdata/lcasten/sli_wgs/HAQER_enrichment/data/', ph, '.hg38.bed')

    tmp %>%
        # filter(! symbol %in% constrained$symbol) %>% ## remove CRE influencing constrained genes
        inner_join(chrom_limits) %>%
        # mutate(chromStart = ifelse(chromStart - 10000 >= 0, chromStart - 10000, 1),
            #    chromEnd = ifelse(chromEnd + 10000 <= chrom_size, chromEnd + 10000, chrom_size)) %>% 
        mutate(chr = as.numeric(str_remove_all(`#chrom`, 'chr'))) %>% 
        arrange(chr, chromStart, chromEnd) %>%
        select(`#chrom`, chromStart, chromEnd) %>% 
        write_tsv(outf, col_names = FALSE)

    ##   
    outf2 = str_c('/wdata/lcasten/sli_wgs/HAQER_enrichment/data/', set_name, '.non_overlapping.hg38.bed')
    cmd <- str_c('bedtools merge -i ', outf, ' > ', outf2)
    system(cmd)

    ## run enrichment: overlapEnrichments method         elements1.lift elements2.lift          noGap.lift                                    out.txt
    outres = str_c("/wdata/lcasten/sli_wgs/HAQER_enrichment/data/", set_name, '.txt')
    outres_har <- str_c("/wdata/lcasten/sli_wgs/HAQER_enrichment/data/", set_name, '.HAR_enrichment.txt')
    ## HAQER enrichment
    cmd <- str_c('~/go/bin/overlapEnrichments normalApproximate ', outf2, ' ', outbed, ' /wdata/lcasten/tools/ref_data/hg38/hg38.chrom.sizes.bed ', outres)
    system(cmd)
    ## HAR enrichment
    cmd2 <- str_c('~/go/bin/overlapEnrichments normalApproximate ', outf2, ' ', outbed_har, ' /wdata/lcasten/tools/ref_data/hg38/hg38.chrom.sizes.bed ', outres_har)
    system(cmd2)
    ## print results
    message('Results HAQERs:')
    system(str_c('cat ', outres, ' | cut -f 9,10,11'))
    message('Results HARs:')
    system(str_c('cat ', outres_har, ' | cut -f 9,10,11'))
}


########################################
## enrichment for gtex bulk eQTLs
########################################
files = list.files('/wdata/lcasten/tools/ref_data/gtex/GTEx_Analysis_v8_eQTL', pattern = 'signif_variant_gene_pairs.', full.names = TRUE)
files = files[str_detect(files, 'Brain_')]
files

##
outbed = '/wdata/lcasten/tools/ref_data/human_accelerated_regions/HAQERs/Mangan-Cell2022/HAQER.hg38.sorted_autosomes_non_overlapping.bed'
outbed_har <- '/wdata/lcasten/tools/ref_data/human_accelerated_regions/HAQERs/Mangan-Cell2022/other_regions/har.hg38_sorted.autosomes_non_overlapping.bed'

##
for(f in files) {
    ph <- basename(f)
    ph = str_split(ph, pattern = '[.]', simplify = TRUE)[,1]
    set_name = str_c('GTEx_v8_bulk_eQTL.', ph)
    tmp <- read_table(f, show_col_types = FALSE) %>% 
            mutate(`#chrom` = str_split(variant_id, pattern = '_', simplify = TRUE)[,1],
                   chromEnd = as.numeric(str_split(variant_id, pattern = '_', simplify = TRUE)[,2]),
                   chromStart = chromEnd - 1) %>%
            filter(str_detect(`#chrom`, '[0-9]')) %>% 
            select(`#chrom`, chromStart, chromEnd)
   
    cat('\n\n\n\n')
    message('============================================')
    message('Set: ', set_name)
    message('There are ', nrow(tmp), ' eQTLs for this tissue')

    ##
    outf <- str_c('/wdata/lcasten/sli_wgs/HAQER_enrichment/data/', ph, '.hg38.bed')

    tmp %>%
        # filter(! symbol %in% constrained$symbol) %>% ## remove CRE influencing constrained genes
        inner_join(chrom_limits) %>%
        # mutate(chromStart = ifelse(chromStart - 10000 >= 0, chromStart - 10000, 1),
            #    chromEnd = ifelse(chromEnd + 10000 <= chrom_size, chromEnd + 10000, chrom_size)) %>% 
        mutate(chr = as.numeric(str_remove_all(`#chrom`, 'chr'))) %>% 
        arrange(chr, chromStart, chromEnd) %>%
        select(`#chrom`, chromStart, chromEnd) %>% 
        write_tsv(outf, col_names = FALSE)

    ##   
    outf2 = str_c('/wdata/lcasten/sli_wgs/HAQER_enrichment/data/', set_name, '.non_overlapping.hg38.bed')
    cmd <- str_c('bedtools merge -i ', outf, ' > ', outf2)
    system(cmd)

    ## run enrichment: overlapEnrichments method         elements1.lift elements2.lift          noGap.lift                                    out.txt
    outres = str_c("/wdata/lcasten/sli_wgs/HAQER_enrichment/data/", set_name, '.txt')
    outres_har <- str_c("/wdata/lcasten/sli_wgs/HAQER_enrichment/data/", set_name, '.HAR_enrichment.txt')
    ## HAQER enrichment
    cmd <- str_c('~/go/bin/overlapEnrichments normalApproximate ', outf2, ' ', outbed, ' /wdata/lcasten/tools/ref_data/hg38/hg38.chrom.sizes.bed ', outres)
    system(cmd)
    ## HAR enrichment
    cmd2 <- str_c('~/go/bin/overlapEnrichments normalApproximate ', outf2, ' ', outbed_har, ' /wdata/lcasten/tools/ref_data/hg38/hg38.chrom.sizes.bed ', outres_har)
    system(cmd2)
    ## print results
    message('Results HAQERs:')
    system(str_c('cat ', outres, ' | cut -f 9,10,11'))
    message('Results HARs:')
    system(str_c('cat ', outres_har, ' | cut -f 9,10,11'))
}

#################################
## run enrichment analysis
#################################
files = list.files('/wdata/lcasten/tools/ref_data/psychENCODE2/brain-cCREs', pattern = '.bed$', full.names = TRUE)
files

##
outbed = '/wdata/lcasten/tools/ref_data/human_accelerated_regions/HAQERs/Mangan-Cell2022/HAQER.hg38.sorted_autosomes_non_overlapping.bed'
outbed_har <- '/wdata/lcasten/tools/ref_data/human_accelerated_regions/HAQERs/Mangan-Cell2022/other_regions/har.hg38_sorted.autosomes_non_overlapping.bed'

##
for(f in files) {
    ph <- basename(f)
    ph = str_remove_all(ph, pattern = '.bed')
    set_name = ph
    tmp <- read_tsv(f, col_names = FALSE)  %>% 
            filter(str_detect(X1, '[0-9]')) %>% 
            rename(`#chrom` = X1, chromStart = X2, chromEnd = X3)
   
    cat('\n\n\n\n')
    message('============================================')
    message('Set: ', set_name)
    ##
    outf <- str_c('/wdata/lcasten/sli_wgs/HAQER_enrichment/data/', ph, '.hg38.bed')

    tmp %>%
        # filter(! symbol %in% constrained$symbol) %>% ## remove CRE influencing constrained genes
        inner_join(chrom_limits) %>%
        # mutate(chromStart = ifelse(chromStart - 10000 >= 0, chromStart - 10000, 1),
            #    chromEnd = ifelse(chromEnd + 10000 <= chrom_size, chromEnd + 10000, chrom_size)) %>% 
        mutate(chr = as.numeric(str_remove_all(`#chrom`, 'chr'))) %>% 
        arrange(chr, chromStart, chromEnd) %>%
        select(`#chrom`, chromStart, chromEnd) %>% 
        write_tsv(outf, col_names = FALSE)

    ##   
    outf2 = str_c('/wdata/lcasten/sli_wgs/HAQER_enrichment/data/', set_name, '.non_overlapping.hg38.bed')
    cmd <- str_c('bedtools merge -i ', outf, ' > ', outf2)
    system(cmd)

    ## run enrichment: overlapEnrichments method         elements1.lift elements2.lift          noGap.lift                                    out.txt
    outres = str_c("/wdata/lcasten/sli_wgs/HAQER_enrichment/data/", set_name, '.txt')
    outres_har <- str_c("/wdata/lcasten/sli_wgs/HAQER_enrichment/data/", set_name, '.HAR_enrichment.txt')
    ## HAQER enrichment
    cmd <- str_c('~/go/bin/overlapEnrichments normalApproximate ', outf2, ' ', outbed, ' /wdata/lcasten/tools/ref_data/hg38/hg38.chrom.sizes.bed ', outres)
    system(cmd)
    ## HAR enrichment
    cmd2 <- str_c('~/go/bin/overlapEnrichments normalApproximate ', outf2, ' ', outbed_har, ' /wdata/lcasten/tools/ref_data/hg38/hg38.chrom.sizes.bed ', outres_har)
    system(cmd2)
    ## print results
    message('Results HAQERs:')
    system(str_c('cat ', outres, ' | cut -f 9,10,11'))
    message('Results HARs:')
    system(str_c('cat ', outres_har, ' | cut -f 9,10,11'))
}


################

#################################
## run enrichment analysis for DEGs
#################################
files = list.files('/wdata/lcasten/tools/ref_data/psychENCODE2/DEGs', pattern = '.csv$', full.names = TRUE)
files

##
outbed = '/wdata/lcasten/tools/ref_data/human_accelerated_regions/HAQERs/Mangan-Cell2022/HAQER.hg38.sorted_autosomes_non_overlapping.bed'
outbed_har <- '/wdata/lcasten/tools/ref_data/human_accelerated_regions/HAQERs/Mangan-Cell2022/other_regions/har.hg38_sorted.autosomes_non_overlapping.bed'

##
for(f in files) {
    ph <- basename(f)
    ph = str_remove_all(ph, pattern = 'combined.csv')
    set_name = ph
    tmp <- read_csv(f) %>% 
            filter(padj < 0.05) %>% 
            distinct(gene)
   
    cat('\n\n\n\n')
    message('============================================')
    message('Set: ', set_name)
    ##
    outf <- str_c('/wdata/lcasten/sli_wgs/HAQER_enrichment/data/', ph, '.hg38.bed')

    gene_coords %>% 
        filter(symbol %in% tmp$gene) %>%
        # filter(! symbol %in% constrained$symbol) %>% ## remove  constrained genes
        inner_join(chrom_limits) %>%
        mutate(chromStart = ifelse(chromStart - 35000 >= 0, chromStart - 35000, 1),
               chromEnd = ifelse(chromEnd + 10000 <= chrom_size, chromEnd + 10000, chrom_size)) %>% 
        mutate(chr = as.numeric(str_remove_all(`#chrom`, 'chr'))) %>% 
        arrange(chr, chromStart, chromEnd) %>%
        select(`#chrom`, chromStart, chromEnd) %>% 
        write_tsv(outf, col_names = FALSE)

    ##   
    outf2 = str_c('/wdata/lcasten/sli_wgs/HAQER_enrichment/data/', set_name, '.non_overlapping.hg38.bed')
    cmd <- str_c('bedtools merge -i ', outf, ' > ', outf2)
    system(cmd)

    ## run enrichment: overlapEnrichments method         elements1.lift elements2.lift          noGap.lift                                    out.txt
    outres = str_c("/wdata/lcasten/sli_wgs/HAQER_enrichment/data/", set_name, '.txt')
    outres_har <- str_c("/wdata/lcasten/sli_wgs/HAQER_enrichment/data/", set_name, '.HAR_enrichment.txt')
    ## HAQER enrichment
    cmd <- str_c('~/go/bin/overlapEnrichments normalApproximate ', outf2, ' ', outbed, ' /wdata/lcasten/tools/ref_data/hg38/hg38.chrom.sizes.bed ', outres)
    system(cmd)
    ## HAR enrichment
    cmd2 <- str_c('~/go/bin/overlapEnrichments normalApproximate ', outf2, ' ', outbed_har, ' /wdata/lcasten/tools/ref_data/hg38/hg38.chrom.sizes.bed ', outres_har)
    system(cmd2)
    ## print results
    message('Results HAQERs:')
    system(str_c('cat ', outres, ' | cut -f 9,10,11'))
    message('Results HARs:')
    system(str_c('cat ', outres_har, ' | cut -f 9,10,11'))
}





#################################
## run enrichment analysis for human specific CRE's
#################################
cre <- read_tsv('/wdata/lcasten/tools/ref_data/scATAC_brain_atlas_Li-Science-2024/NIHMS1959369-supplement-Supplementary_Tables__S1-S29_/Supplementary Files/Table S6_List of cCREs in bed format', col_names = FALSE)
cre
cre_annot <- read_tsv('/wdata/lcasten/tools/ref_data/scATAC_brain_atlas_Li-Science-2024/NIHMS1959369-supplement-Supplementary_Tables__S1-S29_/Supplementary Files/Table S20 - List of different categories of cCREs.txt')
cre_annot 
unique(cre_annot$subclass)
cre_annot %>% 
    filter(subclass == 'MSN') %>% 
    group_by(category) %>% 
    count()

hs_cre <- cre %>% 
    rename(cCRE = X4) %>% 
    rename(`#chrom` = X1, chromStart = X2, chromEnd = X3) %>% 
    inner_join(cre_annot) %>% 
    filter(category == 'humanSpec')

##
outbed = '/wdata/lcasten/tools/ref_data/human_accelerated_regions/HAQERs/Mangan-Cell2022/HAQER.hg38.sorted_autosomes_non_overlapping.bed'
outbed_har <- '/wdata/lcasten/tools/ref_data/human_accelerated_regions/HAQERs/Mangan-Cell2022/other_regions/har.hg38_sorted.autosomes_non_overlapping.bed'

##
for(cl in unique(hs_cre$subclass)) {
    ph <- str_c(cl, '_humanSpec')
    set_name = ph
    tmp <- hs_cre %>% 
        filter(subclass == cl)
   
    cat('\n\n\n\n')
    message('============================================')
    message('Set: ', set_name)
    ##
    outf <- str_c('/wdata/lcasten/sli_wgs/HAQER_enrichment/data/', ph, '.hg38.bed')

    tmp %>% 
        # filter(symbol %in% tmp$gene) %>%
        # filter(! symbol %in% constrained$symbol) %>% ## remove  constrained genes
        # inner_join(chrom_limits) %>%
        # mutate(chromStart = ifelse(chromStart - 35000 >= 0, chromStart - 35000, 1),
            #    chromEnd = ifelse(chromEnd + 10000 <= chrom_size, chromEnd + 10000, chrom_size)) %>% 
        mutate(chr = as.numeric(str_remove_all(`#chrom`, 'chr'))) %>% 
        drop_na() %>%
        arrange(chr, chromStart, chromEnd) %>%
        select(`#chrom`, chromStart, chromEnd) %>% 
        write_tsv(outf, col_names = FALSE)

    ##   
    outf2 = str_c('/wdata/lcasten/sli_wgs/HAQER_enrichment/data/', set_name, '.non_overlapping.hg38.bed')
    cmd <- str_c('bedtools merge -i ', outf, ' > ', outf2)
    system(cmd)

    ## run enrichment: overlapEnrichments method         elements1.lift elements2.lift          noGap.lift                                    out.txt
    outres = str_c("/wdata/lcasten/sli_wgs/HAQER_enrichment/data/", set_name, '.txt')
    outres_har <- str_c("/wdata/lcasten/sli_wgs/HAQER_enrichment/data/", set_name, '.HAR_enrichment.txt')
    ## HAQER enrichment
    cmd <- str_c('~/go/bin/overlapEnrichments normalApproximate ', outf2, ' ', outbed, ' /wdata/lcasten/tools/ref_data/hg38/hg38.chrom.sizes.bed ', outres)
    system(cmd)
    ## HAR enrichment
    cmd2 <- str_c('~/go/bin/overlapEnrichments normalApproximate ', outf2, ' ', outbed_har, ' /wdata/lcasten/tools/ref_data/hg38/hg38.chrom.sizes.bed ', outres_har)
    system(cmd2)
    ## print results
    message('Results HAQERs:')
    system(str_c('cat ', outres, ' | cut -f 9,10,11'))
    message('Results HARs:')
    system(str_c('cat ', outres_har, ' | cut -f 9,10,11'))
}


############################################
## repeat for other types
unique(cre_annot$category)
cons_cre <- cre %>% 
    rename(cCRE = X4) %>% 
    rename(`#chrom` = X1, chromStart = X2, chromEnd = X3) %>% 
    inner_join(cre_annot) %>% 
    filter(category == 'CA_Cons')

##
outbed = '/wdata/lcasten/tools/ref_data/human_accelerated_regions/HAQERs/Mangan-Cell2022/HAQER.hg38.sorted_autosomes_non_overlapping.bed'
outbed_har <- '/wdata/lcasten/tools/ref_data/human_accelerated_regions/HAQERs/Mangan-Cell2022/other_regions/har.hg38_sorted.autosomes_non_overlapping.bed'

##
for(cl in unique(cons_cre$subclass)) {
    ph <- str_c(cl, '_CA_Cons')
    set_name = ph
    tmp <- cons_cre %>% 
        filter(subclass == cl)
   
    cat('\n\n\n\n')
    message('============================================')
    message('Set: ', set_name)
    ##
    outf <- str_c('/wdata/lcasten/sli_wgs/HAQER_enrichment/data/', ph, '.hg38.bed')

    tmp %>% 
        # filter(symbol %in% tmp$gene) %>%
        # filter(! symbol %in% constrained$symbol) %>% ## remove  constrained genes
        # inner_join(chrom_limits) %>%
        # mutate(chromStart = ifelse(chromStart - 35000 >= 0, chromStart - 35000, 1),
            #    chromEnd = ifelse(chromEnd + 10000 <= chrom_size, chromEnd + 10000, chrom_size)) %>% 
        mutate(chr = as.numeric(str_remove_all(`#chrom`, 'chr'))) %>% 
        drop_na() %>%
        arrange(chr, chromStart, chromEnd) %>%
        select(`#chrom`, chromStart, chromEnd) %>% 
        write_tsv(outf, col_names = FALSE)

    ##   
    outf2 = str_c('/wdata/lcasten/sli_wgs/HAQER_enrichment/data/', set_name, '.non_overlapping.hg38.bed')
    cmd <- str_c('bedtools merge -i ', outf, ' > ', outf2)
    system(cmd)

    ## run enrichment: overlapEnrichments method         elements1.lift elements2.lift          noGap.lift                                    out.txt
    outres = str_c("/wdata/lcasten/sli_wgs/HAQER_enrichment/data/", set_name, '.txt')
    outres_har <- str_c("/wdata/lcasten/sli_wgs/HAQER_enrichment/data/", set_name, '.HAR_enrichment.txt')
    ## HAQER enrichment
    cmd <- str_c('~/go/bin/overlapEnrichments normalApproximate ', outf2, ' ', outbed, ' /wdata/lcasten/tools/ref_data/hg38/hg38.chrom.sizes.bed ', outres)
    system(cmd)
    ## HAR enrichment
    cmd2 <- str_c('~/go/bin/overlapEnrichments normalApproximate ', outf2, ' ', outbed_har, ' /wdata/lcasten/tools/ref_data/hg38/hg38.chrom.sizes.bed ', outres_har)
    system(cmd2)
    ## print results
    message('Results HAQERs:')
    system(str_c('cat ', outres, ' | cut -f 9,10,11'))
    message('Results HARs:')
    system(str_c('cat ', outres_har, ' | cut -f 9,10,11'))
}

############################################
## repeat for other types: CA diverged
unique(cre_annot$category)
div_cre <- cre %>% 
    rename(cCRE = X4) %>% 
    rename(`#chrom` = X1, chromStart = X2, chromEnd = X3) %>% 
    inner_join(cre_annot) %>% 
    filter(category == 'CA_Diver')

##
outbed = '/wdata/lcasten/tools/ref_data/human_accelerated_regions/HAQERs/Mangan-Cell2022/HAQER.hg38.sorted_autosomes_non_overlapping.bed'
outbed_har <- '/wdata/lcasten/tools/ref_data/human_accelerated_regions/HAQERs/Mangan-Cell2022/other_regions/har.hg38_sorted.autosomes_non_overlapping.bed'

##
for(cl in unique(div_cre$subclass)) {
    ph <- str_c(cl, '_CA_Diver')
    set_name = ph
    tmp <- div_cre %>% 
        filter(subclass == cl)
   
    cat('\n\n\n\n')
    message('============================================')
    message('Set: ', set_name)
    ##
    outf <- str_c('/wdata/lcasten/sli_wgs/HAQER_enrichment/data/', ph, '.hg38.bed')

    tmp %>% 
        # filter(symbol %in% tmp$gene) %>%
        # filter(! symbol %in% constrained$symbol) %>% ## remove  constrained genes
        # inner_join(chrom_limits) %>%
        # mutate(chromStart = ifelse(chromStart - 35000 >= 0, chromStart - 35000, 1),
            #    chromEnd = ifelse(chromEnd + 10000 <= chrom_size, chromEnd + 10000, chrom_size)) %>% 
        mutate(chr = as.numeric(str_remove_all(`#chrom`, 'chr'))) %>% 
        drop_na() %>%
        arrange(chr, chromStart, chromEnd) %>%
        select(`#chrom`, chromStart, chromEnd) %>% 
        write_tsv(outf, col_names = FALSE)

    ##   
    outf2 = str_c('/wdata/lcasten/sli_wgs/HAQER_enrichment/data/', set_name, '.non_overlapping.hg38.bed')
    cmd <- str_c('bedtools merge -i ', outf, ' > ', outf2)
    system(cmd)

    ## run enrichment: overlapEnrichments method         elements1.lift elements2.lift          noGap.lift                                    out.txt
    outres = str_c("/wdata/lcasten/sli_wgs/HAQER_enrichment/data/", set_name, '.txt')
    outres_har <- str_c("/wdata/lcasten/sli_wgs/HAQER_enrichment/data/", set_name, '.HAR_enrichment.txt')
    ## HAQER enrichment
    cmd <- str_c('~/go/bin/overlapEnrichments normalApproximate ', outf2, ' ', outbed, ' /wdata/lcasten/tools/ref_data/hg38/hg38.chrom.sizes.bed ', outres)
    system(cmd)
    ## HAR enrichment
    cmd2 <- str_c('~/go/bin/overlapEnrichments normalApproximate ', outf2, ' ', outbed_har, ' /wdata/lcasten/tools/ref_data/hg38/hg38.chrom.sizes.bed ', outres_har)
    system(cmd2)
    ## print results
    message('Results HAQERs:')
    system(str_c('cat ', outres, ' | cut -f 9,10,11'))
    message('Results HARs:')
    system(str_c('cat ', outres_har, ' | cut -f 9,10,11'))
}


