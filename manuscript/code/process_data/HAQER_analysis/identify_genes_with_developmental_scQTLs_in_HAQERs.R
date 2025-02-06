library(tidyverse)

## gene list
gene_info <- read_tsv('/wdata/lcasten/tools/ref_data/hg19/gene_map.tsv')
    filter(str_detect(`#chrom`, '[0-9]'))

## haqer bed file
haqer_loc <- read_tsv('/wdata/lcasten/sli_wgs/prs/gene_sets/HAQERs.bed', col_names = FALSE)
names(haqer_loc) <- c('chromosome', 'start', 'end')

## HAQERs w/ SNP
haqer_annot <- read_csv('/wdata/lcasten/sli_wgs/prs/pathway_prs/HAQER_PGS_manual_calculation_annotations.csv') %>% 
    rename(start = haqer_start, end = haqer_end) %>% 
    mutate(type = 'HAQER') %>% 
    select(chromosome, start, end, haqer) %>% 
    distinct() %>% 
    mutate(chromosome = str_c('chr', chromosome))

## subset HAQERs to those with CP-GWAS SNPs
haqer_loc <- haqer_loc %>% 
    inner_join(haqer_annot) %>% 
    select(-haqer)

### identify genes regulated by HAQERs during neurodevelopment
files = list.files('/wdata/lcasten/tools/ref_data/single_cell_QTL-Jerber-NatureGenetics-2021/eqtl_summary_stats_filtered', pattern = '.txt$', full.names = TRUE)

qtl_list = list()

for(f in files) {
    ph = basename(f)
    ph = str_split(ph, pattern = '[.]qtl_', simplify = TRUE)[,1]
    ##
    genes <- read_table(f, show_col_types = FALSE) %>% 
        select(feature_id, snp_chromosome, snp_position) %>% 
        distinct() %>% 
        rename(chromosome = snp_chromosome, end = snp_position) %>% 
        mutate(start = end - 1) %>% 
        relocate(chromosome, start, end)

    ##
    goi <- haqer_loc %>% 
        mutate(chromosome = as.numeric(str_remove_all(chromosome, 'chr'))) %>%
        fuzzyjoin::genome_inner_join(genes) %>% 
        distinct(feature_id, .keep_all = TRUE) %>% 
        mutate(gs = ph)
    qtl_list[[ph]] <- goi
    # gene_info %>% 
    #     filter(ensembl_gene_id %in% goi$feature_id) %>% 
    #     select(symbol) %>% 
    #     unlist() %>% 
    #     unname()
    outf <- str_c('/wdata/lcasten/sli_wgs/HAQER_enrichment/data/neurodevelopmental_scQTL_genes/', ph, '.genelist.txt')
    goi %>% 
        select(feature_id) %>% 
        write_tsv(outf, col_names = FALSE)
}

bind_rows(qtl_list) %>% 
    distinct(feature_id) %>% 
    write_tsv('/wdata/lcasten/sli_wgs/HAQER_enrichment/data/neurodevelopmental_scQTL_genes/all_genes.genelist.txt', col_names = FALSE)

bind_rows(qtl_list) %>% 
    relocate(gs) %>% 
    mutate(gs = str_split(gs, pattern = '[.]', simplify = TRUE)[,1]) %>% 
    distinct(gs, feature_id) %>% 
    group_by(feature_id) %>%
    count() %>% 
    filter(n > 1) %>% 
    ungroup() %>% 
    select(feature_id) %>% 
    write_tsv('/wdata/lcasten/sli_wgs/HAQER_enrichment/data/neurodevelopmental_scQTL_genes/all_genes_at_multiple_timepoints.genelist.txt', col_names = FALSE)

bind_rows(qtl_list) %>% 
    relocate(gs) %>% 
    mutate(gs = str_split(gs, pattern = '[.]', simplify = TRUE)[,1]) %>% 
    distinct(gs, feature_id) %>% 
    group_by(feature_id) %>%
    count() %>% 
    filter(n == 3) %>% 
    ungroup() %>% 
    select(feature_id) %>% 
    write_tsv('/wdata/lcasten/sli_wgs/HAQER_enrichment/data/neurodevelopmental_scQTL_genes/all_genes_at_all_timepoints.genelist.txt', col_names = FALSE)


######
## run genes through FUMA
fuma_res <- read_tsv('/wdata/lcasten/sli_wgs/HAQER_enrichment/data/neurodevelopmental_scQTL_genes/FUMA_gene2func566270/GS.txt')
unique(fuma_res$Category)

fuma_res %>% 
    filter(Category == 'GWAScatalog') %>%
    select(2, p, adjP) %>% 
    as.data.frame()

fuma_res %>% 
    filter(Category == 'Reactome') %>%
    select(2, p, adjP) 

fuma_res %>% 
    filter(Category == 'Wikipathways') %>%
    select(2, p, adjP) 

fuma_res %>% 
    filter(Category == 'KEGG') %>%
    select(2, p, adjP) 

fuma_res %>% 
    filter(Category == 'GO_bp') %>%
    select(2, p, adjP)  %>% 
    # select(1) %>% as.data.frame()
    head(n = 20)

fuma_res %>% 
    filter(Category == 'GO_cc') %>%
    select(2, p, adjP)  %>% 
    # select(1) %>% as.data.frame()
    head(n = 20)