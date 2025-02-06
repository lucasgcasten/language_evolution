##
library(tidyverse)

#########################
## read in data
#########################
## sample info
md <- read_csv('/wdata/lcasten/tools/ref_data/neanderthal_ancestry_through_time-Iasi-Science2024/doi_10_5061_dryad_zw3r228gg__v20241126/Meta_Data_individuals.csv')
md
md2 <- md %>% 
    select(sample = sample_name, superpopulation, ML_BP_Mean)
md2 %>% 
    filter(ML_BP_Mean > 0)
unique(md2$superpopulation)

## regions in humans with neanderthal introgression
nean_anc <- read_csv('/wdata/lcasten/tools/ref_data/neanderthal_ancestry_through_time-Iasi-Science2024/doi_10_5061_dryad_zw3r228gg__v20241126/Neandertal_segments_matching_references_Shared_map.csv')
names(nean_anc)
hist(nean_anc$score)
nean_anc %>% 
    filter(score >= 4.5) %>%
    select(chrom, start, end)

nean_anc2 <- nean_anc %>% 
    filter(score >= 4.5) %>%
    select(1:4, sample, matches('prop_match'))
table(nean_anc$sample)

## neanderthal regions from supplemental of Iasi paper
readxl::excel_sheets('/wdata/lcasten/tools/ref_data/neanderthal_ancestry_through_time-Iasi-Science2024/science.adq3010_supplementary_tables.xlsx')
nean_anc3 <- readxl::read_xlsx('/wdata/lcasten/tools/ref_data/neanderthal_ancestry_through_time-Iasi-Science2024/science.adq3010_supplementary_tables.xlsx', sheet = 'table S26 (All high frequency N', skip = 1)
nean_anc3 <- nean_anc3 %>% 
    select(`#chrom` = chr, chromStart = `Start_pos (bp)`, chromEnd = `End_pos (bp)`) %>% 
    mutate(type = 'supplemental')
nean_anc4 <- readxl::read_xlsx('/wdata/lcasten/tools/ref_data/neanderthal_ancestry_through_time-Iasi-Science2024/science.adq3010_supplementary_tables.xlsx', sheet = 'table S27 (Overlap with perviou', skip = 1)
nean_anc4 <- nean_anc4 %>% 
    select(`#chrom` = chr, chromStart = start, chromEnd = end) %>% 
    mutate(type = 'supplemental_rep') %>% 
    drop_na() %>% 
    mutate(`#chrom` = as.numeric(`#chrom`)) %>% 
    arrange(`#chrom`, chromStart)

## merge sample info and introgression data to identify regions in ancient and modern europeans
anc_intro <- nean_anc2 %>% 
    inner_join(md2) %>% 
    # filter(superpopulation %in% c('WEurAs', 'EUR')) %>% 
    filter(ML_BP_Mean > 0) %>% 
    rename(`#chrom` = chrom, chromStart = start, chromEnd = end) %>% 
    mutate(type = 'ancient')
        
modern_intro <- nean_anc2 %>% 
    inner_join(md2) %>% 
    # filter(superpopulation %in% c('WEurAs', 'EUR')) %>% 
    filter(ML_BP_Mean == 0) %>% 
    rename(`#chrom` = chrom, chromStart = start, chromEnd = end) %>% 
    mutate(type = 'modern')

all_eur <- nean_anc2 %>% 
    inner_join(md2) %>% 
    # filter(superpopulation %in% c('WEurAs', 'EUR')) %>% 
    rename(`#chrom` = chrom, chromStart = start, chromEnd = end) %>% 
    mutate(type = 'all_eur')

mean(anc_intro$chromEnd - anc_intro$chromStart)
mean(modern_intro$chromEnd - modern_intro$chromStart)
mean(nean_anc3$chromEnd - nean_anc3$chromStart)

#################
## look at other lists of introgressed variants
##################
# introg2 <- read_tsv('/wdata/lcasten/tools/ref_data/archaic_hominini/neanderthal_1000genomes_reversions/Rinker_Suppl.File1/Rinker2018_Suppl_File_1/NDA_introgressed_EUR.bed', col_names = FALSE)

# introg_snp <- introg2 %>% 
#     distinct() %>% 
#     mutate(`#chrom` = as.numeric(str_remove_all(X1, pattern = 'chr'))) %>%
#     select(`#chrom`, chromStart = X2, chromEnd = X3) %>% 
#     mutate(type = 'NDA_introgressed_SNPs_1000Geur_Rinker2018') %>% 
#     arrange(`#chrom`, chromStart)

# introg_hap <- introg2 %>% 
#     select(X4) %>% 
#     distinct() %>% 
#     mutate(`#chrom` = str_split(X4, pattern = '_', simplify = TRUE)[,1],
#            `#chrom` = as.numeric(str_remove_all(`#chrom`, pattern = 'chr')),
#            chromStart = str_split(X4, pattern = '_', simplify = TRUE)[,2],
#            chromEnd = str_split(X4, pattern = '_', simplify = TRUE)[,3]) %>%
#     mutate(chromStart = as.numeric(chromStart),
#            chromEnd = as.numeric(chromEnd)) %>%
#     select(`#chrom`, chromStart, chromEnd) %>% 
#     arrange(`#chrom`, chromStart) %>%
#     mutate(type = 'NDA_introgressed_haplotypes_1000Geur_Rinker2018')

introg2 <- read_tsv('/wdata/lcasten/tools/ref_data/archaic_hominini/neanderthal_1000genomes_reversions/Rinker_Suppl.File1/Rinker2018_Suppl_File_1/all_introgressed_EUR.bed')

introg_snp <- introg2 %>% 
    distinct() %>% 
    mutate(`#chrom` = as.numeric(str_remove_all(`#chr1`, pattern = 'chr'))) %>%
    select(`#chrom`, chromStart = `POS-1`, chromEnd = POS) %>% 
    mutate(type = 'all_introgressed_SNPs_Rinker2018') %>% 
    arrange(`#chrom`, chromStart)

introg_hap <- introg2 %>% 
    select(X4 = V16_LD_haplotype) %>% 
    distinct() %>% 
    mutate(`#chrom` = str_split(X4, pattern = '_', simplify = TRUE)[,1],
           `#chrom` = as.numeric(str_remove_all(`#chrom`, pattern = 'chr')),
           chromStart = str_split(X4, pattern = '_', simplify = TRUE)[,2],
           chromEnd = str_split(X4, pattern = '_', simplify = TRUE)[,3]) %>%
    mutate(chromStart = as.numeric(chromStart),
           chromEnd = as.numeric(chromEnd)) %>%
    select(`#chrom`, chromStart, chromEnd) %>% 
    arrange(`#chrom`, chromStart) %>%
    mutate(type = 'introgressed_haplotypes_Rinker2018')

##
ng <- readxl::read_xlsx('/wdata/lcasten/tools/ref_data/human_neurogenesis-Torre-Unieta-Cell-2018/1-s2.0-S0092867417314940-mmc1.xlsx') %>% 
    filter(stat < 0) %>%
    filter(padj < 0.05) %>% 
    select(`#chrom` = Chr, chromStart = peakstart, chromEnd = peakend) %>% 
    mutate(type = 'human_neurogenesis-Torre-Unieta2018') %>% 
    distinct() %>% 
    arrange(`#chrom`, chromStart)

intro_list <- list(ng, anc_intro, modern_intro, all_eur, nean_anc3, nean_anc4, introg_snp, introg_hap)
str(intro_list, 2)
intro_list[[1]]

######################################
## try enrichment analysis in HAQERs
######################################
##
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

##
for(i in 1:length(intro_list)) {
    tmp <- intro_list[[i]]

    ph <- unique(tmp$type)
    set_name = str_c('neanderthal_introgression.', ph)

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


##############
## intersect introgressed neanderthal variation with HAQERs and make BED files for ES-PGS
haqer <- read_tsv('/wdata/lcasten/sli_wgs/prs/gene_sets/HAQERs.bed', col_names = FALSE)
haqer
# list.dirs('/wdata/lcasten/sli_wgs/prs', recursive = FALSE)
intro_reg <- all_eur %>% 
    select(chromosome = `#chrom`, start = chromStart, end = chromEnd) %>% 
    distinct() 
sum(intro_reg$end - intro_reg$start)

## full HAQERs with neanderthal only chunks
haq_nean_intro <- haqer %>% 
    rename(chromosome = X1, start = X2, end = X3) %>% 
    mutate(chromosome = as.numeric(str_remove_all(chromosome, pattern = 'chr'))) %>%
    fuzzyjoin::genome_inner_join(intro_reg) %>% 
    select(1:3) %>% 
    distinct()

## neanderthal only chunks within HAQERs
nean_intro_haq <- haqer %>% 
    rename(chromosome = X1, start = X2, end = X3) %>% 
    mutate(chromosome = as.numeric(str_remove_all(chromosome, pattern = 'chr'))) %>%
    fuzzyjoin::genome_inner_join(intro_reg) %>% 
    mutate(start = ifelse(start.x > start.y, start.x, start.y),
           end = ifelse(end.x < end.y, end.x, end.y)) %>% 
    select(chromosome = chromosome.x, start, end) %>% 
    distinct() %>% 
    arrange(chromosome, start, end)

intro_reg %>% 
    arrange(chromosome, start, end) %>% 
    mutate(chromosome = str_c('chr', chromosome)) %>%
    write_tsv('/wdata/lcasten/sli_wgs/prs/gene_sets/neanderthal_introgression_Iasi2024.bed', col_names = FALSE)

nean_intro_haq %>% 
    mutate(chromosome = str_c('chr', chromosome)) %>%
    write_tsv('/wdata/lcasten/sli_wgs/prs/gene_sets/neanderthal_introgression_within_HAQERs_Iasi2024.bed', col_names = FALSE)

haq_nean_intro %>% 
    select(chromosome = chromosome.x, start = start.x, end = end.x) %>% 
    arrange(chromosome, start, end) %>% 
    mutate(chromosome = str_c('chr', chromosome)) %>%
    write_tsv('/wdata/lcasten/sli_wgs/prs/gene_sets/HAQERs_with_neanderthal_introgression_Iasi2024.bed', col_names = FALSE)

## =============================================================
## make complementary BED file for introgression evo sets 
## =============================================================
d <- '/Dedicated/jmichaelson-wdata/lcasten/sli_wgs/prs/gene_sets'
chr_sizes <- read_table('/wdata/lcasten/tools/ref_data/hg19/hg19.chrom.sizes', col_names = FALSE) %>% 
    filter(X1 %in% str_c('chr', 1:22))
# chr_sizes %>% 
#     mutate(chr = as.numeric(str_remove_all(X1, pattern = 'chr'))) %>% 
#     arrange(chr) %>%
#     select(-chr) %>% 
#     write_delim(str_c(d, '/chromosome_sizes.txt'), delim = '\t', col_names = FALSE)

## make complement bed file of anno
anno_files <- list.files(d, pattern = 'bed$', full.names = TRUE)
anno_files <- anno_files[str_detect(anno_files, pattern = 'complement', negate = TRUE)]
anno_files <- anno_files[str_detect(anno_files, pattern = 'neanderthal_introgression')]
anno_files
geno <- str_c(d, '/chromosome_sizes.txt')

for(f in anno_files){
    cat(sprintf('\n\n\n\n\n'))
    message('=========================================')
    anno_path <- f
    message(basename(anno_path))
    cf <- str_c(d, '/complement_', basename(anno_path))

    bedtools_cmd <- str_c("bedtools complement -i ", anno_path, " -g ", geno, " > ", cf)
    message(bedtools_cmd)
    system(bedtools_cmd)
    # read_tsv(cf, col_names = FALSE)
}

comp_files <- list.files(d, pattern = 'bed$', full.names = TRUE)
comp_files <- comp_files[str_detect(comp_files, pattern = 'complement')]
comp_files <- comp_files[str_detect(comp_files, pattern = 'HAQER')]
comp_files

############
## plot
files_intro = list.files('/wdata/lcasten/sli_wgs/HAQER_enrichment/data', pattern = 'neanderthal_introgression.all_eur')
files_intro = files_intro[str_detect(files_intro, '.txt$')]
files_intro

res_list <- list()

for(i in 1:length(files_intro)) {
    set <- basename(files_intro[i])
    set = str_remove_all(set, '.txt')
    type = ifelse(str_detect(set, 'HAR_enrichment'), 'HAR', 'HAQER')

    res_list[[i]] <- read_table(str_c('/wdata/lcasten/sli_wgs/HAQER_enrichment/data/', files_intro[i]), show_col_types = FALSE) %>% 
        mutate(ph = set, type = type) %>% 
        relocate(ph, type, EnrichPValue, DepletePValue)
}

#################
introg_enrich = bind_rows(res_list)

##################
p_enrich_introg <- introg_enrich %>% 
    # filter(type == 'HAQER') %>% 
    arrange(EnrichPValue) %>% 
    select(ph, type, matches('PValue')) %>% 
    drop_na() %>%
    mutate(p = ifelse(EnrichPValue < DepletePValue, EnrichPValue, DepletePValue),
           stat = -log10(p),
           stat = ifelse(EnrichPValue < DepletePValue, stat, stat * -1)) %>%
    arrange(type, EnrichPValue) %>% 
    mutate(type = factor(type, levels = c('HAR', 'HAQER'))) %>%
    mutate(type_fill = case_when(p < 0.05 & type == 'HAQER' ~ 'HAQER_sig',
                                 p > 0.05 & type == 'HAQER' ~ 'HAQER_ns',
                                 p < 0.05 & type == 'HAR' ~ 'HAR_sig',
                                 p > 0.05 & type == 'HAR' ~ 'HAR_ns')) %>%
    mutate(sig = ifelse(p < 0.05, 1, 0.8)) %>%
    # mutate(EnrichPValue = ifelse(EnrichPValue >= 0.9 & str_detect(evo_type_clean, 'Human'), 0.9, EnrichPValue),
    #        EnrichPValue = ifelse(EnrichPValue >= 0.1 & str_detect(evo_type_clean, 'Human', negate = TRUE), 0.1, EnrichPValue)) %>%
    arrange(p) %>% # filter(str_detect(evo_type_clean, 'Conserv') & type == 'HAQER')
    ggplot(aes(x = -log10(EnrichPValue), y = type, fill = type)) + # , alpha = sig
    geom_bar(stat = 'identity', position = 'dodge', width = 0.3) +
    xlab('-log10(enrichment p-value)') +
    ylab(NULL) +
    geom_vline(xintercept = -log10(0.05), color = 'grey', linetype = 'dashed', size = 1.075) +
    # geom_vline(xintercept = -1 * -log10(0.05), color = 'grey', linetype = 'dashed', size = 1.075) +
    # facet_grid(cols = vars(evo_type_clean), scales = 'free') + # , space = 'free'
    scale_fill_manual(values = c('goldenrod', 'purple2')) +
    labs(fill = NULL) +
    theme_classic() +
    labs(fill = NULL) +
    theme(legend.position = 'none',
          axis.text = element_text(size = 14),
          axis.title = element_text(size = 16),
          strip.text = element_text(size = 16),
          legend.text = element_text(size = 16),
          plot.title = element_text(size = 18, hjust = 0.5),
          panel.border = element_rect(colour = "black", fill=NA)) +
    ggtitle('Neanderthal introgressed regions')


##############################################
##############################################
## compute neanderthal ES-PGS 
##############################################
##############################################
