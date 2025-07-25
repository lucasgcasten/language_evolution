library(tidverse)

## Konopka 2023 evo, human specific brain expression
d <- '/Dedicated/jmichaelson-wdata/lcasten/sli_wgs/prs/gene_sets'
gene_info <- read_tsv('/wdata/lcasten/tools/ref_data/hg19/hg19_hgnc.tsv') %>% 
    select(1:3, symbol)
hs_genes <- readxl::read_xlsx('/wdata/lcasten/tools/ref_data/evolution/human_brain_cellular_complexity_evo-Caglayan-Nature2023/41586_2023_6338_MOESM5_ESM.xlsx')
table(hs_genes$Evolution)
hs_genes <- hs_genes %>%
    filter(Evolution == 'Human_Specific')
table(hs_genes$CellType)
hs_genes <- hs_genes %>% 
    mutate(CellType = str_replace_all(CellType, pattern = '[-]$', replacement = 'negative'),
           CellType = str_replace_all(CellType, pattern = '[+]$', replacement = 'positive'))
hs_genes_unique <- hs_genes %>% # filter(Gene == 'FOXP2')
    distinct(Gene)
vars = unique(hs_genes$CellType)
for(i in vars) {
    tmp <- hs_genes %>% 
        filter(CellType == i) %>% 
        rename(gene = Gene)
    tmp2 <- gene_info %>% 
        filter(symbol %in% tmp$gene) %>% 
        filter(`#chrom` %in% str_c('chr', 1:22)) %>%
        mutate(chromStart = chromStart - 35000,
               chromEnd = chromEnd + 10000) %>%
        mutate(chromStart = ifelse(chromStart <= 1, 1, chromStart))
    outf <- str_c(d, '/human_specific_evolution_brain_expression.', i, '.bed')
    tmp2 %>% 
        write_tsv(outf, col_names = FALSE)
}

tmp2 <- gene_info %>% 
    filter(symbol %in% hs_genes_unique$Gene) %>% 
    filter(`#chrom` %in% str_c('chr', 1:22)) %>%
    mutate(chromStart = chromStart - 35000,
            chromEnd = chromEnd + 10000) %>%
    mutate(chromStart = ifelse(chromStart <= 1, 1, chromStart))
outf <- str_c(d, '/human_specific_evolution_brain_expression.all_cell_types.bed')
tmp2 %>% 
    write_tsv(outf, col_names = FALSE)

## gene set for excitatory neurons interesting from FOXP2 perspective
foxp2_genes <- hs_genes %>% 
    filter(CellType %in% c('L5-6_THEMIS_1', 'L4-6_RORB_2')) %>% 
    distinct(Gene)
tmp2 <- gene_info %>% 
    filter(symbol %in% foxp2_genes$Gene) %>% 
    filter(`#chrom` %in% str_c('chr', 1:22)) %>%
    mutate(chromStart = chromStart - 35000,
            chromEnd = chromEnd + 10000) %>%
    mutate(chromStart = ifelse(chromStart <= 1, 1, chromStart)) %>% 
    select(-symbol)
outf <- str_c(d, '/human_specific_evolution_brain_expression_FOXP2_DEG.bed')
tmp2 %>% 
    write_tsv(outf, col_names = FALSE)

non_foxp2_genes <- hs_genes %>% 
    filter(! CellType %in% c('L5-6_THEMIS_1', 'L4-6_RORB_2')) %>% 
    filter(str_detect(CellType, pattern = '^L[0-9]')) %>%
    filter(! Gene %in% foxp2_genes$Gene) %>%
    distinct(Gene)
tmp2 <- gene_info %>% 
    filter(symbol %in% non_foxp2_genes$Gene) %>% 
    filter(`#chrom` %in% str_c('chr', 1:22)) %>%
    mutate(chromStart = chromStart - 35000,
            chromEnd = chromEnd + 10000) %>%
    mutate(chromStart = ifelse(chromStart <= 1, 1, chromStart)) %>% 
    select(-symbol)
outf <- str_c(d, '/human_specific_evolution_brain_expression_non_FOXP2_excitatory.bed')
tmp2 %>% 
    write_tsv(outf, col_names = FALSE)

## gene set for excitatory neurons interesting from FOXP2 + downstream targets perspective
foxp2_target_genes <- hs_genes %>% 
    filter(CellType %in% c('L5-6_THEMIS_1', 'L5-6_THEMIS_2', 'L4-6_RORB_2', 'L3-5_RORB_1', 'L3-5_RORB_2')) %>% 
    distinct(Gene)
tmp2 <- gene_info %>% 
    filter(symbol %in% foxp2_target_genes$Gene) %>% 
    filter(`#chrom` %in% str_c('chr', 1:22)) %>%
    mutate(chromStart = chromStart - 35000,
            chromEnd = chromEnd + 10000) %>%
    mutate(chromStart = ifelse(chromStart <= 1, 1, chromStart)) %>% 
    select(-symbol)
outf <- str_c(d, '/human_specific_evolution_brain_expression_FOXP2_targets.bed')
tmp2 %>% 
    write_tsv(outf, col_names = FALSE)

## gene set for excitatory neurons interesting from downstream targets of FOXPw only
foxp2_target_genes_only <- hs_genes %>% 
    filter(CellType %in% c('L4-6_RORB_2', 'L3-5_RORB_1', 'L3-5_RORB_2')) %>% 
    distinct(Gene)
tmp2 <- gene_info %>% 
    filter(symbol %in% foxp2_target_genes_only$Gene) %>% 
    filter(`#chrom` %in% str_c('chr', 1:22)) %>%
    mutate(chromStart = chromStart - 35000,
            chromEnd = chromEnd + 10000) %>%
    mutate(chromStart = ifelse(chromStart <= 1, 1, chromStart)) %>% 
    select(-symbol)
outf <- str_c(d, '/human_specific_evolution_brain_expression_FOXP2_targets_only.bed')
tmp2 %>% 
    write_tsv(outf, col_names = FALSE)


## updated list of HS variants overlapping brain CRE
tmp <- readxl::read_xlsx('/wdata/lcasten/tools/ref_data/evolution/human_brain_cellular_complexity_evo-Caglayan-Nature2023/41586_2023_6338_MOESM10_ESM.xlsx')
tmp %>% 
    filter(Is_HS == 'HS') %>%
    select(Modern_Chr, Modern_Pos) %>% 
    distinct() %>% 
    mutate(id = str_c(Modern_Chr, ':', Modern_Pos, '-', Modern_Pos)) %>% 
    select(id) %>%
    write_tsv('/wdata/lcasten/tools/ref_data/evolution/human_brain_cellular_complexity_evo-Caglayan-Nature2023/human_specific_variants_in_brain_CRE.hg38.txt', col_names = FALSE)
read_tsv('/wdata/lcasten/tools/ref_data/evolution/human_brain_cellular_complexity_evo-Caglayan-Nature2023/human_specific_variants_in_brain_CRE.hg19.bed', col_names = FALSE) %>% 
    mutate(chr = str_split(X1, pattern = '[:]', simplify = TRUE)[,1],
           start = str_split(X1, pattern = '[:]', simplify = TRUE)[,2],
           start = as.numeric(str_split(start, pattern = '[-]', simplify = TRUE)[,1]) - 1,
           stop = start + 1) %>% 
    select(chr, start, stop) %>% 
    mutate(chr2 = as.numeric(str_remove_all(chr, pattern = 'chr'))) %>% 
    arrange(chr2, start) %>% 
    select(-chr2) %>%
    write_tsv(str_c(d, '/human_specific_variants_in_brain_CRE.hg19.bed'), col_names = FALSE)

## HS chromatin
tmp <- readxl::read_xlsx('/wdata/lcasten/tools/ref_data/evolution/human_brain_cellular_complexity_evo-Caglayan-Nature2023/41586_2023_6338_MOESM6_ESM.xlsx') %>% 
    filter(Evolution == 'Human_Specific')

tmp
table(tmp$CellType)
tmp %>% 
    select(Gene, CellType) %>% 
    mutate(chr = str_split(Gene, pattern = '[:]', simplify = TRUE)[,1],
           Gene = str_split(Gene, pattern = '[:]', simplify = TRUE)[,2],
           start = str_split(Gene, pattern = '[-]', simplify = TRUE)[,1],
           start = as.numeric(start) - 1,
           end = str_split(Gene, pattern = '[-]', simplify = TRUE)[,2]) %>% 
    select(chr, start, end, CellType) %>% 
    write_tsv('/wdata/lcasten/tools/ref_data/evolution/human_brain_cellular_complexity_evo-Caglayan-Nature2023/human_specific_brain_CRE.hg38.bed', col_names = FALSE)

read_tsv('/wdata/lcasten/tools/ref_data/evolution/human_brain_cellular_complexity_evo-Caglayan-Nature2023/human_specific_brain_CRE.hg19.bed', col_names = FALSE) %>% 
    mutate(chr2 = as.numeric(str_remove_all(X1, pattern = 'chr'))) %>% 
    drop_na() %>%
    arrange(X4, chr2, X2, X3) %>% 
    select(-chr2) %>%
    write_tsv(str_c(d, '/human_specific_brain_CRE_sorted.hg19.bed'), col_names = FALSE)
bcre <- read_tsv(str_c(d, '/human_specific_brain_CRE_sorted.hg19.bed'), col_names = FALSE)
unique(bcre$X1)
ct <- unique(bcre$X4)
ct
for(c in ct) {
    tmp <- bcre %>% 
        filter(X4 == c)
    c2 <- str_replace_all(c, pattern = '[-]$', replacement = 'negative')
    c2 <- str_replace_all(c2, pattern = '[+]$', replacement = 'positive')
    message('There are ', nrow(tmp), ' brain-CREs for ', c2)
    outf <- str_c(d, '/human_specific_brain_CRE_', c2, '.bed')
    tmp %>% 
        select(-X4) %>%
        write_tsv(outf, col_names = FALSE)
}

## add DARs to DEGs 
tmp2 
foxp2_targets_only_dars <- bcre %>% 
    filter(X4 %in% c('L4-6_RORB_2', 'L3-5_RORB_1', 'L3-5_RORB_2')) %>% 
    select(-X4) %>% 
    distinct()
names(foxp2_targets_only_dars) <- names(tmp2)
foxp2_targets_only_all <- bind_rows(tmp2, foxp2_targets_only_dars) %>% 
    mutate(chr = as.numeric(str_remove_all(`#chrom`, 'chr'))) %>%
    arrange(chr, chromStart, chromEnd) %>% 
    select(-chr)

# make string to paste in PRSet
s <- c()
s2 <- c()
cs <- c()
cs2 <- c()
i = 0
for(c in ct) {
    i = i + 1
    tmp <- bcre %>% 
        filter(X4 == c)
    c2 <- str_replace_all(c, pattern = '[-]$', replacement = 'negative')
    c2 <- str_replace_all(c2, pattern = '[+]$', replacement = 'positive')
    s[i] <- str_c(d, '/human_specific_brain_CRE_', c2, '.bed')
    s2[i] <- str_c('human_specific_brain_CRE_', c2)

    cs[i] <- str_c(d, '/complement_human_specific_brain_CRE_', c2, '.bed')
    cs2[i] <- str_c('complement_human_specific_brain_CRE_', c2)
}
s
s2
s_cre <- str_c(s, ':', s2)
cs_cre <- str_c(cs, ':', cs2)
s_cre_str <- str_c(s_cre, collapse = ',')
cs_cre_str <- str_c(cs_cre, collapse = ',')

## primate specific information
psi <- readxl::read_xlsx('/wdata/lcasten/tools/ref_data/evolution/primate_specific_information_regions-Wei2023FrontBioinf/Table1.XLSX', skip = 2, guess_max = 20000)
psi %>% 
    mutate(chr = str_c('chr', chr)) %>%
    arrange(chr, start, end) %>% 
    select(chr, start, end) %>%
    drop_na() %>%
    write_tsv('/wdata/lcasten/sli_wgs/prs/gene_sets/primate_specific_information.bed', col_names = FALSE)

##
linar <- read_tsv('/wdata/lcasten/tools/ref_data/evolution/primate_lineage_accelerated_regions-Bi2023SciAdv/reformatted_lineages_with_humans.hg19.bed', col_names = FALSE) %>% 
    mutate(chr = as.numeric(str_remove_all(X1, pattern = 'chr'))) %>% 
    drop_na() %>% 
    arrange(X4, chr, X2)
sp = unique(linar$X4)
for(s in sp) {
    linar %>% 
        filter(X4 == s) %>% 
        select(X1, X2, X3) %>% 
        write_tsv(str_c('/wdata/lcasten/sli_wgs/prs/gene_sets/LinAR_', s, '.bed'), col_names = FALSE)
}

## get top 5% of SDS
sds <- read_table('/wdata/lcasten/tools/ref_data/singleton_density_scores_Science2016/SDS_UK10K_n3195_release_Sep_19_2016.tab.gz')
sds <- sds %>% 
    arrange(desc(abs(SDS)))

# remove MHC (chr6:28,477,797-33,448,354) and lactase regions (chr2:136,545,415–136,594,750)
sds <- sds %>% 
    mutate(drop = case_when(CHR == 2 & POS <= 136594750 & POS >= 136545415 ~ TRUE,
                            CHR == 6 & POS >= 28477797 & POS <= 33448354 ~ TRUE,
                            TRUE ~ FALSE)) %>% 
    filter(drop == FALSE)
.05 * nrow(sds)

sds_top5pct <- sds %>% 
    mutate(rank = rank(abs(SDS))) %>% 
    filter(rank >= .95 * max(rank))

sds_top5pct %>% 
    mutate(start = POS - 1) %>%
    select(CHR, start, end = POS) %>% 
    arrange(CHR, start) %>% 
    mutate(CHR = str_c('chr', CHR)) %>%
    write_tsv('/wdata/lcasten/sli_wgs/prs/gene_sets/singleton_density_score_top5pct.bed', col_names = FALSE)


## get primate UCE regions (need to convert to hg19)
## liftOver format: (e.g. "chr4 100000 100001", 0-based)  ------ or the format of the position box ("chr4:100,001-100,001", 1-based)
primate_uce <- readxl::read_xlsx('/wdata/lcasten/tools/ref_data/primate_conservation_Nature-Kuderna2024/41586_2023_6798_MOESM3_ESM.xlsx', sheet = 5)
primate_uce %>% 
    mutate(chr = str_remove_all(chromosome, pattern = 'chr'),
           chr = as.numeric(chr)) %>% 
    drop_na() %>%
    arrange(chr, start, stop) %>% 
    select(-chr) %>% 
    write_tsv('/wdata/lcasten/tools/ref_data/primate_conservation_Nature-Kuderna2024/primate_UCE.hg38.bed', col_names = FALSE)
## liftOver on UCSC to hg19
read_tsv('/wdata/lcasten/tools/ref_data/primate_conservation_Nature-Kuderna2024/primate_UCE.hg19.bed', col_names = FALSE) %>% 
    mutate(chr = str_remove_all(X1, pattern = 'chr'),
           chr = as.numeric(chr)) %>% 
    drop_na() %>%
    arrange(chr, X2, X3) %>% 
    select(-chr) %>% 
    filter(X2 < X3) %>%
    write_tsv('/wdata/lcasten/sli_wgs/prs/gene_sets/primate_UCE.hg19.bed', col_names = FALSE)


## primate exon conservation
primate_exon <- readxl::read_xlsx('/wdata/lcasten/tools/ref_data/primate_conservation_Nature-Kuderna2024/41586_2023_6798_MOESM3_ESM.xlsx', sheet = 2)
primate_exon %>% 
    mutate(chr = str_remove_all(chromosome, pattern = 'chr'),
           chr = as.numeric(chr)) %>% 
    drop_na() %>%
    arrange(chr, start, end) %>% 
    select(-chr) %>% 
    filter(start < end) %>%
    # filter(is_primate_specific_constrained_gene == 1)
    write_tsv('/wdata/lcasten/tools/ref_data/primate_conservation_Nature-Kuderna2024/primate_conserved_exons.hg38.bed', col_names = FALSE)
read_tsv('/wdata/lcasten/tools/ref_data/primate_conservation_Nature-Kuderna2024/primate_conserved_exons.hg19.bed', col_names = FALSE) %>% 
    mutate(chr = str_remove_all(X1, pattern = 'chr'),
           chr = as.numeric(chr)) %>% 
    drop_na() %>%
    arrange(chr, X2, X3) %>% 
    select(-chr) %>% 
    select(1:3) %>%
    filter(X2 < X3) %>%
    write_tsv('/wdata/lcasten/sli_wgs/prs/gene_sets/primate_conserved_exons.hg19.bed', col_names = FALSE)


## get conserved regions across distant evo time
vertebrate_cons <- read_table('/wdata/lcasten/tools/ref_data/evolution/phastConsElements46way.txt.gz', col_names = FALSE)
primate_cons <- read_table('/wdata/lcasten/tools/ref_data/evolution/phastConsElements46wayPrimates.txt.gz', col_names = FALSE)
mammal_cons <- read_table('/wdata/lcasten/tools/ref_data/evolution/phastConsElements46wayPlacental.UCSC.txt.gz', col_names = FALSE)

## example BED file for PRSice
read_tsv('/wdata/lcasten/sli_wgs/prs/gene_sets/HAQERs.bed', col_names = FALSE)

vertebrate_cons %>% 
    mutate(chr = str_remove_all(X2, pattern = 'chr'),
            chr = as.numeric(chr)) %>% 
    drop_na() %>%
    arrange(chr, X3, X4) %>% 
    select(-chr) %>% 
    select(X2, X3, X4) %>% 
    write_tsv('/wdata/lcasten/sli_wgs/prs/gene_sets/phastConsElements_vertebrates.160-400_mya.bed', col_names = FALSE)

mammal_cons %>% 
    anti_join(select(vertebrate_cons, X2, X3, X4)) %>% 
    mutate(chr = str_remove_all(X2, pattern = 'chr'),
            chr = as.numeric(chr)) %>% 
    drop_na() %>%
    arrange(chr, X3, X4) %>% 
    select(-chr) %>% 
    select(X2, X3, X4) %>% 
    write_tsv('/wdata/lcasten/sli_wgs/prs/gene_sets/phastConsElements_placental_mammals.100_mya.bed', col_names = FALSE)

primate_cons %>%
    anti_join(select(vertebrate_cons, X2, X3, X4)) %>% 
    anti_join(select(mammal_cons, X2, X3, X4)) %>% 
    mutate(chr = str_remove_all(X2, pattern = 'chr'),
            chr = as.numeric(chr)) %>% 
    drop_na() %>%
    arrange(chr, X3, X4) %>% 
    select(-chr) %>% 
    select(X2, X3, X4) %>% 
    write_tsv('/wdata/lcasten/sli_wgs/prs/gene_sets/phastConsElements_primates.65_mya.bed', col_names = FALSE)

##

## =============================================================
## make complementary BED file for evo sets 
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
anno_files <- anno_files[str_detect(anno_files, pattern = 'HAQERs_under_1700')]
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