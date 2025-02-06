library(tidyverse)


##
fac <- read_csv('/wdata/common/SLI_WGS/public/phenotype/factors/pheno_factors_resid.csv') %>% 
    rename(IID = sample)
fac_long <- fac %>% 
    pivot_longer(cols = matches('Factor'))

##
item_scores <- read_csv('/wdata/common/SLI_WGS/public/phenotype/pheno_imputed.csv')
item_scores_long <- item_scores %>% 
    rename(IID = core_id) %>% 
    inner_join(fac) %>%
    pivot_longer(cols = matches('^g|^Factor'))

##
wd_raw <- read_csv('/wdata/lcasten/sli_wgs/prs/pathway_prs/HAQER_PGS_manual_calculation_raw.csv')
wd_raw2 <- read_csv('/wdata/lcasten/sli_wgs/prs/pathway_prs/HAR_PGS_manual_calculation_raw.csv') %>% 
    rename(background_pgs_HAR = background_pgs)

# wd_raw = wd_raw %>% 
#     inner_join(wd_raw2)

## get QTL files info (i.e., HAQER -> scQTL cell-type map)
files = list.files('/wdata/lcasten/tools/ref_data/single_cell_QTL-Jerber-NatureGenetics-2021/eqtl_summary_stats_filtered', full.names = TRUE)
files = files[str_detect(files, '_treated', negate = TRUE)]
files

bed_list <- list()
for(f in files) {
    bed_list[[basename(f)]] <- read_tsv(f, show_col_types = FALSE) %>% 
        select(snp_chromosome, snp_position) %>%
        mutate(end = snp_position,
               snp_position = snp_position - 1) %>% 
        rename(chromosome = snp_chromosome, start = snp_position) %>%
        mutate(bed = basename(f))
}

annos <- bind_rows(bed_list) %>% 
    mutate(bed = str_split(bed, pattern = '[.]qtl', simplify = TRUE)[,1],
           bed = str_remove_all(bed, '.untreated'))
annos2 <- annos %>%
    mutate(bed = str_split(bed, pattern = '[.]', simplify = TRUE)[,1]) %>% 
    arrange(bed, chromosome, start) # %>% 
    # distinct(chromosome, start, .keep_all = TRUE)
annos3 <- annos %>% 
    mutate(bed = 'all')
annos <- bind_rows(annos, annos2,annos3)
unique(annos$bed)
table(annos$bed)


## intersect HAQERs with neurogenesis data
haqer_annot <- read_csv('/wdata/lcasten/sli_wgs/prs/pathway_prs/HAQER_PGS_manual_calculation_annotations.csv') %>% 
    rename(start = haqer_start, end = haqer_end) %>% 
    mutate(type = 'HAQER') %>% 
    select(chromosome, start, end, haqer) %>% 
    distinct()

## intersection
haqer_intersection <- haqer_annot %>% 
    fuzzyjoin::genome_inner_join(annos)

table(haqer_intersection$bed)

all_haqer_with_scqtl <- haqer_intersection %>% 
    distinct(haqer, bed) %>% 
    filter(bed == 'all')

## annotations to run through GREAT enrichment tool
haqer_annot %>%
    select(1:3) %>%
    mutate(chromosome = str_c('chr', chromosome)) %>%
    write_tsv('/wdata/lcasten/sli_wgs/prs/pathway_prs/HAQERs_with_SNPs.background.bed', col_names = FALSE)
haqer_annot %>% 
    filter(haqer %in% all_haqer_with_scqtl$haqer) %>% 
    select(1:3) %>%
    mutate(chromosome = str_c('chr', chromosome)) %>%
    write_tsv('/wdata/lcasten/sli_wgs/prs/pathway_prs/HAQERs_with_SNPs_and_ND-scQTL.bed', col_names = FALSE)
read_tsv('/wdata/lcasten/sli_wgs/prs/gene_sets/HAQERs.bed', col_names = FALSE) %>% 
    write_tsv('/wdata/lcasten/sli_wgs/prs/pathway_prs/HAQERs.background.bed', col_names = FALSE)

##
haqer_intersection2 <- haqer_intersection %>% 
    rename(haqer_pgs = bed) %>%
    select(name = haqer, haqer_pgs)

## compute PGS for those annot
ng_pgs <- wd_raw %>% 
    select(1:6, any_of(haqer_intersection$haqer)) %>% 
    pivot_longer(cols = matches('^HAQER')) %>% 
    inner_join(haqer_intersection2) %>% 
    distinct(IID, haqer_pgs, name, .keep_all = TRUE) %>% 
    # filter(IID == 'HG00096' & haqer_pgs == 'all') %>% distinct(name) ## check that there are only 207 HAQERs being used per individual (and that each HAQER is only being used once)
    group_by(IID, pc1, pc2, pc3, pc4, pc5, haqer_pgs) %>% 
    summarise(haqer_pgs_sum = sum(value))

#################3
## stats
ng_pgs %>% 
    filter(haqer_pgs == 'all') %>% 
    ggplot(aes(x = haqer_pgs_sum)) +
    geom_histogram()

ng_res <- ng_pgs %>% 
    inner_join(fac_long) %>% 
    inner_join(select(wd_raw, IID, background_pgs)) %>%
    group_by(name, haqer_pgs) %>% 
    mutate(pgs_resid = scale(resid(lm(haqer_pgs_sum ~ pc1 + pc2 + pc3 + pc4 + pc5)))[,1]) %>%
    # do(res = broom::tidy(lm(value ~ haqer_pgs_sum + pc1 + pc2 + pc3 + pc4 + pc5, data = .))) %>% 
    do(res = broom::tidy(cor.test(.$value, .$pgs_resid))) %>% 
    unnest(res) %>% 
    # filter(str_detect(term, 'pgs')) %>% 
    arrange(p.value) %>% 
    filter(name %in% str_c('Factor', 1:3))

p_ht <- ng_res %>% 
    filter(haqer_pgs %in% c('D11', 'D30', 'D52', 'all')) %>% 
    mutate(sig = ifelse(p.value < 0.05, '*', '')) %>% 
    mutate(timepoint = case_when(haqer_pgs == 'all' ~ 'HAQERs with ND-scQTL at any time point',
                                 haqer_pgs == 'D11' ~ 'HAQERs with ND-scQTL at day 11',
                                 haqer_pgs == 'D30' ~ 'HAQERs with ND-scQTL at day 30',
                                 haqer_pgs == 'D52' ~ 'HAQERs with ND-scQTL at day 52'),
           timepoint = factor(timepoint, levels = rev(c('HAQERs with ND-scQTL at day 11', 'HAQERs with ND-scQTL at day 30', 'HAQERs with ND-scQTL at day 52', 'HAQERs with ND-scQTL at any time point')))) %>%
    mutate(name = str_remove_all(name, pattern = 'actor')) %>%
    ggplot(aes(x = name, y = timepoint, fill = estimate)) +
    geom_tile() +
    scale_fill_gradient2(low = 'dodgerblue', mid = 'white', high = 'chocolate1', midpoint = 0) +
    geom_text(aes(label = sig), size = 10, check_overlap = TRUE, hjust = 0.5, vjust = 0.5) +
    labs(fill = "Correlation:") +
    xlab('Factor') +
    ylab('HAQER CP-PGS') +
    theme(legend.key.width = unit(1, 'cm'))

unique(ng_res$haqer_pgs)
# ng_res %>% 
#     filter(term != 'background_pgs') %>% 
#     filter(p.value < 0.05)
ng_res %>% 
    filter(p.value < 0.05)

## scatterplots of PGS using all HAQERs overlapping with neurodevelopmental scQTLs 
p_sc <- ng_pgs %>% 
    inner_join(fac_long) %>% 
    inner_join(select(wd_raw, IID, background_pgs)) %>%
    group_by(name, haqer_pgs) %>% 
    mutate(pgs_resid = scale(resid(lm(haqer_pgs_sum ~ pc1 + pc2 + pc3 + pc4 + pc5)))[,1]) %>% 
    filter(haqer_pgs == 'all') %>% 
    filter(name %in% str_c('Factor', 1:3)) %>% 
    inner_join(ng_res) %>% 
    mutate(lab = str_c('r = ', round(estimate, digits = 2), ', p = ', formatC(p.value, digits = 2))) %>% 
    mutate(sig = ifelse(p.value < 0.05, TRUE, FALSE)) %>%
    ggplot(aes(x = pgs_resid, y = value)) +
    geom_point(aes(color = sig)) +
    geom_smooth(method = 'lm', size = 1.5) +
    facet_wrap(~ name) +
    geom_text(aes(x = -1.5, y = 2.5, label = lab), check_overlap = TRUE, size = 5) +
    xlab('CP-PGS in HAQERs with ND-scQTL') +
    ylab('Factor score') +
    scale_color_manual(values = c('grey80', 'black')) +
    theme(legend.position = 'none')

