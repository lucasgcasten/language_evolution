##
library(tidyverse)

#################################
## run enrichment analysis
#################################
## hg38 chromosome info
chrom_limits = read_table('manuscript/supplemental_materials/hg38.chrom.sizes', col_names = FALSE)
names(chrom_limits) = c('#chrom', 'chrom_size')

## annotation coords
outbed <- 'manuscript/supplemental_materials/HAQER_conservative_set.hg38.bed'
outbed_har <- 'manuscript/supplemental_materials/HAR.hg38.sorted_autosomes.bed'
outbed_rand <- 'manuscript/supplemental_materials/RAND.hg38.bed'
outbed_uce <- 'manuscript/supplemental_materials/UCE.hg38.bed'

## human brain cCRE coords from Supp table 6 from Li 2023 Science paper
cres <- read_table("/wdata/lcasten/tools/ref_data/CATlas_v1/Li_2023_HS/science.adf7044_tables_s1_to_s29/Supplementary Tables/Table S6 - List of cCREs in bed format/Table S6 – List of cCREs in bed format", col_names = FALSE)

## human brain cCRE classification from Supp table 20 from Li 2023 Science paper
cre_class <- read_table("/wdata/lcasten/tools/ref_data/CATlas_v1/Li_2023_HS/science.adf7044_tables_s1_to_s29/Supplementary Tables/Table S20 - List of different categories of cCREs.txt") 
cre_hs <- cre_class[cre_class$category == 'humanSpec',]
cre_cons <- cre_class[cre_class$category == 'CA_Cons',]
####################################################################
## compute enrichment of HAQERs around human-specific cCREs
####################################################################
for(ct in unique(cre_hs$subclass)) {
    ph <- ct
    set_name = str_c('hs_ca.', ph)
    outf <- str_c('manuscript/supplemental_materials/stats/chromatin_accessibility_evo/', ph, '.hg38.bed')
    cat('\n\n\n\n')
    message('============================================')
    message('Set: ', set_name)
    ## make BED file coords
    tmp <- cres %>% 
        filter(X4 %in% cre_hs$cCRE[cre_hs$subclass == ct])
    tmp %>% 
        rename(`#chrom` = X1, chromStart = X2, chromEnd = X3) %>% 
        select(`#chrom`, chromStart, chromEnd) %>%
        inner_join(chrom_limits) %>%
        mutate(chr = as.numeric(str_remove_all(`#chrom`, 'chr'))) %>% 
        arrange(chr, chromStart, chromEnd) %>%
        select(`#chrom`, chromStart, chromEnd) %>% 
        mutate(chromStart = pmax(0,chromStart - 100000), ## add flank around cCREs
               chromEnd = chromEnd + 100000) %>%
        write_tsv(outf, col_names = FALSE)

    ## merge overlapping regions for enrichment analysis
    outf2 = str_c('manuscript/supplemental_materials/stats/chromatin_accessibility_evo/', set_name, '.non_overlapping.hg38.bed')
    cmd <- str_c('bedtools merge -i ', outf, ' > ', outf2)
    system(cmd)

    ## run enrichment: overlapEnrichments method elements1.lift elements2.lift noGap.lift out.txt
    outres = str_c("manuscript/supplemental_materials/stats/chromatin_accessibility_evo/", set_name, '.HAQER_enrichment.txt')
    outres_har <- str_c("manuscript/supplemental_materials/stats/chromatin_accessibility_evo/", set_name, '.HAR_enrichment.txt')
    outres_rand <- str_c("manuscript/supplemental_materials/stats/chromatin_accessibility_evo/", set_name, '.RAND_enrichment.txt')
    outres_uce <- str_c("manuscript/supplemental_materials/stats/chromatin_accessibility_evo/", set_name, '.UCE_enrichment.txt')
    ## HAQER enrichment
    cmd <- str_c('~/go/bin/overlapEnrichments normalApproximate ', outf2, ' ', outbed, ' manuscript/supplemental_materials/hg38.chrom.sizes.bed ', outres)
    system(cmd)
    ## HAR enrichment
    cmd2 <- str_c('~/go/bin/overlapEnrichments normalApproximate ', outf2, ' ', outbed_har, ' manuscript/supplemental_materials/hg38.chrom.sizes.bed ', outres_har)
    system(cmd2)
    ## RAND enrichment
    cmd2 <- str_c('~/go/bin/overlapEnrichments normalApproximate ', outf2, ' ', outbed_rand, ' manuscript/supplemental_materials/hg38.chrom.sizes.bed ', outres_rand)
    system(cmd2)
    ## UCE enrichment
    cmd2 <- str_c('~/go/bin/overlapEnrichments normalApproximate ', outf2, ' ', outbed_uce, ' manuscript/supplemental_materials/hg38.chrom.sizes.bed ', outres_uce)
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
}

####################################################################
## compute enrichment of HAQERs around human-mouse CONSERVED cCREs
####################################################################
for(ct in unique(cre_cons$subclass)) {
    ph <- ct
    set_name = str_c('cons_ca.', ph)
    outf <- str_c('manuscript/supplemental_materials/stats/chromatin_accessibility_evo/', ph, '.hg38.bed')
    cat('\n\n\n\n')
    message('============================================')
    message('Set: ', set_name)
    ## make BED file coords
    tmp <- cres %>% 
        filter(X4 %in% cre_cons$cCRE[cre_cons$subclass == ct])
    tmp %>% 
        rename(`#chrom` = X1, chromStart = X2, chromEnd = X3) %>% 
        select(`#chrom`, chromStart, chromEnd) %>%
        inner_join(chrom_limits) %>%
        mutate(chr = as.numeric(str_remove_all(`#chrom`, 'chr'))) %>% 
        arrange(chr, chromStart, chromEnd) %>%
        select(`#chrom`, chromStart, chromEnd) %>% 
        mutate(chromStart = pmax(0,chromStart - 100000),
               chromEnd = chromEnd + 100000) %>%
        write_tsv(outf, col_names = FALSE)

    ## merge overlapping regions for enrichment analysis
    outf2 = str_c('manuscript/supplemental_materials/stats/chromatin_accessibility_evo/', set_name, '.non_overlapping.hg38.bed')
    cmd <- str_c('bedtools merge -i ', outf, ' > ', outf2)
    system(cmd)

    ## run enrichment: overlapEnrichments method elements1.lift elements2.lift noGap.lift out.txt
    outres = str_c("manuscript/supplemental_materials/stats/chromatin_accessibility_evo/", set_name, '.HAQER_enrichment.txt')
    outres_har <- str_c("manuscript/supplemental_materials/stats/chromatin_accessibility_evo/", set_name, '.HAR_enrichment.txt')
    outres_rand <- str_c("manuscript/supplemental_materials/stats/chromatin_accessibility_evo/", set_name, '.RAND_enrichment.txt')
    outres_uce <- str_c("manuscript/supplemental_materials/stats/chromatin_accessibility_evo/", set_name, '.UCE_enrichment.txt')
    ## HAQER enrichment
    cmd <- str_c('~/go/bin/overlapEnrichments normalApproximate ', outf2, ' ', outbed, ' manuscript/supplemental_materials/hg38.chrom.sizes.bed ', outres)
    system(cmd)
    ## HAR enrichment
    cmd2 <- str_c('~/go/bin/overlapEnrichments normalApproximate ', outf2, ' ', outbed_har, ' manuscript/supplemental_materials/hg38.chrom.sizes.bed ', outres_har)
    system(cmd2)
    ## RAND enrichment
    cmd2 <- str_c('~/go/bin/overlapEnrichments normalApproximate ', outf2, ' ', outbed_rand, ' manuscript/supplemental_materials/hg38.chrom.sizes.bed ', outres_rand)
    system(cmd2)
    ## UCE enrichment
    cmd2 <- str_c('~/go/bin/overlapEnrichments normalApproximate ', outf2, ' ', outbed_uce, ' manuscript/supplemental_materials/hg38.chrom.sizes.bed ', outres_uce)
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
}

################################
## gather and reformat results
################################
## human-specific
files <- list.files('manuscript/supplemental_materials/stats/chromatin_accessibility_evo/', pattern = 'hs_ca', full.names = TRUE)

res_list_hs = list()
for(f in files) {
    res_list_hs[[basename(f)]] <- read_table(f, show_col_types = FALSE) %>% 
        mutate(chromatin_accessibility_set = basename(Filename1),
            evo_annot = basename(Filename2)) %>% 
        relocate(evo_annot, chromatin_accessibility_set) %>% 
        select(-matches('Filename')) %>% 
        mutate(evo_annot = str_split(evo_annot, pattern = '[.]', simplify = TRUE)[,1],
               chromatin_accessibility_set = str_remove_all(chromatin_accessibility_set, pattern = 'hs_ca.|_enrichment.txt'),
               chromatin_accessibility_set = str_split(chromatin_accessibility_set, pattern = '[.]', simplify = TRUE)[,1],
               chromatin_accessibility_set = str_remove_all(chromatin_accessibility_set, pattern = '.non_overlapping')
               ) %>% 
        select(evo_annot, chromatin_accessibility_set, enrichment_method = `#Method`, n_elements_evo_annot = LenElements2, n_elements_chromatin_accessibility_set = LenElements1, n_overlapping_elements = OverlapCount, n_expected_overlaps = ExpectedOverlap, enrichment = Enrichment, enrichment_p = EnrichPValue)
}

hs_enr <- bind_rows(res_list_hs) %>% 
    mutate(type = 'human_specific_chromatin_accessibility') %>% 
    mutate(evo_annot = case_when(str_detect(evo_annot, 'HAQER') ~ 'HAQER',
                                 TRUE ~ evo_annot)) %>%
    relocate(type, .after = chromatin_accessibility_set)

## conserved chromatin accessibility
files <- list.files('manuscript/supplemental_materials/stats/chromatin_accessibility_evo/', pattern = 'cons_ca', full.names = TRUE)

res_list_cons = list()
for(f in files) {
    res_list_cons[[basename(f)]] <- read_table(f, show_col_types = FALSE) %>% 
        mutate(chromatin_accessibility_set = basename(Filename1),
            evo_annot = basename(Filename2)) %>% 
        relocate(evo_annot, chromatin_accessibility_set) %>% 
        select(-matches('Filename')) %>% 
        mutate(evo_annot = str_split(evo_annot, pattern = '[.]', simplify = TRUE)[,1],
               chromatin_accessibility_set = str_remove_all(chromatin_accessibility_set, pattern = 'hs_ca.|_enrichment.txt|cons_ca.'),
               chromatin_accessibility_set = str_split(chromatin_accessibility_set, pattern = '[.]', simplify = TRUE)[,1],
               chromatin_accessibility_set = str_remove_all(chromatin_accessibility_set, pattern = '.non_overlapping')
               ) %>% 
        select(evo_annot, chromatin_accessibility_set, enrichment_method = `#Method`, n_elements_evo_annot = LenElements2, n_elements_chromatin_accessibility_set = LenElements1, n_overlapping_elements = OverlapCount, n_expected_overlaps = ExpectedOverlap, enrichment = Enrichment, enrichment_p = EnrichPValue)
}

cons_enr <- bind_rows(res_list_cons) %>% 
    mutate(type = 'conserved_chromatin_accessibility') %>% 
    mutate(evo_annot = case_when(str_detect(evo_annot, 'HAQER') ~ 'HAQER',
                                 TRUE ~ evo_annot)) %>%
    relocate(type, .after = chromatin_accessibility_set)

################################
## merge results and save
################################
bind_rows(hs_enr, cons_enr) %>% 
    filter(evo_annot != 'UCE') %>% 
    mutate(chromatin_accessibility_set = ifelse(chromatin_accessibility_set == 'ITL23', 'ITL2-3', chromatin_accessibility_set),
           chromatin_accessibility_set = str_replace_all(chromatin_accessibility_set, 'ITL', 'ITL ')) %>%
    drop_na() %>%
    mutate(class = case_when(str_detect(chromatin_accessibility_set, 'FOXP2|LAMP|MSN|PVALB|SST|VIP') ~ 'GABA+',
                             str_detect(chromatin_accessibility_set, '^CT|^IT|L6|NP') ~ 'vGlut+',
                             str_detect(chromatin_accessibility_set, 'ASCT|MGC|OGC|OPC') ~ 'Non-neuronal'),
           class = factor(class, levels = c('GABA+', 'vGlut+', 'Non-neuronal'))) %>% 
    relocate(class, .after = type) %>%
    write_csv('manuscript/supplemental_materials/stats/HAQER_chromatin_accessibility_evo_Li2023_enrichment_stats.csv')

######################
## make figures
######################
cl <- c("#762776", "#e04468", "#dcc699")

hs_df <- hs_enr %>% 
    mutate(chromatin_accessibility_set = ifelse(chromatin_accessibility_set == 'ITL23', 'ITL2-3', chromatin_accessibility_set),
           chromatin_accessibility_set = str_replace_all(chromatin_accessibility_set, 'ITL', 'ITL ')) %>%
    mutate(evo_annot = factor(evo_annot, levels = c('RAND', 'HAR', 'HAQER'))) %>%
    arrange(desc(evo_annot), desc(enrichment_p)) %>%
    mutate(chromatin_accessibility_set = factor(chromatin_accessibility_set, levels = unique(chromatin_accessibility_set))) %>% 
    drop_na() %>%
    mutate(enrichment_p = ifelse(enrichment_p > 0.65, .65, enrichment_p)) %>% ## add a small constant value so they show up on the plot
    mutate(class = case_when(str_detect(chromatin_accessibility_set, 'FOXP2|LAMP|MSN|PVALB|SST|VIP') ~ 'GABA+',
                             str_detect(chromatin_accessibility_set, '^CT|^IT|L6|NP') ~ 'vGlut+',
                             str_detect(chromatin_accessibility_set, 'ASCT|MGC|OGC|OPC') ~ 'Non-neuronal'),
           class = factor(class, levels = c('GABA+', 'vGlut+', 'Non-neuronal')))

p_hs <- hs_df %>% 
    ggplot(aes(x = -log10(enrichment_p), y = chromatin_accessibility_set)) +
    geom_bar(stat = 'identity', aes(fill = evo_annot), position = 'dodge') +
    geom_vline(xintercept = -log10(0.05), color = 'red', linetype = 'dashed', size = 1.075) +
    ggforce::facet_col(facets = vars(class), 
                       scales = "free_y", 
                       space = "free") +
    labs(fill = NULL) +
    theme_classic() +
    theme(axis.text = element_text(size = 12),
          axis.title = element_text(size = 14),
          strip.text = element_text(size = 14),
          legend.text = element_text(size = 12),
          plot.title = element_text(size = 14, hjust = .5),
          legend.position = 'bottom') +
    scale_fill_manual(values = rev(cl)) +
    xlab('-log10(enrichment p-value)') +
    ylab('Cell-type') +
    ggtitle('Human-specific cCREs')

p_cons <- cons_enr %>% 
    mutate(chromatin_accessibility_set = ifelse(chromatin_accessibility_set == 'ITL23', 'ITL2-3', chromatin_accessibility_set),
          chromatin_accessibility_set = str_replace_all(chromatin_accessibility_set, 'ITL', 'ITL ')) %>%
    mutate(evo_annot = factor(evo_annot, levels = c('RAND', 'HAR', 'HAQER'))) %>%
    arrange(desc(evo_annot), desc(enrichment_p)) %>%
    mutate(chromatin_accessibility_set = factor(chromatin_accessibility_set, levels = levels(hs_df$chromatin_accessibility_set))) %>% 
    drop_na() %>%
    mutate(enrichment_p = ifelse(enrichment_p > 0.65, .65, enrichment_p)) %>% ## add a small constant value so they show up on the plot
    mutate(class = case_when(str_detect(chromatin_accessibility_set, 'FOXP2|LAMP|MSN|PVALB|SST|VIP') ~ 'GABA+',
                             str_detect(chromatin_accessibility_set, '^CT|^IT|L6|NP') ~ 'vGlut+',
                             str_detect(chromatin_accessibility_set, 'ASCT|MGC|OGC|OPC') ~ 'Non-neuronal'),
           class = factor(class, levels = c('GABA+', 'vGlut+', 'Non-neuronal'))) %>% 
    ggplot(aes(x = -log10(enrichment_p), y = chromatin_accessibility_set)) +
    geom_bar(stat = 'identity', aes(fill = evo_annot), position = 'dodge') +
    geom_vline(xintercept = -log10(0.05), color = 'red', linetype = 'dashed', size = 1.075) +
    ggforce::facet_col(facets = vars(class), 
                       scales = "free_y", 
                       space = "free") +
    labs(fill = NULL) +
    theme_classic() +
    theme(axis.text = element_text(size = 12),
          axis.title = element_text(size = 14),
          strip.text = element_text(size = 14),
          legend.text = element_text(size = 12),
          plot.title = element_text(size = 14, hjust = .5),
          legend.position = 'bottom') +
    scale_fill_manual(values = rev(cl)) +
    xlab('-log10(enrichment p-value)') +
    ylab('Cell-type') +
    ggtitle('Conserved cCREs')

##
library(patchwork)
p_merge <- p_hs + p_cons + plot_layout(guides = "collect") & theme(legend.position = 'bottom')
ggsave(p_merge,
       filename = 'manuscript/figures/HAQER_evolutionary_chromatin_accessibility_enrichment.png',
       dpi = 300, bg = 'white', device = 'png', 
       units = 'in', width = 11, height = 8)

## save plot object
p_merge %>% 
    write_rds('manuscript/figures/R_plot_objects/HAQER_evolutionary_chromatin_accessibility_enrichment.rds')