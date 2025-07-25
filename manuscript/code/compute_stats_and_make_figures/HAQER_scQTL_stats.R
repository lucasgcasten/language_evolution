##
library(tidyverse)

#################################
## run enrichment analysis
#################################
## hg19 chromosome info
chrom_limits = read_table('manuscript/supplemental_materials/hg19.chrom.sizes', col_names = FALSE)
names(chrom_limits) = c('#chrom', 'chrom_size')

## annotation coords
outbed <- 'manuscript/supplemental_materials/HAQER.hg19.sorted_autosomes_non_overlapping.bed'
outbed_har <- 'manuscript/supplemental_materials/HAR.hg19.sorted_autosomes_non_overlapping.bed'
outbed_rand <- 'manuscript/supplemental_materials/RAND.hg19.sorted_autosomes.bed'
outbed_uce <- 'manuscript/supplemental_materials/UCE.hg19.sorted_autosomes.bed'

## neurodevelopmental scQTL files (filtered to snpXgene paires with p < 5e-04)
files = list.files('manuscript/supplemental_materials/single_cell_QTL-Jerber-NatureGenetics-2021', pattern = '.txt', full.names = TRUE)

##
for(f in files) {
    ph <- basename(f)
    ph = str_remove_all(ph, pattern = '.txt')
    set_name = str_c('scQTL.', ph)
    outf <- str_c('manuscript/supplemental_materials/stats/HAQER_scQTL_enrichment/', ph, '.hg19.bed')
    cat('\n\n\n\n')
    message('============================================')
    message('Set: ', set_name)
    ## make BED file coords
    tmp <- read_tsv(f) 
    tmp <- tmp %>% 
        mutate(chromEnd = snp_position) %>%
        mutate(snp_position = snp_position - 1) %>%
        rename(`#chrom` = snp_chromosome, chromStart = snp_position) %>% 
        relocate(`#chrom`, chromStart, chromEnd)
    tmp %>%
        mutate(`#chrom` = str_c('chr', `#chrom`)) %>%
        inner_join(chrom_limits) %>%
        mutate(chr = as.numeric(str_remove_all(`#chrom`, 'chr'))) %>% 
        arrange(chr, chromStart, chromEnd) %>%
        select(`#chrom`, chromStart, chromEnd) %>% 
        write_tsv(outf, col_names = FALSE)

    ## merge overlapping regions for enrichment analysis
    outf2 = str_c('manuscript/supplemental_materials/stats/HAQER_scQTL_enrichment/', set_name, '.non_overlapping.hg19.bed')
    cmd <- str_c('bedtools merge -i ', outf, ' > ', outf2)
    system(cmd)

    ## run enrichment: overlapEnrichments method elements1.lift elements2.lift noGap.lift out.txt
    outres = str_c("manuscript/supplemental_materials/stats/HAQER_scQTL_enrichment/", set_name, '.HAQER_enrichment.txt')
    outres_har <- str_c("manuscript/supplemental_materials/stats/HAQER_scQTL_enrichment/", set_name, '.HAR_enrichment.txt')
    outres_rand <- str_c("manuscript/supplemental_materials/stats/HAQER_scQTL_enrichment/", set_name, '.RAND_enrichment.txt')
    outres_uce <- str_c("manuscript/supplemental_materials/stats/HAQER_scQTL_enrichment/", set_name, '.UCE_enrichment.txt')
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
}

############################################################
## post mortem adult brain scQTLs from psychENCODE 2 (hg38)
############################################################
chrom_limits <- read_table('manuscript/supplemental_materials/hg38.chrom.sizes', col_names = FALSE)
names(chrom_limits) = c('#chrom', 'chrom_size')
files = list.files('/wdata/lcasten/tools/ref_data/psychENCODE2/cell_type_eQTL', pattern = '.dat$', full.names = TRUE)

outbed_rand <- 'manuscript/supplemental_materials/RAND.hg38.bed'
outbed = 'manuscript/supplemental_materials/HAQER.hg38.sorted_autosomes.bed'
outbed_har <- 'manuscript/supplemental_materials/HAR.hg38.sorted_autosomes.bed'
outbed_uce <- 'manuscript/supplemental_materials/UCE.hg38.bed'


##
for(f in files) {
    ph <- basename(f)
    ph = str_remove_all(ph, pattern = '.dat')
    set_name = str_c('adult_scQTL.', ph)
    tmp <- read_table(f, col_names = FALSE)  %>% 
            filter(str_detect(X9, '[0-9]')) %>% 
            select(9:11) %>%
            rename(`#chrom` = X9, chromStart = X10, chromEnd = X11) %>% 
            mutate(chromStart = ifelse(chromStart == chromEnd, chromStart - 1, chromStart))
   
    cat('\n\n\n\n')
    message('============================================')
    message('Set: ', set_name)
    ##
    outf <- str_c('manuscript/supplemental_materials/stats/HAQER_scQTL_enrichment/', ph, '.hg38.bed')

    tmp %>%
        inner_join(chrom_limits) %>%
        mutate(chr = as.numeric(str_remove_all(`#chrom`, 'chr'))) %>% 
        arrange(chr, chromStart, chromEnd) %>%
        select(`#chrom`, chromStart, chromEnd) %>% 
        write_tsv(outf, col_names = FALSE)
    ## merge overlapping regions for enrichment analysis
    outf2 = str_c('manuscript/supplemental_materials/stats/HAQER_scQTL_enrichment/', set_name, '.non_overlapping.hg38.bed')
    cmd <- str_c('bedtools merge -i ', outf, ' > ', outf2)
    system(cmd)

    ## run enrichment: overlapEnrichments method elements1.lift elements2.lift noGap.lift out.txt
    outres = str_c("manuscript/supplemental_materials/stats/HAQER_scQTL_enrichment/", set_name, '.HAQER_enrichment.txt')
    outres_har <- str_c("manuscript/supplemental_materials/stats/HAQER_scQTL_enrichment/", set_name, '.HAR_enrichment.txt')
    outres_rand <- str_c("manuscript/supplemental_materials/stats/HAQER_scQTL_enrichment/", set_name, '.RAND_enrichment.txt')
    outres_uce <- str_c("manuscript/supplemental_materials/stats/HAQER_scQTL_enrichment/", set_name, '.UCE_enrichment.txt')
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
## neurodevelopmental scQTLs
files <- list.files('manuscript/supplemental_materials/stats/HAQER_scQTL_enrichment', pattern = 'D[0-9]', full.names = TRUE)

res_list_pre = list()
for(f in files) {
    res_list_pre[[basename(f)]] <- read_table(f) %>% 
        mutate(scQTL_set = basename(Filename1),
            evo_annot = basename(Filename2)) %>% 
        relocate(evo_annot, scQTL_set) %>% 
        select(-matches('Filename')) %>% 
        mutate(evo_annot = str_split(evo_annot, pattern = '[.]', simplify = TRUE)[,1],
               scQTL_set = str_remove_all(scQTL_set, pattern = '^scQTL[.]D'),
               day = str_split(scQTL_set, pattern = '[.]', simplify = TRUE)[,1],
               scQTL_set = str_split(scQTL_set, pattern = '[.]', simplify = TRUE)[,2],
               cell_type = case_when(scQTL_set == 'Epen1' ~ 'Ependymal',
                                     scQTL_set == 'DA' ~ 'Dopaminergic',
                                     scQTL_set == 'FPP' ~ 'Floor plate progenitors',
                                     scQTL_set == 'Sert' ~ 'Serotinergic',
                                     scQTL_set == 'Astro' ~ 'Astrocytes',
                                     scQTL_set == 'pseudobulk' ~ 'Any (pseudobulk)',
                                     scQTL_set == 'P_FPP' ~ 'Proliferating floor plate progenitors'),
               scQTL_set = str_c('day', day, '.', cell_type)) %>% 
        select(evo_annot, scQTL_set, enrichment_method = `#Method`, n_elements_evo_annot = LenElements2, n_elements_scQTL_set = LenElements1, n_overlapping_elements = OverlapCount, n_expected_overlaps = ExpectedOverlap, enrichment = Enrichment, enrichment_p = EnrichPValue)
}

pre_scqtl <- bind_rows(res_list_pre) %>% 
    mutate(scQTL_type = 'developing_iPSC_neurons.Jerber-NatGen2021') %>% 
    relocate(scQTL_type, .after = scQTL_set)

## adult scQTLs
files <- list.files('manuscript/supplemental_materials/stats/HAQER_scQTL_enrichment', pattern = 'adult', full.names = TRUE)

res_list_adult = list()
for(f in files) {
    res_list_adult[[basename(f)]] <- read_table(f) %>% 
        mutate(scQTL_set = basename(Filename1),
            evo_annot = basename(Filename2)) %>% 
        relocate(evo_annot, scQTL_set) %>% 
        select(-matches('Filename')) %>% 
        mutate(evo_annot = str_split(evo_annot, pattern = '[.]', simplify = TRUE)[,1],
            scQTL_set = str_split(scQTL_set, pattern = '[.]', simplify = TRUE)[,2],
            scQTL_set = str_remove_all(scQTL_set, pattern = '_sig_QTLs')) %>% 
        select(evo_annot, scQTL_set, enrichment_method = `#Method`, n_elements_evo_annot = LenElements2, n_elements_scQTL_set = LenElements1, n_overlapping_elements = OverlapCount, n_expected_overlaps = ExpectedOverlap, enrichment = Enrichment, enrichment_p = EnrichPValue)
}

adult_scqtl <- bind_rows(res_list_adult) %>% 
    mutate(scQTL_type = 'adult_brain_post_mortem.Emani-Science2024') %>% 
    relocate(scQTL_type, .after = scQTL_set)

#########
## merge results
bind_rows(pre_scqtl, adult_scqtl) %>% 
    filter(evo_annot != 'UCE') %>% 
    write_csv('manuscript/supplemental_materials/stats/HAQER_scQTL_enrichment_stats.csv')


######################
## make figures
######################
cl <- c("#762776", "#e04468", "#dcc699")

p_prenatal <- pre_scqtl %>% 
    mutate(day = str_split(scQTL_set, pattern = '[.]', simplify = TRUE)[,1],
           cell_type = str_split(scQTL_set, pattern = '[.]', simplify = TRUE)[,2]) %>% 
    mutate(evo_annot = factor(evo_annot, levels = c('RAND', 'HAR', 'HAQER'))) %>%
    mutate(day = str_replace_all(day, pattern = 'day', replacement = 'Postdifferention day ')) %>%
    arrange(evo_annot, enrichment_p) %>%
    mutate(scQTL_set = factor(scQTL_set, levels = unique(scQTL_set))) %>% 
    # mutate(cell_type = factor(cell_type, levels = unique(cell_type))) %>%
    drop_na() %>%
    mutate(enrichment_p = ifelse(enrichment_p > 0.65, .65, enrichment_p)) %>% ## add a small constant value so they show up on the plot
    ggplot(aes(x = -log10(enrichment_p), y = cell_type)) +
    geom_bar(stat = 'identity', aes(fill = evo_annot), position = 'dodge') +
    geom_vline(xintercept = -log10(0.05), color = 'red', linetype = 'dashed', size = 1.075) +
    ggforce::facet_col(facets = vars(day), 
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
    ggtitle('Prenatal scQTLs')

p_post <- adult_scqtl %>% 
    mutate(evo_annot = factor(evo_annot, levels = c('RAND', 'HAR', 'HAQER'))) %>%
    arrange(evo_annot, desc(enrichment_p)) %>%
    mutate(scQTL_set = str_replace_all(scQTL_set, pattern = '__', replacement = ' ')) %>%
    mutate(scQTL_set = factor(scQTL_set, levels = unique(scQTL_set))) %>% 
    drop_na() %>%
    mutate(enrichment_p = ifelse(enrichment_p > 0.65, .65, enrichment_p)) %>% ## add a small constant value so they show up on the plot
    ggplot(aes(x = -log10(enrichment_p), y = scQTL_set)) +
    geom_bar(stat = 'identity', aes(fill = evo_annot), position = 'dodge') +
    geom_vline(xintercept = -log10(0.05), color = 'red', linetype = 'dashed', size = 1.075) +
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
    ggtitle('Postnatal scQTLs')

##
library(patchwork)
p_merge <- p_prenatal + p_post + plot_layout(guides = "collect") & theme(legend.position = 'bottom')
ggsave(p_merge,
       filename = 'manuscript/figures/HAQER_scQTL_enrichment.png',
       dpi = 300, bg = 'white', device = 'png', 
       units = 'in', width = 11, height = 8)

## save plot object
p_merge %>% 
    write_rds('manuscript/figures/R_plot_objects/HAQER_scQTL_enrichment.rds')