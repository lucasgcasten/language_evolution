library(tidyverse)

##
pc <- read_csv('/wdata/common/SLI_WGS/public/gene_lists/entrez_ensembl_symbol_map.csv') %>% 
    filter(source == 'Ensembl' & str_detect(biotype, pattern = 'protein_coding'))
##
files = list.files('/wdata/lcasten/sli_wgs/regenie/data/rare_var_results', pattern = 'regenie$', full.names = TRUE)
files = files[str_detect(files, 'sum_out')]
files = files[str_detect(files, 'genesets')]
files = files[str_detect(files, 'VEP_HIGH_or_CADD20', negate = TRUE)]

res_list = list()
i = 0
for(f in files){
    i = i + 1
    ph <- basename(f)
    ph = str_remove_all(ph, pattern = '.regenie|max_out_|sum_out_|.VEP_HIGH_or_CADD20_|.any_variant_|sum_out|max_out|genesets.')
    res_list[[i]] <- read_table(f, skip = 1) %>% 
        mutate(pheno = ph) %>% 
        relocate(pheno, ALLELE1)
}

res <- bind_rows(res_list) %>% 
    mutate(p = 10^(-LOG10P)) %>% 
    relocate(p)  %>% 
    arrange(p)

unique(res$ALLELE1)
unique(res$pheno)
unique(res$TEST)

res2 <- res %>% 
    relocate(ID, p, pheno, ALLELE1, A1FREQ, BETA) %>% 
    mutate(ID = str_split(ID, pattern = '[.]mask', simplify =TRUE)[,1]) # %>% 
    # filter(str_detect(ALLELE1, pattern = 'HIGH|MODERATE')) %>% 
    # filter(ID %in% pc$symbol)
res2 %>% 
    group_by(pheno) %>%
    mutate(padj = p.adjust(p, method = 'fdr')) %>% 
    filter(padj < 0.05)

bonferroni_cutoff = 0.05 / length(unique(res2$ID))
res2  %>% 
    filter(pheno == 'Factor1') %>% 
    distinct(ID, .keep_all = TRUE) %>%
    head(20)

res2  %>% 
    filter(pheno == 'Factor5') %>% 
    distinct(ID, .keep_all = TRUE) %>%
    head(20)

res %>% filter(TEST == 'ADD-ACATO')

## make plots
unique(res2$ALLELE1)
res2 %>% 
    filter(ALLELE1 == 'maskHIGHorMODERATE.0.01') %>%
    rename(POS = GENPOS, P = p) %>% 
    distinct() %>% 
    filter(TEST == 'ADD-ACATO') %>% 
    select(ID, P, pheno) %>% 
    filter(P < 0.05) %>% 
    arrange(pheno, P)

res3 <- res2 %>% 
    filter(ALLELE1 == 'maskHIGHorMODERATE.singleton') %>%
    rename(POS = GENPOS, P = p) %>% 
    distinct() %>% 
    filter(TEST == 'ADD')
res3 %>%
    select(1:3, BETA, TEST) #
res3 %>% filter(P < 0.05)
phenos <- unique(res3$pheno)
dir.create('/wdata/lcasten/sli_wgs/regenie/data/rare_var_results/figures')

unique(res3$ID)[str_detect(unique(res3$ID),'math')]
unique(res3$ID)
gene_set_map <- data.frame(ID = c('migration.covar_TDI', 'cognitive_performance_SSGAC_2018', 'Thickness', 'edu_years_cognitive_skills_Demange_2021', 'high_pLI', 'autism', 'Migration_EA_NEfactor', 'schizophrenia_PGC_2021', 'SurfArea', 'bipolar_disorder_I_and_II_PGC_2021', 'migration', 'low_LOEUF', 'ADHD2022_iPSYCH_deCODE_PGC', 'edu_years_non_cognitive_skills_Demange_2021', 'EA4', 'Migration_EA_Efactor', 'adult_cog_rare_edu', 'epilepsy', 'DDG2P', 'adult_cog_rare_vnr', 'epilepsy_nafe', 'epilepsy_gge', 'ASD_sfari', 'epilepsy_dee', 'adult_cog_rare_reaction', 'bip_rare', 'scz_rare', 'asd_rare',
                                  'neurotransmitters', 'synapse', 'brainCRE', 'functional_connectivity_Ventral.Attention', 'adhd_rare', 'functional_connectivity_Limbic', 'diff_methylated_chimp', 'functional_connectivity_Default', 'human_brain_expr_diverged_non_neur', 'HAR', 'depression', 'functional_connectivity_Frontoparietal', 'functional_connectivity_Somatomotor', 'language_connectivity', 'human_brain_expr_diverged_exc_neur', 'functional_connectivity_Dorsal.Attention', 'human_brain_expr_diverged_inh_neur', 'oligo_myelin', 'neuron', 'functional_connectivity_global', 'microglia', 'stimulus_response', 'diff_methylated_archaic', 'astrocyte', 'HAQER', 'functional_connectivity_Visual', 'left_right_brain_asymmetry', 'handedness_RvsL', 'diff_methylated_AMH',
                                  'glutamate', 'ieg', 'income_Hill2019', 'age_first_birth_Mills2021', 'action_potential', 'gaba', 'action_potential_har', 'rhythm', 'sleep_duration_Jansen2019', 'axon', 'neuroticism_GPC2', 'deafness', 'neuroestimator', 'extraversion_GPC2', 'language', 'insomnia_Jansen2019', 'serotonin', 'TDI_Neale2019', 'dyslexia',
                                  'KDM5B_CHEA_targets', 'KDM5A_CHEA_targets', 'Frontal_Cortex_BA9', 'Amygdala', 'FOXP2_CHEA_targets', 'Putamen_basal_ganglia', 'Anterior_cingulate_cortex_BA24', 'Cortex', 'Cerebellar_Hemisphere', 'FOXP1_CHEA_targets', 'Substantia_nigra', 'Hippocampus', 'language_merged', 'Caudate_basal_ganglia', 'Nucleus_accumbens_basal_ganglia', 'Cerebellum', 'language_gwas_catalog', 'Spinal_cord_cervical_c_1', 'antisocial_Tielbeek2022', 'executive_functioning_Hatoum2023', 'math', 'stuttering', 'memory', 'FOXP2_DAR_subset', 'FOXP2_DAR', 
                                  'jaw', 'larynx', 'tongue', 'ear', 'neural_crest_expr_human_less_than_chimp', 'neural_crest_expr_human_greater_than_chimp', 'neanderthal_selective_sweep_bottom5pct', 'pharynx', 'trachea', 'human_specific_missense_genes', 'vocal_cord', 'epiglottis', 'height_omim', 'forecasd_asd'),
           gs_clean = c('Migration (TDI)', 'Cognitive performance', 'Cortical thickness', 'Educational attainment cognitive skills', 'Constrained genes (pLI > 0.9)', 'Autism', 'Migration non-education factor', 'Schizophrenia', 'Cortical surface area', 'Bipolar', 'Migration', 'Low LOEUF (LOEUF < 0.35)', 'ADHD', 'Educational attainment non-cognitive skills', 'Educational attainment', 'Migration education factor', 'Educational attainment (UKBB rare var.)', 'Epilepsy (Rare var.)', 'Dev. disorders (DDG2P)', 'Verbal num. reasoning (UKBB rare var.)', 'Epilepsy NAFE (rare var.)', 'Epilepsy GGE (rare var.)', 'Autism (SFARI)', 'Epilepsy DEE (rare var.)', 'Reaction time (UKBB rare var.)', 'Bipolar (rare var.)', 'Schizophrenia (rare var.)', 'Autism (rare var.)',
                                  'Neurotransmitters', 'Synapse', 'Brain-CRE', 'FC ventral attention', 'ADHD (rare var.)', 'FC limbic', 'Human-chimp divergence', 'FC default mode', 'Div. brain expression: non-neuronal', 'Human accelerated regions', 'Depression', 'FC frontoparietal', 'FC somatomotor', 'Language network connectivity', 'Div. brain expression: exc-neuron', 'FC dorsal attention', 'Div. brain expression: inh-neuron', 'Oligodendrocytes/myelination', 'Neuron', 'FC global', 'Microglia', 'Stimulus response', 'Archaic hominin DMG', 'Astrocytes', 'HAQER', 'FC visual', 'Left-right brain asymmetry', 'Left handed', 'Modern human div. (DMG)',
                                  'Glutamate', 'Immediate early response genes', 'Income', 'Age of first birth', 'Action potential', 'GABA', 'Action potential speed/dendrites (HARs)', 'Rhythm', 'Sleep duration', 'Axon', 'Neuroticism', 'Deafness', 'Neuroestimator', 'Extraversion', 'Language impairment', 'Insomnia', 'Serotonin', 'TDI', 'Dyslexia',
                                  'KDM5B targets', 'KDM5A', 'Frontal cortex (BA9)', 'Amygdala', 'FOXP2 targets', 'Putamen-basal ganglia', 'ACC (BA24)', 'Cortex', 'Cerebellar hemisphere', 'FOXP1', 'Substantia nigra', 'Hippocampus', 'Language (merged)', 'Caudate-basal ganglia', 'Nucleus accumbens-basal ganglia', 'Cerebellum', 'Language (gwas catalog)', 'Spinal cord', 'Antisocial behavior', 'Executive functioning', 'Math ability', 'Stuttering', 'Memory', 'FOXP2 DARs (subset)', 'FOXP2 DARs',
                                  'Jaw', 'Larynx', 'Tongue', 'Ear', 'Neural crest cell DEGs: human < chimp', 'Neural crest cell DEGs: human > chimp', 'Neanderthal selective sweep', 'Pharynx', 'Trachea', 'Human specific missense genes', 'Vocal cord', 'Epiglottis', 'Height', 'Autism (forecASD)')) %>% 
           mutate(gs = case_when(str_detect(gs_clean, pattern = 'ADHD|Autism|Depress|Schizo|Bipolar|Epilepsy|DD') ~ 'Neuropsychiatric',
                                 str_detect(gs_clean, pattern = 'Migration') ~ 'Migration',
                                 str_detect(gs_clean, pattern = 'FC|connectiv|Cortical|Language network|asymm') ~ 'Brain MRI',
                                #  str_detect(gs_clean, pattern = 'Education|Cognitive|UKBB|Rhyth|Dyslex|Language|Executive|Math|Stuttering|Memory') ~ 'Cognitive',
                                 str_detect(gs_clean, pattern = 'Insomn|Sleep|Income|Age of first|Neuroticis|Extravers|TDI|Anti|Education|Cognitive|UKBB|Rhyth|Dyslex|Language|Executive|Math|Stuttering|Memory') ~ 'Cognitive/Behavioral/SES',
                                #  str_detect(gs_clean, pattern = 'pLI|LOEUF|chimp') ~ 'Evo/Constraint',
                                #  str_detect(gs_clean, pattern = 'KDM5|FOXP') ~ 'TF targets',
                                 str_detect(gs_clean, pattern = 'KDM5|FOXP|pLI|LOEUF|chimp|Height|Left |Deafn|hand') ~ 'Misc.',
                                 str_detect(gs_clean, pattern = 'Neurotransmitter|Synapse|Brain-CRE|Oligo|^Neuron$|Microglia|Astroc|Stimulus|Action|GABA|Glut|Seroton|Immediate|Axon|Neuroestimat') ~ 'Neural',
                                 str_detect(gs_clean, pattern = 'Frontal cortex|Amyg|Putam|ACC |^Cortex$|Cerebell|Substantia|Hippo|Caudate|Nucleus acc|Spinal') ~ 'Brain tissues',
                                 str_detect(gs_clean, pattern = 'Jaw|Laryn|Tong|Ear|Phary|Trache|Vocal|Epiglot') ~ 'Other tissues',
                                 str_detect(gs_clean, pattern = 'DMG|brain expression|Human accel|HAQER|Neander|Neural crest|Human specif') ~ 'Evolution'
                                 )) %>% 
            distinct()
nrow(gene_set_map)
table(gene_set_map$gs)


##################################################################
p <- res3 %>% 
    mutate(sig = case_when(P < 0.05 ~ TRUE, 
                           TRUE ~ FALSE)) %>%
    inner_join(gene_set_map) %>%
    arrange(desc(gs_clean)) %>%
    # filter(gs %in% c('Cognitive', 'Developmental/Psychiatric', 'Behavioral/SES')) %>%
    filter(str_detect(gs_clean, pattern = 'birth|rare|Rare|SFARI|Immediate|subset|cognitive skill|glottis|Jaw|Tongue|Pharynx|Trache|DAR|DMG|specific|Migration|Neural crest|estimator|brain expr|Math|Memory|duration|ction|Stimul|accel|HAQER|FOXP1|KDM5|Neander|Ear|Laryn|Vocal|Seroton|Astro|GABA|Glutamat|catalog|merged|LOEUF', negate = TRUE)) %>%
    mutate(gs_clean = factor(gs_clean, levels = unique(gs_clean))) %>%
    filter(str_detect(gs_clean, pattern = 'Epilepsy [A-Z]', negate = TRUE)) %>%
    mutate(pheno = str_remove_all(pheno, pattern = 'actor')) %>%
    filter(pheno %in% c('F1', 'F2', 'F3')) %>%
    mutate(pheno = case_when(pheno == 'F1' ~ 'Core language (F1)',
                             pheno == 'F2' ~ 'Receptive language (F2)',
                             pheno == 'F3' ~ 'NVIQ (F3)'),
           pheno = factor(pheno, levels = c('Core language (F1)','Receptive language (F2)', 'NVIQ (F3)'))) %>% 
    mutate(gs = factor(gs, levels = c('Cognitive/Behavioral/SES',  'Neuropsychiatric', 'Brain tissues', 'Brain MRI', 'Neural', 'Misc.'))) %>%
    mutate(xmin = BETA - 1.96 * SE, 
           xmax = BETA + 1.96 * SE) %>%
    # mutate(xmin = ifelse(xmin < -0.0825, -0.0825, xmin),
    #        xmax = ifelse(xmax > 0.02, 0.02, xmax)) %>%
    ggplot(aes(x = BETA, y = gs_clean, color = sig)) +
    geom_point(size = 3) +
    geom_linerange(aes(xmin = xmin, xmax = xmax), size = 1.05) +
    facet_grid(cols = vars(pheno), rows = vars(gs), scales = 'free_y', space = 'free_y') +
    geom_vline(xintercept = 0, color = 'red', linetype = 'dashed', size = 1.075) +
    scale_color_manual(values = c('grey85', 'black')) +
    xlab('Regression beta for damaging singleton variants (95% CI)') +
    ylab('Gene set') +
    # ggtitle('REGENIE rare variant association results: HIGH or MODERATE impact singleton variants') +
    theme_classic() +
    theme(legend.position = 'none',
          panel.border = element_rect(color="black", size=1, linetype="solid", fill = NA),
          strip.background = element_rect(colour = "black", fill="white", size = 1),
          strip.placement = 'outside',
          strip.text.x = element_text(size = 12),
          strip.text.y = element_text(size = 12),
          axis.text = element_text(size = 10),
          axis.title = element_text(14)) +
    coord_cartesian(xlim = c(-.08, .03))


ggsave(p, filename = '/wdata/lcasten/sli_wgs/paper_figures/fig3_rare_variant_burden_genesets.png', device = 'png', dpi = 300, units = 'in', width = 12, height = 16)

##
hits <- res3 %>% 
    filter(str_detect(pheno, '1|2|3')) %>%
    filter(P < 0.05) %>% 
    distinct(ID)
p_poster <- res3 %>% 
    mutate(sig = case_when(P < 0.05 ~ TRUE, 
                           TRUE ~ FALSE)) %>%
    inner_join(gene_set_map) %>%
    arrange(desc(gs_clean)) %>%
    # filter(gs %in% c('Cognitive', 'Developmental/Psychiatric', 'Behavioral/SES')) %>%
    filter(str_detect(gs_clean, pattern = 'birth|rare|Rare|SFARI|Immediate|subset|cognitive skill|glottis|Jaw|Tongue|Pharynx|Trache|DAR|DMG|specific|Migration|Neural crest|estimator|brain expr|Math|Memory|duration|ction|Stimul|accel|HAQER|FOXP1|KDM5|Neander|Ear|Laryn|Vocal|Seroton|Astro|GABA|Glutamat|catalog|merged|LOEUF', negate = TRUE)) %>%
    mutate(gs_clean = factor(gs_clean, levels = unique(gs_clean))) %>%
    filter(str_detect(gs_clean, pattern = 'Epilepsy [A-Z]', negate = TRUE)) %>%
    mutate(pheno = str_remove_all(pheno, pattern = 'actor')) %>%
    filter(pheno %in% c('F1', 'F2', 'F3')) %>%
    mutate(pheno = case_when(pheno == 'F1' ~ 'Core language (F1)',
                             pheno == 'F2' ~ 'Receptive language (F2)',
                             pheno == 'F3' ~ 'NVIQ (F3)'),
           pheno = factor(pheno, levels = c('Core language (F1)','Receptive language (F2)', 'NVIQ (F3)'))) %>% 
    mutate(gs = str_replace_all(gs, 'Cognitive/Behavioral/SES', 'Cog.'),
           gs = str_replace_all(gs, 'Brain MRI|Brain tissues|Cog.', 'Neural'),
           gs = str_replace_all(gs, 'Misc.', 'Evo')) %>%
    mutate(gs = factor(gs, levels = c('Cognitive/Behavioral/SES',  'Neuropsychiatric', 'Brain tissues', 'Brain MRI', 'Neural', 'Evo'))) %>%
    mutate(xmin = BETA - 1.96 * SE, 
           xmax = BETA + 1.96 * SE) %>%
    filter(ID %in% hits$ID | str_detect(gs_clean, 'Autism (fore)|ADHD|Schizo|Bipolar')) %>%
    # mutate(xmin = ifelse(xmin < -0.0825, -0.0825, xmin),
    #        xmax = ifelse(xmax > 0.02, 0.02, xmax)) %>%
    ggplot(aes(x = BETA, y = gs_clean, color = sig)) +
    geom_point(size = 3) +
    geom_linerange(aes(xmin = xmin, xmax = xmax), size = 1.05) +
    facet_grid(cols = vars(pheno), rows = vars(gs), scales = 'free_y', space = 'free_y') +
    geom_vline(xintercept = 0, color = 'red', linetype = 'dashed', size = 1.075) +
    scale_color_manual(values = c('grey85', 'black')) +
    xlab('Regression beta for damaging singleton variants (95% CI)') +
    ylab('Gene set') +
    # ggtitle('REGENIE rare variant association results: HIGH or MODERATE impact singleton variants') +
    theme_classic() +
    theme(legend.position = 'none',
          panel.border = element_rect(color="black", size=1, linetype="solid", fill = NA),
          strip.background = element_rect(colour = "black", fill="white", size = 1),
          strip.placement = 'outside',
          strip.text.x = element_text(size = 12),
          strip.text.y = element_text(size = 12),
          axis.text = element_text(size = 10),
          axis.title = element_text(14)) +
    coord_cartesian(xlim = c(-.08, .02)) +
    scale_x_continuous(breaks = c(-0.075, -0.05, -0.025, 0))
ggsave(p_poster, filename = '/wdata/lcasten/sli_wgs/paper_figures/rare_variant_burden_genesets_subset.png', device = 'png', dpi = 300, units = 'in', width = 12, height = 10)


##############################
## make plot for all factors
p2 <- res3 %>% 
    mutate(sig = case_when(P < 0.05 ~ TRUE, 
                           TRUE ~ FALSE)) %>%
    inner_join(gene_set_map) %>%
    arrange(desc(gs_clean)) %>%
    # filter(gs %in% c('Cognitive', 'Developmental/Psychiatric', 'Behavioral/SES')) %>%
    filter(str_detect(gs_clean, pattern = 'birth|rare|Rare|SFARI|Immediate|subset|cognitive skill|glottis|Jaw|Tongue|Pharynx|Trache|DAR|DMG|specific|Migration|Neural crest|estimator|brain expr|Math|Memory|duration|ction|Stimul|accel|HAQER|FOXP1|KDM5|Neander|Ear|Laryn|Vocal|Seroton|Astro|GABA|Glutamat|catalog|merged|LOEUF', negate = TRUE)) %>%
    mutate(gs_clean = factor(gs_clean, levels = unique(gs_clean))) %>%
    filter(str_detect(gs_clean, pattern = 'Epilepsy [A-Z]', negate = TRUE)) %>%
    mutate(pheno = str_remove_all(pheno, pattern = 'actor')) %>%
    mutate(gs = factor(gs, levels = c('Cognitive/Behavioral/SES',  'Neuropsychiatric', 'Brain tissues', 'Brain MRI', 'Neural', 'Misc.'))) %>%
    mutate(xmin = BETA - 1.96 * SE, 
           xmax = BETA + 1.96 * SE) %>%
    mutate(xmin = ifelse(xmin < -0.0825, -0.0825, xmin),
           xmax = ifelse(xmax > 0.02, 0.02, xmax)) %>%
    ggplot(aes(x = BETA, y = gs_clean, color = sig)) +
    geom_point(size = 3) +
    geom_linerange(aes(xmin = xmin, xmax = xmax), size = 1.05) +
    facet_grid(cols = vars(pheno), rows = vars(gs), scales = 'free_y', space = 'free_y') +
    geom_vline(xintercept = 0, color = 'red', linetype = 'dashed', size = 1.075) +
    scale_color_manual(values = c('grey85', 'black')) +
    xlab('Regression beta for damaging singleton variants (95% CI)') +
    ylab('Gene set') +
    # ggtitle('REGENIE rare variant association results: HIGH or MODERATE impact singleton variants') +
    theme_classic() +
    theme(legend.position = 'none',
          panel.border = element_rect(color="black", size=1, linetype="solid", fill = NA),
          strip.background = element_rect(colour = "black", fill="white", size = 1),
          strip.placement = 'outside',
          strip.text.x = element_text(size = 12),
          strip.text.y = element_text(size = 8)) +
    coord_cartesian(xlim = c(-.0775, .0075))

ggsave(p2, filename = '/wdata/lcasten/sli_wgs/paper_figures/supp_fig_rare_variant_burden_genesets_all_factors.png', device = 'png', dpi = 300, units = 'in', width = 16, height = 16)
