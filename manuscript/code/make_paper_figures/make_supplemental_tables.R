library(tidyverse)
library(rio)

## --------------------------------------------
## gather data
## --------------------------------------------
###################
## EpiSLI
###################
## pheno
ph <- read_csv('manuscript/supplemental_materials/EpiSLI_pheno_data.csv') %>% 
    select(1, matches('^F'))
## es-pgs
es_pgs <- read_csv('manuscript/supplemental_materials/EpiSLI_ES-PGS_data.csv')
## CBCL x Factor corr res
cbcl_res <- read_csv('manuscript/supplemental_materials/stats/EpiSLI_factor_CBCL_correlations.csv')
## LDpred2 PGS x Factor corr res
ldpred2_res <- read_csv('manuscript/supplemental_materials/stats/EpiSLI_factor_PGS_correlations.csv') %>% 
    mutate(n = parameter + 2) %>% 
    mutate(pgs_name = clean_name) %>% 
    relocate(factor, pgs_name, type) %>% 
    select(-c(name, clean_name)) %>% 
    rename(correlation_coefficient = estimate)
## ES-PGS x Factor corr res
es_pgs_res <- read_csv('manuscript/supplemental_materials/stats/EpiSLI_factor_ES-PGS_results.csv') %>% 
    filter(factor %in% str_c('F', 1:3))

#####################
## SPARK
#####################
## es-pgs results
spark_es_pgs <- read_csv('manuscript/supplemental_materials/stats/SPARK_ES-PGS_HAQER_validations.csv')
## reversion results
spark_rev <- read_csv('manuscript/supplemental_materials/stats/SPARK_rare_reversion_results.csv')

#####################
## ABCD
#####################
## es-pgs results
abcd_es_pgs <- read_csv('manuscript/supplemental_materials/stats/ABCD_HAQER_ES-PGS_results.csv')
## c-section family analysis
abcd_es_pgs_csec <- read_csv('manuscript/supplemental_materials/stats/ABCD_c-section_HAQER_ES-PGS_results.csv')

#####################
## AADR
#####################
## data
aadr <- read_csv('manuscript/supplemental_materials/AADR_data.csv')
aadr_imp <- read_csv('manuscript/supplemental_materials/AADR_imputed_data_no_neanderthals.csv')
## polygenic selection
aadr_sel <- read_csv('manuscript/supplemental_materials/stats/AADR_HAQER_ES-PGS_polygenic_selection_results.csv')
aadr_sel_imp <- read_csv('manuscript/supplemental_materials/stats/AADR_imputed_no_neanderthals_HAQER_ES-PGS_polygenic_selection_results.csv')
aadr_sel_gathered <- bind_rows(aadr_sel, aadr_sel_imp)
## polygenic correlation
aadr_pgs_corr <- read_csv('manuscript/supplemental_materials/stats/AADR_ES-PGS_polygenic_correlation_results.csv')

#####################
## Cross species
#####################
## data
cs_dat <- read_csv('manuscript/supplemental_materials/cross_species_vocal_learning_HAQER_similarity.csv')
## phylolm regression results
# read_csv('manuscript/supplemental_materials/stats/cross_species_HAQER_results.csv')

#####################
## HAQERs
#####################
## scQTL enrichment
scqtl_enr <- read_csv('manuscript/supplemental_materials/stats/HAQER_scQTL_enrichment_stats.csv')
## EpiSLI TFBS selection
tfbs_sel <- read_csv('manuscript/supplemental_materials/stats/TFBS_reversion_core_language_selection_results.csv')
## EpiSLI TFBS selection + language ability convergence
tfbs_conv <- read_csv('manuscript/supplemental_materials/stats/TFBS_TF_family_binding_convergence_results.csv')

####################################################
## gather data and write to huge excel spreadsheet
####################################################
out <- "manuscript/supplemental_materials/supplemental_tables.xlsx" ## path to gathered table
export(list(SupplementaryTable1 = ph, 
            SupplementaryTable2 = cbcl_res, 
            SupplementaryTable3 = ldpred2_res, 
            SupplementaryTable4 = es_pgs, 
            SupplementaryTable5 = es_pgs_res, 
            SupplementaryTable6 = spark_es_pgs, 
            SupplementaryTable7 = spark_rev, 
            SupplementaryTable8 = abcd_es_pgs, 
            SupplementaryTable9 = abcd_es_pgs_csec, 
            SupplementaryTable10 = aadr, 
            SupplementaryTable11 = aadr_imp, 
            SupplementaryTable12 = aadr_sel_gathered, 
            SupplementaryTable13 = aadr_pgs_corr, 
            SupplementaryTable14 = cs_dat, 
            SupplementaryTable15 = scqtl_enr, 
            SupplementaryTable16 = tfbs_sel, 
            SupplementaryTable17 = tfbs_conv),
        out)