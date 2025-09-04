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
es_pgs_res <- read_csv('manuscript/supplemental_materials/stats/EpiSLI_factor_ES-PGS_results.csv')

#####################
## SPARK
#####################
## es-pgs results
spark_es_pgs <- read_csv('manuscript/supplemental_materials/stats/SPARK_ES-PGS_HAQER_validations.csv')
spark_self_dx_es_pgs <- read_csv("manuscript/supplemental_materials/stats/SPARK_ES-PGS_HAQER_validations_self_reported_language_diagnosis.csv")[,1:8]
## reversion results
spark_rev <- read_csv('manuscript/supplemental_materials/stats/SPARK_rare_reversion_results.csv')

#####################
## ABCD
#####################
## es-pgs results
abcd_es_pgs <- read_csv('manuscript/supplemental_materials/stats/ABCD_HAQER_ES-PGS_results.csv')

#####################
## AADR + archaic
#####################
## polygenic selection
aadr_sel <- read_csv('manuscript/supplemental_materials/stats/AADR_HAQER_ES-PGS_polygenic_selection_results.csv')
## archaic human ES-PGS comparisons with 1000 Genomes
nean_pgs <- read_csv("manuscript/supplemental_materials/stats/AADR_HAQER_ES-PGS_polygenic_group_comparison_results.csv")

#####################
## Cross species
#####################
## data
cs_dat <- read_csv('manuscript/supplemental_materials/cross_species_vocal_learning_HAQER_similarity.csv') %>% 
    select(1:4, brain_mass_g, neonate_body_mass_g_PanTHERIA, adult_body_mass_g_PanTHERIA, matches("HAQER_sequence_sim"))
## phylolm regression results
cs_res <- read_csv('manuscript/supplemental_materials/stats/cross_species_HAQER_results.csv')

#####################
## HAQERs
#####################
## scQTL enrichment
scqtl_enr <- read_csv('manuscript/supplemental_materials/stats/HAQER_scQTL_enrichment_stats.csv')
## birth head circ GWAS loci enrichment
bhc_enr <- read_csv("manuscript/supplemental_materials/stats/HAQER_birth_head_circ_gwas_Vogelezang2022_enrichment_stats.csv")
## mammalian vocal learning regions
vl_enr <- read_csv("manuscript/supplemental_materials/stats/HAQER_vocal_learning_Wirthlin2024_enrichment_stats.csv")

## EpiSLI TFBS selection
tfbs_sel <- read_csv('manuscript/supplemental_materials/stats/TFBS_reversion_core_language_selection_results.csv')
## EpiSLI TFBS selection + language ability convergence
tfbs_conv <- read_csv('manuscript/supplemental_materials/stats/TFBS_TF_family_binding_convergence_results.csv')

#####################
## EpiSLI WGS stats
#####################
## raw table vals to make into df
# 	        Mean	0.00%	25.00%	50.00%	75.00%	100.00%
# Coverage	31.95	22	29	31	35	43
# Insert-size	384.5	260.5	361.8	395.2	411.3	526.7
episli_wgs_cov <- data.frame("Mean" = c(31.95, 384.5), 
            "p0" = c(22, 260.5), 
            "p25" = c(29, 361.8), 
            "p50" = c(31, 395.2), 
            "p75" = c(35, 411.3), 
            "p100" = c(43, 526.7))
rownames(episli_wgs_cov) <- c("Coverage", "Insert-size")

####################################################
## gather data and write to huge excel spreadsheet
####################################################
out <- "manuscript/supplemental_materials/supplemental_tables.xlsx" ## path to gathered table
tb <- data.frame(table = str_c("SupplementaryTable", 1:17),
                 description = c("Correlation results from EpiSLI cohort factor scores with CBCL scales",
                                 'Correlation results from EpiSLI cohort factor scores with genome wide polygenic scores',
                                 'ES-PGS results from EpiSLI cohort factor scores',
                                 'Validation of HAQER CP-PGS in larger SPARK cohort',
                                 'Validation of HAQER CP-PGS in SPARK research match cohort',
                                 'SPARK WGS rare reversion analysis results',
                                 'ABCD HAQER CP-PGS validation analysis and birth related traits',
                                 'Polygenic selection results from ancient human data from the AADR cohort',
                                 'ES-PGS comparison results of AADR groups to 1,000 Genomes Europeans',
                                 'Cross-species HAQER-like similarity and phenotypic data',
                                 'Phylogenetic regression results from the cross-species analysis',
                                 'Enrichment of pre/postnatal scQTLs in HAQERs',
                                 'Enrichment of birth head circumference GWAS loci in HAQERs',
                                 'Enrichment of mammalian vocal learning enhancer regions in HAQERs',
                                 "Results from hominin-gained TFBS analysis in relation to EpiSLI core language (F1)",
                                 "Enrichment results of TF families with both hominin-gained TF motif integrity and positive motif integrity-language associations",
                                 "EpiSLI whole genome sequencing coverage statistics"
                                 ))
export(list(README = tb,
            SupplementaryTable1 = cbcl_res, 
            SupplementaryTable2 = ldpred2_res, 
            SupplementaryTable3 = es_pgs_res, 
            SupplementaryTable4 = spark_es_pgs, 
            SupplementaryTable5 = spark_self_dx_es_pgs, 
            SupplementaryTable6 = spark_rev, 
            SupplementaryTable7 = abcd_es_pgs, 
            SupplementaryTable8 = aadr_sel, 
            SupplementaryTable9 = nean_pgs, 
            SupplementaryTable10 = cs_dat, 
            SupplementaryTable11 = cs_res, 
            SupplementaryTable12 = scqtl_enr, 
            SupplementaryTable13 = bhc_enr, 
            SupplementaryTable14 = vl_enr, 
            SupplementaryTable15 = tfbs_sel, 
            SupplementaryTable16 = tfbs_conv,
            SupplementaryTable17 = episli_wgs_cov),
        out)