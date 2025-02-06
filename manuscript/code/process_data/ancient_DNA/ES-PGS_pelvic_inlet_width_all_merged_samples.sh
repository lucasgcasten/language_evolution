#!/bin/bash

## bash /Dedicated/jmichaelson-wdata/lcasten/sli_wgs/ancient_DNA/code/ES-PGS_pelvic_inlet_width_all_merged_samples.sh

eval "$(conda shell.bash hook)" ## run when calling from bash

conda activate RERCONVERGE

cd /Dedicated/jmichaelson-wdata/lcasten/tools/PRSice

## pelvic inlet width
Rscript PRSice.R --prsice PRSice_linux \
    --base /Dedicated/jmichaelson-sdata/gwas_summary_stats/anthropomorphic/pelvis/Xu_biorxiv2024/pelvic_inlet_width_lmm_combo.tab \
    --a1 ALLELE1 \
    --a2 ALLELE0 \
    --pvalue P_BOLT_LMM_INF \
    --target /Dedicated/jmichaelson-wdata/lcasten/sli_wgs/ancient_DNA/data/merged \
    --ld-keep /Dedicated/jmichaelson-wdata/lcasten/sli_wgs/prs/1000_genomes_EUR_merge/LD_keep_samples.txt \
    --binary-target T \
    --bar-levels 1 --fastscore --no-regress --thread 54 \
    --bed /Dedicated/jmichaelson-wdata/lcasten/sli_wgs/prs/gene_sets/HAQERs.bed:HAQER,/Dedicated/jmichaelson-wdata/lcasten/sli_wgs/prs/gene_sets/complement_HAQERs.bed:complement_HAQER \
    --print-snp \
    --out /Dedicated/jmichaelson-wdata/lcasten/sli_wgs/prs/pathway_prs/1000Genomes_EpiSLI_AADR_neanderthal_merged.human_evolution_complement.female_pelvic_inlet_width

## both sex pelvic inlet width (much larger sample)
Rscript PRSice.R --prsice PRSice_linux \
    --base /Dedicated/jmichaelson-sdata/gwas_summary_stats/anthropomorphic/pelvis/Xu_biorxiv2024/both_sex_pelvic_inlet_width_lmm_combo.tab \
    --a1 ALLELE1 \
    --a2 ALLELE0 \
    --pvalue P_BOLT_LMM_INF \
    --target /Dedicated/jmichaelson-wdata/lcasten/sli_wgs/ancient_DNA/data/merged \
    --ld-keep /Dedicated/jmichaelson-wdata/lcasten/sli_wgs/prs/1000_genomes_EUR_merge/LD_keep_samples.txt \
    --binary-target T \
    --bar-levels 1 --fastscore --no-regress --thread 54 \
    --bed /Dedicated/jmichaelson-wdata/lcasten/sli_wgs/prs/gene_sets/HAQERs.bed:HAQER,/Dedicated/jmichaelson-wdata/lcasten/sli_wgs/prs/gene_sets/complement_HAQERs.bed:complement_HAQER \
    --print-snp \
    --out /Dedicated/jmichaelson-wdata/lcasten/sli_wgs/prs/pathway_prs/1000Genomes_EpiSLI_AADR_neanderthal_merged.human_evolution_complement.pelvic_inlet_width

## waist circumference adj for BMI in females (much larger GWAS than preprint)
# zcat /sdata/gwas_summary_stats/anthropomorphic/pelvis/female_waist_circumference_BMI_adj.GCST90020029.tsv.gz | head
Rscript PRSice.R --prsice PRSice_linux \
    --base /Dedicated/jmichaelson-sdata/gwas_summary_stats/anthropomorphic/pelvis/female_waist_circumference_BMI_adj.GCST90020029.tsv.gz \
    --a1 effect_allele \
    --a2 other_allele \
    --pvalue p_value \
    --snp rs_id \
    --target /Dedicated/jmichaelson-wdata/lcasten/sli_wgs/ancient_DNA/data/merged \
    --ld-keep /Dedicated/jmichaelson-wdata/lcasten/sli_wgs/prs/1000_genomes_EUR_merge/LD_keep_samples.txt \
    --binary-target T \
    --bar-levels 1 --fastscore --no-regress --thread 54 \
    --bed /Dedicated/jmichaelson-wdata/lcasten/sli_wgs/prs/gene_sets/HAQERs.bed:HAQER,/Dedicated/jmichaelson-wdata/lcasten/sli_wgs/prs/gene_sets/complement_HAQERs.bed:complement_HAQER \
    --print-snp \
    --extract /Dedicated/jmichaelson-wdata/lcasten/sli_wgs/prs/pathway_prs/1000Genomes_EpiSLI_AADR_neanderthal_merged.human_evolution_complement.female_waist_circumference_bmi_adj.valid \
    --out /Dedicated/jmichaelson-wdata/lcasten/sli_wgs/prs/pathway_prs/1000Genomes_EpiSLI_AADR_neanderthal_merged.human_evolution_complement.female_waist_circumference_bmi_adj

## pelvic width
Rscript PRSice.R --prsice PRSice_linux \
    --base /Dedicated/jmichaelson-sdata/gwas_summary_stats/anthropomorphic/pelvis/Xu_biorxiv2024/pelvic_width_lmm_combo.tab \
    --a1 ALLELE1 \
    --a2 ALLELE0 \
    --pvalue P_BOLT_LMM_INF \
    --target /Dedicated/jmichaelson-wdata/lcasten/sli_wgs/ancient_DNA/data/merged \
    --ld-keep /Dedicated/jmichaelson-wdata/lcasten/sli_wgs/prs/1000_genomes_EUR_merge/LD_keep_samples.txt \
    --binary-target T \
    --bar-levels 1 --fastscore --no-regress --thread 54 \
    --bed /Dedicated/jmichaelson-wdata/lcasten/sli_wgs/prs/gene_sets/HAQERs.bed:HAQER,/Dedicated/jmichaelson-wdata/lcasten/sli_wgs/prs/gene_sets/complement_HAQERs.bed:complement_HAQER \
    --print-snp \
    --out /Dedicated/jmichaelson-wdata/lcasten/sli_wgs/prs/pathway_prs/1000Genomes_EpiSLI_AADR_neanderthal_merged.human_evolution_complement.female_pelvic_width

## birth weight
Rscript PRSice.R --prsice PRSice_linux \
    --base /Dedicated/jmichaelson-sdata/gwas_summary_stats/anthropomorphic/birth_early_development/EGG_consortium/Fetal_BW_European_meta.NG2019.txt \
    --a1 ea \
    --a2 nea \
    --pvalue p \
    --snp rsid \
    --target /Dedicated/jmichaelson-wdata/lcasten/sli_wgs/ancient_DNA/data/merged \
    --ld-keep /Dedicated/jmichaelson-wdata/lcasten/sli_wgs/prs/1000_genomes_EUR_merge/LD_keep_samples.txt \
    --binary-target T \
    --bar-levels 1 --fastscore --no-regress --thread 54 \
    --bed /Dedicated/jmichaelson-wdata/lcasten/sli_wgs/prs/gene_sets/HAQERs.bed:HAQER,/Dedicated/jmichaelson-wdata/lcasten/sli_wgs/prs/gene_sets/complement_HAQERs.bed:complement_HAQER \
    --print-snp \
    --extract /Dedicated/jmichaelson-wdata/lcasten/sli_wgs/prs/pathway_prs/1000Genomes_EpiSLI_AADR_neanderthal_merged.human_evolution_complement.fetal_birth_weight.valid \
    --out /Dedicated/jmichaelson-wdata/lcasten/sli_wgs/prs/pathway_prs/1000Genomes_EpiSLI_AADR_neanderthal_merged.human_evolution_complement.fetal_birth_weight


## placental weight
Rscript PRSice.R --prsice PRSice_linux \
    --base /Dedicated/jmichaelson-sdata/gwas_summary_stats/anthropomorphic/birth_early_development/EGG_consortium/placental_weight_reformatted.txt \
    --snp snp \
    --a1 A1 \
    --a2 A0 \
    --pvalue p \
    --target /Dedicated/jmichaelson-wdata/lcasten/sli_wgs/ancient_DNA/data/merged \
    --ld-keep /Dedicated/jmichaelson-wdata/lcasten/sli_wgs/prs/1000_genomes_EUR_merge/LD_keep_samples.txt \
    --binary-target T \
    --bar-levels 1 --fastscore --no-regress --thread 54 \
    --bed /Dedicated/jmichaelson-wdata/lcasten/sli_wgs/prs/gene_sets/HAQERs.bed:HAQER,/Dedicated/jmichaelson-wdata/lcasten/sli_wgs/prs/gene_sets/complement_HAQERs.bed:complement_HAQER \
    --print-snp \
    --extract /Dedicated/jmichaelson-wdata/lcasten/sli_wgs/prs/pathway_prs/1000Genomes_EpiSLI_AADR_neanderthal_merged.human_evolution_complement.placental_weight.valid \
    --out /Dedicated/jmichaelson-wdata/lcasten/sli_wgs/prs/pathway_prs/1000Genomes_EpiSLI_AADR_neanderthal_merged.human_evolution_complement.placental_weight


## birth head circumference
# read_table('/Dedicated/jmichaelson-sdata/gwas_summary_stats/anthropomorphic/birth_early_development/EGG_consortium/HC_BIRTH.filteredDf0.txt') %>% 
#     rename(snp = ID, p = P.value, chr = CHR, pos = POS) %>% 
#     mutate(Allele1 = toupper(Allele1), Allele2 = toupper(Allele2)) %>% 
#     drop_na(chr) %>% 
#     select(snp, chr, pos, A1 = Allele1, A0 = Allele2, beta = Effect, se = StdErr, p, N = TotalSampleSize, Freq1) %>% 
#     filter(Freq1 >= 0.01) %>% 
#     filter(snp != '.') %>%
#     arrange(chr, pos) %>%
#     write_tsv('/Dedicated/jmichaelson-sdata/gwas_summary_stats/anthropomorphic/birth_early_development/EGG_consortium/head_circumference_at_birth_reformatted.txt')
Rscript PRSice.R --prsice PRSice_linux \
    --base /Dedicated/jmichaelson-sdata/gwas_summary_stats/anthropomorphic/birth_early_development/EGG_consortium/head_circumference_at_birth_reformatted.txt \
    --a1 A1 \
    --a2 A0 \
    --pvalue p \
    --snp snp \
    --target /Dedicated/jmichaelson-wdata/lcasten/sli_wgs/ancient_DNA/data/merged \
    --ld-keep /Dedicated/jmichaelson-wdata/lcasten/sli_wgs/prs/1000_genomes_EUR_merge/LD_keep_samples.txt \
    --binary-target T \
    --bar-levels 1 --fastscore --no-regress --thread 54 \
    --bed /Dedicated/jmichaelson-wdata/lcasten/sli_wgs/prs/gene_sets/HAQERs.bed:HAQER,/Dedicated/jmichaelson-wdata/lcasten/sli_wgs/prs/gene_sets/complement_HAQERs.bed:complement_HAQER \
    --print-snp \
    --out /Dedicated/jmichaelson-wdata/lcasten/sli_wgs/prs/pathway_prs/1000Genomes_EpiSLI_AADR_neanderthal_merged.human_evolution_complement.birth_head_circumference


## infant head circumference
# read_table('/Dedicated/jmichaelson-sdata/gwas_summary_stats/anthropomorphic/birth_early_development/EGG_consortium/HC_INFANT.filteredDf0.txt') %>% 
#     rename(snp = ID, p = P.value, chr = CHR, pos = POS) %>% 
#     mutate(Allele1 = toupper(Allele1), Allele2 = toupper(Allele2)) %>% 
#     drop_na(chr) %>% 
#     select(snp, chr, pos, A1 = Allele1, A0 = Allele2, beta = Effect, se = StdErr, p, N = TotalSampleSize, Freq1) %>% 
#     filter(Freq1 >= 0.01) %>% 
#     filter(snp != '.') %>%
#     arrange(chr, pos) %>%
#     write_tsv('/Dedicated/jmichaelson-sdata/gwas_summary_stats/anthropomorphic/birth_early_development/EGG_consortium/head_circumference_infancy_reformatted.txt')
Rscript PRSice.R --prsice PRSice_linux \
    --base /Dedicated/jmichaelson-sdata/gwas_summary_stats/anthropomorphic/birth_early_development/EGG_consortium/head_circumference_infancy_reformatted.txt \
    --a1 A1 \
    --a2 A0 \
    --pvalue p \
    --snp snp \
    --target /Dedicated/jmichaelson-wdata/lcasten/sli_wgs/ancient_DNA/data/merged \
    --ld-keep /Dedicated/jmichaelson-wdata/lcasten/sli_wgs/prs/1000_genomes_EUR_merge/LD_keep_samples.txt \
    --binary-target T \
    --bar-levels 1 --fastscore --no-regress --thread 54 \
    --bed /Dedicated/jmichaelson-wdata/lcasten/sli_wgs/prs/gene_sets/HAQERs.bed:HAQER,/Dedicated/jmichaelson-wdata/lcasten/sli_wgs/prs/gene_sets/complement_HAQERs.bed:complement_HAQER \
    --print-snp \
    --out /Dedicated/jmichaelson-wdata/lcasten/sli_wgs/prs/pathway_prs/1000Genomes_EpiSLI_AADR_neanderthal_merged.human_evolution_complement.infant_head_circumference



## waist circumference (BMI adj)
# tmp <-  read_table('/Dedicated/jmichaelson-sdata/gwas_summary_stats/anthropomorphic/GIANT_2015_WCadjBMI_COMBINED_EUR.no_header.txt', skip = 8, col_names = FALSE)
# names(tmp) <- c('snp', 'chr', 'pos', 'A1', 'A2', 'af', 'beta', 'se', 'p', 'N') 
# tmp %>% 
#     drop_na(chr) %>% 
#     filter(af >= 0.01) %>% 
#     filter(snp != '.') %>%
#     arrange(chr, pos) %>%
#     write_tsv('/Dedicated/jmichaelson-sdata/gwas_summary_stats/anthropomorphic/GIANT_2015_HIPadjBMI_COMBINED_EUR_reformatted.txt')
# Rscript PRSice.R --prsice PRSice_linux \
#     --base /Dedicated/jmichaelson-sdata/gwas_summary_stats/anthropomorphic/GIANT_2015_HIPadjBMI_COMBINED_EUR.txt.gz \
#     --a1 A1 \
#     --a2 A0 \
#     --pvalue p \
#     --snp snp \
#     --target /Dedicated/jmichaelson-wdata/lcasten/sli_wgs/ancient_DNA/data/merged \
#     --ld-keep /Dedicated/jmichaelson-wdata/lcasten/sli_wgs/prs/1000_genomes_EUR_merge/LD_keep_samples.txt \
#     --binary-target T \
#     --bar-levels 1 --fastscore --no-regress --thread 54 \
#     --bed /Dedicated/jmichaelson-wdata/lcasten/sli_wgs/prs/gene_sets/HAQERs.bed:HAQER,/Dedicated/jmichaelson-wdata/lcasten/sli_wgs/prs/gene_sets/complement_HAQERs.bed:complement_HAQER \
#     --print-snp \
#     --out /Dedicated/jmichaelson-wdata/lcasten/sli_wgs/prs/pathway_prs/1000Genomes_EpiSLI_AADR_neanderthal_merged.human_evolution_complement.infant_head_circumference


########
## gestational duration
# head /Dedicated/jmichaelson-sdata/gwas_summary_stats/birth/gestation_duration-Liu-NatComm2019/Fetal_gest_duration_NComms2019.txt
# read_table('/Dedicated/jmichaelson-sdata/gwas_summary_stats/birth/gestation_duration-Liu-NatComm2019/Fetal_gest_duration_NComms2019.txt') %>% 
#     select(-SNP) %>% 
#     relocate(snp = Rsid) %>% 
#     janitor::clean_names() %>% 
#     rename(beta = effect, se = std_err, a1 = effect_allele, a0 = non_effect_allele) %>% 
#     select(-het_p_val) %>% 
#     write_tsv('/Dedicated/jmichaelson-sdata/gwas_summary_stats/birth/gestation_duration-Liu-NatComm2019/Fetal_gest_duration_NComms2019.reformatted.txt')
head /Dedicated/jmichaelson-sdata/gwas_summary_stats/birth/gestation_duration-Liu-NatComm2019/Fetal_gest_duration_NComms2019.reformatted.txt

Rscript PRSice.R --prsice PRSice_linux \
    --base /Dedicated/jmichaelson-sdata/gwas_summary_stats/birth/gestation_duration-Liu-NatComm2019/Fetal_gest_duration_NComms2019.reformatted.txt \
    --a1 a1 \
    --a2 a0 \
    --snp snp \
    --pvalue p \
    --target /Dedicated/jmichaelson-wdata/lcasten/sli_wgs/ancient_DNA/data/merged \
    --ld-keep /Dedicated/jmichaelson-wdata/lcasten/sli_wgs/prs/1000_genomes_EUR_merge/LD_keep_samples.txt \
    --binary-target T \
    --bar-levels 1 --fastscore --no-regress --thread 54 \
    --bed /Dedicated/jmichaelson-wdata/lcasten/sli_wgs/prs/gene_sets/HAQERs.bed:HAQER,/Dedicated/jmichaelson-wdata/lcasten/sli_wgs/prs/gene_sets/complement_HAQERs.bed:complement_HAQER \
    --print-snp \
    --extract /Dedicated/jmichaelson-wdata/lcasten/sli_wgs/prs/pathway_prs/1000Genomes_EpiSLI_AADR_neanderthal_merged.human_evolution_complement.gestational_duration.valid \
    --out /Dedicated/jmichaelson-wdata/lcasten/sli_wgs/prs/pathway_prs/1000Genomes_EpiSLI_AADR_neanderthal_merged.human_evolution_complement.gestational_duration


###########################3
## genlang phenos
###########################
ls /Dedicated/jmichaelson-sdata/gwas_summary_stats/genlang/pnas_2022_reading_gwas
########
## NVIQ
head /Dedicated/jmichaelson-sdata/gwas_summary_stats/genlang/pnas_2022_reading_gwas/METAANALYSIS_PIQ_Z_combined_STERR_GCON_1.tbl

# read_table('/Dedicated/jmichaelson-sdata/gwas_summary_stats/genlang/pnas_2022_reading_gwas/METAANALYSIS_PIQ_Z_combined_STERR_GCON_1.tbl') %>% 
#     relocate() %>% 
#     janitor::clean_names() %>% 
#     rename(snp = marker_name, beta = effect, se = std_err, a1 = allele1, a0 = allele2, p = p_value) %>% 
#     write_tsv('/Dedicated/jmichaelson-sdata/gwas_summary_stats/genlang/pnas_2022_reading_gwas/METAANALYSIS_PIQ_Z_combined_STERR_GCON_1.reformatted.txt')

Rscript PRSice.R --prsice PRSice_linux \
    --base /Dedicated/jmichaelson-sdata/gwas_summary_stats/genlang/pnas_2022_reading_gwas/METAANALYSIS_PIQ_Z_combined_STERR_GCON_1.reformatted.txt \
    --a1 a1 \
    --a2 a0 \
    --snp snp \
    --pvalue p \
    --target /Dedicated/jmichaelson-wdata/lcasten/sli_wgs/ancient_DNA/data/merged \
    --ld-keep /Dedicated/jmichaelson-wdata/lcasten/sli_wgs/prs/1000_genomes_EUR_merge/LD_keep_samples.txt \
    --binary-target T \
    --bar-levels 1 --fastscore --no-regress --thread 54 \
    --bed /Dedicated/jmichaelson-wdata/lcasten/sli_wgs/prs/gene_sets/HAQERs.bed:HAQER,/Dedicated/jmichaelson-wdata/lcasten/sli_wgs/prs/gene_sets/complement_HAQERs.bed:complement_HAQER \
    --print-snp \
    --out /Dedicated/jmichaelson-wdata/lcasten/sli_wgs/prs/pathway_prs/1000Genomes_EpiSLI_AADR_neanderthal_merged.human_evolution_complement.genlang_NVIQ

########
## NREAD
ls /Dedicated/jmichaelson-sdata/gwas_summary_stats/genlang/pnas_2022_reading_gwas
head /Dedicated/jmichaelson-sdata/gwas_summary_stats/genlang/pnas_2022_reading_gwas/METAANALYSIS_NREAD_RT_EUR_combined_STERR_GCON_1.tbl

# read_table('/Dedicated/jmichaelson-sdata/gwas_summary_stats/genlang/pnas_2022_reading_gwas/METAANALYSIS_NREAD_RT_EUR_combined_STERR_GCON_1.tbl') %>% 
#     relocate() %>% 
#     janitor::clean_names() %>% 
#     rename(snp = marker_name, beta = effect, se = std_err, a1 = allele1, a0 = allele2, p = p_value) %>% 
#     write_tsv('/Dedicated/jmichaelson-sdata/gwas_summary_stats/genlang/pnas_2022_reading_gwas/METAANALYSIS_NREAD_RT_EUR_combined_STERR_GCON_1.reformatted.txt')

Rscript PRSice.R --prsice PRSice_linux \
    --base /Dedicated/jmichaelson-sdata/gwas_summary_stats/genlang/pnas_2022_reading_gwas/METAANALYSIS_NREAD_RT_EUR_combined_STERR_GCON_1.reformatted.txt \
    --a1 a1 \
    --a2 a0 \
    --snp snp \
    --pvalue p \
    --target /Dedicated/jmichaelson-wdata/lcasten/sli_wgs/ancient_DNA/data/merged \
    --ld-keep /Dedicated/jmichaelson-wdata/lcasten/sli_wgs/prs/1000_genomes_EUR_merge/LD_keep_samples.txt \
    --binary-target T \
    --bar-levels 1 --fastscore --no-regress --thread 54 \
    --bed /Dedicated/jmichaelson-wdata/lcasten/sli_wgs/prs/gene_sets/HAQERs.bed:HAQER,/Dedicated/jmichaelson-wdata/lcasten/sli_wgs/prs/gene_sets/complement_HAQERs.bed:complement_HAQER \
    --print-snp \
    --out /Dedicated/jmichaelson-wdata/lcasten/sli_wgs/prs/pathway_prs/1000Genomes_EpiSLI_AADR_neanderthal_merged.human_evolution_complement.genlang_NREAD

########
## PA
ls /Dedicated/jmichaelson-sdata/gwas_summary_stats/genlang/pnas_2022_reading_gwas
head /Dedicated/jmichaelson-sdata/gwas_summary_stats/genlang/pnas_2022_reading_gwas/METAANALYSIS_PA_RT_EUR_combined_STERR_GCON_1.tbl

# read_table('/Dedicated/jmichaelson-sdata/gwas_summary_stats/genlang/pnas_2022_reading_gwas/METAANALYSIS_PA_RT_EUR_combined_STERR_GCON_1.tbl') %>% 
#     relocate() %>% 
#     janitor::clean_names() %>% 
#     rename(snp = marker_name, beta = effect, se = std_err, a1 = allele1, a0 = allele2, p = p_value) %>% 
#     write_tsv('/Dedicated/jmichaelson-sdata/gwas_summary_stats/genlang/pnas_2022_reading_gwas/METAANALYSIS_PA_RT_EUR_combined_STERR_GCON_1.reformatted.txt')

Rscript PRSice.R --prsice PRSice_linux \
    --base /Dedicated/jmichaelson-sdata/gwas_summary_stats/genlang/pnas_2022_reading_gwas/METAANALYSIS_PA_RT_EUR_combined_STERR_GCON_1.reformatted.txt \
    --a1 a1 \
    --a2 a0 \
    --snp snp \
    --pvalue p \
    --target /Dedicated/jmichaelson-wdata/lcasten/sli_wgs/ancient_DNA/data/merged \
    --ld-keep /Dedicated/jmichaelson-wdata/lcasten/sli_wgs/prs/1000_genomes_EUR_merge/LD_keep_samples.txt \
    --binary-target T \
    --bar-levels 1 --fastscore --no-regress --thread 54 \
    --bed /Dedicated/jmichaelson-wdata/lcasten/sli_wgs/prs/gene_sets/HAQERs.bed:HAQER,/Dedicated/jmichaelson-wdata/lcasten/sli_wgs/prs/gene_sets/complement_HAQERs.bed:complement_HAQER \
    --print-snp \
    --out /Dedicated/jmichaelson-wdata/lcasten/sli_wgs/prs/pathway_prs/1000Genomes_EpiSLI_AADR_neanderthal_merged.human_evolution_complement.genlang_PA


########
## NREP
ls /Dedicated/jmichaelson-sdata/gwas_summary_stats/genlang/pnas_2022_reading_gwas
head /Dedicated/jmichaelson-sdata/gwas_summary_stats/genlang/pnas_2022_reading_gwas/METAANALYSIS_RANDOM__NREP_Z_EUR_combined_STERR_GCON_1.tbl

# read_table('/Dedicated/jmichaelson-sdata/gwas_summary_stats/genlang/pnas_2022_reading_gwas/METAANALYSIS_RANDOM__NREP_Z_EUR_combined_STERR_GCON_1.tbl') %>% 
#     relocate() %>% 
#     janitor::clean_names() %>% 
#     rename(snp = marker_name, beta = effect, se = std_err, a1 = allele1, a0 = allele2, p = pvalue) %>% 
#     write_tsv('/Dedicated/jmichaelson-sdata/gwas_summary_stats/genlang/pnas_2022_reading_gwas/METAANALYSIS_NREP_RT_EUR_combined_STERR_GCON_1.reformatted.txt')

Rscript PRSice.R --prsice PRSice_linux \
    --base /Dedicated/jmichaelson-sdata/gwas_summary_stats/genlang/pnas_2022_reading_gwas/METAANALYSIS_NREP_RT_EUR_combined_STERR_GCON_1.reformatted.txt \
    --a1 a1 \
    --a2 a0 \
    --snp snp \
    --pvalue p \
    --target /Dedicated/jmichaelson-wdata/lcasten/sli_wgs/ancient_DNA/data/merged \
    --ld-keep /Dedicated/jmichaelson-wdata/lcasten/sli_wgs/prs/1000_genomes_EUR_merge/LD_keep_samples.txt \
    --binary-target T \
    --bar-levels 1 --fastscore --no-regress --thread 54 \
    --bed /Dedicated/jmichaelson-wdata/lcasten/sli_wgs/prs/gene_sets/HAQERs.bed:HAQER,/Dedicated/jmichaelson-wdata/lcasten/sli_wgs/prs/gene_sets/complement_HAQERs.bed:complement_HAQER \
    --print-snp \
    --out /Dedicated/jmichaelson-wdata/lcasten/sli_wgs/prs/pathway_prs/1000Genomes_EpiSLI_AADR_neanderthal_merged.human_evolution_complement.genlang_NREP


########
## SP
ls /Dedicated/jmichaelson-sdata/gwas_summary_stats/genlang/pnas_2022_reading_gwas
head /Dedicated/jmichaelson-sdata/gwas_summary_stats/genlang/pnas_2022_reading_gwas/METAANALYSIS_SP_RT_EUR_combined_STERR_GCON_1.tbl

# read_table('/Dedicated/jmichaelson-sdata/gwas_summary_stats/genlang/pnas_2022_reading_gwas/METAANALYSIS_SP_RT_EUR_combined_STERR_GCON_1.tbl') %>% 
#     relocate() %>% 
#     janitor::clean_names() %>% 
#     rename(snp = marker_name, beta = effect, se = std_err, a1 = allele1, a0 = allele2, p = p_value) %>% 
#     write_tsv('/Dedicated/jmichaelson-sdata/gwas_summary_stats/genlang/pnas_2022_reading_gwas/METAANALYSIS_SP_RT_EUR_combined_STERR_GCON_1.reformatted.txt')

Rscript PRSice.R --prsice PRSice_linux \
    --base /Dedicated/jmichaelson-sdata/gwas_summary_stats/genlang/pnas_2022_reading_gwas/METAANALYSIS_SP_RT_EUR_combined_STERR_GCON_1.reformatted.txt \
    --a1 a1 \
    --a2 a0 \
    --snp snp \
    --pvalue p \
    --target /Dedicated/jmichaelson-wdata/lcasten/sli_wgs/ancient_DNA/data/merged \
    --ld-keep /Dedicated/jmichaelson-wdata/lcasten/sli_wgs/prs/1000_genomes_EUR_merge/LD_keep_samples.txt \
    --binary-target T \
    --bar-levels 1 --fastscore --no-regress --thread 54 \
    --bed /Dedicated/jmichaelson-wdata/lcasten/sli_wgs/prs/gene_sets/HAQERs.bed:HAQER,/Dedicated/jmichaelson-wdata/lcasten/sli_wgs/prs/gene_sets/complement_HAQERs.bed:complement_HAQER \
    --print-snp \
    --out /Dedicated/jmichaelson-wdata/lcasten/sli_wgs/prs/pathway_prs/1000Genomes_EpiSLI_AADR_neanderthal_merged.human_evolution_complement.genlang_SP



########
## WR
ls /Dedicated/jmichaelson-sdata/gwas_summary_stats/genlang/pnas_2022_reading_gwas
head /Dedicated/jmichaelson-sdata/gwas_summary_stats/genlang/pnas_2022_reading_gwas/METAANALYSIS_WR_Z_EUR_combined_STERR_GCON_1.tbl

# read_table('/Dedicated/jmichaelson-sdata/gwas_summary_stats/genlang/pnas_2022_reading_gwas/METAANALYSIS_WR_Z_EUR_combined_STERR_GCON_1.tbl') %>% 
#     relocate() %>% 
#     janitor::clean_names() %>% 
#     rename(snp = marker_name, beta = effect, se = std_err, a1 = allele1, a0 = allele2, p = p_value) %>% 
#     write_tsv('/Dedicated/jmichaelson-sdata/gwas_summary_stats/genlang/pnas_2022_reading_gwas/METAANALYSIS_WR_Z_EUR_combined_STERR_GCON_1.reformatted.txt')

Rscript PRSice.R --prsice PRSice_linux \
    --base /Dedicated/jmichaelson-sdata/gwas_summary_stats/genlang/pnas_2022_reading_gwas/METAANALYSIS_WR_Z_EUR_combined_STERR_GCON_1.reformatted.txt \
    --a1 a1 \
    --a2 a0 \
    --snp snp \
    --pvalue p \
    --target /Dedicated/jmichaelson-wdata/lcasten/sli_wgs/ancient_DNA/data/merged \
    --ld-keep /Dedicated/jmichaelson-wdata/lcasten/sli_wgs/prs/1000_genomes_EUR_merge/LD_keep_samples.txt \
    --binary-target T \
    --bar-levels 1 --fastscore --no-regress --thread 54 \
    --bed /Dedicated/jmichaelson-wdata/lcasten/sli_wgs/prs/gene_sets/HAQERs.bed:HAQER,/Dedicated/jmichaelson-wdata/lcasten/sli_wgs/prs/gene_sets/complement_HAQERs.bed:complement_HAQER \
    --print-snp \
    --out /Dedicated/jmichaelson-wdata/lcasten/sli_wgs/prs/pathway_prs/1000Genomes_EpiSLI_AADR_neanderthal_merged.human_evolution_complement.genlang_WR


#######################
## adult brain volume meta analysis
# ss <- read_table('/sdata/gwas_summary_stats/cncr-nl/brainvol-intelligence_nagel_2020/meta_analysis_BV_Jansenetal_2020.sumstats.txt.gz') %>% 
#     janitor::clean_names() %>% 
#     mutate(beta = z / sqrt(2 * eaf_eur_ukb * (1 - eaf_eur_ukb) * (n + z^2)),
#            se = 1 / sqrt(2 * eaf_eur_ukb * (1 - eaf_eur_ukb) * (n + z^2)))
# mean(ss$beta)
# sd(ss$beta)
# hist(ss$beta)

# ss %>%
#     select(-snp) %>% 
#     rename(a0 = a2, af = eaf_eur_ukb) %>%
#     write_tsv('/sdata/gwas_summary_stats/cncr-nl/brainvol-intelligence_nagel_2020/meta_analysis_BV_Jansenetal_2020.sumstats.reformatted.txt')

Rscript PRSice.R --prsice PRSice_linux \
    --base /Dedicated/jmichaelson-sdata/gwas_summary_stats/cncr-nl/brainvol-intelligence_nagel_2020/meta_analysis_BV_Jansenetal_2020.sumstats.reformatted.txt \
    --a1 a1 \
    --a2 a0 \
    --snp rsid \
    --pvalue p \
    --target /Dedicated/jmichaelson-wdata/lcasten/sli_wgs/ancient_DNA/data/merged \
    --ld-keep /Dedicated/jmichaelson-wdata/lcasten/sli_wgs/prs/1000_genomes_EUR_merge/LD_keep_samples.txt \
    --binary-target T \
    --bar-levels 1 --fastscore --no-regress --thread 54 \
    --bed /Dedicated/jmichaelson-wdata/lcasten/sli_wgs/prs/gene_sets/HAQERs.bed:HAQER,/Dedicated/jmichaelson-wdata/lcasten/sli_wgs/prs/gene_sets/complement_HAQERs.bed:complement_HAQER \
    --print-snp \
    --out /Dedicated/jmichaelson-wdata/lcasten/sli_wgs/prs/pathway_prs/1000Genomes_EpiSLI_AADR_neanderthal_merged.human_evolution_complement.intracranial_volume_adult_meta


#####################
## brain surface area
zcat /Dedicated/jmichaelson-sdata/gwas_summary_stats/ENIGMA_brain_imaging/global_measures/raw/ENIGMA3_mixed_se_wo_Mean_Full_SurfArea_20190429.txt.gz | head
Rscript PRSice.R --prsice PRSice_linux \
    --base /Dedicated/jmichaelson-sdata/gwas_summary_stats/ENIGMA_brain_imaging/global_measures/raw/ENIGMA3_mixed_se_wo_Mean_Full_SurfArea_20190429_ZStatAdded.tsv \
    --a1 A1 \
    --a2 A2 \
    --snp SNP \
    --pvalue P \
    --target /Dedicated/jmichaelson-wdata/lcasten/sli_wgs/ancient_DNA/data/merged \
    --ld-keep /Dedicated/jmichaelson-wdata/lcasten/sli_wgs/prs/1000_genomes_EUR_merge/LD_keep_samples.txt \
    --binary-target T \
    --bar-levels 1 --fastscore --no-regress --thread 54 \
    --bed /Dedicated/jmichaelson-wdata/lcasten/sli_wgs/prs/gene_sets/HAQERs.bed:HAQER,/Dedicated/jmichaelson-wdata/lcasten/sli_wgs/prs/gene_sets/complement_HAQERs.bed:complement_HAQER \
    --print-snp \
    --out /Dedicated/jmichaelson-wdata/lcasten/sli_wgs/prs/pathway_prs/1000Genomes_EpiSLI_AADR_neanderthal_merged.human_evolution_complement.brain_surface_area

####################3
## loop over global brain imaging measures
cd /Dedicated/jmichaelson-sdata/gwas_summary_stats/brain_imaging_Warrier_NatureGenetics_2023/cleaned
files=(*.tsv)
cd /Dedicated/jmichaelson-wdata/lcasten/tools/PRSice
echo ${files[@]}

head /Dedicated/jmichaelson-sdata/gwas_summary_stats/brain_imaging_Warrier_NatureGenetics_2023/cleaned/CT_meta.tsv

for f in ${files[@]}
do 
    ## brain imagin pheno
    ph="brain_${f%_meta.tsv}"
    ## dry run (will identify duplicated SNPs)
    Rscript PRSice.R --prsice PRSice_linux \
        --base /Dedicated/jmichaelson-sdata/gwas_summary_stats/brain_imaging_Warrier_NatureGenetics_2023/cleaned/${f} \
        --chr-id c:L-aBd \
        --a1 A1 \
        --a2 A2 \
        --pvalue P \
        --target /Dedicated/jmichaelson-wdata/lcasten/sli_wgs/ancient_DNA/data/merged \
        --ld-keep /Dedicated/jmichaelson-wdata/lcasten/sli_wgs/prs/1000_genomes_EUR_merge/LD_keep_samples.txt \
        --binary-target T \
        --bar-levels 1 --fastscore --no-regress --thread 54 \
        --bed /Dedicated/jmichaelson-wdata/lcasten/sli_wgs/prs/gene_sets/HAQERs.bed:HAQER,/Dedicated/jmichaelson-wdata/lcasten/sli_wgs/prs/gene_sets/complement_HAQERs.bed:complement_HAQER \
        --print-snp \
        --out /Dedicated/jmichaelson-wdata/lcasten/sli_wgs/prs/pathway_prs/1000Genomes_EpiSLI_AADR_neanderthal_merged.human_evolution_complement.${ph}
    ## compute PRS w/o duplicated SNPs
    Rscript PRSice.R --prsice PRSice_linux \
        --base /Dedicated/jmichaelson-sdata/gwas_summary_stats/brain_imaging_Warrier_NatureGenetics_2023/cleaned/${f} \
        --chr-id c:L-aBd \
        --a1 A1 \
        --a2 A2 \
        --pvalue P \
        --target /Dedicated/jmichaelson-wdata/lcasten/sli_wgs/ancient_DNA/data/merged \
        --ld-keep /Dedicated/jmichaelson-wdata/lcasten/sli_wgs/prs/1000_genomes_EUR_merge/LD_keep_samples.txt \
        --binary-target T \
        --bar-levels 1 --fastscore --no-regress --thread 54 \
        --bed /Dedicated/jmichaelson-wdata/lcasten/sli_wgs/prs/gene_sets/HAQERs.bed:HAQER,/Dedicated/jmichaelson-wdata/lcasten/sli_wgs/prs/gene_sets/complement_HAQERs.bed:complement_HAQER \
        --print-snp \
        --extract /Dedicated/jmichaelson-wdata/lcasten/sli_wgs/prs/pathway_prs/1000Genomes_EpiSLI_AADR_neanderthal_merged.human_evolution_complement.${ph}.valid \
        --out /Dedicated/jmichaelson-wdata/lcasten/sli_wgs/prs/pathway_prs/1000Genomes_EpiSLI_AADR_neanderthal_merged.human_evolution_complement.${ph}    
done



####################3
## loop over subcortical vols
## run in R
### fix weird sumstats
# weird_ss = list.files('/sdata/gwas_summary_stats/ENIGMA_brain_imaging/subcortical_unrestricted_use', full.names = T)
# weird_ss = weird_ss[str_detect(weird_ss, '_unrestricted_NG')]
# weird_ss

# for (f in weird_ss) {
#   cat('\n\n\n\n\n')
#   message('====================================================')
#   message(str_c('reformatting ', f))
#   message('====================================================')

#   ss = read_table(file = f, col_types = c('ccccddddcd'))
  
#   ss = ss %>%
#     mutate(CHR = str_split(MarkerName, pattern = ':', simplify = T)[,1],
#            BP = str_split(MarkerName, pattern = ':', simplify = T)[,2]) %>%
#     relocate(CHR, BP) %>%
#     mutate(A1 = toupper(A1),
#            A2 = toupper(A2)
#     ) %>%
#     rename(Z = Zscore) %>%
#     mutate(BETA = Z / sqrt(2 * FreqA1 * (1 - FreqA1) * (N + Z^2)),
#            SE = 1 / sqrt(2 * FreqA1 * (1 - FreqA1) * (N + Z^2)),
#     )
  
#   outf = str_c('/sdata/gwas_summary_stats/ENIGMA_brain_imaging/subcortical_unrestricted_use/processed/', basename(f))
#   outf = str_sub(outf, start = 1, end = -4)
#   outf = str_c(outf, '.txt')
#   message(str_c('writing reformatted sumstats to: ', outf))
  
#   ss %>%
#     write_tsv(file = outf)
#   system(str_c('bgzip -f ', outf))
# }

## run 
cd /Dedicated/jmichaelson-sdata/gwas_summary_stats/ENIGMA_brain_imaging/subcortical_unrestricted_use/processed
ls
files=(*.txt.gz)
cd /Dedicated/jmichaelson-wdata/lcasten/tools/PRSice
echo ${files[@]}

zcat /Dedicated/jmichaelson-sdata/gwas_summary_stats/ENIGMA_brain_imaging/subcortical_unrestricted_use/processed/accumbens_eur_z_ldsc_unrestricted_NG05SEP19.txt.gz | head

for f in ${files[@]}
do 
    ## brain imaging pheno
    ph="brain_${f%_eur_z_ldsc_unrestricted_NG05SEP19.txt.gz}"
    ## dry run (will identify duplicated SNPs)
    Rscript PRSice.R --prsice PRSice_linux \
        --base /Dedicated/jmichaelson-sdata/gwas_summary_stats/ENIGMA_brain_imaging/subcortical_unrestricted_use/processed/${f} \
        --a1 A1 \
        --a2 A2 \
        --pvalue P \
        --snp rsid \
        --target /Dedicated/jmichaelson-wdata/lcasten/sli_wgs/ancient_DNA/data/merged \
        --ld-keep /Dedicated/jmichaelson-wdata/lcasten/sli_wgs/prs/1000_genomes_EUR_merge/LD_keep_samples.txt \
        --binary-target T \
        --bar-levels 1 --fastscore --no-regress --thread 54 \
        --bed /Dedicated/jmichaelson-wdata/lcasten/sli_wgs/prs/gene_sets/HAQERs.bed:HAQER,/Dedicated/jmichaelson-wdata/lcasten/sli_wgs/prs/gene_sets/complement_HAQERs.bed:complement_HAQER \
        --print-snp \
        --out /Dedicated/jmichaelson-wdata/lcasten/sli_wgs/prs/pathway_prs/1000Genomes_EpiSLI_AADR_neanderthal_merged.human_evolution_complement.${ph}
    ## compute PRS w/o duplicated SNPs
    Rscript PRSice.R --prsice PRSice_linux \
        --base /Dedicated/jmichaelson-sdata/gwas_summary_stats/ENIGMA_brain_imaging/subcortical_unrestricted_use/processed/${f} \
        --a1 A1 \
        --a2 A2 \
        --pvalue P \
        --snp rsid \
        --target /Dedicated/jmichaelson-wdata/lcasten/sli_wgs/ancient_DNA/data/merged \
        --ld-keep /Dedicated/jmichaelson-wdata/lcasten/sli_wgs/prs/1000_genomes_EUR_merge/LD_keep_samples.txt \
        --binary-target T \
        --bar-levels 1 --fastscore --no-regress --thread 54 \
        --bed /Dedicated/jmichaelson-wdata/lcasten/sli_wgs/prs/gene_sets/HAQERs.bed:HAQER,/Dedicated/jmichaelson-wdata/lcasten/sli_wgs/prs/gene_sets/complement_HAQERs.bed:complement_HAQER \
        --print-snp \
        --extract /Dedicated/jmichaelson-wdata/lcasten/sli_wgs/prs/pathway_prs/1000Genomes_EpiSLI_AADR_neanderthal_merged.human_evolution_complement.${ph}.valid \
        --out /Dedicated/jmichaelson-wdata/lcasten/sli_wgs/prs/pathway_prs/1000Genomes_EpiSLI_AADR_neanderthal_merged.human_evolution_complement.${ph}    
done
##


#########################
## reproductive success traits
cd /Dedicated/jmichaelson-sdata/gwas_summary_stats/birth/reproductive_success-Mathieson-NatHumBeh-2023/Mathieson_2023_02_07
ls

zcat CL_chr1to22.txt.gz | head

files=(*22.txt.gz)
cd /Dedicated/jmichaelson-wdata/lcasten/tools/PRSice
echo ${files[@]}

for f in ${files[@]}
do 
    ##  pheno
    ph="reproductive_${f%_chr1to22.txt.gz}"
    ## dry run (will identify duplicated SNPs)
    Rscript PRSice.R --prsice PRSice_linux \
        --base /Dedicated/jmichaelson-sdata/gwas_summary_stats/birth/reproductive_success-Mathieson-NatHumBeh-2023/Mathieson_2023_02_07/${f} \
        --a1 Allele1 \
        --a2 Allele2 \
        --pvalue Pvalue \
        --snp MARKER \
        --target /Dedicated/jmichaelson-wdata/lcasten/sli_wgs/ancient_DNA/data/merged \
        --ld-keep /Dedicated/jmichaelson-wdata/lcasten/sli_wgs/prs/1000_genomes_EUR_merge/LD_keep_samples.txt \
        --binary-target T \
        --bar-levels 1 --fastscore --no-regress --thread 54 \
        --bed /Dedicated/jmichaelson-wdata/lcasten/sli_wgs/prs/gene_sets/HAQERs.bed:HAQER,/Dedicated/jmichaelson-wdata/lcasten/sli_wgs/prs/gene_sets/complement_HAQERs.bed:complement_HAQER \
        --print-snp \
        --out /Dedicated/jmichaelson-wdata/lcasten/sli_wgs/prs/pathway_prs/1000Genomes_EpiSLI_AADR_neanderthal_merged.human_evolution_complement.${ph}
    # ## compute PRS w/o duplicated SNPs
    Rscript PRSice.R --prsice PRSice_linux \
        --base /Dedicated/jmichaelson-sdata/gwas_summary_stats/birth/reproductive_success-Mathieson-NatHumBeh-2023/Mathieson_2023_02_07/${f} \
        --a1 Allele1 \
        --a2 Allele2 \
        --pvalue Pvalue \
        --snp MARKER \
        --target /Dedicated/jmichaelson-wdata/lcasten/sli_wgs/ancient_DNA/data/merged \
        --ld-keep /Dedicated/jmichaelson-wdata/lcasten/sli_wgs/prs/1000_genomes_EUR_merge/LD_keep_samples.txt \
        --binary-target T \
        --bar-levels 1 --fastscore --no-regress --thread 54 \
        --bed /Dedicated/jmichaelson-wdata/lcasten/sli_wgs/prs/gene_sets/HAQERs.bed:HAQER,/Dedicated/jmichaelson-wdata/lcasten/sli_wgs/prs/gene_sets/complement_HAQERs.bed:complement_HAQER \
        --print-snp \
        --extract /Dedicated/jmichaelson-wdata/lcasten/sli_wgs/prs/pathway_prs/1000Genomes_EpiSLI_AADR_neanderthal_merged.human_evolution_complement.${ph}.valid \
        --out /Dedicated/jmichaelson-wdata/lcasten/sli_wgs/prs/pathway_prs/1000Genomes_EpiSLI_AADR_neanderthal_merged.human_evolution_complement.${ph} 
done

#########################
## reproductive success traits (women vs men)
files=("/Dedicated/jmichaelson-wdata/trthomas/abcd_spark_merge/polygenic_scores/sumstats/Social_Science_Genetic_Association_Consortium/Polygenic_Index_Repository/NEBwomen/NEBwomen1_single_gwide_sumstats.txt" "/Dedicated/jmichaelson-wdata/trthomas/abcd_spark_merge/polygenic_scores/sumstats/Social_Science_Genetic_Association_Consortium/Polygenic_Index_Repository/NEBmen/NEBmen1_single_gwide_sumstats.txt")
head ${files[0]}

cd /Dedicated/jmichaelson-wdata/lcasten/tools/PRSice

for f in ${files[@]}
do 
    ##  pheno
    ph="reproductive_$(basename ${f%1_single_gwide_sumstats.txt})"
    printf "\n\n\n\n"
    echo "==========================================="
    echo "Working on: $ph"
    echo "==========================================="
    ## dry run (will identify duplicated SNPs)
    Rscript PRSice.R --prsice PRSice_linux \
        --base ${f} \
        --a1 EFFECT_ALLELE \
        --a2 OTHER_ALLELE \
        --pvalue PVALUE \
        --snp SNPID \
        --target /Dedicated/jmichaelson-wdata/lcasten/sli_wgs/ancient_DNA/data/merged \
        --ld-keep /Dedicated/jmichaelson-wdata/lcasten/sli_wgs/prs/1000_genomes_EUR_merge/LD_keep_samples.txt \
        --binary-target T \
        --bar-levels 1 --fastscore --no-regress --thread 54 \
        --bed /Dedicated/jmichaelson-wdata/lcasten/sli_wgs/prs/gene_sets/HAQERs.bed:HAQER,/Dedicated/jmichaelson-wdata/lcasten/sli_wgs/prs/gene_sets/complement_HAQERs.bed:complement_HAQER \
        --print-snp \
        --out /Dedicated/jmichaelson-wdata/lcasten/sli_wgs/prs/pathway_prs/1000Genomes_EpiSLI_AADR_neanderthal_merged.human_evolution_complement.${ph}
done

############################
## vocal pitch
f="/Dedicated/jmichaelson-sdata/gwas_summary_stats/deCODE/vocal_gisladottir_SciAdvances2023/voicepitch2022.gz"
ph="vocal_pitch"
zcat $f | head
Rscript PRSice.R --prsice PRSice_linux \
    --base ${f} \
    --a1 EA \
    --a2 OA \
    --pvalue P \
    --snp rsID \
    --target /Dedicated/jmichaelson-wdata/lcasten/sli_wgs/ancient_DNA/data/merged \
    --ld-keep /Dedicated/jmichaelson-wdata/lcasten/sli_wgs/prs/1000_genomes_EUR_merge/LD_keep_samples.txt \
    --binary-target T \
    --bar-levels 1 --fastscore --no-regress --thread 54 \
    --bed /Dedicated/jmichaelson-wdata/lcasten/sli_wgs/prs/gene_sets/HAQERs.bed:HAQER,/Dedicated/jmichaelson-wdata/lcasten/sli_wgs/prs/gene_sets/complement_HAQERs.bed:complement_HAQER \
    --print-snp \
    --extract /Dedicated/jmichaelson-wdata/lcasten/sli_wgs/prs/pathway_prs/1000Genomes_EpiSLI_AADR_neanderthal_merged.human_evolution_complement.vocal_pitch.valid \
    --out /Dedicated/jmichaelson-wdata/lcasten/sli_wgs/prs/pathway_prs/1000Genomes_EpiSLI_AADR_neanderthal_merged.human_evolution_complement.${ph}


##########################
## waist circumference UKBB
ph="waist_circumference_GCST90302888"
Rscript PRSice.R --prsice PRSice_linux \
    --base /Dedicated/jmichaelson-sdata/gwas_summary_stats/anthropomorphic/GCST90302888.tsv.gz \
    --a1 effect_allele \
    --a2 other_allele \
    --pvalue p_value \
    --snp VARIANT_id \
    --target /Dedicated/jmichaelson-wdata/lcasten/sli_wgs/ancient_DNA/data/merged \
    --ld-keep /Dedicated/jmichaelson-wdata/lcasten/sli_wgs/prs/1000_genomes_EUR_merge/LD_keep_samples.txt \
    --binary-target T \
    --bar-levels 1 --fastscore --no-regress --thread 54 \
    --bed /Dedicated/jmichaelson-wdata/lcasten/sli_wgs/prs/gene_sets/HAQERs.bed:HAQER,/Dedicated/jmichaelson-wdata/lcasten/sli_wgs/prs/gene_sets/complement_HAQERs.bed:complement_HAQER \
    --print-snp \
    --extract /Dedicated/jmichaelson-wdata/lcasten/sli_wgs/prs/pathway_prs/1000Genomes_EpiSLI_AADR_neanderthal_merged.human_evolution_complement.waist_circumference_GCST90302888.valid \
    --out /Dedicated/jmichaelson-wdata/lcasten/sli_wgs/prs/pathway_prs/1000Genomes_EpiSLI_AADR_neanderthal_merged.human_evolution_complement.${ph}