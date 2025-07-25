#!/bin/bash

## ---------------------------------------------------------------------
## script to QC merged AADR v54, 1000 Genomes Eur, EpiSLI dataset
## ---------------------------------------------------------------------

## prune SNPs
/Dedicated/jmichaelson-wdata/lcasten/tools/plink \
    --bfile /Dedicated/jmichaelson-wdata/lcasten/sli_wgs/ancient_DNA/data/merged \
    --indep-pairwise 1000 1 0.05 \
    --maf .05 \
    --keep /Dedicated/jmichaelson-wdata/lcasten/sli_wgs/prs/1000_genomes_EUR_merge/LD_keep_samples.txt \
    --out /Dedicated/jmichaelson-wdata/lcasten/sli_wgs/ancient_DNA/data/merged.1000Genomes_pruned_SNP_for_GCTA

## identify samples with lower missingness
/Dedicated/jmichaelson-wdata/lcasten/tools/plink \
    --bfile /Dedicated/jmichaelson-wdata/lcasten/sli_wgs/ancient_DNA/data/merged \
    --extract /Dedicated/jmichaelson-wdata/lcasten/sli_wgs/ancient_DNA/data/merged.1000Genomes_pruned_SNP_for_GCTA.prune.in \
    --mind 0.5 \
    --make-just-fam \
    --out /Dedicated/jmichaelson-wdata/lcasten/sli_wgs/ancient_DNA/data/merged.1000Genomes_pruned_SNP_for_GCTA.samples_missing_under_50pct

## identify SNPs with lower missingness
/Dedicated/jmichaelson-wdata/lcasten/tools/plink \
    --bfile /Dedicated/jmichaelson-wdata/lcasten/sli_wgs/ancient_DNA/data/merged \
    --extract /Dedicated/jmichaelson-wdata/lcasten/sli_wgs/ancient_DNA/data/merged.1000Genomes_pruned_SNP_for_GCTA.prune.in \
    --geno .1 \
    --keep /Dedicated/jmichaelson-wdata/lcasten/sli_wgs/ancient_DNA/data/merged.1000Genomes_pruned_SNP_for_GCTA.samples_missing_under_50pct.fam \
    --make-just-bim \
    --out /Dedicated/jmichaelson-wdata/lcasten/sli_wgs/ancient_DNA/data/merged.1000Genomes_pruned_SNP_for_GCTA.SNPs_missing_under_10pct

## QCd SNP list
awk '{ print $2 }' /Dedicated/jmichaelson-wdata/lcasten/sli_wgs/ancient_DNA/data/merged.1000Genomes_pruned_SNP_for_GCTA.SNPs_missing_under_10pct.bim > /Dedicated/jmichaelson-wdata/lcasten/sli_wgs/ancient_DNA/data/merged.1000Genomes_pruned_SNP_for_GCTA.SNPs_missing_under_10pct.snplist

## GRM with samples / SNPs of interest
/Dedicated/jmichaelson-wdata/lcasten/tools/gcta/gcta-1.94.1-linux-kernel-3-x86_64/gcta-1.94.1 \
    --bfile /Dedicated/jmichaelson-wdata/lcasten/sli_wgs/ancient_DNA/data/merged \
    --thread-num 50 \
    --keep /Dedicated/jmichaelson-wdata/lcasten/sli_wgs/ancient_DNA/data/merged.1000Genomes_pruned_SNP_for_GCTA.samples_missing_under_50pct.fam \
    --extract /Dedicated/jmichaelson-wdata/lcasten/sli_wgs/ancient_DNA/data/merged.1000Genomes_pruned_SNP_for_GCTA.SNPs_missing_under_10pct.snplist \
    --make-grm --out /Dedicated/jmichaelson-wdata/lcasten/sli_wgs/ancient_DNA/data/merged.GCTA_GRM_pruned_SNPs_QC

## get IDs for samples to use in selection analysis
awk '{ print $2 }' /Dedicated/jmichaelson-wdata/lcasten/sli_wgs/ancient_DNA/data/keep_ancient_Eur_sample_age_for_sel.lenient.pheno > kp_sample.txt
grep -f kp_sample.txt /Dedicated/jmichaelson-wdata/lcasten/sli_wgs/ancient_DNA/data/merged.1000Genomes_pruned_SNP_for_GCTA.samples_missing_under_50pct.fam > kp_sample_QC.txt

/Dedicated/jmichaelson-wdata/lcasten/tools/gcta/gcta-1.94.1-linux-kernel-3-x86_64/gcta-1.94.1 \
    --grm /Dedicated/jmichaelson-wdata/lcasten/sli_wgs/ancient_DNA/data/merged.GCTA_GRM_pruned_SNPs_QC \
    --thread-num 50 \
    --keep kp_sample_QC.txt \
    --grm-singleton 0.9 \
    --out /Dedicated/jmichaelson-wdata/lcasten/sli_wgs/ancient_DNA/data/keep_ancient_Eur_sample_age_for_sel.no_dup

ls /Dedicated/jmichaelson-wdata/lcasten/sli_wgs/ancient_DNA/data/keep_ancient_Eur_sample_age_for_sel.no_dup*

## make plain text GR file to read into R
/Dedicated/jmichaelson-wdata/lcasten/tools/gcta/gcta-1.94.1-linux-kernel-3-x86_64/gcta-1.94.1 \
    --grm /Dedicated/jmichaelson-wdata/lcasten/sli_wgs/ancient_DNA/data/merged.GCTA_GRM_pruned_SNPs_QC \
    --thread-num 10 \
    --keep /Dedicated/jmichaelson-wdata/lcasten/sli_wgs/ancient_DNA/data/keep_ancient_Eur_sample_age_for_sel.no_dup.singleton.txt \
    --make-grm-gz --out /Dedicated/jmichaelson-wdata/lcasten/sli_wgs/ancient_DNA/data/merged.GCTA_GRM_pruned_SNPs_QC.no_dup

## polygenic selection in GCTA, should give same results as the "gaston" R implementation
# /Dedicated/jmichaelson-wdata/lcasten/tools/gcta/gcta-1.94.1-linux-kernel-3-x86_64/gcta-1.94.1 \
#     --reml \
#     --grm /Dedicated/jmichaelson-wdata/lcasten/sli_wgs/ancient_DNA/data/merged.GCTA_GRM_pruned_SNPs_QC \
#     --pheno /Dedicated/jmichaelson-wdata/lcasten/sli_wgs/ancient_DNA/data/keep_ancient_Eur_sample_age_for_sel.lenient.pheno \
#     --qcovar /Dedicated/jmichaelson-wdata/lcasten/sli_wgs/ancient_DNA/data/CP_ES_PGS_all.txt \
#     --thread-num 50 \
#     --keep /Dedicated/jmichaelson-wdata/lcasten/sli_wgs/ancient_DNA/data/merged.1000Genomes_pruned_SNP_for_GCTA.samples_missing_under_50pct.fam \
#     --grm-cutoff 0.9 \
#     --reml-no-constrain \
#     --reml-est-fix \
#     --out /Dedicated/jmichaelson-wdata/lcasten/sli_wgs/ancient_DNA/data/GCTA_ES-PGS_selection_results
