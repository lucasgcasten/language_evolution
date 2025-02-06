#!/bin/bash

## =============================================================================================================================================================
## this script will walk you through how to compute ES-PGS for 2 example annotations (HARs and HAQERs)
## this script relies on great publicly available software like "bedtools" and "PRSet"
## to run this with your data just replace the "target" path in PRSet and in the Rscripts (only if your data is QC'd and ready for polygenic score calculations)
## whenever running ES-PGS please check the genome builds match across the annotation and your genetic data
## if the genome builds don't match, try lifting the coordinates with a tool like UCSC's liftOver
## please reach out to Lucas with any questions (Lucas-casten@uiowa.edu or Lucasgcasten@gmail.com)
## =============================================================================================================================================================

## step 1: gather annotations for ES-PGS and identify "complement" regions for the background PGS
arr=("ES-PGS_example/example_data/HAQERs.bed" "ES-PGS_example/example_data/HARs.bed") ## annotations bed files
for i in ${arr[@]} 
do
    ## print progress
    echo "============================================================="
    echo "Identifying background regions for: $(basename ${i})"
    echo "Complement regions will be written out to here: ES-PGS_example/example_data/complement_$(basename ${i})"
    echo "============================================================="
    ## print some basic summary stats about the annotation
    echo "There are $(wc -l < ${i}) regions in this annotation"
    echo "This annotation accounts for $(awk '{sum += ($3 - $2)} END {print sum}' ${i}) bases in the genome"
    ## identify complementary regions
    bedtools complement -i ${i} -g ES-PGS_example/example_data/chromosome_sizes.txt > ES-PGS_example/example_data/complement_$(basename ${i})
    ## basic summary of the complement
    echo "There are $(wc -l < ES-PGS_example/example_data/complement_$(basename ${i})) regions in this complement of this annotation"
    ##
    printf "\n\n\n"
done 

## step 2: compute ES-PGS with PRSet
#### replace the file paths to match your data (genotype + sumstats data)
Rscript /wdata/lcasten/tools/PRSice/PRSice.R --prsice /wdata/lcasten/tools/PRSice/PRSice_linux \
    --base /wdata/lcasten/sli_wgs/polygenic_neanderthal_introgression/sumstats/cogPerf.tsv \
    --target /wdata/lcasten/sli_wgs/prs/1000_genomes_EUR_merge/merged.1000_genomes_EUR.SLI_WGS.qc \
    --ld-keep /wdata/lcasten/sli_wgs/prs/1000_genomes_EUR_merge/LD_keep_samples.txt \
    --binary-target T \
    --maf 0.01 \
    --bar-levels 1 \
    --fastscore --no-regress \
    --thread 10 \
    --bed ES-PGS_example/example_data/HAQERs.bed:HAQERs,data/complement_HAQERs.bed:complement_HAQERs,data/HARs.bed:HARs,data/complement_HARs.bed:complement_HARs \
    --out ES-PGS_example/example_data/ES-PGS_raw

## step 3: adjust PGS for population stratification
Rscript ES-PGS_example/code/Rscripts/correct_PGS.R

## step 4: phenotypic associations with ES-PGS
Rscript ES-PGS_example/code/Rscripts/analyze_PGS.R

## step 5: visualize ES-PGS results
Rscript ES-PGS_example/code/Rscripts/plot_ES-PGS_results.R