#!/bin/bash

module load bcftools
which bcftools

## ==============================================
## subset VCF to HAQERs
## ==============================================
# HAQER="/Dedicated/jmichaelson-wdata/lcasten/tools/ref_data/human_accelerated_regions/HAQERs/Mangan-Cell2022/HAQER_10Kb_flank.no_chr_pre.hg19.bed"
## remove "chr" prefix on chromosome column in bed file
HAQER="/Dedicated/jmichaelson-wdata/lcasten/sli_wgs/prs/gene_sets/HAQERs_v2.hg19.bed"
HAQER_NO_CHR="/Dedicated/jmichaelson-wdata/lcasten/sli_wgs/HAQERs_v2.hg19.no_chr.bed"
sed 's/^chr//' $HAQER > $HAQER_NO_CHR

## merged high quality neanderthal whole genomes
IV="/Dedicated/jmichaelson-wdata/lcasten/tools/ref_data/archaic_hominini/hg19_neanderthal_filtered_vcf/vcf/merged_all_chromosomes.vcf.gz"

## output
OV=${IV%.vcf.gz}.HAQERs_v2.vcf.gz

## subset to HAQERs
bcftools view -R $HAQER_NO_CHR ${IV} --threads 12 -Ou | bcftools +fill-tags -Oz -o $OV -- -t AC,AN,AF,HWE
tabix -p vcf $OV

## check file
echo "Head of subset VCF:"
zcat $OV | grep -v "^#" | head
echo -e "\n\n\n"
echo "Number of loci in VCF:"
zcat $OV | grep -v "^#" | wc -l

## extract alleles
bcftools query -f '%CHROM %POS %REF %ALT %AN %AF\n' $OV > ${IV%.vcf.gz}.HAQERs_v2.alleles.txt
echo -e "\n\n\n"
echo "Head of allele info from VCF:"
head ${IV%.vcf.gz}.HAQERs_v2.alleles.txt
## 