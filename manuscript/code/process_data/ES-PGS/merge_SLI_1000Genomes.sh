#!/bin/bash

##
PLINK="/Dedicated/jmichaelson-wdata/lcasten/tools/plink"
KP="/Dedicated/jmichaelson-wdata/lcasten/sli_wgs/prs/1000_genomes_EUR_merge/overlapping_variants.qc.snplist" ## intersect VCF IDs with phenotype IDs
KP2="/Dedicated/jmichaelson-wdata/lcasten/sli_wgs/prs/1000_genomes_EUR_merge/overlapping_variants.qc_v2.snplist" ## intersect VCF IDs with phenotype IDs

SLI="/Dedicated/jmichaelson-sdata/SLI_WGS/callsets/ensemble/sli.seq.merge.common_variants.no_dups"
KG="/Dedicated/jmichaelson-sdata/1KG/phase3/EUR/common_variants.no_dups"
MERGEKG="/Dedicated/jmichaelson-wdata/lcasten/sli_wgs/prs/1000_genomes_EUR_merge/1000_genomes_EUR.qc"
MERGESLI="/Dedicated/jmichaelson-wdata/lcasten/sli_wgs/prs/1000_genomes_EUR_merge/SLI_WGS.qc"
MERGE="/Dedicated/jmichaelson-wdata/lcasten/sli_wgs/prs/1000_genomes_EUR_merge/merged.1000_genomes_EUR.SLI_WGS.qc"


##
$PLINK --bfile $KG \
  --keep-allele-order \
  --extract $KP \
  --biallelic-only strict \
  --snps-only just-acgt \
  --make-bed --out $MERGEKG

##
$PLINK --bfile $SLI \
  --keep-allele-order \
  --extract $KP \
  --biallelic-only strict \
  --snps-only just-acgt \
  --make-bed --out $MERGESLI

##
grep -F -f $MERGEKG.bim $MERGESLI.bim > $KP2
wc -l $KP2

## will fail because of multiallelics
$PLINK --bfile $MERGEKG \
  --bmerge $MERGESLI \
  --keep-allele-order \
  --extract $KP2 \
  --make-bed --out $MERGE

## drop problematic variants
$PLINK --bfile $MERGESLI \
  --keep-allele-order \
  --exclude /Dedicated/jmichaelson-wdata/lcasten/sli_wgs/prs/1000_genomes_EUR_merge/merged.1000_genomes_EUR.SLI_WGS.qc-merge.missnp \
  --make-bed --out ${MERGESLI}_tmp
$PLINK --bfile $MERGEKG \
  --keep-allele-order \
  --exclude /Dedicated/jmichaelson-wdata/lcasten/sli_wgs/prs/1000_genomes_EUR_merge/merged.1000_genomes_EUR.SLI_WGS.qc-merge.missnp \
  --make-bed --out ${MERGEKG}_tmp
## re-merge (should work this time)
$PLINK --bfile ${MERGEKG}_tmp \
  --bmerge ${MERGESLI}_tmp \
  --keep-allele-order \
  --extract $KP2 \
  --make-bed --out $MERGE

## get dups
$PLINK --bfile $MERGE \
  --keep-allele-order \
  --list-duplicate-vars \
  --out $MERGE.duplicated_vars
wc -l /Dedicated/jmichaelson-wdata/lcasten/sli_wgs/prs/1000_genomes_EUR_merge/merged.1000_genomes_EUR.SLI_WGS.qc.duplicated_vars.dupvar