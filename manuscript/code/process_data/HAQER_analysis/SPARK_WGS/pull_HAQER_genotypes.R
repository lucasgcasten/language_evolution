## <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
## pulling HAQER variants from SPARK WGS GDS files
## <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

## technical overview of GDS format: https://si.biostat.washington.edu/sites/default/files/modules/GDS_intro.pdf
## markdown tutorial of what the GDS data looks like in an R session: https://uw-gac.github.io/SISG_2021/gds-format.html

args = commandArgs(trailingOnly=TRUE)
chr = as.numeric(args[1])
coi <- chr

## packages / dependencies
library(gdsfmt) ## BioConductor
library(SeqArray) ## BioConductor
library(SeqVarTools) ## BioConductor
library(tidyverse)

haqer <- read_tsv('/Dedicated/jmichaelson-wdata/lcasten/tools/ref_data/human_accelerated_regions/HAQERs/Mangan-Cell2022/HAQER.hg38.bed', col_names = FALSE)

## make table with all ancestral alleles
nean <- read_delim('/Dedicated/jmichaelson-wdata/lcasten/tools/ref_data/archaic_hominini/hg19_neanderthal_filtered_vcf/vcf/merged_all_chromosomes.HAQERs_10Kb_flank.alleles.txt', col_names = FALSE, na = c('NA', '', '.', 'NULL', 'null'), delim = ' ') %>% 
    filter(X5 > 0) %>% 
    mutate(nean_allele = case_when(is.na(X4) ~ X3,
                                   is.na(X4) == FALSE & X6 <= 0.5 ~ X3,
                                   is.na(X4) == FALSE & X6 > 0.5 ~ X4)) %>% 
    select(chromosome = X1, pos = X2, nean_allele) %>% 
    distinct(chromosome, pos, .keep_all = TRUE) %>% 
    drop_na() %>% 
    rename(pos_hg19 = pos) %>% 
    distinct(chromosome, pos_hg19, .keep_all = TRUE)
hca <- read_table('/Dedicated/jmichaelson-wdata/lcasten/tools/ref_data/human_accelerated_regions/HAQERs/Mangan-Cell2022/ancestral_alleles.hg38.txt') %>% 
    mutate(HCA_allele = str_remove_all(INFO, pattern = 'AA=')) %>% 
    mutate(chromosome = as.numeric(str_remove_all(`#CHROM`, 'chr'))) %>% 
    select(chromosome, pos = POS, HCA_allele) %>% 
    distinct(chromosome, pos, .keep_all = TRUE)

table(hca$chromosome)

primates <- read_csv('/Dedicated/jmichaelson-wdata/lcasten/tools/ref_data/human_accelerated_regions/HAQERs/Mangan-Cell2022/5wayPrimate/HAQER_alleles/Mangan_MSA_alleles.csv')
primates <- primates %>% 
    mutate(chromosome = as.numeric(str_remove_all(chr, 'chr'))) %>% 
    select(-chr) %>% 
    rename(pos = pos_hg38) %>% 
    distinct(chromosome, pos, .keep_all = TRUE)

anc_alleles <- primates %>% 
    left_join(nean) %>% 
    left_join(hca) %>% 
    relocate(chromosome, position = pos) %>% 
    relocate(hg38_ref_allele = human_allele, neanderthal_allele = nean_allele, HCA_allele, chimp_allele, bonobo_allele, gorilla_allele, orangutan_allele = orangatuan_allele, .after = pos_hg19)

anc_alleles <- anc_alleles %>% 
    filter(chromosome == coi)
rm(nean, hca, primates)
gc()

####################
cat('\n\n\n\n\n')
message('=============================================================')
message('Working on pulling variant annotations from chromosome ', chr)
message('=============================================================')

###
haqer_chr <- haqer %>% 
    mutate(chr = as.numeric(str_remove_all(X1, 'chr'))) %>% 
    filter(chr == coi) %>% 
    arrange(X2, X3) %>% 
    mutate(X2 = X2 - 10000, 
            X3 = X3 + 10000) ## 10Kb flank around HAQERs


## open GDS of interest
gf = paste0('/Dedicated/jmichaelson-wdata/lcasten/spark/annotation/wgs_all/annotated_geno/gds/WGS_chr', chr, '.gds')
genofile = seqOpen(gf) ## open connection to GDS file
# print(genofile) ## look at available annotations

##
# message("A more detailed overview at all of the info contained in the GDS (VCF header + some info):")
# print(seqSummary(genofile))

## grab variants + basic annotations
start_time = Sys.time()
var_id <- seqGetData(genofile, 'variant.id') ## internal variant identifier for the GDS file
var_type <- seqGetData(genofile, "annotation/info/FunctionalAnnotation/genecode_comprehensive_exonic_category") ## type of variant (synonymous, missense, etc.)
var_info <- seqGetData(genofile, "annotation/info/FunctionalAnnotation/genecode_comprehensive_exonic_info") ## gene/transcript info (what is the variant impacting?)
pos <- seqGetData(genofile, "position") ## hg38 variant BP position
ref <- seqGetData(genofile, "$ref") ## ref allele
alt <- seqGetData(genofile, "$alt") ## alt allele
rsid <- seqGetData(genofile, "annotation/info/FunctionalAnnotation/rsid") ## RSID
cadd <- seqGetData(genofile, "annotation/info/FunctionalAnnotation/cadd_phred") ## CADD score
metasvm = seqGetData(genofile, "annotation/info/FunctionalAnnotation/metasvm_pred") ## MetaSVM (classifier for missense variant deleteriousness)
linsight <- seqGetData(genofile, "annotation/info/FunctionalAnnotation/linsight") ## linsight score
fathmm <- seqGetData(genofile, "annotation/info/FunctionalAnnotation/fathmm_xf") ## fathmm_xf
af <- seqGetData(genofile, "annotation/info/AF") ## allele frequency
samples <- seqGetData(genofile, "sample.id") ## SPARK ID's

#########################
## make table of all variants annotations we pulled to quickly identify variants of interest so we can filter to them + pull genotypes
df <- data.frame(var_id = var_id,
        chromosome = rep(chr, times = length(var_id)),
        position = pos,
        ref = ref,
        alt = alt,
        var_info = var_info,
        var_type = var_type,
        rsid = rsid,
        af = af,
        cadd = cadd,
        linsight = linsight,
        fathmm = fathmm,
        metasvm = metasvm
        ) %>% 
    as_tibble()
end_time = Sys.time()

###################################
## subset to variants of interest
###################################
haqer_chr <- haqer_chr %>% 
    select(chromosome = chr, start = X2, end = X3, haqer_id = X4)

haqer_vars_list <- list()
for(i in 1:nrow(haqer_chr)) {
    haqer_vars_list[[i]] <- subset(df, (position >= haqer_chr$start[i]) & (position <= haqer_chr$end[i])) %>% 
        filter(ref != '*' & alt != '*') %>% 
        mutate(haqer_id = haqer_chr$haqer_id[i],
               haqer_start = haqer_chr$start[i],
               haqer_end = haqer_chr$end[i])
}
haqer_vars <- bind_rows(haqer_vars_list) %>% 
    filter(af > 0 & af < 1) %>% 
    distinct(var_id, .keep_all = TRUE)

kp_vars <- haqer_vars$var_id
haqer_filter <- df$var_id %in% kp_vars ## find variants annotated to HAQERs
message('There are ', sum(haqer_filter), ' variants within 10Kb of HAQERs on this chromosome')
df$haqer <- haqer_filter ## add binary filter back to table for posterity (can filter table to these variants to look at CADD scores or allele frequency)
df$geno_id <- paste0(df$chromosome, ':', df$position, ':', df$ref, ':', df$alt) ## make unique variant identifier
seqSetFilter(genofile, verbose = TRUE, variant.sel = haqer_filter) ## apply filter so we only read in what we need (THIS DOES NOT ACTUALLY MANIPULATE THE GDS FILE, JUST SELECTS SPECIFIC VARIANTS/SAMPLES IN THIS R SESSION)
df_haq <- df %>% 
    filter(haqer == TRUE)

######################################
## get ref pop allele frequencies
######################################
message('Extracting reference population allele frequencies from the VCF annotations now...')
tst <- seqSummary(genofile)
tst<- tst$info %>%
    as_tibble() %>%
    filter(ID == 'CSQ')
vep_str <- tst$Description
gnomadAF <- data.frame(split = unlist(strsplit(vep_str, "[|]"))) %>%
    as_tibble() %>% 
    # tidytext::unnest_tokens(input = vep_str, output = split, to_lower = FALSE, token = 'regex', pattern = '[|]') %>%
    mutate(n = seq(1:nrow(.))) %>%
    filter(split %in% c('AF', 'EUR_AF', 'gnomADe_AF', 'gnomADg_AF', 'gnomADg_NFE_AF', 'MAX_AF', 'MAX_AF_POPS'))
kg <- gnomadAF$n[gnomadAF$split == 'AF']
kg_eur <- gnomadAF$n[gnomadAF$split == 'EUR_AF']
gnme <- gnomadAF$n[gnomadAF$split == 'gnomADe_AF']
gnmg <- gnomadAF$n[gnomadAF$split == 'gnomADg_AF']
gnmg_eur <- gnomadAF$n[gnomadAF$split == 'gnomADg_NFE_AF']
max_af <- gnomadAF$n[gnomadAF$split == 'MAX_AF']
max_af_pop <- gnomadAF$n[gnomadAF$split == 'MAX_AF_POPS']
csq <- seqGetData(genofile, "annotation/info/CSQ", .tolist = FALSE)
tmp_df <- df_haq

df_af <- data.frame(geno_id = rep(str_c(tmp_df$chromosome, ':', tmp_df$position, ':', tmp_df$ref, ':', tmp_df$alt), times = unlist(csq$length)),
        csq = unlist(csq$data))  %>%
        as_tibble() %>%
        distinct(geno_id, .keep_all = TRUE) %>%
        mutate(gnomADe = str_split(csq, pattern = '[|]', simplify = TRUE)[,gnme],
            gnomADg = str_split(csq, pattern = '[|]', simplify = TRUE)[,gnmg],
            gnomADg_NFE_AF = str_split(csq, pattern = '[|]', simplify = TRUE)[,gnmg_eur],
            thousand_genomes = str_split(csq, pattern = '[|]', simplify = TRUE)[,kg],
            thousand_genomes_eur = str_split(csq, pattern = '[|]', simplify = TRUE)[,kg_eur],
            max_gnomAD_af = str_split(csq, pattern = '[|]', simplify = TRUE)[,max_af],
            max_gnomAD_af_pop = str_split(csq, pattern = '[|]', simplify = TRUE)[,max_af_pop]) %>% 
        mutate(gnomADe = as.numeric(gnomADe),
            gnomADg = as.numeric(gnomADg),
            gnomADg_NFE_AF = as.numeric(gnomADg_NFE_AF),
            thousand_genomes = as.numeric(thousand_genomes),
            thousand_genomes_eur = as.numeric(thousand_genomes_eur),
            max_gnomAD_af = as.numeric(max_gnomAD_af)) %>% 
        select(-csq)

rm(csq)
gc()

###################################
## get genotypes
###################################
start_geno <- Sys.time()
gd = seqGetData(genofile, var.name = '$dosage_alt') ## grab genotypes in sample x variant matrix (0 = homozygous reference allele, 1 = heterozygous, 2 = homozygous alternate allele)
seqClose(genofile) ## close connection to GDS file (not necessary, but it's a good habit)
end_geno <- Sys.time()

## add sample names + variant IDs to genotypes and look at them
row.names(gd) <- samples
colnames(gd) <- df$geno_id[df$haqer == TRUE]
gd <- gd[, !duplicated(colnames(gd))]
gd_var_id <- colnames(gd)
message('The HAQER genotype table looks like this:')
print(gd[1:5,1:3])

## print summary stats about HAQER variants
message('There are ', length(var_id), ' total variants on chromosome ', chr)
message('There are ', sum(haqer_filter), ' HAQER variants on this chromosome (max MAF = ', round(max(df$af[df$haqer == TRUE] * 100), digits = 3), '%)')
message('It took ', round(difftime(end_time, start_time), digits = 2), ' seconds to pull the variant info into a dataframe')
message('Genotype table dimensions: samples (rows) = ', nrow(gd), ', variants (cols) = ', ncol(gd))
message('It took ', round(difftime(end_geno, start_geno), digits = 2), ' seconds to read in all of the genotypes for variants of interest')

gd[is.na(gd)] <- 0
# gd[1:5,1:5]
gd_t <- data.table::transpose(as.data.frame(gd))
message('Here, just finished transposing the genotype DF')

rm(gd)
gc()

# dim(gd_t)
# gd_t[1:5,1:5]
colnames(gd_t) <- samples
row.names(gd_t) <- gd_var_id
# length(gd_var_id)
# dim(gd_t)
gd_t <- gd_t %>% 
    rownames_to_column(var = 'geno_id')  %>% 
    as_tibble()

message('Here, about to merge all data together')
wd <- df_haq %>% 
    inner_join(df_af) %>%
    left_join(anc_alleles) %>%
    inner_join(gd_t) %>% 
    filter(af > 0 & af < 1) %>%
    distinct(chromosome, position, ref, alt, .keep_all = TRUE) %>% 
    distinct(geno_id, .keep_all = TRUE)

message('Final dataset has: ', nrow(wd), ' variants')
## remove tmp data
# rm(wd, gd_t, gd, df_af, df_haq, csq, tmp_df, df, haqer_chr)

## save to file
outf <- str_c('/Dedicated/jmichaelson-wdata/lcasten/sli_wgs/HAQER_variant_associations/data/SPARK_WGS/chr', chr, '_HAQER_10Kb_flank.hg38.rds')
wd %>% 
    write_rds(outf)

## done, can now use table for phenotype correlations or save it for future use