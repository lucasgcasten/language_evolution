## ------------------------------------------------------------------------------------------------
## This script will compute LDPred2 PGS
## ------------------------------------------------------------------------------------------------

#### set up ####
library(tidyverse)
library(runonce)
library(remotes)
library(bigsnpr)
options(bigstatsr.check.parallel.blas = FALSE)
options(default.nproc.blas = NULL)
library(data.table)
library(magrittr)

#### get arguments ####
args = commandArgs(trailingOnly = TRUE)

bed_filepath = args[1]
core_count = as.numeric(args[2])

## set output directory ##
my_wd <- c('/Dedicated/jmichaelson-wdata/lcasten/sli_wgs/prs/LDPred2-inf-v2/')
dir.create(my_wd, recursive = TRUE)
setwd(my_wd)

##
message('##########################')
message('Input PLINK files: ')
message(bed_filepath)
message('Directory where everything will be output is here:')
message(my_wd)
message('Using ', core_count, ' CPU cores for computations...')

#### get LDPred2 function ####
source('/Dedicated/jmichaelson-wdata/lcasten/functions/LDPred2_argon_hm3plus.R')

map_ld = read_rds('/Dedicated/jmichaelson-wdata/lcasten/tools/LDPred2/map_hm3_plus.rds')
map_ldref = read_rds('/Dedicated/jmichaelson-wdata/lcasten/tools/LDPred2/map_hm3_plus.rds')

## subset to overlapping SNPs
bim_snp <- read_tsv(file = str_replace_all(bed_filepath, pattern = '[.]bed', replacement = '.bim'), col_names = F) %>%
  rename(chr = X1, pos = X4) %>%
  select(chr, pos, major = X6, minor = X5) %>%
  filter(nchar(major) == 1 & nchar(minor) == 1)

map_ld <- map_ld %>%
  inner_join(bim_snp) %>%
  filter(major == a0 | major == a1) %>%
  filter(minor == a0 | minor == a1) %>%
  select(-c(major, minor))

map_ldref <- map_ldref %>%
  inner_join(bim_snp) %>%
  filter(major == a0 | major == a1) %>%
  filter(minor == a0 | minor == a1) %>%
  select(-c(major, minor))

# message(str_c('Before dropping duplicate positions/overlapping with HM3+ there were ', nrow(bim_snp), ' variants in the BIM file'))
# bim_snp <- bim_snp %>%
#   inner_join(map_ld) %>%
#   filter(major == a0 | major == a1) %>%
#   filter(minor == a0 | minor == a1) %>%
#   distinct(chr, pos, .keep_all = T)
# message(str_c('After dropping problematic SNPs there are ', nrow(bim_snp), ' variants in the BIM file'))


##
######### STEP 0: get geno data ready
# cat("\n\n\n\n")
# message("===========================================")
# message('Converting BED to RDS file now')
# message("===========================================")
# snp_readBed(bed_filepath)
cat('\n\n\n')
message('Finished conversion')
rds_path <- str_replace_all(bed_filepath, pattern = '[.]bed', replacement = '.rds')

#####################################################
## compute PGS for the traits of interest
#####################################################
## NDC gwas from Huang et al. 2024 Nature paper
ph <- "rare_neurodevelopmental_condition"
ncas <- 10015
ncon <- 22937
n = round(4 / (1 / ncas + 1 / ncon), digits = 0) ## effective sample size

## compute PRS
ss <- read_tsv('/Dedicated/jmichaelson-sdata/gwas_summary_stats/neurodevelopmental_conditions-Huang-Nature2024/supplementary_data3_GWAS_meta.txt.gz') %>% 
        mutate(N = n)
pheno = ph
output_path <- str_c(my_wd, pheno, "/")
str_c("the final output path is: ", output_path)
dir.create(output_path)
setwd(output_path)
final_sumstats <- ss %>% 
    select(rsid = SNP,
            chr = CHR,
            pos = BP, # check the genome build in either the paper or look at dbSNP (https://ncbi.nlm.nih.gov/snp/)
            a1 = effect_allele, # a1 is the effect allele
            a0 = other_allele,
            n_eff = N,
            beta,
            beta_se = SE,
            p = P) %>%
    drop_na() %>%
    filter(beta_se > 0) %>%
    inner_join(select(map_ldref, rsid, chr, pos)) %>%
    relocate(rsid, chr, pos) %>%
    inner_join(select(bim_snp, chr, pos)) %>% 
    distinct(rsid, chr, pos, .keep_all = TRUE) ## remove any duplicate SNPs in case there are any (likely multiallelic sites)
try(calculate_ldpred_pgs(sumstats = final_sumstats, n_core = core_count,  plink_rds = rds_path, sd_y = 0, pheno_name = pheno, build = 'hg19'))


###############################
## addiction GWAS
ph <- "addiction_Hatoum2023"
n = 1025550 ## sample size reported in paper abstract

## compute PRS
ss <- read_tsv('/Dedicated/jmichaelson-sdata/gwas_summary_stats/PGC/pgc-sub2023/Hatoum2023AddictionEuropean.txt') %>% 
        mutate(N = n)

pheno = ph
output_path <- str_c(my_wd, pheno, "/")
str_c("the final output path is: ", output_path)
dir.create(output_path)
setwd(output_path)
final_sumstats <- ss %>% 
    mutate(se  = abs(beta / qnorm(p / 2))) %>% 
    select(rsid = snp,
            chr = chr,
            pos = bp, # check the genome build in either the paper or look at dbSNP (https://ncbi.nlm.nih.gov/snp/)
            a1 = a1, # a1 is the effect allele
            a0 = a2,
            n_eff = N,
            beta,
            beta_se = se,
            p = p) %>%
    drop_na() %>%
    filter(beta_se > 0) %>%
    inner_join(select(map_ldref, rsid, chr, pos)) %>%
    relocate(rsid, chr, pos) %>%
    inner_join(select(bim_snp, chr, pos)) %>% 
    distinct(rsid, chr, pos, .keep_all = TRUE) ## remove any duplicate SNPs in case there are any (likely multiallelic sites)
try(calculate_ldpred_pgs(sumstats = final_sumstats, n_core = core_count,  plink_rds = rds_path, sd_y = 0, pheno_name = pheno, build = 'hg19'))



###########################################
## epilepsy GWAS'
## GWAS catalog accessions: GCST90271608-GCST90271620
files = list.files('/sdata/gwas_summary_stats/epilepsy/ILAE-NatureGenetics2023/final_sumstats', pattern = '.tbl$')
files = files[str_detect(files, pattern = 'Cauc', negate = TRUE)]
files

files <- data.frame(file = c('ILAE3_CAE_final.tbl', 'ILAE3_focal_HS_final.tbl', 'ILAE3_focal_lesion_negative_final.tbl', 'ILAE3_focal_other_lesion_final.tbl', 'ILAE3_GTCS_final.tbl', 'ILAE3_JAE_final.tbl', 'ILAE3_JME_final.tbl', 'ILAE3_TRANS_all_epilepsy_final.tbl', 'ILAE3_TRANS_focal_epilepsy_final.tbl', 'ILAE3_TRANS_GGE_final.tbl'),
           pheno = c('epilepsy_subtype_childhood_absence', 'epilepsy_subtype_hippocampal_sclerosis', 'epilepsy_subtype_focal_lesion_negative', 'epilepsy_subtype_other_lesion', 'epilepsy_subtype_generalized_tonic_clonic_seizures', 'epilepsy_subtype_juvenile_abscence', 'epilepsy_subtype_juvenile_myoclonic', 'epilepsy', 'epilepsy_subtype_focal_epilepsy', 'epilepsy_subtype_genetic_generalized'))

for(i in 1:nrow(files)) {
    ## 
    cat('\n\n\n\n')
    message('==================================================================')
    message('Working on: ', files$pheno[i], ' (', i, ' / ', nrow(files), ')')
    message('==================================================================')
    
    ##
    pheno <- files$pheno[i]

    ##
    ss <- read_table(str_c('/Dedicated/jmichaelson-sdata/gwas_summary_stats/epilepsy/ILAE-NatureGenetics2023/final_sumstats/', files$file[i]))

    ##
    output_path <- str_c(my_wd, pheno, "/")
    str_c("the final output path is: ", output_path)
    dir.create(output_path)
    setwd(output_path)

    ##
    final_sumstats <- ss %>% 
        filter(Freq1 >= 0.01) %>%
        select(rsid = MarkerName,
                chr = CHR,
                pos = BP, # check the genome build in either the paper or look at dbSNP (https://ncbi.nlm.nih.gov/snp/)
                a1 = Allele1, # a1 is the effect allele
                a0 = Allele2,
                n_eff = Effective_N,
                beta = Beta,
                beta_se = SE,
                p = `P-value`) %>%
        drop_na() %>%
        filter(beta_se > 0) %>%
        mutate(a1 = toupper(a1),
               a0 = toupper(a0)) %>%
        arrange(chr, pos) %>%
        filter(is.finite(beta) == TRUE & is.finite(beta_se) == TRUE) %>%
        inner_join(select(map_ldref, rsid, chr, pos)) %>%
        relocate(rsid, chr, pos) %>%
        inner_join(select(bim_snp, chr, pos)) %>% 
        distinct(rsid, chr, pos, .keep_all = TRUE) ## remove any duplicate SNPs in case there are any (likely multiallelic sites)

    try(calculate_ldpred_pgs(sumstats = final_sumstats, n_core = core_count,  plink_rds = rds_path, sd_y = 0, pheno_name = pheno, build = 'hg19'))
}