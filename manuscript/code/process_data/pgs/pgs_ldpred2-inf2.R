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

message(str_c('Before dropping duplicate positions/overlapping with HM3+ there were ', nrow(bim_snp), ' variants in the BIM file'))
bim_snp <- bim_snp %>%
  inner_join(map_ld) %>%
  filter(major == a0 | major == a1) %>%
  filter(minor == a0 | minor == a1) %>%
  distinct(chr, pos, .keep_all = T)
message(str_c('After dropping problematic SNPs there are ', nrow(bim_snp), ' variants in the BIM file'))


##
######### STEP 0: get geno data ready
# snp_readBed(bed_filepath)
rds_path <- str_replace_all(bed_filepath, pattern = '[.]bed', replacement = '.rds')


##########################################
################################ 
## brain connectivity
sumstats_files = list.files(path = '/Dedicated/jmichaelson-sdata/gwas_summary_stats/brain_connectivity/Tissink_eNeuro_2023', 
                            pattern = '[.]txt[.]gz', full.names = T, recursive = T)

## file paths
ss_map = data.frame(filepath = sumstats_files, 
                            pheno = basename(sumstats_files)
                            ) %>%
  mutate(pheno = str_remove_all(pheno, pattern = 'Tissinketal_eNeuro2023_|_EUR_sumstats[.]txt[.]gz'),
         pheno = str_replace_all(string = pheno, pattern = '^SC_', replacement = 'structuralConnectivity_Yeo7networkParcellation_'),
         pheno = str_replace_all(string = pheno, pattern = '^FC_', replacement = 'functionalConnectivity_Yeo7networkParcellation_'),
         )

#### calculate pgs ####
for (i in 1:nrow(ss_map)) {
  message('')
  message('')
  message('######################################')
  message('Working on pheno ', i, '/' , nrow(ss_map))
  
  pheno = str_c(ss_map$pheno[i])
  output_path <- str_c(my_wd, pheno, "/")
  str_c("the final output path is: ", output_path)
  system(str_c('rm -r ', output_path))
  
  ## check if PGS folder already exists, skip if it does
  if (pheno %in% list.dirs(path = str_sub(my_wd, start = 1, end = -2), recursive = F, full.names = F)) {
    message('PGS already computed, moving to next one...')
  }
  
  ##
  if (! pheno %in% list.dirs(path = str_sub(my_wd, start = 1, end = -2), recursive = F, full.names = F)) {
    message('PGS has not been computed')
    message('Calculating PGS for: ', pheno)
    ss <- read_tsv(as.character(ss_map$filepath[i]))
    
    str_c("number of rows in the raw summary statistics file is: ", nrow(ss))
    head(ss)
    
    ### 3: give the output path
    output_path <- str_c(my_wd, pheno, "/")
    str_c("the final output path is: ", output_path)
    system(str_c('rm -r ', output_path))
    dir.create(output_path)
    setwd(output_path)
    
    ##### 4: modify/create/filter the summary statistics columns to have: chr, pos, a1 (effect allele), a0, beta, beta_se, p, n_eff
    message('before GWAS sumstats QC there were ', nrow(ss), ' variants')
    final_sumstats <- ss %>%
      mutate(chromosome = as.numeric(chromosome)) %>%
      drop_na(chromosome) %>%
      filter(effect_allele_frequency >= .01) %>%
      filter(nchar(effect_allele) == 1 & nchar(other_allele) == 1)
    
    ##
    final_sumstats$Neff = final_sumstats$n
    
    message('after GWAS sumstats QC and before overlapping w/ HapMap3 there are ', nrow(final_sumstats), ' SNPs')
    message('(dropped ', nrow(ss) - nrow(final_sumstats), ' variants)')
    
    final_sumstats = final_sumstats %>% 
      select(chr = chromosome,
             pos = base_pair_location,
             a1 = effect_allele, # a1 is the effect allele
             a0 = other_allele,
             n_eff = Neff,
             beta = beta,
             beta_se = standard_error,
             p = p_value) %>%
    drop_na() %>%
    inner_join(select(bim_snp, chr, pos)) %>%
    inner_join(select(map_ldref, chr, pos))
    
    ## compute PGS
    try(calculate_ldpred_pgs(sumstats = final_sumstats, n_core = core_count,  plink_rds = rds_path, sd_y = 0, pheno_name = pheno, build = 'hg19'))
  }
}

#### household income ####
# ss <- read_table('/Dedicated/jmichaelson-sdata/gwas_summary_stats/UKBIOBANK/Hill2019_income/DS_10283_3441/Household_Income_UKBiobank.txt.gz')
# pheno = '2019_Hill_householdIncome'
# ss = ss %>%
#   drop_na() %>%
#   filter(nchar(Effect_Allele) == 1 & nchar(Non_effect_Allele) == 1)
# #
# # ##
# # str_c("number of rows in the raw summary statistics file is: ", nrow(ss))
# # head(ss)
# #
# ### 3: give the output path
# output_path <- str_c(my_wd, pheno, "/")
# str_c("the final output path is: ", output_path)
# system(str_c('rm -r ', output_path))
# dir.create(output_path)
# setwd(output_path)
# #
# ##### 4: modify/create/filter the summary statistics columns to have: chr, pos, a1 (effect allele), a0, beta, beta_se, p, n_eff
# final_sumstats <- ss %>%
#   mutate(N = 286301,
#          se_beta = Standard_Error_of_Beta) %>%
#   select(
#     chr = Chr,
#     pos = BPos,
#     a1 = Effect_Allele, # a1 is the effect allele
#     a0 = Non_effect_Allele,
#     n_eff = N,
#     beta = Beta,
#     beta_se = se_beta,
#     p = P) %>%
#   filter(a1 != a0) %>%
#   drop_na() %>%
#   inner_join(select(bim_snp, chr, pos)) %>%
#   inner_join(select(map_ldref, chr, pos))
# #
# #
# #
# try(calculate_ldpred_pgs(sumstats = final_sumstats,  n_core = core_count, plink_rds = rds_path, sd_y = 0, pheno_name = pheno, build = 'hg19'))
# ## done