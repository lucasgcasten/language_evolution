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
my_wd <- c('/Dedicated/jmichaelson-wdata/lcasten/sli_wgs/prs/LDPred2-inf-test/')
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
source('/Dedicated/jmichaelson-wdata/lcasten/functions/LDPred2_argon_hm3plus_v2.R')

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
# ## -------------------------------------------
# #### compute cog perf pgs ####
# ## -------------------------------------------
ss <- read_tsv('/Dedicated/jmichaelson-wdata/trthomas/abcd_spark_merge/polygenic_scores/sumstats/Social_Science_Genetic_Association_Consortium/GWAS_CP_all.txt')
pheno = '2018_ssgac_cognitive_performance'
#
##
str_c("number of rows in the raw summary statistics file is: ", nrow(ss))
#
ss$N <- 257828
#
### 3: give the output path
output_path <- str_c(my_wd, pheno, "/")
str_c("the final output path is: ", output_path)
system(str_c('rm -r ', output_path))
dir.create(output_path)
setwd(output_path)
#
##### 4: modify/create/filter the summary statistics columns to have: chr, pos, a1 (effect allele), a0, beta, beta_se, p, n_eff
final_sumstats <- ss %>%
  select(chr = CHR,
         pos = POS,
         a1 = A1, # a1 is the effect allele
         a0 = A2,
         n_eff = N,
         beta = Beta,
         beta_se = SE,
         p = Pval) %>%
  filter(a1 != a0) %>%
  drop_na() %>%
  inner_join(select(bim_snp, chr, pos)) %>%
  inner_join(select(map_ldref, chr, pos))
#
#
##
try(calculate_ldpred_pgs(sumstats = final_sumstats, n_core = core_count,  plink_rds = rds_path, sd_y = 1, pheno_name = pheno, build = 'hg19'))
# #
# # ## -------------------------------------------
# # #### compute ea pgs ####
# # ## -------------------------------------------
# # # n = 766345
system('cat /Dedicated/jmichaelson-sdata/gwas_summary_stats/ssgac/2022/ReadMe_EA4.txt')
ss <- read_tsv('/Dedicated/jmichaelson-sdata/gwas_summary_stats/ssgac/2022/EA4_additive_excl_23andMe.txt.gz') %>%
  select(-matches('_unadj$'))
pheno = '2022_ssgac_educational_attainment'
#
#
str_c("number of rows in the raw summary statistics file is: ", nrow(ss))
head(ss)
#
ss$N <- 765283
#
### 3: give the output path
output_path <- str_c(my_wd, pheno, "/")
str_c("the final output path is: ", output_path)
system(str_c('rm -r ', output_path))
dir.create(output_path)
setwd(output_path)
#
##### 4: modify/create/filter the summary statistics columns to have: chr, pos, a1 (effect allele), a0, beta, beta_se, p, n_eff
final_sumstats <- ss %>%
  select(chr = Chr,
         pos = BP,
         a1 = Effect_allele, # a1 is the effect allele
         a0 = Other_allele,
         n_eff = N,
         beta = Beta,
         beta_se = SE,
         p = P) %>%
  filter(a1 != a0) %>%
  drop_na() %>%
  inner_join(select(bim_snp, chr, pos)) %>%
  inner_join(select(map_ldref, chr, pos))
#
#
##
try(calculate_ldpred_pgs(sumstats = final_sumstats,  n_core = core_count, plink_rds = rds_path, sd_y = 1, pheno_name = pheno, build = 'hg19'))
# #
# #
# # ## -------------------------------------------
# # #### compute scz pgs ####
# # ## -------------------------------------------
ss <- read_tsv("/Dedicated/jmichaelson-wdata/trthomas/abcd_spark_merge/polygenic_scores/sumstats/Psychiatric_Genomics_Consortium/schizophrenia-PGC-2021")
pheno = '2021_pgc_schizophrenia'
#
#
##
str_c("number of rows in the raw summary statistics file is: ", nrow(ss))
head(ss)
#
cases = 69369
controls = 236642
ss$N <- 4 / (1 / cases + 1 / controls)
#
### 3: give the output path
output_path <- str_c(my_wd, pheno, "/")
str_c("the final output path is: ", output_path)
system(str_c('rm -r ', output_path))
dir.create(output_path)
setwd(output_path)
#
##### 4: modify/create/filter the summary statistics columns to have: chr, pos, a1 (effect allele), a0, beta, beta_se, p, n_eff
final_sumstats <- ss %>%
  mutate(beta = log(OR),
         se_beta = abs(beta / qnorm(P / 2))
  )
#
print(
  str_c('head of OR SE',
        head(final_sumstats$SE))
)
print(
  str_c('head of beta SE',
        head(final_sumstats$se_beta))
)
#
final_sumstats = final_sumstats %>%
  select(chr = CHR,
         pos = BP,
         a1 = A1, # a1 is the effect allele
         a0 = A2,
         n_eff = N,
         beta,
         beta_se = se_beta,
         p = P) %>%
  filter(a1 != a0) %>%
  drop_na() %>%
  inner_join(select(bim_snp, chr, pos)) %>%
  inner_join(select(map_ldref, chr, pos))
#
#
##
try(calculate_ldpred_pgs(sumstats = final_sumstats,  n_core = core_count, plink_rds = rds_path, sd_y = 2, pheno_name = pheno, build = 'hg19'))
# #
