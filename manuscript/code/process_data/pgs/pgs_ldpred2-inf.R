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
ss <- read_table('/Dedicated/jmichaelson-sdata/gwas_summary_stats/UKBIOBANK/Hill2019_income/DS_10283_3441/Household_Income_UKBiobank.txt.gz')
pheno = '2019_Hill_householdIncome'
ss = ss %>%
  drop_na() %>%
  filter(nchar(Effect_Allele) == 1 & nchar(Non_effect_Allele) == 1)
#
# ##
# str_c("number of rows in the raw summary statistics file is: ", nrow(ss))
# head(ss)
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
  mutate(N = 286301,
         se_beta = Standard_Error_of_Beta) %>%
  select(
    chr = Chr,
    pos = BPos,
    a1 = Effect_Allele, # a1 is the effect allele
    a0 = Non_effect_Allele,
    n_eff = N,
    beta = Beta,
    beta_se = se_beta,
    p = P) %>%
  filter(a1 != a0) %>%
  drop_na() %>%
  inner_join(select(bim_snp, chr, pos)) %>%
  inner_join(select(map_ldref, chr, pos))
#
#
#
try(calculate_ldpred_pgs(sumstats = final_sumstats,  n_core = core_count, plink_rds = rds_path, sd_y = 1, pheno_name = pheno, build = 'hg19'))
#




## -------------------
#### new measures ####
## -------------------
## loneliness stuff
file_names = list.files('/Dedicated/jmichaelson-sdata/gwas_summary_stats/loneliness_Day_2018_NatComs/Upload', pattern = '_Imputed.txt.gz', full.names = T)
file_names = as.character(file_names)
files = data.frame(file = file_names,
                   file_base = basename(file_names),
                   ncas = c(80134, 124047, 66259, 135060),
                   ncon = c(364890, 328255, 386043, 317242)
                   ) %>%
  mutate(file = as.character(file)) %>%
  mutate(pheno = case_when(str_detect(string = file_base, pattern = 'Loneliness_') == T ~ 'Loneliness',
                           str_detect(string = file_base, pattern = 'Pub_') == T ~'Regular_Pub_Attendance',
                           str_detect(string = file_base, pattern = 'Rel_EurRel') == T ~ 'Regular_Religious_Service_Attendance',
                           str_detect(string = file_base, pattern = 'Sport') == T ~ 'Regular_Gym_Attendance')
                   ) %>%
  mutate(N = 4 / (1 / ncas + 1 / ncon))
files


message('##########################')
message('Starting weights now...')

for (i in 1:nrow(files)) {
  ss <- read_table(files$file[i])

  ss$N = files$N[i]

  pheno = files$pheno[i]

  ##
  str_c("number of rows in the raw summary statistics file is: ", nrow(ss))

  ### 3: give the output path
  output_path <- str_c(my_wd, pheno, "/")
  str_c("the final output path is: ", output_path)
  system(str_c('rm -r ', output_path))
  dir.create(output_path)
  setwd(output_path)
#
  ##### 4: modify/create/filter the summary statistics columns to have: chr, pos, a1 (effect allele), a0, beta, beta_se, p, n_eff
  final_sumstats <- ss %>%
    filter(INFO >= 0.8) %>%
    filter(A1FREQ >= 0.01) %>%
    select(
      chr = CHR,
      pos = BP,
      a1 = ALLELE1, # a1 is the effect allele
      a0 = ALLELE0,
      n_eff = N,
      beta = BETA,
      beta_se = SE,
      p = P_BOLT_LMM) %>%
    filter(a1 != a0) %>%
    drop_na() %>%
    inner_join(select(bim_snp, chr, pos)) %>%
    inner_join(select(map_ldref, chr, pos))
#
  ##
  try(calculate_ldpred_pgs(sumstats = final_sumstats,  n_core = core_count, plink_rds = rds_path, sd_y = 0, pheno_name = pheno, build = 'hg19'))
}
#
#


sumstats = '/Dedicated/jmichaelson-sdata/gwas_summary_stats/deCODE/vocal_gisladottir_SciAdvances2023/reformatted_hg19_HapMap3_deCode_VoiceSpeech_MedianF0_MeasureReading.csv'
ss = read_csv(sumstats)

#### calculate pgs ####
pheno = "vocal_pitch_median_f0"

str_c("number of rows in the raw summary statistics file is: ", nrow(ss))
head(ss)

### 3: give the output path
output_path <- str_c(my_wd, pheno, "/")
str_c("the final output path is: ", output_path)
dir.create(output_path)
setwd(output_path)

##### 4: modify/create/filter the summary statistics columns to have: chr, pos, a1 (effect allele), a0, beta, beta_se, p, n_eff
final_sumstats <- ss
##
final_sumstats$Neff = 12901

final_sumstats = final_sumstats %>% 
  select(chr,
         pos,
         a1, # a1 is the effect allele
         a0,
         n_eff = Neff,
         beta,
         beta_se = se,
         p = P)


## compute PGS
try(calculate_ldpred_pgs(sumstats = final_sumstats, n_core = core_count,  plink_rds = rds_path, sd_y = 1, pheno_name = pheno, build = 'hg19'))

# #### -------------------------
## CTG lab brain imaging stuff
file_names = list.files('/Dedicated/jmichaelson-wdata/lcasten/misc/gwas_summary_stats/ctg_sumstats_posthuma/imaging_processed', pattern = 'Tissink', full.names = T)
file_names = as.character(file_names)
file_names
files = data.frame(file = file_names) %>%
  mutate(file = as.character(file)) %>%
  mutate(pheno = case_when(str_detect(basename(file), 'cerebellarvolume') == T ~ 'cerebellar_volume',
                           str_detect(basename(file), 'cerebralvolume') == T ~'cerebral_volume',
                           str_detect(basename(file), 'subcorticalvolume') == T ~ 'subcortical_volume'),
#
  )
files
#
for (i in 1:nrow(files)) {
  ss <- read_table(files$file[i])
#
  pheno = files$pheno[i]
#
  ##
  str_c("number of rows in the raw summary statistics file is: ", nrow(ss))
#
  ### 3: give the output path
  output_path <- str_c(my_wd, pheno, "/")
  str_c("the final output path is: ", output_path)
  system(str_c('rm -r ', output_path))
  dir.create(output_path)
  setwd(output_path)
#
  ##### 4: modify/create/filter the summary statistics columns to have: chr, pos, a1 (effect allele), a0, beta, beta_se, p, n_eff
  max_n = max(ss$N, na.rm = T)
#
  final_sumstats <- ss %>%
    filter(MAF >= .01) %>%
    filter(INFO >= .8) %>%
    filter(N >= .75 * max_n) %>%
    select(
      chr = CHR,
      pos = BP,
      a1 = A1, # a1 is the effect allele
      a0 = A2,
      n_eff = N,
      beta = BETA,
      beta_se = SE,
      p = P) %>%
    drop_na() %>%
    filter(a1 != a0) %>%
    drop_na() %>%
    inner_join(select(bim_snp, chr, pos)) %>%
    inner_join(select(map_ldref, chr, pos))
#
#
  ##
  try(calculate_ldpred_pgs(sumstats = final_sumstats,  n_core = core_count, plink_rds = rds_path, sd_y = 0, pheno_name = pheno, build = 'hg19'))
}
#



#
#### cross disorder
f = '/Dedicated/jmichaelson-wdata/lcasten/misc/gwas_summary_stats/ctg_sumstats_posthuma/All_PsychiatricDisorders_MetaAnalysis.txt.gz'
ss <- read_table(f)
#
pheno = 'CTGlab_2022_Cross_Disorder'
#
##
str_c("number of rows in the raw summary statistics file is: ", nrow(ss))
#
### 3: give the output path
output_path <- str_c(my_wd, pheno, "/")
str_c("the final output path is: ", output_path)
system(str_c('rm -r ', output_path))
dir.create(output_path)
setwd(output_path)
#
##### 4: modify/create/filter the summary statistics columns to have: chr, pos, a1 (effect allele), a0, beta, beta_se, p, n_eff
max_n = max(ss$n_sum, na.rm = T)
#
final_sumstats <- ss %>%
  filter(effect_allele_frequency >= .01) %>%
  filter(n_sum >= .75 * max_n) %>%
  select(
    chr = chromosome,
    pos = base_pair_location,
    a1 = effect_allele, # a1 is the effect allele
    a0 = other_allele,
    n_eff,
    beta = beta,
    beta_se = standard_error,
    p = p_value)%>%
  filter(a1 != a0) %>%
  drop_na() %>%
  inner_join(select(bim_snp, chr, pos)) %>%
  inner_join(select(map_ldref, chr, pos))
#
#
##
try(calculate_ldpred_pgs(sumstats = final_sumstats,  n_core = core_count, plink_rds = rds_path, sd_y = 0, pheno_name = pheno, build = 'hg19'))
#
#
#
#### BMI
f = '/Dedicated/jmichaelson-wdata/lcasten/misc/gwas_summary_stats/physical/2018_giantUKBB_meta_bmi.txt'
ss <- read_table(f)
#
pheno = '2018_giantUKBB_meta_bmi'
#
##
str_c("number of rows in the raw summary statistics file is: ", nrow(ss))
#
### 3: give the output path
output_path <- str_c(my_wd, pheno, "/")
str_c("the final output path is: ", output_path)
system(str_c('rm -r ', output_path))
dir.create(output_path)
setwd(output_path)
#
##### 4: modify/create/filter the summary statistics columns to have: chr, pos, a1 (effect allele), a0, beta, beta_se, p, n_eff
max_n = max(ss$N, na.rm = TRUE)
#
final_sumstats <- ss %>%
  filter(Freq_Tested_Allele >= .01) %>%
  filter(INFO >= .8) %>%
  filter(N > .75 * max_n) %>%
  select(
    chr = CHR,
    pos = POS,
    a1 = Tested_Allele, # a1 is the effect allele
    a0 = Other_Allele,
    n_eff = N,
    beta = BETA,
    beta_se = SE,
    p = P) %>%
  filter(a1 != a0) %>%
  drop_na() %>%
  inner_join(select(bim_snp, chr, pos)) %>%
  inner_join(select(map_ldref, chr, pos))
#
#
##
try(calculate_ldpred_pgs(sumstats = final_sumstats,  n_core = core_count, plink_rds = rds_path, sd_y = 0, pheno_name = pheno, build = 'hg19'))
#
#
# #### type 2 diabetes
f = '/Dedicated/jmichaelson-wdata/lcasten/misc/gwas_summary_stats/t2d/Mahajan.NatGenet2018b.T2Dbmiadj.European.txt'
ss <- read_table(f, col_types = c('cddccddddd'))
#
pheno = 'MahajanNatGenet_2018_T2DbmiAdj'
#
##
str_c("number of rows in the raw summary statistics file is: ", nrow(ss))
#
### 3: give the output path
output_path <- str_c(my_wd, pheno, "/")
str_c("the final output path is: ", output_path)
system(str_c('rm -r ', output_path))
dir.create(output_path)
setwd(output_path)
#
##### 4: modify/create/filter the summary statistics columns to have: chr, pos, a1 (effect allele), a0, beta, beta_se, p, n_eff
max_n = max(ss$Neff, na.rm = T)
#
final_sumstats <- ss %>%
  filter(EAF >= .01) %>%
  #  filter(INFO >= .8) %>%
  filter(Neff > .75 * max_n) %>%
  select(
    chr = Chr,
    pos = Pos,
    a1 = EA, # a1 is the effect allele
    a0 = NEA,
    n_eff = Neff,
    beta = Beta,
    beta_se = SE,
    p = Pvalue) %>%
  filter(a1 != a0) %>%
  drop_na() %>%
  inner_join(select(bim_snp, chr, pos)) %>%
  inner_join(select(map_ldref, chr, pos))
#
#
##
try(calculate_ldpred_pgs(sumstats = final_sumstats,  n_core = core_count, plink_rds = rds_path, sd_y = 0, pheno_name = pheno, build = 'hg19'))
#
#
## ----------------------------------
#### genlang #####
## ----------------
# #
genlang = list.files('/Dedicated/jmichaelson-wdata/lcasten/misc/gwas_summary_stats/genlang/pnas_2022_reading_gwas', full.names = T)
genlang = as.character(genlang)
genlang = genlang[str_detect(genlang, pattern = 'NREP', negate = T)]
# genlang = rev(genlang)
#
for (f in genlang) {
  tmp_pheno = basename(f)
  tmp_pheno = str_remove_all(tmp_pheno, 'METAANALYSIS_')
  tmp_pheno = str_split(tmp_pheno, pattern = '_combined', simplify = T)[,1]
  tmp_pheno = str_remove_all(tmp_pheno, '^RANDOM__|_EUR$')
  tmp_pheno = str_c('2022_genlang_', tmp_pheno)
#
  ss <- read_table(file = f)
#
  ss = ss %>%
    mutate(Allele1 = toupper(Allele1),
           Allele2 = toupper(Allele2)
    ) %>%
    inner_join(map_ld, by = c('MarkerName' = 'rsid')) %>%
    mutate(Effect = case_when(Allele1 == a1 & Allele2 == a0 ~ Effect,
                              Allele1 == a0 & Allele2 == a1 ~ -1 * Effect,
                              TRUE ~ NA_real_)
    ) %>%
    drop_na(Effect) %>%
    select(-c(a0, a1)) %>%
    select(-c(af_UKBB, ld, pos_hg18, pos_hg38))
#
  pheno = tmp_pheno
#
  ##
  str_c("number of rows in the raw summary statistics file is: ", nrow(ss))
  head(ss)
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
    select(
      chr, pos,
      a1 = Allele1, # a1 is the effect allele
      a0 = Allele2,
      n_eff = TotalSampleSize,
      beta = Effect,
      beta_se = StdErr,
      p = `P-value`) %>%
    filter(a1 != a0) %>%
    drop_na() %>%
    inner_join(select(bim_snp, chr, pos)) %>%
    inner_join(select(map_ldref, chr, pos))
#
#
  ##
  try(calculate_ldpred_pgs(sumstats = final_sumstats,  n_core = core_count, plink_rds = rds_path, sd_y = 0, pheno_name = pheno, build = 'hg19'))
}
#
#
## -----------------------
#### SES ####
## -----------------------
#### household income ####
ss <- read_table('/Dedicated/jmichaelson-wdata/trthomas/abcd_spark_merge/polygenic_scores/sumstats/UK_Biobank/Lothian_Birth/SES/Hill2016_UKB_Townsend_summary_results_08112016.txt') %>%
  mutate(Chromosome = as.numeric(Chromosome)) %>%
  drop_na()
#
pheno = '2016_Hill_TownsendIndex'
ss = ss %>%
  filter(nchar(Effect_allele) == 1 & nchar(Other_allele) == 1)
#
##
str_c("number of rows in the raw summary statistics file is: ", nrow(ss))
head(ss)
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
  rename(P = `P-value`) %>%
  mutate(N = 112005,
         se_beta = abs(Beta / qnorm(P / 2))
  ) %>%
  select(
    chr = Chromosome,
    pos = Position,
    a1 = Effect_allele, # a1 is the effect allele
    a0 = Other_allele,
    n_eff = N,
    beta = Beta,
    beta_se = se_beta,
    p = P) %>%
  filter(a1 != a0) %>%
  drop_na() %>%
  inner_join(select(bim_snp, chr, pos)) %>%
  inner_join(select(map_ldref, chr, pos))
#
#
#
try(calculate_ldpred_pgs(sumstats = final_sumstats, n_core = core_count,  plink_rds = rds_path, sd_y = 1, pheno_name = pheno, build = 'hg19'))
#
#


################################################
## repeat for raw Neale sumstats
sumstats <- '/Dedicated/jmichaelson-sdata/gwas_summary_stats/UKBIOBANK/IEU_Neale_TownsendDeprivationIndex_ukb-a-44.rsid.tsv'
pheno = "IEU_Neale_TownsendDeprivationIndex_ukb-a-44"
message('Working on ', pheno)

#### read in sumstats
ss = read_tsv(sumstats)

#### calculate pgs ####
str_c("number of rows in the raw summary statistics file is: ", nrow(ss))

### 3: give the output path
output_path <- str_c(my_wd, pheno, "/")
str_c("the final output path is: ", output_path)
dir.create(output_path)
setwd(output_path)

### normalize column names
ss <- ss %>%
  rename(chr = CHR, pos = POS, a1 = A1, a0 = A2, n_eff = N, beta = BETA, beta_se = SE, p = P) %>%
  mutate(beta = -1 * beta) ## flip betas, based on LDSC genetic correlations - the Neale sumstats are flipped

##
final_sumstats = ss %>% 
  select(chr,
          pos,
          a1, # a1 is the effect allele
          a0,
          n_eff,
          beta,
          beta_se,
          p) %>%
  filter(a1 != a0) %>%
  drop_na() %>%
  inner_join(select(bim_snp, chr, pos)) %>%
  inner_join(select(map_ldref, chr, pos))


## compute PGS
try(calculate_ldpred_pgs(sumstats = final_sumstats, n_core = core_count, plink_rds = rds_path, sd_y = 0, pheno_name = pheno, build = 'hg19'))

rm(ss)
rm(final_sumstats)
rm(sumstats)
rm(pheno)



#### household income ####
ss <- read_table('/Dedicated/jmichaelson-wdata/trthomas/abcd_spark_merge/polygenic_scores/sumstats/UK_Biobank/Lothian_Birth/SES/Hill2016_UKB_Income_summary_results_21112016.txt')
pheno = '2016_Hill_householdIncome'
ss = ss %>%
  drop_na() %>%
  filter(nchar(Effect_allele) == 1 & nchar(Other_allele) == 1)
#
# ##
# str_c("number of rows in the raw summary statistics file is: ", nrow(ss))
# head(ss)
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
  rename(P = `P-value`) %>%
  mutate(N = 96900,
         se_beta = abs(Beta / qnorm(P / 2))
  ) %>%
  select(
    chr = Chromosome,
    pos = Position,
    a1 = Effect_allele, # a1 is the effect allele
    a0 = Other_allele,
    n_eff = N,
    beta = Beta,
    beta_se = se_beta,
    p = P) %>%
  filter(a1 != a0) %>%
  drop_na() %>%
  inner_join(select(bim_snp, chr, pos)) %>%
  inner_join(select(map_ldref, chr, pos))
#
#
#
try(calculate_ldpred_pgs(sumstats = final_sumstats, n_core = core_count,  plink_rds = rds_path, sd_y = 1, pheno_name = pheno, build = 'hg19'))


#
#
## -------------------------------------------------------
## -------------------------------------------------------
## -------------------------------------------------------
## -------------------------------------------------------
#
#### left-handed ####
ss <- read_table('/Dedicated/jmichaelson-wdata/lcasten/misc/gwas_summary_stats/language/ukbb_left_handed_1707_2.v1.0.fastGWA')
pheno = 'UKBB_fastGWA_leftHanded'
#
##
str_c("number of rows in the raw summary statistics file is: ", nrow(ss))
head(ss)
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
  select(
    chr = CHR,
    pos = POS,
    a1 = A1, # a1 is the effect allele
    a0 = A2,
    n_eff = N,
    beta = BETA,
    beta_se = SE,
    p = P)%>%
  filter(a1 != a0) %>%
  drop_na() %>%
  inner_join(select(bim_snp, chr, pos)) %>%
  inner_join(select(map_ldref, chr, pos))
#
#
##
try(calculate_ldpred_pgs(sumstats = final_sumstats,  n_core = core_count, plink_rds = rds_path, sd_y = 0, pheno_name = pheno, build = 'hg19'))
#
#
#
#
#
#
# #### risk pc1 ####
ss <- read_table('/Dedicated/jmichaelson-sdata/gwas_summary_stats/ssgac/2019/2019_ssgac_risk_pc1.tsv')
pheno = '2019_ssgac_risk_pc1'
#
##
str_c("number of rows in the raw summary statistics file is: ", nrow(ss))
head(ss)
#
# cases = 12160
# controls = 13145
# ss$N <- 4 / (1 / cases + 1 / controls)
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
  mutate(N = 315894) %>%
  select(
    chr = CHR,
    pos = POS,
    a1 = A1, # a1 is the effect allele
    a0 = A2,
    n_eff = N,
    beta = BETA,
    beta_se = SE,
    p = P)%>%
  filter(a1 != a0) %>%
  drop_na() %>%
  inner_join(select(bim_snp, chr, pos)) %>%
  inner_join(select(map_ldref, chr, pos))
#
#
##
try(calculate_ldpred_pgs(sumstats = final_sumstats, n_core = core_count,  plink_rds = rds_path, sd_y = 0, pheno_name = pheno, build = 'hg19'))
#
#
#
#### risk tolerance ####
ss <- read_table('/Dedicated/jmichaelson-sdata/gwas_summary_stats/ssgac/2019/2019_ssgac_risk_tolerance.tsv')
pheno = '2019_ssgac_risk_tolerance'
#
##
str_c("number of rows in the raw summary statistics file is: ", nrow(ss))
head(ss)
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
  mutate(N =  431126) %>%
  select(
    chr = CHR,
    pos = POS,
    a1 = A1, # a1 is the effect allele
    a0 = A2,
    n_eff = N,
    beta = BETA,
    beta_se = SE,
    p = P) %>%
  filter(a1 != a0) %>%
  drop_na() %>%
  inner_join(select(bim_snp, chr, pos)) %>%
  inner_join(select(map_ldref, chr, pos))
#
#
##
try(calculate_ldpred_pgs(sumstats = final_sumstats, n_core = core_count,  plink_rds = rds_path, sd_y = 0, pheno_name = pheno, build = 'hg19'))
#
#
#
# #### alcohol abuse ####
f = '/Dedicated/jmichaelson-wdata/trthomas/abcd_spark_merge/polygenic_scores/sumstats/Psychiatric_Genomics_Consortium/alcohol_dependency-PGC-2018.txt'
ss <- read_table(f)
pheno = basename(f)
pheno = str_remove_all(pheno, '.txt')
#
##
str_c("number of rows in the raw summary statistics file is: ", nrow(ss))
head(ss)
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
  ) %>%
  select(
    chr = CHR,
    pos = BP,
    a1 = A1, # a1 is the effect allele
    a0 = A2,
    n_eff = Neff,
    beta = beta,
    beta_se = se_beta,
    p = P) %>%
  filter(a1 != a0) %>%
  drop_na() %>%
  inner_join(select(bim_snp, chr, pos)) %>%
  inner_join(select(map_ldref, chr, pos))
#
#
# ##
try(calculate_ldpred_pgs(sumstats = final_sumstats,  n_core = core_count, plink_rds = rds_path, sd_y = 2, pheno_name = pheno, build = 'hg19'))
#
#
#
#
#
#
#
# #### personality traits ####
personality = list.files('/Dedicated/jmichaelson-sdata/gwas_summary_stats/gpc', pattern = '_filt_reform', full.names = T)
personality = as.character(personality)
personality = personality[str_detect(personality, pattern = 'extraversion|neuro', negate = T)]
personality
#
for (f in personality) {
  tmp_pheno = basename(f)
  tmp_pheno = str_split(tmp_pheno, pattern = '-', simplify = T)[,1]
#
  tmp_pheno = str_c('GPC_', tmp_pheno)
#
  ss <- read_table(file = f)
  # names(ss) = c('rsid','CHR','BP',	'A1',	'A2',	'BETA',	 'SE','PVALUE',	'INFO',	'NCOH',	'MAF')
#
  ss = ss %>%
    # mutate(A1 = toupper(A1),
    #        A2 = toupper(A2)
    # ) %>%'
    rename(chr = CHR, pos = POS) %>%
    filter(INFO >= .8) %>%
    inner_join(map_ld) %>%
    # select(-c(chr, BP)) %>%
    mutate(BETA = case_when(A1 == a1 & A2 == a0 ~ BETA,
                            A1 == a0 & A2 == a1 ~ -1 * BETA,
                            TRUE ~ NA_real_)
    ) %>%
    drop_na(BETA) %>%
    select(-c(a0, a1)) %>%
    select(-c(af_UKBB, ld, pos_hg18, pos_hg38)) %>%
    relocate(pos, .after = chr)
#
  pheno = tmp_pheno
#
  ##
  str_c("number of rows in the raw summary statistics file is: ", nrow(ss))
  head(ss)
#
  ### 3: give the output path
  output_path <- str_c(my_wd, pheno, "/")
  str_c("the final output path is: ", output_path)
  system(str_c('rm -r ', output_path))
  dir.create(output_path)  ### 3: give the output path
  setwd(output_path)
#
  ##### 4: modify/create/filter the summary statistics columns to have: chr, pos, a1 (effect allele), a0, beta, beta_se, p, n_eff
  final_sumstats <- ss %>%
    mutate(N = 17375) %>%
    select(
      chr = chr,
      pos,
      a1 = A1, # a1 is the effect allele
      a0 = A2,
      n_eff = N,
      beta = BETA,
      beta_se = SE,
      p = P)%>%
    filter(a1 != a0) %>%
    drop_na() %>%
    inner_join(select(bim_snp, chr, pos)) %>%
    inner_join(select(map_ldref, chr, pos))
#
#
  ##
  try(calculate_ldpred_pgs(sumstats = final_sumstats, n_core = core_count,  plink_rds = rds_path, sd_y = 0, pheno_name = pheno, build = 'hg19'))
}
#
#
####
#
#### GPC2
personality = list.files('/Dedicated/jmichaelson-sdata/gwas_summary_stats/GPC2', pattern = 'full.txt$', full.names = T)
#
for (f in personality) {
  tmp_pheno = basename(f)
  tmp_pheno = str_split(tmp_pheno, pattern = '[.]', simplify = T)[,2]
#
  tmp_pheno = str_c('GPC2_', tmp_pheno)
#
  ss <- read_table(file = f, col_names = T)
  names(ss) = c('rsid','CHR','BP',	'A1',	'A2',	'BETA',	 'SE','PVALUE',	'NCOH',	'MAF')
#
  ss = ss %>%
    mutate(A1 = toupper(A1),
           A2 = toupper(A2)
    ) %>%
    inner_join(map_ld) %>%
    select(-c(chr, BP)) %>%
    mutate(BETA = case_when(A1 == a1 & A2 == a0 ~ BETA,
                            A1 == a0 & A2 == a1 ~ -1 * BETA,
                            TRUE ~ NA_real_)
    ) %>%
    drop_na(BETA) %>%
    select(-c(a0, a1)) %>%
    select(-c(af_UKBB, ld, pos_hg18, pos_hg38)) %>%
    relocate(pos, .after = CHR)
#
  pheno = tmp_pheno
#
  ##
  str_c("number of rows in the raw summary statistics file is: ", nrow(ss))
  head(ss)
#
  ### 3: give the output path
  output_path <- str_c(my_wd, pheno, "/")
  str_c("the final output path is: ", output_path)
  system(str_c('rm -r ', output_path))
  dir.create(output_path)  ### 3: give the output path
  setwd(output_path)
#
  ##### 4: modify/create/filter the summary statistics columns to have: chr, pos, a1 (effect allele), a0, beta, beta_se, p, n_eff
  final_sumstats <- ss %>%
    mutate(N = 63661) %>%
    select(
      chr = CHR,
      pos,
      a1 = A1, # a1 is the effect allele
      a0 = A2,
      n_eff = N,
      beta = BETA,
      beta_se = SE,
      p = PVALUE) %>%
    filter(a1 != a0) %>%
    drop_na() %>%
    inner_join(select(bim_snp, chr, pos)) %>%
    inner_join(select(map_ldref, chr, pos))
#
#
  ##
  try(calculate_ldpred_pgs(sumstats = final_sumstats,  n_core = core_count, plink_rds = rds_path, sd_y = 0, pheno_name = pheno, build = 'hg19'))
}
#
#
#
#
# #
# # #### height ####
ss <- read_table('/Dedicated/jmichaelson-sdata/gwas_summary_stats/anthropomorphic/height/GIANT_HEIGHT_YENGO_2022_GWAS_SUMMARY_STATS_ALL.gz')
pheno = '2022_Yengo_height'
#
# ss = ss %>%
#
#   mutate(A1 = toupper(EFFECT_ALLELE),
#          A2 = toupper()
#   )
#
##
str_c("number of rows in the raw summary statistics file is: ", nrow(ss))
head(ss)
#
# cases = 12160
# controls = 13145
# ss$N <- 4 / (1 / cases + 1 / controls)
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
  select(
    chr = CHR,
    pos = POS,
    a1 = EFFECT_ALLELE, # a1 is the effect allele
    a0 = OTHER_ALLELE,
    n_eff = N,
    beta = BETA,
    beta_se = SE,
    p = P) %>%
  filter(a1 != a0) %>%
  drop_na() %>%
  inner_join(select(bim_snp, chr, pos)) %>%
  inner_join(select(map_ldref, chr, pos))
#
#
##
try(calculate_ldpred_pgs(sumstats = final_sumstats,  n_core = core_count, plink_rds = rds_path, sd_y = 0, pheno_name = pheno, build = 'hg19'))
#
#
# #### brain imaging measures ####
# # ## -------------------------------------------
# # #### compute global volume pgs ####
# # ## -------------------------------------------
ss <- read_table('/Dedicated/jmichaelson-sdata/gwas_summary_stats/ENIGMA_brain_imaging/global_measures/raw/ENIGMA3_mixed_se_wo_Mean_Full_SurfArea_20190429.txt.gz')
pheno = 'ENIGMA_globalSurfaceArea'
#
ss = ss %>%
  mutate(A1 = toupper(A1),
         A2 = toupper(A2)
  )
#
##
str_c("number of rows in the raw summary statistics file is: ", nrow(ss))
head(ss)
#
# cases = 12160
# controls = 13145
# ss$N <- 4 / (1 / cases + 1 / controls)
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
  select(
    chr = CHR,
    pos = BP,
    a1 = A1, # a1 is the effect allele
    a0 = A2,
    n_eff = N,
    beta = BETA1,
    beta_se = SE,
    p = P) %>%
  filter(a1 != a0) %>%
  drop_na() %>%
  inner_join(select(bim_snp, chr, pos)) %>%
  inner_join(select(map_ldref, chr, pos))
#
#
##
try(calculate_ldpred_pgs(sumstats = final_sumstats, n_core = core_count,  plink_rds = rds_path, sd_y = 0, pheno_name = pheno, build = 'hg19'))
#
#
# # ## -------------------------------------------
# # #### compute global volume pgs ####
# # ## -------------------------------------------
ss <- read_table('/Dedicated/jmichaelson-sdata/gwas_summary_stats/ENIGMA_brain_imaging/global_measures/raw/ENIGMA3_mixed_se_wo_Mean_Full_Thickness_20190429.txt.gz')
pheno = 'ENIGMA_globalThickness'
#
ss = ss %>%
  mutate(A1 = toupper(A1),
         A2 = toupper(A2)
  )
#
##
str_c("number of rows in the raw summary statistics file is: ", nrow(ss))
head(ss)
#
# cases = 12160
# controls = 13145
# ss$N <- 4 / (1 / cases + 1 / controls)
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
  select(
    chr = CHR,
    pos = BP,
    a1 = A1, # a1 is the effect allele
    a0 = A2,
    n_eff = N,
    beta = BETA1,
    beta_se = SE,
    p = P)%>%
  filter(a1 != a0) %>%
  drop_na() %>%
  inner_join(select(bim_snp, chr, pos)) %>%
  inner_join(select(map_ldref, chr, pos))
#
#
##
try(calculate_ldpred_pgs(sumstats = final_sumstats,  n_core = core_count, plink_rds = rds_path, sd_y = 0, pheno_name = pheno, build = 'hg19'))
#
#


## reproductive related
message('')
message('')
message('')
message('')
message('--------------------')
message('Working on reproductive related scores now...')
system('cat /Dedicated/jmichaelson-sdata/gwas_summary_stats/reproductive/Mills-NatureHumanBehavior2021/README.txt')
files = list.files('/Dedicated/jmichaelson-sdata/gwas_summary_stats/reproductive/Mills-NatureHumanBehavior2021',
                   pattern = 'tsv.gz', full.names = TRUE)
fmap <- read_csv('/Dedicated/jmichaelson-sdata/gwas_summary_stats/reproductive/Mills-NatureHumanBehavior2021/gwas_catalog_id_phenotype_map.csv')
for (f in files) {
  id = basename(f)
  id = str_split(id, pattern = '_', simplify = TRUE)[,1]
  
  i = which(fmap$gwas_catalog_id == id)
  
  pheno = fmap$trait[i]
  sample_size = fmap$N[i]
  message('Working on ', pheno)
  
  ss <- read_table(file = f)
  
  ss$N = sample_size
  
  
  ##
  str_c("number of rows in the raw summary statistics file is: ", nrow(ss))
  
  ### 3: give the output path
  output_path <- str_c(my_wd, pheno, "/")
  str_c("the final output path is: ", output_path)
  system(str_c('rm -r ', output_path))
  dir.create(output_path)
  setwd(output_path)
  #
  ##### 4: modify/create/filter the summary statistics columns to have: chr, pos, a1 (effect allele), a0, beta, beta_se, p, n_eff
  final_sumstats <- ss %>%
    select(
      chr = chromosome,
      pos = base_pair_location,
      a1 = effect_allele, # a1 is the effect allele
      a0 = other_allele,
      n_eff = N,
      beta = beta,
      beta_se = standard_error,
      p = p_value) %>%
    filter(a1 != a0) %>%
    drop_na() %>%
    filter(nchar(a1) == 1) %>%
    filter(nchar(a0) == 1) %>%
    inner_join(select(bim_snp, chr, pos)) %>%
    inner_join(select(map_ldref, chr, pos))
  #
  ##
  try(calculate_ldpred_pgs(sumstats = final_sumstats, n_core = core_count,  plink_rds = rds_path, sd_y = 0, pheno_name = pheno, build = 'hg19'))
}



# ## ------------------------------
# #### UKBB cog g ####
# ## ---------------------------
sumstats_filepath <- "/Dedicated/jmichaelson-wdata/trthomas/abcd_spark_merge/polygenic_scores/sumstats/UK_Biobank/Lothian_Birth/cog_traits/TRT_processing_2022/cog_gFactor-UKB-2020"
sumstats <- read_tsv(sumstats_filepath, guess_max = 1000000)
pheno = basename(sumstats_filepath)
head(sumstats)
#
output_path <- str_c(my_wd, pheno, "/")
str_c("the final output path is: ", output_path)
system(str_c('rm -r ', output_path))
dir.create(output_path)
setwd(output_path)
final_sumstats <- sumstats %>%
  select(
    chr = chr,
    pos = pos,
    a1 = a1, # a1 is the effect allele
    a0 = a0,
    n_eff = n_eff,
    beta = beta,
    beta_se = beta_se,
    p = p) %>%
  filter(a1 != a0) %>%
  drop_na() %>%
  inner_join(select(bim_snp, chr, pos)) %>%
  inner_join(select(map_ldref, chr, pos))
#
#
##
try(calculate_ldpred_pgs(sumstats = final_sumstats, n_core = core_count,  plink_rds = rds_path, sd_y = 1, pheno_name = pheno, build = 'hg19'))
#
# #
# #
# ## -----------------------
# #### ukbb verbal numeric reas
# ## -------------------------
sumstats_filepath <- "/Dedicated/jmichaelson-wdata/trthomas/abcd_spark_merge/polygenic_scores/sumstats/UK_Biobank/Lothian_Birth/cog_traits/TRT_processing_2022/cog_verbal_numerical_reasoning-UKB-2020"
sumstats <- read_tsv(sumstats_filepath, guess_max = 1000000)
pheno = basename(sumstats_filepath)
head(sumstats)
#
output_path <- str_c(my_wd, pheno, "/")
str_c("the final output path is: ", output_path)
system(str_c('rm -r ', output_path))
dir.create(output_path)
setwd(output_path)
final_sumstats <- sumstats %>%
  select(
    chr = chr,
    pos = pos,
    a1 = a1, # a1 is the effect allele
    a0 = a0,
    n_eff = n_eff,
    beta = beta,
    beta_se = beta_se,
    p = p) %>%
  filter(a1 != a0) %>%
  drop_na() %>%
  inner_join(select(bim_snp, chr, pos)) %>%
  inner_join(select(map_ldref, chr, pos))
#
#
##
try(calculate_ldpred_pgs(sumstats = final_sumstats,  n_core = core_count, plink_rds = rds_path, sd_y = 1, pheno_name = pheno, build = 'hg19'))
#
#
## -------------
#### matrix reas ####
sumstats_filepath <- "/Dedicated/jmichaelson-wdata/trthomas/abcd_spark_merge/polygenic_scores/sumstats/UK_Biobank/Lothian_Birth/cog_traits/TRT_processing_2022/cog_matrix-UKB-2020"
sumstats <- read_tsv(sumstats_filepath, guess_max = 1000000)
pheno = basename(sumstats_filepath)
head(sumstats)
#
output_path <- str_c(my_wd, pheno, "/")
str_c("the final output path is: ", output_path)
system(str_c('rm -r ', output_path))
dir.create(output_path)
setwd(output_path)
final_sumstats <- sumstats %>%
  select(
    chr = chr,
    pos = pos,
    a1 = a1, # a1 is the effect allele
    a0 = a0,
    n_eff = n_eff,
    beta = beta,
    beta_se = beta_se,
    p = p)%>%
  filter(a1 != a0) %>%
  drop_na() %>%
  inner_join(select(bim_snp, chr, pos)) %>%
  inner_join(select(map_ldref, chr, pos))
#
#
##
try(calculate_ldpred_pgs(sumstats = final_sumstats, n_core = core_count,  plink_rds = rds_path, sd_y = 1, pheno_name = pheno, build = 'hg19'))
#
# ## -------------------
# #### reaction time
sumstats_filepath <- "/Dedicated/jmichaelson-wdata/trthomas/abcd_spark_merge/polygenic_scores/sumstats/UK_Biobank/Lothian_Birth/cog_traits/TRT_processing_2022/cog_reaction_time-UKB-2020"
sumstats <- read_tsv(sumstats_filepath, guess_max = 1000000)
pheno = basename(sumstats_filepath)
head(sumstats)
#
output_path <- str_c(my_wd, pheno, "/")
str_c("the final output path is: ", output_path)
system(str_c('rm -r ', output_path))
dir.create(output_path)
setwd(output_path)
final_sumstats <- sumstats %>%
  select(
    chr = chr,
    pos = pos,
    a1 = a1, # a1 is the effect allele
    a0 = a0,
    n_eff = n_eff,
    beta = beta,
    beta_se = beta_se,
    p = p)%>%
  filter(a1 != a0) %>%
  drop_na() %>%
  inner_join(select(bim_snp, chr, pos)) %>%
  inner_join(select(map_ldref, chr, pos))
#
#
##
try(calculate_ldpred_pgs(sumstats = final_sumstats,  n_core = core_count, plink_rds = rds_path, sd_y = 1, pheno_name = pheno, build = 'hg19'))
#
#
## -------------------------------
#### memory ####
sumstats_filepath <- "/Dedicated/jmichaelson-wdata/trthomas/abcd_spark_merge/polygenic_scores/sumstats/UK_Biobank/Lothian_Birth/cog_traits/TRT_processing_2022/cog_memory-UKB-2020"
sumstats <- read_tsv(sumstats_filepath, guess_max = 1000000)
pheno = basename(sumstats_filepath)
head(sumstats)
#
output_path <- str_c(my_wd, pheno, "/")
str_c("the final output path is: ", output_path)
system(str_c('rm -r ', output_path))
dir.create(output_path)
setwd(output_path)
final_sumstats <- sumstats %>%
  select(
    chr = chr,
    pos = pos,
    a1 = a1, # a1 is the effect allele
    a0 = a0,
    n_eff = n_eff,
    beta = beta,
    beta_se = beta_se,
    p = p)%>%
  filter(a1 != a0) %>%
  drop_na() %>%
  inner_join(select(bim_snp, chr, pos)) %>%
  inner_join(select(map_ldref, chr, pos))
#
#
##
try(calculate_ldpred_pgs(sumstats = final_sumstats, n_core = core_count,  plink_rds = rds_path, sd_y = 1, pheno_name = pheno, build = 'hg19'))
#
# ## ----------------------
# #### adhd ####
ss_file = "/Dedicated/jmichaelson-wdata/trthomas/abcd_spark_merge/polygenic_scores/sumstats/Psychiatric_Genomics_Consortium/ADHD-PGC-2019"
ss <- read_tsv(ss_file)
# ss = read_tsv('/wdata/trthomas/abcd_spark_merge/polygenic_scores/sumstats/Psychiatric_Genomics_Consortium/ADHD-PGC-2019')
names(ss)
pheno = basename(ss_file)
#
##
str_c("number of rows in the raw summary statistics file is: ", nrow(ss))
head(ss)
#
cases <- 20183 # input number of cases
controls <- 35191 #input number of controls
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
print('head of OR SE')
print(
  head(final_sumstats$SE))
#
print('head of beta SE')
print(
  head(final_sumstats$se_beta))
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
try(calculate_ldpred_pgs(sumstats = final_sumstats, n_core = core_count,  plink_rds = rds_path, sd_y = 2, pheno_name = pheno, build = 'hg19'))
#
#
#
#
#
#
#
#
#
# #
# # ## -------------------------------------------
# # #### compute T1D pgs ####
# # ## -------------------------------------------
t1d <- read_tsv('/Dedicated/jmichaelson-wdata/lcasten/misc/gwas_summary_stats/imsgc_ms/2021_T1D.tsv')
pheno = '2021_t1d_chiou'
#
##
str_c("number of rows in the raw summary statistics file is: ", nrow(t1d))
head(t1d)
#
cases = 18942
controls = 501638
t1d$N <- 4 / (1 / cases + 1 / controls)
#
### 3: give the output path
output_path <- str_c(my_wd, pheno, "/")
str_c("the final output path is: ", output_path)
dir.create(output_path)
setwd(output_path)
#
##### 4: modify/create/filter the summary statistics columns to have: chr, pos, a1 (effect allele), a0, beta, beta_se, p, n_eff
final_sumstats <- t1d %>%
  select(
    chr = CHR,
    pos = POS,
    a1 = A1, # a1 is the effect allele
    a0 = A2,
    n_eff = N,
    beta = BETA,
    beta_se = SE,
    p = P) %>%
  filter(a1 != a0) %>%
  drop_na() %>%
  inner_join(select(bim_snp, chr, pos)) %>%
  inner_join(select(map_ldref, chr, pos))
#
#
##
try(calculate_ldpred_pgs(sumstats = final_sumstats, n_core = core_count,  plink_rds = rds_path, sd_y = 2, pheno_name = pheno, build = 'hg19'))
#
#
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
# #
# #
# # ### 2: give the path to the summary statistics input file
sumstats <- read_tsv('/Dedicated/jmichaelson-wdata/lcasten/misc/gwas_summary_stats/pgc_bd3/pgc-bip2021-all/pgc-bip2021-all.tsv', guess_max = 10000000)
pheno_name = '2021_PGC_bipolar_I_II'
# pheno_name = str_remove(pheno_name, ".tsv")
str_c("number of rows in the raw summary statistics file is: ", nrow(sumstats))
head(sumstats)
#
cases = 41917
controls = 371549
sumstats$N <- 4 / (1 / cases + 1 / controls)
#
### 3: give the output path
output_path <- str_c(my_wd, pheno_name, "/")
str_c("the final output path is: ", output_path)
system(str_c('rm -r ', output_path))
dir.create(output_path)
setwd(output_path)
#
##### 4: modify/create/filter the summary statistics columns to have: chr, pos, a1 (effect allele), a0, beta, beta_se, p, n_eff
final_sumstats <- sumstats %>%
  mutate(N =  4 / (1 / NCAS + 1 / NCON)) %>%
  filter(IMPINFO >= 0.8 & N >= .75 * max(N)) %>%
  select(
    chr = CHR,
    pos = POS,
    a1 = A1, # a1 is the effect allele
    a0 = A2,
    n_eff = N,
    beta = BETA,
    beta_se = SE,
    p = PVAL) %>%
  filter(a1 != a0) %>%
  drop_na() %>%
  inner_join(select(bim_snp, chr, pos)) %>%
  inner_join(select(map_ldref, chr, pos))
#
head(final_sumstats)
#
### calculate bipolar prs ####
pheno = pheno_name
try(calculate_ldpred_pgs(sumstats = final_sumstats,  n_core = core_count, plink_rds = rds_path, sd_y = 2, pheno_name = pheno, build = 'hg19'))
#
# # ## -------------------------------------------
# # #### compute asd pgs ####
# # ## -------------------------------------------
ss <- read_tsv("/Dedicated/jmichaelson-wdata/trthomas/abcd_spark_merge/polygenic_scores/sumstats/Psychiatric_Genomics_Consortium/autism-PGC-2019")
pheno = '2019_pgc_autism'
#
#
##
str_c("number of rows in the raw summary statistics file is: ", nrow(ss))
head(ss)
#
cases <- 18381 # input number of cases
controls <- 27969 #input number of controls
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
print('head of OR SE')
print(
  head(final_sumstats$SE))
#
print('head of beta SE')
print(
  head(final_sumstats$se_beta))
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
try(calculate_ldpred_pgs(sumstats = final_sumstats, n_core = core_count,  plink_rds = rds_path, sd_y = 2, pheno_name = pheno, build = 'hg19'))
#
# # ## -------------------------------------------
# # #### compute IBD pgs ####
# # ## -------------------------------------------
ss <- read_tsv('/Dedicated/jmichaelson-sdata/gwas_summary_stats/intestinal/2017_delange_ibd.tsv')
pheno = '2017_delange_ibd'
#
##
str_c("number of rows in the raw summary statistics file is: ", nrow(ss))
head(ss)
#
cases = 12160
controls = 13145
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
  select(
    chr = CHR,
    pos = POS,
    a1 = A1, # a1 is the effect allele
    a0 = A2,
    n_eff = N,
    beta = BETA,
    beta_se = SE,
    p = P)%>%
  filter(a1 != a0) %>%
  drop_na() %>%
  inner_join(select(bim_snp, chr, pos)) %>%
  inner_join(select(map_ldref, chr, pos))
#
#
##
try(calculate_ldpred_pgs(sumstats = final_sumstats,  n_core = core_count, plink_rds = rds_path, sd_y = 2, pheno_name = pheno, build = 'hg19'))
#
#
#
####################
#### depression ####
####################
sumstats_filepath <- "/Dedicated/jmichaelson-wdata/trthomas/abcd_spark_merge/polygenic_scores/sumstats/Psychiatric_Genomics_Consortium/depression-PGC-2019"
sumstats <- read_delim(sumstats_filepath, guess_max = 1000000, delim = " ")
pheno = basename(sumstats_filepath)
head(sumstats)
#
### 4: need effective N. For quantitative traits, it is just total N. For binary traits, calculate using N cases and N controls
n_case <- NA
n_control <- NA
n_eff <- NA
n_case <- 170756 #input number of cases
n_control <- 329443 #input number of controls
n_eff <- round(4 / (1 / n_case + 1 / n_control), digits = 0)
str_c("effective N is: ", n_eff)
#
sumstats <- sumstats %>%
  mutate(n_eff = n_eff) %>%
  mutate(A1 = toupper(A1)) %>%
  mutate(A2 = toupper(A2))
#
sumstats <- sumstats %>%
  select(
    rsid = MarkerName,
    a1 = A1, # a1 is the effect allele
    a0 = A2,
    n_eff,
    beta = LogOR,
    beta_se = StdErrLogOR,
    p = P)
#
str_c("Chromsome and position not provided in this file, so going to use the map_ldref to pull them (for only the ones used in PGS)")
#
# map_ldref <- readRDS("/Dedicated/jmichaelson-wdata/trthomas/abcd_spark_merge/polygenic_scores/ldpred2/ld_ref/map.rds")
# map_ldref <- map_ldref %>%
#   select(chr, pos, a0, a1, rsid, af_UKBB, ld)
#
final_sumstats <- sumstats %>%
  inner_join(map_ldref %>% select(rsid, chr, pos), by = "rsid") %>%
  select(chr, pos, a1, a0, n_eff, beta, beta_se, p) %>%
  filter(a1 != a0) %>%
  drop_na() %>%
  inner_join(select(bim_snp, chr, pos)) %>%
  inner_join(select(map_ldref, chr, pos))
#
#
### 3: give the output path
output_path <- str_c(my_wd, pheno, "/")
str_c("the final output path is: ", output_path)
system(str_c('rm -r ', output_path))
dir.create(output_path)
setwd(output_path)
#
##
try(calculate_ldpred_pgs(sumstats = final_sumstats, n_core = core_count,  plink_rds = rds_path, sd_y = 2, pheno_name = pheno, build = 'hg19'))
#
#
# # ###############################
# # #### cannabis use disorder ####
# # ###############################
sumstats_filepath <- "/Dedicated/jmichaelson-wdata/trthomas/abcd_spark_merge/polygenic_scores/sumstats/Psychiatric_Genomics_Consortium/cannabis_use_disorder-PGC-2020"
sumstats <- read_table2(sumstats_filepath, guess_max = 1000000)
pheno = basename(sumstats_filepath)
# head(sumstats)
#
# ### 4: need effective N. For quantitative traits, it is just total N. For binary traits, calculate using N cases and N controls
n_case <- 14080
n_control <- 343726
n_eff <- round(4 / (1 / n_case + 1 / n_control), digits = 0)
str_c("effective N is: ", n_eff)
#
sumstats <- sumstats %>%
  mutate(n_eff = n_eff)
#
final_sumstats <- sumstats %>%
  select(
    chr = CHR,
    pos = BP,
    a1 = A1, # a1 is the effect allele
    a0 = A2,
    n_eff,
    beta = Beta,
    beta_se = SE,
    p = P) %>%
  filter(a1 != a0) %>%
  drop_na() %>%
  inner_join(select(bim_snp, chr, pos)) %>%
  inner_join(select(map_ldref, chr, pos))
#
# head(final_sumstats)
#
# ### 3: give the output path
output_path <- str_c(my_wd, pheno, "/")
str_c("the final output path is: ", output_path)
system(str_c('rm -r ', output_path))
dir.create(output_path)
setwd(output_path)
#
try(calculate_ldpred_pgs(sumstats = final_sumstats,  n_core = core_count, plink_rds = rds_path, sd_y = 2, pheno_name = pheno, build = 'hg19'))
#
#
#
# ##################
# #### Insomnia ####
# ##################
sumstats_filepath <- "/Dedicated/jmichaelson-wdata/lcasten/misc/gwas_summary_stats/ctg_sumstats_posthuma/insomnia/insomnia_ukb2b_EUR_sumstats_20190311_with_chrX_mac_100.txt.gz"
sumstats <- read_table2(sumstats_filepath, guess_max = 1000000)
pheno = 'CTGlab_2019_insomnia'
head(sumstats)
#
sumstats <- sumstats %>%
  mutate(n_eff = NMISS) %>%
  mutate(beta = log(OR),
         se_beta = abs(beta / qnorm(P / 2))
  )
#
sumstats = sumstats %>%
  filter(MAF >= 0.01 & INFO_UKB >= .8)
#
final_sumstats <- sumstats %>%
  select(
    chr = CHR,
    pos = BP,
    a1 = A1, # a1 is the effect allele
    a0 = A2,
    n_eff,
    beta = beta,
    beta_se = se_beta,
    p = P) %>%
  filter(a1 != a0) %>%
  drop_na() %>%
  inner_join(select(bim_snp, chr, pos)) %>%
  inner_join(select(map_ldref, chr, pos))
#
head(final_sumstats)
#
### 3: give the output path
output_path <- str_c(my_wd, pheno, "/")
str_c("the final output path is: ", output_path)
system(str_c('rm -r ', output_path))
dir.create(output_path)
setwd(output_path)
#
##
try(calculate_ldpred_pgs(sumstats = final_sumstats,  n_core = core_count, plink_rds = rds_path, sd_y = 2, pheno_name = pheno, build = 'hg19'))
# #
# #
# # ###################
# # #### longevity ####
# # ###################
sumstats_filepath <- "/Dedicated/jmichaelson-wdata/lcasten/misc/gwas_summary_stats/longevity_Deelen2019NatComm/Results_90th_percentile.txt"
sumstats <- read_table2(sumstats_filepath, guess_max = 1000000)
pheno = 'Deleen2019NatComm_longevity90thPercentile'
# head(sumstats)
#
sumstats = sumstats %>%
  mutate(EA = toupper(EA),
         NEA = toupper(NEA)
  )
#
final_sumstats <- sumstats %>%
  select(
    chr = Chr,
    pos = Position,
    a1 = EA, # a1 is the effect allele
    a0 = NEA,
    n_eff = Effective_N,
    beta = Beta,
    beta_se = SE,
    p = `P-value`
  ) %>%
  filter(a1 != a0) %>%
  drop_na() %>%
  inner_join(select(bim_snp, chr, pos)) %>%
  inner_join(select(map_ldref, chr, pos))
#
#
head(final_sumstats)
#
### 3: give the output path
output_path <- str_c(my_wd, pheno, "/")
str_c("the final output path is: ", output_path)
system(str_c('rm -r ', output_path))
dir.create(output_path)
setwd(output_path)
#
##
try(calculate_ldpred_pgs(sumstats = final_sumstats,  n_core = core_count, plink_rds = rds_path, sd_y = 2, pheno_name = pheno, build = 'hg19'))
#
#
# #####################
# #### PTSD - 2019 ####
# #####################
sumstats_filepath <- "/Dedicated/jmichaelson-sdata/gwas_summary_stats/PGC/tanner_archive/2019/2019_pgc_ptsd_eur_freeze2_overall.results.gz"
sumstats <- read_table2(sumstats_filepath, guess_max = 1000000)
pheno = 'PTSD-PGC-2019'
head(sumstats)
#
sumstats = sumstats %>%
  mutate(beta = log(OR),
         se_beta = abs(beta / qnorm(P / 2))
  ) %>%
  filter(INFO >= .8)
#
final_sumstats <- sumstats %>%
  select(
    chr = CHR,
    pos = BP,
    a1 = A1, # a1 is the effect allele
    a0 = A2,
    n_eff = Neff,
    beta = beta,
    beta_se = se_beta,
    p = P
  ) %>%
  filter(a1 != a0) %>%
  drop_na() %>%
  inner_join(select(bim_snp, chr, pos)) %>%
  inner_join(select(map_ldref, chr, pos))
#
#
head(final_sumstats)
#
### 3: give the output path
output_path <- str_c(my_wd, pheno, "/")
str_c("the final output path is: ", output_path)
system(str_c('rm -r ', output_path))
dir.create(output_path)
setwd(output_path)
#
##
try(calculate_ldpred_pgs(sumstats = final_sumstats,  n_core = core_count, plink_rds = rds_path, sd_y = 2, pheno_name = pheno, build = 'hg19'))
#
#
#
#
#
# #####################################################
# #### new PGS ####
# #####################################################
#
# #############################
# #### antisocial behavior ####
sumstats = '/Dedicated/jmichaelson-sdata/gwas_summary_stats/antisocial_behavior/Tielbeek_MolecularPsychiatry_2022/broadABC2022_Final_CombinedSex.TBL.gz'
ss = read_table(sumstats) %>%
  filter(effect_allele_frequency >= .01) %>%
  mutate(effect_allele = toupper(effect_allele),
         other_allele = toupper(other_allele)
  )
#
#### calculate pgs ####
pheno = "broad_antisocial_behavior"
#
### 3: give the output path
output_path <- str_c(my_wd, pheno, "/")
str_c("the final output path is: ", output_path)
dir.create(output_path)
setwd(output_path)
#
# ##### 4: modify/create/filter the summary statistics columns to have: chr, pos, a1 (effect allele), a0, beta, beta_se, p, n_eff
final_sumstats <- ss %>%
  filter(effect_allele_frequency >= .01) %>%
  filter(nchar(effect_allele) == 1 & nchar(other_allele) == 1)
#
##
final_sumstats$Neff = final_sumstats$n
#
final_sumstats = final_sumstats %>%
  select(chr = chromosome,
         pos = base_pair_location,
         a1 = effect_allele, # a1 is the effect allele
         a0 = other_allele,
         n_eff = Neff,
         beta = beta,
         beta_se = standard_error,
         p = p_value) %>%
  filter(a1 != a0) %>%
  drop_na() %>%
  inner_join(select(bim_snp, chr, pos)) %>%
  inner_join(select(map_ldref, chr, pos))
#
#
#
## compute PGS
try(calculate_ldpred_pgs(sumstats = final_sumstats,  n_core = core_count, plink_rds = rds_path, sd_y = 0, pheno_name = pheno, build = 'hg19'))
#
#
#
####################################################
################## iPSYCH GWAS' ####################
####################################################
#############################
#### ADHD late diagnosis ####
sumstats = '/Dedicated/jmichaelson-sdata/gwas_summary_stats/iPsych/2022/ADHD_subtypes/adulthood.sumstats.gz'
ss = read_table(file = sumstats)
#
pheno = "ADHD_diagnosed_late"
#
### 3: give the output path
output_path <- str_c(my_wd, pheno, "/")
str_c("the final output path is: ", output_path)
dir.create(output_path)
setwd(output_path)
#
##### 4: modify/create/filter the summary statistics columns to have: chr, pos, a1 (effect allele), a0, beta, beta_se, p, n_eff
final_sumstats <- ss %>%
  mutate(CHR = as.numeric(CHR)) %>%
  drop_na(CHR) %>%
  filter(nchar(A1) == 1 & nchar(A2) == 1) %>%
  filter(str_detect(A1, pattern = 'I|D', negate = T),
         str_detect(A2, pattern = 'I|D', negate = T)
  )
#
#
final_sumstats = final_sumstats %>%
  select(chr = CHR,
         pos = BP,
         a1 = A1, # a1 is the effect allele
         a0 = A2,
         n_eff = N,
         beta = BETA,
         beta_se = SE,
         p = P) %>%
  filter(a1 != a0) %>%
  drop_na() %>%
  inner_join(select(bim_snp, chr, pos)) %>%
  inner_join(select(map_ldref, chr, pos))
#
#
#
## compute PGS
try(calculate_ldpred_pgs(sumstats = final_sumstats,  n_core = core_count, plink_rds = rds_path, sd_y = 0, pheno_name = pheno, build = 'hg19'))
#
######################
#### schoolgrades ####
sg_map = data.frame(pheno = c('schoolE1.better_overall_school_performance',
                              'schoolE2.better_language_than_math',
                              'schoolE3.better_oral_than_written_exams',
                              'schoolE4.better_first_language_than_secondary'
)
)
#
for (i in 1:4) {
  sumstats = str_c('/Dedicated/jmichaelson-sdata/gwas_summary_stats/iPsych/2023/schoolgrades/schoolgrades.E', i, '.sumstats.gz')
  ss = read_table(file = sumstats)
#
  pheno = sg_map$pheno[i]
#
  ### 3: give the output path
  output_path <- str_c(my_wd, pheno, "/")
  str_c("the final output path is: ", output_path)
  dir.create(output_path)
  setwd(output_path)
#
  ##### 4: modify/create/filter the summary statistics columns to have: chr, pos, a1 (effect allele), a0, beta, beta_se, p, n_eff
  final_sumstats <- ss %>%
    mutate(CHR = as.numeric(CHR)) %>%
    drop_na(CHR) %>%
    filter(INFO >= .8) %>%
    filter(nchar(A1) == 1 & nchar(A2) == 1) %>%
    filter(str_detect(A1, pattern = 'I|D', negate = T),
           str_detect(A2, pattern = 'I|D', negate = T)
    )
#
#
  final_sumstats = final_sumstats %>%
    select(chr = CHR,
           pos = BP,
           a1 = A1, # a1 is the effect allele
           a0 = A2,
           n_eff = N,
           beta = BETA,
           beta_se = SE,
           p = P)%>%
    filter(a1 != a0) %>%
    drop_na() %>%
    inner_join(select(bim_snp, chr, pos)) %>%
    inner_join(select(map_ldref, chr, pos))
#
#
#
  ## compute PGS
  try(calculate_ldpred_pgs(sumstats = final_sumstats,  n_core = core_count, plink_rds = rds_path, sd_y = 0, pheno_name = pheno, build = 'hg19'))
}






## -------------------
#### new measures ####
## -------------------
## additional PGC stuff
anorexia <- '/Dedicated/jmichaelson-sdata/gwas_summary_stats/PGC/tanner_archive/2019/raw/anorexia_nervosa/pgcAN2.2019-07.vcf.tsv'
tmp <- read_lines(file = anorexia)
to_skip <- min(which(str_detect(tmp, pattern = '^#', negate = T)))
ss <- read_tsv(anorexia, skip = to_skip - 1) %>%
  mutate(N = 4 / (1 / NCAS + 1 / NCON))
#
pheno = '2019-PGC-anorexia'
#
##
str_c("number of rows in the raw summary statistics file is: ", nrow(ss))
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
  filter(IMPINFO >= 0.8) %>%
  select(
    chr = CHROM,
    pos = POS,
    a1 = ALT, # a1 is the effect allele
    a0 = REF,
    n_eff = N,
    beta = BETA,
    beta_se = SE,
    p = PVAL) %>%
  filter(a1 != a0) %>%
  drop_na() %>%
  inner_join(select(bim_snp, chr, pos)) %>%
  inner_join(select(map_ldref, chr, pos))
#
##
try(calculate_ldpred_pgs(sumstats = final_sumstats,  n_core = core_count, plink_rds = rds_path, sd_y = 0, pheno_name = pheno, build = 'hg19'))
rm(ss)
rm(sumstats)
rm(final_sumstats)
rm(pheno)
#
################################
anxiety <- '/Dedicated/jmichaelson-sdata/gwas_summary_stats/iPsych/2019/anxiety_stress_disorders/daner_woautism_ad_sd8-sd6_cleaned.gz'
ss <- read_tsv(anxiety)
ss <- ss %>%
  filter(A1 %in% c('A', 'C', 'T', 'G') & A2 %in% c('A', 'C', 'T', 'G')) %>%
  mutate(beta = log(OR),
         se_beta = abs(beta / qnorm(P / 2))
  )
NCAS = 12665
NCON =  19225
ss$N = 4 / (1 / ss$Nca + 1 / ss$Nco)
ss
pheno = '2019-iPSYCH-anxiety_stress_disorder'
#
##
str_c("number of rows in the raw summary statistics file is: ", nrow(ss))
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
  filter(INFO >= 0.8) %>%
  select(
    chr = CHR,
    pos = BP,
    a1 = A1, # a1 is the effect allele
    a0 = A2,
    n_eff = N,
    beta = beta,
    beta_se = se_beta,
    p = P) %>%
  filter(a1 != a0) %>%
  drop_na() %>%
  inner_join(select(bim_snp, chr, pos)) %>%
  inner_join(select(map_ldref, chr, pos))
#
##
try(calculate_ldpred_pgs(sumstats = final_sumstats,  n_core = core_count, plink_rds = rds_path, sd_y = 0, pheno_name = pheno, build = 'hg19'))
rm(ss)
rm(sumstats)
rm(final_sumstats)
rm(pheno)
#
#
################################
## OCD
ocd <- '/Dedicated/jmichaelson-wdata/trthomas/abcd_spark_merge/polygenic_scores/sumstats/Psychiatric_Genomics_Consortium/OCD-PGC-2018'
ss <- read_tsv(ocd)
ss <- ss %>%
  filter(A1 %in% c('A', 'C', 'T', 'G') & A2 %in% c('A', 'C', 'T', 'G')) %>%
  mutate(beta = log(OR),
         se_beta = abs(beta / qnorm(P / 2))
  )
NCAS = 2688
NCON =  7037
ss$N = 4 / (1 / NCAS + 1 / NCON)
ss
pheno = '2018-PGC-OCD'
#
##
str_c("number of rows in the raw summary statistics file is: ", nrow(ss))
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
  filter(INFO >= 0.8) %>%
  select(
    chr = CHR,
    pos = BP,
    a1 = A1, # a1 is the effect allele
    a0 = A2,
    n_eff = N,
    beta = beta,
    beta_se = se_beta,
    p = P) %>%
  filter(a1 != a0) %>%
  drop_na() %>%
  inner_join(select(bim_snp, chr, pos)) %>%
  inner_join(select(map_ldref, chr, pos))
#
##
try(calculate_ldpred_pgs(sumstats = final_sumstats,  n_core = core_count, plink_rds = rds_path, sd_y = 0, pheno_name = pheno, build = 'hg19'))
rm(ss)
rm(sumstats)
rm(final_sumstats)
rm(pheno)
#
#
##########################################################
# EQ
# GWAS: https://www.nature.com/articles/mp2017122#Sec22
ss <- read_table('https://static-content.springer.com/esm/art%3A10.1038%2Fmp.2017.122/MediaObjects/41380_2018_BFmp2017122_MOESM84_ESM.txt') %>%
  arrange(P) %>%
  filter(Allele1 %in% tolower(c('A', 'C', 'T', 'G')) & Allele2 %in% tolower(c('A', 'C', 'T', 'G'))) %>%
  mutate(N = 88056) %>%
  mutate(Allele1 = toupper(Allele1),
         Allele2 = toupper(Allele2))
pheno = '2017-Warrier-MolecularPsychiatry-cognitive_empathy'
#
##
str_c("number of rows in the raw summary statistics file is: ", nrow(ss))
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
  select(
    chr = CHR,
    pos = BP,
    a1 = Allele1, # a1 is the effect allele
    a0 = Allele2,
    n_eff = N,
    beta = Effect,
    beta_se = StdErr,
    p = P) %>%
  filter(a1 != a0) %>%
  drop_na() %>%
  inner_join(select(bim_snp, chr, pos)) %>%
  inner_join(select(map_ldref, chr, pos))
#
##
try(calculate_ldpred_pgs(sumstats = final_sumstats,  n_core = core_count, plink_rds = rds_path, sd_y = 0, pheno_name = pheno, build = 'hg19'))
rm(ss)
rm(sumstats)
rm(final_sumstats)
rm(pheno)
#
#
###################################################################
# vocal pitch
sumstats = '/Dedicated/jmichaelson-sdata/gwas_summary_stats/deCODE/vocal_gisladottir_SciAdvances2023/reformatted_hg19_HapMap3_plus_deCode_VoiceSpeech_MedianF0_MeasureReading.csv'
ss = read_csv(sumstats)
#
#### calculate pgs ####
pheno = "vocal_pitch_median_f0"
#
str_c("number of rows in the raw summary statistics file is: ", nrow(ss))
head(ss)
#
### 3: give the output path
output_path <- str_c(my_wd, pheno, "/")
str_c("the final output path is: ", output_path)
dir.create(output_path)
setwd(output_path)
#
##### 4: modify/create/filter the summary statistics columns to have: chr, pos, a1 (effect allele), a0, beta, beta_se, p, n_eff
final_sumstats <- ss
##
final_sumstats$Neff = 12901
#
final_sumstats = final_sumstats %>%
  select(chr,
         pos,
         a1, # a1 is the effect allele
         a0,
         n_eff = Neff,
         beta,
         beta_se = se,
         p = P)
#
##
try(calculate_ldpred_pgs(sumstats = final_sumstats,  n_core = core_count, plink_rds = rds_path, sd_y = 0, pheno_name = pheno, build = 'hg19'))
rm(ss)
rm(sumstats)
rm(final_sumstats)
rm(pheno)
#
#
###################################################################
# executive functioning
ef <- '/Dedicated/jmichaelson-sdata/gwas_summary_stats/gwas_catalog/hatoum2023-BiologicalPsychiarty-executiveFunctioning/GCST90162547_buildGRCh37.tsv'
ss <- read_tsv(ef)
#
ss <- ss %>%
  mutate(N = 427037) %>%
  filter(effect_allele_frequency >= 0.01) %>%
  filter(INFO >= 0.8)
ss
pheno = 'hatoum2023-executive_functioning'
#
##
str_c("number of rows in the raw summary statistics file is: ", nrow(ss))
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
  select(
    chr = chromosome,
    pos = base_pair_location,
    a1 = effect_allele, # a1 is the effect allele
    a0 = other_allele,
    n_eff = N,
    beta = beta,
    beta_se = standard_error,
    p = p_value) %>%
  filter(a1 != a0) %>%
  drop_na() %>%
  inner_join(select(bim_snp, chr, pos)) %>%
  inner_join(select(map_ldref, chr, pos))
#
##
try(calculate_ldpred_pgs(sumstats = final_sumstats,  n_core = core_count, plink_rds = rds_path, sd_y = 0, pheno_name = pheno, build = 'hg19'))
rm(ss)
rm(sumstats)
rm(final_sumstats)
rm(pheno)
#
###################################################################
# Physical activity
ac <- '/Dedicated/jmichaelson-sdata/gwas_summary_stats/gwas_catalog/wang2022-NatureGenetics-physicalActivity/GCST90104341_buildGRCh37.tsv'
ss <- read_tsv(ac)
#
ss <- ss %>%
  filter(Freq1 >= 0.01) %>%
  mutate(Allele1 = toupper(Allele1),
         Allele2 = toupper(Allele2)) %>%
  filter(nchar(Allele1) == 1) %>%
  filter(nchar(Allele2) == 1) %>%
  filter(Allele1 %in% c("A", "C", "T", "G") & Allele2 %in% c("A", "C", "T", "G")) %>%
  mutate(N = 608595) %>%
  filter(N >= 0.7 * max(SampleSize))
ss
pheno = 'wang2022-physical_activity_during_leisure_time'
#
##
str_c("number of rows in the raw summary statistics file is: ", nrow(ss))
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
  select(
    chr = chromosome,
    pos = base_pair_location,
    a1 = Allele1, # a1 is the effect allele
    a0 = Allele2,
    n_eff = SampleSize,
    beta = BETA,
    beta_se = SE,
    p = p_value) %>%
  filter(a1 != a0) %>%
  drop_na() %>%
  inner_join(select(bim_snp, chr, pos)) %>%
  inner_join(select(map_ldref, chr, pos))
#
##
try(calculate_ldpred_pgs(sumstats = final_sumstats,  n_core = core_count, plink_rds = rds_path, sd_y = 0, pheno_name = pheno, build = 'hg19'))
rm(ss)
rm(sumstats)
rm(final_sumstats)
rm(pheno)
#
#
#
#
#
#
#
#
## Alzheimer's
f = '/Dedicated/jmichaelson-sdata/gwas_summary_stats/alzheimer/kunke2019-NatureGenetics-alzheimers/Kunkle_etal_Stage1_results.txt?file=1'
pheno = 'kunke2019_alzheimers'
message('Working on ', pheno)
#
ss <- read_table(f, show_col_types = F)
#
ncas = 21982
ncon = 41944
ss$N <- 4 / (1 / ncas + 1 / ncon)
ss <- ss %>%
  filter(nchar(Effect_allele) == 1) %>%
  filter(nchar(Non_Effect_allele) == 1)
#
##
str_c("number of rows in the raw summary statistics file is: ", nrow(ss))
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
  select(
    chr = Chromosome,
    pos = Position,
    a1 = Effect_allele, # a1 is the effect allele
    a0 = Non_Effect_allele,
    n_eff = N,
    beta = Beta,
    beta_se = SE,
    p = Pvalue) %>%
  filter(a1 != a0) %>%
  drop_na() %>%
  inner_join(select(bim_snp, chr, pos)) %>%
  inner_join(select(map_ldref, chr, pos))
#
##
try(calculate_ldpred_pgs(sumstats = final_sumstats, n_core = core_count,  plink_rds = rds_path, sd_y = 0, pheno_name = pheno, build = 'hg19'))
rm(ss)
rm(final_sumstats)
rm(pheno)
#
#########################################################
## cognitive resilience
f = '/Dedicated/jmichaelson-sdata/gwas_summary_stats/UKBIOBANK/fitzgerald2022-genes-cognitive_resilience/all.res.txt'
pheno = 'fitzgerald2022_cognitive_resilience'
message('Working on ', pheno)
#
ss <- read_table(f, show_col_types = F)
#
ss$N <- 330097
ss <- ss %>%
  filter(nchar(A1) == 1) %>%
  filter(nchar(A2) == 1) %>%
  filter(MAF >= 0.01)
#
##
str_c("number of rows in the raw summary statistics file is: ", nrow(ss))
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
  select(
    chr = CHR,
    pos = POS,
    a1 = A1, # a1 is the effect allele
    a0 = A2,
    n_eff = N,
    beta = Beta,
    beta_se = SE,
    p = P) %>%
  filter(a1 != a0) %>%
  drop_na() %>%
  inner_join(select(bim_snp, chr, pos)) %>%
  inner_join(select(map_ldref, chr, pos))
#
##
try(calculate_ldpred_pgs(sumstats = final_sumstats, n_core = core_count,  plink_rds = rds_path, sd_y = 0, pheno_name = pheno, build = 'hg19'))
rm(ss)
rm(final_sumstats)
rm(pheno)
#
#
#########################################################
## lean muscle mass
f = '/Dedicated/jmichaelson-sdata/gwas_summary_stats/gwas_catalog/pei2020-CommunicationBiology-appendicularLeanMass/GCST90000025_GRCh37.tsv'
pheno = 'pei2020_appendicularLeanMass'
message('Working on ', pheno)
#
ss <- read_table(f, show_col_types = F)
#
ss$N <- 450243
ss <- ss %>%
  filter(nchar(effect_allele) == 1) %>%
  filter(nchar(other_allele) == 1) %>%
  filter(effect_allele_frequency >= 0.01) %>%
  mutate(effect_allele = toupper(effect_allele),
         other_allele = toupper(other_allele)
  )
#
##
str_c("number of rows in the raw summary statistics file is: ", nrow(ss))
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
  select(
    chr = chromosome,
    pos = base_pair_location,
    a1 = effect_allele, # a1 is the effect allele
    a0 = other_allele,
    n_eff = N,
    beta = beta,
    beta_se = standard_error,
    p = p_value) %>%
  filter(a1 != a0) %>%
  drop_na() %>%
  inner_join(select(bim_snp, chr, pos)) %>%
  inner_join(select(map_ldref, chr, pos))
#
##
try(calculate_ldpred_pgs(sumstats = final_sumstats,  n_core = core_count, plink_rds = rds_path, sd_y = 0, pheno_name = pheno, build = 'hg19'))
rm(ss)
rm(final_sumstats)
rm(pheno)
#
#
#####################################################################
## epigenetic aging
files = list.files('/Dedicated/jmichaelson-sdata/gwas_summary_stats/epigenetic_age/mccartney2021-genomeBiology-epigenetic_clocks',
                   pattern = 'txt$', full.names = T)
for (f in files) {
  pheno = basename(f)
  pheno = str_split(pheno, pattern = '_', simplify = T)[,1]
  pheno = str_c('epigenetic_', pheno)
  message('Working on ', pheno)
#
  ss <- read_table(f, show_col_types = F)
#
  ss <- ss %>%
    filter(nchar(A1) == 1) %>%
    filter(nchar(A2) == 1) %>%
    filter(Freq1 >= 0.01) %>%
    mutate(A1 = toupper(A1),
           A2 = toupper(A2)
    ) %>%
    filter(N >= .5 * max(N))
#
  ##
  str_c("number of rows in the raw summary statistics file is: ", nrow(ss))
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
    select(
      chr = chr,
      pos = bp,
      a1 = A1, # a1 is the effect allele
      a0 = A2,
      n_eff = N,
      beta = Effect,
      beta_se = SE,
      p = P) %>%
    filter(a1 != a0) %>%
    drop_na() %>%
    inner_join(select(bim_snp, chr, pos)) %>%
    inner_join(select(map_ldref, chr, pos))
#
  ##
  try(calculate_ldpred_pgs(sumstats = final_sumstats,  n_core = core_count, plink_rds = rds_path, sd_y = 0, pheno_name = pheno, build = 'hg19'))
  rm(ss)
  rm(final_sumstats)
  rm(pheno)
}
#
#
## Alzheimer's
f = '/Dedicated/jmichaelson-wdata/lcasten/misc/gwas_summary_stats/ctg_sumstats_posthuma/alzheimer/PGCALZ2ExcludingUKBand23andME_METALInverseVariance_MetaAnalysis.txt'
pheno = 'PGC2_exclUKBB_excl23andMe_alzheimers'
message('Working on ', pheno)
#
ss <- read_table(f, show_col_types = F)
#
ss <- ss %>%
  filter(nchar(effect_allele) == 1) %>%
  filter(nchar(other_allele) == 1) %>%
  filter(effect_allele_frequency >= 0.01) %>%
  # filter(Build == 'GRCh37') %>%
  filter(N >= 0.7 * max(Neffective))
#
##
str_c("number of rows in the raw summary statistics file is: ", nrow(ss))
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
  select(
    chr = chromosome,
    pos = base_pair_location,
    a1 = effect_allele, # a1 is the effect allele
    a0 = other_allele,
    n_eff = Neffective,
    beta = beta,
    beta_se = standard_error,
    p = p_value) %>%
  filter(a1 != a0) %>%
  drop_na() %>%
  inner_join(select(bim_snp, chr, pos)) %>%
  inner_join(select(map_ldref, chr, pos))
#
##
try(calculate_ldpred_pgs(sumstats = final_sumstats, n_core = core_count,  plink_rds = rds_path, sd_y = 0, pheno_name = pheno, build = 'hg19'))
rm(ss)
rm(final_sumstats)
rm(pheno)
#
#
##
#
#### get files ####
brain_files = list.files('/Dedicated/jmichaelson-sdata/gwas_summary_stats/brain_imaging_Warrier_NatureGenetics_2023',
                         pattern = '_meta.txt', full.names = T)
files <- c(brain_files)
#
files = files[str_detect(files, pattern = '_Warrier_')]
#
phenos <- c(basename(files))
phenos = case_when(str_detect(phenos, pattern = 'CT_meta') == T ~ 'cortical_thickness.Warrier2023',
                   str_detect(phenos, pattern = 'FA_meta') == T ~ 'fractional_anisotropy.Warrier2023',
                   str_detect(phenos, pattern = 'Foldingindex_meta') == T ~ 'folding_index.Warrier2023',
                   str_detect(phenos, pattern = 'Gaussiancurvature_meta') == T ~ 'gaussian_curvature.Warrier2023',
                   str_detect(phenos, pattern = 'ICVF_meta') == T ~ 'intracellular_volume_fraction.Warrier2023',
                   str_detect(phenos, pattern = 'Intrinsic_meta') == T ~ 'intrinsic_curvature_index.Warrier2023',
                   str_detect(phenos, pattern = 'ISOVF_meta') == T ~ 'isotropic_volume_fraction.Warrier2023',
                   str_detect(phenos, pattern = 'LGI_meta') == T ~ 'local_gyrification_index.Warrier2023',
                   str_detect(phenos, pattern = 'MD_meta') == T ~ 'mean_diffusivity.Warrier2023',
                   str_detect(phenos, pattern = 'Meancurvature_meta') == T ~ 'mean_curvature.Warrier2023',
                   str_detect(phenos, pattern = 'OD_meta') == T ~ 'orientation_diffusion_index.Warrier2023',
                   str_detect(phenos, pattern = 'SA_meta') == T ~ 'surface_area.Warrier2023',
                   str_detect(phenos, pattern = 'Volume_meta') == T ~ 'volume.Warrier2023',
                   TRUE ~ phenos
)
#
#
#################################################
#### loop over sumstats files to compute PGS ####
message('Working on polygenic score calculations now...')
message('Files: ')
message(str_c(files, collapse = ', '))
source('/Dedicated/jmichaelson-wdata/lcasten/functions/LDPred2_argon_hm3plus.R')
for (sumstats in files) {
  message('')
  message('')
  message('------------------------------------------')
  i <- which(sumstats == files)
  pheno = phenos[i]
  message('Working on ', pheno, ' (', i, ' / ', length(files), ')')
#
  #### read in sumstats
  ss = read_table(sumstats)
#
  #### calculate pgs ####
  str_c("number of rows in the raw summary statistics file is: ", nrow(ss))
#
  ### 3: give the output path
  output_path <- str_c(my_wd, pheno, "/")
  str_c("the final output path is: ", output_path)
  dir.create(output_path)
  setwd(output_path)
#
  ### normalize column names
  if (str_detect(pheno, pattern = 'log10') == T) {
    ss <- ss %>%
      rename(chr = CHR, pos = POS, a1 = A1, a0 = A2, n_eff = N, beta = BETA, beta_se = SE, p = P)
  }
#
  if (str_detect(pheno, pattern = 'log10') == F) {
    ss <- ss %>%
      #  N =1, sample size = 31,797 (UKB only). If N = 2, sample size = 36,663
      mutate(N = case_when(N == 1 ~ 31797,
                           N == 2 ~ 36663)
      ) %>%
      rename(chr = CHR, pos = BP, a1 = A1, a0 = A2, n_eff = N, beta = BETA, beta_se = SE, p = P)
  }
#
  ##
  final_sumstats = ss %>%
    select(chr,
           pos,
           a1, # a1 is the effect allele
           a0,
           n_eff,
           beta,
           beta_se,
           p) %>%
    filter(a1 != a0) %>%
    drop_na() %>%
    inner_join(select(bim_snp, chr, pos)) %>%
    inner_join(select(map_ldref, chr, pos))
#
#
  ## compute PGS
  try(calculate_ldpred_pgs(sumstats = final_sumstats, n_core = core_count, plink_rds = rds_path, sd_y = 0, pheno_name = pheno, build = 'hg19'))
#
  rm(ss)
  rm(final_sumstats)
  rm(sumstats)
  rm(pheno)
}
#
#
# #####################################################
# #### new PGS ####
# #####################################################
## new adhd
adhd <- '/Dedicated/jmichaelson-sdata/gwas_summary_stats/PGC/pgc-adhd2022/ADHD2022_iPSYCH_deCODE_PGC.meta.gz'
ss <- read_delim(adhd, delim = ' ')
ss
ss <- ss %>%
  mutate(CHR = as.numeric(CHR)) %>%
  drop_na(CHR) %>%
  mutate(beta = log(OR),
         se_beta = abs(beta / qnorm(P / 2))) %>%
  filter(FRQ_A_38691 >= 0.001 & FRQ_U_186843 >= 0.001) %>%
  filter(nchar(A1) == 1) %>%
  filter(nchar(A2) == 1) %>%
  filter(A1 %in% c("A", "C", "T", "G") & A2 %in% c("A", "C", "T", "G")) %>%
  mutate(N = 4 / (1 / Nca + 1 / Nco)) %>%
  filter(N >= 0.5 * max(N))
ss
pheno = 'PGC2022-ADHD'
#
##
str_c("number of rows in the raw summary statistics file is: ", nrow(ss))
#
### 3: give the output path
output_path <- str_c(my_wd, pheno, "/")
str_c("the final output path is: ", output_path)
system(str_c('rm -r ', output_path))
dir.create(output_path)
setwd(output_path)
#
##### 4: modify/create/filter the summary statistics columns to have: chr, pos, a1 (effect a), a0, beta, beta_se, p, n_eff
final_sumstats <- ss %>%
  select(
    chr = CHR,
    pos = BP,
    a1 = A1, # a1 is the effect a
    a0 = A2,
    n_eff = N,
    beta = beta,
    beta_se = se_beta,
    p = P) %>%
  filter(a1 != a0) %>%
  drop_na() %>%
  inner_join(select(bim_snp, chr, pos)) %>%
  inner_join(select(map_ldref, chr, pos))
#
##
try(calculate_ldpred_pgs(sumstats = final_sumstats,  n_core = core_count, plink_rds = rds_path, sd_y = 0, pheno_name = pheno, build = 'hg19'))
rm(ss)
rm(sumstats)
rm(final_sumstats)
rm(pheno)
#


##########################

#####################
##### tobacco use
####################
f <- '/Dedicated/jmichaelson-sdata/gwas_summary_stats/UKBIOBANK/Watanabe2019_NatureGenetics_GWAS_atlas/tobacco_use/f.22506.0.0_res.EUR.sumstats.MACfilt.txt.gz'
ss <- read_delim(f, delim = '\t') 
ss
ss <- ss %>%
  mutate(CHR = as.numeric(CHR)) %>%
  drop_na(CHR)

ss <- ss %>%
  filter(nchar(A1) == 1) %>%
  filter(nchar(A2) == 1) %>%
  filter(A1 %in% c("A", "C", "T", "G") & A2 %in% c("A", "C", "T", "G")) %>%
  # mutate(N = 4 / (1 / NCAS + 1 / NCON),
  #        N = round(N, digits = 0)) %>%
  filter(MAF >= 0.001)
ss
pheno = 'GWASatlas2019-tobacco_use'

##
str_c("number of rows in the raw summary statistics file is: ", nrow(ss))

### 3: give the output path
output_path <- str_c(my_wd, pheno, "/")
str_c("the final output path is: ", output_path)
system(str_c('rm -r ', output_path))
dir.create(output_path)
setwd(output_path)

##### 4: modify/create/filter the summary statistics columns to have: chr, pos, a1 (effect a), a0, beta, beta_se, p, n_eff
final_sumstats <- ss %>%
  select(
    chr = CHR,
    pos = BP,
    a1 = A1, # a1 is the effect a
    a0 = A2,
    n_eff = NMISS,
    beta = BETA,
    beta_se = SE,
    p = P) %>%
  filter(a1 != a0) %>%
  drop_na() %>%
  inner_join(select(bim_snp, chr, pos)) %>%
  inner_join(select(map_ldref, chr, pos))

##
try(calculate_ldpred_pgs(sumstats = final_sumstats,  n_core = core_count, plink_rds = rds_path, sd_y = 0, pheno_name = pheno, build = 'hg19'))
rm(ss)
rm(sumstats)
rm(final_sumstats)
rm(pheno)

##########
## NHB
f = '/Dedicated/jmichaelson-wdata/trthomas/abcd_spark_merge/polygenic_scores/sumstats/UK_Biobank/non_heterosexual_behavior-UKB-2019'
ss <- read_delim(f, delim = '\t') 
ss
ss <- ss %>%
  mutate(CHR = as.numeric(CHR)) %>%
  drop_na(CHR)
NCON = 477522
NCAS = 26827
NEFF = 4 / (1 / NCAS + 1 / NCON)
NEFF = round(NEFF, digits = 0)
NEFF
ss <- ss %>%
  filter(nchar(A1) == 1) %>%
  filter(nchar(A2) == 1) %>%
  filter(A1 %in% c("A", "C", "T", "G") & A2 %in% c("A", "C", "T", "G")) %>%
  mutate(N = NEFF
         ) %>%
  filter(A1FREQ >= 0.001) %>%
  filter(INFO >= 0.8)
ss
pheno = 'UKBB2019-non_heterosexual_behavior'

##
str_c("number of rows in the raw summary statistics file is: ", nrow(ss))

### 3: give the output path
output_path <- str_c(my_wd, pheno, "/")
str_c("the final output path is: ", output_path)
system(str_c('rm -r ', output_path))
dir.create(output_path)
setwd(output_path)

##### 4: modify/create/filter the summary statistics columns to have: chr, pos, a1 (effect a), a0, beta, beta_se, p, n_eff
final_sumstats <- ss %>%
  select(
    chr = CHR,
    pos = BP,
    a1 = A1, # a1 is the effect a
    a0 = A2,
    n_eff = N,
    beta = BETA,
    beta_se = SE,
    p = P_BOLT_LMM_INF) %>%
  filter(a1 != a0) %>%
  drop_na() %>%
  inner_join(select(bim_snp, chr, pos)) %>%
  inner_join(select(map_ldref, chr, pos))

##
try(calculate_ldpred_pgs(sumstats = final_sumstats,  n_core = core_count, plink_rds = rds_path, sd_y = 0, pheno_name = pheno, build = 'hg19'))
rm(ss)
rm(sumstats)
rm(final_sumstats)
rm(pheno)


# #####################################################
# #### new PGS ####
# #####################################################
## new adhd
adhd <- '/Dedicated/jmichaelson-sdata/gwas_summary_stats/PGC/pgc-adhd2022/ADHD2022_iPSYCH_deCODE_PGC.meta.gz'
ss <- read_delim(adhd, delim = ' ')
ss
ss <- ss %>%
  mutate(CHR = as.numeric(CHR)) %>%
  drop_na(CHR) %>%
  mutate(beta = log(OR),
         se_beta = abs(beta / qnorm(P / 2))) %>%
  filter(FRQ_A_38691 >= 0.001 & FRQ_U_186843 >= 0.001) %>%
  filter(nchar(A1) == 1) %>%
  filter(nchar(A2) == 1) %>%
  filter(A1 %in% c("A", "C", "T", "G") & A2 %in% c("A", "C", "T", "G")) %>%
  mutate(N = 4 / (1 / Nca + 1 / Nco)) %>%
  filter(N >= 0.5 * max(N))
ss
pheno = 'PGC2022-ADHD'
#
##
str_c("number of rows in the raw summary statistics file is: ", nrow(ss))
#
### 3: give the output path
output_path <- str_c(my_wd, pheno, "/")
str_c("the final output path is: ", output_path)
system(str_c('rm -r ', output_path))
dir.create(output_path)
setwd(output_path)
#
##### 4: modify/create/filter the summary statistics columns to have: chr, pos, a1 (effect a), a0, beta, beta_se, p, n_eff
final_sumstats <- ss %>%
  select(
    chr = CHR,
    pos = BP,
    a1 = A1, # a1 is the effect a
    a0 = A2,
    n_eff = N,
    beta = beta,
    beta_se = se_beta,
    p = P) %>%
  filter(a1 != a0) %>%
  drop_na() %>%
  inner_join(select(bim_snp, chr, pos)) %>%
  inner_join(select(map_ldref, chr, pos))
#
##
try(calculate_ldpred_pgs(sumstats = final_sumstats,  n_core = core_count, plink_rds = rds_path, sd_y = 0, pheno_name = pheno, build = 'hg19'))
rm(ss)
rm(sumstats)
rm(final_sumstats)
rm(pheno)
#
## tourette
ts <- '/Dedicated/jmichaelson-sdata/gwas_summary_stats/PGC/pgc-ts2019/TS_Oct2018.gz'
ss <- read_delim(ts, delim = ' ')
ss
Nca = 4819
Nco = 9488
NEFF = 4 / (1 / Nca + 1 / Nco)
NEFF = round(NEFF, digits = 0)
ss <- ss %>%
  mutate(CHR = as.numeric(CHR)) %>%
  drop_na(CHR) %>%
  mutate(beta = log(OR),
         se_beta = abs(beta / qnorm(P / 2))) %>%
  filter(nchar(A1) == 1) %>%
  filter(nchar(A2) == 1) %>%
  filter(A1 %in% c("A", "C", "T", "G") & A2 %in% c("A", "C", "T", "G")) %>%
  mutate(N = NEFF)
ss
pheno = 'PGC2019-TouretteSyndrome'
#
##
str_c("number of rows in the raw summary statistics file is: ", nrow(ss))
#
### 3: give the output path
output_path <- str_c(my_wd, pheno, "/")
str_c("the final output path is: ", output_path)
system(str_c('rm -r ', output_path))
dir.create(output_path)
setwd(output_path)
#
##### 4: modify/create/filter the summary statistics columns to have: chr, pos, a1 (effect a), a0, beta, beta_se, p, n_eff
final_sumstats <- ss %>%
  select(
    chr = CHR,
    pos = BP,
    a1 = A1, # a1 is the effect a
    a0 = A2,
    n_eff = N,
    beta = beta,
    beta_se = se_beta,
    p = P) %>%
  filter(a1 != a0) %>%
  drop_na() %>%
  inner_join(select(bim_snp, chr, pos)) %>%
  inner_join(select(map_ldref, chr, pos))
#
##
try(calculate_ldpred_pgs(sumstats = final_sumstats,  n_core = core_count, plink_rds = rds_path, sd_y = 0, pheno_name = pheno, build = 'hg19'))
rm(ss)
rm(sumstats)
rm(final_sumstats)
rm(pheno)

##########################
## hoarding
hd <- '/Dedicated/jmichaelson-sdata/gwas_summary_stats/PGC/pgc-hoarding2022/39399566'
vcf = read_lines(hd)
to_skip = max(which(str_detect(vcf, pattern = '^#')))
ss <- read_delim(hd, skip = to_skip, delim = '\t')
ss <- ss %>%
  mutate(CHROM = as.numeric(CHROM)) %>%
  drop_na(CHROM)
#
ss <- ss %>%
  filter(nchar(A1) == 1) %>%
  filter(nchar(A2) == 1) %>%
  filter(A1 %in% c("A", "C", "T", "G") & A2 %in% c("A", "C", "T", "G")) %>%
  mutate(N = 4 / (1 / NCAS + 1 / NCON),
         N = round(N, digits = 0)) %>%
  filter(IMPINFO >= 0.8) %>%
  filter(FCAS >= 0.001 & FCON >= 0.001)
ss
pheno = 'PGC2022-Hoarding'
#
##
str_c("number of rows in the raw summary statistics file is: ", nrow(ss))
#
### 3: give the output path
output_path <- str_c(my_wd, pheno, "/")
str_c("the final output path is: ", output_path)
system(str_c('rm -r ', output_path))
dir.create(output_path)
setwd(output_path)
#
##### 4: modify/create/filter the summary statistics columns to have: chr, pos, a1 (effect a), a0, beta, beta_se, p, n_eff
final_sumstats <- ss %>%
  select(
    chr = CHROM,
    pos = POS,
    a1 = A1, # a1 is the effect a
    a0 = A2,
    n_eff = N,
    beta = BETA,
    beta_se = SE,
    p = PVAL) %>%
  filter(a1 != a0) %>%
  drop_na() %>%
  inner_join(select(bim_snp, chr, pos)) %>%
  inner_join(select(map_ldref, chr, pos))

##
try(calculate_ldpred_pgs(sumstats = final_sumstats,  n_core = core_count, plink_rds = rds_path, sd_y = 0, pheno_name = pheno, build = 'hg19'))
rm(ss)
rm(sumstats)
rm(final_sumstats)
rm(pheno)


##################################
#### PGI scores ####
##################################
d = '/Dedicated/jmichaelson-wdata/trthomas/abcd_spark_merge/polygenic_scores/sumstats/Social_Science_Genetic_Association_Consortium/Polygenic_Index_Repository'
files = list.files(path = d,
                   recursive = TRUE, pattern = 'gwide')
files

for (f in files) {
  pheno = str_c('PGI-', basename(f))
  message('Working on: ', pheno)
  
  f = str_c(d, '/', f)
  ss <- read_delim(f, delim = '\t') 
  ss
  ss <- ss %>%
    mutate(CHR = as.numeric(CHR)) %>%
    drop_na(CHR)
  
  ss <- ss %>%
    filter(nchar(EFFECT_ALLELE) == 1) %>%
    filter(nchar(OTHER_ALLELE) == 1) %>%
    filter(EFFECT_ALLELE %in% c("A", "C", "T", "G") & OTHER_ALLELE %in% c("A", "C", "T", "G"))
  ss
  
  ##
  str_c("number of rows in the raw summary statistics file is: ", nrow(ss))
  
  ### 3: give the output path
  output_path <- str_c(my_wd, pheno, "/")
  str_c("the final output path is: ", output_path)
  system(str_c('rm -r ', output_path))
  dir.create(output_path)
  setwd(output_path)
  
  ##### 4: modify/create/filter the summary statistics columns to have: chr, pos, a1 (effect a), a0, beta, beta_se, p, n_eff
  final_sumstats <- ss %>%
    select(
      chr = CHR,
      pos = BP,
      a1 = EFFECT_ALLELE, # a1 is the effect a
      a0 = OTHER_ALLELE,
      n_eff = GWASeqN,
      beta = BETA,
      beta_se = SE,
      p = PVALUE) %>%
    filter(a1 != a0) %>%
    drop_na() %>%
    inner_join(select(bim_snp, chr, pos)) %>%
    inner_join(select(map_ldref, chr, pos))
  
  ##
  try(calculate_ldpred_pgs(sumstats = final_sumstats,  n_core = core_count, plink_rds = rds_path, sd_y = 0, pheno_name = pheno, build = 'hg19'))
  rm(ss)
  rm(sumstats)
  rm(final_sumstats)
  rm(pheno)
}


##################################
#### diet stuff ####
##################################
d = '/Dedicated/jmichaelson-sdata/gwas_summary_stats/ssgac/2020'
files = list.files(path = d,
                   recursive = TRUE, pattern = 'txt.gz$')
files

for (f in files) {
  pheno = str_remove_all(f, pattern = '.txt.gz')
  message('Working on: ', pheno)
  
  f = str_c(d, '/', f)
  ss <- read_delim(f, delim = '\t') 
  ss
  ss <- ss %>%
    mutate(CHR = as.numeric(CHR)) %>%
    drop_na(CHR)
  
  ss <- ss %>%
    filter(nchar(A1) == 1) %>%
    filter(nchar(A2) == 1) %>%
    mutate(A1 = toupper(A1),
           A2 = toupper(A2)) %>%
    filter(A1 %in% c("A", "C", "T", "G") & A2 %in% c("A", "C", "T", "G"))
  ss
  
  ##
  str_c("number of rows in the raw summary statistics file is: ", nrow(ss))
  
  ### 3: give the output path
  output_path <- str_c(my_wd, pheno, "/")
  str_c("the final output path is: ", output_path)
  system(str_c('rm -r ', output_path))
  dir.create(output_path)
  setwd(output_path)
  
  ##### 4: modify/create/filter the summary statistics columns to have: chr, pos, a1 (effect a), a0, beta, beta_se, p, n_eff
  final_sumstats <- ss %>%
    select(
      chr = CHR,
      pos = POS,
      a1 = A1, # a1 is the effect a
      a0 = A2,
      n_eff = N,
      beta = Beta,
      beta_se = SE,
      p = Pval) %>%
    filter(a1 != a0) %>%
    drop_na() %>%
    inner_join(select(bim_snp, chr, pos)) %>%
    inner_join(select(map_ldref, chr, pos))
  
  ##
  try(calculate_ldpred_pgs(sumstats = final_sumstats,  n_core = core_count, plink_rds = rds_path, sd_y = 0, pheno_name = pheno, build = 'hg19'))
  rm(ss)
  rm(sumstats)
  rm(final_sumstats)
  rm(pheno)
}

###########################
#### EA gwas by subtraction
## non-cognitive 
pheno = 'educational_attainment_non_cognitive-Demange2021'
message('Working on: ', pheno)

f = c('/Dedicated/jmichaelson-sdata/gwas_summary_stats/gwas_catalog/demange2021-NatureGenetics-eduAttainmentGWASBySubstraction/EA_non_cognitive/GCST90011874_buildGRCh37.tsv.gz')
ss <- read_delim(f, delim = '\t') 
ss
ss <- ss %>%
  mutate(chromosome = as.numeric(chromosome)) %>%
  drop_na(chromosome)

ss <- ss %>%
  filter(nchar(effect_allele) == 1) %>%
  filter(nchar(other_allele) == 1) %>%
  filter(effect_allele %in% c("A", "C", "T", "G") & other_allele %in% c("A", "C", "T", "G")) %>%
  filter(MAF >= 0.001)
ss$N = 510795

##
str_c("number of rows in the raw summary statistics file is: ", nrow(ss))

### 3: give the output path
output_path <- str_c(my_wd, pheno, "/")
str_c("the final output path is: ", output_path)
system(str_c('rm -r ', output_path))
dir.create(output_path)
setwd(output_path)

##### 4: modify/create/filter the summary statistics columns to have: chr, pos, a1 (effect a), a0, beta, beta_se, p, n_eff
final_sumstats <- ss %>%
  select(
    chr = chromosome,
    pos = base_pair_location,
    a1 = effect_allele, # a1 is the effect a
    a0 = other_allele,
    n_eff = N,
    beta = est,
    beta_se = standard_error,
    p = p_value) %>%
  filter(a1 != a0) %>%
  drop_na() %>%
  inner_join(select(bim_snp, chr, pos)) %>%
  inner_join(select(map_ldref, chr, pos))

##
try(calculate_ldpred_pgs(sumstats = final_sumstats,  n_core = core_count, plink_rds = rds_path, sd_y = 0, pheno_name = pheno, build = 'hg19'))
rm(ss)
rm(sumstats)
rm(final_sumstats)
rm(pheno)


## cognitive 
pheno = 'educational_attainment_cognitive-Demange2021'
message('Working on: ', pheno)

f = c('/Dedicated/jmichaelson-sdata/gwas_summary_stats/gwas_catalog/demange2021-NatureGenetics-eduAttainmentGWASBySubstraction/EA_cognitive/GCST90011875_buildGRCh37.tsv.gz')
ss <- read_delim(f, delim = '\t') 
ss
ss <- ss %>%
  mutate(chromosome = as.numeric(chromosome)) %>%
  drop_na(chromosome)

ss <- ss %>%
  filter(nchar(effect_allele) == 1) %>%
  filter(nchar(other_allele) == 1) %>%
  filter(effect_allele %in% c("A", "C", "T", "G") & other_allele %in% c("A", "C", "T", "G")) %>%
  filter(MAF >= 0.001)
ss$N = 257700

##
str_c("number of rows in the raw summary statistics file is: ", nrow(ss))

### 3: give the output path
output_path <- str_c(my_wd, pheno, "/")
str_c("the final output path is: ", output_path)
system(str_c('rm -r ', output_path))
dir.create(output_path)
setwd(output_path)

##### 4: modify/create/filter the summary statistics columns to have: chr, pos, a1 (effect a), a0, beta, beta_se, p, n_eff
final_sumstats <- ss %>%
  select(
    chr = chromosome,
    pos = base_pair_location,
    a1 = effect_allele, # a1 is the effect a
    a0 = other_allele,
    n_eff = N,
    beta = est,
    beta_se = standard_error,
    p = p_value) %>%
  filter(a1 != a0) %>%
  drop_na() %>%
  inner_join(select(bim_snp, chr, pos)) %>%
  inner_join(select(map_ldref, chr, pos))

##
try(calculate_ldpred_pgs(sumstats = final_sumstats,  n_core = core_count, plink_rds = rds_path, sd_y = 0, pheno_name = pheno, build = 'hg19'))
rm(ss)
rm(sumstats)
rm(final_sumstats)
rm(pheno)

## done