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
sumstats_files = list.files(path = '/Dedicated/jmichaelson-sdata/gwas_summary_stats/ENIGMA_brain_imaging/subcortical_unrestricted_use', 
                            pattern = '[.]gz$', full.names = T)

## file paths
ss_map = data.frame(filepath = sumstats_files, 
                            pheno = basename(sumstats_files)
                            ) %>%
  mutate(pheno = str_replace_all(pheno, pattern = '_eur_z_ldsc_unrestricted_NG05SEP19.gz', replacement = '_volume'))

#### calculate pgs ####
for (i in 2:nrow(ss_map)) {
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
    ss <- read_delim(as.character(ss_map$filepath[i]), delim = ' ', guess_max = 8000000)
    
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
    #   mutate(chromosome = str_split(MarkerName, pattern = '[:]', simplify = TRUE)[,1],
    #          position = str_split(MarkerName, pattern = '[:]', simplify = TRUE)[,2]) %>%
    #   mutate(chromosome = as.numeric(chromosome),
    #          position = as.numeric(position)) %>%
    #   drop_na(chromosome) %>%
      mutate(A2 = toupper(A2),
             A1 = toupper(A1)) %>%
      filter(nchar(A1) == 1 & nchar(A2) == 1) %>% 
      filter(A1 %in% c('A', 'C', 'T', 'G') & A2 %in% c('A', 'C', 'T', 'G')) %>%
      mutate(Beta = Zscore / sqrt(2 * FreqA1 * (1 - FreqA1) * (N + Zscore^2)),
             SE = 1 / sqrt(2 * FreqA1 * (1 - FreqA1) * (N + Zscore^2)))
    
    ##
    final_sumstats$Neff = final_sumstats$N
    
    message('after GWAS sumstats QC and before overlapping w/ HapMap3 there are ', nrow(final_sumstats), ' SNPs')
    message('(dropped ', nrow(ss) - nrow(final_sumstats), ' variants)')
    
    final_sumstats = final_sumstats %>% 
      select(rsid, 
            #  chr = chromosome,
            #  pos = position,
             a1 = A1, # a1 is the effect allele
             a0 = A2,
             n_eff = Neff,
             beta = Beta,
             beta_se = SE,
             p = P) %>%
    drop_na() %>%
    inner_join(select(map_ldref, rsid, chr, pos)) %>%
    relocate(chr, pos, .after = rsid) %>%
    inner_join(select(bim_snp, chr, pos))
    
    ## compute PGS
    try(calculate_ldpred_pgs(sumstats = final_sumstats, n_core = core_count,  plink_rds = rds_path, sd_y = 0, pheno_name = pheno, build = 'hg19'))
  }
}



#####################################
#####################################
#####################################
## personality
sumstats_files = list.files(path = '/Dedicated/jmichaelson-sdata/gwas_summary_stats/personality/BIG5-Gupta2024-NatHumBeh', 
                            pattern = '_sumstat_file$', full.names = T)

## file paths
ss_map = data.frame(filepath = sumstats_files, 
                            pheno = basename(sumstats_files)
                            ) %>%
  mutate(pheno = str_replace_all(pheno, pattern = '_MVP_GPC1_sumstat_file|_MVP_UKB_sumstat_file', replacement = '_BIG5'))

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
    ss <- read_delim(as.character(ss_map$filepath[i]), delim = '\t')
    
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
      mutate(CHR = as.numeric(CHR),
             BP = as.numeric(BP)) %>%
      drop_na(CHR) %>%
      filter(nchar(A1) == 1 & nchar(A2) == 1) %>% 
      filter(A1 %in% c('A', 'C', 'T', 'G') & A2 %in% c('A', 'C', 'T', 'G'))
    
    ##
    final_sumstats$Neff = final_sumstats$N
    
    message('after GWAS sumstats QC and before overlapping w/ HapMap3 there are ', nrow(final_sumstats), ' SNPs')
    message('(dropped ', nrow(ss) - nrow(final_sumstats), ' variants)')
    
    final_sumstats = final_sumstats %>% 
      select(rsid = SNP, 
             chr = CHR,
             pos = BP,
             a1 = A1, # a1 is the effect allele
             a0 = A2,
             n_eff = Neff,
             beta = BETA,
             beta_se = SE,
             p = P) %>%
    drop_na() %>%
    inner_join(select(bim_snp, chr, pos)) %>%
    inner_join(select(map_ldref, rsid, chr, pos))
    
    ## compute PGS
    try(calculate_ldpred_pgs(sumstats = final_sumstats, n_core = core_count,  plink_rds = rds_path, sd_y = 0, pheno_name = pheno, build = 'hg19'))
  }
}