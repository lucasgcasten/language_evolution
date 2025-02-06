## =========================================================================
## this script will take the raw PGS and adjust them for population stratification
## specifically, it will adjust for the main effects of the first 5 genetic PCs based on 1000 Genomes Europeans
## it will output a csv with all of the adjusted PGS which will then be ready for ES-PGS analysis
## =========================================================================

################################################
## set up: load packages and read in data
################################################

## load packages
library(tidyverse)

message('Adjusting PGS for genetic PCs now...')

## read in genetic PCs
pc <- read_csv('ES-PGS_example/example_data/PCA_results.merged.1000_genomes_EUR.SLI_WGS.qc.demo.csv')
names(pc)[1] <- 'IID'

## 1000 genomes sample list
kg <- pc$IID[pc$population == '1000 Genomes']

## read in raw ES-PGS (PRSet output)
pgs <- read_table(file = 'ES-PGS_example/example_data/ES-PGS_raw.all_score',
                  show_col_types = FALSE)
prsice <- read_table(file = 'ES-PGS_example/example_data/ES-PGS_raw.prsice',
                     show_col_types = FALSE)

################################################
## adjust for PCs: loop over each PGS and adjust for PCs based on 1000 Genomes
################################################
names(pgs)[-c(1:2)] <- str_c('cp.', names(pgs)[-c(1:2)])

## merge PGS with genetic PCs
tmp <- pgs %>%
    inner_join(pc)

##
kg <- tmp %>%
    filter(population == '1000 Genomes')
dg <- tmp %>%
    filter(population != '1000 Genomes')
pgs_list_tmp <- list()
##
for(j in 3:ncol(pgs)){
    nm <- colnames(kg)[j]
    kg$PGS <- unname(unlist(kg[,colnames(kg)[j]]))
    dg$PGS <- unname(unlist(dg[,colnames(kg)[j]]))
    mod <- lm(PGS ~ pc1 + pc2 + pc3 + pc4 + pc5, 
                data = kg)
    ##
    preds_kg <- predict(mod, newdata = kg)
    preds_all <- predict(mod, newdata = dg)
    ##
    resids_kg = kg$PGS - preds_kg 
    resids_dg <- dg$PGS - preds_all
    ##
    kg_mean = mean(resids_kg)
    kg_sd = sd(resids_kg)
    z_kg <- (resids_kg - kg_mean) / kg_sd
    z = (resids_dg - kg_mean) / kg_sd
    dg$pgs_pc_corrected = z
    kg$pgs_pc_corrected = z_kg

    tmp2 <- bind_rows(dg, kg)
    tmp2 <- tmp2 %>%
        mutate(pgs_name = nm) %>%
        select(IID, pgs_name, PGS, pgs_pc_corrected)
    pgs_list_tmp[[j]] <- tmp2
}

## reformat data
pgs_all <- bind_rows(pgs_list_tmp)

pgs_gw <- pgs_all %>% 
    filter(str_detect(pgs_name, pattern = 'Base_')) %>% 
    select(IID, pgs_genome_wide = pgs_pc_corrected)

dat <- pgs_gw %>%
  inner_join(pgs_all) %>%
  inner_join(select(pc, IID, population)) %>%
  rename(pgs_raw = PGS) %>%
  mutate(cohort = ifelse(population == '1000 Genomes', '1000 Genomes (Europeans)', 'EpiSLI')) %>%
  relocate(cohort, .after = IID) %>%
  select(-population)

## make df w/ cols for complement + anno scores
datc <- dat %>% 
  filter(str_detect(pgs_name, 'complement_')) %>% 
  mutate(pgs_name = str_remove_all(pgs_name, pattern = 'complement_')) 
names(datc)[5:6] <- str_c(names(datc)[5:6], '_complement')

data <- dat %>% 
  filter(str_detect(pgs_name, 'complement_|Base_', negate = TRUE)) %>% 
  inner_join(datc) %>% 
  select(-matches('pgs_raw'))

################################################
## save data
################################################
data %>%
  write_csv('ES-PGS_example/example_data/gathered_es-pgs_corrected.cogPerf.csv')
message("ES-PGS adjusted for population stratification")
message('Saved analysis ready PGS to: ES-PGS_example/example_data/gathered_es-pgs_corrected.cogPerf.csv')