library(tidyverse)

##
f <- "/Dedicated/jmichaelson-wdata/lcasten/sli_wgs/prs/pathway_prs/SPARK_ABCD.cogPerf.HAQER.pathways.all_score"
pgs <- read_table(f)
dim(pgs)

trait = 'cogPerf_HAQER'
prsice <- read_table(file = str_replace_all(f, pattern = 'all_score', replacement = 'prsice'),
                show_col_types = FALSE)
to_drop <- prsice %>% 
    filter(Num_SNP < 10)

to_keep <- prsice %>% 
    filter(Num_SNP >= 10) %>% 
    filter(Set != 'Base')

pgs <- pgs %>% 
    select(IID, any_of(str_c(to_keep$Set, '_1')))
names(pgs) <- str_replace_all(names(pgs), pattern = '-', replacement = '_')

################
pc <- read_csv('/Dedicated/jmichaelson-wdata/lcasten/spark/prs/HapMap3_plus/PCA/raw_KING_pca_results.csv')
names(pc)[1] <- 'IID'
pca <- pc %>% 
    select(IID, str_c('pc', 1:5))

##
unrel <- read_table("/Dedicated/jmichaelson-wdata/lcasten/spark/prs/pathway/SPARK_ABCD_unrelated_europeans_for_LD.txt", col_names = FALSE)


#################
tmp <- pgs %>%
    inner_join(pca)

##
kg_samp <- unrel

kg <- tmp %>%
    filter(IID %in% kg_samp$X2)
dg <- tmp %>%
    filter(! IID %in% kg_samp$X2)
pgs_list_tmp <- list()
##
for(j in 2:ncol(pgs)){
    if(j %% 100 == 0) {
        message(j, '/', ncol(pgs))
    }
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
pgs_corrected <- bind_rows(pgs_list_tmp)

##
pgs_corrected %>% 
    write_csv('/Dedicated/jmichaelson-wdata/lcasten/sli_wgs/prs/pathway_prs/SPARK_ABCD.cogPerf.HAQER_pathway.csv')

