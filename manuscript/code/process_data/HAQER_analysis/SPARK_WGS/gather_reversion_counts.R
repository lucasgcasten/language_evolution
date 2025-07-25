
rand_files <- list.files('/Dedicated/jmichaelson-wdata/lcasten/sli_wgs/HAQER_variant_associations/data/SPARK_WGS', pattern = 'RAND_10Kb_flank.hg38.reversion_count', full.names = TRUE)
rand_files 
rand_list <- list()
for(f in rand_files) {
    rand_list[[basename(f)]] <- read_csv(f, show_col_types = FALSE)
}

rand_burden <- bind_rows(rand_list) %>% 
    group_by(sample) %>% 
    summarise(rand_rare_variant_count = sum(haqer_rare_variant_count),
              rand_rare_variant_reversion_count = sum(haqer_rare_variant_reversion_count)) %>% 
    rename(spid = sample)


har_files <- list.files('/Dedicated/jmichaelson-wdata/lcasten/sli_wgs/HAQER_variant_associations/data/SPARK_WGS', pattern = 'HAR_10Kb_flank.hg38.reversion_count', full.names = TRUE)
har_files 
har_list <- list()
for(f in har_files) {
    har_list[[basename(f)]] <- read_csv(f, show_col_types = FALSE)
}

har_burden <- bind_rows(har_list) %>% 
    group_by(sample) %>% 
    summarise(har_rare_variant_count = sum(haqer_rare_variant_count),
              har_rare_variant_reversion_count = sum(haqer_rare_variant_reversion_count)) %>% 
    rename(spid = sample)

##
haqer_burden <- read_csv('/wdata/lcasten/sli_wgs/HAQER_variant_associations/data/SPARK_WGS/HAQER_10Kb_flank_rare_variant_burden_counts.csv') 

rare_burden <- haqer_burden %>% 
    inner_join(har_burden) %>% 
    inner_join(rand_burden)

rare_burden %>% 
    write_csv('/wdata/lcasten/sli_wgs/HAQER_variant_associations/data/SPARK_WGS/HAQER_HAR_RAND_10Kb_flank_rare_variant_burden_counts.csv') 

cor(rare_burden[,-1])
