library(tidyverse)

##
bim1 <- read_table('/Dedicated/jmichaelson-sdata/SLI_WGS/callsets/ensemble/sli.seq.merge.common_variants.no_dups.bim', col_names = FALSE)
bim2 <- read_table('/Dedicated/jmichaelson-sdata/1KG/phase3/EUR/common_variants.no_dups.bim', col_names = FALSE)

##
length(intersect(bim1$X2, bim2$X2))

##
kp <- bim1 %>% 
    inner_join(bim2) %>% 
    filter(X2 != '.') %>%
    distinct(X2, .keep_all = TRUE)

kp %>% 
    write_delim('/wdata/lcasten/sli_wgs/prs/1000_genomes_EUR_merge/overlapping_variants.qc.bim', delim = ' ', col_names = FALSE)
kp %>% 
    select(X2) %>%
    write_delim('/wdata/lcasten/sli_wgs/prs/1000_genomes_EUR_merge/overlapping_variants.qc.snplist', delim = ' ', col_names = FALSE)

