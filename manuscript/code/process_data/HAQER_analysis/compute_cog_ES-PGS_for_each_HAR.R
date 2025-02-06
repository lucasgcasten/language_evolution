library(tidyverse)

## get SNPs
snp_map <- read_table('/wdata/lcasten/sli_wgs/prs/pathway_prs/cogPerf.human_evolution_complement.snp')
har_snp <- snp_map %>% 
    filter(HAR == 1) %>% 
    select(1:4)
table(har_snp$CHR)

## subset PLINK file to those SNPs
har_snp %>% 
    arrange(CHR, BP) %>%
    select(SNP) %>% 
    write_tsv('/wdata/lcasten/sli_wgs/prs/pathway_prs/HAR_EpiSLI_1000Genomes_SNPs.txt', col_names = FALSE)

cmd <- c('/wdata/lcasten/tools/plink --bfile /wdata/lcasten/sli_wgs/prs/1000_genomes_EUR_merge/merged.1000_genomes_EUR.SLI_WGS.qc --extract /wdata/lcasten/sli_wgs/prs/pathway_prs/HAR_EpiSLI_1000Genomes_SNPs.txt --keep-allele-order --make-bed --out /wdata/lcasten/sli_wgs/prs/pathway_prs/HAR_EpiSLI_1000Genomes_independent_SNPs')
cmd
system(cmd)

## read in genotypes and reformat
bim <- read_table('/wdata/lcasten/sli_wgs/prs/pathway_prs/HAR_EpiSLI_1000Genomes_independent_SNPs.bim', col_names = FALSE)
names(bim) <- c('chromosome', 'SNP', 'CM', 'BP', 'ALT', 'REF')
bim <- bim %>% 
    mutate(end = BP) %>%
    select(chromosome, start = BP, end, SNP, REF, ALT)

geno <- BEDMatrix::BEDMatrix("/wdata/lcasten/sli_wgs/prs/pathway_prs/HAR_EpiSLI_1000Genomes_independent_SNPs.bed")
geno <- geno %>% 
    as.matrix() %>% 
    as.data.frame() %>% 
    rownames_to_column(var = 'IID') %>% 
    mutate(IID = str_split(IID, pattern = '_', simplify = TRUE)[,2]) %>% 
    as_tibble()
colnames(geno) <- str_split(colnames(geno), pattern = '_', simplify = TRUE)[,1]
geno
geno_vcf <- data.table::transpose(geno[,-1])
colnames(geno_vcf) <- geno$IID 
geno_vcf <- geno_vcf %>% 
    as.data.frame() %>% 
    as_tibble() %>% 
    mutate(SNP = colnames(geno)[-1]) %>% 
    relocate(SNP)
geno_vcf

## read in HARs
haq <- read_tsv('/wdata/lcasten/tools/ref_data/human_accelerated_regions/pollard_lab/nchaes_merged_hg19.bed', col_names = FALSE)
haq <- haq %>% 
    select(chromosome = X1, start = X2, end = X3, har = X4) %>% 
    mutate(chromosome = as.numeric(str_remove_all(chromosome, 'chr'))) %>% 
    drop_na(chromosome) %>% 
    arrange(chromosome, start, end) %>% 
    mutate(har = str_replace_all(har, pattern = '[;]', replacement = '__'))

## map HARs to EpiSLI variants
var_map <- bim %>% 
    fuzzyjoin::genome_inner_join(haq) %>% 
    select(chromosome = chromosome.x, BP = start.x, SNP, REF, ALT, har, har_start = start.y, har_end = end.y)
table(var_map$har)
tmp <- var_map %>% 
    group_by(har, har_start, har_end) %>% 
    count() %>% 
    mutate(har_size = har_end - har_start) %>% 
    ungroup()
cor.test(tmp$n, tmp$har_size)
tmp %>% 
    arrange(desc(har_size))

## read in sumstats
cog <- read_tsv('/wdata/lcasten/sli_wgs/polygenic_neanderthal_introgression/sumstats/cogPerf.tsv') %>% 
    rename(chromosome = CHR)
annot <- var_map %>% 
    inner_join(cog) %>% 
    mutate(BETA_matching_EpiSLI = case_when(ALT == A2 ~ -1 * BETA,
                                            TRUE ~ BETA))
annot %>% filter(ALT != A1) %>% nrow(.) ## number of flipped alleles

annot %>% 
    write_csv('/wdata/lcasten/sli_wgs/prs/pathway_prs/HAR_PGS_manual_calculation_annotations.csv')


## merge genotypes with annotations/weights
wd <- annot %>% 
    inner_join(geno_vcf)

## compute "PGS" for each HAR
wd_score <- wd %>% 
    pivot_longer(cols = matches('sample|HG|NA')) %>% 
    group_by(name) %>% 
    mutate(maf = mean(value, na.rm = TRUE)) %>%
    ungroup() %>%
    mutate(value = ifelse(is.na(value), 0, value)) %>% ## "impute" missing genos with MAF (same as PRSice)
    mutate(score = value * BETA_matching_EpiSLI) %>% 
    group_by(name, har) %>% 
    summarise(pgs_raw = sum(score)) %>% 
    select(IID = name, har, pgs_raw) %>% 
    ungroup()
wd_score
wd_score_sum = wd_score %>% 
    group_by(IID) %>% 
    summarise(pgs_raw = sum(pgs_raw)) %>% 
    ungroup()
hist(wd_score_sum$pgs_raw)
wd_score_wide <- wd_score %>% 
    group_by(har) %>% 
    mutate(pgs_raw = scale(pgs_raw)[,1]) %>%
    pivot_wider(id_cols= IID, names_from = har, values_from = pgs_raw)

## check the sum score is equivalent to actual PRSet scores and get background PGS
prset <- read_table('/wdata/lcasten/sli_wgs/prs/pathway_prs/cogPerf.human_evolution_complement.all_score') %>% 
    select(IID, background_pgs = Base_1, HAR_1)
tmp <- prset %>% 
    inner_join(wd_score_sum)
cor.test(tmp$pgs_raw, tmp$HAR_1) ## prset does some weird scaling?? but the correlation is essentially perfect (any difference likely due to rounding errors)

## save to file
pc <- read_csv('/wdata/lcasten/sli_wgs/PCA/PCA_results.merged.1000_genomes_EUR.SLI_WGS.qc.demo.csv') %>% 
    select(IID = sample.ID, pc1, pc2, pc3, pc4, pc5)

pc %>% 
    inner_join(select(prset, IID, background_pgs)) %>%
    inner_join(wd_score_sum) %>% 
    mutate(background_pgs = scale(background_pgs)[,1],
           pgs_raw = scale(pgs_raw)[,1]) %>%
    rename(HAR_pgs_sum = pgs_raw) %>% 
    inner_join(wd_score_wide) %>% 
    select(-ANC985) %>%
    write_csv('/wdata/lcasten/sli_wgs/prs/pathway_prs/HAR_PGS_manual_calculation_raw.csv')
