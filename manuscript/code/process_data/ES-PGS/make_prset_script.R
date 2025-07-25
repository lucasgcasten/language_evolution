library(tidyverse)

## get sumstats
sdir <- '/Dedicated/jmichaelson-wdata/lcasten/sli_wgs/polygenic_neanderthal_introgression/sumstats'
ss <- list.files(sdir, full.names = TRUE)
ss

## plink files
bfile = '/Dedicated/jmichaelson-wdata/lcasten/sli_wgs/prs/1000_genomes_EUR_merge/merged.1000_genomes_EUR.SLI_WGS.qc'
## neurotransmitter gmt filepath
gmt <- '/Dedicated/jmichaelson-wdata/lcasten/sli_wgs/prs/gene_sets/custom_gene_sets.gmt'
## get samples for LD (1000 Genomes Eur)
fam <- read_table('/Dedicated/jmichaelson-wdata/lcasten/sli_wgs/prs/1000_genomes_EUR_merge/merged.1000_genomes_EUR.SLI_WGS.qc.fam', col_names = FALSE)
fam %>% 
    filter(str_detect(X2, pattern = 'sample', negate = TRUE)) %>% 
    select(1:2) %>% 
    write_delim('/Dedicated/jmichaelson-wdata/lcasten/sli_wgs/prs/1000_genomes_EUR_merge/LD_keep_samples.txt', delim = ' ', col_names = FALSE)

## make prsset command
## example:
# Rscript PRSice.R \
#     --prsice ./bin/PRSice \
#     --base TOY_BASE_GWAS.assoc \
#     --target TOY_TARGET_DATA \
#     --binary-target T \
#     --thread 4 \
#     --gtf gene.gtf \
#     --msigdb set.txt \
#     --multi-plot 10
prsice = "Rscript PRSice.R"
tc = 54

gtf = '/Dedicated/jmichaelson-wdata/lcasten/tools/ref_data/hg19/Homo_sapiens.GRCh37.87.gtf'
outd = "/Dedicated/jmichaelson-wdata/lcasten/sli_wgs/prs/pathway_prs/"

# list.files('/wdata/lcasten/sli_wgs/prs/gene_sets', pattern = 'bed$')
bedstring = "/Dedicated/jmichaelson-wdata/lcasten/sli_wgs/prs/gene_sets/brainCREs.bed:brainCREs,/Dedicated/jmichaelson-wdata/lcasten/sli_wgs/prs/gene_sets/neanderthal_selective_sweep.bed:NeanderthalSelectiveSweep,/Dedicated/jmichaelson-wdata/lcasten/sli_wgs/prs/gene_sets/HAQERs.bed:HAQER,/Dedicated/jmichaelson-wdata/lcasten/sli_wgs/prs/gene_sets/HARs.bed:HAR,/Dedicated/jmichaelson-wdata/lcasten/sli_wgs/prs/gene_sets/human_brain_expression_divergence.Div_in_ExN.bed:diverged_ExN,/Dedicated/jmichaelson-wdata/lcasten/sli_wgs/prs/gene_sets/human_brain_expression_divergence.Div_in_InN.bed:diverged_InN,/Dedicated/jmichaelson-wdata/lcasten/sli_wgs/prs/gene_sets/human_brain_expression_divergence.Div_in_NonN.bed:diverged_NonN,/Dedicated/jmichaelson-wdata/lcasten/sli_wgs/prs/gene_sets/introgressed_SNPs.bed:archaic_introgression,/Dedicated/jmichaelson-wdata/lcasten/sli_wgs/prs/gene_sets/neanderthal_introgressed_SNPs.bed:neanderthal_introgression"
dir.create(outd)

##
cmd_list = list()
i = 0
for(f in ss){
    cat(sprintf('\n\n\n\n\n'))
    message('===============================================')
    message('working on sumstats: ', basename(f))
    i = i + 1
    trait <- basename(f)
    trait = str_remove_all(trait, pattern = '.tsv')

    ## evo stuff
    cmd <- str_c(prsice, '--prsice PRSice_linux --base', f, 
              '--target', bfile, '--ld-keep /Dedicated/jmichaelson-wdata/lcasten/sli_wgs/prs/1000_genomes_EUR_merge/LD_keep_samples.txt --binary-target T',
              '--bar-levels 0.00000005,0.0005,0.05,0.2,1 --fastscore --no-regress',
              '--thread', tc, '--gtf', gtf, '--bed', bedstring, 
              '--out', paste0(outd, trait, '.human_evolution'), sep = ' ')
    ## neurotransmitters
    cmd2 <- str_c(prsice, '--prsice PRSice_linux --base', f, 
              '--target', bfile, '--ld-keep /Dedicated/jmichaelson-wdata/lcasten/sli_wgs/prs/1000_genomes_EUR_merge/LD_keep_samples.txt --binary-target T',
              '--bar-levels 0.00000005,0.0005,0.05,0.2,1 --fastscore --wind-5 35K --wind-3 10K --no-regress',
              '--thread', tc, '--gtf', gtf, '--msigdb', gmt, 
              '--out', paste0(outd, trait, '.custom_gene_sets'), sep = ' ')
    cmds <- c(cmd, cmd2) 
    cmd_list[[i]] <- cmds
}

##
bash = '#!/bin/bash'
dir <- 'cd /Dedicated/jmichaelson-wdata/lcasten/tools/PRSice'
script = c(bash, dir, unlist(cmd_list))
message('PRSset script that will be run:')
print(script)
script_out = '/Dedicated/jmichaelson-wdata/lcasten/sli_wgs/polygenic_neanderthal_introgression/code/prset_v3.sh'
write_lines(script, file = script_out)

message('==================================')
message('to run pathway PRS script call:')
message(str_c('bash ', script_out))
system(str_c('head -n 4 ', script_out))
# system(str_c('bash ', script_out))
