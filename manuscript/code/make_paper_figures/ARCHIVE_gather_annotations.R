library(tidyverse)

################3
to_copy <- list.files('/wdata/lcasten/sli_wgs/prs/gene_sets', pattern = 'bed$', full.names = TRUE)
to_copy <- to_copy[str_detect(to_copy, pattern = 'human_specific|human_brain_expression_divergence|brainCRE|phastCons|Iasi|_under_|primate_conserved_exons|primate_specific_information|introgressed_SNPs', negate = TRUE)]

for(f in to_copy) {
    fbase <- basename(f)
    outf <- str_c('manuscript/supplemental_materials/ES-PGS_annotations/', fbase)
    file.copy(from = f, to = outf)
}
