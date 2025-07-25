library(tidyverse)

message('Making figures now...')

####################################
## read in ES-PGS model results
####################################
res <- read_csv('ES-PGS_example/example_data/ES-PGS_lm_comparison_results.csv')

######################################
## make performance improvement figure
######################################
## initialize objects to store figures
fig_list_perf <- list()
fig_list_beta <- list()

## make figures and save to lists
for(ph in unique(res$phenotype)) {
    cat('\n\n')
    message('Making figures for ', ph)
    ## subset data to results from single phenotype
    tmp <- res %>% 
        filter(phenotype == ph) %>% 
        mutate(relative_rsq_gain_per_1000indSNP = ifelse(relative_rsq_gain_per_1000indSNP >= 1, 1, relative_rsq_gain_per_1000indSNP)) ## cap R-squared gain if > 100%
    ## make bar plot with relative performance increase compared to background PGS
    fig_perf <- tmp %>% 
        ggplot(aes(x = ES_PGS_model, y = relative_rsq_gain_per_1000indSNP, fill = ES_PGS_model)) +
        geom_bar(stat = 'identity', width = 0.5) +
        scale_y_continuous(labels = scales::percent, limits = c(0, 1)) +
        xlab('ES-PGS annotation') +
        ylab('Relative R-squared gain\nper 1000 SNPs in annotation') +
        theme_classic() +
        theme(legend.position = 'none',
              axis.text = element_text(size = 12),
              axis.title = element_text(size = 16, face = 'bold'),
              plot.title = element_text(size = 16, face = 'bold', hjust = 0.5)) +
        ggtitle(str_c('Phenotype: ', ph))
    fig_list_perf[[ph]] <- fig_perf
    ## make forest plot with ES-PGS betas
    fig_beta <- tmp %>% 
        ggplot(aes(x = ES_PGS_beta, y = ES_PGS_model)) +
        geom_point(size = 3) +
        geom_linerange(aes(xmin = ES_PGS_beta - 1.96 * ES_PGS_std_err, xmax = ES_PGS_beta + 1.96 * ES_PGS_std_err), size = 1.2) +
        geom_vline(xintercept = 0, linetype = 'dashed', size = 1.075, color = 'red') +
        xlab('ES-PGS annotation beta') +
        ylab('ES-PGS annotation') +
        theme_classic() +
        theme(legend.position = 'none',
              axis.text = element_text(size = 12),
              axis.title = element_text(size = 16, face = 'bold'),
              plot.title = element_text(size = 16, face = 'bold', hjust = 0.5)) +
        ggtitle(str_c('Phenotype: ', ph))
    fig_list_beta[[ph]] <- fig_beta
}


## save performance figures to file
pdf("ES-PGS_example/example_figures/ES-PGS_model_performance.pdf", onefile = TRUE)
for (i in seq(length(fig_list_perf))) {
    print(fig_list_perf[[i]])
}
dev.off()

## save forest plots to file
pdf("ES-PGS_example/example_figures/ES-PGS_model_betas.pdf", onefile = TRUE)
for (i in seq(length(fig_list_beta))) {
    print(fig_list_beta[[i]])
}
dev.off()

message('Saved figures to: "ES-PGS_example/example_figures" directory')