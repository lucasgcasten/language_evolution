## -------------------------------------
## polygenic selection
## -------------------------------------
######################
## load packages
######################
library(tidyverse)
library(data.table)
library(Matrix)
library(gaston)

######################
## load GRM
######################
k <- readRDS("manuscript/supplemental_materials/AADR_GCTA_GRM.rds")

#######################################
## Load phenotype and PGS
#######################################
data <- fread("manuscript/supplemental_materials/AADR_data.csv", header = TRUE)

## make sure everyone with pheno and PGS data is in the GRM
data <- data[data$IID %in% rownames(k), ]

## add intercept column for later (lmm.aireml won't warn you and you'll get crazy inflated selection estimates)
data$int <- 1 # add col of 1's for the intercept 

## log transform, flip, and scale sample age so we can interpret the selection coefficient
data$age <- -1 * log10(data$sample_age)

## double check IDs are in same order between data and GRM
identical(data$IID, colnames(k)) ## identical order and elements of GRM columns with our data
identical(data$IID, rownames(k)) ## identical order and elements of GRM rows with our data

#################################################################
## ES-PGS selection test: accounting for relatedness
#################################################################
## fit the GREML LMM: age ~ background_pgs + es_pgs (this will take a few minutes...)
lmm_mod <- lmm.aireml(Y = scale(data$age)[,1], 
                      X = as.matrix(data[,c("int", "cp_pgs.HAQER", "cp_pgs.background", "cp_pgs.matched")]), 
                      K = k)

## extract betas (intercept, cp_pgs.HAQER, cp_pgs.background)
betas <- lmm_mod$BLUP_beta

## extract SEs (intercept, cp_pgs.HAQER, cp_pgs.background)
var_beta <- lmm_mod$varbeta
se_betas <- sqrt(diag(var_beta))

## get t-statistics and p-vals
t_stats <- betas / se_betas
pvals <- 2 * pt(-abs(t_stats), df = nrow(data) - length(betas))

## gather stats
sel_results <- data.frame(
    y = rep('sample_age', times = length(betas)),
    x = c("intercept", "cp_pgs.HAQER", "cp_pgs.background", "cp_pgs.matched"),
    beta = betas,
    se = se_betas,
    t_stat = t_stats,
    pval = pvals
)

## print result table
print(sel_results)

## save results
write.csv(sel_results[sel_results$x != 'intercept',], 
          file = "manuscript/supplemental_materials/stats/AADR_HAQER_ES-PGS_polygenic_selection_results.csv", 
          quote = FALSE, row.names = FALSE)

###############################
## polygenic selection figure
###############################
library(tidyverse)
## get GRM residual age value and convert to a "raw age" so we can better visualize what LMM is testing
mn <- mean(data$age)
sdev <- sd(data$age)
data$age_resid_scaled <- (scale(data$age)[,1] - lmm_mod$BLUP_omega) * sdev + mn
data$sample_age_GRM_adj <- round(10^(-1 * data$age_resid_scaled))

## labs for plot
sel_results <- sel_results %>% 
    mutate(lab = str_c('selection coef. = ', round(beta, digits = 3), ', p-val = ', formatC(pval, digits = 2)))
sel_haq_lab <- sel_results$lab[sel_results$x == 'cp_pgs.HAQER']
sel_bg_lab <- sel_results$lab[sel_results$x == 'cp_pgs.background']

## plot
p_sel <- data %>% 
    pivot_longer(cols = matches('cp_pgs')) %>% 
    mutate(pgs = case_when(name == 'cp_pgs.HAQER' ~ 'HAQERs',
                           name == 'cp_pgs.background' ~ 'Background')) %>%
    mutate(sel_haq_lab = sel_haq_lab,
           sel_bg_lab = sel_bg_lab) %>%
    ggplot(aes(x = -1 * round(sample_age_GRM_adj), y = value, color = pgs)) +
    geom_smooth(method = 'lm', size = 1.5, alpha = .2) +
    geom_text(aes(x = -8000, y = -3, label = sel_bg_lab), check_overlap = TRUE, size = 6, color = 'grey50') +
    geom_text(aes(x = -8000, y = 1, label = sel_haq_lab), check_overlap = TRUE, size = 6, color = '#762776') +
    xlab('Years ago') +
    ylab('CP-PGS') +
    scale_color_manual(values = c('grey70', "#762776")) +
    theme_classic() +
    labs(color = NULL) +
    theme(axis.text = element_text(size = 16),
          axis.title = element_text(size = 18),
          strip.text = element_text(size = 18),
          legend.text = element_text(size = 16),
          plot.title = element_text(size = 18, hjust = .5),
          legend.position = 'bottom') +
    scale_x_continuous(breaks = c(seq(-20000, -5000, 5000), -100), labels = c('20,000', '15,000', '10,000', '5,000', '100')) +
    labs(color = 'CP-PGS:')
p_sel %>% 
    ggsave(filename = 'manuscript/figures/AADR_HAQER_polygenic_selection.png',
           device = 'png', dpi = 300, bg = 'white',
           units = 'in', width = 5, height = 5)


## ---------------------------------------------------
## neanderthal ES-PGS
## ---------------------------------------------------
cl <- c("#762776", "#e04468", "#dcc699") ## color palette (HAQER, HAR, RAND)
## read in PGS data
df <- read_csv('manuscript/supplemental_materials/neanderthal_1000Genomes_raw_ES-PGS_data.csv') %>% 
    filter(type != "EpiSLI") %>% ## drop EpiSLI samples 
    mutate(cp_pgs.background = scale(cp_pgs.background)[,1],
           cp_pgs.HAQER = scale(cp_pgs.HAQER)[,1])

## make density plots with neanderthals vs 1000 Genomes Europeans
nean_id <- c('AltaiNea', 'Chagyrskaya-Phalanx', 'DenisovaPinky', 'Vindija33.19')
tmp <- df %>% 
    filter(IID %in% nean_id) %>% 
    mutate(type = 'Neanderthals and Denisovan')
kg_df <- df %>% 
    filter(type == '1000 Genomes')

p_nean_haq <- kg_df %>% 
    mutate(type = 'Modern Europeans') %>%
    ggplot(aes(x = cp_pgs.HAQER)) +
    geom_density(aes(fill = type), alpha = 0.7) +
    xlab('HAQER CP-PGS') +
    ylab('Density') +
    theme_classic() +
    theme(axis.text = element_text(size = 12),
          axis.title = element_text(size = 14),
          strip.text = element_text(size = 14),
          legend.text = element_text(size = 12),
          legend.position = 'bottom') +
    geom_vline(data = tmp, aes(xintercept = cp_pgs.HAQER, color = type), size = 1, linetype = 'dashed') +
    scale_color_manual(name = NULL, values = c("Neanderthals and Denisovan" = "black")) +
    scale_fill_manual(name = NULL, values = c('Modern Europeans' = cl[1])) +
    guides(fill = guide_legend(order = 1, override.aes = list(alpha = 0.7)),
           color = guide_legend(order = 2, override.aes = list(size = 1.1, linetype = 'solid')))

p_nean_bg <- kg_df %>% 
    mutate(type = 'Modern Europeans') %>%
    ggplot(aes(x = cp_pgs.background)) +
    geom_density(aes(fill = type), alpha = 0.7) +
    xlab('Background CP-PGS') +
    ylab('Density') +
    theme_classic() +
    theme(axis.text = element_text(size = 12),
          axis.title = element_text(size = 14),
          strip.text = element_text(size = 14),
          legend.text = element_text(size = 12),
          legend.position = 'bottom') +
    geom_vline(data = tmp, aes(xintercept = cp_pgs.background, color = type), size = 1, linetype = 'dashed') +
    scale_color_manual(name = NULL, values = c("Neanderthals and Denisovan" = "black")) +
    scale_fill_manual(name = NULL, values = c('Modern Europeans' = 'grey75')) +
    guides(fill = guide_legend(order = 1, override.aes = list(alpha = 0.7)),
           color = guide_legend(order = 2, override.aes = list(size = 1.1, linetype = 'solid')))

## save density plots
p_nean_bg %>% 
    ggsave(filename = 'manuscript/figures/AADR_neanderthal_background_ES-PGS.png',
           device = 'png', dpi = 300, bg = 'white',
           units = 'in', width = 5, height = 5)
p_nean_haq %>% 
    ggsave(filename = 'manuscript/figures/AADR_neanderthal_HAQER_ES-PGS.png',
           device = 'png', dpi = 300, bg = 'white',
           units = 'in', width = 5, height = 5)


###################################
## save plot objects to merge later
p_nean_bg %>% 
    write_rds('manuscript/figures/R_plot_objects/AADR_background-CP-PGS_nean_dist.rds')
p_nean_haq %>% 
    write_rds('manuscript/figures/R_plot_objects/AADR_HAQER-CP-PGS_nean_dist.rds')
p_sel %>% 
    write_rds('manuscript/figures/R_plot_objects/AADR_HAQER-CP-PGS_selection.rds')
