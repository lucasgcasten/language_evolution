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
# data$sex_unknown = ifelse(data$molecular_sex == 'U', 1, 0)
# data$sex_female <- ifelse(data$molecular_sex != 'U' & data$molecular_sex != 'M', 1, 0)
lmm_mod <- lmm.aireml(Y = scale(data$age)[,1], 
                      X = as.matrix(data[,c("int", "cp_pgs.HAQER", "cp_pgs.background", "cp_pgs.matched")]), ## , "sex_unknown", "sex_female", "lat", "long"
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
    x = c("intercept", "cp_pgs.HAQER", "cp_pgs.background", "cp_pgs.matched"), ## , "sex_unknown", "sex_female", "lat", "long"
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
    drop_na(pgs) %>%                           
    mutate(sel_haq_lab = sel_haq_lab,
           sel_bg_lab = sel_bg_lab) %>%
    ggplot(aes(x = -1 * sample_age, y = value, color = pgs)) +
    geom_smooth(method = 'lm', size = 2, alpha = .2) +
    geom_text(aes(x = -8000, y = -2.1, label = sel_bg_lab), check_overlap = TRUE, size = 7, color = 'grey50') +
    geom_text(aes(x = -8000, y = .5, label = sel_haq_lab), check_overlap = TRUE, size = 7, color = '#762776') +
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


p_sel2 <- data %>% 
    pivot_longer(cols = matches('cp_pgs')) %>% 
    mutate(pgs = case_when(name == 'cp_pgs.HAQER' ~ 'HAQERs',
                           name == 'cp_pgs.background' ~ 'Background')) %>%
    drop_na(pgs) %>%                           
    mutate(sel_haq_lab = sel_haq_lab,
           sel_bg_lab = sel_bg_lab) %>%
    ggplot(aes(x = -1 * sample_age, y = value, color = pgs)) +
    geom_smooth(size = 2, alpha = .2) +
    geom_text(aes(x = -8000, y = -2.1, label = sel_bg_lab), check_overlap = TRUE, size = 7, color = 'grey50') +
    geom_text(aes(x = -8000, y = .5, label = sel_haq_lab), check_overlap = TRUE, size = 7, color = '#762776') +
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
p_sel2 %>% 
    ggsave(filename = 'manuscript/figures/AADR_HAQER_polygenic_selection2.png',
           device = 'png', dpi = 300, bg = 'white',
           units = 'in', width = 5, height = 5)

## ---------------------------------------------------
## neanderthal ES-PGS
## ---------------------------------------------------
cl <- c("#762776", "#e04468", "#dcc699") ## color palette (HAQER, HAR, RAND)

## read in PGS data
df <- read_csv("manuscript/supplemental_materials/AADR_1000_genomes_ES-PGS_data.csv") %>% 
    ## add each group's N for figure
    mutate(type = factor(type, levels = rev(c('Modern Europeans', 'Ancient Europeans', 'Neanderthals and Denisovans')))) %>% 
    arrange(type) %>% 
    group_by(type) %>% 
    mutate(n = n(),
           type = str_c(type, '\n(N = ', prettyNum(n, big.mark = ','), ')')) %>% 
    ungroup() %>% 
    mutate(type = ifelse(type == 'Neanderthals and Denisovans\n(N = 10)', 'Neanderthals and\nDenisovans (N = 10)', type)) %>%
    mutate(type = factor(type, levels = unique(type)))

## compare raw PGS across groups (ex. neanderthals vs modern humans) to look for trends
res_haq <- broom::tidy(pairwise.t.test(x = df$cp_pgs.HAQER, g = df$type, p.adjustment.method = 'none')) %>% 
    mutate(x = 'cp_pgs.HAQER') %>% 
    relocate(x)
res_bg <- broom::tidy(pairwise.t.test(x = df$cp_pgs.background, g = df$type)) %>% 
    mutate(x = 'cp_pgs.background') %>% 
    relocate(x)
res_mt <- broom::tidy(pairwise.t.test(x = df$cp_pgs.matched, g = df$type)) %>% 
    mutate(x = 'cp_pgs.matched') %>% 
    relocate(x)

## gather PGS group comparison stats and save
pgs_sumstats <- df %>% 
    pivot_longer(cols = matches('cp_pgs')) %>% 
    group_by(name, type, n) %>% 
    summarise(mean_pgs = mean(value),
              median_pgs = median(value),
              sd_pgs = sd(value)) %>% 
    ungroup() %>% 
    rename(x = name)

tmp1 <- pgs_sumstats  %>% 
    rename(group1 = type) %>% 
    mutate(group1 = str_split(group1, pattern = '[\n]', simplify = TRUE)[,1])
names(tmp1)[-c(1:2)] <- str_c("group1_", names(tmp1)[-c(1:2)])

tmp2 <- pgs_sumstats  %>% 
    rename(group2 = type) %>% 
    mutate(group2 = str_split(group2, pattern = '[\n]', simplify = TRUE)[,1]) %>% 
    mutate(group2 = ifelse(group2 == 'Neanderthals and', 'Neanderthals and Denisovans', group2))
names(tmp2)[-c(1:2)] <- str_c("group2_", names(tmp2)[-c(1:2)])

res_comp <- bind_rows(res_haq, res_bg, res_mt) %>% 
    mutate(group1 = str_split(group1, pattern = '[\n]', simplify = TRUE)[,1],
           group2 = str_split(group2, pattern = '[\n]', simplify = TRUE)[,1]) %>% 
    mutate(group2 = ifelse(group2 == 'Neanderthals and', 'Neanderthals and Denisovans', group2))

res_comp %>% 
    inner_join(tmp1) %>% 
    inner_join(tmp2) %>%
    write_csv("manuscript/supplemental_materials/stats/AADR_HAQER_ES-PGS_polygenic_group_comparison_results.csv")

## make 3 boxplots comparing neanderthals (AADR) vs ancient europeans (AADR) vs modern europeans (1000 Genomes) across 3 PGS' (HAQER, background, and matched)
### HAQER pgs
p_nean_haq <- df %>% 
    ggplot(aes(x = type, y = cp_pgs.HAQER)) +
    geom_violin(size = 1.5) +
    geom_boxplot(aes(fill = type), size = 1.5, width = .3, alpha = .75) +
    geom_hline(yintercept = mean(df$cp_pgs.HAQER), size = 1.075, color = 'red', linetype = 'dashed') +
    xlab(NULL) +
    ylab('HAQER CP-PGS') +
    theme_classic() +
    theme(axis.text = element_text(size = 12),
          axis.title = element_text(size = 14),
          strip.text = element_text(size = 14),
          legend.text = element_text(size = 12),
          legend.position = 'none') +
    scale_fill_manual(name = NULL, values = viridis::viridis(3, option = "D")) +
    ## add significance comparing neanderthals vs ancient europeans
    geom_linerange(aes(xmin = 1, xmax = 2, y = 3.8), size = 1.075) +
    geom_text(aes(x = 1.5, y = 3.85, label = "*"), size = 7, check_overlap = TRUE) + ## check t-test results to get significance level
    ## add significance comparing neanderthals vs modern europeans
    geom_linerange(aes(xmin = 1, xmax = 3, y = 4.05), size = 1.075) +
    geom_text(aes(x = 2, y = 4.1, label = "*"), size = 7, check_overlap = TRUE) ## check t-test results to get significance level

### background PGS
p_bg_haq <- df %>% 
    ggplot(aes(x = type, y = cp_pgs.background)) +
    geom_violin(size = 1.5, width = 1.25) +
    geom_boxplot(aes(fill = type), size = 1.1, width = .2, alpha = .75) +
    geom_hline(yintercept = mean(df$cp_pgs.background), size = 1.075, color = 'red', linetype = 'dashed') +
    xlab(NULL) +
    ylab('Background CP-PGS') +
    theme_classic() +
    theme(axis.text = element_text(size = 12),
          axis.title = element_text(size = 14),
          strip.text = element_text(size = 14),
          legend.text = element_text(size = 12),
          legend.position = 'none') +
    scale_fill_manual(name = NULL, values = viridis::viridis(3, option = "D")) +
    ## add significance comparing neanderthals vs ancient europeans
    geom_linerange(aes(xmin = 1, xmax = 2, y = 3.95), size = 1.075) +
    geom_text(aes(x = 1.5, y = 4, label = "***"), size = 7, check_overlap = TRUE) + ## check t-test results to get significance level
    ## add significance comparing neanderthals vs modern europeans
    geom_linerange(aes(xmin = 1, xmax = 3, y = 4.25), size = 1.075) +
    geom_text(aes(x = 2, y = 4.3, label = "***"), size = 7, check_overlap = TRUE) + ## check t-test results to get significance level
    ## add significance comparing ancient europeans vs modern europeans
    geom_linerange(aes(xmin = 2, xmax = 3, y = 4.55), size = 1.075) +
    geom_text(aes(x = 2.5, y = 4.6, label = "***"), size = 7, check_overlap = TRUE) ## check t-test results to get significance level

### matched control PGS
p_mt_haq <- df %>% 
    ggplot(aes(x = type, y = cp_pgs.matched)) +
    geom_violin(size = 1.5) +
    geom_boxplot(aes(fill = type), size = 1.5, width = .2, alpha = .75) +
    geom_hline(yintercept = mean(df$cp_pgs.matched), size = 1.075, color = 'red', linetype = 'dashed') +
    xlab(NULL) +
    ylab('Matched CP-PGS') +
    theme_classic() +
    theme(axis.text = element_text(size = 12),
          axis.title = element_text(size = 14),
          strip.text = element_text(size = 14),
          legend.text = element_text(size = 12),
          legend.position = 'none') +
    scale_fill_manual(name = NULL, values = viridis::viridis(3, option = "D"))

## save plots
p_bg_haq %>% 
    ggsave(filename = 'manuscript/figures/AADR_neanderthal_background_ES-PGS.png',
           device = 'png', dpi = 300, bg = 'white',
           units = 'in', width = 6, height = 6)
p_nean_haq %>% 
    ggsave(filename = 'manuscript/figures/AADR_neanderthal_HAQER_ES-PGS.png',
           device = 'png', dpi = 300, bg = 'white',
           units = 'in', width = 6, height = 6)
p_mt_haq %>% 
    ggsave(filename = 'manuscript/figures/AADR_neanderthal_matched_ES-PGS.png',
           device = 'png', dpi = 300, bg = 'white',
           units = 'in', width = 6, height = 6)

## -------------------------------------------------
## make set of ES-PGS boxplots for larger figure
## -------------------------------------------------
p_nean_haq_fig5 <- df %>% 
    ggplot(aes(x = type, y = cp_pgs.HAQER)) +
    geom_violin(size = 1.5) +
    geom_boxplot(aes(fill = type), size = 1.5, width = .3, alpha = .75) +
    geom_hline(yintercept = mean(df$cp_pgs.HAQER), size = 1.075, color = 'red', linetype = 'dashed') +
    xlab(NULL) +
    ylab('HAQER CP-PGS') +
    scale_fill_manual(name = NULL, values = c('slategray2', viridis::viridis(3, option = "D")[2:3])) +
    theme_classic() +
    theme(axis.text = element_text(size = 12),
          axis.text.x = element_blank(),
          axis.title = element_text(size = 14),
          strip.text = element_text(size = 14),
          legend.text = element_text(size = 20),
          legend.position = 'bottom') +
    ## add significance comparing neanderthals vs ancient europeans
    geom_linerange(aes(xmin = 1, xmax = 2, y = 4.15), size = 1.075) +
    geom_text(aes(x = 1.5, y = 4.2, label = "*"), size = 7, check_overlap = TRUE) + ## check t-test results to get significance level
    ## add significance comparing neanderthals vs modern europeans
    geom_linerange(aes(xmin = 1, xmax = 3, y = 3.8), size = 1.075) +
    geom_text(aes(x = 2, y = 3.85, label = "*"), size = 7, check_overlap = TRUE) ## check t-test results to get significance level

### background PGS
p_bg_haq_fig5 <- df %>% 
    ggplot(aes(x = type, y = cp_pgs.background)) +
    geom_violin(size = 1.5, width = 1.25) +
    geom_boxplot(aes(fill = type), size = 1.5, width = .2, alpha = .75) +
    geom_hline(yintercept = mean(df$cp_pgs.background), size = 1.075, color = 'red', linetype = 'dashed') +
    xlab(NULL) +
    ylab('Background CP-PGS') +
    scale_fill_manual(name = NULL, values = c('slategray2', viridis::viridis(3, option = "D")[2:3])) +
    theme_classic() +
    theme(axis.text = element_text(size = 12),
          axis.text.x = element_blank(),
          axis.title = element_text(size = 14),
          strip.text = element_text(size = 14),
          legend.text = element_text(size = 20),
          legend.position = 'bottom') +
    ## add significance comparing neanderthals vs ancient europeans
    geom_linerange(aes(xmin = 1, xmax = 2, y = 3.95), size = 1.075) +
    geom_text(aes(x = 1.5, y = 4, label = "***"), size = 7, check_overlap = TRUE) + ## check t-test results to get significance level
    ## add significance comparing neanderthals vs modern europeans
    geom_linerange(aes(xmin = 1, xmax = 3, y = 4.3), size = 1.075) +
    geom_text(aes(x = 2, y = 4.35, label = "***"), size = 7, check_overlap = TRUE) + ## check t-test results to get significance level
    ## add significance comparing ancient europeans vs modern europeans
    geom_linerange(aes(xmin = 2, xmax = 3, y = 4.65), size = 1.075) +
    geom_text(aes(x = 2.5, y = 4.7, label = "***"), size = 7, check_overlap = TRUE) ## check t-test results to get significance level

### matched control PGS
p_mt_haq_fig5 <- df %>% 
    ggplot(aes(x = type, y = cp_pgs.matched)) +
    geom_violin(size = 1.5) +
    geom_boxplot(aes(fill = type), size = 1.5, width = .2, alpha = .75) +
    geom_hline(yintercept = mean(df$cp_pgs.matched), size = 1.075, color = 'red', linetype = 'dashed') +
    xlab(NULL) +
    ylab('Matched CP-PGS') +
    scale_fill_manual(name = NULL, values = c('slategray2', viridis::viridis(3, option = "D")[2:3])) +
    theme_classic() +
    theme(axis.text = element_text(size = 12),
          axis.text.x = element_blank(),
          axis.title = element_text(size = 14),
          strip.text = element_text(size = 14),
          legend.text = element_text(size = 20),
          legend.position = 'bottom')

library(patchwork)
fig5_neand <- (p_nean_haq_fig5 + p_bg_haq_fig5 + p_mt_haq_fig5) + plot_layout(guides = 'collect') & theme(legend.position = 'bottom')

###################################
## save plot objects to merge later
p_bg_haq %>% 
    write_rds('manuscript/figures/R_plot_objects/AADR_background-CP-PGS_nean_dist.rds')
p_nean_haq %>% 
    write_rds('manuscript/figures/R_plot_objects/AADR_HAQER-CP-PGS_nean_dist.rds')
p_mt_haq %>% 
    write_rds('manuscript/figures/R_plot_objects/AADR_matched-CP-PGS_nean_dist.rds')    
p_sel %>% 
    write_rds('manuscript/figures/R_plot_objects/AADR_HAQER-CP-PGS_selection.rds')
fig5_neand %>% 
    write_rds("manuscript/figures/R_plot_objects/AADR_HAQER-CP-PGS_nean_dist_all_merged.rds")
