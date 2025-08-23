## load packages
library(tidyverse)
library(ape)
library(phylolm)

## read in phylogeny downloaded from UCSC (https://hgdownload.soe.ucsc.edu/goldenPath/hg38/cactus447way/hg38.447way.nh.txt)
tr <- read.tree('manuscript/supplemental_materials/hg38.447way.nh.txt')

## read in data
cs_dat <- read_csv('manuscript/supplemental_materials/cross_species_vocal_learning_HAQER_similarity.csv')

## make birth:adult weight ratio (proxy for obstetric dilemma score)
cs_dat$obstetric_risk = scale(log10(cs_dat$neonate_body_mass_g_PanTHERIA) / log10(cs_dat$adult_body_mass_g_PanTHERIA))[,1]

## z-scale the sequence similarity scores
cs_dat$HAQER_sequence_similarity_scaled <- scale(cs_dat$HAQER_sequence_similarity)[,1]
cs_dat$HAR_sequence_similarity_scaled <- scale(cs_dat$HAR_sequence_similarity)[,1]

## reformat data to work with phylolm
cs_dat <- as.data.frame(cs_dat)
rownames(cs_dat) <- cs_dat$species_name

## ------------------------------------------------------------------------------------------
## phylogenetic logistic regression predicting vocal learning using HAQER sequence identity
## ------------------------------------------------------------------------------------------
## fit phylogenetic regression model
fit = phyloglm(vocal_learner_binary ~ HAQER_sequence_similarity_scaled, phy = tr, data = cs_dat, method = 'logistic_IG10')
summary(fit)

## pull stats for figure
b <- fit$coefficients['HAQER_sequence_similarity_scaled']
b_plot <- unname(round(b, digits = 2))
p <- summary(fit)$coefficients['HAQER_sequence_similarity_scaled','p.value']
# p_plot <- formatC(p, digits = 2)
p_plot <- format(round(p, 7), nsmall = 3)
# p_vl_lab <- str_c('beta = ', round(b, digits = 2), '\np = ', formatC(p, digits = 2))
p_vl_lab <- bquote(beta == .(b_plot) ~ " " ~ italic(p) == .(p_plot))

b_part <- bquote(beta == .(b_plot))
p_parts <- strsplit(format(p, scientific = TRUE, digits = 2), "e")[[1]]
base <- as.numeric(p_parts[1])
exponent <- as.numeric(p_parts[2])
p_text <- bquote(italic(p) == .(round(base, 2)) %*% 10^.(exponent))
# Combine all expressions
stats_expr_vl <- as.expression(bquote(.(b_part) ~ "" ~ .(p_text)))

## gather stats
n <- nrow(cs_dat)
n_voc <- nrow(cs_dat[cs_dat$vocal_learner_binary == 1,])
n_nonvoc <- n - n_voc

res_vl <- as.data.frame(summary(fit)$coefficients) %>% 
    rownames_to_column(var = 'x') %>% 
    mutate(y = 'vocal_learner_binary') %>% 
    as_tibble() %>% 
    relocate(y) %>% 
    rename(beta = Estimate, std.error = StdErr, statistic = z.value) %>% 
    filter(x == 'HAQER_sequence_similarity_scaled') %>% 
    mutate(n_species = n,
           n_vocal_learning_species = n_voc,
           n_non_vocal_learning_species = n_nonvoc)

## make figure comparing vocal learners and non-vocal learners for HAQER-like seq sim
set.seed(65666)
cl <- sample(RColorBrewer::brewer.pal(name="Dark2",n=3))
p_vl <- cs_dat %>%
    group_by(vocal_learner_binary) %>% 
    mutate(n = n()) %>%
    ungroup() %>%
    mutate(vlab = as.logical(vocal_learner_binary),
            vlab = str_c(vlab, '\nN = ', n, ' species')
            ) %>%
    mutate(common_name = case_when(str_detect(common_name, pattern = 'Egyptian fruit bat| red bat|ottlenose dolphin|Brazilian gu|tree shrew|Sumatran rhino|Cuvier|Orca|Killer|dog|Domestic cat') ~ common_name,
                                   TRUE ~ NA_character_)) %>%
    ggplot(aes(x = vlab, y = HAQER_sequence_similarity_scaled)) +
    geom_violin(linewidth = 0, aes(fill = vlab), alpha = 0.4) +
    geom_point(aes(color = vlab), size = .7) +
    ggrepel::geom_text_repel(aes(label = common_name, color = vlab), size = 3.25, max.overlaps = 15) +
    xlab('Vocal learner (Wirthlin, et al.)') +
    ylab('HAQER-like sequence score') +
    scale_fill_manual(values = cl[1:2]) +
    scale_color_manual(values = cl[1:2]) +
    theme_classic() +
    theme(axis.text = element_text(size = 12),
          axis.title = element_text(size = 14),
          strip.text = element_text(size = 14),
          legend.text = element_text(size = 12),
          plot.title = element_text(size = 14, hjust = .5),
          legend.position = 'none') +
    geom_text(aes(x = 1.5, y = 1.85), label = stats_expr_vl, 
              check_overlap = TRUE, size = 5, parse = TRUE)

## fig comparing HAQER and HAR like sequences
p_comp <- cs_dat %>% 
    select(-matches('scaled')) %>%
    pivot_longer(cols = matches('sequence_similarity')) %>% 
    mutate(name = case_when(name == 'HAQER_sequence_similarity' ~ 'HAQER-like',
                            name == 'HAR_sequence_similarity' ~ 'HAR-like')) %>%
    ggplot(aes(x = name, y = value / 100)) +
    geom_violin(linewidth = 0, aes(fill = as.logical(vocal_learner_binary)), alpha = 0.5) +
    xlab('Sequence class') +
    ylab('Sequence similarity') +
    scale_y_continuous(labels = scales::percent) +
    theme_classic() +
    theme(axis.text = element_text(size = 12),
          axis.title = element_text(size = 14),
          strip.text = element_text(size = 14),
          legend.text = element_text(size = 12),
          plot.title = element_text(size = 14, hjust = .5),
          legend.position = 'bottom') +
    labs(fill = 'Vocal learner (Wirthlin, et al.):') +
    scale_fill_manual(values = cl[1:2])

## ------------------------------------------------------------------------------------------
## phylogenetic regression predicting log(brain size) using HAQER sequence identity
## ------------------------------------------------------------------------------------------
## fit phylogenetic regression model
fit = phylolm(log10(brain_mass_g) ~ HAQER_sequence_similarity_scaled, phy = tr, data = cs_dat)
summary(fit)

## pull stats for figure
b <- fit$coefficients['HAQER_sequence_similarity_scaled']
b_plot <- unname(round(b, digits = 2))
p <- summary(fit)$coefficients['HAQER_sequence_similarity_scaled','p.value']
p_plot <- format(round(p, 7), nsmall = 3)
p_bs_lab <- bquote(beta == .(b_plot) ~ " " ~ italic(p) == .(p_plot))

b_part <- bquote(beta == .(b_plot))
p_parts <- strsplit(format(p, scientific = TRUE, digits = 2), "e")[[1]]
base <- as.numeric(p_parts[1])
exponent <- as.numeric(p_parts[2])
# p_text <- bquote(italic(p) == .(round(base, 2)) %*% 10^.(exponent))
p_text = bquote(italic(p) == .(formatC(p,2)))
# Combine all expressions
stats_expr <- as.expression(bquote(.(b_part) ~ "" ~ .(p_text)))

## gather stats
n_bs <- sum(is.na(cs_dat$brain_mass_g) == FALSE)

res_bs <- as.data.frame(summary(fit)$coefficients) %>% 
    rownames_to_column(var = 'x') %>% 
    mutate(y = 'log10(brain_mass_g)') %>% 
    as_tibble() %>% 
    relocate(y) %>% 
    rename(beta = Estimate, std.error = StdErr, statistic = t.value) %>% 
    filter(x == 'HAQER_sequence_similarity_scaled') %>% 
    mutate(n_species = n_bs,
           n_vocal_learning_species = NA,
           n_non_vocal_learning_species = NA)

## make figure for PGLS stats of brain size ~ HAQER similarity
p_bs <- cs_dat %>%
    mutate(common_name = case_when(str_detect(common_name, pattern = 'Egyptian fruit bat| red bat|ottlenose dolphin|Brazilian gu|tree shrew|Sumatran rhino|Cuvier|Orca|Killer|dog|Domestic cat') ~ common_name,
                                   TRUE ~ NA_character_)) %>%
    ggplot(aes(x = log10(brain_mass_g), y = HAQER_sequence_similarity_scaled)) +
    # geom_point(aes(color = as.logical(vocal_learner_binary))) +
    geom_point() +
    geom_smooth(method = 'lm', size = 1.5) +
    # ggrepel::geom_text_repel(aes(label = common_name), size = 3.25, max.overlaps = 3) +
    xlab('Log10(brain size)') +
    ylab('HAQER-like sequence score') +
    scale_color_manual(values = cl[1:2]) +
    theme_classic() +
    theme(axis.text = element_text(size = 12),
          axis.title = element_text(size = 14),
          strip.text = element_text(size = 14),
          legend.text = element_text(size = 12),
          plot.title = element_text(size = 14, hjust = .5),
          legend.position = 'bottom') +
    # labs(color = 'Vocal learner (Wirthlin et al.):') +          
    geom_text(aes(x = 1.525, y = 1.85), label = stats_expr, 
              check_overlap = TRUE, size = 5, parse = TRUE)


## ------------------------------------------------------------------------------------------
## phylogenetic regression predicting birth:adult weight ratio using HAQER sequence identity
## ------------------------------------------------------------------------------------------
## fit phylogenetic regression model
fit = phylolm(obstetric_risk ~ HAQER_sequence_similarity_scaled, phy = tr, data = cs_dat)
summary(fit)

## pull stats for figure
b <- fit$coefficients['HAQER_sequence_similarity_scaled']
b_plot <- unname(round(b, digits = 2))
p <- summary(fit)$coefficients['HAQER_sequence_similarity_scaled','p.value']
p_plot <- format(round(p, 7), nsmall = 3)
p_wt_lab <- bquote(beta == .(b_plot) ~ " " ~ italic(p) == .(p_plot))

b_part <- bquote(beta == .(b_plot))
p_parts <- strsplit(format(p, scientific = TRUE, digits = 2), "e")[[1]]
base <- as.numeric(p_parts[1])
exponent <- as.numeric(p_parts[2])
# p_text <- bquote(italic(p) == .(round(base, 2)) %*% 10^.(exponent))
p_text = bquote(italic(p) == .(formatC(p,2)))
# Combine all expressions
stats_expr <- as.expression(bquote(.(b_part) ~ "" ~ .(p_text)))

## gather stats
n_wt <- sum(is.na(cs_dat$obstetric_risk) == FALSE)

res_wt <- as.data.frame(summary(fit)$coefficients) %>% 
    rownames_to_column(var = 'x') %>% 
    mutate(y = 'birth_weight_to_adult_weight_ratio') %>% 
    as_tibble() %>% 
    relocate(y) %>% 
    rename(beta = Estimate, std.error = StdErr, statistic = t.value) %>% 
    filter(x == 'HAQER_sequence_similarity_scaled') %>% 
    mutate(n_species = n_wt,
           n_vocal_learning_species = NA,
           n_non_vocal_learning_species = NA)

## make figure for PGLS of birth:adult weight ratio ~ HAQER similarity
p_wt <- cs_dat %>%
    mutate(common_name = case_when(str_detect(common_name, pattern = 'Egyptian fruit bat| red bat|ottlenose dolphin|Brazilian gu|tree shrew|Sumatran rhino|Cuvier|Orca|Killer|dog|Domestic cat') ~ common_name,
                                   TRUE ~ NA_character_)) %>%
    ggplot(aes(x = obstetric_risk, y = HAQER_sequence_similarity_scaled)) +
    geom_point() +
    geom_smooth(method = 'lm', size = 1.5) +
    xlab('Obstetric risk proxy\n(birth:adult weight ratio)') +
    ylab('HAQER-like sequence score') +
    scale_color_manual(values = cl[1:2]) +
    theme_classic() +
    theme(axis.text = element_text(size = 12),
          axis.title = element_text(size = 14),
          strip.text = element_text(size = 14),
          legend.text = element_text(size = 12),
          plot.title = element_text(size = 14, hjust = .5),
          legend.position = 'bottom') +
    geom_text(aes(x = -1, y = 1.85), label = stats_expr, 
              check_overlap = TRUE, size = 5, parse = TRUE)

## ---------------------------------
## save all stats
## ---------------------------------
bind_rows(res_vl, res_bs, res_wt) %>% 
    write_csv('manuscript/supplemental_materials/stats/cross_species_HAQER_results.csv')

## ---------------------------------
## save figures
## ---------------------------------
p_vl %>% 
    ggsave(filename = 'manuscript/figures/cross_species_vocal_learners_HAQER-like_no_primates_authors.png',
           device = 'png', dpi = 300, bg = 'white',
           units = 'in', width = 5, height = 5)
p_bs %>% 
    ggsave(filename = 'manuscript/figures/cross_species_brain_size_HAQER-like_no_primates_authors.png',
           device = 'png', dpi = 300, bg = 'white',
           units = 'in', width = 5, height = 5)           
p_wt %>% 
    ggsave(filename = 'manuscript/figures/cross_species_birth_weight_ratio_HAQER_and_HAR-like_no_primates_authors.png',
           device = 'png', dpi = 300, bg = 'white',
           units = 'in', width = 5, height = 5)
p_comp %>% 
    ggsave(filename = 'manuscript/figures/cross_species_vocal_learners_HAQER_and_HAR-like_no_primates_authors.png',
           device = 'png', dpi = 300, bg = 'white',
           units = 'in', width = 5, height = 5)

## save figure object
p_vl %>% 
    write_rds('manuscript/figures/R_plot_objects/vocal_learners_HAQER-like_no_primates.rds')
p_comp %>% 
    write_rds('manuscript/figures/R_plot_objects/vocal_learners_HAQER_and_HAR-like_no_primates.rds')
p_bs %>% 
    write_rds('manuscript/figures/R_plot_objects/brain_size_HAQER-like_no_primates.rds')
p_wt %>% 
    write_rds('manuscript/figures/R_plot_objects/birth_weight_ratio_HAQER-like_no_primates.rds')    