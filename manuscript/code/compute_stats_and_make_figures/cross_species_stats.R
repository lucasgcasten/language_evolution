## load packages
library(tidyverse)
library(ape)
library(phylolm)

## read in phylogeny downloaded from UCSC (https://hgdownload.soe.ucsc.edu/goldenPath/hg38/cactus447way/hg38.447way.nh.txt)
tr <- read.tree('manuscript/supplemental_materials/hg38.447way.nh.txt')

## read in data
cs_dat <- read_csv('manuscript/supplemental_materials/cross_species_vocal_learning_HAQER_similarity.csv')

## z-scale the sequence similarity scores
cs_dat$HAQER_sequence_similarity_scaled <- scale(cs_dat$HAQER_sequence_similarity)[,1]
cs_dat$HAR_sequence_similarity_scaled <- scale(cs_dat$HAR_sequence_similarity)[,1]

## reformat data to work with phylolm
cs_dat <- as.data.frame(cs_dat)
rownames(cs_dat) <- cs_dat$species_name

############################################################################
############################################################################
############## DELETE
cs_dat$bw_ratio <- scale(log10(cs_dat$birth_weight_g) / log10(cs_dat$adult_body_mass_g))[,1]

fit_bw_vl <- phylolm(bw_ratio ~ as.factor(vocal_learner_binary), phy = tr, data = cs_dat)
fit_bw_seq <- phylolm(bw_ratio ~ HAQER_sequence_similarity_scaled + HAR_sequence_similarity_scaled, phy = tr, data = cs_dat)
n_bw <- sum(is.na(cs_dat$bw_ratio) == FALSE)

## HAQER score vs prenatal growth score figure
## pull stats for figure
b <- fit_bw_seq$coefficients['HAQER_sequence_similarity_scaled']
b_plot <- unname(round(b, digits = 2))
p <- summary(fit_bw_seq)$coefficients['HAQER_sequence_similarity_scaled','p.value']
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
stats_expr_bw <- as.expression(bquote(.(b_part) ~ "" ~ .(p_text) ~ "" ~ .(str_c("N = ", n_bw, ' species'))))

## pull stats for prenatal development ~ vl model
b <- fit_bw_vl$coefficients['as.factor(vocal_learner_binary)1']
b_plot <- unname(round(b, digits = 2))
p <- summary(fit_bw_vl)$coefficients['as.factor(vocal_learner_binary)1','p.value']
# p_plot <- formatC(p, digits = 2)
p_plot <- format(round(p, 7), nsmall = 3)
# p_vl_lab <- str_c('beta = ', round(b, digits = 2), '\np = ', formatC(p, digits = 2))
p_vl_lab <- bquote(beta == .(b_plot) ~ " " ~ italic(p) == .(p_plot))

b_part <- bquote(beta == .(b_plot))
p_parts <- strsplit(format(p, scientific = TRUE, digits = 2), "e")[[1]]
base <- as.numeric(p_parts[1])
exponent <- as.numeric(p_parts[2])
p_text <- bquote(italic(p) == .(round(p, 3)))
# Combine all expressions
stats_expr_bw <- as.expression(bquote(.(b_part) ~ "" ~ .(p_text) ~ ""))


cl <- c("#762776", "#e04468", "#dcc699")
p_bw_haq_seq <- cs_dat %>% 
    ggplot(aes(x = HAQER_sequence_similarity_scaled, y = bw_ratio)) +
    geom_point() +
    geom_smooth(method = 'lm', color = "#762776", size = 1.5) +
    xlab("HAQER-like score") +
    ylab("Prenatal vs postnatal growth score") +
    theme_classic() +
    theme(axis.text = element_text(size = 12),
          axis.title = element_text(size = 14),
          strip.text = element_text(size = 14),
          legend.text = element_text(size = 12),
          plot.title = element_text(size = 14, hjust = .5),
          legend.position = 'none') +
    geom_text(aes(x =-.5, y = 2.1), label = stats_expr, 
              check_overlap = TRUE, size = 5, parse = TRUE)

n_voc_bw <- nrow(cs_dat[cs_dat$vocal_learner_binary == 1 & is.na(cs_dat$bw_ratio) == FALSE,])
n_nonvoc_bw <- n_bw - n_voc
set.seed(65666)
cl <- sample(RColorBrewer::brewer.pal(name="Dark2",n=3))
p_bw_vl <- cs_dat %>% 
    drop_na(bw_ratio) %>%
    group_by(vocal_learner_binary) %>% 
    mutate(n = n()) %>%
    ungroup() %>%
    mutate(vlab = as.logical(vocal_learner_binary),
            vlab = str_c(vlab, '\nN = ', n, ' species')
            ) %>%
    mutate(common_name = case_when(str_detect(common_name, pattern = 'Egyptian fruit bat| red bat|ottlenose dolphin|Brazilian gu|tree shrew|Sumatran rhino|Cuvier|Orca|Killer|dog|Domestic cat') ~ common_name,
                                   TRUE ~ NA_character_)) %>%
    ggplot(aes(x = vlab, y = bw_ratio)) +
    geom_violin(linewidth = 0, aes(fill = vlab), alpha = 0.4) +
    geom_point(aes(color = vlab), size = .7) +
    ggrepel::geom_text_repel(aes(label = common_name, color = vlab), size = 3.25, max.overlaps = 15) +
    xlab('Vocal learner (Wirthlin, et al.)') +
    ylab('Prenatal vs postnatal growth score') +
    scale_fill_manual(values = cl[1:2]) +
    scale_color_manual(values = cl[1:2]) +
    theme_classic() +
    theme(axis.text = element_text(size = 12),
          axis.title = element_text(size = 14),
          strip.text = element_text(size = 14),
          legend.text = element_text(size = 12),
          plot.title = element_text(size = 14, hjust = .5),
          legend.position = 'none') +
    geom_text(aes(x = 1.5, y = 1.7), label = stats_expr_bw, 
              check_overlap = TRUE, size = 5, parse = TRUE)
library(patchwork)
p_bw_merged <- p_bw_vl + p_bw_haq_seq + plot_annotation(tag_levels = 'a') & 
        theme(plot.tag = element_text(size = 28, face = 'bold'),
               axis.text = element_text(size = 14),
               axis.title = element_text(size = 16),
               legend.text = element_text(size = 20),
               legend.title = element_text(size = 20),
               strip.text = element_text(size = 20),
               plot.title = element_text(size = 20, face = 'bold'))
p_bw_merged %>% 
    ggsave(filename = "/wdata/lcasten/sli_wgs/cross_species_prenatal.png", device = 'png', dpi = 300, bg = 'white',
           width = 12, height = 6, units = 'in')
############## DELETE
############################################################################
############################################################################


## fit phylogenetic regression model
fit = phyloglm(vocal_learner_binary ~ HAQER_sequence_similarity_scaled + HAR_sequence_similarity_scaled, phy = tr, data = cs_dat, method = 'logistic_IG10')
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
stats_expr <- as.expression(bquote(.(b_part) ~ "" ~ .(p_text)))
 


## gather stats and save
n <- nrow(cs_dat)
n_voc <- nrow(cs_dat[cs_dat$vocal_learner_binary == 1,])
n_nonvoc <- n - n_voc

as.data.frame(summary(fit)$coefficients) %>% 
    rownames_to_column(var = 'x') %>% 
    mutate(y = 'vocal_learner_binary') %>% 
    as_tibble() %>% 
    relocate(y) %>% 
    rename(beta = Estimate, std.error = StdErr, statistic = z.value) %>% 
    filter(x == 'HAQER_sequence_similarity_scaled') %>% 
    mutate(n_species = n,
           n_vocal_learning_species = n_voc,
           n_non_vocal_learning_species = n_nonvoc) %>% 
    write_csv('manuscript/supplemental_materials/stats/cross_species_HAQER_results.csv')

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
    geom_text(aes(x = 1.5, y = 2.1), label = stats_expr, 
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

## save figure
p_vl %>% 
    ggsave(filename = 'manuscript/figures/cross_species_vocal_learners_HAQER-like_no_primates_authors.png',
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
