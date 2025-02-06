library(tidyverse)

## helper function for correlations in long dataframes
long_corr <- function(x, y) {
    corr_res <- x %>%
        inner_join(y) %>% 
        group_by(factor, name) %>% 
        do(res = broom::tidy(cor.test(.$factor_val, .$value))) %>% 
        unnest(res) %>% 
        arrange(factor, p.value) %>% 
        group_by(factor) %>% 
        mutate(fdr = p.adjust(p.value, method = 'fdr')) %>%  ## FDR correction
        ungroup() %>% 
        relocate(fdr, .after = p.value) %>% 
        mutate(sig = case_when(fdr < 0.05 ~ '**',
                               fdr >= 0.05 & p.value < 0.05 ~ '*',
                               p.value >= 0.05 ~ '')
            )
    return(corr_res)
}

#######################
## read in EpiSLI data
#######################
ph <- read_csv('manuscript/supplemental_materials/EpiSLI_pheno_data.csv')
pgs <- read_csv('manuscript/supplemental_materials/EpiSLI_LDpred2-PGS_data.csv')
es_pgs <- read_csv('manuscript/supplemental_materials/EpiSLI_ES-PGS_data.csv')
fact_dat <- read_rds("manuscript/supplemental_materials/EpiSLI_factor_analysis_data.rds")$coef[[7]] ## read in 7 factor solution loadings

#######################################
## association between factors + CBCL
#######################################
## subset to factor scores and reformat to "long" df
fc <- ph %>% 
    select(IID, matches('^F[1-7]')) %>% 
    pivot_longer(cols = -IID, names_to = 'factor', values_to = 'factor_val')
## subset to CBCL scores and reformat to "long" df
cbcl <- ph %>% 
    select(IID, matches('^cbcl_')) %>% 
    pivot_longer(cols = -IID)
## merge "long" data to correlate all possible Factor x CBCL subscales 
cbcl_cor <- long_corr(x = fc, y = cbcl)

## save correlation sumstats
cbcl_cor %>% 
    write_csv('manuscript/supplemental_materials/stats/EpiSLI_factor_CBCL_correlations.csv')

## make heatmap of Pearson's r
p_cbcl <- cbcl_cor %>%
    mutate(clean_name = str_remove_all(name, pattern = '^cbcl_'),
           clean_name = str_replace_all(clean_name, pattern = '_', replacement = ' ')) %>%
    arrange(factor) %>% 
    mutate(factor = factor(factor, levels = unique(factor))) %>%
    arrange(factor, estimate) %>% 
    mutate(clean_name = str_replace_all(clean_name, pattern = 'anxiousdepressed', replacement = 'anxious / depressed'),
            clean_name = str_replace_all(clean_name, pattern = 'externalizing', replacement = 'externalizing problems'),
            clean_name = str_replace_all(clean_name, pattern = 'internalizing', replacement = 'internalizing problems'),
            clean_name = str_replace_all(clean_name, pattern = 'attention scale', replacement = 'attention problems'),
            clean_name = str_replace_all(clean_name, pattern = ' behavior scale', replacement = ' behavior'),
            clean_name = str_replace_all(clean_name, pattern = ' complaints scale', replacement = ' complaints'),
            clean_name = str_replace_all(clean_name, pattern = ' problems scale', replacement = ' problems'),
            clean_name = str_replace_all(clean_name, pattern = ' compentence', replacement = ' competence')
            ) %>%
    mutate(clean_name = factor(clean_name, levels = unique(clean_name)))  %>%
    ggplot(aes(x = factor, y = clean_name, fill = estimate)) +
    geom_tile() +
    geom_text(aes(label = sig), size = 5, hjust=0.5) +
    scale_fill_gradient2(low = 'dodgerblue', mid = 'white', high = 'chocolate1', midpoint = 0) +
    xlab(NULL) +
    ylab('CBCL subscale') +
    labs(fill = 'Correlation:') +
    theme_minimal() +  
    scale_x_discrete(expand = c(0,0), position = 'top') + 
    theme(
          panel.spacing = unit(-1, 'pt'), 
          panel.grid = element_blank(), 
          panel.border = element_rect(color = 'grey30', size =1/3, fill = NA),
          strip.placement = 'outside',
          strip.text.y.left = element_text(angle = 0, hjust = 1),
          axis.text.x = element_text(angle = 0, hjust = 1/2, vjust = 0, face = 'bold'),
          axis.title = element_blank(),
          legend.key.width = unit(1, "cm"),
          legend.position = 'bottom'
          )
p_cbcl %>% 
    ggsave(filename = 'manuscript/figures/EpiSLI_factor_CBCL_heatmap.png', 
           device = 'png', dpi = 300, bg = 'white', 
           units = 'in', width = 5, height = 7)


##################################################
## association between factors + genome wide PGS
##################################################
## subset to PGS of interest
pgs_str <- "functional|BIG5|^epilepsy$|insomn|cogntive|perf|PGC|schiz|bip|ADHD|leftHand|height|autism|neurodev|insomn|Cognitiv|trauma|osychiatric|antisoc|aggress|addiction|educat|pgc|execut|vocal"
pgs2 <- pgs %>% 
    select(IID, name = pgs_name, value = pgs_pc_corrected) %>% 
    filter(str_detect(name, pattern = pgs_str)) %>% 
    filter(! name %in% c('ADHD_diagnosed_late', 'schoolE1.better_overall_school_performance', 'PGC2022-ADHD', 'height_Yengo'))
## run correlations
pgs_cor <- long_corr(x = fc, y = pgs2)
## add clean PGS names and categories
pgs_map <- read_csv('manuscript/supplemental_materials/PGS_name_map.csv')
pgs_cor <- pgs_map %>% 
    inner_join(pgs_cor)
## save sumstats
pgs_cor %>% 
    write_csv('manuscript/supplemental_materials/stats/EpiSLI_factor_PGS_correlations.csv')
## heatmap of factor x PGS correlations
p_pgs <- pgs_cor %>% 
    mutate(type = ifelse(type == 'Cognitive', 'Cog.', type)) %>%
    mutate(type = factor(type, levels= c('Cog.', 'Neuropsychiatric', 'Personality', 'Brain MRI', 'Other'))) %>%
    filter(type != 'Brain MRI') %>% ## remove MRI from figure since it is too busy and doesn't offer much anyway
    filter(! type %in% c('Alcohol dependency', 'Cannabis use disorder', 'Antisocial behavior')) %>% ## remove some redundant PGS from figure
    arrange(factor, estimate) %>%
    arrange(desc(clean_name)) %>% 
    mutate(clean_name = factor(clean_name, levels = unique(clean_name))) %>%
    ggplot(aes(x = factor, y = clean_name, fill = estimate)) +
    geom_tile() +
    geom_text(aes(label = sig), vjust = 0.7, hjust = 0.5, size = 5) +
    facet_grid(rows = vars(type), scales = 'free_y', space = 'free_y') +
    scale_fill_gradient2(low = 'dodgerblue', mid = 'white', high = 'chocolate1', midpoint = 0) +
    xlab(NULL) +
    ylab('Polygenic score') +
    labs(fill = 'Correlation:') +
    scale_x_discrete(expand = c(0,0), position = 'top') + 
    theme_classic() +  
    theme(legend.key.height = unit(1, "cm"),
          panel.border = element_rect(color = "black", fill = NA, size = 1),
          strip.background = element_rect(colour = "black", fill="white", size = 1),
          strip.placement = 'outside', 
          legend.position = 'right',
          axis.text = element_text(size = 12),
          strip.text = element_text(size = 14),
          axis.title = element_text(size = 14),
          legend.text = element_text(size = 12),
          legend.title = element_text(size = 14))

p_pgs %>% 
    ggsave(filename = 'manuscript/figures/EpiSLI_factor_PGS_heatmap.png', 
           device = 'png', dpi = 300, bg = 'white', 
           units = 'in', width = 8, height = 14)


#########################
## ES-PGS analysis
#########################
## ES-PGS modelling
coi <- names(es_pgs)[str_detect(names(es_pgs), pattern = 'pgs') & str_detect(names(es_pgs), 'complement', negate = TRUE)]
tmp <- fc %>% 
  inner_join(es_pgs)

iter = 0
res_list = list()
for(evo in coi){
  cat(sprintf('\n\n\n\n'))
  message(evo)
  for(i in 1:7){
      message('Factor ', i)
      iter = iter + 1
      tmp2 <- tmp %>% 
        filter(factor == str_c('F', i)) %>% 
        select(IID, factor_val, matches(str_c(str_remove_all(evo, pattern = 'cp_pgs.'), '$')))
      bdat <- tmp2 %>% 
        select(IID, factor_val, matches('complement'))

      baseline <- lm(factor_val ~ . - IID, data = bdat)
      baseline_plus_anno <- lm(factor_val ~ . - IID, data = tmp2)
      baseline_rsq = summary(baseline)$r.squared
      baseline_plus_anno_rsq = summary(baseline_plus_anno)$r.squared
      baseline_plus_anno_coef <- broom::tidy(baseline_plus_anno) %>% 
        filter(term != '(Intercept)' & str_detect(term, 'complement', negate = TRUE))

      res <- broom::tidy(anova(baseline, baseline_plus_anno)) %>% 
        drop_na() %>% 
        mutate(factor = str_c('F', i)) %>% 
        mutate(evo = str_remove_all(evo, 'cp_pgs.')) %>% 
        relocate(factor, evo) %>% 
        mutate(baseline_rsq = baseline_rsq,
               baseline_plus_anno_rsq = baseline_plus_anno_rsq,
               annotation_beta = baseline_plus_anno_coef$estimate,
               annotation_std_err = baseline_plus_anno_coef$std.error,
               annotation_pval = baseline_plus_anno_coef$p.value)
      res_list[[iter]] <- res
  }
}

## get SNP info
snp_counts <- read_csv('manuscript/supplemental_materials/EpiSLI_ES-PGS_PRSet_independent_SNP_counts.csv')

## reformat ES-PGS results and merge with number of SNPs used in calculation
es_pgs_res <- bind_rows(res_list) %>%
  arrange(p.value) %>% 
  select(-term) %>% 
  rename(model = evo, p.value_model_comparison = p.value) %>%
  inner_join(snp_counts) %>% 
  arrange(factor, model) %>%
  relocate(complement_PGS_n_snp, .before = annotation_PGS_n_snp) %>%
  mutate(rsq_gain_per_1000indSNP = (baseline_plus_anno_rsq - baseline_rsq) / (annotation_PGS_n_snp / 1000))
es_pgs_res %>% 
  write_csv('manuscript/supplemental_materials/stats/EpiSLI_factor_ES-PGS_results.csv')

## make figures for ES-PGS
p_es_pgs_forest <- es_pgs_res %>% 
    group_by(factor) %>%
    filter(factor == 'F1') %>% 
    mutate(mod_clean = case_when(model == 'LinAR_Catarrhini' ~ 'Catarrhini',
                                 model == 'LinAR_Hominidae' ~ 'Great ape acceleration',
                                 model == 'LinAR_Homininae' ~ 'Homininae',
                                 model == 'LinAR_Hominoidea' ~ 'Hominoidea',
                                 model == 'LinAR_Simiformes' ~ 'Simiformes',
                                 model == 'NeanderthalSelectiveSweep' ~ 'Neanderthal deserts',
                                 model == 'consPrimates_UCE' ~ 'Primate UCEs',
                                 model == 'human_chimp_div_DMG' ~ 'Human-chimp divergence',
                                 model == 'human_singleton_density_score_top5pct' ~ 'Recent selection',
                                 model == 'HAQER' ~ 'HAQERs',
                                 model == 'HAR' ~ 'HARs',
                                 TRUE ~ NA_character_)) %>% 
    drop_na(mod_clean) %>%
    mutate(mod_clean = factor(mod_clean, levels = c('Primate UCEs', 'Simiformes', 'Catarrhini', 'Hominoidea', 'Great ape acceleration', 'Homininae', 'Human-chimp divergence', 'HAQERs','HARs',  'Neanderthal deserts', 'Recent selection'))) %>%
    filter(mod_clean %in% c('Primate UCEs', 'Great ape acceleration', 'Human-chimp divergence', 'HAQERs','HARs', 'Neanderthal deserts', 'Recent selection')) %>%
    mutate(sig = ifelse(p.value_model_comparison < 0.05, TRUE, FALSE)) %>%
    ggplot(aes(x = mod_clean, y = annotation_beta, color = sig)) +
    geom_linerange(aes(ymin = annotation_beta - 1.96 * annotation_std_err, ymax = annotation_beta + 1.96 * annotation_std_err), size = 1.2) +
    geom_point(size = 3.5, aes(shape = sig)) +
    scale_shape_manual(values = c(1, 16)) +
    geom_hline(yintercept = 0, color = 'red', linetype = 'dashed', size = 1.075) +
    xlab('ES-PGS model') +
    ylab('Effect on core language (F1)') +
    scale_color_manual(values = c('grey75', 'black')) +
    theme_classic() +  
    theme(axis.text = element_text(size = 12),
         axis.text.x = element_text(angle = 15, hjust = 1),
          axis.title = element_text(size = 14),
          legend.text = element_text(size = 12),
          legend.title = element_text(size = 14),
          legend.position = 'none') +
    guides(shape = 'none')
p_es_pgs_forest %>%
    ggsave(filename = 'manuscript/figures/EpiSLI_factor_ES-PGS_forest.png', 
           device = 'png', dpi = 300, bg = 'white', 
           units = 'in', width = 8, height = 7)

## scatterplots of residualized F1-F3 and HAQERs
wd <- tmp %>% 
    select(IID, factor, factor_val, matches('HAQER')) %>% 
    filter(factor %in% str_c('F', 1:3)) %>% 
    group_by(factor) %>% 
    mutate(factor_val_resid = resid(lm(factor_val ~ cp_pgs.complement_HAQER)),
           factor_val_resid = scale(factor_val_resid)[,1])
es_pgs_cor <- wd %>% 
    group_by(factor) %>% 
    do(res = broom::tidy(cor.test(.$factor_val_resid, .$cp_pgs.HAQER))) %>% 
    unnest(res) %>% 
    mutate(lab = str_c('r = ', round(estimate, digits = 2), ', p-val = ', formatC(p.value, digits = 2)))
p_f1 <- wd %>% 
    inner_join(es_pgs_cor) %>%    
    filter(factor == 'F1') %>%
    ggplot(aes(x = factor_val_resid, y = cp_pgs.HAQER)) +
    geom_point(size = 0.7) +
    geom_smooth(method = 'lm', size = 1.5, color = 'black') +
    xlab('Core language (F1) adjusted\nfor background CP-PGS') +
    ylab('HAQER CP-PGS') +
    geom_text(aes(x = -1.5, y = 2.85, label = lab), size = 4, check_overlap = TRUE) +
    theme_classic() +  
    theme(axis.text = element_text(size = 12),
          axis.title = element_text(size = 14),
          legend.text = element_text(size = 12),
          legend.title = element_text(size = 14))
p_f2 <- wd %>% 
    inner_join(es_pgs_cor) %>%    
    filter(factor == 'F2') %>%
    ggplot(aes(x = factor_val_resid, y = cp_pgs.HAQER)) +
    geom_point(size = 0.7) +
    geom_smooth(method = 'lm', size = 1.5, color = 'black') +
    xlab('Receptive language (F2) adjusted\nfor background CP-PGS') +
    ylab('HAQER CP-PGS') +
    geom_text(aes(x = -1.5, y = 2.85, label = lab), size = 4, check_overlap = TRUE) +
    theme_classic() +  
    theme(axis.text = element_text(size = 12),
          axis.title = element_text(size = 14),
          legend.text = element_text(size = 12),
          legend.title = element_text(size = 14))
p_f3 <- wd %>% 
    inner_join(es_pgs_cor) %>%    
    filter(factor == 'F3') %>%
    ggplot(aes(x = factor_val_resid, y = cp_pgs.HAQER)) +
    geom_point(size = 0.7, color = 'grey75') +
    geom_smooth(method = 'lm', size = 1.5, color = 'grey60') +
    xlab('Nonverbal IQ (F3) adjusted\nfor background CP-PGS') +
    ylab('HAQER CP-PGS') +
    geom_text(aes(x = -1.5, y = 2.85, label = lab), size = 4, check_overlap = TRUE) +
    theme_classic() +  
    theme(axis.text = element_text(size = 12),
          axis.title = element_text(size = 14),
          legend.text = element_text(size = 12),
          legend.title = element_text(size = 14))
p_f1 %>% 
    ggsave(filename = 'manuscript/figures/EpiSLI_factor1_HAQER_ES-PGS.png', 
           device = 'png', dpi = 300, bg = 'white', 
           units = 'in', width = 5, height = 5)
p_f2 %>% 
    ggsave(filename = 'manuscript/figures/EpiSLI_factor2_HAQER_ES-PGS.png', 
           device = 'png', dpi = 300, bg = 'white', 
           units = 'in', width = 5, height = 5)
p_f3 %>% 
    ggsave(filename = 'manuscript/figures/EpiSLI_factor3_HAQER_ES-PGS.png', 
           device = 'png', dpi = 300, bg = 'white', 
           units = 'in', width = 5, height = 5)

###############################
## make factor loadings fig
###############################
## loading heatmap
plot_a <- fact_dat %>%
  as_tibble(rownames = 'assessment_name') %>%
  separate(assessment_name, 
    into = c('grade', 'category', 'assessment_name'), 
    extra = 'merge') %>%
  mutate(assessment_name = str_replace_all(assessment_name, pattern = 'vocab', replacement = 'vocabulary')) %>%
  mutate(assessment_name = str_replace_all(assessment_name, pattern = 'iq', replacement = 'IQ')) %>%
  mutate(assessment_name = str_replace_all(assessment_name, pattern = '_', replacement = ' ')) %>% 
  mutate(assessment_name = as_factor(assessment_name) %>% 
  fct_reorder2(assessment_name, category, first2)) %>%
  gather(factor_name, factor_loading, matches("Factor\\d")) %>%
  mutate(factor_name = str_remove(factor_name, 'actor')) %>%
  ggplot(aes(y = grade, x = factor_name, fill = factor_loading)) + 
  geom_tile() + 
  ggh4x::facet_nested(
    rows = vars(assessment_name),
    scales = 'free',
    space = 'free',
    nest_line = TRUE,
    switch = 'y',
    resect = unit(1/4, 'lines')
    ) +
  theme_minimal() +  
  scale_x_discrete(expand = c(0,0), position = 'top') + 
  scale_y_discrete(expand = c(0,0), position = 'right') + 
  scale_fill_fermenter(palette = 'RdGy', 
    limits = c(-1, 1), 
    breaks = c(-5:-1,-1/4,1/4, 1:5)/5, 
    labels = c('-1', '', '-0.6', '', '-0.2', '', '',  '0.2',  '',  '0.6', '',  '1'),
    name = 'Loadings:', 
    guide = guide_colorbar(ticks = FALSE)) +  # barwidth = 1/2, barheight = 12, 
  theme(
    panel.spacing = unit(-1, 'pt'), 
    panel.grid = element_blank(), 
    panel.border = element_rect(color = 'grey30', size =1/3, fill = NA), 
    strip.placement = 'outside', 
    strip.text.y.left = element_text(angle = 0, hjust = 1),
    axis.text = element_text(size = 12),
    axis.text.x = element_text(angle = 0, hjust = 1/2, vjust = 0), 
    axis.title = element_blank(),
    legend.position = 'bottom',
    legend.key.width = unit(1.25, "cm")
    )
plot_a %>% 
    ggsave(filename = 'manuscript/figures/EpiSLI_factor_loadings.png', 
           device = 'png', dpi = 300, bg = 'white', 
           units = 'in', width = 5, height = 8)  

## table with human readable names
table_a <- tribble(
  ~factor_name, ~interpretation, 
  'F1', 'Core language',
  'F2', 'Receptive language', 
  'F3', 'Nonverbal IQ', 
  'F4', 'Early language', 
  'F5', 'Talkativeness', 
  'F6', 'Following directions', 
  'F7', 'Vocabulary'
  ) %>% 
  gridExtra::tableGrob(
    rows = NULL,
    cols = NULL, 
    theme = gridExtra::ttheme_minimal(
      core = list(
        bg_params = list(fill = c("grey99", "grey96"), lwd = 1.5, col = "white"),
        fg_params=list(hjust=0, x=0)),
      colhead = list(fg_params=list(hjust=0, x=0))
      )
    )
pdf('manuscript/figures/EpiSLI_factor_description.pdf')
plot(table_a)
dev.off()

###############################################
## figure with factor intercorrelations and distributions
library(patchwork)
library(GGally)
library(gridExtra)
# Create the ggmatrix
ggpairs(
  data,
  lower = list(continuous = function(data, mapping, ...) {
    ggally_text(data = data, mapping = mapping, label = paste(mapping$x, "vs", mapping$y), ...)
  }),
  upper = list(continuous = "blank"),     # Optionally, leave the upper triangle blank
  diag = list(continuous = "blankDiag")   # Optionally, leave the diagonal blank
)

### pairs plot
lower_text <- function(data, mapping, ...) {
  lower_dat <- tab
  ggplot(data, mapping) + 
    geom_text(aes(label = paste(..x.., ..y.., sep = ", ")), ...) +
    theme_void()
}

smooth_fn <- function(data, mapping, pts=list(), smt=list(), ...){
              ggplot(data = data, mapping = mapping, ...) + 
                         do.call(geom_point, pts) +
                         do.call(geom_smooth, smt) +
    theme_minimal() + 
    theme(
      panel.grid.minor = element_blank(), 
      panel.grid.major = element_line(size = 1/4)
      )
}

dens_fn <- function(data, mapping, dens = list(), ...){
              ggplot(data = data, mapping = mapping, ...) + 
                         do.call(geom_density, dens) +
    theme_void()
                 }

cor_fn <- function(data, mapping, method, symbol, ...){
  x <- eval_data_col(data, mapping$x)
  y <- eval_data_col(data, mapping$y)

  corr <- cor(x, y, method=method, use='complete.obs')

  colFn <- colorRampPalette(
    RColorBrewer::brewer.pal(name = 'RdBu', n = 9) %>%
      rev(),
    interpolate ='spline')
  fill <- colFn(100)[findInterval(corr, seq(-.2, .2, length = 100))]

  ggally_text(
    label = paste(symbol, as.character(round(corr, 2))), 
    mapping = aes(),
    xP = 0.5, yP = 0.5,
    color = 'grey10',
    ...) + 
    theme_void() +
    theme(panel.background = element_rect(fill = fill))
}

plot_b <- ph %>%
  select(matches('^F[0-9]')) %>%
  GGally::ggpairs(
    progress = F ,
    lower = 'blank',
    # lower = list(
    #   continuous = wrap(smooth_fn, 
    #     pts = list(alpha = .15, shape = 21, size = 1/2, color = 'grey30'), 
    #     smt = list(color = scales::muted('green', l = 40, c = 50), 
    #       alpha = 0, method = 'lm', size = 2/3))
    #   ), 
    upper = list(
      continuous = wrap(cor_fn, method = 'p', symbol = 'r = ', size = 3)
      ), 
    diag = list(
      continuous = wrap(dens_fn,
        mapping = aes(color = name),
        dens = list())
    )
    ) + 
  # theme_classic() + 
  theme(
    axis.text = element_blank(), 
    axis.ticks = element_blank(), 
    #panel.background = element_rect(fill = 'grey99', color = 'grey95'), 
    strip.background = element_blank(), 
    panel.spacing = unit(1/4, 'lines'),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 14)
  )

plot_b


#####################
## save plot objects
plot_a %>% 
  write_rds('manuscript/figures/R_plot_objects/EpiSLI_factor_loadings.rds')
plot_b %>% 
  write_rds('manuscript/figures/R_plot_objects/EpiSLI_factor_correlations.rds')
table_a %>% 
  write_rds('manuscript/figures/R_plot_objects/EpiSLI_factor_names.rds')
p_cbcl %>% 
  write_rds('manuscript/figures/R_plot_objects/EpiSLI_factor_CBCL_correlations.rds')
p_pgs %>% 
  write_rds('manuscript/figures/R_plot_objects/EpiSLI_factor_PGS_correlations.rds')
p_es_pgs_forest %>% 
  write_rds('manuscript/figures/R_plot_objects/EpiSLI_factor_ES-PGS_forest.rds')
p_f1 %>% 
  write_rds('manuscript/figures/R_plot_objects/EpiSLI_factor1_HAQER-CP-PGS.rds')
p_f2 %>% 
  write_rds('manuscript/figures/R_plot_objects/EpiSLI_factor2_HAQER-CP-PGS.rds')
p_f3 %>% 
  write_rds('manuscript/figures/R_plot_objects/EpiSLI_factor3_HAQER-CP-PGS.rds')
