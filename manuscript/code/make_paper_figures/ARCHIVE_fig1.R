
library(tidyverse)
library(tidymodels)
library(patchwork)
library(GGally)
library(gridExtra)

setwd('/wdata/common/SLI_WGS/public')
# prs_dat <- read_csv(here("prs/prs_n_threshold_data.csv"))
pheno_dat <- read_csv("phenotype/factors/pheno_factors.csv")
pheno_resid_dat <- read_csv("phenotype/factors/pheno_factors_resid.csv")
fact_dat <- read_rds("phenotype/factors/factors_1_to_10.rds")$coef[[7]]



## Figure 1 Drafts
### Factor loading
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
    guide = guide_colorbar(ticks = FALSE)) &  # barwidth = 1/2, barheight = 12, 
  theme(
    panel.spacing = unit(-1, 'pt'), 
    panel.grid = element_blank(), 
    panel.border = element_rect(color = 'grey30', size =1/3, fill = NA), 
    strip.placement = 'outside', 
    strip.text.y.left = element_text(angle = 0, hjust = 1),
    axis.text.x = element_text(angle = 0, hjust = 1/2, vjust = 0, face = 'bold'), 
    axis.title = element_blank(),
    legend.position = 'bottom',
    legend.key.width = unit(1.75, "cm")
    )


table_a <- tribble(
  ~factor_name, ~interpretation, 
  'F1', 'Verbal memory',
  'F2', 'Receptive language', 
  'F3', 'Nonverbal VIQ', 
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
plot(table_a)
tmp <- tribble(~Factor, ~Description, 
  'F1', 'Core language',
  'F2', 'Receptive language', 
  'F3', 'Nonverbal VIQ', 
  'F4', 'Early language', 
  'F5', 'Talkativeness', 
  'F6', 'Following directions', 
  'F7', 'Vocabulary')
tab <- tableGrob(tmp,  rows = NULL)
plot(tab)



data <- mtcars[, 1:4]

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

plot_b <- pheno_dat %>%
  rename_all(~str_remove(.x, 'actor')) %>%
  select(-sample, -id) %>%
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
  theme(
    axis.text = element_blank(), 
    axis.ticks = element_blank(), 
    #panel.background = element_rect(fill = 'grey99', color = 'grey95'), 
    strip.background = element_blank(), 
    panel.spacing = unit(1/4, 'lines')
  )

plot_b

####################
## behavioral problems
beh_cor <- read_csv('/wdata/lcasten/sli_wgs/behavior/cbcl_correlation_heatmap_ggplot_data.csv')
plot_c <- beh_cor %>% 
  arrange(desc(name)) %>% 
  mutate(name = factor(name, levels = unique(name))) %>%
  arrange(desc(name), estimate) %>% 
  mutate(clean_name = str_replace_all(clean_name, pattern = 'anxiousdepressed', replacement = 'anxious/depressed')) %>%
  mutate(clean_name = factor(clean_name, levels = unique(clean_name)))  %>%
  ggplot(aes(y = name, x = clean_name, fill = estimate)) +
  geom_tile() +
  geom_text(aes(label = sig), size = 8) +
  scale_fill_gradient2(low = 'dodgerblue', mid = 'white', high = 'chocolate1', midpoint = 0) +
  ylab('Factor') +
  xlab('CBCL subscale') +
  labs(fill = 'Correlation:') +
  theme(axis.ticks = element_blank(), 
        #panel.background = element_rect(fill = 'grey99', color = 'grey95'), 
        strip.background = element_blank(), 
        panel.background = element_blank(),
        panel.spacing = unit(1/4, 'lines'),
        legend.key.size = unit(0.75, "cm"),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = 'bottom')

unique(beh_cor$clean_name)
plot_c <- beh_cor %>% 
  arrange(name) %>% 
  mutate(name = factor(name, levels = unique(name))) %>%
  arrange(name, estimate) %>% 
  mutate(clean_name = str_replace_all(clean_name, pattern = 'anxiousdepressed', replacement = 'anxious/depressed'),
         clean_name = str_replace_all(clean_name, pattern = 'externalizing', replacement = 'externalizing problems'),
         clean_name = str_replace_all(clean_name, pattern = 'internalizing', replacement = 'internalizing problems'),
         clean_name = str_replace_all(clean_name, pattern = 'attention scale', replacement = 'attention problems'),
         clean_name = str_replace_all(clean_name, pattern = ' behavior scale', replacement = ' behavior'),
         clean_name = str_replace_all(clean_name, pattern = ' complaints scale', replacement = ' complaints'),
         clean_name = str_replace_all(clean_name, pattern = ' problems scale', replacement = ' problems')
         ) %>%
  mutate(clean_name = factor(clean_name, levels = unique(clean_name)))  %>%
  ggplot(aes(x = name, y = clean_name, fill = estimate)) +
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


plot_c

## assemble
layout <- c(
  area(1, l = 1, r = 20), area(1, l = 24, r = 24),  area(1, l = 30, r = 100), 
  area(2, l = 1, r = 20, b =3), area(2, l = 21, r = 100, b = 3)
  )

layout <- c(
  area(1, l = 1, r = 20), area(1, l = 24, r = 24),  area(1, l = 30, r = 100), 
  area(2, l = 21, r = 100, b = 2)
  )

layout <- "AAAAABBBBB
CCCCCDDDDD
"

mytheme <- gridExtra::ttheme_default(
    core = list(fg_params=list(cex = 1.75)),
    colhead = list(fg_params=list(cex = 1.75)),
    rowhead = list(fg_params=list(cex = 1.75)))
tab2 = tableGrob(tmp,  rows = NULL, theme = mytheme)
# plot(tab2)

## merge together with patchwork
# plot_fig1_v2 <- (plot_a +  plot_c) / GGally::ggmatrix_gtable(plot_b) +  plot_layout(widths = c(2, 2, 1), nrow = 2) + plot_annotation(tag_levels = 'A') & theme(plot.tag = element_text(size = 12))
library(GGally)
plot_fig1_v2 <- plot_a +  plot_c + wrap_elements(tab2) + ggmatrix_gtable(plot_b) +  plot_layout(design = layout) + plot_annotation(tag_levels = 'A') & theme(plot.tag = element_text(size = 12))
# plot_fig1_v2
# dir.create('/wdata/lcasten/sli_wgs/paper_figures')
ggsave(plot_fig1_v2, filename = '/wdata/lcasten/sli_wgs/paper_figures/fig1_fact_plus_v2.png', width = 12, height = 12, dpi = 300)

plot(tab2)
ggsave(plot_a, filename = '/wdata/lcasten/sli_wgs/paper_figures/subplots/fig1A.png', width = 6, height = 7, dpi = 300, bg = 'white')
ggsave(plot_c, filename = '/wdata/lcasten/sli_wgs/paper_figures/subplots/fig1B.png', width = 6, height = 7, dpi = 300, bg = 'white')








# plot_prs_plus <- wrap_plots(
#   list(
#     plot_a, 
#     guide_area(),
#     GGally::ggmatrix_gtable(plot_b), 
#     plot_c
# ),
#    design = layout,
#   heights = c(6,5), 
#   guides = 'collect', 
#   tag_level = 'new'
# )   +
#   plot_annotation(tag_levels = 'A', title = 'PRSs Positive') & 
#   theme(
#     axis.title = element_blank()
#     ) 
# plot_prs_plus