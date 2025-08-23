library(tidyverse)
library(patchwork)

## read in ABCD imaging data and ES-PGS
abcd_df <- read_csv('manuscript/supplemental_materials/ABCD_data.csv')

######################################################
## ABCD verbal learning and cognitive score analysis
######################################################
abcd_df_cog <- abcd_df %>% 
    select(IID, ravlt_age, sex, ravlt_total_learning, nihtbx_cryst_fc, nihtbx_fluidcomp_fc, nihtbx_reading_fc, nihtbx_picvocab_fc, matches("HAQER_v2"), matches('pc'))

abcd_cog_espgs_stats <- abcd_df_cog %>%
    rename(age = ravlt_age) %>%
    pivot_longer(cols = matches("ravlt|nihtbx")) %>% 
    drop_na() %>%
    group_by(name) %>% 
    mutate(n = n()) %>% 
    group_by(name, n) %>%
    do(res = broom::tidy(lm(value ~ age + as.factor(sex) + pc1 + pc2 + pc3 + pc4 + pc5 + cp_pgs.HAQER_v2 + cp_pgs.background_HAQER_v2 + cp_pgs.matched_control_HAQER_v2, data = .))) %>%
    unnest(res) %>% 
    filter(term == 'cp_pgs.HAQER_v2') %>% 
    arrange(p.value) %>% 
    rename(x = term, beta = estimate, pheno = name) %>%
    relocate(n, .after = p.value)

#################################################
## ABCD cephalo-pelvic disproportion analysis
#################################################
## compute cephalo-pelvic disproportion using residualized ICV (obstetric risk proxy)
abcd_df_cpd <- abcd_df %>% 
    select(IID, mri_age, sex, smri_vol_scs_intracranialv, bmi, anthroheightcalc, anthroweight1lb, anthro_waist_cm, cp_pgs.HAQER_v2, cp_pgs.background_HAQER_v2, cp_pgs.matched_control_HAQER_v2, str_c('pc', 1:5), matches('fsqc'), csection, birth_weight, mom_age_birth, bio_mom, preeclampsia_toxemia, premature_weeks) %>% 
    drop_na(smri_vol_scs_intracranialv, fsqc_qc) %>% 
    group_by(sex) %>%
    mutate(smri_vol_scs_intracranialv2 = qnorm((rank(smri_vol_scs_intracranialv,na.last="keep")-0.5)/sum(!is.na(smri_vol_scs_intracranialv))),
           resid_vol = scale(resid(lm(smri_vol_scs_intracranialv2 ~ mri_age + pc1 + pc2 + pc3 + pc4 + pc5 + anthro_waist_cm + anthroheightcalc + anthroweight1lb + fsqc_nrev + fsqc_revdisp + fsqc_qc + fsqc_qu_motion + fsqc_qu_pialover + fsqc_qu_wmunder + fsqc_qu_inhomogeneity + fsqc_qu_artifact)))[,1],
           obstetric_risk = ifelse(resid_vol > 1.96, 1, 0)) %>% 
    drop_na()

## sample size
n_FALSE <- sum(abcd_df_cpd$obstetric_risk == 0)
n_TRUE <- sum(abcd_df_cpd$obstetric_risk == 1)


## run ES-PGS on cephalo-pelvic disproportion
abcd_cpd_es_pgs_stats <- broom::tidy(glm(obstetric_risk ~ cp_pgs.HAQER_v2 + cp_pgs.background_HAQER_v2 + cp_pgs.matched_control_HAQER_v2, data = abcd_df_cpd)) %>% 
    filter(term != '(Intercept)') %>% 
    rename(x = term, beta = estimate) %>%
    mutate(pheno = 'cephalopelvic_disproportion',
           n = n_FALSE + n_TRUE,
           n_TRUE = n_TRUE,
           n_FALSE = n_FALSE) %>% 
    relocate(pheno)

## make labels for figure
res_obs <- abcd_cpd_es_pgs_stats %>% 
    mutate(lab = str_c('Beta = ', round(beta, digits = 2), ', p-val = ', formatC(p.value, digits = 2)))

p_haq <- abcd_df_cpd %>% 
    mutate(obstetric_risk = case_when(obstetric_risk == 1 ~ 'High risk',
                                      obstetric_risk == 0 ~ 'Typical',
                                      TRUE ~ NA_character_)) %>%
    group_by(obstetric_risk) %>% 
    mutate(n = n()) %>% 
    ungroup() %>% 
    mutate(obstetric_risk = str_c(obstetric_risk, '\nN = ', prettyNum(n, big.mark = ','))) %>% 
    arrange(desc(obstetric_risk)) %>%
    mutate(obstetric_risk = factor(obstetric_risk, levels = unique(obstetric_risk))) %>%
    mutate(lab = res_obs$lab[res_obs$x == 'cp_pgs.HAQER_v2']) %>%
    ggplot(aes(x = obstetric_risk, y = cp_pgs.HAQER_v2)) +
    geom_violin(size = 1.4) +
    geom_boxplot(aes(fill = obstetric_risk), width = .3, alpha = .8, size = 1.4) +
    xlab("Cephalopelvic disproportion") +
    ylab("HAQER CP-PGS") +
    geom_hline(yintercept = 0, color = 'red', linetype = 'dashed', size = 1.075) +
    scale_fill_manual(values = c('grey70', 'chocolate1')) +
    theme_classic(base_size = 18) +
    theme(legend.position = 'none') +
    geom_text(aes(x = 1.5, y = 3.75, label = lab), check_overlap = TRUE, size = 5)

p_bg <- abcd_df_cpd %>% 
    mutate(obstetric_risk = case_when(obstetric_risk == 1 ~ 'High risk',
                                      obstetric_risk == 0 ~ 'Typical',
                                      TRUE ~ NA_character_)) %>%
    group_by(obstetric_risk) %>% 
    mutate(n = n()) %>% 
    ungroup() %>% 
    mutate(obstetric_risk = str_c(obstetric_risk, '\nN = ', prettyNum(n, big.mark = ','))) %>% 
    arrange(desc(obstetric_risk)) %>%
    mutate(obstetric_risk = factor(obstetric_risk, levels = unique(obstetric_risk))) %>%
    ggplot(aes(x = obstetric_risk, y = cp_pgs.background_HAQER_v2)) +
    geom_violin(size = 1.4) +
    geom_boxplot(aes(fill = obstetric_risk), width = .3, alpha = .8, size = 1.4) +
    xlab("Cephalopelvic disproportion") +
    ylab("Background CP-PGS") +
    geom_hline(yintercept = 0, color = 'red', linetype = 'dashed', size = 1.075) +
    scale_fill_manual(values = c('grey70', 'chocolate1')) +
    theme_classic(base_size = 18) +
    theme(legend.position = 'none') +
    geom_text(aes(x = 1.5, y = 3.75, label = res_obs$lab[res_obs$x == 'cp_pgs.background_HAQER_v2']), check_overlap = TRUE, size = 5)

p_con <- abcd_df_cpd %>% 
    mutate(obstetric_risk = case_when(obstetric_risk == 1 ~ 'High risk',
                                      obstetric_risk == 0 ~ 'Typical',
                                      TRUE ~ NA_character_)) %>%
    group_by(obstetric_risk) %>% 
    mutate(n = n()) %>% 
    ungroup() %>% 
    mutate(obstetric_risk = str_c(obstetric_risk, '\nN = ', prettyNum(n, big.mark = ','))) %>% 
    arrange(desc(obstetric_risk)) %>%
    mutate(obstetric_risk = factor(obstetric_risk, levels = unique(obstetric_risk))) %>%
    ggplot(aes(x = obstetric_risk, y = cp_pgs.matched_control_HAQER_v2)) +
    geom_violin(size = 1.4) +
    geom_boxplot(aes(fill = obstetric_risk), width = .3, alpha = .8, size = 1.4) +
    xlab("Cephalopelvic disproportion") +
    ylab("Matched CP-PGS") +
    geom_hline(yintercept = 0, color = 'red', linetype = 'dashed', size = 1.075) +
    scale_fill_manual(values = c('grey70', 'chocolate1')) +
    theme_classic(base_size = 18) +
    theme(legend.position = 'none') +
    geom_text(aes(x = 1.5, y = 3.75, label = res_obs$lab[res_obs$x == 'cp_pgs.matched_control_HAQER_v2']), check_overlap = TRUE, size = 5)


####################################################################
## link cephalopelvic disproportion to c-sections and birth weight
####################################################################
## csection ~ cephalopelvic disproportion + ...
n_FALSE = sum(abcd_df_cpd$csection == 0)
n_TRUE = sum(abcd_df_cpd$csection == 1)

abcd_csection_cpd_stats <- broom::tidy(glm(csection ~ as.factor(obstetric_risk) + pc1 + pc2 + pc3 + pc4 + pc5 + as.factor(sex) + mri_age + premature_weeks + as.factor(preeclampsia_toxemia) + as.factor(bio_mom), data = abcd_df_cpd, family = 'binomial')) %>% 
    filter(term == 'as.factor(obstetric_risk)1') %>% 
    rename(beta = estimate, x = term) %>%
    mutate(x = 'cephalopelvic_disproportion',
           pheno = 'csection',
           n = n_FALSE + n_TRUE,
           n_TRUE = n_TRUE,
           n_FALSE = n_FALSE) %>% 
    relocate(pheno)

## birth weight ~ cephalopelvic disproportion + ...
n <- nrow(abcd_df_cpd)
abcd_bw_cpd_stats <- broom::tidy(lm(birth_weight ~ as.factor(obstetric_risk) + pc1 + pc2 + pc3 + pc4 + pc5 + as.factor(sex) + mri_age + premature_weeks + as.factor(preeclampsia_toxemia) + as.factor(bio_mom) + as.factor(csection), data = abcd_df_cpd)) %>%
    filter(term == 'as.factor(obstetric_risk)1') %>% 
    rename(beta = estimate, x = term) %>%
    mutate(x = 'cephalopelvic_disproportion',
           pheno = 'birth_weight',
           n = n) %>% 
    relocate(pheno)

## save figures
p_haq %>% 
    ggsave(filename = 'manuscript/figures/ABCD_cephalopelvic_disproportion_HAQER_ES-PGS.png',
           device = 'png', dpi = 300, bg = 'white', 
           units = 'in', width = 6, height = 6)

p_all <- p_haq + p_bg + p_con
p_all %>% 
    ggsave(filename = 'manuscript/figures/ABCD_cephalopelvic_disproportion_HAQER_ES-PGS_comparison.png',
           device = 'png', dpi = 300, bg = 'white', 
           units = 'in', width = 18, height = 6)

#######################################
## gather stats
#######################################
## merge stats results and save to file
bind_rows(abcd_cog_espgs_stats, abcd_cpd_es_pgs_stats, abcd_csection_cpd_stats, abcd_bw_cpd_stats) %>% 
    write_csv('manuscript/supplemental_materials/stats/ABCD_HAQER_ES-PGS_results.csv')

#######################################
## save figure objects
#######################################
p_haq %>% 
    write_rds('manuscript/figures/R_plot_objects/ABCD_cephalopelvic_disproportion_HAQER_ES-PGS.rds')
p_all  %>% 
    write_rds('manuscript/figures/R_plot_objects/ABCD_cephalopelvic_disproportion_HAQER_ES-PGS_comparison.rds')  
