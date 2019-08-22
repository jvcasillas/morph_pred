# Growth curve analyisis ------------------------------------------------------
#
# - Question 1: Are the groups different from each other in when they begin
#   to fixate on the target?
#     - test 3 groups at each level of 'condition'
#     - hypothesis: SS has steeper slope for both conditions
# - Question 2: W/in groups, is the there a difference between
#   oxytone/paroxytone items?
#     - test oxytone vs. paroxytone for each group
#     - hypothesis: steeper slope/earlier break in oxytone condition
#
# -----------------------------------------------------------------------------

source(here::here("scripts", "02_load_data.R"))

mem_data <- read_csv(here("data", "raw", "dur_stress_demographics.csv"))
mem_data$id <- toupper(mem_data$id)
mem_data <- mem_data %>%
  rename(participant = "id") %>%
  select(., -group)


phon_data <- read_csv(here("data", "raw", "phonotactic_frequency.csv"))
phon_data <- phon_data %>%
  select(., -coda)

stress50 <- stress50 %>%
  filter(., group %in% c("la", "int", "ss")) %>%
  left_join(mem_data, by = "participant") %>%
  left_join(phon_data, by = "target")

stress50$wm <- as.numeric(stress50$wm)
glimpse(stress50)

# Get path to saved models
gca_mods_path  <- here("models", "stress", "s3_adv_int_nat", "eye_track", "gca")

# Load models as lists
load(paste0(gca_mods_path, "/full_mods.Rdata"))
load(paste0(gca_mods_path, "/full_mods_lang_learn.Rdata"))
load(paste0(gca_mods_path, "/ind_mods_phon.Rdata"))
load(paste0(gca_mods_path, "/nested_model_comparisons_wm.Rdata"))
load(paste0(gca_mods_path, "/model_preds_phon.Rdata"))

# Store objects in global env
list2env(full_mods, globalenv())
list2env(full_mods_lang_learn, globalenv())
list2env(ind_mods_phon, globalenv())
list2env(nested_model_comparisons_wm, globalenv())
list2env(model_preds_wm, globalenv())

# -----------------------------------------------------------------------------
stress_gc_subset <- stress50 %>%
  filter(., group %in% c('int', 'la', 'ss'),
         !participant %in% c("LA04",
                             "LA06", "LA09", "LA14", "LA15", "LA19"),
         time_zero >= -4 & time_zero <= 12) %>%
  mutate(., group = fct_relevel(group, "ss", "la", "int"),
         condition_sum = if_else(condition == "stressed", 1, -1),
         coda_sum = if_else(coda == 0, 1, -1),
         wm_std = (wm - mean(wm, na.rm = T)) / sd(wm, na.rm = T),
         phon_std = (phon_prob - mean(phon_prob)) / sd(phon_prob),
         biphon_std =(biphon_prob - mean(biphon_prob)) / sd(biphon_prob)) %>%
  poly_add_columns(., time_zero, degree = 3, prefix = "ot")

# -----------------------------------------------------------------------------



## Adding phon to full model from BLC <gca_full_mod_int_3>
gca_phon_mod_int_0 <- update(gca_full_mod_int_3, . ~ . + phon_std)
gca_phon_mod_int_1 <- update(gca_phon_mod_int_0, . ~ . + ot1:phon_std)
gca_phon_mod_int_2 <- update(gca_phon_mod_int_1, . ~ . + ot2:phon_std)
gca_phon_mod_int_3 <- update(gca_phon_mod_int_2, . ~ . + ot3:phon_std)

full_phon_anova <-
  anova(gca_phon_mod_int_0, gca_phon_mod_int_1, gca_phon_mod_int_2, gca_phon_mod_int_3)

# Nothing significant here

## Interactions
# phon x Group
gca_phon_group_mod_int_0 <- update(gca_full_mod_int_3, . ~ . + phon_std:group)
gca_phon_group_mod_int_1 <- update(gca_phon_group_mod_int_0, . ~ . + ot1:phon_std:group)
gca_phon_group_mod_int_2 <- update(gca_phon_group_mod_int_1, . ~ . + ot2:phon_std:group)
gca_phon_group_mod_int_3 <- update(gca_phon_group_mod_int_2, . ~ . + ot3:phon_std:group)

full_phon_group_anova <-
  anova(gca_phon_group_mod_int_0, gca_phon_group_mod_int_1, gca_phon_group_mod_int_2, gca_phon_group_mod_int_3)


#                          Df    AIC    BIC logLik deviance   Chisq Chi Df Pr(>Chisq)
# gca_phon_group_mod_int_0 67 104299 104825 -52083   104165
# gca_phon_group_mod_int_1 70 104303 104853 -52082   104163  1.6792      3    0.64156
# gca_phon_group_mod_int_2 73 104299 104872 -52076   104153 10.5332      3    0.01454 *
# gca_phon_group_mod_int_3 76 104302 104899 -52075   104150  2.5327      3    0.46941

# phon x Stress

gca_phon_stress_mod_int_0 <- update(gca_phon_group_mod_int_2, . ~ . + phon_std:condition_sum)
gca_phon_stress_mod_int_1 <- update(gca_phon_stress_mod_int_0, . ~ . + ot1:phon_std:condition_sum)
gca_phon_stress_mod_int_2 <- update(gca_phon_stress_mod_int_1, . ~ . + ot2:phon_std:condition_sum)
gca_phon_stress_mod_int_3 <- update(gca_phon_stress_mod_int_2, . ~ . + ot3:phon_std:condition_sum)

full_phon_stress_anova <-
  anova(gca_phon_stress_mod_int_0, gca_phon_stress_mod_int_1, gca_phon_stress_mod_int_2, gca_phon_stress_mod_int_3)

#                           Df    AIC    BIC logLik deviance  Chisq Chi Df Pr(>Chisq)
# gca_phon_stress_mod_int_0 74 104300 104881 -52076   104152
# gca_phon_stress_mod_int_1 75 104299 104888 -52074   104149 3.6243      1    0.05694 .
# gca_phon_stress_mod_int_2 76 104301 104897 -52074   104149 0.0866      1    0.76849
# gca_phon_stress_mod_int_3 77 104297 104901 -52071   104143 5.9802      1    0.01447 *


# phon x Coda
gca_phon_coda_mod_int_0 <- update(gca_phon_stress_mod_int_3, . ~ . + phon_std:coda_sum)
gca_phon_coda_mod_int_1 <- update(gca_phon_coda_mod_int_0, . ~ . + ot1:phon_std:coda_sum)
gca_phon_coda_mod_int_2 <- update(gca_phon_coda_mod_int_1, . ~ . + ot2:phon_std:coda_sum)
gca_phon_coda_mod_int_3 <- update(gca_phon_coda_mod_int_2, . ~ . + ot3:phon_std:coda_sum)

full_phon_coda_anova <-
  anova(gca_phon_coda_mod_int_0, gca_phon_coda_mod_int_1, gca_phon_coda_mod_int_2, gca_phon_coda_mod_int_3)
# Nothing significant here



# phon x coda x group
gca_phon_coda_grp_mod_int_0 <- update(gca_phon_stress_mod_int_3, . ~ . + phon_std:coda_sum:group)
gca_phon_coda_grp_mod_int_1 <- update(gca_phon_coda_grp_mod_int_0, . ~ . + ot1:phon_std:coda_sum:group)
gca_phon_coda_grp_mod_int_2 <- update(gca_phon_coda_grp_mod_int_1, . ~ . + ot2:phon_std:coda_sum:group)
gca_phon_coda_grp_mod_int_3 <- update(gca_phon_coda_grp_mod_int_2, . ~ . + ot3:phon_std:coda_sum:group)

full_phon_coda_grp_anova <-
  anova(gca_phon_coda_grp_mod_int_0, gca_phon_coda_grp_mod_int_1, gca_phon_coda_grp_mod_int_2, gca_phon_coda_grp_mod_int_3)
# Nothing significant here

# phon x stress x group
gca_phon_condition_sum_grp_mod_int_0 <- update(gca_phon_stress_mod_int_3, . ~ . + phon_std:condition_sum:group)
gca_phon_condition_sum_grp_mod_int_1 <- update(gca_phon_condition_sum_grp_mod_int_0, . ~ . + ot1:phon_std:condition_sum:group)
gca_phon_condition_sum_grp_mod_int_2 <- update(gca_phon_condition_sum_grp_mod_int_1, . ~ . + ot2:phon_std:condition_sum:group)
gca_phon_condition_sum_grp_mod_int_3 <- update(gca_phon_condition_sum_grp_mod_int_2, . ~ . + ot3:phon_std:condition_sum:group)

full_phon_stress_grp_anova <-
  anova(gca_phon_condition_sum_grp_mod_int_0,
        gca_phon_condition_sum_grp_mod_int_1,
        gca_phon_condition_sum_grp_mod_int_2,
        gca_phon_condition_sum_grp_mod_int_3)

#                                      Df    AIC    BIC logLik deviance  Chisq Chi Df Pr(>Chisq)
# gca_phon_condition_sum_grp_mod_int_0 79 104298 104918 -52070   104140
# gca_phon_condition_sum_grp_mod_int_1 81 104293 104929 -52065   104131 9.4322      2    0.00895 **
# gca_phon_condition_sum_grp_mod_int_2 83 104289 104941 -52062   104123 7.4317      2    0.02433 *
# gca_phon_condition_sum_grp_mod_int_3 85 104293 104960 -52061   104123 0.6314      2    0.72927


# phon x stres x group
# Relevel for pairwise comparisons
stress_gc_subset %<>% mutate(., group = fct_relevel(group, "int"))
gca_full_mod_phon_relevel <- update(gca_phon_condition_sum_grp_mod_int_2)




# Model predictions for plottting ---------------------------------------------

# Create design dataframe for predictions
new_dat_phon_minus1 <- stress_gc_subset %>%
  dplyr::select(group, time_zero, ot1:ot3, coda_sum, condition_sum, phon_std) %>%
  distinct(.) %>%
  mutate(phon_std = -1)
new_dat_phon_0 <- stress_gc_subset %>%
  dplyr::select(group, time_zero, ot1:ot3, coda_sum, condition_sum, phon_std) %>%
  distinct(.) %>%
  mutate(phon_std = 0)
new_dat_phon_1 <- stress_gc_subset %>%
  dplyr::select(group, time_zero, ot1:ot3, coda_sum, condition_sum, phon_std) %>%
  distinct(.) %>%
  mutate(phon_std = 1)

# Get model predictions and SE
fits_all_phon <-
  bind_rows(
    predictSE(gca_phon_condition_sum_grp_mod_int_2, new_dat_phon_minus1) %>%
    as_tibble %>%
    bind_cols(new_dat_phon_minus1),
    predictSE(gca_phon_condition_sum_grp_mod_int_2, new_dat_phon_0) %>%
    as_tibble %>%
    bind_cols(new_dat_phon_0),
    predictSE(gca_phon_condition_sum_grp_mod_int_2, new_dat_phon_1) %>%
    as_tibble %>%
    bind_cols(new_dat_phon_1)
  ) %>%
  rename(se = se.fit) %>%
  mutate(ymin = fit - se, ymax = fit + se,
         group = fct_recode(group, M = "ss", NIN = "la", IN = "int"))

# Filter preds at target offset
target_offset_phon_preds <- filter(fits_all_phon, time_zero == 4) %>%
  select(group, coda = coda_sum, cond = condition_sum, phon = phon_std,
         elog = fit, elog_lb = ymin, elog_ub = ymax) %>%
  mutate(prob = plogis(elog),
         prob_lb = plogis(elog_lb),
         prob_ub = plogis(elog_ub)) %>%
  arrange(group)

# -----------------------------------------------------------------------------









# Save models -----------------------------------------------------------------

if(F) {

  # Save full models with WM

  full_mods_phon <- mget(c(
    "gca_phon_mod_int_0", "gca_phon_mod_int_1", "gca_phon_mod_int_2", "gca_phon_mod_int_3",
    "gca_phon_group_mod_int_0", "gca_phon_group_mod_int_1", "gca_phon_group_mod_int_2", "gca_phon_group_mod_int_3",
    "gca_phon_stress_mod_int_0", "gca_phon_stress_mod_int_1", "gca_phon_stress_mod_int_2", "gca_phon_stress_mod_int_3",
    "gca_phon_coda_mod_int_0", "gca_phon_coda_mod_int_1", "gca_phon_coda_mod_int_2", "gca_phon_coda_mod_int_3",
    "gca_phon_coda_grp_mod_int_0", "gca_phon_coda_grp_mod_int_1", "gca_phon_coda_grp_mod_int_2", "gca_phon_coda_grp_mod_int_3",
    "gca_phon_condition_sum_grp_mod_int_0", "gca_phon_condition_sum_grp_mod_int_1", "gca_phon_condition_sum_grp_mod_int_2", "gca_phon_condition_sum_grp_mod_int_3",
    "gca_full_mod_phon_relevel"))

  save(full_mods_phon,
       file = here("models", "stress", "s3_adv_int_nat", "eye_track", "gca",
                   "full_mods_phon.Rdata"))




   # Save anova model comparisons
  nested_model_comparisons_phon <-
    mget(c("full_phon_anova", "full_phon_group_anova", "full_phon_stress_anova", "full_phon_coda_anova",
      "full_phon_coda_grp_anova", "full_phon_stress_grp_anova"))




  save(nested_model_comparisons_phon,
       file = here("models", "stress", "s3_adv_int_nat", "eye_track", "gca",
                   "nested_model_comparisons_phon.Rdata"))





  # Save models predictions
  model_preds_phon <- mget(c("fits_all_phon", "target_offset_phon_preds"))

  save(model_preds_phon,
       file = here("models", "stress", "s3_adv_int_nat", "eye_track", "gca",
                   "model_preds_phon.Rdata"))

}

# -----------------------------------------------------------------------------
