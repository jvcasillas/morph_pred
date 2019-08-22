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
#
# Consider updating here to reflect the specific questions of this paper
#
# -----------------------------------------------------------------------------





# Load data and models --------------------------------------------------------

# Load data + add other variables to the framework ()
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
load(paste0(gca_mods_path, "/ind_mods_wm.Rdata"))
load(paste0(gca_mods_path, "/nested_model_comparisons_wm.Rdata"))
load(paste0(gca_mods_path, "/model_preds_wm.Rdata"))

# Store objects in global env
list2env(full_mods, globalenv())
list2env(full_mods_lang_learn, globalenv())
list2env(ind_mods_wm, globalenv())
list2env(nested_model_comparisons_wm, globalenv())
list2env(model_preds_wm, globalenv())

# -----------------------------------------------------------------------------







# Data prep -------------------------------------------------------------------

# - subset using time course
#    - We need to reduce the time course to a relevant time window that
#      that includes enough of the trajectory from before and after the
#      target syllable onset
#    - Importantly, we need to make sure that the adjusted time course
#      is centered at 200ms after the offset of the first syllable
#    - This is because the orthogonal polynomials center the time course,
#      thus the parameter estimates on the intercept and the linear slope
#      are calculated for the midpoint (0).
#    - This has an added bonus of assessing group differences at the mid
#      point (200ms after target syllable offset), which will corroborate
#      the results from the GLMMs.
#    - We can select the appropriate time course subset by selecting the
#      target syllable offset, bin 4 (200ms / 50 = 4), and keeping an
#      equal number of bins on each side:
#                     8 7 6 5 4 3 2 1 X 1 2 3 4 5 6 7 8
#                                     ^
#                     center of time course (bin 4)
#
#
# Number of bins:     1  2  3  4 5 6 7 8 9 10 11 12 13 14 15 16 17
# Actual bin number: -4 -3 -2 -1 0 1 2 3 4  5  6  7  8  9 10 11 12

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


## Adding WM to full model from BLC <gca_full_mod_int_3>
gca_wm_mod_int_0 <- update(gca_full_mod_int_3, . ~ . + wm_std)
gca_wm_mod_int_1 <- update(gca_wm_mod_int_0, . ~ . + ot1:wm_std)
gca_wm_mod_int_2 <- update(gca_wm_mod_int_1, . ~ . + ot2:wm_std)
gca_wm_mod_int_3 <- update(gca_wm_mod_int_2, . ~ . + ot3:wm_std)

full_wm_anova <-
  anova(gca_wm_mod_int_0, gca_wm_mod_int_1, gca_wm_mod_int_2, gca_wm_mod_int_3)

# Nothing significant here

## Interactions
# WM x Group
gca_wm_group_mod_int_0 <- update(gca_full_mod_int_3, . ~ . + wm_std:group)
gca_wm_group_mod_int_1 <- update(gca_wm_group_mod_int_0, . ~ . + ot1:wm_std:group)
gca_wm_group_mod_int_2 <- update(gca_wm_group_mod_int_1, . ~ . + ot2:wm_std:group)
gca_wm_group_mod_int_3 <- update(gca_wm_group_mod_int_2, . ~ . + ot3:wm_std:group)

full_wm_group_anova <-
  anova(gca_wm_group_mod_int_0, gca_wm_group_mod_int_1, gca_wm_group_mod_int_2, gca_wm_group_mod_int_3)
# Nothing significant here

# WM x Stress

gca_wm_stress_mod_int_0 <- update(gca_full_mod_int_3, . ~ . + wm_std:condition_sum)
gca_wm_stress_mod_int_1 <- update(gca_wm_stress_mod_int_0, . ~ . + ot1:wm_std:condition_sum)
gca_wm_stress_mod_int_2 <- update(gca_wm_stress_mod_int_1, . ~ . + ot2:wm_std:condition_sum)
gca_wm_stress_mod_int_3 <- update(gca_wm_stress_mod_int_2, . ~ . + ot3:wm_std:condition_sum)

full_wm_stress_anova <-
  anova(gca_wm_stress_mod_int_0, gca_wm_stress_mod_int_1, gca_wm_stress_mod_int_2, gca_wm_stress_mod_int_3)

#                         Df    AIC    BIC logLik deviance  Chisq Chi Df Pr(>Chisq)
# gca_wm_stress_mod_int_0 65 101471 101980 -50671   101341
# gca_wm_stress_mod_int_1 66 101468 101984 -50668   101336 5.4179      1    0.01993 *
# gca_wm_stress_mod_int_2 67 101464 101988 -50665   101330 5.7151      1    0.01682 *
# gca_wm_stress_mod_int_3 68 101466 101998 -50665   101330 0.0840      1    0.77193

# WM x Coda
gca_wm_coda_mod_int_0 <- update(gca_wm_stress_mod_int_2, . ~ . + wm_std:coda_sum)
gca_wm_coda_mod_int_1 <- update(gca_wm_coda_mod_int_0, . ~ . + ot1:wm_std:coda_sum)
gca_wm_coda_mod_int_2 <- update(gca_wm_coda_mod_int_1, . ~ . + ot2:wm_std:coda_sum)
gca_wm_coda_mod_int_3 <- update(gca_wm_coda_mod_int_2, . ~ . + ot3:wm_std:coda_sum)

full_wm_coda_anova <-
  anova(gca_wm_coda_mod_int_0, gca_wm_coda_mod_int_1, gca_wm_coda_mod_int_2, gca_wm_coda_mod_int_3)

#                       Df    AIC    BIC logLik deviance  Chisq Chi Df Pr(>Chisq)
# gca_wm_coda_mod_int_0 68 101466 101998 -50665   101330
# gca_wm_coda_mod_int_1 69 101464 102004 -50663   101326 4.1858      1    0.04076 *
# gca_wm_coda_mod_int_2 70 101466 102013 -50663   101326 0.4654      1    0.49509
# gca_wm_coda_mod_int_3 71 101465 102020 -50661   101323 2.8083      1    0.09378 .



# WM x coda x group
gca_wm_coda_grp_mod_int_0 <- update(gca_wm_coda_mod_int_1, . ~ . + wm_std:coda_sum:group)
gca_wm_coda_grp_mod_int_1 <- update(gca_wm_coda_grp_mod_int_0, . ~ . + ot1:wm_std:coda_sum:group)
gca_wm_coda_grp_mod_int_2 <- update(gca_wm_coda_grp_mod_int_1, . ~ . + ot2:wm_std:coda_sum:group)
gca_wm_coda_grp_mod_int_3 <- update(gca_wm_coda_grp_mod_int_2, . ~ . + ot3:wm_std:coda_sum:group)

full_wm_coda_grp_anova <-
  anova(gca_wm_coda_grp_mod_int_0, gca_wm_coda_grp_mod_int_1, gca_wm_coda_grp_mod_int_2, gca_wm_coda_grp_mod_int_3)
# Nothing significant here

# WM x stress x group
gca_wm_condition_sum_grp_mod_int_0 <- update(gca_wm_coda_mod_int_1, . ~ . + wm_std:condition_sum:group)
gca_wm_condition_sum_grp_mod_int_1 <- update(gca_wm_condition_sum_grp_mod_int_0, . ~ . + ot1:wm_std:condition_sum:group)
gca_wm_condition_sum_grp_mod_int_2 <- update(gca_wm_condition_sum_grp_mod_int_1, . ~ . + ot2:wm_std:condition_sum:group)
gca_wm_condition_sum_grp_mod_int_3 <- update(gca_wm_condition_sum_grp_mod_int_2, . ~ . + ot3:wm_std:condition_sum:group)

full_wm_stress_grp_anova <-
  anova(gca_wm_condition_sum_grp_mod_int_0,
        gca_wm_condition_sum_grp_mod_int_1,
        gca_wm_condition_sum_grp_mod_int_2,
        gca_wm_condition_sum_grp_mod_int_3)

#                                    Df    AIC    BIC logLik deviance  Chisq Chi Df Pr(>Chisq)
# gca_wm_condition_sum_grp_mod_int_0 71 101462 102017 -50660   101320
# gca_wm_condition_sum_grp_mod_int_1 73 101460 102031 -50657   101314 6.0313      2    0.04902 *
# gca_wm_condition_sum_grp_mod_int_2 75 101457 102044 -50653   101307 7.1231      2    0.02840 *
# gca_wm_condition_sum_grp_mod_int_3 78 101456 102066 -50650   101300 6.7480      3    0.08038 .

# WM x stres x group
# Relevel for pairwise comparisons
stress_gc_subset %<>% mutate(., group = fct_relevel(group, "int"))
gca_full_mod_wm_relevel <- update(gca_wm_condition_sum_grp_mod_int_2)




# Model predictions for plottting ---------------------------------------------

# Create design dataframe for predictions
new_dat_wm_minus1 <- stress_gc_subset %>%
  dplyr::select(group, time_zero, ot1:ot3, coda_sum, condition_sum, wm_std) %>%
  distinct(.) %>%
  mutate(wm_std = -1)
new_dat_wm_0 <- stress_gc_subset %>%
  dplyr::select(group, time_zero, ot1:ot3, coda_sum, condition_sum, wm_std) %>%
  distinct(.) %>%
  mutate(wm_std = 0)
new_dat_wm_1 <- stress_gc_subset %>%
  dplyr::select(group, time_zero, ot1:ot3, coda_sum, condition_sum, wm_std) %>%
  distinct(.) %>%
  mutate(wm_std = 1)

# Get model predictions and SE
fits_all_wm <-
  bind_rows(
    predictSE(gca_wm_condition_sum_grp_mod_int_2, new_dat_wm_minus1) %>%
      as_tibble %>%
      bind_cols(new_dat_wm_minus1),
    predictSE(gca_wm_condition_sum_grp_mod_int_2, new_dat_wm_0) %>%
      as_tibble %>%
      bind_cols(new_dat_wm_0),
    predictSE(gca_wm_condition_sum_grp_mod_int_2, new_dat_wm_1) %>%
      as_tibble %>%
      bind_cols(new_dat_wm_1)
  ) %>%
  rename(se = se.fit) %>%
  mutate(ymin = fit - se, ymax = fit + se,
         group = fct_recode(group, M = "ss", NIN = "la", IN = "int"))

# Filter preds at target offset
target_offset_wm_preds <- filter(fits_all_wm, time_zero == 4) %>%
  select(group, coda = coda_sum, cond = condition_sum, wm = wm_std,
         elog = fit, elog_lb = ymin, elog_ub = ymax) %>%
  mutate(prob = plogis(elog),
         prob_lb = plogis(elog_lb),
         prob_ub = plogis(elog_ub)) %>%
  arrange(group)

# -----------------------------------------------------------------------------









# Save models -----------------------------------------------------------------

if(F) {

  # Save full models with WM

  full_mods_lang_learn <- mget(c(
    "gca_wm_mod_int_0", "gca_wm_mod_int_1", "gca_wm_mod_int_2", "gca_wm_mod_int_3",
    "gca_wm_group_mod_int_0", "gca_wm_group_mod_int_1", "gca_wm_group_mod_int_2", "gca_wm_group_mod_int_3",
    "gca_wm_stress_mod_int_0", "gca_wm_stress_mod_int_1", "gca_wm_stress_mod_int_2", "gca_wm_stress_mod_int_3",
    "gca_wm_coda_mod_int_0", "gca_wm_coda_mod_int_1", "gca_wm_coda_mod_int_2", "gca_wm_coda_mod_int_3",
    "gca_wm_coda_grp_mod_int_0", "gca_wm_coda_grp_mod_int_1", "gca_wm_coda_grp_mod_int_2", "gca_wm_coda_grp_mod_int_3",
    "gca_wm_condition_sum_grp_mod_int_0", "gca_wm_condition_sum_grp_mod_int_1", "gca_wm_condition_sum_grp_mod_int_2", "gca_wm_condition_sum_grp_mod_int_3",
    "gca_full_mod_wm_relevel"))

  save(full_mods_lang_learn,
       file = here("models", "stress", "s3_adv_int_nat", "eye_track", "gca",
                   "full_mods_lang_learn.Rdata"))




  # Save anova model comparisons
  nested_model_comparisons_wm <-
    mget(c("full_wm_anova", "full_wm_group_anova", "full_wm_stress_anova", "full_wm_coda_anova",
           "full_wm_coda_grp_anova", "full_wm_stress_grp_anova"))




  save(nested_model_comparisons_wm,
       file = here("models", "stress", "s3_adv_int_nat", "eye_track", "gca",
                   "nested_model_comparisons_wm.Rdata"))





  # Save models predictions
  model_preds_wm <- mget(c("fits_all_wm", "target_offset_wm_preds"))

  save(model_preds_wm,
       file = here("models", "stress", "s3_adv_int_nat", "eye_track", "gca",
                   "model_preds_wm.Rdata"))
