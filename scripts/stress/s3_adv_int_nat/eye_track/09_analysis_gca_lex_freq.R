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

# READ:
# There are not interesting effects here
# This has been tested using the frecuencies from LEXESP

source(here::here("scripts", "02_load_data.R"))

mem_data <- read_csv(here("data", "raw", "dur_stress_demographics.csv"))
mem_data$id <- toupper(mem_data$id)
mem_data <- mem_data %>%
  rename(participant = "id") %>%
  select(., -group)


phon_data <- read_csv(here("data", "raw", "phonotactic_frequency.csv"))
phon_data <- phon_data %>%
  select(., -coda)

freq_data <- read_excel(here("data", "raw", "lex_freq_stress.xlsx"))

stress50 <- stress50 %>%
  filter(., group %in% c("la", "int", "ss")) %>%
  left_join(mem_data, by = "participant") %>%
  left_join(phon_data, by = "target") %>%
  left_join(freq_data, by = "target") %>%
  rename(., lex_prob = "freq_nim")


glimpse(stress50)

# Get path to saved models
gca_mods_path  <- here("models", "stress", "s3_adv_int_nat", "eye_track", "gca")

# Load models as lists
load(paste0(gca_mods_path, "/full_mods.Rdata"))
load(paste0(gca_mods_path, "/ind_mods_lex.Rdata"))
load(paste0(gca_mods_path, "/nested_model_comparisons_lex.Rdata"))
load(paste0(gca_mods_path, "/model_preds_lex.Rdata"))

# Store objects in global env
list2env(full_mods, globalenv())
list2env(ind_mods_lex, globalenv())
list2env(nested_model_comparisons_lex, globalenv())
list2env(model_preds_lex, globalenv())

# -----------------------------------------------------------------------------
stress_gc_subset <- stress50 %>%
  filter(., group %in% c('int', 'la', 'ss'),
         !participant %in% c("LA04",
                             "LA06", "LA09", "LA14", "LA15", "LA19"),
         time_zero >= -4 & time_zero <= 12) %>%
  mutate(., group = fct_relevel(group, "ss", "la", "int"),
         condition_sum = if_else(condition == "stressed", 1, -1),
         coda_sum = if_else(coda == 0, 1, -1),
         lex_std = (lex_prob - mean (lex_prob)) / sd(lex_prob)) %>%
  poly_add_columns(., time_zero, degree = 3, prefix = "ot")

# -----------------------------------------------------------------------------



## Adding lex to full model from BLC <gca_full_mod_int_3>
gca_lex_mod_int_0 <- update(gca_full_mod_int_3, . ~ . + lex_std)
gca_lex_mod_int_1 <- update(gca_lex_mod_int_0, . ~ . + ot1:lex_std)
gca_lex_mod_int_2 <- update(gca_lex_mod_int_1, . ~ . + ot2:lex_std)
gca_lex_mod_int_3 <- update(gca_lex_mod_int_2, . ~ . + ot3:lex_std)

full_lex_anova <-
  anova(gca_lex_mod_int_0, gca_lex_mod_int_1, gca_lex_mod_int_2, gca_lex_mod_int_3)

# Nothing significant here

## Interactions
# lex x Group
gca_lex_group_mod_int_0 <- update(gca_full_mod_int_3, . ~ . + lex_std:group)
gca_lex_group_mod_int_1 <- update(gca_lex_group_mod_int_0, . ~ . + ot1:lex_std:group)
gca_lex_group_mod_int_2 <- update(gca_lex_group_mod_int_1, . ~ . + ot2:lex_std:group)
gca_lex_group_mod_int_3 <- update(gca_lex_group_mod_int_2, . ~ . + ot3:lex_std:group)

full_lex_group_anova <-
  anova(gca_lex_group_mod_int_0, gca_lex_group_mod_int_1, gca_lex_group_mod_int_2, gca_lex_group_mod_int_3)
# Nothing significant here


# lex x Stress

gca_lex_stress_mod_int_0 <- update(gca_full_mod_int_3, . ~ . + lex_std:condition_sum)
gca_lex_stress_mod_int_1 <- update(gca_lex_stress_mod_int_0, . ~ . + ot1:lex_std:condition_sum)
gca_lex_stress_mod_int_2 <- update(gca_lex_stress_mod_int_1, . ~ . + ot2:lex_std:condition_sum)
gca_lex_stress_mod_int_3 <- update(gca_lex_stress_mod_int_2, . ~ . + ot3:lex_std:condition_sum)

full_lex_stress_anova <-
  anova(gca_lex_stress_mod_int_0, gca_lex_stress_mod_int_1, gca_lex_stress_mod_int_2, gca_lex_stress_mod_int_3)
# Nothing significant here



# lex x Coda
gca_lex_coda_mod_int_0 <- update(gca_full_mod_int_3, . ~ . + lex_std:coda_sum)
gca_lex_coda_mod_int_1 <- update(gca_lex_coda_mod_int_0, . ~ . + ot1:lex_std:coda_sum)
gca_lex_coda_mod_int_2 <- update(gca_lex_coda_mod_int_1, . ~ . + ot2:lex_std:coda_sum)
gca_lex_coda_mod_int_3 <- update(gca_lex_coda_mod_int_2, . ~ . + ot3:lex_std:coda_sum)

full_lex_coda_anova <-
  anova(gca_lex_coda_mod_int_0, gca_lex_coda_mod_int_1, gca_lex_coda_mod_int_2, gca_lex_coda_mod_int_3)
# Nothing significant here



# lex x coda x group
gca_lex_coda_grp_mod_int_0 <- update(gca_full_mod_int_3, . ~ . + lex_std:coda_sum:group)
gca_lex_coda_grp_mod_int_1 <- update(gca_lex_coda_grp_mod_int_0, . ~ . + ot1:lex_std:coda_sum:group)
gca_lex_coda_grp_mod_int_2 <- update(gca_lex_coda_grp_mod_int_1, . ~ . + ot2:lex_std:coda_sum:group)
gca_lex_coda_grp_mod_int_3 <- update(gca_lex_coda_grp_mod_int_2, . ~ . + ot3:lex_std:coda_sum:group)

full_lex_coda_grp_anova <-
  anova(gca_lex_coda_grp_mod_int_0, gca_lex_coda_grp_mod_int_1, gca_lex_coda_grp_mod_int_2, gca_lex_coda_grp_mod_int_3)

#                            Df    AIC    BIC logLik deviance   Chisq Chi Df Pr(>Chisq)
# gca_lex_coda_grp_mod_int_0 67 104289 104815 -52078   104155
# gca_lex_coda_grp_mod_int_1 70 104281 104830 -52070   104141 14.4648      3   0.002336 **
# gca_lex_coda_grp_mod_int_2 73 104287 104860 -52070   104141  0.1006      3   0.991769
# gca_lex_coda_grp_mod_int_3 76 104290 104887 -52069   104138  2.7945      3   0.424411


# lex x stress x group
gca_lex_condition_sum_grp_mod_int_0 <- update(gca_lex_coda_grp_mod_int_1, . ~ . + lex_std:condition_sum:group)
gca_lex_condition_sum_grp_mod_int_1 <- update(gca_lex_condition_sum_grp_mod_int_0, . ~ . + ot1:lex_std:condition_sum:group)
gca_lex_condition_sum_grp_mod_int_2 <- update(gca_lex_condition_sum_grp_mod_int_1, . ~ . + ot2:lex_std:condition_sum:group)
gca_lex_condition_sum_grp_mod_int_3 <- update(gca_lex_condition_sum_grp_mod_int_2, . ~ . + ot3:lex_std:condition_sum:group)

full_lex_stress_grp_anova <-
  anova(gca_lex_condition_sum_grp_mod_int_0,
        gca_lex_condition_sum_grp_mod_int_1,
        gca_lex_condition_sum_grp_mod_int_2,
        gca_lex_condition_sum_grp_mod_int_3)
# Nothing significant here





# Model predictions for plottting ---------------------------------------------

# Create design dataframe for predictions
new_dat_lex_minus1 <- stress_gc_subset %>%
  dplyr::select(group, time_zero, ot1:ot3, coda_sum, condition_sum, lex_std) %>%
  distinct(.) %>%
  mutate(lex_std = -1)
new_dat_lex_0 <- stress_gc_subset %>%
  dplyr::select(group, time_zero, ot1:ot3, coda_sum, condition_sum, lex_std) %>%
  distinct(.) %>%
  mutate(lex_std = 0)
new_dat_lex_1 <- stress_gc_subset %>%
  dplyr::select(group, time_zero, ot1:ot3, coda_sum, condition_sum, lex_std) %>%
  distinct(.) %>%
  mutate(lex_std = 1)

# Get model predictions and SE
fits_all_lex <-
  bind_rows(
    predictSE(gca_lex_coda_grp_mod_int_1, new_dat_lex_minus1) %>%
      as_tibble %>%
      bind_cols(new_dat_lex_minus1),
    predictSE(gca_lex_coda_grp_mod_int_1, new_dat_lex_0) %>%
      as_tibble %>%
      bind_cols(new_dat_lex_0),
    predictSE(gca_lex_coda_grp_mod_int_1, new_dat_lex_1) %>%
      as_tibble %>%
      bind_cols(new_dat_lex_1)
  ) %>%
  rename(se = se.fit) %>%
  mutate(ymin = fit - se, ymax = fit + se,
         group = fct_recode(group, M = "ss", NIN = "la", IN = "int"))

# Filter preds at target offset
target_offset_lex_preds <- filter(fits_all_lex, time_zero == 4) %>%
  select(group, coda = coda_sum, cond = condition_sum, lex = lex_std,
         elog = fit, elog_lb = ymin, elog_ub = ymax) %>%
  mutate(prob = plogis(elog),
         prob_lb = plogis(elog_lb),
         prob_ub = plogis(elog_ub)) %>%
  arrange(group)

# -----------------------------------------------------------------------------









# Save models -----------------------------------------------------------------

if(F) {

  # Save full models with WM

  full_mods_lex <- mget(c(
    "gca_lex_mod_int_0", "gca_lex_mod_int_1", "gca_lex_mod_int_2", "gca_lex_mod_int_3",
    "gca_lex_group_mod_int_0", "gca_lex_group_mod_int_1", "gca_lex_group_mod_int_2", "gca_lex_group_mod_int_3",
    "gca_lex_stress_mod_int_0", "gca_lex_stress_mod_int_1", "gca_lex_stress_mod_int_2", "gca_lex_stress_mod_int_3",
    "gca_lex_coda_mod_int_0", "gca_lex_coda_mod_int_1", "gca_lex_coda_mod_int_2", "gca_lex_coda_mod_int_3",
    "gca_lex_coda_grp_mod_int_0", "gca_lex_coda_grp_mod_int_1", "gca_lex_coda_grp_mod_int_2", "gca_lex_coda_grp_mod_int_3",
    "gca_lex_condition_sum_grp_mod_int_0", "gca_lex_condition_sum_grp_mod_int_1", "gca_lex_condition_sum_grp_mod_int_2", "gca_lex_condition_sum_grp_mod_int_3"))

  save(full_mods_lex,
       file = here("models", "stress", "s3_adv_int_nat", "eye_track", "gca",
                   "full_mods_lex.Rdata"))




  # Save anova model comparisons
  nested_model_comparisons_lex <-
    mget(c("full_lex_anova", "full_lex_group_anova", "full_lex_stress_anova", "full_lex_coda_anova",
           "full_lex_coda_grp_anova", "full_lex_stress_grp_anova"))




  save(nested_model_comparisons_lex,
       file = here("models", "stress", "s3_adv_int_nat", "eye_track", "gca",
                   "nested_model_comparisons_lex.Rdata"))





  # Save models predictions
  model_preds_lex <- mget(c("fits_all_lex", "target_offset_lex_preds"))

  save(model_preds_lex,
       file = here("models", "stress", "s3_adv_int_nat", "eye_track", "gca",
                   "model_preds_lex.Rdata"))

}

# -----------------------------------------------------------------------------
