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





# Load data and models --------------------------------------------------------

# Load data + add other variables to the framework ()
source(here::here("scripts", "02_load_data.R"))

mem_data <- read_csv(here("data", "raw", "dur_stress_demographics.csv"))
mem_data$id <- toupper(mem_data$id)
mem_data <- mem_data %>%
  rename(participant = "id") %>%
  select(., -group)


stress50 <- stress50 %>%
  filter(., group %in% c("la", "int", "ss")) %>%
  left_join(mem_data, by = "participant")

stress50$wm <- as.numeric(stress50$wm)


# Get path to saved models
gca_mods_path  <- here("models", "stress", "s3_adv_int_nat", "eye_track", "gca")

# Load models as lists
load(paste0(gca_mods_path, "/ind_mods.Rdata"))
load(paste0(gca_mods_path, "/full_mods.Rdata"))
load(paste0(gca_mods_path, "/nested_model_comparisons.Rdata"))
load(paste0(gca_mods_path, "/model_preds.Rdata"))

# Store objects in global env
list2env(ind_mods, globalenv())
list2env(full_mods, globalenv())
list2env(nested_model_comparisons, globalenv())
list2env(model_preds, globalenv())

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
         wm_std = (wm - mean(wm, na.rm = T) / sd(wm, na.rm = T))) %>%
  poly_add_columns(., time_zero, degree = 3, prefix = "ot")

# -----------------------------------------------------------------------------






# Random effects structure ----------------------------------------------------

# Build up random effects to test time terms
if(F){

  mod_ot1 <-
    lmer(eLog ~ 1 + ot1 +
           (1 + coda_sum + condition_sum + ot1 | participant),
         control = lmerControl(optimizer = 'bobyqa'),
         data = stress_gc_subset, weights = 1/wts, REML = F)

  mod_ot2 <-
    update(mod_ot1, . ~ . -(1 + coda_sum + condition_sum + ot1 | participant) +
             ot2 + (1 + coda_sum + condition_sum + ot1 + ot2 | participant))

  mod_ot3 <-
    update(mod_ot2, . ~ . -(1 + coda_sum + condition_sum + ot1 + ot2 | participant) +
             ot3 + (1 + coda_sum + condition_sum + ot1 + ot2 + ot3 | participant))

  mod_ot4 <- update(mod_ot3, . ~ . + (1 | target))

  anova(mod_ot1, mod_ot2, mod_ot3, mod_ot4)

}

# -----------------------------------------------------------------------------





# Individual models -----------------------------------------------------------

#
# only ss
#


gca_mod_ss_base <-
  lmer(eLog ~ 1 + (ot1 + ot2 + ot3) +
         (1 + coda_sum + condition_sum + (ot1 + ot2 + ot3) | participant) +
         (1 + ot1 + ot2 + ot3 | target),
       control = lmerControl(optimizer = 'bobyqa'), REML = F,
       data = filter(stress_gc_subset, group == "ss"))

# add coda effect to intercept, linear slope, quadratic, and cubic time terms
gca_mod_ss_coda_0 <- update(gca_mod_ss_base,   . ~ . + coda_sum)
gca_mod_ss_coda_1 <- update(gca_mod_ss_coda_0, . ~ . + ot1:coda_sum)
gca_mod_ss_coda_2 <- update(gca_mod_ss_coda_1, . ~ . + ot2:coda_sum)
gca_mod_ss_coda_3 <- update(gca_mod_ss_coda_2, . ~ . + ot3:coda_sum)

ss_coda_anova <-
  anova(gca_mod_ss_base, gca_mod_ss_coda_0, gca_mod_ss_coda_1,
        gca_mod_ss_coda_2, gca_mod_ss_coda_3)

# add coda effect to intercept, linear slope, quadratic, and cubic time terms
gca_mod_ss_cond_0 <- update(gca_mod_ss_coda_3,   . ~ . + condition_sum)
gca_mod_ss_cond_1 <- update(gca_mod_ss_cond_0,   . ~ . + ot1:condition_sum)
gca_mod_ss_cond_2 <- update(gca_mod_ss_cond_1,   . ~ . + ot2:condition_sum)
gca_mod_ss_cond_3 <- update(gca_mod_ss_cond_2,   . ~ . + ot3:condition_sum)

ss_cond_anova <-
  anova(gca_mod_ss_coda_3, gca_mod_ss_cond_0, gca_mod_ss_cond_1,
        gca_mod_ss_cond_2, gca_mod_ss_cond_3)

# add coda x cond int to intercept, linear slope, quadratic, and cubic terms
gca_mod_ss_int_0 <- update(gca_mod_ss_cond_3, . ~ . + coda_sum:condition_sum)
gca_mod_ss_int_1 <- update(gca_mod_ss_int_0,  . ~ . + ot1:coda_sum:condition_sum)
gca_mod_ss_int_2 <- update(gca_mod_ss_int_1,  . ~ . + ot2:coda_sum:condition_sum)
gca_mod_ss_int_3 <- update(gca_mod_ss_int_2,  . ~ . + ot3:coda_sum:condition_sum)

ss_int_anova <-
  anova(gca_mod_ss_cond_3, gca_mod_ss_int_0, gca_mod_ss_int_1,
        gca_mod_ss_int_2, gca_mod_ss_int_3)

# add wm effect to intercept, linear slope, quadratic, and cubic time terms
gca_mod_ss_wm_0 <- update(gca_mod_ss_base,   . ~ . + wm_std)
gca_mod_ss_wm_1 <- update(gca_mod_ss_wm_0,   . ~ . + ot1:wm_std)
gca_mod_ss_wm_2 <- update(gca_mod_ss_wm_1,   . ~ . + ot2:wm_std)
gca_mod_ss_wm_3 <- update(gca_mod_ss_wm_2,   . ~ . + ot3:wm_std)

anova(gca_mod_ss_wm_0, gca_mod_ss_wm_1, gca_mod_ss_wm_2, gca_mod_ss_wm_3)


#
# only la
#

gca_mod_la_base <-
  lmer(eLog ~ 1 + (ot1 + ot2 + ot3) +
         (1 + coda_sum + condition_sum + (ot1 + ot2 + ot3) | participant) +
         (1 + ot1 + ot2 + ot3 | target),
       control = lmerControl(optimizer = 'bobyqa'), REML = F,
       data = filter(stress_gc_subset, group == "la"))

# add coda effect to intercept, linear slope, quadratic, and cubic time terms
gca_mod_la_coda_0 <- update(gca_mod_la_base,   . ~ . + coda_sum)
gca_mod_la_coda_1 <- update(gca_mod_la_coda_0, . ~ . + ot1:coda_sum)
gca_mod_la_coda_2 <- update(gca_mod_la_coda_1, . ~ . + ot2:coda_sum)
gca_mod_la_coda_3 <- update(gca_mod_la_coda_2, . ~ . + ot3:coda_sum)

la_coda_anova <-
  anova(gca_mod_la_base, gca_mod_la_coda_0, gca_mod_la_coda_1,
        gca_mod_la_coda_2, gca_mod_la_coda_3)

# add coda effect to intercept, linear slope, quadratic, and cubic time terms
gca_mod_la_cond_0 <- update(gca_mod_la_coda_3,   . ~ . + condition_sum)
gca_mod_la_cond_1 <- update(gca_mod_la_cond_0,   . ~ . + ot1:condition_sum)
gca_mod_la_cond_2 <- update(gca_mod_la_cond_1,   . ~ . + ot2:condition_sum)
gca_mod_la_cond_3 <- update(gca_mod_la_cond_2,   . ~ . + ot3:condition_sum)

la_cond_anova <-
  anova(gca_mod_la_coda_3, gca_mod_la_cond_0, gca_mod_la_cond_1,
        gca_mod_la_cond_2, gca_mod_la_cond_3)

# add coda x cond int to intercept, linear slope, quadratic, and cubic terms
gca_mod_la_int_0 <- update(gca_mod_la_cond_3, . ~ . + coda_sum:condition_sum)
gca_mod_la_int_1 <- update(gca_mod_la_int_0,  . ~ . + ot1:coda_sum:condition_sum)
gca_mod_la_int_2 <- update(gca_mod_la_int_1,  . ~ . + ot2:coda_sum:condition_sum)
gca_mod_la_int_3 <- update(gca_mod_la_int_2,  . ~ . + ot3:coda_sum:condition_sum)

la_int_anova <-
  anova(gca_mod_la_cond_3, gca_mod_la_int_0, gca_mod_la_int_1,
        gca_mod_la_int_2, gca_mod_la_int_3)

# add wm effect to intercept, linear slope, quadratic, and cubic time terms
gca_mod_la_wm_0 <- update(gca_mod_la_base,   . ~ . + wm_std)
gca_mod_la_wm_1 <- update(gca_mod_la_wm_0,   . ~ . + ot1:wm_std)
gca_mod_la_wm_2 <- update(gca_mod_la_wm_1,   . ~ . + ot2:wm_std)
gca_mod_la_wm_3 <- update(gca_mod_la_wm_2,   . ~ . + ot3:wm_std) # GAVE A WEIRD WARNING MESSAGE

anova(gca_mod_la_wm_0, gca_mod_la_wm_1, gca_mod_la_wm_2, gca_mod_la_wm_3)



#
# only int
#

gca_mod_int_base <-
  lmer(eLog ~ 1 + (ot1 + ot2 + ot3) +
         (1 + coda_sum + condition_sum + (ot1 + ot2 + ot3) | participant) +
         (1 + ot1 + ot2 + ot3 | target),
       control = lmerControl(optimizer = 'bobyqa'), REML = F,
       data = filter(stress_gc_subset, group == "int"))

# add coda effect to intercept, linear slope, quadratic, and cubic time terms
gca_mod_int_coda_0 <- update(gca_mod_int_base,   . ~ . + coda_sum)
gca_mod_int_coda_1 <- update(gca_mod_int_coda_0, . ~ . + ot1:coda_sum)
gca_mod_int_coda_2 <- update(gca_mod_int_coda_1, . ~ . + ot2:coda_sum)
gca_mod_int_coda_3 <- update(gca_mod_int_coda_2, . ~ . + ot3:coda_sum)

int_coda_anova <-
  anova(gca_mod_int_base, gca_mod_int_coda_0, gca_mod_int_coda_1,
        gca_mod_int_coda_2, gca_mod_int_coda_3)

# add coda effect to intercept, linear slope, quadratic, and cubic time terms
gca_mod_int_cond_0 <- update(gca_mod_int_coda_3,   . ~ . + condition_sum)
gca_mod_int_cond_1 <- update(gca_mod_int_cond_0,   . ~ . + ot1:condition_sum)
gca_mod_int_cond_2 <- update(gca_mod_int_cond_1,   . ~ . + ot2:condition_sum)
gca_mod_int_cond_3 <- update(gca_mod_int_cond_2,   . ~ . + ot3:condition_sum)

int_cond_anova <-
  anova(gca_mod_int_coda_3, gca_mod_int_cond_0, gca_mod_int_cond_1,
        gca_mod_int_cond_2, gca_mod_int_cond_3)

# add coda x cond int to intercept, linear slope, quadratic, and cubic terms
gca_mod_int_int_0 <- update(gca_mod_int_cond_3, . ~ . + coda_sum:condition_sum)
gca_mod_int_int_1 <- update(gca_mod_int_int_0,  . ~ . + ot1:coda_sum:condition_sum)
gca_mod_int_int_2 <- update(gca_mod_int_int_1,  . ~ . + ot2:coda_sum:condition_sum)
gca_mod_int_int_3 <- update(gca_mod_int_int_2,  . ~ . + ot3:coda_sum:condition_sum)

int_int_anova <-
  anova(gca_mod_int_cond_3, gca_mod_int_int_0, gca_mod_int_int_1,
        gca_mod_int_int_2, gca_mod_int_int_3)

# Check age

stress_int_age_subset <- stress_gc_subset %>%
  filter(., group == "int") %>%
  left_join(.,
            read_csv(here("data", "raw", "dur_stress_demographics.csv")) %>%
              filter(group == "in", id != "IN17") %>%
              mutate(., group = fct_recode(group, int = "in"),
                     participant = id), by = "participant") %>%
  mutate(., age_std = (age - mean(age)) / sd(age))


gca_mod_int_age <-
  lmer(eLog ~ 1 + (ot1 + ot2 + ot3) * coda_sum * condition_sum +
         (1 + coda_sum + condition_sum + (ot1 + ot2 + ot3) | participant) +
         (1 + ot1 + ot2 + ot3 | target),
       control = lmerControl(optimizer = 'bobyqa'), REML = F,
       data = stress_int_age_subset)

gca_mod_int_age_0 <- update(gca_mod_int_age, . ~ . + age_std)
gca_mod_int_age_1 <- update(gca_mod_int_age_0, . ~ . + ot1:age_std)
gca_mod_int_age_2 <- update(gca_mod_int_age_1, . ~ . + ot2:age_std)
gca_mod_int_age_3 <- update(gca_mod_int_age_2, . ~ . + ot3:age_std)



int_age_anova <-
  anova(gca_mod_int_age, gca_mod_int_age_0, gca_mod_int_age_1,
        gca_mod_int_age_2,  gca_mod_int_age_3)


# add wm effect to intercept, linear slope, quadratic, and cubic time terms


gca_mod_int_base_wm <-
  lmer(eLog ~ 1 + (ot1 + ot2 + ot3) +
         (1 + ot1 + ot2 + ot3 | participant) +
         (1 + ot1 + ot2 + ot3 | target),
       control = lmerControl(optimizer = 'bobyqa'), REML = F,
       data = filter(stress_gc_subset, group == "int"))


gca_mod_int_wm_0 <- update(gca_mod_int_base_wm,   . ~ . + wm_std)
gca_mod_int_wm_1 <- update(gca_mod_int_wm_0,   . ~ . + ot1:wm_std)
gca_mod_int_wm_2 <- update(gca_mod_int_wm_1,   . ~ . + ot2:wm_std)
gca_mod_int_wm_3 <- update(gca_mod_int_wm_2,   . ~ . + ot3:wm_std)

anova(gca_mod_int_wm_0, gca_mod_int_wm_1, gca_mod_int_wm_2, gca_mod_int_wm_3)
# -----------------------------------------------------------------------------






# Full model ------------------------------------------------------------------

if(F){
  # Base model
  gca_full_mod_base <-
    lmer(eLog ~ 1 + (ot1 + ot2 + ot3) * coda_sum * condition_sum +
           (1 + coda_sum + condition_sum + (ot1 + ot2 + ot3) | participant) +
           (1 + ot1 + ot2 + ot3 | target),
         control = lmerControl(optimizer = 'bobyqa',
                               optCtrl = list(maxfun = 2e4)),
         data = stress_gc_subset, REML = F)

  # add group effect to intercept, linear slope, quadratic, and cubic time terms
  gca_full_mod_group_0 <- update(gca_full_mod_base,    . ~ . + group)
  gca_full_mod_group_1 <- update(gca_full_mod_group_0, . ~ . + ot1:group)
  gca_full_mod_group_2 <- update(gca_full_mod_group_1, . ~ . + ot2:group)
  gca_full_mod_group_3 <- update(gca_full_mod_group_2, . ~ . + ot3:group)

  full_group_anova <-
    anova(gca_full_mod_base, gca_full_mod_group_0, gca_full_mod_group_1,
          gca_full_mod_group_2, gca_full_mod_group_3)

  # add 3-way int to intercept, linear slope, quadratic, and cubic time terms
  gca_full_mod_int_0 <- update(gca_full_mod_group_3, . ~ . + coda_sum:condition_sum:group)
  gca_full_mod_int_1 <- update(gca_full_mod_int_0,   . ~ . + ot1:coda_sum:condition_sum:group)
  gca_full_mod_int_2 <- update(gca_full_mod_int_1,   . ~ . + ot2:coda_sum:condition_sum:group)
  gca_full_mod_int_3 <- update(gca_full_mod_int_2,   . ~ . + ot3:coda_sum:condition_sum:group)

  full_int_anova <-
    anova(gca_full_mod_group_3, gca_full_mod_int_0, gca_full_mod_int_1,
          gca_full_mod_int_2, gca_full_mod_int_3)

  # Relevel for pairwise comparisons
  stress_gc_subset %<>% mutate(., group = fct_relevel(group, "int"))
  gca_full_mod_int_relevel <- update(gca_full_mod_int_3)

# WM and group interaction
  gca_mod_full_wm_0 <- update(gca_full_mod_base,   . ~ . + wm_std)
  gca_mod_full_wm_1 <- update(gca_mod_full_wm_0,   . ~ . + ot1:wm_std)
  gca_mod_full_wm_2 <- update(gca_mod_full_wm_1,   . ~ . + ot2:wm_std)
  gca_mod_full_wm_3 <- update(gca_mod_full_wm_2,   . ~ . + ot3:wm_std)

  anova(gca_mod_full_wm_0, gca_mod_full_wm_1, gca_mod_full_wm_2, gca_mod_full_wm_3)

}

# -----------------------------------------------------------------------------





# Model predictions for plottting ---------------------------------------------

# Create design dataframe for predictions
new_dat_all <- stress_gc_subset %>%
  dplyr::select(group, time_zero, ot1:ot3, coda_sum, condition_sum) %>%
  distinct

# Get model predictions and SE
fits_all <- predictSE(gca_full_mod_int_3, new_dat_all) %>%
  as_tibble %>%
  bind_cols(new_dat_all) %>%
  rename(se = se.fit) %>%
  mutate(ymin = fit - se, ymax = fit + se,
         group = fct_recode(group, M = "ss", NIN = "la", IN = "int"))

# Filter preds at target offset
target_offset_preds <- filter(fits_all, time_zero == 4) %>%
  select(group, coda = coda_sum, cond = condition_sum,
         elog = fit, elog_lb = ymin, elog_ub = ymax) %>%
  mutate(prob = plogis(elog),
         prob_lb = plogis(elog_lb),
         prob_ub = plogis(elog_ub)) %>%
  arrange(group)

# -----------------------------------------------------------------------------









# Save models -----------------------------------------------------------------

if(F) {
  # Build model names programatically
  mod_type <- "gca_mod_"
  mod_spec <- c("_base", "_coda_0", "_coda_1", "_coda_2", "_coda_3", "_cond_0",
                "_cond_1", "_cond_2", "_cond_3", "_int_0", "_int_1", "_int_2",
                "_int_3")

  # Store ind models in list
  ind_mods <- mget(c(paste0(mod_type, "ss", mod_spec),
                     paste0(mod_type, "la", mod_spec),
                     paste0(mod_type, "int", mod_spec)))

  # Add age models from interpreters
  ind_mods$gca_mod_int_age_0 <- gca_mod_int_age_0
  ind_mods$gca_mod_int_age_1 <- gca_mod_int_age_1
  ind_mods$gca_mod_int_age_2 <- gca_mod_int_age_2
  ind_mods$gca_mod_int_age_3 <- gca_mod_int_age_3

  save(ind_mods,
       file = here("models", "stress", "s3_adv_int_nat", "eye_track", "gca",
                   "ind_mods.Rdata"))

  # Store full (ot1, ot2, ot3, group, coda, cond) models in list
  full_mods <- mget(c(
    "gca_full_mod_base", "gca_full_mod_group_0", "gca_full_mod_group_1",
    "gca_full_mod_group_2", "gca_full_mod_group_3", "gca_full_mod_int_0",
    "gca_full_mod_int_1", "gca_full_mod_int_2", "gca_full_mod_int_3",
    "gca_full_mod_int_relevel"))

  save(full_mods,
       file = here("models", "stress", "s3_adv_int_nat", "eye_track", "gca",
                   "full_mods.Rdata"))

  # Save anova model comparisons
  nested_model_comparisons <-
    mget(c("ss_coda_anova", "ss_cond_anova", "ss_int_anova",
           "la_coda_anova", "la_cond_anova", "la_int_anova",
           "int_coda_anova", "int_cond_anova", "int_int_anova",
           "int_age_anova", "full_group_anova", "full_int_anova"))

  save(nested_model_comparisons,
       file = here("models", "stress", "s3_adv_int_nat", "eye_track", "gca",
                   "nested_model_comparisons.Rdata"))

  # Save models predictions
  model_preds <- mget(c("fits_all", "target_offset_preds"))

  save(model_preds,
       file = here("models", "stress", "s3_adv_int_nat", "eye_track", "gca",
                   "model_preds.Rdata"))

}

# -----------------------------------------------------------------------------

