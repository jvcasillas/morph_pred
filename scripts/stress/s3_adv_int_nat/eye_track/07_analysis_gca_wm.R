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
         wm_std = (wm - mean(wm, na.rm = T)) / sd(wm, na.rm = T),
         phon_std = (phon_prob - mean(phon_prob)) / sd(phon_prob),
         biphon_std =(biphon_prob - mean(biphon_prob)) / sd(biphon_prob)) %>%
  poly_add_columns(., time_zero, degree = 3, prefix = "ot")

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



# add wm effect to intercept, linear slope, quadratic, and cubic time terms
gca_mod_ss_wm_0 <- update(gca_mod_ss_base,   . ~ . + wm_std) # WM as the only fixed effect
gca_mod_ss_wm_1 <- update(gca_mod_ss_wm_0,   . ~ . + ot1:wm_std)
gca_mod_ss_wm_2 <- update(gca_mod_ss_wm_1,   . ~ . + ot2:wm_std)
gca_mod_ss_wm_3 <- update(gca_mod_ss_wm_2,   . ~ . + ot3:wm_std)


#                 Df   AIC   BIC logLik deviance  Chisq Chi Df Pr(>Chisq)
# gca_mod_ss_wm_0 37 33446 33695 -16686    33372
# gca_mod_ss_wm_1 38 33448 33704 -16686    33372 0.0347      1    0.85226
# gca_mod_ss_wm_2 39 33448 33710 -16685    33370 2.8118      1    0.09357 .
# gca_mod_ss_wm_3 40 33450 33719 -16685    33370 0.0895      1    0.76483

# add wm effect to intercept, linear slope, quadratic, and cubic time terms
gca_mod_ss_wm_0_all <- update(gca_mod_ss_cond_3,   . ~ . + wm_std) # Stress + coda + wm as fixed effects
gca_mod_ss_wm_1_all <- update(gca_mod_ss_wm_0_all,   . ~ . + ot1:wm_std)
gca_mod_ss_wm_2_all <- update(gca_mod_ss_wm_1_all,   . ~ . + ot2:wm_std)
gca_mod_ss_wm_3_all <- update(gca_mod_ss_wm_2_all,   . ~ . + ot3:wm_std)

anova(gca_mod_ss_wm_0_all, gca_mod_ss_wm_1_all, gca_mod_ss_wm_2_all, gca_mod_ss_wm_3_all)

#                     Df   AIC   BIC logLik deviance  Chisq Chi Df Pr(>Chisq)
# gca_mod_ss_wm_0_all 45 33448 33751 -16679    33358
# gca_mod_ss_wm_1_all 46 33450 33760 -16679    33358 0.0405      1    0.84042
# gca_mod_ss_wm_2_all 47 33449 33765 -16678    33355 2.7777      1    0.09558 .
# gca_mod_ss_wm_3_all 48 33451 33774 -16678    33355 0.0738      1    0.78594

# -----------------------------------------------------------------------------

###########################
# Interactions:           #
# 1. Ox -> WM:coda        #
# 2. Parox -> WM:coda     #
# 3. CV -> WM:stress      #
# 4. CVC -> WM:stress     #
###########################

# 1. For Paroxytone targets, check interaction of WM:Coda

stress_gc_parox <- stress_gc_subset %>%
  filter(., condition == "stressed")


gca_mod_ss_parox_base <-
  lmer(eLog ~ 1 + (ot1 + ot2 + ot3) + wm_std + coda_sum +
         (1 + ot1 + ot2 + ot3 | participant) +
         (1 + ot1 + ot2 + ot3 | target),
       control = lmerControl(optimizer = 'bobyqa'), REML = F,
       data = filter(stress_gc_parox, group == "ss"))

gca_mod_ss_parox_int_0 <- update(gca_mod_ss_parox_base, . ~ . + wm_std:coda_sum)
gca_mod_ss_parox_int_1 <- update(gca_mod_ss_parox_int_0, . ~ . + ot1:wm_std:coda_sum)
gca_mod_ss_parox_int_2 <- update(gca_mod_ss_parox_int_1, . ~ . + ot2:wm_std:coda_sum)
gca_mod_ss_parox_int_3 <- update(gca_mod_ss_parox_int_2, . ~ . + ot3:wm_std:coda_sum)

anova(gca_mod_ss_parox_int_0, gca_mod_ss_parox_int_1, gca_mod_ss_parox_int_2, gca_mod_ss_parox_int_3)

#                     Df   AIC   BIC  logLik deviance  Chisq Chi Df Pr(>Chisq)
# gca_mod_parox_int_0 28 16488 16657 -8216.0    16432
# gca_mod_parox_int_1 29 16483 16658 -8212.4    16425 7.1896      1   0.007333 **
# gca_mod_parox_int_2 30 16484 16665 -8212.3    16424 0.2518      1   0.615830
# gca_mod_parox_int_3 31 16485 16672 -8211.6    16423 1.3908      1   0.238273




# 2. For Oxytone targets, check interaction of WM:Coda

stress_gc_ox <- stress_gc_subset %>%
  filter(., condition == "unstressed")

gca_mod_ss_ox_base <-
  lmer(eLog ~ 1 + (ot1 + ot2 + ot3) + wm_std + coda_sum +
         (1 + coda_sum + (ot1 + ot2 + ot3) | participant) +
         (1 + ot1 + ot2 + ot3 | target),
       control = lmerControl(optimizer = 'bobyqa'), REML = F,
       data = filter(stress_gc_ox, group == "ss"))

gca_mod_ss_ox_int_0 <- update(gca_mod_ss_ox_base, . ~ . + wm_std:coda_sum)
gca_mod_ss_ox_int_1 <- update(gca_mod_ss_ox_int_0, . ~ . + ot1:wm_std:coda_sum)
gca_mod_ss_ox_int_2 <- update(gca_mod_ss_ox_int_1, . ~ . + ot2:wm_std:coda_sum)
gca_mod_ss_ox_int_3 <- update(gca_mod_ss_ox_int_2, . ~ . + ot3:wm_std:coda_sum)

anova(gca_mod_ss_ox_int_0, gca_mod_ss_ox_int_1, gca_mod_ss_ox_int_2, gca_mod_ss_ox_int_3)

#                  Df   AIC   BIC  logLik deviance   Chisq Chi Df Pr(>Chisq)
# gca_mod_ox_int_0 33 16701 16901 -8317.6    16635
# gca_mod_ox_int_1 34 16686 16892 -8309.1    16618 17.0229      1  3.693e-05 ***
# gca_mod_ox_int_2 35 16688 16900 -8309.0    16618  0.1863      1     0.6660
# gca_mod_ox_int_3 36 16689 16907 -8308.5    16617  0.9879      1     0.3203



# 3. For CV targets, check interaction of WM:stress

stress_gc_cv <- stress_gc_subset %>%
  filter(., coda == "0")

gca_mod_ss_cv_base <-
  lmer(eLog ~ 1 + (ot1 + ot2 + ot3) + wm_std + condition_sum +
         (1 + ot1 + ot2 + ot3 | participant) +
         (1 + ot1 + ot2 + ot3 | target),
       control = lmerControl(optimizer = 'bobyqa'), REML = F,
       data = filter(stress_gc_cv, group == "ss"))

gca_mod_ss_cv_int_0 <- update(gca_mod_ss_cv_base, . ~ . + wm_std:condition_sum)
gca_mod_ss_cv_int_1 <- update(gca_mod_ss_cv_int_0, . ~ . + ot1:wm_std:condition_sum)
gca_mod_ss_cv_int_2 <- update(gca_mod_ss_cv_int_1, . ~ . + ot2:wm_std:condition_sum)
gca_mod_ss_cv_int_3 <- update(gca_mod_ss_cv_int_2, . ~ . + ot3:wm_std:condition_sum)

anova(gca_mod_ss_cv_int_0, gca_mod_ss_cv_int_1, gca_mod_ss_cv_int_2, gca_mod_ss_cv_int_3)

#                  Df   AIC   BIC  logLik deviance   Chisq Chi Df Pr(>Chisq)
# gca_mod_cv_int_0 28 14138 14302 -7040.9    14082
# gca_mod_cv_int_1 29 14121 14291 -7031.5    14063 18.7382      1    1.5e-05 ***
# gca_mod_cv_int_2 30 14120 14296 -7029.9    14060  3.3016      1    0.06921 .
# gca_mod_cv_int_3 31 14122 14303 -7029.9    14060  0.0109      1    0.91699


# 4. For CVC targets, check interaction of WM:stress

stress_gc_cvc <- stress_gc_subset %>%
  filter(., coda == "1")

gca_mod_ss_cvc_base <-
  lmer(eLog ~ 1 + (ot1 + ot2 + ot3) + wm_std + condition_sum +
         (1 + ot1 + ot2 + ot3 | participant) +
         (1 + ot1 + ot2 + ot3 | target),
       control = lmerControl(optimizer = 'bobyqa'), REML = F,
       data = filter(stress_gc_cvc, group == "ss"))

gca_mod_ss_cvc_int_0 <- update(gca_mod_ss_cvc_base, . ~ . + wm_std:condition_sum)
gca_mod_ss_cvc_int_1 <- update(gca_mod_ss_cvc_int_0, . ~ . + ot1:wm_std:condition_sum)
gca_mod_ss_cvc_int_2 <- update(gca_mod_ss_cvc_int_1, . ~ . + ot2:wm_std:condition_sum)
gca_mod_ss_cvc_int_3 <- update(gca_mod_ss_cvc_int_2, . ~ . + ot3:wm_std:condition_sum)

anova(gca_mod_ss_cvc_int_0, gca_mod_ss_cvc_int_1, gca_mod_ss_cvc_int_2, gca_mod_ss_cvc_int_3)

#                   Df   AIC   BIC  logLik deviance  Chisq Chi Df Pr(>Chisq)
# gca_mod_cvc_int_0 28 19545 19719 -9744.7    19489
# gca_mod_cvc_int_1 29 19542 19721 -9742.0    19484 5.4244      1   0.019857 *
# gca_mod_cvc_int_2 30 19534 19720 -9737.1    19474 9.6971      1   0.001846 **
# gca_mod_cvc_int_3 31 19536 19728 -9737.1    19474 0.0740      1   0.785545



# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------




#
# only la
#

gca_mod_la_base <-
  lmer(eLog ~ 1 + (ot1 + ot2 + ot3) +
         (1 + coda_sum + condition_sum + (ot1 + ot2 + ot3) | participant) +
         (1 + ot1 + ot2 + ot3 | target),
       control = lmerControl(optimizer = 'bobyqa'), REML = F,
       data = filter(stress_gc_subset, group == "la"))


# add wm effect to intercept, linear slope, quadratic, and cubic time terms
gca_mod_la_wm_0 <- update(gca_mod_la_base,   . ~ . + wm_std) # only WM as fixed effect
gca_mod_la_wm_1 <- update(gca_mod_la_wm_0,   . ~ . + ot1:wm_std)
gca_mod_la_wm_2 <- update(gca_mod_la_wm_1,   . ~ . + ot2:wm_std)
gca_mod_la_wm_3 <- update(gca_mod_la_wm_2,   . ~ . + ot3:wm_std)

anova(gca_mod_la_wm_0, gca_mod_la_wm_1, gca_mod_la_wm_2, gca_mod_la_wm_3)

#                 Df   AIC   BIC logLik deviance  Chisq Chi Df Pr(>Chisq)
# gca_mod_la_wm_0 37 34943 35193 -17434    34869
# gca_mod_la_wm_1 38 34945 35201 -17434    34869 0.1543      1     0.6945
# gca_mod_la_wm_2 39 34947 35210 -17434    34869 0.2409      1     0.6236
# gca_mod_la_wm_3 40 34948 35219 -17434    34868 0.0716      1     0.7890

# -----------------------------------------------------------------------------

###########################
# Interactions:           #
# 1. Ox -> WM:coda        #
# 2. Parox -> WM:coda     #
# 3. CV -> WM:stress      #
# 4. CVC -> WM:stress     #
###########################

# 1. For Paroxytone targets, check interaction of WM:Coda


gca_mod_la_parox_base <-
  lmer(eLog ~ 1 + (ot1 + ot2 + ot3) + wm_std + coda_sum +
         (1 + ot1 + ot2 + ot3 | participant) +
         (1 + ot1 + ot2 + ot3 | target),
       control = lmerControl(optimizer = 'bobyqa'), REML = F,
       data = filter(stress_gc_parox, group == "la"))

gca_mod_la_parox_int_0 <- update(gca_mod_la_parox_base, . ~ . + wm_std:coda_sum)
gca_mod_la_parox_int_1 <- update(gca_mod_la_parox_int_0, . ~ . + ot1:wm_std:coda_sum)
gca_mod_la_parox_int_2 <- update(gca_mod_la_parox_int_1, . ~ . + ot2:wm_std:coda_sum)
gca_mod_la_parox_int_3 <- update(gca_mod_la_parox_int_2, . ~ . + ot3:wm_std:coda_sum)

anova(gca_mod_la_parox_int_0, gca_mod_la_parox_int_1, gca_mod_la_parox_int_2, gca_mod_la_parox_int_3)

#                        Df   AIC   BIC  logLik deviance  Chisq Chi Df Pr(>Chisq)
# gca_mod_la_parox_int_0 28 17316 17486 -8630.2    17260
# gca_mod_la_parox_int_1 29 17318 17494 -8630.1    17260 0.0844      1    0.77136
# gca_mod_la_parox_int_2 30 17317 17499 -8628.5    17257 3.1663      1    0.07517 .
# gca_mod_la_parox_int_3 31 17319 17507 -8628.5    17257 0.0448      1    0.83229




# 2. For Oxytone targets, check interaction of WM:Coda

gca_mod_la_ox_base <-
  lmer(eLog ~ 1 + (ot1 + ot2 + ot3) + wm_std + coda_sum +
         (1 + coda_sum + (ot1 + ot2 + ot3) | participant) +
         (1 + ot1 + ot2 + ot3 | target),
       control = lmerControl(optimizer = 'bobyqa'), REML = F,
       data = filter(stress_gc_ox, group == "la"))

gca_mod_la_ox_int_0 <- update(gca_mod_la_ox_base, . ~ . + wm_std:coda_sum)
gca_mod_la_ox_int_1 <- update(gca_mod_la_ox_int_0, . ~ . + ot1:wm_std:coda_sum)
gca_mod_la_ox_int_2 <- update(gca_mod_la_ox_int_1, . ~ . + ot2:wm_std:coda_sum)
gca_mod_la_ox_int_3 <- update(gca_mod_la_ox_int_2, . ~ . + ot3:wm_std:coda_sum)

anova(gca_mod_la_ox_int_0, gca_mod_la_ox_int_1, gca_mod_la_ox_int_2, gca_mod_la_ox_int_3)

#                     Df   AIC   BIC  logLik deviance  Chisq Chi Df Pr(>Chisq)
# gca_mod_la_ox_int_0 33 17488 17688 -8711.1    17422
# gca_mod_la_ox_int_1 34 17490 17696 -8711.0    17422 0.2418      1     0.6229
# gca_mod_la_ox_int_2 35 17492 17704 -8711.0    17422 0.1043      1     0.7467
# gca_mod_la_ox_int_3 36 17493 17712 -8710.7    17421 0.4888      1     0.4844


# 3. For CV targets, check interaction of WM:stress

gca_mod_la_cv_base <-
  lmer(eLog ~ 1 + (ot1 + ot2 + ot3) + wm_std + condition_sum +
         (1 + ot1 + ot2 + ot3 | participant) +
         (1 + ot1 + ot2 + ot3 | target),
       control = lmerControl(optimizer = 'bobyqa'), REML = F,
       data = filter(stress_gc_cv, group == "la"))

gca_mod_la_cv_int_0 <- update(gca_mod_la_cv_base, . ~ . + wm_std:condition_sum)
gca_mod_la_cv_int_1 <- update(gca_mod_la_cv_int_0, . ~ . + ot1:wm_std:condition_sum)
gca_mod_la_cv_int_2 <- update(gca_mod_la_cv_int_1, . ~ . + ot2:wm_std:condition_sum)
gca_mod_la_cv_int_3 <- update(gca_mod_la_cv_int_2, . ~ . + ot3:wm_std:condition_sum)

anova(gca_mod_la_cv_int_0, gca_mod_la_cv_int_1, gca_mod_la_cv_int_2, gca_mod_la_cv_int_3)

#                     Df   AIC   BIC  logLik deviance  Chisq Chi Df Pr(>Chisq)
# gca_mod_la_cv_int_0 28 14541 14705 -7242.6    14485
# gca_mod_la_cv_int_1 29 14542 14712 -7241.9    14484 1.2715      1    0.25948
# gca_mod_la_cv_int_2 30 14544 14720 -7241.8    14484 0.2549      1    0.61367
# gca_mod_la_cv_int_3 31 14542 14724 -7240.1    14480 3.3669      1    0.06652 .


# 4. For CVC targets, check interaction of WM:stress


gca_mod_la_cvc_base <-
  lmer(eLog ~ 1 + (ot1 + ot2 + ot3) + wm_std + condition_sum +
         (1 + ot1 + ot2 + ot3 | participant) +
         (1 + ot1 + ot2 + ot3 | target),
       control = lmerControl(optimizer = 'bobyqa'), REML = F,
       data = filter(stress_gc_cvc, group == "la"))

gca_mod_la_cvc_int_0 <- update(gca_mod_la_cvc_base, . ~ . + wm_std:condition_sum)
gca_mod_la_cvc_int_1 <- update(gca_mod_la_cvc_int_0, . ~ . + ot1:wm_std:condition_sum)
gca_mod_la_cvc_int_2 <- update(gca_mod_la_cvc_int_1, . ~ . + ot2:wm_std:condition_sum)
gca_mod_la_cvc_int_3 <- update(gca_mod_la_cvc_int_2, . ~ . + ot3:wm_std:condition_sum)

anova(gca_mod_la_cvc_int_0, gca_mod_la_cvc_int_1, gca_mod_la_cvc_int_2, gca_mod_la_cvc_int_3)

#                      Df   AIC   BIC logLik deviance  Chisq Chi Df Pr(>Chisq)
# gca_mod_la_cvc_int_0 28 20645 20819 -10294    20589
# gca_mod_la_cvc_int_1 29 20646 20826 -10294    20588 1.0949      1    0.29539
# gca_mod_la_cvc_int_2 30 20645 20832 -10293    20585 2.7327      1    0.09831 .
# gca_mod_la_cvc_int_3 31 20646 20839 -10292    20584 1.2244      1    0.26849



# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------



#
# only int
#

gca_mod_int_base <-
  lmer(eLog ~ 1 + (ot1 + ot2 + ot3) +
         (1 + coda_sum + condition_sum + (ot1 + ot2 + ot3) | participant) +
         (1 + ot1 + ot2 + ot3 | target),
       control = lmerControl(optimizer = 'bobyqa'), REML = F,
       data = filter(stress_gc_subset, group == "int"))



# add wm effect to intercept, linear slope, quadratic, and cubic time terms

gca_mod_int_wm_0 <- update(gca_mod_int_base,   . ~ . + wm_std) # only wm as fixed effect
gca_mod_int_wm_1 <- update(gca_mod_int_wm_0,   . ~ . + ot1:wm_std)
gca_mod_int_wm_2 <- update(gca_mod_int_wm_1,   . ~ . + ot2:wm_std)
gca_mod_int_wm_3 <- update(gca_mod_int_wm_2,   . ~ . + ot3:wm_std)

anova(gca_mod_int_wm_0, gca_mod_int_wm_1, gca_mod_int_wm_2, gca_mod_int_wm_3)

#                   Df   AIC   BIC logLik deviance  Chisq Chi Df Pr(>Chisq)
# gca_mod_int_wm_0 37 32834 33081 -16380    32760
# gca_mod_int_wm_1 38 32835 33089 -16380    32759 0.4050      1     0.5245
# gca_mod_int_wm_2 39 32837 33098 -16380    32759 0.2694      1     0.6037
# gca_mod_int_wm_3 40 32838 33105 -16379    32758 1.3810      1     0.2399


# -----------------------------------------------------------------------------

###########################
# Interactions:           #
# 1. Ox -> WM:coda        #
# 2. Parox -> WM:coda     #
# 3. CV -> WM:stress      #
# 4. CVC -> WM:stress     #
###########################

# 1. For Paroxytone targets, check interaction of WM:Coda


gca_mod_int_parox_base <-
  lmer(eLog ~ 1 + (ot1 + ot2 + ot3) + wm_std + coda_sum +
         (1 + ot1 + ot2 + ot3 | participant) +
         (1 + ot1 + ot2 + ot3 | target),
       control = lmerControl(optimizer = 'bobyqa'), REML = F,
       data = filter(stress_gc_parox, group == "int"))

gca_mod_int_parox_int_0 <- update(gca_mod_int_parox_base, . ~ . + wm_std:coda_sum)
gca_mod_int_parox_int_1 <- update(gca_mod_int_parox_int_0, . ~ . + ot1:wm_std:coda_sum)
gca_mod_int_parox_int_2 <- update(gca_mod_int_parox_int_1, . ~ . + ot2:wm_std:coda_sum)
gca_mod_int_parox_int_3 <- update(gca_mod_int_parox_int_2, . ~ . + ot3:wm_std:coda_sum)

anova(gca_mod_int_parox_int_0, gca_mod_int_parox_int_1, gca_mod_int_parox_int_2, gca_mod_int_parox_int_3)

#                         Df   AIC   BIC  logLik deviance  Chisq Chi Df Pr(>Chisq)
# gca_mod_int_parox_int_0 28 16408 16575 -8175.9    16352
# gca_mod_int_parox_int_1 29 16408 16582 -8175.1    16350 1.6402      1     0.2003
# gca_mod_int_parox_int_2 30 16409 16589 -8174.5    16349 1.1957      1     0.2742
# gca_mod_int_parox_int_3 31 16411 16596 -8174.4    16349 0.2093      1     0.6473


# 2. For Oxytone targets, check interaction of WM:Coda

gca_mod_int_ox_base <-
  lmer(eLog ~ 1 + (ot1 + ot2 + ot3) + wm_std + coda_sum +
         (1 + coda_sum + (ot1 + ot2 + ot3) | participant) +
         (1 + ot1 + ot2 + ot3 | target),
       control = lmerControl(optimizer = 'bobyqa'), REML = F,
       data = filter(stress_gc_ox, group == "int"))

gca_mod_int_ox_int_0 <- update(gca_mod_int_ox_base, . ~ . + wm_std:coda_sum)
gca_mod_int_ox_int_1 <- update(gca_mod_int_ox_int_0, . ~ . + ot1:wm_std:coda_sum)
gca_mod_int_ox_int_2 <- update(gca_mod_int_ox_int_1, . ~ . + ot2:wm_std:coda_sum)
gca_mod_int_ox_int_3 <- update(gca_mod_int_ox_int_2, . ~ . + ot3:wm_std:coda_sum)

anova(gca_mod_int_ox_int_0, gca_mod_int_ox_int_1, gca_mod_int_ox_int_2, gca_mod_int_ox_int_3)

#                     Df   AIC   BIC  logLik deviance  Chisq Chi Df Pr(>Chisq)
# gca_mod_int_ox_int_0 33 16338 16536 -8135.8    16272
# gca_mod_int_ox_int_1 34 16339 16543 -8135.3    16271 0.9950      1   0.318522
# gca_mod_int_ox_int_2 35 16340 16550 -8134.9    16270 0.7917      1   0.373574
# gca_mod_int_ox_int_3 36 16333 16549 -8130.5    16261 8.8345      1   0.002956 **


# 3. For CV targets, check interaction of WM:stress

gca_mod_int_cv_base <-
  lmer(eLog ~ 1 + (ot1 + ot2 + ot3) + wm_std + condition_sum +
         (1 + ot1 + ot2 + ot3 | participant) +
         (1 + ot1 + ot2 + ot3 | target),
       control = lmerControl(optimizer = 'bobyqa'), REML = F,
       data = filter(stress_gc_cv, group == "int"))

gca_mod_int_cv_int_0 <- update(gca_mod_int_cv_base, . ~ . + wm_std:condition_sum)
gca_mod_int_cv_int_1 <- update(gca_mod_int_cv_int_0, . ~ . + ot1:wm_std:condition_sum)
gca_mod_int_cv_int_2 <- update(gca_mod_int_cv_int_1, . ~ . + ot2:wm_std:condition_sum)
gca_mod_int_cv_int_3 <- update(gca_mod_int_cv_int_2, . ~ . + ot3:wm_std:condition_sum)

anova(gca_mod_int_cv_int_0, gca_mod_int_cv_int_1, gca_mod_int_cv_int_2, gca_mod_int_cv_int_3)

#                     Df   AIC   BIC  logLik deviance  Chisq Chi Df Pr(>Chisq)
# gca_mod_int_cv_int_0 28 13286 13448 -6614.8    13230
# gca_mod_int_cv_int_1 29 13276 13444 -6608.9    13218 11.875      1  0.0005688 ***
# gca_mod_int_cv_int_2 30 13278 13451 -6608.8    13218  0.088      1  0.7667430
# gca_mod_int_cv_int_3 31 13268 13448 -6603.1    13206 11.355      1  0.0007523 ***


# 4. For CVC targets, check interaction of WM:stress


gca_mod_int_cvc_base <-
  lmer(eLog ~ 1 + (ot1 + ot2 + ot3) + wm_std + condition_sum +
         (1 + ot1 + ot2 + ot3 | participant) +
         (1 + ot1 + ot2 + ot3 | target),
       control = lmerControl(optimizer = 'bobyqa'), REML = F,
       data = filter(stress_gc_cvc, group == "int"))

gca_mod_int_cvc_int_0 <- update(gca_mod_int_cvc_base, . ~ . + wm_std:condition_sum)
gca_mod_int_cvc_int_1 <- update(gca_mod_int_cvc_int_0, . ~ . + ot1:wm_std:condition_sum)
gca_mod_int_cvc_int_2 <- update(gca_mod_int_cvc_int_1, . ~ . + ot2:wm_std:condition_sum)
gca_mod_int_cvc_int_3 <- update(gca_mod_int_cvc_int_2, . ~ . + ot3:wm_std:condition_sum)

anova(gca_mod_int_cvc_int_0, gca_mod_int_cvc_int_1, gca_mod_int_cvc_int_2, gca_mod_int_cvc_int_3)

#                       Df   AIC   BIC  logLik deviance  Chisq Chi Df Pr(>Chisq)
# gca_mod_int_cvc_int_0 28 19671 19843 -9807.4    19615
# gca_mod_int_cvc_int_1 29 19670 19849 -9805.9    19612 2.8406      1    0.09191 .
# gca_mod_int_cvc_int_2 30 19668 19854 -9804.2    19608 3.3936      1    0.06545 .
# gca_mod_int_cvc_int_3 31 19670 19862 -9804.2    19608 0.0102      1    0.91949



# -----------------------------------------------------------------------------











# Model predictions for plottting ---------------------------------------------

# Create design dataframe for predictions
new_dat_all <- stress_gc_subset %>%
  dplyr::select(group, time_zero, ot1:ot3, coda_sum, condition_sum, wm_std) %>%
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
  mod_spec <- c("_base", "_parox_base", "_parox_int_0", "_parox_int_1", "_parox_int_2", "_parox_int_3",
                "_ox_int_0", "_ox_int_1", "_ox_int_2", "_ox_int_3", "_cv_int_0", "_cv_int_1", "_cv_int_2",
                "_cv_int_3", "_cvc_int_0", "_cvc_int_1", "_cvc_int_2", "_cvc_int_3")

  # Store ind models in list
  ind_mods_wm <- mget(c(paste0(mod_type, "ss", mod_spec),
                     paste0(mod_type, "la", mod_spec),
                     paste0(mod_type, "int", mod_spec)))

  # Add age models from interpreters

  ind_mods_wm$gca_mod_ss_parox_int_0 <- gca_mod_ss_parox_int_0
  ind_mods_wm$gca_mod_ss_parox_int_1 <- gca_mod_ss_parox_int_1
  ind_mods_wm$gca_mod_ss_parox_int_2 <- gca_mod_ss_parox_int_2
  ind_mods_wm$gca_mod_ss_parox_int_3 <- gca_mod_ss_parox_int_3

  ind_mods_wm$gca_mod_ss_ox_int_0 <- gca_mod_ss_ox_int_0
  ind_mods_wm$gca_mod_ss_ox_int_1 <- gca_mod_ss_ox_int_1
  ind_mods_wm$gca_mod_ss_ox_int_2 <- gca_mod_ss_ox_int_2
  ind_mods_wm$gca_mod_ss_ox_int_3 <- gca_mod_ss_ox_int_3

  ind_mods_wm$gca_mod_ss_cv_int_0 <- gca_mod_ss_cv_int_0
  ind_mods_wm$gca_mod_ss_cv_int_1 <- gca_mod_ss_cv_int_1
  ind_mods_wm$gca_mod_ss_cv_int_2 <- gca_mod_ss_cv_int_2
  ind_mods_wm$gca_mod_ss_cv_int_3 <- gca_mod_ss_cv_int_3

  ind_mods_wm$gca_mod_ss_cvc_int_0 <- gca_mod_ss_cvc_int_0
  ind_mods_wm$gca_mod_ss_cvc_int_1 <- gca_mod_ss_cvc_int_1
  ind_mods_wm$gca_mod_ss_cvc_int_2 <- gca_mod_ss_cvc_int_2
  ind_mods_wm$gca_mod_ss_cvc_int_3 <- gca_mod_ss_cvc_int_3

  ind_mods_wm$gca_mod_la_parox_int_0 <- gca_mod_la_parox_int_0
  ind_mods_wm$gca_mod_la_parox_int_1 <- gca_mod_la_parox_int_1
  ind_mods_wm$gca_mod_la_parox_int_2 <- gca_mod_la_parox_int_2
  ind_mods_wm$gca_mod_la_parox_int_3 <- gca_mod_la_parox_int_3

  ind_mods_wm$gca_mod_la_ox_int_0 <- gca_mod_la_ox_int_0
  ind_mods_wm$gca_mod_la_ox_int_1 <- gca_mod_la_ox_int_1
  ind_mods_wm$gca_mod_la_ox_int_2 <- gca_mod_la_ox_int_2
  ind_mods_wm$gca_mod_la_ox_int_3 <- gca_mod_la_ox_int_3

  ind_mods_wm$gca_mod_la_cv_int_0 <- gca_mod_la_cv_int_0
  ind_mods_wm$gca_mod_la_cv_int_1 <- gca_mod_la_cv_int_1
  ind_mods_wm$gca_mod_la_cv_int_2 <- gca_mod_la_cv_int_2
  ind_mods_wm$gca_mod_la_cv_int_3 <- gca_mod_la_cv_int_3

  ind_mods_wm$gca_mod_la_cvc_int_0 <- gca_mod_la_cvc_int_0
  ind_mods_wm$gca_mod_la_cvc_int_1 <- gca_mod_la_cvc_int_1
  ind_mods_wm$gca_mod_la_cvc_int_2 <- gca_mod_la_cvc_int_2
  ind_mods_wm$gca_mod_la_cvc_int_3 <- gca_mod_la_cvc_int_3

  ind_mods_wm$gca_mod_int_parox_int_0 <- gca_mod_int_parox_int_0
  ind_mods_wm$gca_mod_int_parox_int_1 <- gca_mod_int_parox_int_1
  ind_mods_wm$gca_mod_int_parox_int_2 <- gca_mod_int_parox_int_2
  ind_mods_wm$gca_mod_int_parox_int_3 <- gca_mod_int_parox_int_3

  ind_mods_wm$gca_mod_int_ox_int_0 <- gca_mod_int_ox_int_0
  ind_mods_wm$gca_mod_int_ox_int_1 <- gca_mod_int_ox_int_1
  ind_mods_wm$gca_mod_int_ox_int_2 <- gca_mod_int_ox_int_2
  ind_mods_wm$gca_mod_int_ox_int_3 <- gca_mod_int_ox_int_3

  ind_mods_wm$gca_mod_int_cv_int_0 <- gca_mod_int_cv_int_0
  ind_mods_wm$gca_mod_int_cv_int_1 <- gca_mod_int_cv_int_1
  ind_mods_wm$gca_mod_int_cv_int_2 <- gca_mod_int_cv_int_2
  ind_mods_wm$gca_mod_int_cv_int_3 <- gca_mod_int_cv_int_3

  ind_mods_wm$gca_mod_int_cvc_int_0 <- gca_mod_int_cvc_int_0
  ind_mods_wm$gca_mod_int_cvc_int_1 <- gca_mod_int_cvc_int_1
  ind_mods_wm$gca_mod_int_cvc_int_2 <- gca_mod_int_cvc_int_2
  ind_mods_wm$gca_mod_int_cvc_int_3 <- gca_mod_int_cvc_int_3


  save(ind_mods_wm,
       file = here("models", "stress", "s3_adv_int_nat", "eye_track", "gca",
                   "ind_mods_wm.Rdata"))

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

