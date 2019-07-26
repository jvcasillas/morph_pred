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
load(paste0(gca_mods_path, "/ind_mods_phon.Rdata"))
load(paste0(gca_mods_path, "/nested_model_comparisons_phon.Rdata"))
load(paste0(gca_mods_path, "/model_preds.Rdata"))

# Store objects in global env
list2env(ind_mods_phon, globalenv())
list2env(nested_model_comparisons_phon, globalenv())
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



# add phonotactic frequency effect to intercept, linear slope, quadratic, and cubic time terms
gca_mod_ss_phon_0 <- update(gca_mod_ss_base,   . ~ . + phon_std) # WM as the only fixed effect
gca_mod_ss_phon_1 <- update(gca_mod_ss_phon_0,   . ~ . + ot1:phon_std)
gca_mod_ss_phon_2 <- update(gca_mod_ss_phon_1,   . ~ . + ot2:phon_std)
gca_mod_ss_phon_3 <- update(gca_mod_ss_phon_2,   . ~ . + ot3:phon_std)

ss_phon_anova <-
anova(gca_mod_ss_phon_0, gca_mod_ss_phon_1, gca_mod_ss_phon_2, gca_mod_ss_phon_3)

#                   Df   AIC   BIC logLik deviance  Chisq Chi Df Pr(>Chisq)
# gca_mod_ss_phon_0 37 34957 35208 -17442    34883
# gca_mod_ss_phon_1 38 34959 35216 -17442    34883 0.3016      1     0.5829
# gca_mod_ss_phon_2 39 34959 35223 -17440    34881 2.2483      1     0.1338
# gca_mod_ss_phon_3 40 34961 35231 -17440    34881 0.2734      1     0.6011

# add phonotactic frequency effect to intercept, linear slope, quadratic, and cubic time terms

gca_mod_ss_cond_3 <-
  lmer(eLog ~ 1 + (ot1 + ot2 + ot3) + condition_sum + coda_sum +
         (1 + coda_sum + condition_sum + (ot1 + ot2 + ot3) | participant) +
         (1 + ot1 + ot2 + ot3 | target),
       control = lmerControl(optimizer = 'bobyqa'), REML = F,
       data = filter(stress_gc_subset, group == "ss"))


gca_mod_ss_phon_0_all <- update(gca_mod_ss_cond_3,   . ~ . + phon_std) # Stress + coda + phon as fixed effects
gca_mod_ss_phon_1_all <- update(gca_mod_ss_phon_0_all,   . ~ . + ot1:phon_std)
gca_mod_ss_phon_2_all <- update(gca_mod_ss_phon_1_all,   . ~ . + ot2:phon_std)
gca_mod_ss_phon_3_all <- update(gca_mod_ss_phon_2_all,   . ~ . + ot3:phon_std)

ss_phon_anova_all <-
anova(gca_mod_ss_phon_0_all, gca_mod_ss_phon_1_all, gca_mod_ss_phon_2_all, gca_mod_ss_phon_3_all)

#                       Df   AIC   BIC logLik deviance  Chisq Chi Df Pr(>Chisq)
# gca_mod_ss_phon_0_all 39 34959 35224 -17441    34881
# gca_mod_ss_phon_1_all 40 34961 35232 -17441    34881 0.3052      1     0.5806
# gca_mod_ss_phon_2_all 41 34961 35239 -17440    34879 2.2295      1     0.1354
# gca_mod_ss_phon_3_all 42 34963 35247 -17439    34879 0.2729      1     0.6014

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


gca_mod_ss_parox_phon_base <-
  lmer(eLog ~ 1 + (ot1 + ot2 + ot3) + phon_std + coda_sum +
         (1 + ot1 + ot2 + ot3 | participant) +
         (1 + ot1 + ot2 + ot3 | target),
       control = lmerControl(optimizer = 'bobyqa'), REML = F,
       data = filter(stress_gc_parox, group == "ss"))

gca_mod_ss_parox_phon_int_0 <- update(gca_mod_ss_parox_phon_base, . ~ . + phon_std:coda_sum)
gca_mod_ss_parox_phon_int_1 <- update(gca_mod_ss_parox_phon_int_0, . ~ . + ot1:phon_std:coda_sum)
gca_mod_ss_parox_phon_int_2 <- update(gca_mod_ss_parox_phon_int_1, . ~ . + ot2:phon_std:coda_sum)
gca_mod_ss_parox_phon_int_3 <- update(gca_mod_ss_parox_phon_int_2, . ~ . + ot3:phon_std:coda_sum)

ss_parox_phon_anova <-
anova(gca_mod_ss_parox_phon_int_0, gca_mod_ss_parox_phon_int_1, gca_mod_ss_parox_phon_int_2, gca_mod_ss_parox_phon_int_3)

#                          Df   AIC   BIC  logLik deviance  Chisq Chi Df Pr(>Chisq)
# gca_mod_parox_phon_int_0 28 17313 17483 -8628.7    17257
# gca_mod_parox_phon_int_1 29 17313 17489 -8627.7    17255 1.9022      1     0.1678
# gca_mod_parox_phon_int_2 30 17314 17496 -8627.3    17254 0.9288      1     0.3352
# gca_mod_parox_phon_int_3 31 17316 17504 -8627.0    17254 0.5766      1     0.4477




# 2. For Oxytone targets, check interaction of WM:Coda

stress_gc_ox <- stress_gc_subset %>%
  filter(., condition == "unstressed")

gca_mod_ss_ox_phon_base <-
  lmer(eLog ~ 1 + (ot1 + ot2 + ot3) + phon_std + coda_sum +
         (1 + coda_sum + (ot1 + ot2 + ot3) | participant) +
         (1 + ot1 + ot2 + ot3 | target),
       control = lmerControl(optimizer = 'bobyqa'), REML = F,
       data = filter(stress_gc_ox, group == "ss"))

gca_mod_ss_ox_phon_int_0 <- update(gca_mod_ss_ox_phon_base, . ~ . + phon_std:coda_sum)
gca_mod_ss_ox_phon_int_1 <- update(gca_mod_ss_ox_phon_int_0, . ~ . + ot1:phon_std:coda_sum)
gca_mod_ss_ox_phon_int_2 <- update(gca_mod_ss_ox_phon_int_1, . ~ . + ot2:phon_std:coda_sum)
gca_mod_ss_ox_phon_int_3 <- update(gca_mod_ss_ox_phon_int_2, . ~ . + ot3:phon_std:coda_sum)

ss_ox_phon_anova <-
anova(gca_mod_ss_ox_phon_int_0, gca_mod_ss_ox_phon_int_1, gca_mod_ss_ox_phon_int_2, gca_mod_ss_ox_phon_int_3)

#                       Df   AIC   BIC  logLik deviance  Chisq Chi Df Pr(>Chisq)
# gca_mod_ox_phon_int_0 33 17431 17632 -8682.3    17365
# gca_mod_ox_phon_int_1 34 17432 17639 -8682.2    17364 0.2535      1     0.6146
# gca_mod_ox_phon_int_2 35 17434 17647 -8682.1    17364 0.2611      1     0.6094
# gca_mod_ox_phon_int_3 36 17436 17655 -8682.0    17364 0.0104      1     0.9186



# 3. For CV targets, check interaction of WM:stress

stress_gc_cv <- stress_gc_subset %>%
  filter(., coda == "0")

gca_mod_ss_cv_phon_base <-
  lmer(eLog ~ 1 + (ot1 + ot2 + ot3) + phon_std + condition_sum +
         (1 + ot1 + ot2 + ot3 | participant) +
         (1 + ot1 + ot2 + ot3 | target),
       control = lmerControl(optimizer = 'bobyqa'), REML = F,
       data = filter(stress_gc_cv, group == "ss"))

gca_mod_ss_cv_phon_int_0 <- update(gca_mod_ss_cv_phon_base, . ~ . + phon_std:condition_sum)
gca_mod_ss_cv_phon_int_1 <- update(gca_mod_ss_cv_phon_int_0, . ~ . + ot1:phon_std:condition_sum)
gca_mod_ss_cv_phon_int_2 <- update(gca_mod_ss_cv_phon_int_1, . ~ . + ot2:phon_std:condition_sum)
gca_mod_ss_cv_phon_int_3 <- update(gca_mod_ss_cv_phon_int_2, . ~ . + ot3:phon_std:condition_sum)

ss_cv_phon_anova <-
anova(gca_mod_ss_cv_phon_int_0, gca_mod_ss_cv_phon_int_1, gca_mod_ss_cv_phon_int_2, gca_mod_ss_cv_phon_int_3)

#                       Df   AIC   BIC  logLik deviance  Chisq Chi Df Pr(>Chisq)
# gca_mod_cv_phon_int_0 28 14718 14883 -7331.0    14662
# gca_mod_cv_phon_int_1 29 14720 14891 -7331.0    14662 0.0003      1     0.9867
# gca_mod_cv_phon_int_2 30 14721 14898 -7330.7    14661 0.5535      1     0.4569
# gca_mod_cv_phon_int_3 31 14723 14906 -7330.7    14661 0.0091      1     0.9241


# 4. For CVC targets, check interaction of WM:stress

stress_gc_cvc <- stress_gc_subset %>%
  filter(., coda == "1")

gca_mod_ss_cvc_phon_base <-
  lmer(eLog ~ 1 + (ot1 + ot2 + ot3) + phon_std + condition_sum +
         (1 + ot1 + ot2 + ot3 | participant) +
         (1 + ot1 + ot2 + ot3 | target),
       control = lmerControl(optimizer = 'bobyqa'), REML = F,
       data = filter(stress_gc_cvc, group == "ss"))

gca_mod_ss_cvc_phon_int_0 <- update(gca_mod_ss_cvc_phon_base, . ~ . + phon_std:condition_sum)
gca_mod_ss_cvc_phon_int_1 <- update(gca_mod_ss_cvc_phon_int_0, . ~ . + ot1:phon_std:condition_sum)
gca_mod_ss_cvc_phon_int_2 <- update(gca_mod_ss_cvc_phon_int_1, . ~ . + ot2:phon_std:condition_sum)
gca_mod_ss_cvc_phon_int_3 <- update(gca_mod_ss_cvc_phon_int_2, . ~ . + ot3:phon_std:condition_sum)

ss_cvc_phon_anova <-
anova(gca_mod_ss_cvc_phon_int_0, gca_mod_ss_cvc_phon_int_1, gca_mod_ss_cvc_phon_int_2, gca_mod_ss_cvc_phon_int_3)

#                        Df   AIC   BIC logLik deviance  Chisq Chi Df Pr(>Chisq)
# gca_mod_cvc_phon_int_0 28 20563 20737 -10253    20507
# gca_mod_cvc_phon_int_1 29 20564 20745 -10253    20506 0.8293      1     0.3625
# gca_mod_cvc_phon_int_2 30 20566 20753 -10253    20506 0.0301      1     0.8623
# gca_mod_cvc_phon_int_3 31 20568 20761 -10253    20506 0.0527      1     0.8185



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


# add phon effect to intercept, linear slope, quadratic, and cubic time terms
gca_mod_la_phon_0 <- update(gca_mod_la_base,   . ~ . + phon_std) # only WM as fixed effect
gca_mod_la_phon_1 <- update(gca_mod_la_phon_0,   . ~ . + ot1:phon_std)
gca_mod_la_phon_2 <- update(gca_mod_la_phon_1,   . ~ . + ot2:phon_std)
gca_mod_la_phon_3 <- update(gca_mod_la_phon_2,   . ~ . + ot3:phon_std)


la_phon_anova <-
anova(gca_mod_la_phon_0, gca_mod_la_phon_1, gca_mod_la_phon_2, gca_mod_la_phon_3)

#                   Df   AIC   BIC logLik deviance  Chisq Chi Df Pr(>Chisq)
# gca_mod_la_phon_0 37 36297 36548 -18111    36223
# gca_mod_la_phon_1 38 36299 36557 -18111    36223 0.0226      1     0.8804
# gca_mod_la_phon_2 39 36301 36566 -18111    36223 0.0746      1     0.7848
# gca_mod_la_phon_3 40 36302 36574 -18111    36222 0.5824      1     0.4454


# -----------------------------------------------------------------------------

###########################
# Interactions:           #
# 1. Ox -> WM:coda        #
# 2. Parox -> WM:coda     #
# 3. CV -> WM:stress      #
# 4. CVC -> WM:stress     #
###########################

# 1. For Paroxytone targets, check interaction of WM:Coda


gca_mod_la_parox_phon_base <-
  lmer(eLog ~ 1 + (ot1 + ot2 + ot3) + phon_std + coda_sum +
         (1 + ot1 + ot2 + ot3 | participant) +
         (1 + ot1 + ot2 + ot3 | target),
       control = lmerControl(optimizer = 'bobyqa'), REML = F,
       data = filter(stress_gc_parox, group == "la"))

gca_mod_la_parox_phon_int_0 <- update(gca_mod_la_parox_phon_base, . ~ . + phon_std:coda_sum)
gca_mod_la_parox_phon_int_1 <- update(gca_mod_la_parox_phon_int_0, . ~ . + ot1:phon_std:coda_sum)
gca_mod_la_parox_phon_int_2 <- update(gca_mod_la_parox_phon_int_1, . ~ . + ot2:phon_std:coda_sum)
gca_mod_la_parox_phon_int_3 <- update(gca_mod_la_parox_phon_int_2, . ~ . + ot3:phon_std:coda_sum)

la_parox_phon_anova <-
anova(gca_mod_la_parox_phon_int_0, gca_mod_la_parox_phon_int_1, gca_mod_la_parox_phon_int_2, gca_mod_la_parox_phon_int_3)

#                        Df   AIC   BIC  logLik deviance  Chisq Chi Df Pr(>Chisq)
# gca_mod_la_parox_int_0 28 17921 18091 -8932.4    17865
# gca_mod_la_parox_int_1 29 17922 18098 -8931.8    17864 1.2003      1    0.27327
# gca_mod_la_parox_int_2 30 17924 18106 -8931.8    17864 0.0663      1    0.79686
# gca_mod_la_parox_int_3 31 17921 18110 -8929.6    17859 4.3906      1    0.03614 *



# 2. For Oxytone targets, check interaction of WM:Coda

gca_mod_la_ox_phon_base <-
  lmer(eLog ~ 1 + (ot1 + ot2 + ot3) + phon_std + coda_sum +
         (1 + coda_sum + (ot1 + ot2 + ot3) | participant) +
         (1 + ot1 + ot2 + ot3 | target),
       control = lmerControl(optimizer = 'bobyqa'), REML = F,
       data = filter(stress_gc_ox, group == "la"))

gca_mod_la_ox_phon_int_0 <- update(gca_mod_la_ox_phon_base, . ~ . + phon_std:coda_sum)
gca_mod_la_ox_phon_int_1 <- update(gca_mod_la_ox_phon_int_0, . ~ . + ot1:phon_std:coda_sum)
gca_mod_la_ox_phon_int_2 <- update(gca_mod_la_ox_phon_int_1, . ~ . + ot2:phon_std:coda_sum)
gca_mod_la_ox_phon_int_3 <- update(gca_mod_la_ox_phon_int_2, . ~ . + ot3:phon_std:coda_sum)

la_ox_phon_anova <-
anova(gca_mod_la_ox_phon_int_0, gca_mod_la_ox_phon_int_1, gca_mod_la_ox_phon_int_2, gca_mod_la_ox_phon_int_3)

#                           Df   AIC   BIC  logLik deviance  Chisq Chi Df Pr(>Chisq)
# gca_mod_la_ox_phon_int_0 33 18224 18425 -9078.8    18158
# gca_mod_la_ox_phon_int_1 34 18226 18433 -9078.8    18158 0.0117      1     0.9137
# gca_mod_la_ox_phon_int_2 35 18225 18439 -9077.5    18155 2.5856      1     0.1078
# gca_mod_la_ox_phon_int_3 36 18227 18447 -9077.4    18155 0.0700      1     0.7914


# 3. For CV targets, check interaction of WM:stress

gca_mod_la_cv_phon_base <-
  lmer(eLog ~ 1 + (ot1 + ot2 + ot3) + phon_std + condition_sum +
         (1 + ot1 + ot2 + ot3 | participant) +
         (1 + ot1 + ot2 + ot3 | target),
       control = lmerControl(optimizer = 'bobyqa'), REML = F,
       data = filter(stress_gc_cv, group == "la"))

gca_mod_la_cv_phon_int_0 <- update(gca_mod_la_cv_phon_base, . ~ . + phon_std:condition_sum)
gca_mod_la_cv_phon_int_1 <- update(gca_mod_la_cv_phon_int_0, . ~ . + ot1:phon_std:condition_sum)
gca_mod_la_cv_phon_int_2 <- update(gca_mod_la_cv_phon_int_1, . ~ . + ot2:phon_std:condition_sum)
gca_mod_la_cv_phon_int_3 <- update(gca_mod_la_cv_phon_int_2, . ~ . + ot3:phon_std:condition_sum)

la_cv_phon_anova <-
anova(gca_mod_la_cv_phon_int_0, gca_mod_la_cv_phon_int_1, gca_mod_la_cv_phon_int_2, gca_mod_la_cv_phon_int_3)

#                          Df   AIC   BIC  logLik deviance  Chisq Chi Df Pr(>Chisq)
# gca_mod_la_cv_phon_int_0 28 15122 15288 -7533.2    15066
# gca_mod_la_cv_phon_int_1 29 15124 15296 -7533.2    15066 0.0239      1     0.8772
# gca_mod_la_cv_phon_int_2 30 15126 15303 -7533.2    15066 0.0135      1     0.9074
# gca_mod_la_cv_phon_int_3 31 15128 15311 -7533.0    15066 0.3551      1     0.5512



# 4. For CVC targets, check interaction of WM:stress


gca_mod_la_cvc_phon_base <-
  lmer(eLog ~ 1 + (ot1 + ot2 + ot3) + phon_std + condition_sum +
         (1 + ot1 + ot2 + ot3 | participant) +
         (1 + ot1 + ot2 + ot3 | target),
       control = lmerControl(optimizer = 'bobyqa'), REML = F,
       data = filter(stress_gc_cvc, group == "la"))

gca_mod_la_cvc_phon_int_0 <- update(gca_mod_la_cvc_phon_base, . ~ . + phon_std:condition_sum)
gca_mod_la_cvc_phon_int_1 <- update(gca_mod_la_cvc_phon_int_0, . ~ . + ot1:phon_std:condition_sum)
gca_mod_la_cvc_phon_int_2 <- update(gca_mod_la_cvc_phon_int_1, . ~ . + ot2:phon_std:condition_sum)
gca_mod_la_cvc_phon_int_3 <- update(gca_mod_la_cvc_phon_int_2, . ~ . + ot3:phon_std:condition_sum)

la_cvc_phon_anova <-
anova(gca_mod_la_cvc_phon_int_0, gca_mod_la_cvc_phon_int_1, gca_mod_la_cvc_phon_int_2, gca_mod_la_cvc_phon_int_3)

#                           Df   AIC   BIC logLik deviance  Chisq Chi Df Pr(>Chisq)
# gca_mod_la_cvc_phon_int_0 28 21465 21641 -10705    21409
# gca_mod_la_cvc_phon_int_1 29 21467 21649 -10704    21409 0.2659      1    0.60609
# gca_mod_la_cvc_phon_int_2 30 21467 21655 -10704    21407 1.9438      1    0.16326
# gca_mod_la_cvc_phon_int_3 31 21465 21659 -10701    21403 4.3817      1    0.03633 *


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



# add phon effect to intercept, linear slope, quadratic, and cubic time terms

gca_mod_int_phon_0 <- update(gca_mod_int_base,   . ~ . + phon_std) # only phon as fixed effect
# This model gives the following message
# Warning message:
#   In optwrap(optimizer, devfun, getStart(start, rho$lower, rho$pp),  :
#                convergence code 1 from bobyqa: bobyqa -- maximum number of function evaluations exceeded
gca_mod_int_phon_1 <- update(gca_mod_int_phon_0,   . ~ . + ot1:phon_std)
gca_mod_int_phon_2 <- update(gca_mod_int_phon_1,   . ~ . + ot2:phon_std)
gca_mod_int_phon_3 <- update(gca_mod_int_phon_2,   . ~ . + ot3:phon_std)

int_phon_anova <-
anova(gca_mod_int_phon_0, gca_mod_int_phon_1, gca_mod_int_phon_2, gca_mod_int_phon_3)

#                    Df   AIC   BIC logLik deviance  Chisq Chi Df Pr(>Chisq)
# gca_mod_int_phon_0 37 32835 33082 -16380    32761
# gca_mod_int_phon_1 38 32837 33091 -16380    32761 0.0375      1     0.8464
# gca_mod_int_phon_2 39 32838 33099 -16380    32760 0.6448      1     0.4220
# gca_mod_int_phon_3 40 32840 33107 -16380    32760 0.4125      1     0.5207


# -----------------------------------------------------------------------------

###########################
# Interactions:           #
# 1. Ox -> WM:coda        #
# 2. Parox -> WM:coda     #
# 3. CV -> WM:stress      #
# 4. CVC -> WM:stress     #
###########################

# 1. For Paroxytone targets, check interaction of WM:Coda


gca_mod_int_parox_phon_base <-
  lmer(eLog ~ 1 + (ot1 + ot2 + ot3) + phon_std + coda_sum +
         (1 + ot1 + ot2 + ot3 | participant) +
         (1 + ot1 + ot2 + ot3 | target),
       control = lmerControl(optimizer = 'bobyqa'), REML = F,
       data = filter(stress_gc_parox, group == "int"))

gca_mod_int_parox_phon_int_0 <- update(gca_mod_int_parox_phon_base, . ~ . + phon_std:coda_sum)
gca_mod_int_parox_phon_int_1 <- update(gca_mod_int_parox_phon_int_0, . ~ . + ot1:phon_std:coda_sum)
gca_mod_int_parox_phon_int_2 <- update(gca_mod_int_parox_phon_int_1, . ~ . + ot2:phon_std:coda_sum)
gca_mod_int_parox_phon_int_3 <- update(gca_mod_int_parox_phon_int_2, . ~ . + ot3:phon_std:coda_sum)

int_parox_phon_anova <-
anova(gca_mod_int_parox_phon_int_0, gca_mod_int_parox_phon_int_1, gca_mod_int_parox_phon_int_2, gca_mod_int_parox_phon_int_3)

#                              Df   AIC   BIC  logLik deviance  Chisq Chi Df Pr(>Chisq)
# gca_mod_int_parox_phon_int_0 28 16406 16573 -8174.8    16350
# gca_mod_int_parox_phon_int_1 29 16402 16575 -8171.8    16344 6.0244      1    0.01411 *
# gca_mod_int_parox_phon_int_2 30 16404 16583 -8171.8    16344 0.0743      1    0.78524
# gca_mod_int_parox_phon_int_3 31 16404 16590 -8171.0    16342 1.5594      1    0.21175

# 2. For Oxytone targets, check interaction of WM:Coda

gca_mod_int_ox_phon_base <-
  lmer(eLog ~ 1 + (ot1 + ot2 + ot3) + phon_std + coda_sum +
         (1 + coda_sum + (ot1 + ot2 + ot3) | participant) +
         (1 + ot1 + ot2 + ot3 | target),
       control = lmerControl(optimizer = 'bobyqa'), REML = F,
       data = filter(stress_gc_ox, group == "int"))

gca_mod_int_ox_phon_int_0 <- update(gca_mod_int_ox_phon_base, . ~ . + phon_std:coda_sum)
gca_mod_int_ox_phon_int_1 <- update(gca_mod_int_ox_phon_int_0, . ~ . + ot1:phon_std:coda_sum)
gca_mod_int_ox_phon_int_2 <- update(gca_mod_int_ox_phon_int_1, . ~ . + ot2:phon_std:coda_sum)
gca_mod_int_ox_phon_int_3 <- update(gca_mod_int_ox_phon_int_2, . ~ . + ot3:phon_std:coda_sum)

int_ox_phon_anova <-
anova(gca_mod_int_ox_phon_int_0, gca_mod_int_ox_phon_int_1, gca_mod_int_ox_phon_int_2, gca_mod_int_ox_phon_int_3)

#                           Df   AIC   BIC  logLik deviance  Chisq Chi Df Pr(>Chisq)
# gca_mod_int_ox_phon_int_0 33 16336 16534 -8135.1    16270
# gca_mod_int_ox_phon_int_1 34 16338 16542 -8135.1    16270 0.0001      1     0.9913
# gca_mod_int_ox_phon_int_2 35 16340 16550 -8135.1    16270 0.0010      1     0.9749
# gca_mod_int_ox_phon_int_3 36 16342 16558 -8134.9    16270 0.2937      1     0.5878


# 3. For CV targets, check interaction of WM:stress

gca_mod_int_cv_phon_base <-
  lmer(eLog ~ 1 + (ot1 + ot2 + ot3) + phon_std + condition_sum +
         (1 + ot1 + ot2 + ot3 | participant) +
         (1 + ot1 + ot2 + ot3 | target),
       control = lmerControl(optimizer = 'bobyqa'), REML = F,
       data = filter(stress_gc_cv, group == "int"))

gca_mod_int_cv_phon_int_0 <- update(gca_mod_int_cv_phon_base, . ~ . + phon_std:condition_sum)
gca_mod_int_cv_phon_int_1 <- update(gca_mod_int_cv_phon_int_0, . ~ . + ot1:phon_std:condition_sum)
gca_mod_int_cv_phon_int_2 <- update(gca_mod_int_cv_phon_int_1, . ~ . + ot2:phon_std:condition_sum)
gca_mod_int_cv_phon_int_3 <- update(gca_mod_int_cv_phon_int_2, . ~ . + ot3:phon_std:condition_sum)

int_cv_phon_anova <-
anova(gca_mod_int_cv_phon_int_0, gca_mod_int_cv_phon_int_1, gca_mod_int_cv_phon_int_2, gca_mod_int_cv_phon_int_3)

#                           Df   AIC   BIC  logLik deviance  Chisq Chi Df Pr(>Chisq)
# gca_mod_int_cv_phon_int_0 28 13305 13467 -6624.6    13249
# gca_mod_int_cv_phon_int_1 29 13305 13473 -6623.4    13247 2.2904      1     0.1302
# gca_mod_int_cv_phon_int_2 30 13305 13478 -6622.4    13245 1.9980      1     0.1575
# gca_mod_int_cv_phon_int_3 31 13307 13486 -6622.3    13245 0.2126      1     0.6447


# 4. For CVC targets, check interaction of WM:stress


gca_mod_int_cvc_phon_base <-
  lmer(eLog ~ 1 + (ot1 + ot2 + ot3) + phon_std + condition_sum +
         (1 + ot1 + ot2 + ot3 | participant) +
         (1 + ot1 + ot2 + ot3 | target),
       control = lmerControl(optimizer = 'bobyqa'), REML = F,
       data = filter(stress_gc_cvc, group == "int"))

gca_mod_int_cvc_phon_int_0 <- update(gca_mod_int_cvc_phon_base, . ~ . + phon_std:condition_sum)
gca_mod_int_cvc_phon_int_1 <- update(gca_mod_int_cvc_phon_int_0, . ~ . + ot1:phon_std:condition_sum)
gca_mod_int_cvc_phon_int_2 <- update(gca_mod_int_cvc_phon_int_1, . ~ . + ot2:phon_std:condition_sum)
gca_mod_int_cvc_phon_int_3 <- update(gca_mod_int_cvc_phon_int_2, . ~ . + ot3:phon_std:condition_sum)

int_cvc_phon_anova <-
anova(gca_mod_int_cvc_phon_int_0, gca_mod_int_cvc_phon_int_1, gca_mod_int_cvc_phon_int_2, gca_mod_int_cvc_phon_int_3)

#                            Df   AIC   BIC  logLik deviance  Chisq Chi Df Pr(>Chisq)
# gca_mod_int_cvc_phon_int_0 28 19676 19849 -9810.0    19620
# gca_mod_int_cvc_phon_int_1 29 19676 19854 -9808.8    19618 2.5484      1    0.11041
# gca_mod_int_cvc_phon_int_2 30 19677 19862 -9808.7    19617 0.1070      1    0.74354
# gca_mod_int_cvc_phon_int_3 31 19676 19867 -9806.8    19614 3.8081      1    0.05101 .



# -----------------------------------------------------------------------------







# Model predictions for plottting ---------------------------------------------

# Create design dataframe for predictions
new_dat_all <- stress_gc_subset %>%
  dplyr::select(group, time_zero, ot1:ot3, coda_sum, condition_sum, wm_std, phon_std, biphon_std) %>%
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
  mod_spec <- c("_base", "_parox_phon_base", "_parox_phon_int_0", "_parox_phon_int_1",
              "_parox_phon_int_2", "_parox_phon_int_3", "_ox_phon_int_0", "_ox_phon_int_1",
              "_ox_phon_int_2", "_ox_phon_int_3", "_cv_phon_int_0", "_cv_phon_int_1",
              "_cv_phon_int_2", "_cv_phon_int_3", "_cvc_phon_int_0", "_cvc_phon_int_1",
              "_cvc_phon_int_2", "_cvc_phon_int_3")

  # Store ind models in list
  ind_mods_phon <- mget(c(paste0(mod_type, "ss", mod_spec),
                        paste0(mod_type, "la", mod_spec),
                        paste0(mod_type, "int", mod_spec)))


  save(ind_mods_phon,
       file = here("models", "stress", "s3_adv_int_nat", "eye_track", "gca",
                   "ind_mods_phon.Rdata"))





  # Save anova model comparisons
  nested_model_comparisons_phon <-
    mget(c("ss_phon_anova", "ss_phon_anova_all", "ss_parox_phon_anova", "ss_ox_phon_anova",
      "ss_cv_phon_anova", "ss_cvc_phon_anova", "la_phon_anova", "la_parox_phon_anova",
      "la_ox_phon_anova", "la_cv_phon_anova", "la_cvc_phon_anova", "int_phon_anova",
      "int_parox_phon_anova", "int_ox_phon_anova", "int_cv_phon_anova", "int_cvc_phon_anova"))




  save(nested_model_comparisons_phon,
       file = here("models", "stress", "s3_adv_int_nat", "eye_track", "gca",
                   "nested_model_comparisons_phon.Rdata"))

  # Save models predictions
  model_preds <- mget(c("fits_all", "target_offset_preds"))

  save(model_preds,
       file = here("models", "stress", "s3_adv_int_nat", "eye_track", "gca",
                   "model_preds.Rdata"))

}

# -----------------------------------------------------------------------------

