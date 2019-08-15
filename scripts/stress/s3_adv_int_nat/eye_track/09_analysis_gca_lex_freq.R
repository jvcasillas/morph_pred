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
  select(., -coda, -freq)

freq_data <- read_excel(here("data", "raw", "lex_freq_stress.xlsx"))



stress50 <- stress50 %>%
  filter(., group %in% c("la", "int", "ss")) %>%
  left_join(mem_data, by = "participant") %>%
  left_join(phon_data, by = "target") %>%
  left_join(freq_data, by = "target")

stress50$wm <- as.numeric(stress50$wm)
glimpse(stress50)

# Get path to saved models
gca_mods_path  <- here("models", "stress", "s3_adv_int_nat", "eye_track", "gca")

# Load models as lists
load(paste0(gca_mods_path, "/ind_mods.Rdata"))
load(paste0(gca_mods_path, "/nested_model_comparisons_lex.Rdata"))
load(paste0(gca_mods_path, "/model_preds.Rdata"))

# Store objects in global env
list2env(ind_mods, globalenv())
list2env(nested_model_comparisons_lex, globalenv())
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
         biphon_std =(biphon_prob - mean(biphon_prob)) / sd(biphon_prob),
         nim_std = (freq_nim - mean(freq_nim)) / sd(freq_nim),
         davies_std = (freq_davies - mean(freq_davies)) / sd(freq_davies)) %>%
  poly_add_columns(., time_zero, degree = 3, prefix = "ot")

# -----------------------------------------------------------------------------


# -----------------------------------------------------------------------------

## Adding lexical frequency to full model
gca_lex_mod_int_0 <- update(gca_full_mod_int_3, . ~ . + nim_std)
gca_lex_mod_int_1 <- update(gca_lex_mod_int_0, . ~ . + ot1:nim_std)
gca_lex_mod_int_2 <- update(gca_lex_mod_int_1, . ~ . + ot2:nim_std)
gca_lex_mod_int_3 <- update(gca_lex_mod_int_2, . ~ . + ot3:nim_std)

full_lex_anova <-
  anova(gca_lex_mod_int_0, gca_lex_mod_int_1, gca_lex_mod_int_2, gca_lex_mod_int_3)

## Interactions
# lex freq x Group
gca_lex_group_mod_int_0 <- update(gca_full_mod_int_3, . ~ . + nim_std:group)
gca_lex_group_mod_int_1 <- update(gca_lex_group_mod_int_0, . ~ . + ot1:nim_std:group)
gca_lex_group_mod_int_2 <- update(gca_lex_group_mod_int_1, . ~ . + ot2:nim_std:group)
gca_lex_group_mod_int_3 <- update(gca_lex_group_mod_int_2, . ~ . + ot3:nim_std:group)

full_lex_group_anova <-
  anova(gca_lex_group_mod_int_0, gca_lex_group_mod_int_1, gca_lex_group_mod_int_2, gca_lex_group_mod_int_3)
# Nothing significant here

# lex freq x Stress

gca_lex_stress_mod_int_0 <- update(gca_full_mod_int_3, . ~ . + nim_std:condition_sum)
gca_lex_stress_mod_int_1 <- update(gca_lex_stress_mod_int_0, . ~ . + ot1:nim_std:condition_sum)
gca_lex_stress_mod_int_2 <- update(gca_lex_stress_mod_int_1, . ~ . + ot2:nim_std:condition_sum)
gca_lex_stress_mod_int_3 <- update(gca_lex_stress_mod_int_2, . ~ . + ot3:nim_std:condition_sum)

full_lex_stress_anova <-
  anova(gca_lex_stress_mod_int_0, gca_lex_stress_mod_int_1, gca_lex_stress_mod_int_2, gca_lex_stress_mod_int_3)



# lex freq x Coda
gca_lex_coda_mod_int_0 <- update(gca_full_mod_int_3, . ~ . + nim_std:coda_sum)
gca_lex_coda_mod_int_1 <- update(gca_lex_coda_mod_int_0, . ~ . + ot1:nim_std:coda_sum)
gca_lex_coda_mod_int_2 <- update(gca_lex_coda_mod_int_1, . ~ . + ot2:nim_std:coda_sum)
gca_lex_coda_mod_int_3 <- update(gca_lex_coda_mod_int_2, . ~ . + ot3:nim_std:coda_sum)

full_lex_coda_anova <-
  anova(gca_lex_coda_mod_int_0, gca_lex_coda_mod_int_1, gca_lex_coda_mod_int_2, gca_lex_coda_mod_int_3)







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



# add lex_freq_nim effect to intercept, linear slope, quadratic, and cubic time terms
gca_mod_ss_nim_0 <- update(gca_mod_ss_base,   . ~ . + nim_std) # lex freq as the only fixed effect
gca_mod_ss_nim_1 <- update(gca_mod_ss_nim_0,   . ~ . + ot1:nim_std)
gca_mod_ss_nim_2 <- update(gca_mod_ss_nim_1,   . ~ . + ot2:nim_std)
gca_mod_ss_nim_3 <- update(gca_mod_ss_nim_2,   . ~ . + ot3:nim_std)

ss_lex_anova <-
anova(gca_mod_ss_nim_0, gca_mod_ss_nim_1, gca_mod_ss_nim_2, gca_mod_ss_nim_3)

#                  Df   AIC   BIC logLik deviance  Chisq Chi Df Pr(>Chisq)
# gca_mod_ss_nim_0 37 34957 35208 -17442    34883
# gca_mod_ss_nim_1 38 34959 35217 -17442    34883 0.0148      1     0.9032
# gca_mod_ss_nim_2 39 34960 35224 -17441    34882 1.2800      1     0.2579
# gca_mod_ss_nim_3 40 34962 35233 -17441    34882 0.3033      1     0.5818


# add wm effect to intercept, linear slope, quadratic, and cubic time terms

gca_mod_ss_cond_3 <-
  lmer(eLog ~ 1 + (ot1 + ot2 + ot3) + condition_sum + coda_sum +
         (1 + coda_sum + condition_sum + (ot1 + ot2 + ot3) | participant) +
         (1 + ot1 + ot2 + ot3 | target),
       control = lmerControl(optimizer = 'bobyqa'), REML = F,
       data = filter(stress_gc_subset, group == "ss"))

gca_mod_ss_nim_0_all <- update(gca_mod_ss_cond_3,   . ~ . + nim_std) # Stress + coda + nim as fixed effects
gca_mod_ss_nim_1_all <- update(gca_mod_ss_nim_0_all,   . ~ . + ot1:nim_std)
gca_mod_ss_nim_2_all <- update(gca_mod_ss_nim_1_all,   . ~ . + ot2:nim_std)
gca_mod_ss_nim_3_all <- update(gca_mod_ss_nim_2_all,   . ~ . + ot3:nim_std)

ss_lex_anova_all <-
anova(gca_mod_ss_nim_0_all, gca_mod_ss_nim_1_all, gca_mod_ss_nim_2_all, gca_mod_ss_nim_3_all)

#                      Df   AIC   BIC logLik deviance  Chisq Chi Df Pr(>Chisq)
# gca_mod_ss_nim_0_all 45 34958 35262 -17434    34868
# gca_mod_ss_nim_1_all 46 34960 35271 -17434    34868 0.0356      1     0.8504
# gca_mod_ss_nim_2_all 47 34957 35275 -17432    34863 4.4597      1     0.0347 *
# gca_mod_ss_nim_3_all 48 34959 35284 -17432    34863 0.0000      1     0.9955



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
  lmer(eLog ~ 1 + (ot1 + ot2 + ot3) + nim_std + coda_sum +
         (1 + ot1 + ot2 + ot3 | participant) +
         (1 + ot1 + ot2 + ot3 | target),
       control = lmerControl(optimizer = 'bobyqa'), REML = F,
       data = filter(stress_gc_parox, group == "ss"))

gca_mod_ss_parox_int_0 <- update(gca_mod_ss_parox_base, . ~ . + nim_std:coda_sum)
gca_mod_ss_parox_int_1 <- update(gca_mod_ss_parox_int_0, . ~ . + ot1:nim_std:coda_sum)
gca_mod_ss_parox_int_2 <- update(gca_mod_ss_parox_int_1, . ~ . + ot2:nim_std:coda_sum)
gca_mod_ss_parox_int_3 <- update(gca_mod_ss_parox_int_2, . ~ . + ot3:nim_std:coda_sum)

ss_parox_lex_anova <-
anova(gca_mod_ss_parox_int_0, gca_mod_ss_parox_int_1, gca_mod_ss_parox_int_2, gca_mod_ss_parox_int_3)

#                     Df   AIC   BIC  logLik deviance  Chisq Chi Df Pr(>Chisq)
# gca_mod_parox_int_0 28 17308 17478 -8625.9    17252
# gca_mod_parox_int_1 29 17308 17484 -8625.2    17250 1.4600      1     0.2269
# gca_mod_parox_int_2 30 17310 17492 -8625.0    17250 0.3269      1     0.5675
# gca_mod_parox_int_3 31 17310 17498 -8624.2    17248 1.5676      1     0.2106



# 2. For Oxytone targets, check interaction of WM:Coda

stress_gc_ox <- stress_gc_subset %>%
  filter(., condition == "unstressed")

gca_mod_ss_ox_base <-
  lmer(eLog ~ 1 + (ot1 + ot2 + ot3) + nim_std + coda_sum +
         (1 + coda_sum + (ot1 + ot2 + ot3) | participant) +
         (1 + ot1 + ot2 + ot3 | target),
       control = lmerControl(optimizer = 'bobyqa'), REML = F,
       data = filter(stress_gc_ox, group == "ss"))

gca_mod_ss_ox_int_0 <- update(gca_mod_ss_ox_base, . ~ . + nim_std:coda_sum)
gca_mod_ss_ox_int_1 <- update(gca_mod_ss_ox_int_0, . ~ . + ot1:nim_std:coda_sum)
gca_mod_ss_ox_int_2 <- update(gca_mod_ss_ox_int_1, . ~ . + ot2:nim_std:coda_sum)
gca_mod_ss_ox_int_3 <- update(gca_mod_ss_ox_int_2, . ~ . + ot3:nim_std:coda_sum)

ss_ox_lex_anova <-
anova(gca_mod_ss_ox_int_0, gca_mod_ss_ox_int_1, gca_mod_ss_ox_int_2, gca_mod_ss_ox_int_3)

#                  Df   AIC   BIC  logLik deviance  Chisq Chi Df Pr(>Chisq)
# gca_mod_ox_int_0 33 17431 17632 -8682.5    17365
# gca_mod_ox_int_1 34 17432 17639 -8682.1    17364 0.8414      1     0.3590
# gca_mod_ox_int_2 35 17433 17646 -8681.5    17363 1.1953      1     0.2743
# gca_mod_ox_int_3 36 17434 17653 -8681.0    17362 1.0336      1     0.3093



# 3. For CV targets, check interaction of WM:stress

stress_gc_cv <- stress_gc_subset %>%
  filter(., coda == "0")

gca_mod_ss_cv_base <-
  lmer(eLog ~ 1 + (ot1 + ot2 + ot3) + nim_std + condition_sum +
         (1 + ot1 + ot2 + ot3 | participant) +
         (1 + ot1 + ot2 + ot3 | target),
       control = lmerControl(optimizer = 'bobyqa'), REML = F,
       data = filter(stress_gc_cv, group == "ss"))

gca_mod_ss_cv_int_0 <- update(gca_mod_ss_cv_base, . ~ . + nim_std:condition_sum)
gca_mod_ss_cv_int_1 <- update(gca_mod_ss_cv_int_0, . ~ . + ot1:nim_std:condition_sum)
gca_mod_ss_cv_int_2 <- update(gca_mod_ss_cv_int_1, . ~ . + ot2:nim_std:condition_sum)
gca_mod_ss_cv_int_3 <- update(gca_mod_ss_cv_int_2, . ~ . + ot3:nim_std:condition_sum)

ss_cv_lex_anova <-
anova(gca_mod_ss_cv_int_0, gca_mod_ss_cv_int_1, gca_mod_ss_cv_int_2, gca_mod_ss_cv_int_3)

#                  Df   AIC   BIC  logLik deviance  Chisq Chi Df Pr(>Chisq)
# gca_mod_cv_int_0 28 14721 14886 -7332.7    14665
# gca_mod_cv_int_1 29 14720 14891 -7330.8    14662 3.7509      1    0.05278 .
# gca_mod_cv_int_2 30 14722 14898 -7330.8    14662 0.0310      1    0.86033
# gca_mod_cv_int_3 31 14723 14906 -7330.7    14661 0.2999      1    0.58394

# 4. For CVC targets, check interaction of WM:stress

stress_gc_cvc <- stress_gc_subset %>%
  filter(., coda == "1")

gca_mod_ss_cvc_base <-
  lmer(eLog ~ 1 + (ot1 + ot2 + ot3) + nim_std + condition_sum +
         (1 + ot1 + ot2 + ot3 | participant) +
         (1 + ot1 + ot2 + ot3 | target),
       control = lmerControl(optimizer = 'bobyqa'), REML = F,
       data = filter(stress_gc_cvc, group == "ss"))

gca_mod_ss_cvc_int_0 <- update(gca_mod_ss_cvc_base, . ~ . + nim_std:condition_sum)
gca_mod_ss_cvc_int_1 <- update(gca_mod_ss_cvc_int_0, . ~ . + ot1:nim_std:condition_sum)
gca_mod_ss_cvc_int_2 <- update(gca_mod_ss_cvc_int_1, . ~ . + ot2:nim_std:condition_sum)
gca_mod_ss_cvc_int_3 <- update(gca_mod_ss_cvc_int_2, . ~ . + ot3:nim_std:condition_sum)

ss_cvc_lex_anova <-
anova(gca_mod_ss_cvc_int_0, gca_mod_ss_cvc_int_1, gca_mod_ss_cvc_int_2, gca_mod_ss_cvc_int_3)

#                   Df   AIC   BIC logLik deviance  Chisq Chi Df Pr(>Chisq)
# gca_mod_cvc_int_0 28 20558 20732 -10251    20502
# gca_mod_cvc_int_1 29 20559 20740 -10250    20501 0.4040      1     0.5251
# gca_mod_cvc_int_2 30 20559 20746 -10250    20499 2.1764      1     0.1401
# gca_mod_cvc_int_3 31 20561 20754 -10249    20499 0.0930      1     0.7605


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


# add nim effect to intercept, linear slope, quadratic, and cubic time terms
gca_mod_la_nim_0 <- update(gca_mod_la_base,   . ~ . + nim_std) # only WM as fixed effect
gca_mod_la_nim_1 <- update(gca_mod_la_nim_0,   . ~ . + ot1:nim_std)
gca_mod_la_nim_2 <- update(gca_mod_la_nim_1,   . ~ . + ot2:nim_std)
gca_mod_la_nim_3 <- update(gca_mod_la_nim_2,   . ~ . + ot3:nim_std)

la_lex_anova <-
anova(gca_mod_la_nim_0, gca_mod_la_nim_1, gca_mod_la_nim_2, gca_mod_la_nim_3)

#                   Df   AIC   BIC logLik deviance  Chisq Chi Df Pr(>Chisq)
# gca_mod_cvc_int_0 28 20558 20732 -10251    20502
# gca_mod_cvc_int_1 29 20559 20740 -10250    20501 0.4040      1     0.5251
# gca_mod_cvc_int_2 30 20559 20746 -10250    20499 2.1764      1     0.1401
# gca_mod_cvc_int_3 31 20561 20754 -10249    20499 0.0930      1     0.7605

-----------------------------------------------------------------------------

###########################
# Interactions:           #
# 1. Ox -> WM:coda        #
# 2. Parox -> WM:coda     #
# 3. CV -> WM:stress      #
# 4. CVC -> WM:stress     #
###########################

# 1. For Paroxytone targets, check interaction of WM:Coda



gca_mod_la_parox_base<-lmer(eLog ~ 1 + (ot1 + ot2 + ot3) + nim_std + coda_sum +
         (1 + coda_sum + (ot1 + ot2 + ot3) | participant) +
         (1 + ot1 + ot2 + ot3 | target),
       control = lmerControl(optimizer = 'bobyqa'), REML = F,
       data = filter(stress_gc_parox, group == "la"))

gca_mod_la_parox_int_0 <- update(gca_mod_la_parox_base, . ~ . + nim_std:coda_sum)
gca_mod_la_parox_int_1 <- update(gca_mod_la_parox_int_0, . ~ . + ot1:nim_std:coda_sum)
gca_mod_la_parox_int_2 <- update(gca_mod_la_parox_int_1, . ~ . + ot2:nim_std:coda_sum)
gca_mod_la_parox_int_3 <- update(gca_mod_la_parox_int_2, . ~ . + ot3:nim_std:coda_sum)

la_parox_lex_anova <-
anova(gca_mod_la_parox_int_0, gca_mod_la_parox_int_1, gca_mod_la_parox_int_2, gca_mod_la_parox_int_3)




# 2. For Oxytone targets, check interaction of WM:Coda

gca_mod_la_ox_base <-
  lmer(eLog ~ 1 + (ot1 + ot2 + ot3) + nim_std + coda_sum +
         (1 + coda_sum + (ot1 + ot2 + ot3) | participant) +
         (1 + ot1 + ot2 + ot3 | target),
       control = lmerControl(optimizer = 'bobyqa'), REML = F,
       data = filter(stress_gc_ox, group == "la"))

gca_mod_la_ox_int_0 <- update(gca_mod_la_ox_base, . ~ . + nim_std:coda_sum)
gca_mod_la_ox_int_1 <- update(gca_mod_la_ox_int_0, . ~ . + ot1:nim_std:coda_sum)
gca_mod_la_ox_int_2 <- update(gca_mod_la_ox_int_1, . ~ . + ot2:nim_std:coda_sum)
gca_mod_la_ox_int_3 <- update(gca_mod_la_ox_int_2, . ~ . + ot3:nim_std:coda_sum)

la_ox_lex_anova <-
anova(gca_mod_la_ox_int_0, gca_mod_la_ox_int_1, gca_mod_la_ox_int_2, gca_mod_la_ox_int_3)

#                     Df   AIC   BIC  logLik deviance  Chisq Chi Df Pr(>Chisq)
# gca_mod_la_ox_int_0 33 18221 18423 -9077.5    18155
# gca_mod_la_ox_int_1 34 18223 18430 -9077.3    18155 0.2408      1     0.6237
# gca_mod_la_ox_int_2 35 18224 18438 -9077.1    18154 0.4427      1     0.5058
# gca_mod_la_ox_int_3 36 18226 18446 -9076.9    18154 0.3529      1     0.5525


# 3. For CV targets, check interaction of WM:stress

gca_mod_la_cv_base <-
  lmer(eLog ~ 1 + (ot1 + ot2 + ot3) + nim_std + condition_sum +
         (1 + ot1 + ot2 + ot3 | participant) +
         (1 + ot1 + ot2 + ot3 | target),
       control = lmerControl(optimizer = 'bobyqa'), REML = F,
       data = filter(stress_gc_cv, group == "la"))

gca_mod_la_cv_int_0 <- update(gca_mod_la_cv_base, . ~ . + nim_std:condition_sum)
gca_mod_la_cv_int_1 <- update(gca_mod_la_cv_int_0, . ~ . + ot1:nim_std:condition_sum)
gca_mod_la_cv_int_2 <- update(gca_mod_la_cv_int_1, . ~ . + ot2:nim_std:condition_sum)
gca_mod_la_cv_int_3 <- update(gca_mod_la_cv_int_2, . ~ . + ot3:nim_std:condition_sum)

la_cv_lex_anova <-
anova(gca_mod_la_cv_int_0, gca_mod_la_cv_int_1, gca_mod_la_cv_int_2, gca_mod_la_cv_int_3)

#                     Df   AIC   BIC  logLik deviance  Chisq Chi Df Pr(>Chisq)
# gca_mod_la_cv_int_0 28 15126 15291 -7534.8    15070
# gca_mod_la_cv_int_1 29 15128 15299 -7534.7    15070 0.0861      1     0.7692
# gca_mod_la_cv_int_2 30 15129 15306 -7534.7    15069 0.0874      1     0.7675
# gca_mod_la_cv_int_3 31 15130 15313 -7534.1    15068 1.1799      1     0.2774


# 4. For CVC targets, check interaction of WM:stress


gca_mod_la_cvc_base <-
  lmer(eLog ~ 1 + (ot1 + ot2 + ot3) + nim_std + condition_sum +
         (1 + ot1 + ot2 + ot3 | participant) +
         (1 + ot1 + ot2 + ot3 | target),
       control = lmerControl(optimizer = 'bobyqa'), REML = F,
       data = filter(stress_gc_cvc, group == "la"))

gca_mod_la_cvc_int_0 <- update(gca_mod_la_cvc_base, . ~ . + nim_std:condition_sum)
gca_mod_la_cvc_int_1 <- update(gca_mod_la_cvc_int_0, . ~ . + ot1:nim_std:condition_sum)
gca_mod_la_cvc_int_2 <- update(gca_mod_la_cvc_int_1, . ~ . + ot2:nim_std:condition_sum)
gca_mod_la_cvc_int_3 <- update(gca_mod_la_cvc_int_2, . ~ . + ot3:nim_std:condition_sum)

la_cvc_lex_anova <-
anova(gca_mod_la_cvc_int_0, gca_mod_la_cvc_int_1, gca_mod_la_cvc_int_2, gca_mod_la_cvc_int_3)

#                       Df   AIC   BIC logLik deviance  Chisq Chi Df Pr(>Chisq)
# gca_mod_la_cvc_int_0 28 21462 21638 -10703    21406
# gca_mod_la_cvc_int_1 29 21462 21643 -10702    21404 2.5946      1     0.1072
# gca_mod_la_cvc_int_2 30 21464 21652 -10702    21404 0.0022      1     0.9628
# gca_mod_la_cvc_int_3 31 21464 21659 -10701    21402 1.3817      1     0.2398



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



# add nim effect to intercept, linear slope, quadratic, and cubic time terms

gca_mod_int_nim_0 <- update(gca_mod_int_base,   . ~ . + nim_std) # only nim as fixed effect
gca_mod_int_nim_1 <- update(gca_mod_int_nim_0,   . ~ . + ot1:nim_std)
gca_mod_int_nim_2 <- update(gca_mod_int_nim_1,   . ~ . + ot2:nim_std)
gca_mod_int_nim_3 <- update(gca_mod_int_nim_2,   . ~ . + ot3:nim_std)

int_lex_anova <-
anova(gca_mod_int_nim_0, gca_mod_int_nim_1, gca_mod_int_nim_2, gca_mod_int_nim_3)

#                   Df   AIC   BIC logLik deviance  Chisq Chi Df Pr(>Chisq)
# gca_mod_int_nim_0 37 32833 33080 -16379    32759
# gca_mod_int_nim_1 38 32834 33089 -16379    32758 0.3690      1     0.5435
# gca_mod_int_nim_2 39 32836 33097 -16379    32758 0.2109      1     0.6461
# gca_mod_int_nim_3 40 32837 33105 -16379    32757 1.0039      1     0.3164


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
  lmer(eLog ~ 1 + (ot1 + ot2 + ot3) + nim_std + coda_sum +
         (1 + ot1 + ot2 + ot3 | participant) +
         (1 + ot1 + ot2 + ot3 | target),
       control = lmerControl(optimizer = 'bobyqa'), REML = F,
       data = filter(stress_gc_parox, group == "int"))

gca_mod_int_parox_int_0 <- update(gca_mod_int_parox_base, . ~ . + nim_std:coda_sum)
gca_mod_int_parox_int_1 <- update(gca_mod_int_parox_int_0, . ~ . + ot1:nim_std:coda_sum)
gca_mod_int_parox_int_2 <- update(gca_mod_int_parox_int_1, . ~ . + ot2:nim_std:coda_sum)
gca_mod_int_parox_int_3 <- update(gca_mod_int_parox_int_2, . ~ . + ot3:nim_std:coda_sum)

int_parox_lex_anova <-
anova(gca_mod_int_parox_int_0, gca_mod_int_parox_int_1, gca_mod_int_parox_int_2, gca_mod_int_parox_int_3)

#                         Df   AIC   BIC  logLik deviance  Chisq Chi Df Pr(>Chisq)
# gca_mod_int_parox_int_0 28 16396 16564 -8170.1    16340
# gca_mod_int_parox_int_1 29 16397 16571 -8169.6    16339 0.9659      1     0.3257
# gca_mod_int_parox_int_2 30 16399 16579 -8169.6    16339 0.0136      1     0.9073
# gca_mod_int_parox_int_3 31 16401 16586 -8169.3    16339 0.6268      1     0.4285

# 2. For Oxytone targets, check interaction of WM:Coda

gca_mod_int_ox_base <-
  lmer(eLog ~ 1 + (ot1 + ot2 + ot3) + nim_std + coda_sum +
         (1 + coda_sum + (ot1 + ot2 + ot3) | participant) +
         (1 + ot1 + ot2 + ot3 | target),
       control = lmerControl(optimizer = 'bobyqa'), REML = F,
       data = filter(stress_gc_ox, group == "int"))

gca_mod_int_ox_int_0 <- update(gca_mod_int_ox_base, . ~ . + nim_std:coda_sum)
gca_mod_int_ox_int_1 <- update(gca_mod_int_ox_int_0, . ~ . + ot1:nim_std:coda_sum)
gca_mod_int_ox_int_2 <- update(gca_mod_int_ox_int_1, . ~ . + ot2:nim_std:coda_sum)
gca_mod_int_ox_int_3 <- update(gca_mod_int_ox_int_2, . ~ . + ot3:nim_std:coda_sum)

int_ox_lex_anova <-
anova(gca_mod_int_ox_int_0, gca_mod_int_ox_int_1, gca_mod_int_ox_int_2, gca_mod_int_ox_int_3)

#                       Df   AIC   BIC  logLik deviance  Chisq Chi Df Pr(>Chisq)
# gca_mod_int_ox_int_0 33 16339 16537 -8136.3    16273
# gca_mod_int_ox_int_1 34 16340 16544 -8135.9    16272 0.7798      1     0.3772
# gca_mod_int_ox_int_2 35 16341 16551 -8135.5    16271 0.8343      1     0.3610
# gca_mod_int_ox_int_3 36 16342 16558 -8134.8    16270 1.3198      1     0.2506


# 3. For CV targets, check interaction of WM:stress

gca_mod_int_cv_base <-
  lmer(eLog ~ 1 + (ot1 + ot2 + ot3) + nim_std + condition_sum +
         (1 + ot1 + ot2 + ot3 | participant) +
         (1 + ot1 + ot2 + ot3 | target),
       control = lmerControl(optimizer = 'bobyqa'), REML = F,
       data = filter(stress_gc_cv, group == "int"))

gca_mod_int_cv_int_0 <- update(gca_mod_int_cv_base, . ~ . + nim_std:condition_sum)
gca_mod_int_cv_int_1 <- update(gca_mod_int_cv_int_0, . ~ . + ot1:nim_std:condition_sum)
gca_mod_int_cv_int_2 <- update(gca_mod_int_cv_int_1, . ~ . + ot2:nim_std:condition_sum)
gca_mod_int_cv_int_3 <- update(gca_mod_int_cv_int_2, . ~ . + ot3:nim_std:condition_sum)

int_cv_lex_anova <-
anova(gca_mod_int_cv_int_0, gca_mod_int_cv_int_1, gca_mod_int_cv_int_2, gca_mod_int_cv_int_3)

#                       Df   AIC   BIC  logLik deviance  Chisq Chi Df Pr(>Chisq)
# gca_mod_int_cv_int_0 28 13304 13466 -6624.1    13248
# gca_mod_int_cv_int_1 29 13303 13471 -6622.6    13245 2.8463      1    0.09159 .
# gca_mod_int_cv_int_2 30 13305 13479 -6622.5    13245 0.3004      1    0.58361
# gca_mod_int_cv_int_3 31 13304 13484 -6621.2    13242 2.5179      1    0.11256


# 4. For CVC targets, check interaction of WM:stress


gca_mod_int_cvc_base <-
  lmer(eLog ~ 1 + (ot1 + ot2 + ot3) + nim_std + condition_sum +
         (1 + ot1 + ot2 + ot3 | participant) +
         (1 + ot1 + ot2 + ot3 | target),
       control = lmerControl(optimizer = 'bobyqa'), REML = F,
       data = filter(stress_gc_cvc, group == "int"))

gca_mod_int_cvc_int_0 <- update(gca_mod_int_cvc_base, . ~ . + nim_std:condition_sum)
gca_mod_int_cvc_int_1 <- update(gca_mod_int_cvc_int_0, . ~ . + ot1:nim_std:condition_sum)
gca_mod_int_cvc_int_2 <- update(gca_mod_int_cvc_int_1, . ~ . + ot2:nim_std:condition_sum)
gca_mod_int_cvc_int_3 <- update(gca_mod_int_cvc_int_2, . ~ . + ot3:nim_std:condition_sum)

int_cvc_lex_anova <-
anova(gca_mod_int_cvc_int_0, gca_mod_int_cvc_int_1, gca_mod_int_cvc_int_2, gca_mod_int_cvc_int_3)

#                       Df   AIC   BIC  logLik deviance  Chisq Chi Df Pr(>Chisq)
# gca_mod_int_cvc_int_0 28 19674 19846 -9808.8    19618
# gca_mod_int_cvc_int_1 29 19675 19854 -9808.7    19617 0.2695      1     0.6037
# gca_mod_int_cvc_int_2 30 19676 19861 -9808.2    19616 1.0605      1     0.3031
# gca_mod_int_cvc_int_3 31 19678 19870 -9808.2    19616 0.0011      1     0.9734



# -----------------------------------------------------------------------------










# Model predictions for plottting ---------------------------------------------

# Create design dataframe for predictions
new_dat_all <- stress_gc_subset %>%
  dplyr::select(group, time_zero, ot1:ot3, coda_sum, condition_sum, nim_std) %>%
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
  ind_mods_lex <- mget(c(paste0(mod_type, "ss", mod_spec),
                          paste0(mod_type, "la", mod_spec),
                          paste0(mod_type, "int", mod_spec)))


  save(ind_mods_lex,
       file = here("models", "stress", "s3_adv_int_nat", "eye_track", "gca",
                   "ind_mods_lex.Rdata"))



  # Save anova model comparisons
  nested_model_comparisons_lex <-
    mget(c("ss_lex_anova", "ss_lex_anova_all", "ss_parox_lex_anova", "ss_ox_lex_anova",
      "ss_cv_lex_anova", "ss_cvc_lex_anova", "la_lex_anova", "la_parox_lex_anova",
      "la_ox_lex_anova", "la_cv_lex_anova", "la_cvc_lex_anova", "int_lex_anova",
      "int_parox_lex_anova", "int_ox_lex_anova", "int_cv_lex_anova", "int_cvc_lex_anova"))




  save(nested_model_comparisons_lex,
       file = here("models", "stress", "s3_adv_int_nat", "eye_track", "gca",
                   "nested_model_comparisons_lex.Rdata"))

  # Save models predictions
  model_preds <- mget(c("fits_all", "target_offset_preds"))

  save(model_preds,
       file = here("models", "stress", "s3_adv_int_nat", "eye_track", "gca",
                   "model_preds.Rdata"))


}

# -----------------------------------------------------------------------------

