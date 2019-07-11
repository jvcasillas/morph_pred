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


#                 Df   AIC   BIC logLik deviance  Chisq Chi Df Pr(>Chisq)
# gca_mod_ss_wm_0 37 33446 33695 -16686    33372
# gca_mod_ss_wm_1 38 33448 33704 -16686    33372 0.0347      1    0.85226
# gca_mod_ss_wm_2 39 33448 33710 -16685    33370 2.8118      1    0.09357 .
# gca_mod_ss_wm_3 40 33450 33719 -16685    33370 0.0895      1    0.76483

# add wm effect to intercept, linear slope, quadratic, and cubic time terms
gca_mod_ss_wm_0_all <- update(gca_mod_ss_cond_3,   . ~ . + wm_std)
gca_mod_ss_wm_1_all <- update(gca_mod_ss_wm_0_all,   . ~ . + ot1:wm_std)
gca_mod_ss_wm_2_all <- update(gca_mod_ss_wm_1_all,   . ~ . + ot2:wm_std)
gca_mod_ss_wm_3_all <- update(gca_mod_ss_wm_2_all,   . ~ . + ot3:wm_std)

anova(gca_mod_ss_wm_0_all, gca_mod_ss_wm_1_all, gca_mod_ss_wm_2_all, gca_mod_ss_wm_3_all)

#                     Df   AIC   BIC logLik deviance  Chisq Chi Df Pr(>Chisq)
# gca_mod_ss_wm_0_all 45 33448 33751 -16679    33358
# gca_mod_ss_wm_1_all 46 33450 33760 -16679    33358 0.0405      1    0.84042
# gca_mod_ss_wm_2_all 47 33449 33765 -16678    33355 2.7777      1    0.09558 .
# gca_mod_ss_wm_3_all 48 33451 33774 -16678    33355 0.0738      1    0.78594

# add wm x coda int to intercept, linear slope, quadratic, and cubic terms
gca_mod_ss_wm_0_all_int <- update(gca_mod_ss_wm_3_all,   . ~ . + wm_std:coda_sum)
gca_mod_ss_wm_1_all_int <- update(gca_mod_ss_wm_0_all,   . ~ . + ot1:wm_std:coda_sum)
gca_mod_ss_wm_2_all_int <- update(gca_mod_ss_wm_1_all,   . ~ . + ot2:wm_std:coda_sum)
gca_mod_ss_wm_3_all_int <- update(gca_mod_ss_wm_2_all,   . ~ . + ot3:wm_std:coda_sum)

anova(gca_mod_ss_wm_0_all_int, gca_mod_ss_wm_1_all_int, gca_mod_ss_wm_2_all_int, gca_mod_ss_wm_3_all_int)


#                         Df   AIC   BIC logLik deviance  Chisq Chi Df Pr(>Chisq)
# gca_mod_ss_wm_1_all_int 46 33449 33758 -16678    33357
# gca_mod_ss_wm_2_all_int 47 33452 33768 -16679    33358 0.0000      1    1.00000
# gca_mod_ss_wm_3_all_int 48 33449 33772 -16676    33353 5.3591      1    0.02061 *
# gca_mod_ss_wm_0_all_int 49 33452 33782 -16677    33354 0.0000      1    1.00000

# check only interaction between working memory and coda (without condition parox/ox)

gca_mod_ss_base_wm_coda_int <-
  lmer(eLog ~ 1 + (ot1 + ot2 + ot3) + wm_std + coda_sum +
         (1 + coda_sum + condition_sum + (ot1 + ot2 + ot3) | participant) +
         (1 + ot1 + ot2 + ot3 | target),
       control = lmerControl(optimizer = 'bobyqa'), REML = F,
       data = filter(stress_gc_subset, group == "ss"))

# add coda effect to intercept, linear slope, quadratic, and cubic time terms
gca_mod_ss_base_wm_coda_int_0 <- update(gca_mod_ss_base_wm_coda_int,   . ~ . + wm_std:coda_sum)
gca_mod_ss_base_wm_coda_int_1 <- update(gca_mod_ss_coda_0, . ~ . + ot1:wm_std:coda_sum)
gca_mod_ss_base_wm_coda_int_2 <- update(gca_mod_ss_coda_1, . ~ . + ot2:wm_std:coda_sum)
gca_mod_ss_base_wm_coda_int_3 <- update(gca_mod_ss_coda_2, . ~ . + ot3:wm_std:coda_sum)

anova(gca_mod_ss_base_wm_coda_int_0, gca_mod_ss_base_wm_coda_int_1,
      gca_mod_ss_base_wm_coda_int_2, gca_mod_ss_base_wm_coda_int_3)

#                               Df   AIC   BIC logLik deviance  Chisq Chi Df Pr(>Chisq)
# gca_mod_ss_base_wm_coda_int_1 38 33448 33704 -16686    33372
# gca_mod_ss_base_wm_coda_int_0 39 33449 33711 -16686    33371 1.1363      1    0.28644
# gca_mod_ss_base_wm_coda_int_2 39 33450 33712 -16686    33372 0.0000      0    1.00000
# gca_mod_ss_base_wm_coda_int_3 40 33447 33716 -16684    33367 4.4172      1    0.03558 *



# check only interaction between working memory and cond parox/ox (without coda)

gca_mod_ss_base_wm_cond_int <-
  lmer(eLog ~ 1 + (ot1 + ot2 + ot3) + wm_std + condition_sum +
         (1 + coda_sum + condition_sum + (ot1 + ot2 + ot3) | participant) +
         (1 + ot1 + ot2 + ot3 | target),
       control = lmerControl(optimizer = 'bobyqa'), REML = F,
       data = filter(stress_gc_subset, group == "ss"))

# add coda effect to intercept, linear slope, quadratic, and cubic time terms
gca_mod_ss_base_wm_cond_int_0 <- update(gca_mod_ss_base_wm_cond_int,   . ~ . + wm_std:condition_sum)
gca_mod_ss_base_wm_cond_int_1 <- update(gca_mod_ss_coda_0, . ~ . + ot1:wm_std:condition_sum)
gca_mod_ss_base_wm_cond_int_2 <- update(gca_mod_ss_coda_1, . ~ . + ot2:wm_std:condition_sum)
gca_mod_ss_base_wm_cond_int_3 <- update(gca_mod_ss_coda_2, . ~ . + ot3:wm_std:condition_sum)

anova(gca_mod_ss_base_wm_cond_int_0, gca_mod_ss_base_wm_cond_int_1,
      gca_mod_ss_base_wm_cond_int_2, gca_mod_ss_base_wm_cond_int_3)

#                               Df   AIC   BIC logLik deviance   Chisq Chi Df Pr(>Chisq)
# gca_mod_ss_base_wm_cond_int_1 38 33449 33704 -16686    33373
# gca_mod_ss_base_wm_cond_int_0 39 33448 33711 -16685    33370  2.3885      1     0.1222
# gca_mod_ss_base_wm_cond_int_2 39 33436 33698 -16679    33358 12.5150      0     <2e-16 ***
# gca_mod_ss_base_wm_cond_int_3 40 33451 33720 -16685    33371  0.0000      1     1.0000

summary(gca_mod_ss_base_wm_cond_int_2)

# add phonotactic freq as a variable

gca_mod_ss_phon_0 <- update(gca_mod_ss_base,   . ~ . + phon_std)
gca_mod_ss_phon_1 <- update(gca_mod_ss_phon_0,   . ~ . + ot1:phon_std)
gca_mod_ss_phon_2 <- update(gca_mod_ss_phon_1,   . ~ . + ot2:phon_std)
gca_mod_ss_phon_3 <- update(gca_mod_ss_phon_2,   . ~ . + ot3:phon_std)

anova(gca_mod_ss_phon_0, gca_mod_ss_phon_1, gca_mod_ss_phon_2, gca_mod_ss_phon_3)

#                   Df   AIC   BIC logLik deviance  Chisq Chi Df Pr(>Chisq)
# gca_mod_ss_phon_0 37 34957 35208 -17442    34883
# gca_mod_ss_phon_1 38 34959 35216 -17442    34883 0.3016      1     0.5829
# gca_mod_ss_phon_2 39 34959 35223 -17440    34881 2.2483      1     0.1338
# gca_mod_ss_phon_3 40 34961 35231 -17440    34881 0.2734      1     0.6011

# add biphon freq as a variable

gca_mod_ss_biphon_0 <- update(gca_mod_ss_base,   . ~ . + biphon_std)
gca_mod_ss_biphon_1 <- update(gca_mod_ss_biphon_0,   . ~ . + ot1:biphon_std)
gca_mod_ss_biphon_2 <- update(gca_mod_ss_biphon_1,   . ~ . + ot2:biphon_std)
gca_mod_ss_biphon_3 <- update(gca_mod_ss_biphon_2,   . ~ . + ot3:biphon_std)

anova(gca_mod_ss_biphon_0, gca_mod_ss_biphon_1, gca_mod_ss_biphon_2, gca_mod_ss_biphon_3)

#                     Df   AIC   BIC logLik deviance  Chisq Chi Df Pr(>Chisq)
# gca_mod_ss_biphon_0 37 34957 35208 -17442    34883
# gca_mod_ss_biphon_1 38 34957 35215 -17441    34881 1.7941      1     0.1804
# gca_mod_ss_biphon_2 39 34959 35223 -17441    34881 0.0053      1     0.9419
# gca_mod_ss_biphon_3 40 34961 35231 -17440    34881 0.7252      1     0.3945



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
gca_mod_la_wm_3 <- update(gca_mod_la_wm_2,   . ~ . + ot3:wm_std)

anova(gca_mod_la_wm_0, gca_mod_la_wm_1, gca_mod_la_wm_2, gca_mod_la_wm_3)

#                 Df   AIC   BIC logLik deviance  Chisq Chi Df Pr(>Chisq)
# gca_mod_la_wm_0 37 34943 35193 -17434    34869
# gca_mod_la_wm_1 38 34945 35201 -17434    34869 0.1543      1     0.6945
# gca_mod_la_wm_2 39 34947 35210 -17434    34869 0.2409      1     0.6236
# gca_mod_la_wm_3 40 34948 35219 -17434    34868 0.0716      1     0.7890

# add phonotactic freq as a variable

gca_mod_la_phon_0 <- update(gca_mod_la_base,   . ~ . + phon_std)
gca_mod_la_phon_1 <- update(gca_mod_la_phon_0,   . ~ . + ot1:phon_std)
gca_mod_la_phon_2 <- update(gca_mod_la_phon_1,   . ~ . + ot2:phon_std)
gca_mod_la_phon_3 <- update(gca_mod_la_phon_2,   . ~ . + ot3:phon_std)

anova(gca_mod_la_phon_0, gca_mod_la_phon_1, gca_mod_la_phon_2, gca_mod_la_phon_3)

#                   Df   AIC   BIC logLik deviance  Chisq Chi Df Pr(>Chisq)
# gca_mod_la_phon_0 37 36297 36548 -18111    36223
# gca_mod_la_phon_1 38 36299 36557 -18111    36223 0.0226      1     0.8804
# gca_mod_la_phon_2 39 36301 36566 -18111    36223 0.0746      1     0.7848
# gca_mod_la_phon_3 40 36302 36574 -18111    36222 0.5824      1     0.4454


# add biphon freq as a variable

gca_mod_la_biphon_0 <- update(gca_mod_la_base,   . ~ . + biphon_std)
gca_mod_la_biphon_1 <- update(gca_mod_la_biphon_0,   . ~ . + ot1:biphon_std)
gca_mod_la_biphon_2 <- update(gca_mod_la_biphon_1,   . ~ . + ot2:biphon_std)
gca_mod_la_biphon_3 <- update(gca_mod_la_biphon_2,   . ~ . + ot3:biphon_std)

anova(gca_mod_la_biphon_0, gca_mod_la_biphon_1, gca_mod_la_biphon_2, gca_mod_la_biphon_3)

#                     Df   AIC   BIC logLik deviance  Chisq Chi Df Pr(>Chisq)
# gca_mod_la_biphon_0 37 36296 36547 -18111    36222
# gca_mod_la_biphon_1 38 36297 36555 -18110    36221 1.0234      1     0.3117
# gca_mod_la_biphon_2 39 36299 36564 -18110    36221 0.0358      1     0.8499
# gca_mod_la_biphon_3 40 36299 36570 -18109    36219 2.2243      1     0.1359

# Check PSTM
stress_la_pstm_subset <- stress_gc_subset %>%
  filter(., group == "la") %>%
  mutate(., pstm_std = (age - mean(pstm)) / sd(pstm))

gca_mod_la_base_pstm <-
  lmer(eLog ~ 1 + (ot1 + ot2 + ot3) +
         (1 + coda_sum + condition_sum + (ot1 + ot2 + ot3) | participant) +
         (1 + ot1 + ot2 + ot3 | target),
       control = lmerControl(optimizer = 'bobyqa'), REML = F,
       data = stress_la_pstm_subset)

gca_mod_la_pstm_0 <- update(gca_mod_la_base_pstm,   . ~ . + pstm_std)
gca_mod_la_pstm_1 <- update(gca_mod_la_pstm_0,   . ~ . + ot1:pstm_std)
gca_mod_la_pstm_2 <- update(gca_mod_la_pstm_1,   . ~ . + ot2:pstm_std)
gca_mod_la_pstm_3 <- update(gca_mod_la_pstm_2,   . ~ . + ot3:pstm_std)

anova(gca_mod_la_pstm_0, gca_mod_la_pstm_1, gca_mod_la_pstm_2, gca_mod_la_pstm_3)

#                   Df   AIC   BIC logLik deviance  Chisq Chi Df Pr(>Chisq)
# gca_mod_la_pstm_0 37 36300 36552 -18113    36226
# gca_mod_la_pstm_1 38 36301 36559 -18112    36225 1.5691      1     0.2103
# gca_mod_la_pstm_2 39 36301 36566 -18111    36223 1.7927      1     0.1806
# gca_mod_la_pstm_3 40 36303 36574 -18111    36223 0.0264      1     0.8708


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


# add wm effect to intercept, linear slope, quadratic, and cubic time terms

gca_mod_int_wm_0 <- update(gca_mod_int_base,   . ~ . + wm_std)
gca_mod_int_wm_1 <- update(gca_mod_int_wm_0,   . ~ . + ot1:wm_std)
gca_mod_int_wm_2 <- update(gca_mod_int_wm_1,   . ~ . + ot2:wm_std)
gca_mod_int_wm_3 <- update(gca_mod_int_wm_2,   . ~ . + ot3:wm_std)

anova(gca_mod_int_wm_0, gca_mod_int_wm_1, gca_mod_int_wm_2, gca_mod_int_wm_3)

#                   Df   AIC   BIC logLik deviance  Chisq Chi Df Pr(>Chisq)
# gca_mod_int_wm_0 37 32834 33081 -16380    32760
# gca_mod_int_wm_1 38 32835 33089 -16380    32759 0.4050      1     0.5245
# gca_mod_int_wm_2 39 32837 33098 -16380    32759 0.2694      1     0.6037
# gca_mod_int_wm_3 40 32838 33105 -16379    32758 1.3810      1     0.2399


# add phonotactic freq as a variable

gca_mod_int_phon_0 <- update(gca_mod_int_base,   . ~ . + phon_std)
# Warning message:
#   In optwrap(optimizer, devfun, getStart(start, rho$lower, rho$pp),  :
#                convergence code 1 from bobyqa: bobyqa -- maximum number of function evaluations exceeded
gca_mod_int_phon_1 <- update(gca_mod_int_phon_0,   . ~ . + ot1:phon_std)
gca_mod_int_phon_2 <- update(gca_mod_int_phon_1,   . ~ . + ot2:phon_std)
gca_mod_int_phon_3 <- update(gca_mod_int_phon_2,   . ~ . + ot3:phon_std)

anova(gca_mod_int_phon_0, gca_mod_int_phon_1, gca_mod_int_phon_2, gca_mod_int_phon_3)

#                    Df   AIC   BIC logLik deviance  Chisq Chi Df Pr(>Chisq)
# gca_mod_int_phon_0 37 32835 33082 -16380    32761
# gca_mod_int_phon_1 38 32837 33091 -16380    32761 0.0375      1     0.8464
# gca_mod_int_phon_2 39 32838 33099 -16380    32760 0.6448      1     0.4220
# gca_mod_int_phon_3 40 32840 33107 -16380    32760 0.4125      1     0.5207



# add biphon freq as a variable

gca_mod_int_biphon_0 <- update(gca_mod_int_base,   . ~ . + biphon_std)
gca_mod_int_biphon_1 <- update(gca_mod_int_biphon_0,   . ~ . + ot1:biphon_std)
gca_mod_int_biphon_2 <- update(gca_mod_int_biphon_1,   . ~ . + ot2:biphon_std)
gca_mod_int_biphon_3 <- update(gca_mod_int_biphon_2,   . ~ . + ot3:biphon_std)

anova(gca_mod_int_biphon_0, gca_mod_int_biphon_1, gca_mod_int_biphon_2, gca_mod_int_biphon_3)

#                      Df   AIC   BIC logLik deviance  Chisq Chi Df Pr(>Chisq)
# gca_mod_int_biphon_0 37 32835 33082 -16380    32761
# gca_mod_int_biphon_1 38 32837 33091 -16380    32761 0.0046      1     0.9460
# gca_mod_int_biphon_2 39 32837 33098 -16380    32759 1.5273      1     0.2165
# gca_mod_int_biphon_3 40 32838 33106 -16379    32758 1.0466      1     0.3063

# Check PSTM
stress_int_pstm_subset <- stress_gc_subset %>%
  filter(., group == "int") %>%
  mutate(., pstm_std = (age - mean(pstm)) / sd(pstm))

gca_mod_int_base_pstm <-
  lmer(eLog ~ 1 + (ot1 + ot2 + ot3) +
         (1 + coda_sum + condition_sum + (ot1 + ot2 + ot3) | participant) +
         (1 + ot1 + ot2 + ot3 | target),
       control = lmerControl(optimizer = 'bobyqa'), REML = F,
       data = stress_int_pstm_subset)

gca_mod_int_pstm_0 <- update(gca_mod_int_base_pstm,   . ~ . + pstm_std)
gca_mod_int_pstm_1 <- update(gca_mod_int_pstm_0,   . ~ . + ot1:pstm_std)
gca_mod_int_pstm_2 <- update(gca_mod_int_pstm_1,   . ~ . + ot2:pstm_std)
gca_mod_int_pstm_3 <- update(gca_mod_int_pstm_2,   . ~ . + ot3:pstm_std)

anova(gca_mod_int_pstm_0, gca_mod_int_pstm_1, gca_mod_int_pstm_2, gca_mod_int_pstm_3)

#                    Df   AIC   BIC logLik deviance  Chisq Chi Df Pr(>Chisq)
# gca_mod_int_pstm_0 37 32835 33082 -16380    32761
# gca_mod_int_pstm_1 38 32837 33091 -16380    32761 0.2361      1     0.6270
# gca_mod_int_pstm_2 39 32837 33098 -16380    32759 1.4470      1     0.2290
# gca_mod_int_pstm_3 40 32839 33107 -16380    32759 0.1740      1     0.6766


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
  gca_mod_full_wm_0 <- update(gca_full_mod_base,   . ~ . + wm_std:group)
  gca_mod_full_wm_1 <- update(gca_mod_full_wm_0,   . ~ . + ot1:wm_std:group)
  gca_mod_full_wm_2 <- update(gca_mod_full_wm_1,   . ~ . + ot2:wm_std:group)
  gca_mod_full_wm_3 <- update(gca_mod_full_wm_2,   . ~ . + ot3:wm_std:group)

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

