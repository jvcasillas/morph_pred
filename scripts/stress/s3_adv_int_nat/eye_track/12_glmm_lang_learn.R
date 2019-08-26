# GLMMs -----------------------------------------------------------------------
#
# - Target fixation as a function of group, stress (condition), and coda
#   at the offset of the first syllable (time_zero == 20)
# - This model builds on the t-test analyses by looking for between group
#   differences in target fixation at the offfset of first syllable
#
# -----------------------------------------------------------------------------




# Load data and models --------------------------------------------------------

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

stress10 <- stress10 %>%
  filter(., group %in% c("la", "int", "ss")) %>%
  left_join(mem_data, by = "participant") %>%
  left_join(phon_data, by = "target") %>%
  left_join(freq_data, by = "target") %>%
  rename(., lex_prob = "freq_nim")

stress10$wm <- as.numeric(stress10$wm)

glimpse(stress10)



prop_0_mod_0     <- readRDS(here("models", "stress", "s3_adv_int_nat",
                                 "eye_track", "glmm", "0_prop_0_mod_0.rds"))
prop_0_mod_group <- readRDS(here("models", "stress", "s3_adv_int_nat",
                                 "eye_track", "glmm", "1_prop_0_mod_group.rds"))
prop_0_mod_cond  <- readRDS(here("models", "stress", "s3_adv_int_nat",
                                 "eye_track", "glmm", "2_prop_0_mod_cond.rds"))
prop_0_mod_coda  <- readRDS(here("models", "stress", "s3_adv_int_nat",
                                 "eye_track", "glmm", "3_prop_0_mod_coda.rds"))
prop_0_mod_int1  <- readRDS(here("models", "stress", "s3_adv_int_nat",
                                 "eye_track", "glmm", "4_prop_0_mod_int1.rds"))
prop_0_mod_int2  <- readRDS(here("models", "stress", "s3_adv_int_nat",
                                 "eye_track", "glmm", "5_prop_0_mod_int2.rds"))
prop_0_mod_int3  <- readRDS(here("models", "stress", "s3_adv_int_nat",
                                 "eye_track", "glmm", "6_prop_0_mod_int3.rds"))
prop_0_mod_full  <- readRDS(here("models", "stress", "s3_adv_int_nat",
                                 "eye_track", "glmm", "7_prop_0_mod_full.rds"))
prop_0_mod_final <- readRDS(here("models", "stress", "s3_adv_int_nat",
                                 "eye_track", "glmm", "8_prop_0_mod_final.rds"))
# -----------------------------------------------------------------------------









# Data prep -------------------------------------------------------------------

# Get subset of int, la, and ss groups
# Filter time course to offset of 1st syllable (time_zero == 20)
# Create sum coded fixed factors (condition and coda)
# Create standarized fixed factors (lex freq, phon freq & wm)
df_stress <- stress10 %>%
  filter(., group %in% c('int', 'la', 'ss'),
         !participant == "LA07", # missing WM
         time_zero == 20) %>%
  mutate(., condition_sum = if_else(condition == "stressed", 1, -1),
         coda_sum = if_else(coda == 1, 1, -1),
         lex_std = (lex_prob - mean (lex_prob)) / sd(lex_prob),
         wm_std = (wm - mean(wm, na.rm = T)) / sd(wm, na.rm = T),
         phon_std = (phon_prob - mean(phon_prob)) / sd(phon_prob))

# -----------------------------------------------------------------------------

View(df_stress)





# Random effects building -----------------------------------------------------

if(F) {

  ranefA <- glmer(cbind(targetCount, 10 - targetCount) ~ 1 +
                           (1 | participant),
                         data = df_stress, family = 'binomial',
                         control = glmerControl(optimizer = 'bobyqa'))

  ranefB <- glmer(cbind(targetCount, 10 - targetCount) ~ 1 +
                           (1 | participant) +
                           (1 | target),
                         data = df_stress, family = 'binomial',
                         control = glmerControl(optimizer = 'bobyqa'))

  anova(ranefA, ranefB, refit = F) # keep intercept for target

  ranefC <- glmer(cbind(targetCount, 10 - targetCount) ~ 1 +
                           (1 + condition_sum | participant) +
                           (1 | target),
                         data = df_stress, family = 'binomial',
                         control = glmerControl(optimizer = 'bobyqa'))

  anova(ranefB, ranefC, refit = F) # keep slope for condition

  ranefD <- glmer(cbind(targetCount, 10 - targetCount) ~ 1 +
                           (1 + condition_sum + coda_sum | participant) +
                           (1 | target),
                         data = df_stress, family = 'binomial',
                         control = glmerControl(optimizer = 'bobyqa'))

  anova(ranefC, ranefD, refit = F) # Keep slope for coda

  ranefE <- glmer(cbind(targetCount, 10 - targetCount) ~ 1 +
                           (1 + condition_sum * coda_sum | participant) +
                           (1 | target),
                         data = df_stress, family = 'binomial',
                         control = glmerControl(optimizer = 'bobyqa'))

  anova(ranefD, ranefE, refit = F) # Keep interaction slope
}

# -----------------------------------------------------------------------------







# Test fixed effects ----------------------------------------------------------

if(F) {
  mod_0 <- glmer(cbind(targetCount, 10 - targetCount) ~ 1 +
                          (1 + condition_sum * coda_sum | participant) +
                          (1 | target),
                        data = df_stress, family = 'binomial',
                        control = glmerControl(optimizer = 'bobyqa'))

  mod_group <- update(mod_0,     . ~ . + group)
  mod_cond  <- update(mod_group, . ~ . + condition_sum)
  mod_coda  <- update(mod_cond,  . ~ . + coda_sum)
  mod_wm    <- update(mod_coda,  . ~ . + wm_std)
  mod_phon  <- update(mod_wm,  . ~ . + phon_std)
  mod_lex  <- update(mod_phon,  . ~ . + lex_std)
  # mod_int1  <- update(mod_coda,  . ~ . + group:coda_sum)
  # mod_int2  <- update(mod_int1,  . ~ . + group:condition_sum)
  # mod_int3  <- update(mod_int2,  . ~ . + coda_sum:condition_sum)
  # mod_full  <- update(mod_int3,  . ~ . + group:coda_sum:condition_sum)


  anova(mod_0, mod_group, test = "Chisq")    # main effect of group
  anova(mod_group, mod_cond, test = "Chisq") # no effect of condition
  anova(mod_cond, mod_coda, test = "Chisq") # main effect of coda
  anova(mod_coda, mod_phon, test = "Chisq") # no effect of phonotactic frequency
  anova(mod_group, mod_lex, test = "Chisq") # no effect of lexical frequency
  anova(mod_coda, mod_wm, test = "Chisq") #
  # anova(mod_coda, mod_int1, test = "Chisq")  # no group x coda interaction
  # anova(mod_coda, mod_int2, test = "Chisq")  # no group condition interaction
  # anova(mod_coda, mod_int3, test = "Chisq")  # no condi x coda interaction
  # anova(mod_coda, mod_full, test = "Chisq")  # no three way interaction


  #    AIC    BIC  logLik deviance  Chisq Chi Df Pr(>Chisq)
  #  11876  11946 -5923.7    11848 12.077      2   0.002386  prop_0_mod_group **
  #  11876  11951 -5923.0    11846 1.4415      1     0.2299  prop_0_mod_cond
  #  11873  11954 -5920.6    11841 6.3035      2    0.04278  prop_0_mod_coda  *
  #  11876  11966 -5919.9    11840 1.4065      2      0.495  group x coda
  #  11879  11980 -5919.5    11839 2.1075      4      0.716  group x condition
  #  11881  11986 -5919.3    11839 2.5777      5     0.7648  cond x coda
  #  11884  12000 -5919.0    11838 3.11        7     0.8746  prop_0_mod_full


  df_stress$group <- factor(df_stress$group, levels = c("ss", "la",  "int"))

  mod_final <- glmer(cbind(targetCount, 10 - targetCount) ~ 1 +
                              group + coda_sum +
                              (1 + condition_sum * coda_sum | participant) +
                              (1 | target),
                            data = df_stress, family = 'binomial',
                            control = glmerControl(optimizer = 'bobyqa'))

}

# -----------------------------------------------------------------------------








# Save models -----------------------------------------------------------------

if(F) {

  saveRDS(mod_0, here("models", "stress", "s3_adv_int_nat",
                             "eye_track", "glmm", "0_mod_0.rds"))
  saveRDS(mod_group, here("models", "stress", "s3_adv_int_nat",
                                 "eye_track", "glmm", "1_mod_group.rds"))
  saveRDS(mod_cond, here("models", "stress", "s3_adv_int_nat",
                                "eye_track", "glmm", "2_mod_cond.rds"))
  saveRDS(mod_coda, here("models", "stress", "s3_adv_int_nat",
                                "eye_track", "glmm", "3_mod_coda.rds"))
  saveRDS(mod_int1, here("models", "stress", "s3_adv_int_nat",
                                "eye_track", "glmm", "4_mod_int1.rds"))
  saveRDS(mod_int2, here("models", "stress", "s3_adv_int_nat",
                                "eye_track", "glmm", "5_mod_int2.rds"))
  saveRDS(mod_int3, here("models", "stress", "s3_adv_int_nat",
                                "eye_track", "glmm", "6_mod_int3.rds"))
  saveRDS(mod_full, here("models", "stress", "s3_adv_int_nat",
                                "eye_track", "glmm", "7_mod_full.rds"))
  saveRDS(mod_final, here("models", "stress", "s3_adv_int_nat",
                                 "eye_track", "glmm", "8_mod_final.rds"))
}

# -----------------------------------------------------------------------------









# Model descriptives ----------------------------------------------------------

MuMIn::r.squaredGLMM(mod_final)
#              R2m  R2c
# theoretical 0.04 0.57
# delta       0.03 0.52

# summary(prop_0_mod_final)
# confint(prop_0_mod_final, method = "Wald")

# Fixed effects:
#               Estimate Std. Error CI-low  CI-high z value Pr(>|z|)
#   (Intercept)   1.1638     0.2469   0.68     1.65   4.714 2.43e-06 ***
#   groupla      -0.8865     0.3021  -1.48    -0.29  -2.934  0.00334 **
#   groupint     -1.0232     0.3120  -1.63    -0.41  -3.279  0.00104 **
#   coda_sum      0.2580     0.1545  -0.04     0.56   1.670  0.09485 .


# Relevel to test int vs la
df_stress$group <- factor(df_stress$group, levels = c("int", "la",  "ss"))

summary(glmer(cbind(targetCount, 10 - targetCount) ~ 1 +
                group + coda_sum +
                (1 + condition_sum * coda_sum | participant) +
                (1 | target),
              data = df_stress, family = 'binomial',
              control = glmerControl(optimizer = 'bobyqa')))

# Fixed effects:
#             Estimate Std. Error z value Pr(>|z|)
# (Intercept)   0.1405     0.2525   0.557  0.57778
# groupla       0.1368     0.3030   0.451  0.65176
# groupss       1.0232     0.3120   3.279  0.00104 **
# coda_sum      0.2580     0.1545   1.670  0.09486 .

# -----------------------------------------------------------------------------
