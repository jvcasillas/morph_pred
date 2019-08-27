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
                                 "eye_track", "glmm", "0_mod_0.rds"))
mod_group <- readRDS(here("models", "stress", "s3_adv_int_nat",
                                 "eye_track", "glmm", "1_mod_group.rds"))
mod_cond  <- readRDS(here("models", "stress", "s3_adv_int_nat",
                                 "eye_track", "glmm", "2_mod_cond.rds"))
mod_coda  <- readRDS(here("models", "stress", "s3_adv_int_nat",
                                 "eye_track", "glmm", "3_mod_coda.rds"))
mod_int1  <- readRDS(here("models", "stress", "s3_adv_int_nat",
                                 "eye_track", "glmm", "4_mod_int1.rds"))
mod_int2  <- readRDS(here("models", "stress", "s3_adv_int_nat",
                                 "eye_track", "glmm", "5_mod_int2.rds"))
mod_int3  <- readRDS(here("models", "stress", "s3_adv_int_nat",
                                 "eye_track", "glmm", "6_mod_int3.rds"))
mod_full  <- readRDS(here("models", "stress", "s3_adv_int_nat",
                                 "eye_track", "glmm", "7_mod_full.rds"))
mod_final <- readRDS(here("models", "stress", "s3_adv_int_nat",
                                 "eye_track", "glmm", "8_mod_final.rds"))
# -----------------------------------------------------------------------------








# Data prep -------------------------------------------------------------------

# Get subset of int, la, and ss groups
# Filter time course to offset of 1st syllable (time_zero == 20)
# Create sum coded fixed factors (condition and coda)
# Create standarized fixed factors (lex freq, phon freq & wm)
df_stress <- stress10 %>%
  filter(., group %in% c('int', 'la', 'ss'),
         !participant %in% c("LA04", "LA06", "LA07", "LA14", "LA19"), # missing WM
         time_zero == 20) %>%
  mutate(., condition_sum = if_else(condition == "stressed", 1, -1),
         coda_sum = if_else(coda == 1, 1, -1),
         lex_std = (lex_prob - mean (lex_prob)) / sd(lex_prob),
         wm_std = (wm - mean(wm, na.rm = T)) / sd(wm, na.rm = T),
         phon_std = (phon_prob - mean(phon_prob)) / sd(phon_prob))

# -----------------------------------------------------------------------------






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







# Test fixed effects  with WM ----------------------------------------------------------

if(F) {
  mod_0 <- glmer(cbind(targetCount, 10 - targetCount) ~ 1 +
                          (1 + condition_sum * coda_sum | participant) +
                          (1 | target),
                        data = df_stress, family = 'binomial',
                        control = glmerControl(optimizer = 'bobyqa'))

  mod_group     <- update(mod_0,     . ~ . + group)
  mod_cond      <- update(mod_group, . ~ . + condition_sum)
  mod_coda      <- update(mod_cond,  . ~ . + coda_sum)
  mod_wm        <- update(mod_coda,  . ~ . + wm_std)
  mod_wm_group  <- update(mod_coda, . ~ . + group:wm_std)
  mod_wm_coda   <- update(mod_coda, . ~ . + coda_sum:wm_std)
  mod_wm_cond   <- update(mod_coda, . ~ . + condition_sum:wm_std)
  mod_wm_triple <- update(mod_coda, . ~ . + condition_sum:coda_sum:wm_std)
  mod_wm_quadr  <- update(mod_coda, . ~ . + condition_sum:coda_sum:wm_std:group)



  anova(mod_0, mod_group, test = "Chisq")    # main effect of group
  anova(mod_group, mod_cond, test = "Chisq") # no effect of condition
  anova(mod_cond, mod_coda, test = "Chisq") # main effect of coda
  anova(mod_coda, mod_wm, test = "Chisq") # no main effect of wm
  anova(mod_coda, mod_wm_group, test = "Chisq")  # no main effect of wm:group
  anova(mod_coda, mod_wm_coda, test = "Chisq") # no main effect of wm:coda
  anova(mod_coda, mod_wm_cond, test = "Chisq") # no main effect of wm:cond
  anova(mod_coda, mod_wm_triple, test = "Chisq") # no main effect of wm:cond:coda
  anova(mod_coda, mod_wm_quadr, test = "Chisq") # no main effect of wm:cond:coda:group



  # df_stress$group <- factor(df_stress$group, levels = c("ss", "la",  "int"))
  #
  # mod_final <- glmer(cbind(targetCount, 10 - targetCount) ~ 1 +
  #                             group + coda_sum +
  #                             (1 + condition_sum * coda_sum | participant) +
  #                             (1 | target),
  #                           data = df_stress, family = 'binomial',
  #                           control = glmerControl(optimizer = 'bobyqa'))

}

# -----------------------------------------------------------------------------

# Test fixed effects  with lexical frequency ----------------------------------------------------------

if(F) {
  mod_0 <- glmer(cbind(targetCount, 10 - targetCount) ~ 1 +
                   (1 + condition_sum * coda_sum | participant) +
                   (1 | target),
                 data = df_stress, family = 'binomial',
                 control = glmerControl(optimizer = 'bobyqa'))

  mod_group     <- update(mod_0,     . ~ . + group)
  mod_cond      <- update(mod_group, . ~ . + condition_sum)
  mod_coda      <- update(mod_cond,  . ~ . + coda_sum)
  mod_lex        <- update(mod_coda,  . ~ . + lex_std)
  mod_lex_group  <- update(mod_coda, . ~ . + group:lex_std)
  mod_lex_coda   <- update(mod_lex_group, . ~ . + coda_sum:lex_std)
  mod_lex_cond   <- update(mod_lex_coda, . ~ . + condition_sum:lex_std)
  mod_lex_triple <- update(mod_lex_cond, . ~ . + condition_sum:coda_sum:lex_std)
  mod_lex_quadr  <- update(mod_lex_triple, . ~ . + condition_sum:coda_sum:lex_std:group)



  anova(mod_0, mod_group, test = "Chisq")    # main effect of group
  anova(mod_group, mod_cond, test = "Chisq") # no effect of condition
  anova(mod_cond, mod_coda, test = "Chisq") # main effect of coda
  anova(mod_coda, mod_lex, test = "Chisq") # no main effect of lex
  anova(mod_lex, mod_lex_group, test = "Chisq")  # main effect of lex:group
  anova(mod_lex_group, mod_lex_coda, test = "Chisq") # no main effect of lex:coda
  anova(mod_lex_group, mod_lex_cond, test = "Chisq") # no main effect of lex:cond
  anova(mod_lex_group, mod_lex_triple, test = "Chisq") # no main effect of lex:cond:coda
  anova(mod_lex_group, mod_lex_quadr, test = "Chisq") # no main effect of lex:cond:coda:group



  df_stress$group <- factor(df_stress$group, levels = c("ss", "la", "int"))

  mod_final_lex <- glmer(cbind(targetCount, 10 - targetCount) ~ 1 +
                              group + coda_sum + lex_std + lex_std:group +
                              (1 + condition_sum * coda_sum | participant) +
                              (1 | target),
                            data = df_stress, family = 'binomial',
                            control = glmerControl(optimizer = 'bobyqa'))
summary(mod_final_lex)

# Fixed effects:
#                  Estimate Std. Error z value Pr(>|z|)
# (Intercept)       1.14538    0.24680   4.641 3.47e-06 ***
# groupla          -0.92824    0.29938  -3.101 0.001931 **
# groupint         -0.98170    0.31261  -3.140 0.001688 **
# coda_sum          0.24017    0.15206   1.579 0.114245
# lex_std           0.15723    0.12952   1.214 0.224772
# groupla:lex_std  -0.31372    0.07309  -4.292 1.77e-05 ***
# groupint:lex_std -0.25168    0.07379  -3.411 0.000647 ***

df_stress %<>% mutate(., group = fct_relevel(group, "int"))
mod_full_mod_int_relevel <- update(mod_final_lex)

summary(mod_full_mod_int_relevel)


# Fixed effects:
#                 Estimate Std. Error z value Pr(>|z|)
# (Intercept)      0.16368    0.25177   0.650 0.515612
# groupss          0.98169    0.31262   3.140 0.001688 **
# groupla          0.05345    0.30142   0.177 0.859246
# coda_sum         0.24016    0.15206   1.579 0.114248
# lex_std         -0.09445    0.13091  -0.721 0.470621
# groupss:lex_std  0.25168    0.07379   3.411 0.000647 ***
# groupla:lex_std -0.06204    0.07196  -0.862 0.388622

}

# -----------------------------------------------------------------------------

# Test fixed effects  with phonotactic frequency ----------------------------------------------------------

df_stress$group <- factor(df_stress$group, levels = c("ss", "la", "int"))

if(F) {
  mod_0 <- glmer(cbind(targetCount, 10 - targetCount) ~ 1 +
                   (1 + condition_sum * coda_sum | participant) +
                   (1 | target),
                 data = df_stress, family = 'binomial',
                 control = glmerControl(optimizer = 'bobyqa'))

  mod_group     <- update(mod_0,     . ~ . + group)
  mod_cond      <- update(mod_group, . ~ . + condition_sum)
  mod_coda      <- update(mod_cond,  . ~ . + coda_sum)
  mod_phon        <- update(mod_coda,  . ~ . + phon_std)
  mod_phon_group  <- update(mod_coda, . ~ . + group:phon_std)
  mod_phon_coda   <- update(mod_phon_group, . ~ . + coda_sum:phon_std)
  mod_phon_cond   <- update(mod_phon_group, . ~ . + condition_sum:phon_std)
  mod_phon_triple <- update(mod_phon_group, . ~ . + condition_sum:coda_sum:phon_std)
  mod_phon_quadr  <- update(mod_phon_group, . ~ . + condition_sum:coda_sum:phon_std:group)



  anova(mod_0, mod_group, test = "Chisq")    # main effect of group
  anova(mod_group, mod_cond, test = "Chisq") # no effect of condition
  anova(mod_cond, mod_coda, test = "Chisq") # main effect of coda
  anova(mod_coda, mod_phon, test = "Chisq") # no main effect of phon
  anova(mod_coda, mod_phon_group, test = "Chisq")  # main effect of phon:group
  anova(mod_phon_group, mod_phon_coda, test = "Chisq") # no main effect of phon:coda
  anova(mod_phon_group, mod_phon_cond, test = "Chisq") # no main effect of phon:cond
  anova(mod_phon_group, mod_phon_triple, test = "Chisq") # no main effect of phon:cond:coda
  anova(mod_phon_group, mod_phon_quadr, test = "Chisq") # main effect of phon:cond:coda:group



  df_stress$group <- factor(df_stress$group, levels = c("ss", "la", "int"))

  mod_final_phon <- glmer(cbind(targetCount, 10 - targetCount) ~ 1 +
                              group + coda_sum + phon_std + group:phon_std +
                               condition_sum:coda_sum:phon_std:group +
                              (1 + condition_sum * coda_sum | participant) +
                              (1 | target),
                            data = df_stress, family = 'binomial',
                            control = glmerControl(optimizer = 'bobyqa'))

  summary(mod_final_phon)

  # Fixed effects:
  #                                           Estimate Std. Error z value Pr(>|z|)
  # (Intercept)                               1.19617    0.24516   4.879 1.07e-06 ***
  # groupla                                  -1.02290    0.29960  -3.414 0.000640 ***
  # groupint                                 -1.09490    0.31091  -3.522 0.000429 ***    coda_sum                                  0.24298    0.15505   1.567 0.117097
  # phon_std                                  0.24965    0.12763   1.956 0.050466 .
  # groupla:phon_std                         -0.15007    0.06470  -2.319 0.020374 *
  # groupint:phon_std                        -0.41318    0.06460  -6.396 1.60e-10 ***
  # groupss:coda_sum:phon_std:condition_sum  -0.06747    0.12453  -0.542 0.587930
  # groupla:coda_sum:phon_std:condition_sum  -0.33358    0.12406  -2.689 0.007169 **
  # groupint:coda_sum:phon_std:condition_sum -0.14202    0.12387  -1.147 0.251575



  df_stress %<>% mutate(., group = fct_relevel(group, "int"))
  mod_full_phon_int_relevel <- update(mod_final_phon)

  summary(mod_full_phon_int_relevel)


  # Fixed effects:
  #                                           Estimate Std. Error z value Pr(>|z|)
  # (Intercept)                               0.10128    0.25035   0.405 0.685797
  # groupss                                   1.09489    0.31088   3.522 0.000428 ***
  # groupla                                   0.07198    0.30181   0.239 0.811483
  # coda_sum                                  0.24298    0.15505   1.567 0.117093
  # phon_std                                 -0.16353    0.12698  -1.288 0.197793
  # groupss:phon_std                          0.41318    0.06460   6.396 1.60e-10 ***
  # groupla:phon_std                          0.26311    0.06317   4.165 3.11e-05 ***
  # groupint:coda_sum:phon_std:condition_sum -0.14202    0.12386  -1.147 0.251567
  # groupss:coda_sum:phon_std:condition_sum  -0.06747    0.12452  -0.542 0.587925
  # groupla:coda_sum:phon_std:condition_sum  -0.33358    0.12406  -2.689 0.007169 **

}

# -----------------------------------------------------------------------------







# Save models -----------------------------------------------------------------

if(F) {
  lang_learn_glmm_mods <- mget(c("mod_0", "mod_group", "mod_cond", "mod_coda",
                                 "mod_wm", "mod_wm_group", "mod_wm_coda", "mod_wm_cond",
                                 "mod_wm_triple", "mod_wm_quadr", "mod_lex", "mod_lex_group", "mod_lex_coda", "mod_lex_cond",
                                 "mod_lex_triple", "mod_lex_quadr", "mod_phon", "mod_phon_group", "mod_phon_coda", "mod_phon_cond",
                                 "mod_phon_triple", "mod_phon_quadr"))

  save(lang_learn_glmm_mods,
        file = here("models", "stress", "s3_adv_int_nat", "eye_track", "gca",
                 "lang_learn_glmm_mods.Rdata"))


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
