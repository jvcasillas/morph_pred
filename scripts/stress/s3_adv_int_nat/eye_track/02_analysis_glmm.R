# GLMMs -----------------------------------------------------------------------
#
# - Target fixation as a function of group, stress (condition), and coda
#   at the offset of the first syllable (time_zero == 20)
# - This model builds on the t-test analyses by looking for between group
#   differences in target fixation at the offfset of first syllable
#
# -----------------------------------------------------------------------------




# Load data -------------------------------------------------------------------

source(here::here("scripts", "01_load_data.R"))

# -----------------------------------------------------------------------------









# Data prep -------------------------------------------------------------------

# Get subset of int, la, and ss groups
# Remove la participants (not sure why)
# Filter time course to offset of 1st syllable (time_zero == 20)
# Create sum coded fixed factors (condition and coda)
df_stress <- stress10 %>%
  filter(., group %in% c('int', 'la', 'ss'),
            !participant %in% c('L01', 'L02', 'L03', 'L04', 'L05',
                                'L06', 'L07', 'L08', 'L09', 'L10',
                                'L15', 'L20', 'L21', 'L22', 'L23',
                                'L26', 'L30', 'L31', 'L33', 'LA04',
                                'LA06', 'LA07', 'LA14'),
            time_zero == 20) %>%
  mutate(., condition_sum = if_else(condition == "stressed", 1, -1),
            coda_sum = if_else(coda == 1, 1, -1))

# -----------------------------------------------------------------------------








# Random effects building -----------------------------------------------------

prop_0_ranefA <- glmer(cbind(targetCount, 10 - targetCount) ~ 1 +
                    (1 | participant),
                    data = df_stress, REML = F, family = 'binomial',
                    control = glmerControl(optimizer = 'bobyqa'))

prop_0_ranefB <- glmer(cbind(targetCount, 10 - targetCount) ~ 1 +
                    (1 | participant) +
                    (1 | target),
                    data = df_stress, family = 'binomial',
                    control = glmerControl(optimizer = 'bobyqa'))

anova(prop_0_ranefA, prop_0_ranefB, refit = F) # keep intercept for target

prop_0_ranefC <- glmer(cbind(targetCount, 10 - targetCount) ~ 1 +
                    (1 + condition_sum | participant) +
                    (1 | target),
                    data = df_stress, family = 'binomial',
                    control = glmerControl(optimizer = 'bobyqa'))

anova(prop_0_ranefB, prop_0_ranefC, refit = F) # keep slope for condition

prop_0_ranefD <- glmer(cbind(targetCount, 10 - targetCount) ~ 1 +
                    (1 + condition_sum + coda_sum | participant) +
                    (1 | target),
                    data = df_stress, family = 'binomial',
                    control = glmerControl(optimizer = 'bobyqa'))

anova(prop_0_ranefC, prop_0_ranefD, refit = F) # Keep slope for coda

prop_0_ranefE <- glmer(cbind(targetCount, 10 - targetCount) ~ 1 +
                    (1 + condition_sum * coda_sum | participant) +
                    (1 | target),
                    data = df_stress, family = 'binomial',
                    control = glmerControl(optimizer = 'bobyqa'))

anova(prop_0_ranefD, prop_0_ranefE, refit = F) # Keep interaction slope

# -----------------------------------------------------------------------------






# Test fixed effects ----------------------------------------------------------

prop_0_mod_0 <- glmer(cbind(targetCount, 10 - targetCount) ~ 1 +
                    (1 + condition_sum * coda_sum | participant) +
                    (1 | target),
                    data = df_stress, family = 'binomial',
                    control = glmerControl(optimizer = 'bobyqa'))

prop_0_mod_group <- update(prop_0_mod_0,     . ~ . + group)
prop_0_mod_cond  <- update(prop_0_mod_group, . ~ . + condition_sum)
prop_0_mod_coda  <- update(prop_0_mod_cond,  . ~ . + coda_sum)
prop_0_mod_int1  <- update(prop_0_mod_coda,  . ~ . + group:coda_sum)
prop_0_mod_int2  <- update(prop_0_mod_int1,  . ~ . + group:condition_sum)
prop_0_mod_int3  <- update(prop_0_mod_int2,  . ~ . + coda_sum:condition_sum)
prop_0_mod_full  <- update(prop_0_mod_int3,  . ~ . + group:coda_sum:condition_sum)


anova(prop_0_mod_0, prop_0_mod_group)    # main effect of group
anova(prop_0_mod_group, prop_0_mod_cond) # no effect of condition
anova(prop_0_mod_group, prop_0_mod_coda) # marginal effect of coda
anova(prop_0_mod_coda, prop_0_mod_int1)  # no group x coda interaction
anova(prop_0_mod_coda, prop_0_mod_int2)  # no group condition interaction
anova(prop_0_mod_coda, prop_0_mod_int3)  # no condi x coda interaction
anova(prop_0_mod_coda, prop_0_mod_full)  # no three way interaction


#    AIC    BIC  logLik deviance  Chisq Chi Df Pr(>Chisq)
#  12132  12203 -6052.2    12104 11.452      2   0.003261  prop_0_mod_group **
#  12133  12209 -6051.6    12103 1.3599      1     0.2436  prop_0_mod_cond
#  12131  12212 -6049.6    12099 5.3169      2    0.07006  prop_0_mod_coda  .
#  12134  12224 -6048.8    12098 1.5352      2     0.4641  group x coda
#  12137  12238 -6048.6    12097 2.0497      4     0.7266  group x condition
#  12138  12244 -6048.2    12096 2.7347      5     0.7408  cond x coda
#  12142  12258 -6047.9    12096 3.3472      7     0.8511  prop_0_mod_full


df_stress$group <- factor(df_stress$group, levels = c("ss", "la",  "int"))

prop_0_mod_final <- glmer(cbind(targetCount, 10 - targetCount) ~ 1 +
                    group + coda_sum +
                    (1 + condition_sum * coda_sum | participant) +
                    (1 | target),
                    data = df_stress, family = 'binomial',
                    control = glmerControl(optimizer = 'bobyqa'))

MuMIn::r.squaredGLMM(prop_0_mod_final)
#              R2m  R2c
# theoretical 0.03 0.57
# delta       0.03 0.52

# summary(prop_0_mod_final)
# confint(prop_0_mod_final, method = "Wald")

# Fixed effects:
#               Estimate Std. Error z value Pr(>|z|)
#   (Intercept)   1.1633     0.2502   4.649 3.33e-06 ***
#   groupla      -0.8907     0.3070  -2.902  0.00371 **
#   groupint     -1.0214     0.3213  -3.179  0.00148 **
#   coda_sum      0.2210     0.1511   1.463  0.14351


# Relevel to test lb vs la
df_stress$group <- factor(df_stress$group, levels = c("int", "la",  "ss"))

summary(glmer(cbind(targetCount, 10 - targetCount) ~ 1 +
        group + coda_sum +
        (1 + condition_sum * coda_sum | participant) +
        (1 | target),
        data = df_stress, family = 'binomial',
        control = glmerControl(optimizer = 'bobyqa')))

# Fixed effects:
#             Estimate Std. Error z value Pr(>|z|)
# (Intercept)   0.1419     0.2552   0.556  0.57822
# groupla       0.1307     0.3063   0.427  0.66958
# groupss       1.0214     0.3212   3.180  0.00147 **
# coda_sum      0.2210     0.1511   1.463  0.14350

# -----------------------------------------------------------------------------
