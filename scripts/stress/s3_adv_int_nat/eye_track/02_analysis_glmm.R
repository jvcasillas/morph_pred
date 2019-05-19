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
            !participant %in% c("L01", "L02", "L03", "L04", "L05",
                                "L06", "L07", "L08", "L09", "L10",
                                "L15", "L20", "L21", "L22", "L23",
                                "L26", "L30", "L31", "L33", "LA04",
                                "LA06", "LA09", "LA14", "LA15", "LA19"),
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
anova(prop_0_mod_group, prop_0_mod_coda) # main effect of coda
anova(prop_0_mod_coda, prop_0_mod_int1)  # no group x coda interaction
anova(prop_0_mod_coda, prop_0_mod_int2)  # no group condition interaction
anova(prop_0_mod_coda, prop_0_mod_int3)  # no condi x coda interaction
anova(prop_0_mod_coda, prop_0_mod_full)  # no three way interaction


#    AIC    BIC  logLik deviance  Chisq Chi Df Pr(>Chisq)
#  11876  11946 -5923.7    11848 12.077      2   0.002386  prop_0_mod_group **
#  11876  11951 -5923.0    11846 1.4415      1     0.2299  prop_0_mod_cond
#  11873  11954 -5920.6    11841 6.3035      2    0.04278  prop_0_mod_coda  .
#  11876  11966 -5919.9    11840 1.4065      2      0.495  group x coda
#  11879  11980 -5919.5    11839 2.1075      4      0.716  group x condition
#  11881  11986 -5919.3    11839 2.5777      5     0.7648  cond x coda
#  11884  12000 -5919.0    11838 3.11        7     0.8746  prop_0_mod_full


df_stress$group <- factor(df_stress$group, levels = c("ss", "la",  "int"))

prop_0_mod_final <- glmer(cbind(targetCount, 10 - targetCount) ~ 1 +
                    group + coda_sum +
                    (1 + condition_sum * coda_sum | participant) +
                    (1 | target),
                    data = df_stress, family = 'binomial',
                    control = glmerControl(optimizer = 'bobyqa'))

MuMIn::r.squaredGLMM(prop_0_mod_final)
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
# (Intercept)   0.1405     0.2525   0.557  0.57778
# groupla       0.1368     0.3030   0.451  0.65176
# groupss       1.0232     0.3120   3.279  0.00104 **
# coda_sum      0.2580     0.1545   1.670  0.09486 .

# -----------------------------------------------------------------------------
