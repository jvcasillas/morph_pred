# Load data -------------------------------------------------------------------

source(here::here("scripts", "01_load_data.R"))

# -----------------------------------------------------------------------------


# read data
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







## @knitr stressGroupCompare


# Load wm data and combine with stress_subset_0 (proportion data)
# in order to add working memory as a covariate
wm_df <- read_csv(here("data", "raw", "wm_all.csv")) %>%
         filter(., !(group %in% c("HS", "L")),
                   !(participant %in% c('L01', 'L02', 'L03', 'L04', 'L05',
                                        'L06', 'L07', 'L08', 'L09', 'L10',
                                        'L15', 'L20', 'L21', 'L22', 'L23',
                                        'L26', 'L30', 'L31', 'L33', 'La04',
                                        'LA06', 'LA07', 'LA14')))

scale_this <- function(x) as.vector(scale(x))








# random effects building

prop_0_ranefA <- glmer(cbind(targetCount, 10 - targetCount) ~ 1 +
                    (1 | participant),
                    data = df_stress, REML = F, family = 'binomial',
                    control = glmerControl(optimizer = 'bobyqa'))

prop_0_ranefB <- glmer(cbind(targetCount, 10 - targetCount) ~ 1 +
                    (1 | participant) + (1 | target),
                    data = df_stress, family = 'binomial',
                    control = glmerControl(optimizer = 'bobyqa'))

anova(prop_0_ranefA, prop_0_ranefB, refit = F) # keep intercept for target

prop_0_ranefC <- glmer(cbind(targetCount, 10 - targetCount) ~ 1 +
                    (1 | participant) + (1 | target) +
                    (0 + condition_sum | participant),
                    data = df_stress, family = 'binomial',
                    control = glmerControl(optimizer = 'bobyqa'))

anova(prop_0_ranefB, prop_0_ranefC, refit = F) # keep uncorrelated slope for condition

prop_0_ranefD <- glmer(cbind(targetCount, 10 - targetCount) ~ 1 +
                    (1 | participant) + (1 | target) +
                    (0 + condition_sum | participant) + (0 + coda_sum | participant),
                    data = df_stress, family = 'binomial',
                    control = glmerControl(optimizer = 'bobyqa'))

anova(prop_0_ranefC, prop_0_ranefD, refit = F) # Keep uncorrelated slope for coda

prop_0_ranefE <- glmer(cbind(targetCount, 10 - targetCount) ~ 1 +
                    (1 | participant) + (1 | target) +
                    (0 + condition_sum * coda_sum | participant),
                    data = df_stress, family = 'binomial',
                    control = glmerControl(optimizer = 'bobyqa'))

anova(prop_0_ranefD, prop_0_ranefE, refit = F) # Keep condition x coda interaction slope

prop_0_ranefF <- glmer(cbind(targetCount, 10 - targetCount) ~ 1 +
                    (1 | participant) + (1 | target) +
                    (1 + condition_sum * coda_sum | participant),
                    data = df_stress, family = 'binomial',
                    control = glmerControl(optimizer = 'bobyqa'))

anova(prop_0_ranefE, prop_0_ranefF, refit = F) # Keep correlated slope



# test fixed effects

prop_0_mod_0 <- glmer(cbind(targetCount, 10 - targetCount) ~ 1 +
                    (1 | participant) + (1 | target) +
                    (1 + condition_sum * coda_sum | participant),
                    data = df_stress, family = 'binomial',
                    control = glmerControl(optimizer = 'bobyqa'))

prop_0_mod_group <- glmer(cbind(targetCount, 10 - targetCount) ~ 1 + group +
                    (1 | participant) + (1 | target) +
                    (1 + condition_sum * coda_sum | participant),
                    data = df_stress, family = 'binomial',
                    control = glmerControl(optimizer = 'bobyqa'))

prop_0_mod_cond <- glmer(cbind(targetCount, 10 - targetCount) ~ 1 + group + condition_sum +
                    (1 | participant) + (1 | target) +
                    (1 + condition_sum * coda_sum | participant),
                    data = df_stress, family = 'binomial',
                    control = glmerControl(optimizer = 'bobyqa'))

prop_0_mod_coda <- glmer(cbind(targetCount, 10 - targetCount) ~ 1 + group + coda_sum +
                    (1 | participant) + (1 | target) +
                    (1 + condition_sum * coda_sum | participant),
                    data = df_stress, family = 'binomial',
                    control = glmerControl(optimizer = 'bobyqa'))

prop_0_mod_int1 <- glmer(cbind(targetCount, 10 - targetCount) ~ 1 + group + coda_sum + group:coda_sum +
                    (1 | participant) + (1 | target) +
                    (1 + condition_sum * coda_sum | participant),
                    data = df_stress, family = 'binomial',
                    control = glmerControl(optimizer = 'bobyqa'))

prop_0_mod_int2 <- glmer(cbind(targetCount, 10 - targetCount) ~ 1 + group + coda_sum + group:condition_sum +
                    (1 | participant) + (1 | target) +
                    (1 + condition_sum * coda_sum | participant),
                    data = df_stress, family = 'binomial',
                    control = glmerControl(optimizer = 'bobyqa'))

prop_0_mod_int3 <- glmer(cbind(targetCount, 10 - targetCount) ~ 1 + group + coda_sum + condition_sum:coda_sum +
                    (1 | participant) + (1 | target) +
                    (1 + condition_sum * coda_sum | participant),
                    data = df_stress, family = 'binomial',
                    control = glmerControl(optimizer = 'bobyqa'))

prop_0_mod_full <- glmer(cbind(targetCount, 10 - targetCount) ~ 1 + group * condition_sum * coda_sum +
                    (1 | participant) + (1 | target) +
                    (1 + condition_sum * coda_sum | participant),
                    data = df_stress, family = 'binomial',
                    control = glmerControl(optimizer = 'bobyqa'))



anova(prop_0_mod_0, prop_0_mod_group) # main effect of group
anova(prop_0_mod_group, prop_0_mod_cond) # no effect of condition
anova(prop_0_mod_group, prop_0_mod_coda) # main effect of coda
anova(prop_0_mod_coda, prop_0_mod_int1) # no group x coda interaction
anova(prop_0_mod_coda, prop_0_mod_int2) # no group condition interaction
anova(prop_0_mod_coda, prop_0_mod_int3) # no condi x coda interaction
anova(prop_0_mod_coda, prop_0_mod_full) # no three way interaction


#    AIC    BIC  logLik deviance  Chisq Chi Df Pr(>Chisq)
# 8419.2 8490.6 -4194.6   8389.2 11.273      2   0.003565  prop_0_mod_group ***
# 8420.6 8496.8 -4194.3   8388.6 0.6178      1     0.4319  prop_0_mod_cond
# 8413.5 8489.7 -4190.7   8381.5 7.7467      1   0.005381  prop_0_mod_coda  **
# 8415.3 8501.0 -4189.7   8379.3 2.1485      2     0.3416  group x coda
# 8412.6 8503.1 -4187.3   8374.6 6.897       3    0.07525  group x condition
# 8412.2 8497.9 -4188.1   8376.2 5.2914      2    0.07096  cond x coda
# 8420.5 8534.8 -4186.3   8372.5 8.9397      8     0.3474  prop_0_mod_full


df_stress$group <- factor(df_stress$group, levels = c("ss", "la",  "int"))

prop_0_mod_final <- glmer(cbind(targetCount, 10 - targetCount) ~ 1 + group + coda_sum +
                    (1 | participant) + (1 | target) +
                    (1 + condition_sum * coda_sum | participant),
                    data = df_stress, family = 'binomial',
                    control = glmerControl(optimizer = 'bobyqa'))

MuMIn::r.squaredGLMM(prop_0_mod_final)

# summary(prop_0_mod_final)
# confint(prop_0_mod_final, method = "Wald")

# Fixed effects:
#               Estimate Std. Error z value Pr(>|z|)
# (Intercept)     1.6025     0.2886   5.553 2.81e-08 ***
# groupla        -1.1382     0.3274  -3.476 0.000509 ***
# groupint       -0.4926     0.4270  -1.153 0.248711
# codaSum1       -0.8177     0.2258  -3.621 0.000294 ***


# Relevel to test lb vs la
df_stress$group <- factor(df_stress$group, levels = c("int", "la",  "ss"))

summary(glmer(cbind(targetCount, 10 - targetCount) ~ 1 + group + coda_sum +
       (1 | participant) + (1 | target) +
       (1 + condition_sum * coda_sum | participant),
       data = df_stress, family = 'binomial',
       control = glmerControl(optimizer = 'bobyqa')))

# Fixed effects:
#             Estimate Std. Error z value Pr(>|z|)
# (Intercept)   1.0620     0.3840   2.765  0.00568 **
# groupla      -0.6426     0.4049  -1.587  0.11249
# groupss       0.4879     0.4244   1.150  0.25028
# codaSum1     -0.6111     0.2170  -2.816  0.00486 **
