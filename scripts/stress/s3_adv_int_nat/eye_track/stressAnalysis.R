#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Morphosyntactic predictability                                              #
# Stress analyses                                                             #
# 11/09/2016                                                                  #
# Script 1                                                                    #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# clean working directory
rm(list = ls(all = TRUE))

## @knitr stressLibs

library(tidyverse); library(broom); library(sjPlot)
library(lme4); library(lmerTest); library(gridExtra); library(cowplot)

## @knitr ignore



# read data
df_stress <- read_csv("./mySources/data/clean/stressBIN10iaClean.csv") %>%
  filter(., corr == 1,
            group %in% c('int', 'la', 'ss'),
            !participant %in% c('L01', 'L02', 'L03', 'L04', 'L05',
                                 'L06', 'L07', 'L08', 'L09', 'L10',
                                 'L15', 'L20', 'L21', 'L22', 'L23',
                                 'L26', 'L30', 'L31', 'L33', 'LA04',
                                 'LA06', 'LA07', 'LA14'))

glimpse(df_stress)
unique(df_stress$group)
unique(df_stress$condition)

predict.glmm <- function(fit, newdata, type = 'gaussian') {

  ## Extract DV:

  DV = strsplit(as.character(fit@call), ' ~ ')[[2]][1]

  ## Append empty DV column:

  newdata = cbind(newdata, rep(0, nrow(newdata)))
  colnames(newdata)[ncol(newdata)] = DV

  ## Extract model matrix:

    mm <- model.matrix(terms(fit), newdata)

    ## Add fitted vals:

    newdata[, DV] <- predict(fit, newdata, re.form = NA)

    ## Extract confidence intervals:

    pvar1 <- diag(mm %*% tcrossprod(vcov(fit), mm))
    newdata$UB <- newdata[, DV] + 1.96 * sqrt(pvar1)
    newdata$LB <- newdata[, DV] - 1.96 * sqrt(pvar1)

    ## Transform if working with a non-gaussian GLMM:

    if (type == 'poisson') {
      newdata$UB <- exp(newdata$UB)
      newdata$LB <- exp(newdata$LB)
      newdata[,DV] <- exp(newdata[,DV])
      }

    if (type == 'binomial') {
      newdata$UB <- plogis(newdata$UB)
      newdata$LB <- plogis(newdata$LB)
      newdata[,DV] <- plogis(newdata[,DV])
      }

    return(newdata)
    }






# calculate lowest max
maxes <- df_stress %>%
  group_by(., target) %>%
  summarize(max = max(binAdj)) %>%
  as.data.frame(.)
binAdjMaxMin <- min(maxes$max)

# calculate highest low
mins <- df_stress %>%
  group_by(., target) %>%
  summarize(min = min(binAdj)) %>%
  as.data.frame(.)
binAdjMinMax <- max(mins$min)

# subset data based on new ranges
df_short <- df_stress %>% filter(., binAdj <= binAdjMaxMin & binAdj >= binAdjMinMax)

# create new adjusted variable that ranges from 1 to max
df_short$binREadj <- (df_short$binAdj - binAdjMinMax) + 1









## @knitr binAdjustments

# Bin adjustments

# Where does target suffix begin in time course
suffixOnsets <- df_stress %>%
  group_by(., target) %>%
  summarize(., sufOnset = unique(word3_suffix) / 10)

# Center time course so that suffix onset = 0
suffixOnsetAdj <- df_stress %>%
  group_by(., target) %>%
  summarize(., sufOnsetAdj = ((unique(word3_suffix) / 10) -
                              (unique(word3_suffix) / 10)))

# Where does the target word begin in the time course?
twOnsets <- df_stress %>%
  group_by(., target) %>%
  summarize(., twOnset = unique(word3_c1v1) / 10)


# Adjust to centered time course
twOnsetAdj <- cbind(suffixOnsets, twOnsets[, 2])
twOnsetAdj <- mutate(twOnsetAdj, diff = sufOnset - twOnset, twOnsetAdj = 0 - diff)

df_stress %>%
  ggplot(., aes(x = bin, y = target)) +
  geom_path(size = 1) +
  geom_point(data = suffixOnsets, aes(x = sufOnset, y = target), color = 'red') +
  geom_point(data = twOnsets, aes(x = twOnset, y = target), color = 'blue') +
  xlab("Time course") + ylab("Items") +
  theme_bw(base_size = 12, base_family = "Times") -> onsetP

df_stress %>%
  ggplot(., aes(x = binAdj, y = target)) +
  geom_path(size = 1) +
  geom_point(data = suffixOnsetAdj, aes(x = sufOnsetAdj, y = target), color = 'red') +
  geom_point(data = twOnsetAdj, aes(x = twOnsetAdj, y = target), color = 'blue') +
  geom_vline(xintercept = binAdjMinMax, color = 'red') +
  geom_vline(xintercept = binAdjMaxMin, color = 'red') +
  scale_y_discrete(position = "right") +
  xlab("Adjusted time course") + ylab("Items") +
  theme_bw(base_size = 12, base_family = "Times") -> onsetAdjP

stressBinPlots <- plot_grid(onsetP, onsetAdjP, ncol = 2)

# ggsave('stressBinPlots.png', plot = stressBinPlots, dpi = 600, device = "png",
#          path = "./mySources/figs/stress/general", height = 6.5, width = 10, unit = "in")



















df_stress_50 %>%
  group_by(group) %>%
  summarize(n = n_distinct(participant))

df_stress_50 <- read_csv("./mySources/data/clean/stressBIN50iaClean_20190513.csv") %>%
  select(., participant, group, target, condition, coda, bin, binAdj,
            binTonsetAlign, binTsuffixAlign, targetCount, distractorCount,
            targetProp, distractorProp, eLog, wts, corr) %>%
  mutate(., group = factor(group, levels = c("ss", "la", "int"))) %>%
  filter(., corr == 1,
            group %in% c('int', 'la', 'ss'),
            !participant %in% c('L01', 'L02', 'L03', 'L04', 'L05',
                                 'L06', 'L07', 'L08', 'L09', 'L10',
                                 'L15', 'L20', 'L21', 'L22', 'L23',
                                 'L26', 'L30', 'L31', 'L33', 'LA04',
                                 'LA06', 'LA07', 'LA14'),
            #binTonsetAlign >= 24 & binTonsetAlign <= 49,
            !target %in% c('cambia', 'cambió'))




condition_names <- c(
                    `stressed` = 'Paroxytone',
                    `unstressed` = 'Oxytone',
                    `0` = 'CV',
                    `1` = 'CVC'
                    )

(df_stress_50 %>%
  na.omit(.) %>%
  filter(., binTonsetAlign >= -2, binTonsetAlign <= 19) %>%
  ggplot(., aes(x = binTonsetAlign, y = targetProp, color = group, shape = group)) +
  facet_grid(. ~ coda, labeller = as_labeller(condition_names)) +
  geom_hline(yintercept = 0.5, color = 'white', size = 2) +
  geom_vline(xintercept = 11, color = 'grey40', lty = 3) +
  geom_vline(xintercept = 5, color = 'grey40', lty = 3) +
  stat_summary(fun.data = mean_cl_boot, geom = 'pointrange',  size = 0.75,
               fun.args = list(conf.int = .95, B = 1000)) +
  stat_summary(fun.y = mean, geom = 'point', color = 'white', alpha = 0.3, size = 1.75) +
  scale_shape_manual(name = "", values = 17:15, labels = c("SS", "LA", "INT")) +
  scale_color_brewer(palette = 'Set1', name = "", labels = c("SS", "LA", "INT")) +
  scale_x_continuous(breaks = c(-4, 1, 6, 11, 16), labels = c("-750", "-500", "-250", "0", "250")) +
  labs(y = 'Proportion of target fixations',
       x = 'Time relative to target syllable offset (ms)',
       caption = "Mean +/- 95% CI") +
  coord_cartesian(ylim = c(0, 1)) +
  annotate("text", x = 11.35, y = 0.02, label = 'Target syllable offset',
           angle = 90, size = 3, hjust = 0) +
  annotate("text", x = 5.35, y = 0.02, label = 'Mean target word onset',
           angle = 90, size = 3, hjust = 0) +
  theme_grey(base_size = 12, base_family = "Times") -> stressP1)

# ggsave('stressP1.png', plot = stressP1, dpi = 600, device = "png",
          # path = "./mySources/figs/stress/s3_adv_int_nat/eye_track",
          # height = 4, width = 9, unit = 'in')




























## @knitr stressTTests


#######################################
# t-tests                             #
# - Question: can they predict after  #
#   hearing @ onset of target suffix? #
# - 'they' refers to each group of    #
#   participants for each type of     #
#   word (paroxytone, oxytone)        #
# - This analysis does not compare    #
#   groups or conditions              #
#######################################

# We want to analyze proportion of target gaze at target onset
# so we need to make a subset of the data that only uses the
# target onset bin (adjusted 200ms for VWP)

glimpse(df_short)
df_short_temp <- df_short

df_short_temp[df_short_temp$group == 'int' & df_short_temp$condition == 'stressed', 'binTsuffixAlign'] <-
  df_short_temp[df_short_temp$group == 'int' & df_short_temp$condition == 'stressed', 'binTsuffixAlign'] - 12

df_short_temp[df_short_temp$group == 'int' & df_short_temp$condition == 'unstressed', 'binTsuffixAlign'] <-
  df_short_temp[df_short_temp$group == 'int' & df_short_temp$condition == 'unstressed', 'binTsuffixAlign'] - 12

# stress_subset_0 <- df_short_temp %>% filter(., binTsuffixAlign == 147, !target %in% c('cambia', 'cambio'))
df_short_temp <- df_short_temp %>% filter(., binTsuffixAlign == 147, !target %in% c('cambia', 'cambió'))




# Quick and dirty mean of target fixations as a function of
# group and condition (stressed, unstressed *1st syllable*)

df_short %>%
  na.omit(.) %>%
  group_by(., group, condition) %>%
  summarise(., meanFix = mean(targetProp))


# We will test this for each group in each condition (stressed, untressed)
# using a one-sided t-test. Specifically, we are testing the
# hypothesis that the proportion of looks is greater than
# chance (50%).
# - H0: u = 0.50
# - Ha: u > 0.50
# The generic code is: t.test(myVector, alternative = "greater", my = 0.33, conf.level = 0.95)

stress_ttest <- stress_subset_0 %>%
  na.omit(.) %>%
  group_by(., group, condition, participant) %>%
  summarise(., meanFix = mean(targetProp)) %>%
  do(tidy(t.test(.$meanFix, alternative = "greater", mu = 0.5, conf.level = 0.99)))

# Convert pvalues from scientific notation
stress_ttest$p.value <- format(stress_ttest$p.value, scientific = F)
stress_ttest$sig <- "N.S."
stress_ttest[stress_ttest$p.value <= 0.05/6, 'sig'] <- "*"

# Print results
print(as.data.frame(stress_ttest[, c(1:7, 11)]))

saveRDS(stress_ttest, "./mySources/reports/stress/stress_int/mods/stress_ttest.rds", compress = "xz")

# group  condition  estimate  statistic       p.value parameter  conf.low  sig
#   int   stressed 0.6344643 3.46115554 0.00357388323         9 0.5248527    *
#   int unstressed 0.6753571 3.37992533 0.00406396949         9 0.5289754    *
#    la   stressed 0.5022928 0.05122216 0.47976999874        26 0.3913462 N.S.
#    la unstressed 0.5841931 2.12785440 0.02149710485        26 0.4861208 N.S.
#    ss   stressed 0.6322078 2.79896749 0.00537726129        21 0.5132880    *
#    ss unstressed 0.7235390 4.76754875 0.00005197742        21 0.6054925    *


# We will plot the models
# This will be almost exactly the same as
# the previous plot, but it will use the
# confidence interval from the test we
# actually conducted

(stress_ttest %>%
  ungroup(.) %>%
  mutate(., group = factor(group, levels = c("la", "int", "ss")),
            condition = factor(condition, levels = c("stressed", "unstressed"))) %>%
  ggplot(., aes(x = group, y = estimate, color = condition,
                group = interaction(group, condition), dodge = condition)) +
    geom_hline(yintercept = 0.5, lty = 3) +
    geom_linerange(aes(ymin = conf.low, ymax = estimate), color = 'grey40',
                       position = position_dodge(width = 0.75), size = 1) +
    geom_point(position = position_dodge(width = 0.75), size = 4) +
    geom_point(position = position_dodge(width = 0.75), size = 2.75, color = 'grey90') +
    ylim(0.3, 1.0) +
    scale_color_brewer(palette = "Set1", name = '', labels = c('Oxytone', 'Paroxytone')) +
    scale_x_discrete(labels = c('LA', 'INT', 'SS')) +
    labs(title = 'Mean fixations at target syllable offset',
         y = 'Target fixations', x = '', caption = 'Mean and lower-bound 99% CI') +
    theme_minimal(base_size = 12, base_family = 'Times') -> stressTargetFixMOD1)

# ggsave('stressTargetFixMOD1.png', plot = stressTargetFixMOD1, dpi = 600, device = "png",
#         path = "./mySources/figs/stress/s3_adv_int_nat/eye_track",
#         height = 4, width = 7, unit = 'in')






















## @knitr stressGroupCompare


# Load wm data and combine with stress_subset_0 (proportion data)
# in order to add working memory as a covariate
wm_df <- read_csv("./mySources/data/raw/wm_all.csv") %>%
         filter(., !(group %in% c("HS", "L")),
                   !(participant %in% c('L01', 'L02', 'L03', 'L04', 'L05',
                                        'L06', 'L07', 'L08', 'L09', 'L10',
                                        'L15', 'L20', 'L21', 'L22', 'L23',
                                        'L26', 'L30', 'L31', 'L33', 'La04',
                                        'LA06', 'LA07', 'LA14')))

scale_this <- function(x) as.vector(scale(x))


# Are groups different from each other?
stress_subset_0_prop <- stress_subset_0 %>%
  na.omit(.) %>%
  select(., group, condition, participant, target, coda,
            targetCount, distractorCount, eLog, wts, targetProp) %>%
  # summarise(., meanFix = mean(targetProp)) %>%
 # ungroup(.) %>%
 # left_join(x = ., y = wm_df[, -1], by = 'participant') %>%
 mutate(., group = factor(group, levels = c("ss", "la", "int")),
           coda = as.factor(coda)) # %>%
 # group_by(., group) %>%
 # mutate(., wmScaled = scale_this(WM))

stress_subset_0_prop$coda <- as.factor(stress_subset_0_prop$coda)
stress_subset_0_prop$codaSum <- C(stress_subset_0_prop$coda, sum)
# contrasts(stress_subset_0_prop$codaSum)

stress_subset_0_prop$condition <- as.factor(stress_subset_0_prop$condition)
stress_subset_0_prop$conditionSum <- C(stress_subset_0_prop$condition, sum)
# contrasts(stress_subset_0_prop$codaSum)



# random effects building

prop_0_ranefA <- glmer(cbind(targetCount, 10 - targetCount) ~ 1 +
                    (1 | participant),
                    data = stress_subset_0_prop, REML = F, family = 'binomial',
                    control = glmerControl(optimizer = 'bobyqa'))

prop_0_ranefB <- glmer(cbind(targetCount, 10 - targetCount) ~ 1 +
                    (1 | participant) + (1 | target),
                    data = stress_subset_0_prop, family = 'binomial',
                    control = glmerControl(optimizer = 'bobyqa'))

anova(prop_0_ranefA, prop_0_ranefB, refit = F) # keep intercept for target

prop_0_ranefC <- glmer(cbind(targetCount, 10 - targetCount) ~ 1 +
                    (1 | participant) + (1 | target) +
                    (0 + conditionSum | participant),
                    data = stress_subset_0_prop, family = 'binomial',
                    control = glmerControl(optimizer = 'bobyqa'))

anova(prop_0_ranefB, prop_0_ranefC, refit = F) # keep uncorrelated slope for condition

prop_0_ranefD <- glmer(cbind(targetCount, 10 - targetCount) ~ 1 +
                    (1 | participant) + (1 | target) +
                    (0 + conditionSum | participant) + (0 + codaSum | participant),
                    data = stress_subset_0_prop, family = 'binomial',
                    control = glmerControl(optimizer = 'bobyqa'))

anova(prop_0_ranefC, prop_0_ranefD, refit = F) # Keep uncorrelated slope for coda

prop_0_ranefE <- glmer(cbind(targetCount, 10 - targetCount) ~ 1 +
                    (1 | participant) + (1 | target) +
                    (0 + conditionSum*codaSum | participant),
                    data = stress_subset_0_prop, family = 'binomial',
                    control = glmerControl(optimizer = 'bobyqa'))

anova(prop_0_ranefD, prop_0_ranefE, refit = F) # Keep condition x coda interaction slope

prop_0_ranefF <- glmer(cbind(targetCount, 10 - targetCount) ~ 1 +
                    (1 | participant) + (1 | target) +
                    (1 + conditionSum*codaSum | participant),
                    data = stress_subset_0_prop, family = 'binomial',
                    control = glmerControl(optimizer = 'bobyqa'))

anova(prop_0_ranefE, prop_0_ranefF, refit = F) # Keep correlated slope



# test fixed effects

prop_0_mod_0 <- glmer(cbind(targetCount, 10 - targetCount) ~ 1 +
                    (1 | participant) + (1 | target) +
                    (1 + conditionSum*codaSum | participant),
                    data = stress_subset_0_prop, family = 'binomial',
                    control = glmerControl(optimizer = 'bobyqa'))

prop_0_mod_group <- glmer(cbind(targetCount, 10 - targetCount) ~ 1 + group +
                    (1 | participant) + (1 | target) +
                    (1 + conditionSum*codaSum | participant),
                    data = stress_subset_0_prop, family = 'binomial',
                    control = glmerControl(optimizer = 'bobyqa'))

prop_0_mod_cond <- glmer(cbind(targetCount, 10 - targetCount) ~ 1 + group + conditionSum +
                    (1 | participant) + (1 | target) +
                    (1 + conditionSum*codaSum | participant),
                    data = stress_subset_0_prop, family = 'binomial',
                    control = glmerControl(optimizer = 'bobyqa'))

prop_0_mod_coda <- glmer(cbind(targetCount, 10 - targetCount) ~ 1 + group + codaSum +
                    (1 | participant) + (1 | target) +
                    (1 + conditionSum*codaSum | participant),
                    data = stress_subset_0_prop, family = 'binomial',
                    control = glmerControl(optimizer = 'bobyqa'))

prop_0_mod_int1 <- glmer(cbind(targetCount, 10 - targetCount) ~ 1 + group + codaSum + group:codaSum +
                    (1 | participant) + (1 | target) +
                    (1 + conditionSum*codaSum | participant),
                    data = stress_subset_0_prop, family = 'binomial',
                    control = glmerControl(optimizer = 'bobyqa'))

prop_0_mod_int2 <- glmer(cbind(targetCount, 10 - targetCount) ~ 1 + group + codaSum + group:conditionSum +
                    (1 | participant) + (1 | target) +
                    (1 + conditionSum*codaSum | participant),
                    data = stress_subset_0_prop, family = 'binomial',
                    control = glmerControl(optimizer = 'bobyqa'))

prop_0_mod_int3 <- glmer(cbind(targetCount, 10 - targetCount) ~ 1 + group + codaSum + conditionSum:codaSum +
                    (1 | participant) + (1 | target) +
                    (1 + conditionSum*codaSum | participant),
                    data = stress_subset_0_prop, family = 'binomial',
                    control = glmerControl(optimizer = 'bobyqa'))

prop_0_mod_full <- glmer(cbind(targetCount, 10 - targetCount) ~ 1 + group * conditionSum * codaSum +
                    (1 | participant) + (1 | target) +
                    (1 + conditionSum*codaSum | participant),
                    data = stress_subset_0_prop, family = 'binomial',
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


stress_subset_0_prop$group <- factor(stress_subset_0_prop$group, levels = c("ss", "la",  "int"))

prop_0_mod_final <- glmer(cbind(targetCount, 10 - targetCount) ~ 1 + group + codaSum + conditionSum +
                    (1 | participant) + (1 | target) +
                    (1 + conditionSum*codaSum | participant),
                    data = stress_subset_0_prop, family = 'binomial',
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
stress_subset_0_prop$group <- factor(stress_subset_0_prop$group, levels = c("int", "la",  "ss"))

summary(glmer(cbind(targetCount, 10 - targetCount) ~ 1 + group + codaSum +
       (1 | participant) + (1 | target) +
       (1 + conditionSum*codaSum | participant),
       data = stress_subset_0_prop, family = 'binomial',
       control = glmerControl(optimizer = 'bobyqa')))

# Fixed effects:
#             Estimate Std. Error z value Pr(>|z|)
# (Intercept)   1.0620     0.3840   2.765  0.00568 **
# groupla      -0.6426     0.4049  -1.587  0.11249
# groupss       0.4879     0.4244   1.150  0.25028
# codaSum1     -0.6111     0.2170  -2.816  0.00486 **




et_ci <- confint(prop_0_mod_final, method = "Wald", level = 0.99) %>%
  as.data.frame(.) %>%
  slice(., 13:16) %>%
  rename(., ciLow = `0.5 %`, ciHi = `99.5 %`)

stressFixModP0 <- broom::tidy(prop_0_mod_final) %>% slice(1:4) %>%
  cbind(., et_ci) %>%
  mutate(., term = recode(term, `(Intercept)` = '(Intercept)',
                                codaSum1 = 'Syllable\nstructure',
                                groupla = 'NIN',
                                groupint = 'IN'),
            term = factor(term, levels = c('IN',
                                           'NIN',
                                           'Syllable\nstructure',
                                           '(Intercept)'))) %>%
  ggplot(., aes(x = estimate, y = term)) +
    geom_vline(xintercept = 0, lty = 3) +
    geom_errorbarh(aes(xmin = ciLow, xmax = ciHi), height = 0.2, size = 0.65) +
    geom_point(size = 3) +
    geom_point(size = 2, color = 'lightgrey') +
    labs(y = 'Term', x = 'Estimate +/- 95% CI') +
    theme_bw(base_size = 15, base_family = 'Times') +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

# ggsave('stressP0.png', plot = stressFixModP0, dpi = 600, device = "png",
#         path = "./mySources/figs/stress/s3_adv_int_nat/eye_track",
#         height = 2.5, width = 7, unit = 'in')





stressFixModsHLS <- stress_subset_0_prop %>%
  mutate(., group = factor(group, levels = c("la", "int",  "ss"))) %>%
  group_by(., group, coda, participant) %>%
  summarise(., meanFix = mean(targetProp)) %>%
  ggplot(., aes(x = group, y = meanFix,
                group = interaction(group, coda), dodge = coda, color = coda)) +
    geom_hline(yintercept = 0.5, lty = 3) +
    stat_summary(fun.data = mean_cl_boot, geom = 'pointrange',
                 position = position_dodge(width = 0.5), color = 'black',
                 size = 0.90, show.legend = FALSE, fun.args = list(conf.int = 0.95)) +
    stat_summary(fun.y = mean, geom = 'point', size = 2.75,
                 position = position_dodge(width = 0.5), show.legend = TRUE) +
    labs(y = 'Target fixations', x = '', caption = 'Mean +/- 99% bootstrap CI.',
         title = "Mean fixations at target syllable offset") +
    coord_cartesian(ylim = c(0.4, 1)) +
    scale_x_discrete(labels = c('LA', 'INT', 'SS')) +
    scale_color_brewer(name = '', palette = 'Set1', labels = c("CV", "CVC")) +
    theme_minimal(base_size = 12, base_family = 'Times')

# ggsave('stressFixModsHLS.png', plot = stressFixModsHLS, dpi = 600, device = "png",
#         path = "./mySources/figs/stress/s3_adv_int_nat/eye_track",
#         height = 4, width = 7, unit = 'in')











# Calculate mean target fixation as a function of group, condition,
# for each participant. We will plot the mean and calculate the
# bootstrapped 95% confidence interval and plot it all.

stressFixModsP1 <- stress_subset_0_prop %>%
  mutate(., group = factor(group, levels = c("la", "int",  "ss"))) %>%
  group_by(., group, condition, participant) %>%
  summarise(., meanFix = mean(targetProp)) %>%
  ggplot(., aes(x = group, y = meanFix, shape = condition,
                group = interaction(group, condition), dodge = condition, color = group)) +
    geom_hline(yintercept = 0.5, lty = 3) +
    stat_summary(fun.data = mean_cl_boot, geom = 'pointrange',
                 position = position_dodge(width = 0.5), color = 'black',
                 size = 0.90, show.legend = FALSE, fun.args = list(conf.int = 0.99)) +
    stat_summary(fun.y = mean, geom = 'point', size = 2.75,
                 position = position_dodge(width = 0.5)) +
    ylim(0, 1) +
    labs(y = '% Correct', x = '', caption = '') +
    scale_x_discrete(labels = c('NIN', 'IN', 'M')) +
    scale_color_manual(name = '', values = c('grey90', 'grey75', 'grey55'), guide = FALSE) +
    scale_shape_manual(name = '', values = c(16, 17), labels = c('Paroxytone', 'Oxytone')) +
    guides(shape = guide_legend(override.aes = list(shape = c(1, 2), color = 'black'))) +
    theme_bw(base_size = 15, base_family = 'Times') +
    theme(legend.position = c(0.26, 0.14),
          legend.box.just = "left",
          legend.background = element_rect(fill = "transparent"),
          legend.key = element_rect(fill = "transparent", colour = "transparent"),
          legend.key.size = unit(0.75, 'lines'),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank())


stressFixModsP2 <- stress_subset_0_prop %>%
  mutate(., group = factor(group, levels = c("la", "int",  "ss"))) %>%
  group_by(., group, coda, participant) %>%
  summarise(., meanFix = mean(targetProp)) %>%
  ggplot(., aes(x = group, y = meanFix, shape = coda,
                group = interaction(group, coda), dodge = coda, color = group)) +
    geom_hline(yintercept = 0.5, lty = 3) +
    stat_summary(fun.data = mean_cl_boot, geom = 'pointrange',
                 position = position_dodge(width = 0.5), color = 'black',
                 size = 0.90, show.legend = FALSE, fun.args = list(conf.int = 0.99)) +
    stat_summary(fun.y = mean, geom = 'point', size = 2.75,
                 position = position_dodge(width = 0.5)) +
    labs(y = '% Correct', x = '', caption = '+/- 99% bootstrap CI.') +
    scale_y_continuous(position = "right", limits = c(0, 1)) +
    scale_x_discrete(labels = c('NIN', 'IN', 'M')) +
    scale_color_manual(name = '', values = c('grey90', 'grey75', 'grey55'), guide = FALSE) +
    scale_shape_manual(name = '', values = c(16, 17), labels = c('CV', 'CVC')) +
    guides(shape = guide_legend(override.aes = list(shape = c(1, 2), color = 'black'))) +
    theme_bw(base_size = 15, base_family = 'Times') +
    theme(legend.position = c(0.16, 0.14),
          legend.box.just = "left",
          legend.background = element_rect(fill = "transparent"),
          legend.key = element_rect(fill = "transparent", colour = "transparent"),
          legend.key.size = unit(0.75, 'lines'),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank())


# arrange all plots together
stress_eyet_plot <- grid.arrange(stressFixModsP1,
                                 stressFixModsP2,
                                   ncol = 2)

# ggsave('stress_eyet_plot.png',
#        plot = stress_eyet_plot, dpi = 600, device = "png",
#       path = "./mySources/figs/stress/s3_adv_int_nat/eye_track",
#        height = 3.5, width = 8.5, units = 'in')















































glimpse(df_stress)
glimpse(df_short_temp)

# calculate proportion of target fixations by sub, by target, by stress cond, by group
# for each landmark in the time course

df_timecourse <- df_short_temp %>%
  select(., participant, group, bin, BIN_START_TIME, BIN_END_TIME, target,
            startsentence:target, condition, coda, targetCount:targetProp, eLog) %>%
  mutate(., bin = bin * 10)

range(df_timecourse$bin)

df_timecourse_startsentence <- df_timecourse %>%
  filter(., startsentence + 1550 >= BIN_START_TIME, startsentence + 1550 <= BIN_END_TIME) %>%
  select(., participant, group, target, condition:eLog) %>%
  mutate(., landmark = 'start_sentence')

df_timecourse_word2_c1v1 <- df_timecourse %>%
  filter(., word2_c1v1 + 1550 >= BIN_START_TIME, word2_c1v1 + 1550 <= BIN_END_TIME) %>%
  select(., participant, group, target, condition:eLog) %>%
  mutate(., landmark = 'word2_c1v1')


df_timecourse_word3_20msafterv1 <- df_timecourse %>%
  filter(., word3_20msafterv1 + 1550 >= BIN_START_TIME, word3_20msafterv1 + 1550 <= BIN_END_TIME) %>%
  select(., participant, group, target, condition:eLog) %>%
  mutate(., landmark = 'word3_20msafterv1')


df_timecourse_word3_c1v1 <- df_timecourse %>%
  filter(., word3_c1v1 + 1550 >= BIN_START_TIME, word3_c1v1 + 1550 <= BIN_END_TIME) %>%
  select(., participant, group, target, condition:eLog) %>%
  mutate(., landmark = 'word3_c1v1')


df_timecourse_word3_c2 <- df_timecourse %>%
  filter(., word3_c2 + 1550 >= BIN_START_TIME, word3_c2 + 1550 <= BIN_END_TIME) %>%
  select(., participant, group, target, condition:eLog) %>%
  mutate(., landmark = 'word3_c2')

df_timecourse_word3_c3 <- df_timecourse %>%
  filter(., word3_c3 + 1550 >= BIN_START_TIME, word3_c3 + 1550 <= BIN_END_TIME)  %>%
  select(., participant, group, target, condition:eLog) %>%
  mutate(., landmark = 'word3_c3')

df_timecourse_word3_suffix <- df_timecourse %>%
  filter(., word3_suffix + 1550 >= BIN_START_TIME, word3_suffix + 1550 <= BIN_END_TIME) %>%
  select(., participant, group, target, condition:eLog) %>%
  mutate(., landmark = 'word3_suffix')

df_timecourse_word4_c1v1 <- df_timecourse %>%
  filter(., word4_c1v1 + 1550 >= BIN_START_TIME, word4_c1v1 + 1550 <= BIN_END_TIME) %>%
  select(., participant, group, target, condition:eLog) %>%
  mutate(., landmark = 'word4_c1v1')

df_timecourse_word5 <- df_timecourse %>%
  filter(., word5 + 1550 >= BIN_START_TIME, word5 + 1550 <= BIN_END_TIME) %>%
  select(., participant, group, target, condition:eLog) %>%
  mutate(., landmark = 'word5')

df_timecourse_end_sentence <- df_timecourse %>%
  filter(., end_sentence + 1550 >= BIN_START_TIME, end_sentence + 1550 <= BIN_END_TIME)  %>%
  select(., participant, group, target, condition:eLog) %>%
  mutate(., landmark = 'end_sentence')

df_landmarks <- do.call("rbind", list(df_timecourse_startsentence,
                      df_timecourse_word2_c1v1,
                      df_timecourse_word3_20msafterv1,
                      df_timecourse_word3_c1v1,
                      df_timecourse_word3_c2,
                      df_timecourse_word3_c3,
                      df_timecourse_word3_suffix,
                      df_timecourse_word4_c1v1,
                      df_timecourse_word5,
                      df_timecourse_end_sentence))

glimpse(df_landmarks)

df_landmarks <- mutate(df_landmarks,
  landmark = factor(landmark, levels = c('start_sentence',
                                         'word2_c1v1',
                                         'word3_c1v1',
                                         'word3_20msafterv1',
                                         'word3_c2',
                                         'word3_c3',
                                         'word3_suffix',
                                         'word4_c1v1',
                                         'word5',
                                         'end_sentence')))




df_landmarks %>%
  na.omit(.) %>%
  filter(., coda == 0,
            landmark %in% c('word3_c1v1', 'word3_20msafterv1', 'word3_suffix', 'word4_c1v1')) %>%
  group_by(., participant, target, group, coda, landmark) %>%
  summarize(., target_fix = mean(targetProp)) %>%
  ggplot(., aes(x = landmark, y = target_fix, shape = group, dodge = group)) +
    geom_hline(yintercept = 0.5, color = 'black') +
    ylim(0, 1) +
    stat_summary(fun.data = mean_cl_boot, geom = 'pointrange', position = position_dodge(width = 0.3)) +
    theme_bw()

df_landmarks %>%
  na.omit(.) %>%
  filter(., coda == 1,
            landmark %in% c('word3_c1v1', 'word3_20msafterv1', 'word3_c2', 'word3_suffix', 'word4_c1v1')) %>%
  group_by(., participant, target, group, coda, landmark) %>%
  summarize(., target_fix = mean(targetProp)) %>%
  ggplot(., aes(x = landmark, y = target_fix, shape = group, dodge = group)) +
    geom_hline(yintercept = 0.5, color = 'black') +
    ylim(0, 1) +
    stat_summary(fun.data = mean_cl_boot, geom = 'pointrange', position = position_dodge(width = 0.3)) +
    theme_bw()





















































## @knitr stressGCA


####################################################
# GROWTH CURVE ANALYSIS                            #
# - Question 1: Are the groups different from      #
#   each other in when they begin to fixate        #
#   on the target?                                 #
#     - test 3 groups at each level of 'condition' #
#     - hypothesis: SS has steeper slope for both  #
#       conditions                                 #
# - Question 2: W/in groups, is the there a        #
#   difference between oxytone/paroxyton           #
#   items?                                         #
#     - test oxytone vs. paroxytone for each group #
#     - hypothesis: steeper slope/earlier break in #
#       oxytone condition                          #
####################################################

# Prep
# - subset using time course
#    - we want to find the earliest target word onset
#    - we will substract 'a little more' and use that
#      as the starting point
#    - we already have this information in the 'twOnsetAdj'
#      dataframe

# print(twOnsetAdj)
# min(twOnsetAdj$twOnsetAdj)

#    - the lowest target word onset is: 'firma' @ -51.3
#    - thus we can use -60 binAdj as the starting point
#      and be sure we are including the entire target word
#    - we can also do a little higher to lighten the models (-50)
stress_gc_subset <- filter(df_stress_50, binTonsetAlign >= 33 &
                                         binTonsetAlign <= 43) %>% as.data.frame


# - Readjust time course
#    - now we will make the time course positive, starting at 1
#    - to do this we will add the lowest 'binAdj' value to each
#      bin, plus 1 (to avoid starting at 0)

# lsrl_gc_subset$binGC <- lsrl_gc_subset$binAdj + 51
stress_gc_subset$binGC <- stress_gc_subset$binTonsetAlign - 32


# - Now we add higher order polynomials for analyses

t <- poly(min(stress_gc_subset$binGC):max(stress_gc_subset$binGC), 4)
stress_gc_subset[, paste('ot', 1:4, sep = "")] <- t[stress_gc_subset$binGC, 1:4]

# glimpse(lsrl_gc_subset)






####################################################
# - Question 1: Are the groups different from      #
#   each other in when they begin to fixate        #
#   on the target?                                 #
#     - test 3 groups at each level of 'condition' #
#     - hypothesis: SS has steeper slope for both  #
#       conditions                                 #
####################################################


# Set SS as reference level
stress_gc_subset$group <- factor(stress_gc_subset$group, levels = c("ss", "la", "int"))

# Set condition as factor and set conding to contrast
stress_gc_subset$condition <- as.factor(stress_gc_subset$condition)

# contrasts(stress_gc_subset$condition)
# contrasts(stress_gc_subset$group)

stress_gc_subset$conditionSum <- C(stress_gc_subset$condition, sum)
# contrasts(stress_gc_subset$conditionSum)

stress_gc_subset$coda <- as.factor(stress_gc_subset$coda)
stress_gc_subset$codaSum <- C(stress_gc_subset$coda, sum)
# contrasts(stress_gc_subset$codaSum)


## @knitr ignore2

# load the models:
gc_mod_base    <- readRDS('./mySources/models/stress/s3_adv_int_nat/eye_track/gc_mod_base.rds')
gc_mod_group_0 <- readRDS('./mySources/models/stress/s3_adv_int_nat/eye_track/gc_mod_group_0.rds')
gc_mod_group_1 <- readRDS('./mySources/models/stress/s3_adv_int_nat/eye_track/gc_mod_group_1.rds')
gc_mod_group_2 <- readRDS('./mySources/models/stress/s3_adv_int_nat/eye_track/gc_mod_group_2.rds')
gc_mod_group_3 <- readRDS('./mySources/models/stress/s3_adv_int_nat/eye_track/gc_mod_group_3.rds')
gc_mod_cond_0  <- readRDS('./mySources/models/stress/s3_adv_int_nat/eye_track/gc_mod_cond_0.rds')
gc_mod_cond_1  <- readRDS('./mySources/models/stress/s3_adv_int_nat/eye_track/gc_mod_cond_1.rds')
gc_mod_cond_2  <- readRDS('./mySources/models/stress/s3_adv_int_nat/eye_track/gc_mod_cond_2.rds')
gc_mod_cond_3  <- readRDS('./mySources/models/stress/s3_adv_int_nat/eye_track/gc_mod_cond_3.rds')
gc_mod_full    <- readRDS('./mySources/models/stress/s3_adv_int_nat/eye_track/gc_mod_full_0.rds')


# random effects structure
mod_ot1 <- lmer(eLog ~ ot1 + (1 + ot1 | participant) +
                (1 + ot1 | participant:codaSum),
                control = lmerControl(optimizer = 'bobyqa'),
                data = stress_gc_subset, weights = 1/wts, REML = F)
mod_ot2 <- lmer(eLog ~ ot1 + ot2 + (1 + (ot1+ot2) | participant) +
                (1 + (ot1+ot2) | participant:codaSum),
                control = lmerControl(optimizer = 'bobyqa'),
                data = stress_gc_subset, weights = 1/wts, REML = F)
mod_ot3 <- lmer(eLog ~ ot1 + ot2 + ot3 + (1 + (ot1+ot2+ot3) | participant) +
                (1 + (ot1+ot2+ot3) | participant:codaSum),
                control = lmerControl(optimizer = 'bobyqa'),
                data = stress_gc_subset, weights = 1/wts, REML = F)

anova(mod_ot1, mod_ot2, mod_ot3) # best fit with quadratic poly

mod_ot2_simp <- lmer(eLog ~ ot1 + ot2 + (1 + (ot1+ot2) | participant),
                control = lmerControl(optimizer = 'bobyqa'),
                data = stress_gc_subset, weights = 1/wts, REML = F)

mod_ot2_cond <- lmer(eLog ~ ot1 + ot2 + (1 + (ot1+ot2) | participant) +
                (1 + (ot1+ot2) | participant:conditionSum),
                control = lmerControl(optimizer = 'bobyqa'),
                data = stress_gc_subset, weights = 1/wts, REML = F)

anova(mod_ot2_simp, mod_ot2_cond) # cond is better than just part

mod_ot2_coda <- lmer(eLog ~ ot1 + ot2 + (1 + (ot1+ot2) | participant) +
                (1 + (ot1+ot2) | participant:codaSum),
                control = lmerControl(optimizer = 'bobyqa'),
                data = stress_gc_subset, weights = 1/wts, REML = F)

anova(mod_ot2_simp, mod_ot2_coda) # coda is better than simp

mod_ot2_max <- lmer(eLog ~ ot1 + ot2 + (1 + (ot1+ot2) | participant) +
                (1 + (ot1+ot2) | participant:codaSum:conditionSum),
                control = lmerControl(optimizer = 'bobyqa'),
                data = stress_gc_subset, weights = 1/wts, REML = F)

anova(mod_ot2_coda, mod_ot2_max) # coda is better than simp








# Base model
if(T){
  gc_mod_base <- lmer(eLog ~ (ot1+ot2) +
                 (1 + (ot1+ot2) | participant) +
                 (1 + (ot1+ot2) | participant:codaSum),
                 control = lmerControl(optimizer = 'bobyqa'),
                 data = stress_gc_subset, weights = 1/wts, REML = F)
  saveRDS(gc_mod_base, file = "./mySources/models/stress/s3_adv_int_nat/eye_track/gc_mod_base.rds", compress = 'xz')
}




# Add group effect on intercept
if(T){
  gc_mod_group_0 <- lmer(eLog ~ (ot1+ot2) + group +
                    (1 + (ot1+ot2) | participant) +
                    (1 + (ot1+ot2) | participant:codaSum),
                    control = lmerControl(optimizer = 'bobyqa'),
                    data = stress_gc_subset, weights = 1/wts, REML = F)
  saveRDS(gc_mod_group_0, file = "./mySources/models/stress/s3_adv_int_nat/eye_track/gc_mod_group_0.rds", compress = 'xz')
}

# Add group effect on slope
if(T){
 gc_mod_group_1 <- lmer(eLog ~ (ot1+ot2) + group +
                   ot1:group +
                   (1 + (ot1+ot2) | participant) +
                   (1 + (ot1+ot2) | participant:codaSum),
                   control = lmerControl(optimizer = 'bobyqa'),
                   data = stress_gc_subset, weights = 1/wts, REML = F)
  saveRDS(gc_mod_group_1, file = "./mySources/models/stress/s3_adv_int_nat/eye_track/gc_mod_group_1.rds", compress = 'xz')
}

# Add group effect on quadratic poly
if(T){
  gc_mod_group_2 <- lmer(eLog ~ (ot1+ot2) + group +
                    ot1:group + ot2:group +
                    (1 + (ot1+ot2) | participant) +
                    (1 + (ot1+ot2) | participant:codaSum),
                    control = lmerControl(optimizer = 'bobyqa'),
                    data = stress_gc_subset, weights = 1/wts, REML = F)
  saveRDS(gc_mod_group_2, file = "./mySources/models/stress/s3_adv_int_nat/eye_track/gc_mod_group_2.rds", compress = 'xz')
}

anova(gc_mod_base,
      gc_mod_group_0,
      gc_mod_group_1,
      gc_mod_group_2, test = 'Chisq')




# full mod
if(T){
gc_mod_full_0 <- lmer(eLog ~ (ot1+ot2+ot3) * group * codaSum +
               (1 + (ot1+ot2+ot3) | participant) +
               # (1 + (ot1+ot2) | participant:codaSum) +
               (1 + (ot1+ot2+ot3) | participant:conditionSum),
               control = lmerControl(optimizer = 'bobyqa'),
               data = stress_gc_subset, weights = 1/wts, REML = F)
  saveRDS(gc_mod_full_0, file = "./mySources/models/stress/s3_adv_int_nat/eye_track/gc_mod_full_0.rds", compress = 'xz')
}

if(T){
gc_mod_full_1 <- lmer(eLog ~ (ot1+ot2) * group * codaSum + conditionSum +
               (1 + (ot1+ot2) | participant) +
               # (1 + (ot1+ot2) | participant:codaSum) +
               (1 + (ot1+ot2) | participant:conditionSum),
               control = lmerControl(optimizer = 'bobyqa'),
               data = stress_gc_subset, weights = 1/wts, REML = F)
  saveRDS(gc_mod_full_1, file = "./mySources/models/stress/s3_adv_int_nat/eye_track/gc_mod_full_1.rds", compress = 'xz')
}

anova(gc_mod_full_0, gc_mod_full_1) # no main effect of cond
# ..1    32 60048 60281 -29992    59984 3.386      1    0.06575 .



# summary(gc_mod_full_0)



# Fixed effects:
#                         Estimate Std. Error         df t value Pr(>|t|)
# (Intercept)              0.90951    0.14791   68.09840   6.149 4.63e-08 ***
# ot1                      0.31167    0.33240 8438.55185   0.938 0.348457
# ot2                     -0.16466    0.23044   54.49578  -0.715 0.477930

# groupla                 -0.38784    0.19936   68.23686  -1.945 0.055843 .
# ot1:groupla             -0.80774    0.50262   67.96785  -1.607 0.112680
# ot2:groupla              1.04344    0.31158   55.32223   3.349 0.001466 **

# groupint                -0.09539    0.26477   68.37041  -0.360 0.719754
# ot1:groupint             0.07018    0.67222   70.37943   0.104 0.917153
# ot2:groupint            -0.20198    0.41851   59.32022  -0.483 0.631149

# codaSum1                 0.07093    0.05074 9317.63184   1.398 0.162166
# ot1:codaSum1             0.27608    0.17746 7863.15846   1.556 0.119827
# ot2:codaSum1             0.06745    0.17106 6075.19780   0.394 0.693379

# groupla:codaSum1        -0.17768    0.06837 9279.03361  -2.599 0.009371 **
# ot1:groupla:codaSum1    -0.33628    0.23924 7949.21744  -1.406 0.159880
# ot2:groupla:codaSum1     0.07746    0.23258 6122.17549   0.333 0.739124

# groupint:codaSum1       -0.31439    0.09270 9243.38468  -3.392 0.000698 ***
# ot1:groupint:codaSum1    3.30410    0.37224   67.43673   8.876 6.05e-13 ***
# ot2:groupint:codaSum1    0.04064    0.31389 7094.47346   0.129 0.896999




# create new df including the fitted model
data.comp <- data.frame(na.omit(stress_gc_subset),
                        GCA_Full = fitted(gc_mod_full_0))
# glimpse(data.comp)


condition_namesGCAMod <- c(
                    `0` = "CV",
                    `1` = "CVC"
                    )

(gca_full <- data.comp %>%
  ggplot(., aes(x = binTonsetAlign, y = eLog, color = group, shape = group)) +
  facet_grid(. ~ coda, labeller = as_labeller(condition_namesGCAMod)) +
  geom_smooth(method = 'gam', formula = y ~ poly(x, 3), se = F, show.legend = FALSE) +
  # stat_summary(aes(y = GCA_Full, color = group), fun.y = mean, geom = 'line', size = 0.4) +
  stat_summary(fun.data = mean_cl_boot, geom = 'pointrange',  size = 0.75,
               fun.args = list(conf.int = .95, B = 1000)) +
  stat_summary(fun.y = mean, geom = 'point', color = 'white', alpha = 0.3, size = 1.75) +
  scale_shape_manual(name = "", values = 17:15, labels = c("M", "NIN", "IN")) +
  labs(x = "Time relative to target syllable offset (ms)", y = "Fixation empirical logit",
       caption = "Mean +/- 95% CI") +
  scale_color_brewer(palette = "Set1", name = "", guide = 'legend',
                     labels = c("M", "NIN", "IN")) +
  scale_x_continuous(breaks = c(33, 38, 43), labels = c("-250", "0", "250")) +
  theme_grey(base_size = 15, base_family = "Times New Roman"))

# ggsave('stressP2.png', plot = gca_full, dpi = 600, device = "png",
#           path = "./mySources/figs/stress/s3_adv_int_nat/eye_track",
#           height = 4, width = 9, unit = "in")













# OLD OUTPUT

# Fixed effects:
#                         Estimate Std. Error         df t value Pr(>|t|)
# (Intercept)              0.89485    0.14747   68.00000   6.068 6.41e-08 ***
# ot1                      3.33309    0.37630   67.00000   8.858 7.25e-13 ***
# ot2                     -0.16281    0.23645   55.00000  -0.689  0.49401
#
# groupla                 -0.37843    0.19871   68.00000  -1.904  0.06107 .
# ot1:groupla             -0.84030    0.50822   67.00000  -1.653  0.10291
# ot2:groupla              1.04784    0.31978   55.00000   3.277  0.00182 **
#
# groupint                -0.11347    0.26369   68.00000  -0.430  0.66832
# ot1:groupint             0.07456    0.67747   69.00000   0.110  0.91268
# ot2:groupint            -0.18355    0.42840   59.00000  -0.428  0.66989
#
# codaSum1                 0.07384    0.05076 9377.00000   1.455  0.14578
# ot1:codaSum1             0.25070    0.17745 8565.00000   1.413  0.15777
# ot2:codaSum1             0.05033    0.17117 6255.00000   0.294  0.76874
#
# groupla:codaSum1        -0.17708    0.06821 9430.00000  -2.596  0.00945 **
# ot1:groupla:codaSum1    -0.33634    0.23883 8657.00000  -1.408  0.15909
# ot2:groupla:codaSum1     0.09856    0.23274 6338.00000   0.423  0.67196
#
# groupint:codaSum1       -0.29203    0.09172 9378.00000  -3.184  0.00146 **
# ot1:groupint:codaSum1    0.32001    0.32918 8746.00000   0.972  0.33102
# ot2:groupint:codaSum1    0.07729    0.31381 7054.00000   0.246  0.80545


# ot2:groupla              1.04784    0.31978   55.00000   3.277  0.00182 **
# groupla:codaSum1        -0.17708    0.06821 9430.00000  -2.596  0.00945 **
# groupint:codaSum1       -0.29203    0.09172 9378.00000  -3.184  0.00146 **




