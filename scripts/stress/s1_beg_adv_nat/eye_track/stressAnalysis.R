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

# Set working directory
# setwd("~/Desktop/morph_pred/")
setwd("~/academia/research/in_progress/morph_pred")


# read data
df_stress <- read_csv("./mySources/data/clean/stressBIN10iaClean.csv") %>%
  filter(., corr == 1,
            group %in% c('lb', 'la', 'ss'), 
            !participant %in% c('L01', 'L02', 'L03', 'L04', 'L05', 
                                 'L06', 'L07', 'L08', 'L09', 'L10', 
                                 'L15', 'L20', 'L21', 'L22', 'L23', 
                                 'L26', 'L30', 'L31', 'L33', 'LA04', 
                                 'LA06', 'LA07', 'LA14'))




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





















df_stress_50 <- read_csv("./mySources/data/clean/stressBIN50iaClean.csv") %>%
  select(., participant, group, target, condition, coda, bin, binAdj, 
            binTonsetAlign, binTsuffixAlign, targetCount, distractorCount, 
            targetProp, distractorProp, eLog, wts, corr) %>%
  mutate(., group = factor(group, levels = c("ss", "la", "lb"))) %>%
  filter(., corr == 1,  
            group %in% c('lb', 'la', 'ss'), 
            !participant %in% c('L01', 'L02', 'L03', 'L04', 'L05', 
                                 'L06', 'L07', 'L08', 'L09', 'L10', 
                                 'L15', 'L20', 'L21', 'L22', 'L23', 
                                 'L26', 'L30', 'L31', 'L33', 'LA04', 
                                 'LA06', 'LA07', 'LA14'), 
            binTonsetAlign >= 24 & binTonsetAlign <= 49, 
            !target %in% c('cambia', 'cambio'))




condition_names <- c(
                    `0` = 'CV', 
                    `1` = 'CVC'
                    )

(df_stress_50 %>%
  na.omit(.) %>% 
  ggplot(., aes(x = binTonsetAlign, y = targetProp, color = group, shape = group)) + 
  facet_grid(coda ~ ., labeller = as_labeller(condition_names)) + 
  geom_hline(yintercept = 0.5, color = 'white', size = 2) + 
  geom_vline(xintercept = 38, color = 'grey40', lty = 3) + 
  geom_vline(xintercept = 32, color = 'grey40', lty = 3) + 
  stat_summary(fun.data = mean_cl_boot, geom = 'pointrange',  size = 0.75, 
               fun.args = list(conf.int = .95, B = 1000)) +
  stat_summary(fun.y = mean, geom = 'point', color = 'white', alpha = 0.3, size = 1.75) + 
  scale_shape_manual(name = "", values = 17:15, labels = c("SS", "LA", "LB")) + 
  scale_color_brewer(palette = 'Set1', name = "", labels = c("SS", "LA", "LB")) + 
  scale_x_continuous(breaks = c(23, 28, 33, 38, 43, 48), labels = c("-750", "-500", "-250", "0", "250", "500")) + 
  labs(y = 'Proportion of target fixations', 
       x = 'Time relative to target syllable offset (ms)', 
       caption = "Mean +/- 95% CI") +
  coord_cartesian(ylim = c(0, 1)) + 
  annotate("text", x = 38.35, y = 0.02, label = 'Target syllable offset', 
           angle = 90, size = 3, hjust = 0) + 
  annotate("text", x = 32.35, y = 0.02, label = 'Mean target word onset', 
           angle = 90, size = 3, hjust = 0) + 
  theme_grey(base_size = 15, base_family = "Times") -> stressP1)

# ggsave('stressP1.png', plot = stressP1, dpi = 600, device = "png", 
#           path = "./mySources/figs/stress/s1_beg_adv_nat/eye_track", 
#           height = 6.75, width = 8.15, unit = "in")




























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
stress_subset_0 <- df_short %>% filter(., binTsuffixAlign == 147, !target %in% c('cambia', 'cambio'))

# Quick and dirty mean of target fixations as a function of 
# group and condition (stressed, unstressed *1st syllable*)
stress_subset_0 %>% 
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
  do(tidy(t.test(.$meanFix, alternative = "greater", mu = 0.5, conf.level = 0.95)))

# Convert pvalues from scientific notation
stress_ttest$p.value <- format(stress_ttest$p.value, scientific = F)
stress_ttest$sig <- "N.S."
stress_ttest[stress_ttest$p.value <= 0.05/6, 'sig'] <- "*"

# Print results
print(as.data.frame(stress_ttest[, c(1:7, 11)]))

# group  condition  estimate   statistic       p.value parameter  conf.low  sig
#    ss   stressed 0.6322078  2.79896749 0.00537726129        21 0.5509294    *
#    ss unstressed 0.7235390  4.76754875 0.00005197742        21 0.6428574    *
#    la   stressed 0.5022928  0.05122216 0.47976999874        26 0.4259471 N.S.
#    la unstressed 0.5841931  2.12785440 0.02149710485        26 0.5167067 N.S.
#    lb   stressed 0.5009398  0.01595298 0.49372369885        18 0.3987797 N.S.
#    lb unstressed 0.4880326 -0.21468852 0.58378816118        18 0.3913704 N.S.



# LB do not predict above chance 
# nor do the LA in the stressed condition

# We will plot the models 
# This will be almost exactly the same as 
# the previous plot, but it will use the 
# confidence interval from the test we 
# actually conducted

(stress_ttest %>%
  ungroup(.) %>%
  mutate(., group = factor(group, levels = c("lb", "la", "ss")), 
            condition = factor(condition, levels = c("stressed", "unstressed"))) %>% 
  ggplot(., aes(x = group, y = estimate, color = condition, 
                group = interaction(group, condition), dodge = condition)) +
    geom_hline(yintercept = 0.5, lty = 3) + 
    geom_linerange(aes(ymin = conf.low, ymax = estimate), color = 'grey40',
                       position = position_dodge(width = 0.75), size = 1) +
    geom_point(position = position_dodge(width = 0.75), size = 4) +
    # ylim(0, 1.0) +
    scale_color_brewer(palette = "Set1", name = '', labels = c('Paroxytone', 'Oxytone')) + 
    scale_x_discrete(labels = c('LB', 'LA', 'SS')) + 
    labs(title = 'Mean fixations at target syllable offset', 
         y = 'Target fixations', x = '', caption = 'Mean and lower-bound 95% CI') + 
    theme_bw(base_size = 12, base_family = 'Times') + 
    theme(legend.position = c(0.15, 0.9), 
          legend.background = element_rect(fill = NA)) -> stressTargetFixMOD1)

# ggsave('stressTargetFixMOD1.png', plot = stressTargetFixMOD1, dpi = 600, device = "png", 
#         path = "./mySources/figs/stress/s1_beg_adv_nat/eye_track", 
#         height = 4, width = 6, unit = 'in')






















## @knitr stressGroupCompare


# Load wm data and combine with stress_subset_0 (proportion data)
# in order to add working memory as a covariate
wm_df <- read_csv("./mySources/data/raw/wm_all.csv") %>% 
         filter(., !(group %in% c("HS", "IN")), 
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
 mutate(., group = factor(group, levels = c("ss", "la", "lb")), 
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
# 9697.2 9770.5 -4833.6   9667.2 17.116      2   0.000192  prop_0_mod_group ***
# 9699.1 9777.3 -4833.6   9667.1 0.0504      1     0.8223  prop_0_mod_cond
# 9691.6 9769.8 -4829.8   9659.6 7.5734      1   0.005924  prop_0_mod_coda  **
# 9693.1 9781.1 -4828.6   9657.1 2.4681      2     0.2911  group x coda
# 9692.3 9785.2 -4827.1   9654.3 5.3049      3     0.1508  group x condition
# 9693.7 9781.7 -4828.9   9657.7 1.8708      2     0.3924  cond x coda
# 9695.3 9812.7 -4823.7   9647.3 12.252      8     0.1403  prop_0_mod_full


stress_subset_0_prop$group <- factor(stress_subset_0_prop$group, levels = c("ss", "la",  "lb"))

prop_0_mod_final <- glmer(cbind(targetCount, 10 - targetCount) ~ 1 + group + codaSum + conditionSum + 
                    (1 | participant) + (1 | target) +  
                    (1 + conditionSum*codaSum | participant),
                    data = stress_subset_0_prop, family = 'binomial', 
                    control = glmerControl(optimizer = 'bobyqa'))

MuMIn::r.squaredGLMM(prop_0_mod_final)

# summary(prop_0_mod_final)
# confint(prop_0_mod_final, method = "Wald")

# Fixed effects:
#             Estimate  2.5%  97.5% Std. Error z value Pr(>|z|)    
# (Intercept)   1.6347  1.08   2.19     0.2851   5.735 9.77e-09 ***
# groupla      -1.2489 -1.92  -0.57     0.3446  -3.624  0.00029 ***
# grouplb      -1.5237 -2.28  -0.77     0.3869  -3.939 8.20e-05 ***
# codaSum1     -0.5303 -0.90  -0.16     0.1899  -2.793  0.00523 ** 


# Relevel to test lb vs la 
stress_subset_0_prop$group <- factor(stress_subset_0_prop$group, levels = c("lb", "la",  "ss"))

summary(glmer(cbind(targetCount, 10 - targetCount) ~ 1 + group + codaSum + 
       (1 | participant) + (1 | target) + 
       (1 + conditionSum*codaSum | participant),
       data = stress_subset_0_prop, family = 'binomial', 
       control = glmerControl(optimizer = 'bobyqa')))

# Fixed effects:
#             Estimate Std. Error z value Pr(>|z|)    
# (Intercept)   0.1110     0.3059   0.363  0.71664    
# groupla       0.2747     0.3714   0.740  0.45951    
# groupss       1.5237     0.3868   3.939 8.19e-05 ***
# codaSum1     -0.5303     0.1899  -2.793  0.00523 ** 




et_ci <- confint(prop_0_mod_final, method = "Wald", level = 0.99) %>% 
  as.data.frame(.) %>% 
  slice(., 13:16) %>% 
  rename(., ciLow = `0.5 %`, ciHi = `99.5 %`)

stressFixModP0 <- broom::tidy(prop_0_mod_final) %>% slice(1:4) %>% 
  cbind(., et_ci) %>% 
  mutate(., term = recode(term, `(Intercept)` = '(Intercept)', 
                                codaSum1 = 'Syllable\nstructure', 
                                groupla = 'LA', 
                                grouplb = 'LB'), 
            term = factor(term, levels = c('LB', 
                                           'LA', 
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
#         path = "./mySources/figs/stress/s1_beg_adv_nat/eye_track", 
#         height = 2.5, width = 7, unit = 'in')






# Calculate mean target fixation as a function of group, condition, 
# for each participant. We will plot the mean and calculate the 
# bootstrapped 95% confidence interval and plot it all. 

stressFixModsP1 <- stress_subset_0_prop %>% 
  mutate(., group = factor(group, levels = c("lb", "la",  "ss"))) %>%
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
    scale_x_discrete(labels = c('LB', 'LA', 'SS')) + 
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
  mutate(., group = factor(group, levels = c("lb", "la",  "ss"))) %>%
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
    scale_x_discrete(labels = c('LB', 'LA', 'SS')) + 
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

ggsave('stressGatingAll.png', 
       plot = stress_eyet_plot, dpi = 600, device = "png", 
      path = "./mySources/figs/stress/s1_beg_adv_nat/eye_track", 
       height = 3.5, width = 8.5, units = 'in')













































# LANDMARK ANALYSES

glimpse(df_stress)

# calculate proportion of target fixations by sub, by target, by stress cond, by group
# for each landmark in the time course

df_timecourse <- df_stress %>% 
  filter(., bin <= 400) %>% 
  select(., participant, group, bin, BIN_START_TIME, BIN_END_TIME, target, 
            startsentence:target, condition, coda, targetCount:targetProp, eLog) %>% 
  mutate(., bin = bin * 10)

# ggplot(df_timecourse, aes(x = bin, y = eLog, color = group)) + 
#   stat_summary(fun.data = mean_se, geom = 'pointrange')


range(df_timecourse$bin)

df_timecourse_startsentence <- df_timecourse %>% 
  filter(., startsentence + 1540 >= BIN_START_TIME, startsentence + 1540 <= BIN_END_TIME) %>% 
  select(., participant, group, target, condition:eLog) %>% 
  mutate(., landmark = 'start_sentence')

df_timecourse_word2_c1v1 <- df_timecourse %>% 
  filter(., word2_c1v1 + 1540 >= BIN_START_TIME, word2_c1v1 + 1540 <= BIN_END_TIME) %>% 
  select(., participant, group, target, condition:eLog) %>% 
  mutate(., landmark = 'word2_c1v1')


df_timecourse_word3_20msafterv1 <- df_timecourse %>% 
  filter(., word3_20msafterv1 + 1500 >= BIN_START_TIME, word3_20msafterv1 + 1500 <= BIN_END_TIME) %>% 
  select(., participant, group, target, condition:eLog) %>% 
  mutate(., landmark = 'word3_20msafterv1')


df_timecourse_word3_c1v1 <- df_timecourse %>% 
  filter(., word3_c1v1 + 1540 >= BIN_START_TIME, word3_c1v1 + 1540 <= BIN_END_TIME) %>% 
  select(., participant, group, target, condition:eLog) %>% 
  mutate(., landmark = 'word3_c1v1')


df_timecourse_word3_c2 <- df_timecourse %>% 
  filter(., word3_c2 + 1540 >= BIN_START_TIME, word3_c2 + 1540 <= BIN_END_TIME) %>% 
  select(., participant, group, target, condition:eLog) %>% 
  mutate(., landmark = 'word3_c2')

df_timecourse_word3_c3 <- df_timecourse %>% 
  filter(., word3_c3 + 1550 >= BIN_START_TIME, word3_c3 + 1550 <= BIN_END_TIME)  %>% 
  select(., participant, group, target, condition:eLog) %>% 
  mutate(., landmark = 'word3_c3')

df_timecourse_word3_suffix <- df_timecourse %>% 
  filter(., word3_suffix + 1540 >= BIN_START_TIME, word3_suffix + 1540 <= BIN_END_TIME) %>% 
  select(., participant, group, target, condition:eLog) %>% 
  mutate(., landmark = 'word3_suffix')

df_timecourse_word4_c1v1 <- df_timecourse %>% 
  filter(., word4_c1v1 + 1540 >= BIN_START_TIME, word4_c1v1 + 1540 <= BIN_END_TIME) %>% 
  select(., participant, group, target, condition:eLog) %>% 
  mutate(., landmark = 'word4_c1v1')

df_timecourse_word5 <- df_timecourse %>% 
  filter(., word5 + 1540 >= BIN_START_TIME, word5 + 1540 <= BIN_END_TIME) %>% 
  select(., participant, group, target, condition:eLog) %>% 
  mutate(., landmark = 'word5')

df_timecourse_end_sentence <- df_timecourse %>% 
  filter(., end_sentence + 1540 >= BIN_START_TIME, end_sentence + 1540 <= BIN_END_TIME)  %>% 
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
            landmark %in% c('word3_c1v1', 'word3_20msafterv1', 'word3_c2', 
                            'word3_suffix', 'word4_c1v1')) %>%
  group_by(., participant, group, condition, landmark) %>% 
  summarize(., target_fix = mean(targetProp)) %>% 
  ggplot(., aes(x = landmark, y = target_fix, shape = group, 
                dodge = group, fill = group)) + 
    geom_hline(yintercept = 0.5, color = 'black', lty = 3) + 
    stat_summary(fun.data = mean_cl_boot, geom = 'pointrange', 
                 position = position_dodge(0.5), size = 1.1) + 
    scale_shape_manual(values = c(21, 22, 23)) + 
    scale_fill_brewer(palette = "Set1") + 
    scale_x_discrete(name = 'Time point', 
                     labels = c('Target word\nonset', '20ms after\nV1', 
                                'Target syllable\noffset', 'Target suffix', 
                                'Following\nword')) + 
    coord_cartesian(ylim = c(0.3, 1)) + 
    theme_minimal()

df_landmarks %>% 
  na.omit(.) %>% 
  filter(., coda == 1, 
            landmark %in% c('word3_c1v1', 'word3_20msafterv1', 'word3_c2', 'word3_c3', 'word3_suffix', 'word4_c1v1')) %>%
  group_by(., participant, target, group, condition, landmark) %>% 
  summarize(., target_fix = mean(targetProp)) %>% 
  ggplot(., aes(x = landmark, y = target_fix, shape = group, 
                dodge = group, fill = group)) + 
    geom_hline(yintercept = 0.5, color = 'black', lty = 3) + 
    stat_summary(fun.data = mean_cl_boot, geom = 'pointrange', 
                 position = position_dodge(0.5), size = 1.1) + 
    scale_shape_manual(values = c(21, 22, 23)) + 
    scale_fill_brewer(palette = "Set1") + 
    scale_x_discrete(name = 'Time point', 
                     labels = c('Target word\nonset', '20ms after\nV1', 'Coda onset', 
                                'Target syllable\noffset', 'Target suffix', 'Following\nword')) + 
    coord_cartesian(ylim = c(0.3, 1)) + 
    theme_minimal()

# 


















































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
stress_gc_subset$group <- factor(stress_gc_subset$group, levels = c("ss", "la", "lb"))

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
gc_mod_base    <- readRDS('./mySources/models/stress/s1_beg_adv_nat/eye_track/gc_mod_base.rds')
gc_mod_group_0 <- readRDS('./mySources/models/stress/s1_beg_adv_nat/eye_track/gc_mod_group_0.rds')
gc_mod_group_1 <- readRDS('./mySources/models/stress/s1_beg_adv_nat/eye_track/gc_mod_group_1.rds')
gc_mod_group_2 <- readRDS('./mySources/models/stress/s1_beg_adv_nat/eye_track/gc_mod_group_2.rds')
gc_mod_group_3 <- readRDS('./mySources/models/stress/s1_beg_adv_nat/eye_track/gc_mod_group_3.rds')
gc_mod_cond_0  <- readRDS('./mySources/models/stress/s1_beg_adv_nat/eye_track/gc_mod_cond_0.rds')
gc_mod_cond_1  <- readRDS('./mySources/models/stress/s1_beg_adv_nat/eye_track/gc_mod_cond_1.rds')
gc_mod_cond_2  <- readRDS('./mySources/models/stress/s1_beg_adv_nat/eye_track/gc_mod_cond_2.rds')
gc_mod_cond_3  <- readRDS('./mySources/models/stress/s1_beg_adv_nat/eye_track/gc_mod_cond_3.rds')
gc_mod_full    <- readRDS('./mySources/models/stress/s1_beg_adv_nat/eye_track/gc_mod_full.rds')


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
  saveRDS(gc_mod_base, file = "./mySources/models/stress/s1_beg_adv_nat/eye_track/gc_mod_base.rds", compress = 'xz')
}




# Add group effect on intercept 
if(T){
  gc_mod_group_0 <- lmer(eLog ~ (ot1+ot2) + group + 
                    (1 + (ot1+ot2) | participant) + 
                    (1 + (ot1+ot2) | participant:codaSum), 
                    control = lmerControl(optimizer = 'bobyqa'), 
                    data = stress_gc_subset, weights = 1/wts, REML = F)
  saveRDS(gc_mod_group_0, file = "./mySources/models/stress/s1_beg_adv_nat/eye_track/gc_mod_group_0.rds", compress = 'xz')
}

# Add group effect on slope
if(T){
 gc_mod_group_1 <- lmer(eLog ~ (ot1+ot2) + group + 
                   ot1:group + 
                   (1 + (ot1+ot2) | participant) + 
                   (1 + (ot1+ot2) | participant:codaSum), 
                   control = lmerControl(optimizer = 'bobyqa'), 
                   data = stress_gc_subset, weights = 1/wts, REML = F)
  saveRDS(gc_mod_group_1, file = "./mySources/models/stress/s1_beg_adv_nat/eye_track/gc_mod_group_1.rds", compress = 'xz')
}

# Add group effect on quadratic poly 
if(T){
  gc_mod_group_2 <- lmer(eLog ~ (ot1+ot2) + group + 
                    ot1:group + ot2:group + 
                    (1 + (ot1+ot2) | participant) + 
                    (1 + (ot1+ot2) | participant:codaSum), 
                    control = lmerControl(optimizer = 'bobyqa'), 
                    data = stress_gc_subset, weights = 1/wts, REML = F)
  saveRDS(gc_mod_group_2, file = "./mySources/models/stress/s1_beg_adv_nat/eye_track/gc_mod_group_2.rds", compress = 'xz')
}

anova(gc_mod_base, 
      gc_mod_group_0, 
      gc_mod_group_1, 
      gc_mod_group_2, test = 'Chisq')




# full mod
if(T){
gc_mod_full_0 <- lmer(eLog ~ (ot1+ot2) * group * codaSum + 
               (1 + (ot1+ot2) | participant) + 
               # (1 + (ot1+ot2) | participant:codaSum) +
               (1 + (ot1+ot2) | participant:conditionSum), 
               control = lmerControl(optimizer = 'bobyqa'), 
               data = stress_gc_subset, weights = 1/wts, REML = F)
  saveRDS(gc_mod_full_0, file = "./mySources/models/stress/s1_beg_adv_nat/eye_track/gc_mod_ful_0.rds", compress = 'xz')
}

if(T){
gc_mod_full_1 <- lmer(eLog ~ (ot1+ot2) * group * codaSum + conditionSum + 
               (1 + (ot1+ot2) | participant) + 
               # (1 + (ot1+ot2) | participant:codaSum) +
               (1 + (ot1+ot2) | participant:conditionSum), 
               control = lmerControl(optimizer = 'bobyqa'), 
               data = stress_gc_subset, weights = 1/wts, REML = F)
  saveRDS(gc_mod_full_1, file = "./mySources/models/stress/s1_beg_adv_nat/eye_track/gc_mod_full_1.rds", compress = 'xz')
}

anova(gc_mod_full_0, gc_mod_full_1) # no main effect of cond
# ..1    32 60048 60281 -29992    59984 3.386      1    0.06575 .



# summary(gc_mod_full_0)

# Random effects:
#  Groups                   Name        Variance Std.Dev. Corr       
#  participant:conditionSum ot1         3.52654  1.8779              
#                           ot2         0.77497  0.8803   0.18       
#  participant:codaSum      (Intercept) 0.19775  0.4447              
#                           ot1         1.32949  1.1530   -0.35      
#                           ot2         0.08958  0.2993   -0.45  0.59
#  participant              (Intercept) 0.24703  0.4970              
#                           ot1         0.25367  0.5037   0.85       
#                           ot2         0.08633  0.2938   0.80  0.37 
#  Residual                             9.04131  3.0069              
# 
# Fixed effects:
#                        Estimate Std. Error         df t value Pr(>|t|)    
# (Intercept)           8.897e-01  1.507e-01  8.500e+01   5.906 6.94e-08 ***
# ot1                   3.295e+00  3.774e-01  8.600e+01   8.730 1.80e-13 ***
# ot2                  -2.040e-01  2.332e-01  6.800e+01  -0.875  0.38466    

# groupla              -3.743e-01  2.030e-01  8.500e+01  -1.844  0.06871 .  
# ot1:groupla          -8.039e-01  5.097e-01  8.600e+01  -1.577  0.11844    
# ot2:groupla           1.073e+00  3.154e-01  6.900e+01   3.402  0.00111 ** 

# grouplb              -5.739e-01  2.218e-01  8.600e+01  -2.587  0.01136 *  
# ot1:grouplb          -1.015e+00  5.581e-01  8.700e+01  -1.819  0.07237 .  
# ot2:grouplb           1.045e+00  3.458e-01  7.100e+01   3.021  0.00350 ** 

# codaSum1              7.244e-02  5.109e-02  1.060e+04   1.418  0.15629    
# ot1:codaSum1          2.509e-01  1.787e-01  9.833e+03   1.404  0.16035    
# ot2:codaSum1          4.669e-02  1.717e-01  6.695e+03   0.272  0.78567    

# groupla:codaSum1     -1.757e-01  6.868e-02  1.068e+04  -2.559  0.01051 *  
# ot1:groupla:codaSum1 -3.413e-01  2.405e-01  9.909e+03  -1.419  0.15579    
# ot2:groupla:codaSum1  1.036e-01  2.335e-01  6.815e+03   0.444  0.65712    

# grouplb:codaSum1     -1.504e-02  7.660e-02  1.060e+04  -0.196  0.84436    
# ot1:grouplb:codaSum1  1.384e-01  2.677e-01  9.608e+03   0.517  0.60529    
# ot2:grouplb:codaSum1  5.830e-02  2.586e-01  7.001e+03   0.225  0.82166    














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
  geom_smooth(method = 'gam', formula = y ~ poly(x, 2), se = F, show.legend = FALSE) + 
  # stat_summary(aes(y = GCA_Full, color = group), fun.y = mean, geom = 'line', size = 0.4) + 
  stat_summary(fun.data = mean_cl_boot, geom = 'pointrange',  size = 0.75, 
               fun.args = list(conf.int = .95, B = 1000)) +
  stat_summary(fun.y = mean, geom = 'point', color = 'white', alpha = 0.3, size = 1.75) + 
  scale_shape_manual(name = "", values = 17:15, labels = c("SS", "LA", "LB")) + 
  labs(x = "Time relative to target syllable offset (ms)", y = "Fixation empirical logit", 
       caption = "Mean +/- 95% CI") + 
  scale_color_brewer(palette = "Set1", name = "", guide = 'legend', 
                     labels = c("SS", "LA", "LB")) + 
  scale_x_continuous(breaks = c(28, 33, 38, 43, 48), labels = c("-500", "-250", "0", "250", "500")) + 
  theme_grey(base_size = 15, base_family = "Times New Roman"))

ggsave('stressP2.png', plot = gca_full, dpi = 600, device = "png", 
          path = "./mySources/figs/stress/s1_beg_adv_nat/eye_track", 
          height = 4.5, width = 8.5, unit = "in")














