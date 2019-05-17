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
            group %in% c('hs', 'la', 'ss'), 
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
  mutate(., group = factor(group, levels = c("ss", "la", "hs"))) %>%
  filter(., corr == 1,  
            group %in% c('hs', 'la', 'ss'), 
            !participant %in% c('L01', 'L02', 'L03', 'L04', 'L05', 
                                 'L06', 'L07', 'L08', 'L09', 'L10', 
                                 'L15', 'L20', 'L21', 'L22', 'L23', 
                                 'L26', 'L30', 'L31', 'L33', 'LA04', 
                                 'LA06', 'LA07', 'LA14'), 
            binTonsetAlign >= 24 & binTonsetAlign <= 49, 
            !target %in% c('cambia', 'cambio'))


df_stress_50[df_stress_50$group == 'hs', 'binTonsetAlign'] <- df_stress_50[df_stress_50$group == 'hs', 'binTonsetAlign'] - 1


condition_names <- c(
                    `stressed` = 'Paroxytone', 
                    `unstressed` = 'Oxytone',
                    `0` = 'CV', 
                    `1` = 'CVC'
                    )

(df_stress_50 %>%
  na.omit(.) %>% 
  filter(., binTonsetAlign >= 25, binTonsetAlign <= 46) %>%
  ggplot(., aes(x = binTonsetAlign, y = targetProp, color = group, shape = group)) + 
  facet_grid(. ~ coda, labeller = as_labeller(condition_names)) + 
  geom_hline(yintercept = 0.5, color = 'white', size = 2) + 
  geom_vline(xintercept = 38, color = 'grey40', lty = 3) + 
  geom_vline(xintercept = 32, color = 'grey40', lty = 3) + 
  stat_summary(fun.data = mean_cl_boot, geom = 'pointrange',  size = 0.75, 
               fun.args = list(conf.int = .95, B = 1000)) +
  stat_summary(fun.y = mean, geom = 'point', color = 'white', alpha = 0.3, size = 1.75) + 
  scale_shape_manual(name = "", values = 17:15, labels = c("SS", "LA", "HS")) + 
  scale_color_brewer(palette = 'Set1', name = "", labels = c("SS", "LA", "HS")) + 
  scale_x_continuous(breaks = c(23, 28, 33, 38, 43), labels = c("-750", "-500", "-250", "0", "250")) + 
  labs(y = 'Proportion of target fixations', 
       x = 'Time relative to target syllable offset (ms)', 
       caption = "Mean +/- 95% CI") +
  coord_cartesian(ylim = c(0, 1)) + 
  annotate("text", x = 38.35, y = 0.02, label = 'Target syllable offset', 
           angle = 90, size = 3, hjust = 0) + 
  annotate("text", x = 32.35, y = 0.02, label = 'Mean target word onset', 
           angle = 90, size = 3, hjust = 0) + 
  theme_grey(base_size = 12, base_family = "Times") -> stressP1)

# ggsave('stressP1.png', plot = stressP1, dpi = 600, device = "png", 
#           path = "./mySources/figs/stress/s2_adv_her_nat/eye_track", 
#           height = 4, width = 9, unit = 'in')




























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

df_short_temp <- df_short
df_short_temp[df_short_temp$group == 'hs' & df_short_temp$condition == 'stressed', 'binTsuffixAlign'] <- df_short_temp[df_short_temp$group == 'hs' & df_short_temp$condition == 'stressed', 'binTsuffixAlign'] - 10
df_short_temp[df_short_temp$group == 'hs' & df_short_temp$condition == 'unstressed', 'binTsuffixAlign'] <- df_short_temp[df_short_temp$group == 'hs' & df_short_temp$condition == 'unstressed', 'binTsuffixAlign'] - 8
stress_subset_0 <- df_short_temp %>% filter(., binTsuffixAlign == 147, !target %in% c('cambia', 'cambio'))




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
  do(tidy(t.test(.$meanFix, alternative = "greater", mu = 0.5, conf.level = 0.99)))

# Convert pvalues from scientific notation
stress_ttest$p.value <- format(stress_ttest$p.value, scientific = F)
stress_ttest$sig <- "N.S."
stress_ttest[stress_ttest$p.value <= 0.05/6, 'sig'] <- "*"

# Print results
print(as.data.frame(stress_ttest[, c(1:7, 11)]))

# group  condition  estimate  statistic        p.value parameter  conf.low  sig
#    hs   stressed 0.6354396 2.89648426 0.003866034457        25 0.5192360    *
#    hs unstressed 0.7092995 5.39546267 0.000006721561        25 0.6128978    *
#    la   stressed 0.5022928 0.05122216 0.479769998735        26 0.3913462 N.S.
#    la unstressed 0.5841931 2.12785440 0.021497104850        26 0.4861208 N.S.
#    ss   stressed 0.6322078 2.79896749 0.005377261285        21 0.5132880    *
#    ss unstressed 0.7235390 4.76754875 0.000051977419        21 0.6054925    *


# We will plot the models 
# This will be almost exactly the same as 
# the previous plot, but it will use the 
# confidence interval from the test we 
# actually conducted

(stress_ttest %>%
  ungroup(.) %>%
  mutate(., group = factor(group, levels = c("la", "hs", "ss")), 
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
    scale_x_discrete(labels = c('LA', 'HS', 'SS')) + 
    labs(title = 'Mean fixations at target syllable offset', 
         y = 'Target fixations', x = '', caption = 'Mean and lower-bound 99% CI') + 
    theme_minimal(base_size = 12, base_family = 'Times') -> stressTargetFixMOD1)

# ggsave('stressTargetFixMOD1.png', plot = stressTargetFixMOD1, dpi = 600, device = "png", 
#         path = "./mySources/figs/stress/s2_adv_her_nat/eye_track", 
#         height = 4, width = 7, unit = 'in')






















## @knitr stressGroupCompare


# Load wm data and combine with stress_subset_0 (proportion data)
# in order to add working memory as a covariate
wm_df <- read_csv("./mySources/data/raw/wm_all.csv") %>% 
         filter(., !(group %in% c("IN", "L")), 
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
 mutate(., group = factor(group, levels = c("ss", "la", "hs")), 
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
anova(prop_0_mod_coda, prop_0_mod_int2) # group condition interaction
anova(prop_0_mod_coda, prop_0_mod_int3) # condi x coda interaction
anova(prop_0_mod_coda, prop_0_mod_full) # no three way interaction


#    AIC    BIC  logLik deviance  Chisq Chi Df Pr(>Chisq) 
#  10462  10537 -5216.1    10432 10.404      2   0.005505  prop_0_mod_group ***
#  10464  10543 -5215.7    10432  0.684      1     0.4082  prop_0_mod_cond
#  10456  10536 -5212.0    10424 8.2164      1   0.004151  prop_0_mod_coda  **
#  10458  10548 -5210.9    10422 2.0495      2     0.3589  group x coda
#  10454  10548 -5207.8    10416 8.2255      3    0.04157  group x condition
#  10453  10543 -5208.7    10417 6.4792      2    0.03918  cond x coda
#  10462  10581 -5206.7    10414 10.457      8     0.2344  prop_0_mod_full


stress_subset_0_prop$group <- factor(stress_subset_0_prop$group, levels = c("ss", "la",  "hs"))

prop_0_mod_final <- glmer(cbind(targetCount, 10 - targetCount) ~ 1 + group + codaSum + group:codaSum + 
                    (1 | participant) + (1 | target) +  
                    (1 + conditionSum*codaSum | participant),
                    data = stress_subset_0_prop, family = 'binomial', 
                    control = glmerControl(optimizer = 'bobyqa'))

MuMIn::r.squaredGLMM(prop_0_mod_final)

# summary(prop_0_mod_final)
# confint(prop_0_mod_final, method = "Wald")

#                  Estimate Std. Error z value Pr(>|z|)    
# (Intercept)        1.7355     0.3317   5.233 1.67e-07 ***
# groupla           -1.4098     0.4237  -3.327 0.000878 ***
# grouphs           -0.3085     0.4288  -0.720 0.471828    
# codaSum1          -0.8238     0.2804  -2.938 0.003300 ** 


# Relevel to test lb vs la 
stress_subset_0_prop$group <- factor(stress_subset_0_prop$group, levels = c("hs", "la",  "ss"))

summary(glmer(cbind(targetCount, 10 - targetCount) ~ 1 + group + codaSum + group:codaSum + 
       (1 | participant) + (1 | target) + 
       (1 + conditionSum*codaSum | participant),
       data = stress_subset_0_prop, family = 'binomial', 
       control = glmerControl(optimizer = 'bobyqa')))

# Fixed effects:
#              Estimate Std. Error z value Pr(>|z|)   
# (Intercept)  1.4270     0.3098    4.606   4.11e-06 ***
# groupla     -1.1013     0.4077   -2.701     0.0069 ** 
# groupss      0.3085     0.4288    0.720     0.4718    
# codaSum1    -0.6111     0.2170   -2.816    0.00486 **




et_ci <- confint(prop_0_mod_final, method = "Wald", level = 0.99) %>% 
  as.data.frame(.) %>% 
  slice(., 13:16) %>% 
  rename(., ciLow = `0.5 %`, ciHi = `99.5 %`)

stressFixModP0 <- broom::tidy(prop_0_mod_final) %>% slice(1:4) %>% 
  cbind(., et_ci) %>% 
  mutate(., term = recode(term, `(Intercept)` = '(Intercept)', 
                                codaSum1 = 'Syllable\nstructure', 
                                groupla = 'LA', 
                                groupint = 'HS'), 
            term = factor(term, levels = c('HS', 
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
#         path = "./mySources/figs/stress/s2_adv_her_nat/eye_track", 
#         height = 2.5, width = 7, unit = 'in')





# stressFixModsHLS <- stress_subset_0_prop %>% 
#   mutate(., group = factor(group, levels = c("la", "hs",  "ss"))) %>%
#   group_by(., group, coda, participant) %>%
#   summarise(., meanFix = mean(targetProp)) %>%
#   ggplot(., aes(x = group, y = meanFix, 
#                 group = interaction(group, coda), dodge = coda, color = coda)) + 
#     geom_hline(yintercept = 0.5, lty = 3) + 
#     stat_summary(fun.data = mean_cl_boot, geom = 'pointrange', 
#                  position = position_dodge(width = 0.5), color = 'black', 
#                  size = 0.90, show.legend = FALSE, fun.args = list(conf.int = 0.95)) + 
#     stat_summary(fun.y = mean, geom = 'point', size = 2.75,
#                  position = position_dodge(width = 0.5), show.legend = TRUE) + 
#     labs(y = 'Target fixations', x = '', caption = 'Mean +/- 99% bootstrap CI.', 
#          title = "Mean fixations at target syllable offset") + 
#     coord_cartesian(ylim = c(0.4, 1)) + 
#     scale_x_discrete(labels = c('LA', 'HS', 'SS')) + 
#     scale_color_brewer(name = '', palette = 'Set1', labels = c("CV", "CVC")) + 
#     theme_minimal(base_size = 12, base_family = 'Times') 

# ggsave('stressFixModsHLS.png', plot = stressFixModsHLS, dpi = 600, device = "png", 
#         path = "./mySources/figs/stress/s2_adv_her_nat/eye_track", 
#         height = 4, width = 7, unit = 'in')











# Calculate mean target fixation as a function of group, condition, 
# for each participant. We will plot the mean and calculate the 
# bootstrapped 95% confidence interval and plot it all. 

stressFixModsP1 <- stress_subset_0_prop %>% 
  mutate(., group = factor(group, levels = c("int", "la",  "ss"))) %>%
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
    scale_x_discrete(labels = c('INT', 'LA', 'SS')) + 
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
  mutate(., group = factor(group, levels = c("int", "la",  "ss"))) %>%
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
    scale_x_discrete(labels = c('INT', 'LA', 'SS')) + 
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

# ggsave('stressGatingAll.png', 
#        plot = stress_eyet_plot, dpi = 600, device = "png", 
#       path = "./mySources/figs/stress/s3_adv_int_nat/eye_track", 
#        height = 3.5, width = 8.5, units = 'in')











































# LANDMARK ANALYSES



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
            landmark %in% c('word3_c1v1', 'word3_20msafterv1', 'word3_c2', 'word3_c3', 
                            'word3_suffix', 'word4_c1v1')) %>%
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
stress_gc_subset$group <- factor(stress_gc_subset$group, levels = c("ss", "la", "hs"))

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
gc_mod_full    <- readRDS('./mySources/models/stress/s3_adv_int_nat/eye_track/gc_mod_full.rds')


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
  saveRDS(gc_mod_base, file = "./mySources/models/stress/s2_adv_her_nat/eye_track/gc_mod_base.rds", compress = 'xz')
}




# Add group effect on intercept 
if(T){
  gc_mod_group_0 <- lmer(eLog ~ (ot1+ot2) + group + 
                    (1 + (ot1+ot2) | participant) + 
                    (1 + (ot1+ot2) | participant:codaSum), 
                    control = lmerControl(optimizer = 'bobyqa'), 
                    data = stress_gc_subset, weights = 1/wts, REML = F)
  saveRDS(gc_mod_group_0, file = "./mySources/models/stress/s2_adv_her_nat/eye_track/gc_mod_group_0.rds", compress = 'xz')
}

# Add group effect on slope
if(T){
 gc_mod_group_1 <- lmer(eLog ~ (ot1+ot2) + group + 
                   ot1:group + 
                   (1 + (ot1+ot2) | participant) + 
                   (1 + (ot1+ot2) | participant:codaSum), 
                   control = lmerControl(optimizer = 'bobyqa'), 
                   data = stress_gc_subset, weights = 1/wts, REML = F)
  saveRDS(gc_mod_group_1, file = "./mySources/models/stress/s2_adv_her_nat/eye_track/gc_mod_group_1.rds", compress = 'xz')
}

# Add group effect on quadratic poly 
if(T){
  gc_mod_group_2 <- lmer(eLog ~ (ot1+ot2) + group + 
                    ot1:group + ot2:group + 
                    (1 + (ot1+ot2) | participant) + 
                    (1 + (ot1+ot2) | participant:codaSum), 
                    control = lmerControl(optimizer = 'bobyqa'), 
                    data = stress_gc_subset, weights = 1/wts, REML = F)
  saveRDS(gc_mod_group_2, file = "./mySources/models/stress/s2_adv_her_nat/eye_track/gc_mod_group_2.rds", compress = 'xz')
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
  saveRDS(gc_mod_full_0, file = "./mySources/models/stress/s2_adv_her_nat/eye_track/gc_mod_ful_0.rds", compress = 'xz')
}

if(T){
gc_mod_full_1 <- lmer(eLog ~ (ot1+ot2) * group * codaSum + conditionSum + 
               (1 + (ot1+ot2) | participant) + 
               # (1 + (ot1+ot2) | participant:codaSum) +
               (1 + (ot1+ot2) | participant:conditionSum), 
               control = lmerControl(optimizer = 'bobyqa'), 
               data = stress_gc_subset, weights = 1/wts, REML = F)
  saveRDS(gc_mod_full_1, file = "./mySources/models/stress/s2_adv_her_nat/eye_track/gc_mod_full_1.rds", compress = 'xz')
}

anova(gc_mod_full_0, gc_mod_full_1) # no main effect of cond
# ..1    32 60048 60281 -29992    59984 3.386      1    0.06575 .



# summary(gc_mod_full_0)

# Fixed effects:
#                        Estimate Std. Error         df t value Pr(>|t|)    
# (Intercept)           9.099e-01  1.561e-01  7.400e+01   5.828 1.34e-07 ***
# ot1                   3.319e+00  3.647e-01  7.700e+01   9.100 7.59e-14 ***
# ot2                  -1.655e-01  2.365e-01  1.040e+02  -0.700  0.48553    
# ot3                  -1.977e-01  1.837e-01  1.210e+02  -1.076  0.28410    
# 
# groupla              -3.887e-01  2.104e-01  7.500e+01  -1.847  0.06865 .  
# ot1:groupla          -8.104e-01  4.925e-01  7.800e+01  -1.645  0.10392    
# ot2:groupla           1.043e+00  3.197e-01  1.060e+02   3.263  0.00149 ** 
# ot3:groupla           1.925e-01  2.467e-01  1.200e+02   0.780  0.43670    
# 
# grouphs              -2.366e-01  2.127e-01  7.500e+01  -1.112  0.26954    
# ot1:grouphs           1.409e-02  5.003e-01  8.000e+01   0.028  0.97760    
# ot2:grouphs           5.808e-01  3.236e-01  1.070e+02   1.795  0.07552 .  
# ot3:grouphs          -5.068e-01  2.560e-01  1.380e+02  -1.979  0.04977 *  
# 
# codaSum1              7.096e-02  5.080e-02  1.120e+04   1.397  0.16250    
# ot1:codaSum1          2.754e-01  1.776e-01  8.429e+03   1.550  0.12112    
# ot2:codaSum1          6.292e-02  1.716e-01  7.111e+03   0.367  0.71382    
# ot3:codaSum1          2.682e-01  1.701e-01  4.316e+03   1.577  0.11485    
# 
# groupla:codaSum1     -1.761e-01  6.848e-02  1.077e+04  -2.571  0.01014 *  
# ot1:groupla:codaSum1 -3.423e-01  2.396e-01  8.369e+03  -1.429  0.15300    
# ot2:groupla:codaSum1  8.477e-02  2.333e-01  7.069e+03   0.363  0.71631    
# ot3:groupla:codaSum1 -1.532e-01  2.298e-01  5.171e+03  -0.666  0.50513    
# 
# grouphs:codaSum1     -8.252e-02  7.105e-02  1.122e+04  -1.161  0.24550    
# ot1:grouphs:codaSum1  1.294e-01  2.511e-01  8.185e+03   0.515  0.60638    
# ot2:grouphs:codaSum1 -3.631e-02  2.374e-01  7.181e+03  -0.153  0.87846    
# ot3:grouphs:codaSum1 -2.498e-01  2.381e-01  4.970e+03  -1.049  0.29413    




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
  scale_shape_manual(name = "", values = 17:15, labels = c("SS", "LA", "HS")) + 
  labs(x = "Time relative to target syllable offset (ms)", y = "Fixation empirical logit", 
       caption = "Mean +/- 95% CI") + 
  scale_color_brewer(palette = "Set1", name = "", guide = 'legend', 
                     labels = c("SS", "LA", "HS")) + 
  scale_x_continuous(breaks = c(28, 33, 38, 43, 48), labels = c("-500", "-250", "0", "250", "500")) + 
  theme_grey(base_size = 15, base_family = "Times New Roman"))

# ggsave('stressP2.png', plot = gca_full, dpi = 600, device = "png", 
#           path = "./mySources/figs/stress/s2_adv_her_nat/eye_track", 
#           height = 4, width = 9, unit = "in")














