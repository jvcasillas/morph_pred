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
            !(participant %in% c('L01', 'L02', 'L03', 'L04', 'L05', 
                                 'L06', 'L07', 'L08', 'L09', 'L10', 
                                 'L15', 'L20', 'L21', 'L22', 'L23', 
                                 'L26', 'L30', 'L31', 'L33', 'LA04', 
                                 'LA06', 'LA07', 'LA14')))









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
            !(participant %in% c('L01', 'L02', 'L03', 'L04', 'L05', 
                                 'L06', 'L07', 'L08', 'L09', 'L10', 
                                 'L15', 'L20', 'L21', 'L22', 'L23', 
                                 'L26', 'L30', 'L31', 'L33', 'LA04', 
                                 'LA06', 'LA07', 'LA14')), 
            binTonsetAlign >= 5 & binTonsetAlign <= 45)




condition_names <- c(
                    `stressed` = "Paroxytone",
                    `unstressed` = "Oxytone"
                    )

(df_stress_50 %>%
  na.omit(.) %>% 
  # filter(., coda == 1) %>% 
  ggplot(., aes(x = binTonsetAlign, y = targetProp, color = group, shape = group)) + 
  facet_grid(. ~ condition, labeller = as_labeller(condition_names)) + 
  geom_hline(yintercept = 0.5, color = 'white', size = 2) + 
  geom_vline(xintercept = 20, color = 'grey60', lty = 3) + 
  geom_vline(xintercept = 24, color = 'grey60', lty = 3) + 
  stat_summary(fun.data = mean_cl_boot, geom = 'pointrange',  size = 0.5, 
               fun.args = list(conf.int = .95, B = 5000)) +
  stat_summary(fun.y = mean, geom = 'point', color = 'white', alpha = 0.3, size = 1.5) + 
  scale_shape_manual(name = "", values = 17:15, labels = c("SS", "LA", "LB")) + 
  scale_color_brewer(palette = 'Set1', name = "", labels = c("SS", "LA", "LB")) + 
  scale_x_continuous(breaks = c(10, 20, 30, 40), labels = c("-500", "0", "500", "1000")) + 
  labs(y = 'Proportion of target fixations', 
       x = 'Time after target onset (ms)', 
       caption = "Mean +/- 95% CI") +
  #coord_cartesian(ylim = c(0, 1)) + 
  theme_grey(base_size = 12, base_family = "Times") -> stressP1)

# ggsave('stressP1.png', plot = stressP1, dpi = 600, device = "png", 
          # path = "./mySources/figs/stress/s1_beg_adv_nat/eye_track", 
          # height = 3.25, width = 10.5, unit = "in")




















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
stress_subset_0 <- df_short %>% filter(., binTsuffixAlign == 147)

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
stress_ttest[stress_ttest$p.value <= 0.05, 'sig'] <- "*"

# Print results
print(as.data.frame(stress_ttest[, c(1:7, 11)]))

#group  condition  estimate    statistic       p.value parameter  conf.low  sig
#   la   stressed 0.4933201 -0.148803174 0.55857145203        26 0.4167535 N.S.
#   la unstressed 0.5839286  2.030959997 0.02630044260        26 0.5134446    *
#   lb   stressed 0.4996867 -0.005911574 0.50232585022        18 0.4077902 N.S.
#   lb unstressed 0.4827694 -0.311696637 0.62057374910        18 0.3869105 N.S.
#   ss   stressed 0.6167749  2.555903872 0.00920282205        21 0.5381571    *
#   ss unstressed 0.7287338  4.696397066 0.00006155838        21 0.6449265    *



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
    ylim(0, 1.0) +
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

# Are groups different from each other?
stress_subset_0_prop <- stress_subset_0 %>%
  na.omit(.) %>% 
  group_by(., group, condition, participant) %>%
  summarise(., meanFix = mean(targetProp)) # %>% 
#  ungroup(.) %>%
#  left_join(x = ., y = wm_df[, -1], by = 'participant') %>% 
#  mutate(., group = factor(group, levels = c("ss", "la", "lb")), 
#            condition = as.factor(condition))

# tested with WM, not significant, so leaving out in 
# order to keep more data

# Set condition factors to sum contrasts
#stress_subset_0_prop$condition <- C(stress_subset_0_prop$condition, sum)

prop_0_mod_0 <- lmer(meanFix ~ 1 + 
                    (1 | participant), 
                    data = stress_subset_0_prop, REML = F, control = lmerControl(optimizer = 'bobyqa'))

prop_0_mod_group <- lmer(meanFix ~ 1 + group + 
                        (1 | participant), 
                        data = stress_subset_0_prop, REML = F, control = lmerControl(optimizer = 'bobyqa'))

prop_0_mod_cond <- lmer(meanFix ~ 1 + group + condition + 
                        (1 | participant), 
                        data = stress_subset_0_prop, REML = F, control = lmerControl(optimizer = 'bobyqa'))

prop_0_mod_full <- lmer(meanFix ~ 1 + group * condition + 
                      (1 | participant), 
                       data = stress_subset_0_prop, REML = F, control = lmerControl(optimizer = 'bobyqa'))

# prop_0_mod_cov <- lmer(meanFix ~ 1 + group * condition + WM + 
#                       (1 | participant), 
#                        data = stress_subset_0_prop, REML = F, control = lmerControl(optimizer = 'bobyqa'))

# prop_0_mod_cov_full <- lmer(meanFix ~ 1 + group * condition * WM + 
#                       (1 | participant), 
#                        data = stress_subset_0_prop, REML = F, control = lmerControl(optimizer = 'bobyqa'))


# anova(prop_0_mod_0, prop_0_mod_group, prop_0_mod_cond, prop_0_mod_full, test = "Chisq")

#        Df     AIC     BIC  logLik deviance   Chisq Chi Df Pr(>Chisq)    
# object  3  1.1555  9.8935  2.4222  -4.8445                              
# ..1     5 -8.8770  5.6863  9.4385 -18.8770 14.0325      2  0.0008972 ***
# ..2     6 -9.9486  7.5273 10.9743 -21.9486  3.0717      1  0.0796670 .  
# ..3     8 -7.9025 15.3988 11.9512 -23.9025  1.9538      2  0.3764707    

# summary(prop_0_mod_group)

# Fixed effects:
#              Estimate Std. Error        df t value Pr(>|t|)    
# (Intercept)   0.53862    0.03072 136.00000  17.533  < 2e-16 ***
# grouplb      -0.04740    0.04780 136.00000  -0.992  0.32318    
# groupss       0.13413    0.04585 136.00000   2.926  0.00403 ** 

# Relevel to test lb vs la 
stress_subset_0_prop$group <- factor(stress_subset_0_prop$group, levels = c("lb", "la",  "ss"))

summary(lmer(meanFix ~ 1 + group + (1 | participant), data = stress_subset_0_prop, REML = F, control = lmerControl(optimizer = 'bobyqa')))
# Learner groups are not different from each other




# Calculate mean target fixation as a function of group, condition, 
# for each participant. We will plot the mean and calculate the 
# bootstrapped 95% confidence interval and plot it all. 
(stress_subset_0 %>%
  na.omit(.) %>%
  group_by(., group, condition, participant) %>%
  summarise(., meanFix = mean(targetProp)) %>% 
  ungroup(.) %>%
  mutate(., group = factor(group, levels = c("lb", "la",  "ss"))) %>%
  ggplot(., aes(x = group, y = meanFix, 
                dodge = condition, color = condition,
                group = interaction(group, condition))) +
    geom_hline(yintercept = 0.5, lty = 3) + 
    stat_summary(fun.data = mean_cl_boot, geom = 'errorbar', 
                 position = position_dodge(width = 0.5), 
                 width = 0.35, color = 'grey40') + 
    stat_summary(fun.y = mean, geom = 'point', size = 4,
                 position = position_dodge(width = 0.5)) + 
    coord_cartesian(ylim = c(0, 1)) + 
    labs(title = 'j', y = '', x = '', caption = 'Mean +/- 95% CI') + 
    scale_x_discrete(labels = c('LB', 'LA', 'SS')) + 
    scale_y_continuous(position = 'right') + 
    scale_color_brewer(palette = "Set1", name = '', guide = F) + 
    theme_bw(base_size = 12, base_family = 'Times') + 
    theme(plot.title = element_text(color = "white")) -> stressTargetFixMOD2)

stressFixMods <- plot_grid(stressTargetFixMOD1, stressTargetFixMOD2, ncol = 2)


# Looks good, save as .png file
# ggsave('stressFixMods.png', plot = stressFixMods, dpi = 600, device = "png", 
#          path = "./mySources/figs/stress/s1_beg_adv_nat/eye_track", 
#          height = 3.25, width = 8.5, unit = "in")

























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
stress_gc_subset <- filter(df_stress_50, binTonsetAlign >= 16 & 
                                         binTonsetAlign <= 45) %>% as.data.frame 


# - Readjust time course 
#    - now we will make the time course positive, starting at 1
#    - to do this we will add the lowest 'binAdj' value to each 
#      bin, plus 1 (to avoid starting at 0)

# lsrl_gc_subset$binGC <- lsrl_gc_subset$binAdj + 51
stress_gc_subset$binGC <- stress_gc_subset$binTonsetAlign - 15


# - Now we add higher order polynomials for analyses

t <- poly(min(stress_gc_subset$binGC):max(stress_gc_subset$binGC), 3)
stress_gc_subset[, paste('ot', 1:3, sep = "")] <- t[stress_gc_subset$binGC, 1:3]

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

stress_gc_subset$groupSum <- C(stress_gc_subset$group, sum)
# contrasts(stress_gc_subset$groupSum)

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


# Base model 
if(T){
  gc_mod_base <- lmer(eLog ~ (ot1+ot2+ot3) + 
                 ((ot1+ot2+ot3) | participant) + 
                 ((ot1+ot2+ot3) | participant:conditionSum),
                 ((ot1+ot2+ot3) | target),
                 control = lmerControl(optimizer = 'bobyqa'), 
                 data = stress_gc_subset, weights = 1/wts, REML = F)
  saveRDS(gc_mod_base, file = "./mySources/models/stress/s1_beg_adv_nat/eye_track/gc_mod_base.rds", compress = 'xz')
}


# Add group effect on intercept 
if(T){
  gc_mod_group_0 <- lmer(eLog ~ (ot1+ot2+ot3) + groupSum + 
                    ((ot1+ot2+ot3) | participant) + 
                    ((ot1+ot2+ot3) | participant:conditionSum),
                    ((ot1+ot2+ot3) | target),
                    control = lmerControl(optimizer = 'bobyqa'), 
                    data = stress_gc_subset, weights = 1/wts, REML = F)
  saveRDS(gc_mod_group_0, file = "./mySources/models/stress/s1_beg_adv_nat/eye_track/gc_mod_group_0.rds", compress = 'xz')
}

# Add group effect on slope
if(T){
 gc_mod_group_1 <- lmer(eLog ~ (ot1+ot2+ot3) + groupSum + 
                   ot1:groupSum + 
                   ((ot1+ot2+ot3) | participant) + 
                   ((ot1+ot2+ot3) | participant:conditionSum),
                   ((ot1+ot2+ot3) | target),
                   control = lmerControl(optimizer = 'bobyqa'), 
                   data = stress_gc_subset, weights = 1/wts, REML = F)
  saveRDS(gc_mod_group_1, file = "./mySources/models/stress/s1_beg_adv_nat/eye_track/gc_mod_group_1.rds", compress = 'xz')
}

# Add group effect on quadratic poly 
if(T){
  gc_mod_group_2 <- lmer(eLog ~ (ot1+ot2+ot3) + groupSum + 
                    ot1:groupSum + ot2:groupSum + 
                    ((ot1+ot2+ot3) | participant) + 
                    ((ot1+ot2+ot3) | participant:conditionSum),
                    ((ot1+ot2+ot3) | target),
                    control = lmerControl(optimizer = 'bobyqa'), 
                    data = stress_gc_subset, weights = 1/wts, REML = F)
  saveRDS(gc_mod_group_2, file = "./mySources/models/stress/s1_beg_adv_nat/eye_track/gc_mod_group_2.rds", compress = 'xz')
}

# Add group effect on cubic poly 
if(T){
  gc_mod_group_3 <- lmer(eLog ~ (ot1+ot2+ot3) + groupSum + 
                    ot1:groupSum + ot2:groupSum + ot3:groupSum + 
                    ((ot1+ot2+ot3) | participant) + 
                    ((ot1+ot2+ot3) | participant:conditionSum),
                    ((ot1+ot2+ot3) | target),
                    control = lmerControl(optimizer = 'bobyqa'), 
                    data = stress_gc_subset, weights = 1/wts, REML = F)
  saveRDS(gc_mod_group_3, file = "./mySources/models/stress/s1_beg_adv_nat/eye_track/gc_mod_group_3.rds", compress = 'xz')
}

# Add condition effect on intercept 
if(T){
  gc_mod_cond_0 <- lmer(eLog ~ (ot1+ot2+ot3) * groupSum + conditionSum + 
                   ((ot1+ot2+ot3) | participant) + 
                   ((ot1+ot2+ot3) | participant:conditionSum),
                   ((ot1+ot2+ot3) | target),
                   control = lmerControl(optimizer = 'bobyqa'), 
                   data = stress_gc_subset, weights = 1/wts, REML = F)
  saveRDS(gc_mod_cond_0, file = "./mySources/models/stress/s1_beg_adv_nat/eye_track/gc_mod_cond_0.rds", compress = 'xz')
}

# Add condition effect on slope 
if(T){
  gc_mod_cond_1 <- lmer(eLog ~ (ot1+ot2+ot3) * groupSum + conditionSum + 
                   ot1:conditionSum + 
                   ((ot1+ot2+ot3) | participant) + 
                   ((ot1+ot2+ot3) | participant:conditionSum),
                   ((ot1+ot2+ot3) | target),
                   control = lmerControl(optimizer = 'bobyqa'), 
                   data = stress_gc_subset, weights = 1/wts, REML = F)
  saveRDS(gc_mod_cond_1, file = "./mySources/models/stress/s1_beg_adv_nat/eye_track/gc_mod_cond_1.rds", compress = 'xz')
}

# Add condition effect on quadratic poly 
if(T){
  gc_mod_cond_2 <- lmer(eLog ~ (ot1+ot2+ot3) * groupSum + conditionSum + 
                   ot1:conditionSum + ot2:conditionSum + 
                   ((ot1+ot2+ot3) | participant) + 
                   ((ot1+ot2+ot3) | participant:conditionSum),
                   ((ot1+ot2+ot3) | target),
                   control = lmerControl(optimizer = 'bobyqa'), 
                   data = stress_gc_subset, weights = 1/wts, REML = F)
  saveRDS(gc_mod_cond_2, file = "./mySources/models/stress/s1_beg_adv_nat/eye_track/gc_mod_cond_2.rds", compress = 'xz')
}

# Add condition effect on cubic poly 
if(T){
  gc_mod_cond_3 <- lmer(eLog ~ (ot1+ot2+ot3) * groupSum + conditionSum + 
                   ot1:conditionSum + ot2:conditionSum + ot3:conditionSum + 
                   ((ot1+ot2+ot3) | participant) + 
                   ((ot1+ot2+ot3) | participant:conditionSum),
                   ((ot1+ot2+ot3) | target),
                   control = lmerControl(optimizer = 'bobyqa'), 
                   data = stress_gc_subset, weights = 1/wts, REML = F)
  saveRDS(gc_mod_cond_3, file = "./mySources/models/stress/s1_beg_adv_nat/eye_track/gc_mod_cond_3.rds", compress = 'xz')
}

# Include all interactions
if(T){
gc_mod_full <- glmer(cbind(targetCount, 50 - targetCount) ~ (ot1+ot2+ot3) * groupSum * conditionSum + 
               ((ot1+ot2+ot3) | participant) + 
               ((ot1+ot2+ot3) | participant:conditionSum) + 
               ((ot1+ot2+ot3) | target),  
               control = glmerControl(optimizer = 'bobyqa'), 
               data = stress_gc_subset, family = 'binomial', REML = F)
  saveRDS(gc_mod_full, file = "./mySources/models/stress/s1_beg_adv_nat/eye_track/gc_mod_full.rds", compress = 'xz')
}

# Model comparison 
# - subsequent models test three time terms: 
#    - the intercept (additive effects)
#    - the linear slope (ot1)
#    - the steepness of the quadratic curvature (ot2)

anova(gc_mod_base, 
      gc_mod_group_0, 
      gc_mod_group_1, 
      gc_mod_group_2, 
      gc_mod_group_3, 
      gc_mod_cond_0, 
      gc_mod_cond_1, 
      gc_mod_cond_2,
      gc_mod_cond_3, 
      gc_mod_full, test = 'Chisq')

# Df    AIC    BIC  logLik deviance   Chisq Chi Df Pr(>Chisq)    
# 25 242262 242478 -121106   242212                               # base model 
# 27 242263 242497 -121105   242209  2.2181      2   0.329865     # add group effect on intercept 
# 29 242267 242518 -121104   242209  0.5739      2   0.750561     # add group effect on slope
# 31 242264 242532 -121101   242202  7.3073      2   0.025896 *   # add group effect on quadratic poly 
# 33 242265 242551 -121100   242199  2.0534      2   0.358187     # add group effect on cubic poly 
# 34 242260 242554 -121096   242192  7.6926      1   0.005545 **  # add cond effect on intercept
# 35 242261 242564 -121095   242191  1.1893      1   0.275474     # add cond effect on slope 
# 36 242234 242546 -121081   242162 28.5858      1  8.964e-08 *** # add cond effect on quadratic poly
# 37 242236 242556 -121081   242162  0.2078      1   0.648520     # add cond effect on cubic poly 
# 45 242162 242552 -121036   242072 89.6219      8  5.549e-16 *** # full model 



# summary(gc_mod_full)
# Fixed effects:
#                                   Estimate Std. Error         df t value  Pr(>|t|)    
# (Intercept)                      5.118e-01  2.436e-02  8.900e+01  21.009   < 2e-16 ***
# ot1                              6.938e-01  1.183e-01  8.600e+01   5.863  8.12e-08 ***
# ot2                              5.388e-01  8.124e-02  9.200e+01   6.632  2.22e-09 ***
# groupla                          3.404e-02  2.532e-02  6.500e+01   1.344   0.18355    
# grouplb                          2.456e-02  3.129e-02  6.300e+01   0.785   0.43543    
# conditionunstressed              5.412e-02  2.335e-02  3.500e+01   2.318   0.02641 *  
# ot1:groupla                     -5.446e-01  1.201e-01  6.900e+01  -4.535  2.37e-05 ***
# ot1:grouplb                     -7.034e-01  1.467e-01  6.400e+01  -4.794  1.01e-05 ***
# ot2:groupla                     -1.173e-01  8.461e-02  8.300e+01  -1.387   0.16927    
# ot2:grouplb                     -3.481e-01  1.015e-01  7.200e+01  -3.430   0.00101 ** 
# ot1:conditionunstressed          2.866e-01  1.223e-01  4.000e+01   2.344   0.02411 *  
# ot2:conditionunstressed         -2.692e-01  8.849e-02  5.200e+01  -3.043   0.00368 ** 
# groupla:conditionunstressed     -6.961e-02  8.374e-03  7.252e+04  -8.313   < 2e-16 ***
# grouplb:conditionunstressed     -1.483e-01  9.160e-03  7.252e+04 -16.186   < 2e-16 ***
# ot1:groupla:conditionunstressed  3.680e-01  6.377e-02  7.253e+04   5.770  7.95e-09 ***
# ot1:grouplb:conditionunstressed  3.157e-01  6.976e-02  7.252e+04   4.526  6.02e-06 ***
# ot2:groupla:conditionunstressed  7.170e-02  6.377e-02  7.254e+04   1.124   0.26086    
# ot2:grouplb:conditionunstressed  3.227e-01  6.976e-02  7.254e+04   4.627  3.72e-06 ***




# create new df including the fitted model 
data.comp <- data.frame(na.omit(stress_gc_subset), 
                        GCA_Full = fitted(gc_mod_full))
# glimpse(data.comp)

condition_namesGCAMod <- c(
                    `ss` = "SS", 
                    `la` = "LA", 
                    `lb` = "LB"
                    )

data.comp %>% 
  ggplot(., aes(x = binGC, y = targetProp, color = condition)) + 
  facet_grid(. ~ group, labeller = as_labeller(condition_namesGCAMod)) + 
  stat_summary(fun.data = mean_se, geom = 'errorbar', 
               show.legend = FALSE, size = 0.1) +
  stat_summary(fun.y = mean, geom = 'point', size = 0.2) + 
  stat_summary(aes(y = GCA_Full, color = condition), fun.y = mean, geom = 'line', size = 0.4) + 
  labs(x = "Time after target onset (ms)", y = "Target fixations") + 
  scale_x_continuous(breaks = c(1, 45), labels = c("200", "1000")) + 
  scale_color_brewer(palette = "Set1", name = "", 
                     labels = c("Oxytone", "Paroxytone")) + 
  theme_bw(base_size = 12, base_family = "Times New Roman") -> stressGCAfullMod

# ggsave('stressGCAfullMod.png', plot = stressGCAfullMod, dpi = 600, device = "png", path = "./mySources/figs/stress/s1_beg_adv_nat/eye_track")
















# RUN GCA WITHOUT LB GROUP 

lsrl_gc_ss_la <- lsrl_gc_subset[lsrl_gc_subset$group != 'lb', ]

lsrl_gc_ss_la$group <- droplevels(lsrl_gc_ss_la$group)

lsrl_gc_ss_la$groupSum <- C(lsrl_gc_ss_la$group, sum)

contrasts(lsrl_gc_ss_la$groupSum)
contrasts(lsrl_gc_ss_la$conditionSum)

# load the models:
# gc_mod2_base    <- readRDS('./mySources/models/stress/s1_beg_adv_nat/eye_track/gc_mod2_base.rds')
# gc_mod2_group_0 <- readRDS('./mySources/models/stress/s1_beg_adv_nat/eye_track/gc_mod2_group_0.rds')
# gc_mod2_group_1 <- readRDS('./mySources/models/stress/s1_beg_adv_nat/eye_track/gc_mod2_group_1.rds')
# gc_mod2_group_2 <- readRDS('./mySources/models/stress/s1_beg_adv_nat/eye_track/gc_mod2_group_2.rds')
# gc_mod2_cond_0  <- readRDS('./mySources/models/stress/s1_beg_adv_nat/eye_track/gc_mod2_cond_0.rds')
# gc_mod2_cond_1  <- readRDS('./mySources/models/stress/s1_beg_adv_nat/eye_track/gc_mod2_cond_1.rds')
# gc_mod2_cond_2  <- readRDS('./mySources/models/stress/s1_beg_adv_nat/eye_track/gc_mod2_cond_2.rds')
# gc_mod2_full    <- readRDS('./mySources/models/stress/s1_beg_adv_nat/eye_track/gc_mod2_full.rds')


# Base model 
if(T){
  gc_mod2_base <- lmer(targetProp ~ (ot1+ot2) + 
                 ((ot1+ot2) | participant) + 
                 ((ot1+ot2) | target),
                 control = lmerControl(optimizer = 'bobyqa'), 
                 data = lsrl_gc_ss_la, REML = F)
  saveRDS(gc_mod2_base, file = "./mySources/models/stress/s1_beg_adv_nat/eye_track/gc_mod2_base.rds", compress = 'xz')
}


# Add group effect on intercept 
if(T){
  gc_mod2_group_0 <- lmer(targetProp ~ (ot1+ot2) + groupSum + 
                    ((ot1+ot2) | participant) + 
                    ((ot1+ot2) | target), 
                    control = lmerControl(optimizer = 'bobyqa'), 
                    data = lsrl_gc_ss_la, REML = F)
  saveRDS(gc_mod2_group_0, file = "./mySources/models/stress/s1_beg_adv_nat/eye_track/gc_mod2_group_0.rds", compress = 'xz')
}

# Add group effect on slope
if(T){
 gc_mod2_group_1 <- lmer(targetProp ~ (ot1+ot2) + groupSum + 
                   ot1:groupSum + 
                   ((ot1+ot2) | participant) + 
                   ((ot1+ot2) | target),
                   control = lmerControl(optimizer = 'bobyqa'), 
                   data = lsrl_gc_ss_la, REML = F)
  saveRDS(gc_mod2_group_1, file = "./mySources/models/stress/s1_beg_adv_nat/eye_track/gc_mod2_group_1.rds", compress = 'xz')
}

# Add group effect on quadratic poly 
if(T){
  gc_mod2_group_2 <- lmer(targetProp ~ (ot1+ot2) + groupSum + 
                    ot1:groupSum + ot2:groupSum + 
                    ((ot1+ot2) | participant) + 
                    ((ot1+ot2) | target),
                    control = lmerControl(optimizer = 'bobyqa'), 
                    data = lsrl_gc_ss_la, REML = F)
  saveRDS(gc_mod2_group_2, file = "./mySources/models/stress/s1_beg_adv_nat/eye_track/gc_mod2_group_2.rds", compress = 'xz')
}


# Add condition effect on intercept 
if(T){
  gc_mod2_cond_0 <- lmer(targetProp ~ (ot1+ot2) * groupSum + conditionSum + 
                   ((ot1+ot2) | participant) + 
                   ((ot1+ot2) | target),
                   control = lmerControl(optimizer = 'bobyqa'), 
                   data = lsrl_gc_ss_la, REML = F)
  saveRDS(gc_mod2_cond_0, file = "./mySources/models/stress/s1_beg_adv_nat/eye_track/gc_mod2_cond_0.rds", compress = 'xz')
}

# Add condition effect on slope 
if(T){
  gc_mod2_cond_1 <- lmer(targetProp ~ (ot1+ot2) * groupSum + conditionSum + 
                   ot1:conditionSum + 
                   ((ot1+ot2) | participant) + 
                   ((ot1+ot2) | target),
                   control = lmerControl(optimizer = 'bobyqa'), 
                   data = lsrl_gc_ss_la, REML = F)
  saveRDS(gc_mod2_cond_1, file = "./mySources/models/stress/s1_beg_adv_nat/eye_track/gc_mod2_cond_1.rds", compress = 'xz')
}

# Add condition effect on quadratic poly 
if(T){
  gc_mod2_cond_2 <- lmer(targetProp ~ (ot1+ot2) * groupSum + conditionSum + 
                   ot1:conditionSum + ot2:conditionSum + 
                   ((ot1+ot2) | participant) + 
                   ((ot1+ot2) | target),
                   control = lmerControl(optimizer = 'bobyqa'), 
                   data = lsrl_gc_ss_la, REML = F)
  saveRDS(gc_mod2_cond_2, file = "./mySources/models/stress/s1_beg_adv_nat/eye_track/gc_mod2_cond_2.rds", compress = 'xz')
}


# Include all interactions
if(T){
gc_mod2_full <- lmer(targetProp ~ (ot1+ot2) * groupSum * conditionSum + 
               ((ot1+ot2) | participant) + 
               ((ot1+ot2) | target),
               control = lmerControl(optimizer = 'bobyqa'), 
               data = lsrl_gc_ss_la, REML = F)
  saveRDS(gc_mod2_full, file = "./mySources/models/stress/s1_beg_adv_nat/eye_track/gc_mod2_full.rds", compress = 'xz')
}

# Model comparison 
# - subsequent models test three time terms: 
#    - the intercept (additive effects)
#    - the linear slope (ot1)
#    - the steepness of the quadratic curvature (ot2)

# anova(gc_mod2_base, 
#       gc_mod2_group_0, 
#       gc_mod2_group_1, 
#       gc_mod2_group_2, 
#       gc_mod2_cond_0, 
#       gc_mod2_cond_1, 
#       gc_mod2_cond_2, 
#       gc_mod2_full, test = 'Chisq')

#        Df   AIC   BIC logLik deviance    Chisq Chi Df Pr(>Chisq)    
# object 16 68773 68915 -34370    68741                                # base model 
# ..1    17 68775 68925 -34370    68741   0.0040      1   0.949859     # add group effect on intercept 
# ..2    18 68768 68928 -34366    68732   8.4043      1   0.003743 **  # add group effect on slope
# ..3    19 68769 68938 -34366    68731   1.1115      1   0.291749     # add group effect on quadratic poly 
# ..4    20 68771 68948 -34365    68731   0.5965      1   0.439899     # add cond effect on intercept
# ..5    21 68765 68951 -34362    68723   7.4548      1   0.006327 **  # add cond effect on slope 
# ..6    22 68760 68955 -34358    68716   6.7032      1   0.009624 **  # add cond effect on quadratic poly
# ..7    25 68656 68878 -34303    68606 110.4731      3  < 2.2e-16 *** # full model 



summary(gc_mod2_full)

# Fixed effects:
#                               Estimate Std. Error         df t value Pr(>|t|)
# (Intercept)                  5.384e-01  1.762e-02  7.000e+01  30.550  < 2e-16  ***
# ot1                          6.557e-01  8.989e-02  6.300e+01   7.294 6.11e-10  ***
# ot2                          3.631e-01  5.487e-02  6.400e+01   6.617 8.90e-09  ***
# groupSum1                    1.649e-04  1.191e-02  4.900e+01   0.014  0.98901     
# conditionSum1               -8.414e-03  1.314e-02  3.200e+01  -0.640  0.52671     
# ot1:groupSum1                1.838e-01  5.766e-02  4.700e+01   3.187  0.00255  ** 
# ot2:groupSum1                4.226e-02  3.982e-02  4.700e+01   1.061  0.29386     
# ot1:conditionSum1           -2.320e-01  7.072e-02  3.200e+01  -3.280  0.00253  ** 
# ot2:conditionSum1            1.150e-01  4.092e-02  3.200e+01   2.809  0.00840  ** 
# groupSum1:conditionSum1     -1.817e-02  2.077e-03  5.180e+04  -8.750  < 2e-16  ***
# ot1:groupSum1:conditionSum1  9.025e-02  1.582e-02  5.180e+04   5.706 1.16e-08  ***
# ot2:groupSum1:conditionSum1  1.923e-02  1.581e-02  5.182e+04   1.216  0.22409     





# create new df including the fitted model 
data_comp2 <- data.frame(na.omit(lsrl_gc_ss_la), 
                        GCA_Full = fitted(gc_mod2_full))
# glimpse(data_comp2)


suffix_area2 <- data.frame(x = 53:58, y = Inf)

condition_namesMod2LSRL <- c(
                    `ss` = "SS", 
                    `la` = "LA"
                    )

data_comp2 %>% 
  ggplot(., aes(x = binGC, y = targetProp, color = condition)) + 
  facet_grid(. ~ group, labeller = as_labeller(condition_namesMod2LSRL)) + 
  geom_area(data = suffix_area2, aes(x = x, y = y), inherit = FALSE, alpha = 0.3, fill = 'lightcyan2') +
  stat_summary(fun.data = mean_se, geom = 'errorbar', 
               show.legend = FALSE, size = 0.1) +
  stat_summary(fun.y = mean, geom = 'point', size = 0.2) + 
  stat_summary(aes(y = GCA_Full, color = condition), fun.y = mean, geom = 'line', size = 0.4) + 
  xlab("Adjusted time course") +
  ylab("Target fixations") +
  coord_cartesian(ylim = c(0.0, 1.0)) + 
  scale_x_continuous(breaks = c(1, 53), labels = c("Approx.\ntarget\nonset", "Target\nsyllable\nonset")) + 
  scale_color_brewer(palette = "Set1", name = "", labels = c("Paroxytone", "Oxytone")) + 
  theme_bw(base_size = 16, base_family = "Times New Roman") -> stressGCAfullMod2

# ggsave('stressGCAfullMod2.png', plot = stressGCAfullMod2, dpi = 600, device = "png", path = "./mySources/figs/stress/s1_beg_adv_nat/eye_track")
