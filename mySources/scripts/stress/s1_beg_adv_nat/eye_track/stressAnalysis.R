#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Morphosyntactic predictability                                              #
# Stress analyses                                                             #
# 11/09/2016                                                                  #
# Script 1                                                                    #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# clean working directory
rm(list = ls(all = TRUE))

## @knitr stressLibs

library(plotly); library(tidyverse); library(broom); library(sjPlot)
library(lme4); library(lmerTest); library(gridExtra); library(cowplot)

## @knitr ignore

# Set working directory
# setwd("~/Desktop/morph_pred/")
setwd("~/academia/research/in_progress/morph_pred")


# read data
# df <- read.csv("./mySources/data/stressBIN5Clean.csv", header = TRUE, quote = "", sep = ',')
# df_stress_temp <- read_csv("./mySources/data/stressBIN5Clean.csv")

df_stress_temp <- read_csv("./mySources/data/stressBIN10Clean.csv")


## @knitr stressSubsets

# Remove particpants (working memory)
remove_partLB <- c('l01', 'l02', 'l03', 'l04', 'l05', 'l06', 'l07', 'l08', 'l09', 'l10', 'l15', 'l20', 'l21', 'l22', 'l23', 'l26', 'l30', 'l31', 'l33')
remove_partLA <- c('la04', 'la06', 'la07', 'la14')

# summary(as.factor(df_stress_temp$participant))


# glimpse(df_stress_temp)

df_stress_temp1 <- df_stress_temp[!(df_stress_temp$participant) %in% remove_partLB, ]
df_stress <- df_stress_temp1[!(df_stress_temp1$participant) %in% remove_partLA, ]

# summary(as.factor(df_stress$participant))

# - for target words without codas the target syllable onset 
#   is 'word3_c2'
# - for words with codas the target syllable onset is 'word3_c3'
# - we add a new column called 'targetSylOnset'
# - initialy it is just a copy of 'word3_c3'
# - we make a vector of words without coda and add the value from 'word3_c2' 
#   to 'targetSylOnset'

df_stress$targetSylOnset <- df_stress$word3_c3
noCodas <- c('bebe', 'bebio', 'llena', 'lleno', 'sube', 'subio', 'come', 'comio', 'saca', 'saco', 'lava', 'lavo', 'graba', 'grabo')
df_stress[df_stress$target %in% noCodas, 'targetSylOnset'] <- df_stress[df_stress$target %in% noCodas, 'word3_c2']

# glimpse(df_stress)




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

# Where does target suffic begin in time course
suffixOnsets <- df_stress %>%
  group_by(., target) %>%
  summarize(., sufOnset = (unique(targetSylOnset) - unique(word2_c1v1)) / 10)

# Center time course so that suffix onset = 0
suffixOnsetAdj <- df_stress %>%
  group_by(., target) %>%
  summarize(., sufOnsetAdj = ((unique(targetSylOnset) - unique(word2_c1v1)) / 10) - 
                        (unique(targetSylOnset) - unique(word2_c1v1)) / 10)

# Where does the target word begin in the time course?
twOnsets <- df_stress %>%
  group_by(., target) %>%
  summarize(., twOnset = (unique(word3_c1v1) - unique(word2_c1v1)) / 10)


# Adjust to centered time course 
twOnsetAdj <- cbind(suffixOnsets, twOnsets[, 2])
twOnsetAdj <- mutate(twOnsetAdj, diff = sufOnset - twOnset, twOnsetAdj = 0 - diff)

df_stress %>%
  ggplot(., aes(x = bin, y = target)) + 
  geom_point(alpha = 0.2, size = 0.2) + 
  geom_point(data = suffixOnsets, aes(x = sufOnset, y = target), color = 'red') + 
  geom_point(data = twOnsets, aes(x = twOnset, y = target), color = 'blue') + 
  xlab("Time course") + ylab("Items") + 
  theme_bw(base_size = 12, base_family = "Times") -> onsetP

df_stress %>%
  ggplot(., aes(x = binAdj, y = target)) + 
  geom_point(alpha = 0.2, size = 0.2) +
  geom_point(data = suffixOnsetAdj, aes(x = sufOnsetAdj, y = target), color = 'red') + 
  geom_point(data = twOnsetAdj, aes(x = twOnsetAdj, y = target), color = 'blue') + 
  geom_vline(xintercept = binAdjMinMax, color = 'red') + 
  geom_vline(xintercept = binAdjMaxMin, color = 'red') + 
  scale_y_discrete(position = "right") +
  xlab("Adjusted time course") + ylab("Items") + 
  theme_bw(base_size = 12, base_family = "Times") -> onsetAdjP

stressBinPlots <- plot_grid(onsetP, onsetAdjP, ncol = 2)

# ggsave('stressBinPlots.png', plot = stressBinPlots, dpi = 600, device = "png", path = "./mySources/figs/stress/general")








## @knitr stressPlotAlldata

# Time sequence for plots 

# Create vector of times for canta 
cantaTimes <- c(100, 341, 1156, 1244, 1420, 1515, 1565, 1699, 1878, 2684)
# Create vector of labels 
cantaSeq   <- c("la", "senora", "c", "a", "n", "t", "a", "la", "cancion", "end")
# add labels as names arg for vector cantaTimes
names(cantaTimes) <- cantaSeq
# adjust bins 
cantaAdj <- (cantaTimes / 10) - (cantaTimes[6] / 10)
# turn in into a dataframe 
canta_df_temp <- data.frame(group = rep(c("lb", "la", "ss"), each = 10), 
                            step  = 1:10, 
                            condition = "stressed", 
                            binAdj     = cantaAdj, 
                            targetProp = 0.05, 
                            text  = names(cantaAdj))





cantoTimes <- c(100, 377, 1137, 1234, 1343, 1415, 1538, 1735, 1929, 2807)
# Create vector of labels 
cantoSeq   <- c("la", "senora", "c", "a", "n", "t", "o", "la", "cancion", "end")
# add labels as names arg for vector cantaTimes
names(cantoTimes) <- cantoSeq
# adjust bins 
cantoAdj <- (cantoTimes / 10) - (cantoTimes[6] / 10)
# turn in into a dataframe 
canto_df_temp <- data.frame(group = rep(c("lb", "la", "ss"), each = 10), 
                            step  = 1:10, 
                            condition = "unstressed", 
                            binAdj     = cantoAdj, 
                            targetProp     = 0.05, 
                            text  = names(cantoAdj))
cantEx <- rbind(canta_df_temp, canto_df_temp)






# Relevel factor
df_short$group <- factor(df_short$group, levels = c("lb", "la", "hs", "ss", "int"))

condition_names <- c(
                    `stressed` = "Paroxytone",
                    `unstressed` = "Oxytone", 
                    `lb` = "LB", 
                    `la` = "LA", 
                    `hs` = "HS", 
                    `ss` = "SS", 
                    `int` = "IN"
                    )

df_short %>% 
  na.omit(.) %>% 
  ggplot(., aes(x = binAdj, y = targetProp)) + 
  facet_grid(group ~ condition, labeller = as_labeller(condition_names)) + 
  geom_vline(xintercept = 0, color = 'grey60') + 
  geom_hline(yintercept = 0.5, color = 'grey60') + 
  stat_summary(fun.data = mean_se, geom = 'errorbar', width = 0.2, size = 0.1,
               show.legend = FALSE) +
  stat_summary(fun.y = mean, geom = 'point', size = 0.2, color = 'darkgreen') + 
  stat_summary(data = na.omit(df_short), aes(x = binAdj, y = distractorProp), fun.data = mean_se, geom = 'errorbar', width = 0.2, size = 0.1) + 
  stat_summary(data = na.omit(df_short), aes(x = binAdj, y = distractorProp), fun.y = mean, geom = 'point', size = 0.2, color = 'red') + 
  ylab('Proportion of fixations') + 
  xlab('Adjusted time course') + 
  theme_bw(base_size = 16, base_family = "Times") -> stressPlotAll






## @knitr stressPlotLSRL

# PLOT FOR LSRL

condition_namesLSRL <- c(
                    `stressed` = "Paroxytone",
                    `unstressed` = "Oxytone", 
                    `lb` = "LB", 
                    `la` = "LA", 
                    `ss` = "SS" 
                    )

lsrl_subset <- df_short %>% filter(., group %in% c('lb', 'la', 'ss'))

lsrl_subset %>%
  na.omit(.) %>% 
  ggplot(., aes(x = binAdj, y = targetProp)) + 
  facet_grid(group ~ condition, labeller = as_labeller(condition_namesLSRL)) + 
  geom_vline(xintercept = 0, color = 'grey60') + 
  geom_hline(yintercept = 0.5, color = 'grey60') + 
  stat_summary(fun.data = mean_cl_boot, geom = 'errorbar', width = 0, size = 0.1,
               show.legend = FALSE, color = 'darkgrey') +
  stat_summary(fun.y = mean, geom = 'point', size = 0.2, color = 'darkgreen') + 
  stat_summary(data = na.omit(lsrl_subset), aes(x = binAdj, y = distractorProp), fun.data = mean_cl_boot, geom = 'errorbar', width = 0, size = 0.1, color = 'darkgrey') + 
  stat_summary(data = na.omit(lsrl_subset), aes(x = binAdj, y = distractorProp), fun.y = mean, geom = 'point', size = 0.2, color = 'red') + 
  ylab('Proportion of fixations') + 
  xlab('Adjusted time course') + 
  ylim(0, 1) + 
  geom_text(data = cantEx[cantEx$step >= 3 & cantEx$step < 8, ], aes(label = text), color = 'blue') + 
  theme_bw(base_size = 16, base_family = "Times") -> stressLSRLp1

# ggsave('stressLSRLp1.png', plot = stressLSRLp1, dpi = 600, device = "png", path = "./mySources/figs/stress/s1_beg_adv_nat/eye_track")























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
# target onset bin (binAdj = 0)
lsrl_subset_0 <- lsrl_subset %>% filter(., binAdj == 0)

# Quick and dirty mean of target fixations as a function of 
# group and condition (stressed, unstressed *1st syllable*)
lsrl_subset_0 %>% 
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

stress_ttest <- lsrl_subset_0 %>%
  na.omit(.) %>% 
  group_by(., group, condition, participant) %>%
  summarise(., meanFix = mean(targetProp)) %>% 
  do(tidy(t.test(.$meanFix, alternative = "greater", mu = 0.5, conf.level = 0.95)))

# Convert pvalues from scientific notation
stress_ttest$p.value <- format(stress_ttest$p.value, scientific = F)
stress_ttest$sig <- "N.S."
stress_ttest[stress_ttest$p.value <= 0.05, 'sig'] <- "*"

# Print results
# print(as.data.frame(stress_ttest[, c(1:7, 11)]))

#   group  condition  estimate statistic          p.value parameter  conf.low   sig
# 1    lb   stressed 0.5358135 0.6104244 0.27699151201698        11 0.4304492  N.S.
# 2    lb unstressed 0.5711806 1.5804339 0.07115660755082        11 0.4902964  N.S.
# 3    la   stressed 0.6555021 3.5018326 0.00087869973568        25 0.5796506     *
# 4    la unstressed 0.6974359 5.1505881 0.00001262759817        25 0.6319583     *
# 5    ss   stressed 0.7109472 4.8578837 0.00003718370061        22 0.6363825     *
# 6    ss unstressed 0.7895186 7.9056402 0.00000003592184        22 0.7266337     *




# Only LB do not predict above chance 

# We will plot the models 
# This will be almost exactly the same as 
# the previous plot, but it will use the 
# confidence interval from the test we 
# actually conducted
stress_ttest$group <- factor(stress_ttest$group, levels = c("lb", "la", "ss"))
stress_ttest$condition <- factor(stress_ttest$condition, levels = c("stressed", "unstressed"))


ggplot(stress_ttest, aes(x = group, y = estimate, color = condition, 
                      group = interaction(group, condition), dodge = condition)) +
    geom_hline(yintercept = 0.5, lty = 3) + 
    geom_linerange(aes(ymin = conf.low, ymax = estimate), color = 'grey40',
                     position = position_dodge(width = 0.75), size = 1) +
    geom_point(position = position_dodge(width = 0.75), size = 4) +
    ylim(0, 1.0) + ylab('Target fixations') + xlab('Group') + 
    ggtitle('Mean target fixations and lower-bound 95% confidence interval') + 
    scale_color_brewer(palette = "Set1", name = '', labels = c('Paroxytone', 'Oxytone')) + 
    scale_x_discrete(labels = c('LB', 'LA', 'SS')) + 
    theme_bw(base_size = 16, base_family = 'Times') -> stressTargetFixMODP3

# Looks good, save as .png file. 
# ggsave('stressTargetFixMODP3.png', plot = stressTargetFixMODP3, dpi = 600, device = "png", path = "./mySources/figs/stress/s1_beg_adv_nat/eye_track")




## @knitr stressGroupCompare


# Are groups different from each other?
lsrl_subset_0_prop <- lsrl_subset_0 %>%
  na.omit(.) %>% 
  group_by(., group, condition, participant) %>%
  summarise(., meanFix = mean(targetProp)) 

# Add working memory covariate
wm_df <- structure(list(participant = structure(c(1L, 2L, 3L, 4L, 5L, 
6L, 7L, 8L, 9L, 10L, 11L, 12L, 1L, 2L, 3L, 4L, 5L, 6L, 7L, 8L, 
9L, 10L, 11L, 12L, 13L, 14L, 15L, 16L, 17L, 18L, 19L, 20L, 21L, 
22L, 23L, 24L, 25L, 26L, 27L, 28L, 29L, 30L, 31L, 32L, 33L, 34L, 
35L, 36L, 37L, 38L, 13L, 14L, 15L, 16L, 17L, 18L, 19L, 20L, 21L, 
22L, 23L, 24L, 25L, 26L, 27L, 28L, 29L, 30L, 31L, 32L, 33L, 34L, 
35L, 36L, 37L, 38L, 39L, 48L, 51L, 54L, 55L, 56L, 57L, 58L, 59L, 
40L, 41L, 42L, 43L, 44L, 45L, 46L, 47L, 49L, 50L, 60L, 61L, 52L, 
53L, 39L, 48L, 51L, 54L, 55L, 56L, 57L, 58L, 59L, 40L, 41L, 42L, 
43L, 44L, 45L, 46L, 47L, 49L, 50L, 60L, 61L, 52L, 53L), .Label = c("l11", 
"l12", "l13", "l16", "l17", "l18", "l19", "l24", "l25", "l27", 
"l28", "l29", "la01", "la02", "la03", "la05", "la08", "la09", 
"la10", "la11", "la12", "la13", "la15", "la16", "la17", "la18", 
"la19", "la20", "la21", "la22", "la23", "la24", "la25", "la26", 
"la27", "la28", "la29", "la30", "s1", "s10", "s11", "s12", "s13", 
"s14", "s16", "s17", "s19", "s2", "s20", "s21", "s3", "s32", 
"s38", "s4", "s5", "s6", "s7", "s8", "s9", "ss10", "ss30"), class = "factor"), 
    wm = c(6L, 8L, 4L, 1L, 9L, 4L, 12L, 5L, 7L, 6L, 6L, 6L, 6L, 
    8L, 4L, 1L, 9L, 4L, 12L, 5L, 7L, 6L, 6L, 6L, 10L, 13L, 1L, 
    11L, 9L, 8L, 10L, 8L, 10L, 9L, 11L, 7L, 8L, 10L, 8L, 8L, 
    10L, 9L, 9L, 10L, 8L, 10L, 11L, 9L, 8L, 10L, 10L, 13L, 1L, 
    11L, 9L, 8L, 10L, 8L, 10L, 9L, 11L, 7L, 8L, 10L, 8L, 8L, 
    10L, 9L, 9L, 10L, 8L, 10L, 11L, 9L, 8L, 10L, 6L, 10L, 11L, 
    6L, 11L, 9L, 10L, 7L, 8L, 12L, 8L, 10L, 11L, 10L, 6L, 12L, 
    8L, 9L, 14L, 8L, 10L, 12L, 6L, 6L, 10L, 11L, 6L, 11L, 9L, 
    10L, 7L, 8L, 12L, 8L, 10L, 11L, 10L, 6L, 12L, 8L, 9L, 14L, 
    8L, 10L, 12L, 6L)), .Names = c("participant", "wm"), class = "data.frame", row.names = c(NA, 
-122L))

lsrl_subset_0_prop$wm <- wm_df$wm



lsrl_subset_0_prop$group <- factor(lsrl_subset_0_prop$group, levels = c("ss", "la", "lb"))

prop_0_mod_0 <- lmer(meanFix ~ 1 + 
                    (1 | participant), 
                    data = lsrl_subset_0_prop, REML = F, control = lmerControl(optimizer = 'bobyqa'))

prop_0_mod_group <- lmer(meanFix ~ 1 + group + 
                        (1 | participant), 
                        data = lsrl_subset_0_prop, REML = F, control = lmerControl(optimizer = 'bobyqa'))

prop_0_mod_cond <- lmer(meanFix ~ 1 + group + condition + 
                        (1 | participant), 
                        data = lsrl_subset_0_prop, REML = F, control = lmerControl(optimizer = 'bobyqa'))

prop_0_mod_full <- lmer(meanFix ~ 1 + group * condition + 
                      (1 | participant), 
                       data = lsrl_subset_0_prop, REML = F, control = lmerControl(optimizer = 'bobyqa'))

prop_0_mod_cov <- lmer(meanFix ~ 1 + group * condition + wm + 
                      (1 | participant), 
                       data = lsrl_subset_0_prop, REML = F, control = lmerControl(optimizer = 'bobyqa'))

prop_0_mod_cov_full <- lmer(meanFix ~ 1 + group * condition * wm + 
                      (1 | participant), 
                       data = lsrl_subset_0_prop, REML = F, control = lmerControl(optimizer = 'bobyqa'))


# anova(prop_0_mod_0, prop_0_mod_group, prop_0_mod_cond, prop_0_mod_full, prop_0_mod_cov, prop_0_mod_cov_full, test = "Chisq")

#        Df     AIC      BIC logLik deviance   Chisq Chi Df Pr(>Chisq)    
# object  3 -30.302 -21.8895 18.151  -36.302                              
# ..1     5 -41.247 -27.2271 25.624  -51.247 14.9456      2  0.0005683 ***
# ..2     6 -41.621 -24.7972 26.811  -53.621  2.3742      1  0.1233552    
# ..3     8 -37.911 -15.4788 26.956  -53.911  0.2896      2  0.8651943    
# ..4     9 -36.484 -11.2482 27.242  -54.484  0.5734      1  0.4489036    
# ..5    14 -34.897   4.3597 31.448  -62.897  8.4122      5  0.1349345    

# summary(prop_0_mod_group)

# Fixed effects:
#              Estimate Std. Error        df t value Pr(>|t|)    
# (Intercept)   0.75023    0.02892 122.00000  25.943  < 2e-16 ***
# groupla      -0.07376    0.03970 122.00000  -1.858 0.065569 .  
# grouplb      -0.19674    0.04939 122.00000  -3.984 0.000116 ***

# Relevel to test lb vs la 
lsrl_subset_0_prop$group <- factor(lsrl_subset_0_prop$group, levels = c("lb", "la",  "ss"))

summary(lmer(meanFix ~ 1 + group + (1 | participant), data = lsrl_subset_0_prop, REML = F, control = lmerControl(optimizer = 'bobyqa')))





# Calculate mean target fixation as a function of group, condition, 
# for each participant. We will plot the mean and calculate the 
# bootstrapped 95% confidence interval and plot it all. 
lsrl_subset_0 %>%
  na.omit(.) %>%
  group_by(., group, condition, participant) %>%
  summarise(., meanFix = mean(targetProp)) %>%
  ggplot(., aes(x = group, y = meanFix, 
                dodge = condition, color = condition,
                group = interaction(group, condition))) +
    geom_hline(yintercept = 0.5, color = "black", size = 0.75, 
               lty = 3) + 
    stat_summary(fun.data = mean_cl_boot, geom = 'errorbar', 
                 position = position_dodge(width = 0.5), 
                 width = 0.35, color = 'grey40') + 
    stat_summary(fun.y = mean, geom = 'point', size = 4,
                 position = position_dodge(width = 0.5)) + 
    coord_cartesian(ylim = c(0, 1)) + ylab('Target fixations') + xlab('Group') + 
    scale_x_discrete(labels = c('LB', 'LA', 'SS')) + 
    ggtitle('Mean target fixations as a function of group and target type') + 
    scale_color_brewer(palette = "Set1", name = '', labels = c('Paroxytone', 'Oxytone')) + 
    theme_bw(base_size = 16, base_family = 'Times') -> stressTargetFixLSRLp2

# Looks good, save as .png file
# ggsave('stressTargetFixPlot.png', plot = stressTargetFixLSRLp2, dpi = 600, device = "png", path = "./mySources/figs/stress/s1_beg_adv_nat/eye_track")


























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
lsrl_gc_subset <- filter(lsrl_subset, binAdj >= -52 & binAdj <= 5)


# - Readjust time course 
#    - now we will make the time course positive, starting at 1
#    - to do this we will add the lowest 'binAdj' value to each 
#      bin, plus 1 (to avoid starting at 0)

# lsrl_gc_subset$binGC <- lsrl_gc_subset$binAdj + 51
lsrl_gc_subset$binGC <- lsrl_gc_subset$binAdj + 53


# - Now we add higher order polynomials for analyses

t <- poly(min(lsrl_gc_subset$binGC):max(lsrl_gc_subset$binGC), 3)
lsrl_gc_subset[, paste('ot', 1:3, sep = "")] <- t[lsrl_gc_subset$binGC, 1:3]

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
lsrl_gc_subset$group <- factor(lsrl_gc_subset$group, levels = c("ss", "la", "lb"))

# Set condition as factor and set conding to contrast 
lsrl_gc_subset$condition <- as.factor(lsrl_gc_subset$condition)

# contrasts(lsrl_gc_subset$condition)
# contrasts(lsrl_gc_subset$group)

lsrl_gc_subset$conditionSum <- C(lsrl_gc_subset$condition, sum)
# contrasts(lsrl_gc_subset$conditionSum)

lsrl_gc_subset$groupSum <- C(lsrl_gc_subset$group, sum)
# contrasts(lsrl_gc_subset$groupSum)

## @knitr ignore2

# load the models:
gc_mod_base    <- readRDS('./mySources/models/stress/s1_beg_adv_nat/eye_track/gc_mod_base.rds')
gc_mod_group_0 <- readRDS('./mySources/models/stress/s1_beg_adv_nat/eye_track/gc_mod_group_0.rds')
gc_mod_group_1 <- readRDS('./mySources/models/stress/s1_beg_adv_nat/eye_track/gc_mod_group_1.rds')
gc_mod_group_2 <- readRDS('./mySources/models/stress/s1_beg_adv_nat/eye_track/gc_mod_group_2.rds')
gc_mod_cond_0  <- readRDS('./mySources/models/stress/s1_beg_adv_nat/eye_track/gc_mod_cond_0.rds')
gc_mod_cond_1  <- readRDS('./mySources/models/stress/s1_beg_adv_nat/eye_track/gc_mod_cond_1.rds')
gc_mod_cond_2  <- readRDS('./mySources/models/stress/s1_beg_adv_nat/eye_track/gc_mod_cond_2.rds')
gc_mod_full    <- readRDS('./mySources/models/stress/s1_beg_adv_nat/eye_track/gc_mod_full.rds')


# Base model 
if(F){
  gc_mod_base <- lmer(targetProp ~ (ot1+ot2) + 
                 ((ot1+ot2) | participant) + 
                 ((ot1+ot2) | target),
                 control = lmerControl(optimizer = 'bobyqa'), 
                 data = lsrl_gc_subset, REML = F)
  saveRDS(gc_mod_base, file = "./mySources/models/stress/s1_beg_adv_nat/eye_track/gc_mod_base.rds", compress = 'xz')
}


# Add group effect on intercept 
if(F){
  gc_mod_group_0 <- lmer(targetProp ~ (ot1+ot2) + group + 
                    ((ot1+ot2) | participant) + 
                    ((ot1+ot2) | target), 
                    control = lmerControl(optimizer = 'bobyqa'), 
                    data = lsrl_gc_subset, REML = F)
  saveRDS(gc_mod_group_0, file = "./mySources/models/stress/s1_beg_adv_nat/eye_track/gc_mod_group_0.rds", compress = 'xz')
}

# Add group effect on slope
if(F){
 gc_mod_group_1 <- lmer(targetProp ~ (ot1+ot2) + group + 
                   ot1:group + 
                   ((ot1+ot2) | participant) + 
                   ((ot1+ot2) | target),
                   control = lmerControl(optimizer = 'bobyqa'), 
                   data = lsrl_gc_subset, REML = F)
  saveRDS(gc_mod_group_1, file = "./mySources/models/stress/s1_beg_adv_nat/eye_track/gc_mod_group_1.rds", compress = 'xz')
}

# Add group effect on quadratic poly 
if(F){
  gc_mod_group_2 <- lmer(targetProp ~ (ot1+ot2) + group + 
                    ot1:group + ot2:group + 
                    ((ot1+ot2) | participant) + 
                    ((ot1+ot2) | target),
                    control = lmerControl(optimizer = 'bobyqa'), 
                    data = lsrl_gc_subset, REML = F)
  saveRDS(gc_mod_group_2, file = "./mySources/models/stress/s1_beg_adv_nat/eye_track/gc_mod_group_2.rds", compress = 'xz')
}


# Add condition effect on intercept 
if(F){
  gc_mod_cond_0 <- lmer(targetProp ~ (ot1+ot2) * group + condition + 
                   ((ot1+ot2) | participant) + 
                   ((ot1+ot2) | target),
                   control = lmerControl(optimizer = 'bobyqa'), 
                   data = lsrl_gc_subset, REML = F)
  saveRDS(gc_mod_cond_0, file = "./mySources/models/stress/s1_beg_adv_nat/eye_track/gc_mod_cond_0.rds", compress = 'xz')
}

# Add condition effect on slope 
if(F){
  gc_mod_cond_1 <- lmer(targetProp ~ (ot1+ot2) * group + condition + 
                   ot1:condition + 
                   ((ot1+ot2) | participant) + 
                   ((ot1+ot2) | target),
                   control = lmerControl(optimizer = 'bobyqa'), 
                   data = lsrl_gc_subset, REML = F)
  saveRDS(gc_mod_cond_1, file = "./mySources/models/stress/s1_beg_adv_nat/eye_track/gc_mod_cond_1.rds", compress = 'xz')
}

# Add condition effect on quadratic poly 
if(F){
  gc_mod_cond_2 <- lmer(targetProp ~ (ot1+ot2) * group + condition + 
                   ot1:condition + ot2:condition + 
                   ((ot1+ot2) | participant) + 
                   ((ot1+ot2) | target),
                   control = lmerControl(optimizer = 'bobyqa'), 
                   data = lsrl_gc_subset, REML = F)
  saveRDS(gc_mod_cond_2, file = "./mySources/models/stress/s1_beg_adv_nat/eye_track/gc_mod_cond_2.rds", compress = 'xz')
}


# Include all interactions
if(F){
gc_mod_full <- lmer(targetProp ~ (ot1+ot2) * group * condition + 
               ((ot1+ot2) | participant) + 
               ((ot1+ot2) | target),
               control = lmerControl(optimizer = 'bobyqa'), 
               data = lsrl_gc_subset, REML = F)
  saveRDS(gc_mod_full, file = "./mySources/models/stress/s1_beg_adv_nat/eye_track/gc_mod_full.rds", compress = 'xz')
}

# Model comparison 
# - subsequent models test three time terms: 
#    - the intercept (additive effects)
#    - the linear slope (ot1)
#    - the steepness of the quadratic curvature (ot2)

# anova(gc_mod_base, 
#       gc_mod_group_0, 
#       gc_mod_group_1, 
#       gc_mod_group_2, 
#       gc_mod_cond_0, 
#       gc_mod_cond_1, 
#       gc_mod_cond_2, 
#       gc_mod_full, test = 'Chisq')

#        Df   AIC   BIC logLik deviance    Chisq Chi Df Pr(>Chisq)    
# object 16 97453 97601 -48711    97421                                # base model 
# ..1    18 97453 97618 -48708    97417   4.6033      2   0.100095     # add group effect on intercept 
# ..2    20 97445 97628 -48702    97405  12.3091      2   0.002124 **  # add group effect on slope
# ..3    22 97445 97647 -48700    97401   3.7757      2   0.151399     # add group effect on quadratic poly 
# ..4    23 97445 97656 -48699    97399   2.0400      1   0.153205     # add cond effect on intercept
# ..5    24 97431 97652 -48692    97383  15.2828      1  9.255e-05 *** # add cond effect on slope 
# ..6    25 97430 97660 -48690    97380   3.3546      1   0.067017 .   # add cond effect on quadratic poly
# ..7    31 97120 97405 -48529    97058 321.7254      6  < 2.2e-16 *** # full model 



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
data.comp <- data.frame(na.omit(lsrl_gc_subset), 
                        GCA_Full = fitted(gc_mod_full))
# glimpse(data.comp)


suffix_area <- data.frame(x = 53:58, y = Inf)

condition_namesModLSRL <- c(
                    `ss` = "SS", 
                    `la` = "LA", 
                    `lb` = "LB" 
                    )

data.comp %>% 
  ggplot(., aes(x = binGC, y = targetProp, color = condition)) + 
  facet_grid(. ~ group, labeller = as_labeller(condition_namesModLSRL)) + 
  geom_area(data = suffix_area, aes(x = x, y = y), inherit = FALSE, alpha = 0.3, fill = 'lightcyan2') +
  stat_summary(fun.data = mean_se, geom = 'errorbar', 
               show.legend = FALSE, size = 0.1) +
  stat_summary(fun.y = mean, geom = 'point', size = 0.2) + 
  stat_summary(aes(y = GCA_Full, color = condition), fun.y = mean, geom = 'line', size = 0.4) + 
  xlab("Adjusted time course") +
  ylab("Target fixations") +
  coord_cartesian(ylim = c(0.0, 1.0)) + 
  scale_x_continuous(breaks = c(1, 53), labels = c("Approx.\ntarget\nonset", "Target\nsyllable\nonset")) + 
  scale_color_brewer(palette = "Set1", name = "", labels = c("Paroxytone", "Oxytone")) + 
  theme_bw(base_size = 16, base_family = "Times New Roman") -> stressGCAfullMod

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
