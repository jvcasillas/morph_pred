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
remove <- c('hs01', 'hs06', 'hs10', 'la05', 'la24')

# summary(as.factor(df_stress_temp$participant))


# glimpse(df_stress_temp)

df_stress <- df_stress_temp[!(df_stress_temp$participant) %in% remove, ]

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
canta_df_temp <- data.frame(group = rep(c("la", "hs", "ss"), each = 10), 
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
canto_df_temp <- data.frame(group = rep(c("la", "hs", "ss"), each = 10), 
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






## @knitr stressPlotHLS

# PLOT FOR HLS

condition_namesHLS <- c(
                    `stressed` = "Paroxytone",
                    `unstressed` = "Oxytone", 
                    `la` = "LL", 
                    `hs` = "EL", 
                    `ss` = "SS" 
                    )

hls_subset <- df_short %>% filter(., group %in% c('la', 'hs', 'ss'))

hls_subset %>%
  na.omit(.) %>% 
  ggplot(., aes(x = binAdj, y = targetProp)) + 
  facet_grid(group ~ condition, labeller = as_labeller(condition_namesHLS)) + 
  geom_vline(xintercept = 0, color = 'grey60') + 
  geom_hline(yintercept = 0.5, color = 'grey60') + 
  stat_summary(fun.data = 'mean_cl_boot', geom = 'errorbar', width = 0, size = 0.1,
               show.legend = FALSE, color = 'black') +
  stat_summary(fun.y = 'mean', geom = 'point', size = 0.2, color = 'darkgreen') + 
  stat_summary(data = na.omit(hls_subset), aes(x = binAdj, y = distractorProp), fun.data = 'mean_cl_boot', geom = 'errorbar', width = 0, size = 0.1, color = 'black') + 
  stat_summary(data = na.omit(hls_subset), aes(x = binAdj, y = distractorProp), fun.y = 'mean', geom = 'point', size = 0.2, color = 'red') + 
  ylab('Proportion of fixations') + 
  xlab('Adjusted time course') + 
  ylim(0, 1) + 
  geom_text(data = cantEx[cantEx$step >= 3 & cantEx$step < 8 & cantEx$group != "lb", ], aes(label = text), color = 'blue') + 
  theme_bw(base_size = 16, base_family = "Times") -> stressHLSp1

# ggsave('stressHLSp1.png', plot = stressHLSp1, dpi = 600, device = "png", path = "./mySources/figs/stress/s2_adv_her_nat/eye_track")























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
hls_subset_0 <- hls_subset %>% filter(., binAdj == 0)

# Quick and dirty mean of target fixations as a function of 
# group and condition (stressed, unstressed *1st syllable*)
hls_subset_0 %>% 
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

stress_ttest <- hls_subset_0 %>%
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

# group  condition  estimate statistic          p.value parameter  conf.low sig
#    la   stressed 0.6198413  2.818709 0.00445832937620        27 0.5474236   *
#    la unstressed 0.7152530  6.201943 0.00000062240234        27 0.6561364   *
#    hs   stressed 0.6483696  4.384343 0.00011797253742        22 0.5902601   *
#    hs unstressed 0.7451087  7.728286 0.00000005212898        22 0.6906481   *
#    ss   stressed 0.7109472  4.857884 0.00003718370061        22 0.6363825   *
#    ss unstressed 0.7895186  7.905640 0.00000003592184        22 0.7266337   *




# All predict above chance 

# We will plot the models 
# This will be almost exactly the same as 
# the previous plot, but it will use the 
# confidence interval from the test we 
# actually conducted
stress_ttest$group <- factor(stress_ttest$group, levels = c("la", "hs", "ss"))
stress_ttest$condition <- factor(stress_ttest$condition, levels = c("stressed", "unstressed"))


ggplot(stress_ttest, aes(x = group, y = estimate, color = condition, 
                      group = interaction(group, condition), dodge = condition)) +
    geom_hline(yintercept = 0.5, lty = 3) + 
    geom_linerange(aes(ymin = conf.low, ymax = estimate), color = 'grey40',
                     position = position_dodge(width = 0.75), size = 1) +
    geom_point(position = position_dodge(width = 0.75), size = 4) +
    # ylim(0, 1.0) + 
    ylab('Target fixations') + xlab('Group') + 
    ggtitle('Mean target fixations and lower-bound 95% confidence interval') + 
    scale_color_brewer(palette = "Set1", name = '', labels = c('Paroxytone', 'Oxytone')) + 
    scale_x_discrete(labels = c('LA', 'HS', 'SS')) + 
    theme_bw(base_size = 16, base_family = 'Times') -> stressTargetFixMODP3

# Looks good, save as .png file. 
# ggsave('stressTargetFixMODP3.png', plot = stressTargetFixMODP3, dpi = 600, device = "png", path = "./mySources/figs/stress/s2_adv_her_nat/eye_track")




## @knitr stressGroupCompare


# Are groups different from each other?
hls_subset_0_glmm <- hls_subset_0 %>%
  na.omit(.) %>% 
  group_by(., group, condition, participant) %>%
  summarise(., meanFix = mean(targetProp)) 

# Add working memory covariate
wm_df <- structure(list(participant = structure(1:74, .Label = c("hs02", 
"hs03", "hs04", "hs05", "hs07", "hs08", "hs09", "hs11", "hs13", 
"hs14", "hs15", "hs17", "hs18", "hs20", "hs21", "hs22", "hs24", 
"hs25", "hs27", "hs28", "hs29", "hs30", "hs31", "la01", "la02", 
"la03", "la04", "la06", "la07", "la08", "la09", "la10", "la11", 
"la12", "la13", "la14", "la15", "la16", "la17", "la18", "la19", 
"la20", "la21", "la22", "la23", "la25", "la26", "la27", "la28", 
"la29", "la30", "s01", "s02", "s03", "s04", "s05", "s06", "s07", 
"s08", "s09", "s10", "s11", "s12", "s13", "s14", "s16", "s17", 
"s19", "s20", "s21", "ss10", "ss30", "ss32", "ss38"), class = "factor"), 
    wm = c(12L, 11L, 5L, 10L, 9L, 10L, 11L, 11L, 12L, 0L, 9L, 
    11L, 12L, 9L, 8L, 11L, 10L, 7L, 9L, 9L, 11L, 9L, 7L, 10L, 
    13L, 1L, 8L, 10L, 9L, 9L, 8L, 10L, 8L, 10L, 9L, 11L, 7L, 
    8L, 10L, 8L, 8L, 10L, 9L, 9L, 8L, 10L, 11L, 9L, 8L, 10L, 
    9L, 6L, 10L, 11L, 6L, 11L, 9L, 10L, 7L, 8L, 8L, 12L, 8L, 
    10L, 11L, 10L, 6L, 12L, 8L, 9L, 14L, 12L, 6L, 10L)), .Names = c("participant", 
"wm"), class = "data.frame", row.names = c(NA, -74L))

# remove unwanted participants
wm <- wm_df[!wm_df$participant %in% remove, ]

# set wm dataframe as a named vector 
wm_vector <- setNames(as.character(wm$wm), wm$participant)

# subset the wm vector using the participant column of df as a 
# way to look up the wm values in the vector and assign the values 
# to a new 'wm' column 
hls_subset_0_glmm$wm <- wm_vector[as.character(hls_subset_0_glmm$participant)]
hls_subset_0_glmm$wm <- as.numeric(hls_subset_0_glmm$wm)





hls_subset_0_glmm$group <- factor(hls_subset_0_glmm$group, levels = c("ss", "la", "hs"))

prop_0_mod_0 <- lmer(meanFix ~ 1 + 
                    (1 | participant), 
                    data = hls_subset_0_glmm, control = lmerControl(optimizer = 'bobyqa'))

prop_0_mod_group <- lmer(meanFix ~ 1 + group + 
                        (1 | participant), 
                        data = hls_subset_0_glmm, control = lmerControl(optimizer = 'bobyqa'))

prop_0_mod_cond <- lmer(meanFix ~ 1 + group + condition + 
                        (1 | participant), 
                        data = hls_subset_0_glmm, control = lmerControl(optimizer = 'bobyqa'))

prop_0_mod_full <- lmer(meanFix ~ 1 + group * condition + 
                      (1 | participant), 
                       data = hls_subset_0_glmm, control = lmerControl(optimizer = 'bobyqa'))

prop_0_mod_cov <- lmer(meanFix ~ 1 + group * condition + wm + 
                      (1 | participant), 
                       data = hls_subset_0_glmm, control = lmerControl(optimizer = 'bobyqa'))

prop_0_mod_cov_full <- lmer(meanFix ~ 1 + group * condition + wm + wm:condition +  
                      (1 | participant), 
                       data = hls_subset_0_glmm, control = lmerControl(optimizer = 'bobyqa'))


# anova(prop_0_mod_0, prop_0_mod_group, prop_0_mod_cond, prop_0_mod_full, prop_0_mod_cov, prop_0_mod_cov_full, test = "Chisq")

#        Df     AIC     BIC logLik deviance  Chisq Chi Df Pr(>Chisq)   
# object  3 -61.838 -52.846 33.919  -67.838                            
# ..1     5 -62.635 -47.649 36.317  -72.635 4.7969      2    0.09086 . 
# ..2     6 -69.365 -51.382 40.683  -81.365 8.7301      1    0.00313 **
# ..3     8 -65.437 -41.460 40.719  -81.437 0.0724      2    0.96447   
# ..4     9 -65.259 -38.284 41.629  -83.259 1.8213      1    0.17716   
# ..5    10 -63.572 -33.600 41.786  -83.572 0.3133      1    0.57564   

# summary(prop_0_mod_cond)

# Fixed effects:
#                      Estimate Std. Error        df t value Pr(>|t|)    
# (Intercept)           0.70494    0.03165 109.98000  22.271  < 2e-16 ***
# groupla              -0.08269    0.03751  71.00000  -2.204  0.03075 *  
# grouphs              -0.05349    0.03931  71.00000  -1.361  0.17785    
# conditionunstressed   0.09059    0.03029  73.00000   2.991  0.00379 ** 

# Relevel to test la vs hs 
hls_subset_0_glmm$group <- factor(hls_subset_0_glmm$group, levels = c("la", "hs",  "ss"))
hls_subset_0_glmm$condition <- factor(hls_subset_0_glmm$condition, levels = c("unstressed",  "stressed"))

summary(lmer(meanFix ~ 1 + group + condition + wm + (1 | participant), data = hls_subset_0_glmm, control = lmerControl(optimizer = 'bobyqa')))

# Fixed effects:
#                    Estimate Std. Error        df t value Pr(>|t|)    
# (Intercept)        0.793369   0.067631 77.380000  11.731  < 2e-16 ***
# grouphs            0.032189   0.037385 70.000000   0.861  0.39217    
# groupss            0.086075   0.037405 70.000000   2.301  0.02436 *  
# conditionstressed -0.090590   0.030288 73.000000  -2.991  0.00379 ** 
# wm                -0.009019   0.006828 70.000000  -1.321  0.19084    




# Calculate mean target fixation as a function of group, condition, 
# for each participant. We will plot the mean and calculate the 
# bootstrapped 95% confidence interval and plot it all. 
hls_subset_0 %>%
  na.omit(.) %>%
  ggplot(., aes(x = group, y = targetProp, 
                dodge = condition, color = condition,
                group = interaction(group, condition))) +
    geom_hline(yintercept = 0.5, color = "black", size = 0.75, 
               lty = 3) + 
    stat_summary(fun.data = 'mean_cl_boot', geom = 'errorbar', 
                 position = position_dodge(width = 0.5), 
                 width = 0.35, color = 'grey40') + 
    stat_summary(fun.y = 'mean', geom = 'point', size = 4,
                 position = position_dodge(width = 0.5)) + 
    # coord_cartesian(ylim = c(0, 1)) + 
    ylab('Target fixations') + xlab('Group') + 
    scale_x_discrete(labels = c('LA', 'HS', 'SS')) + 
    ggtitle('Mean target fixations as a function of group and target type') + 
    scale_color_brewer(palette = "Set1", name = '', labels = c('Paroxytone', 'Oxytone')) + 
    theme_bw(base_size = 16, base_family = 'Times') -> stressTargetFixP2

# Looks good, save as .png file
# ggsave('stressTargetFixP2.png', plot = stressTargetFixP2, dpi = 600, device = "png", path = "./mySources/figs/stress/s2_adv_her_nat/eye_track")


























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
hls_gc_subset <- filter(hls_subset, binAdj >= -52 & binAdj <= 5)


# - Readjust time course 
#    - now we will make the time course positive, starting at 1
#    - to do this we will add the lowest 'binAdj' value to each 
#      bin, plus 1 (to avoid starting at 0)

# lsrl_gc_subset$binGC <- lsrl_gc_subset$binAdj + 51
hls_gc_subset$binGC <- hls_gc_subset$binAdj + 53


# - Now we add higher order polynomials for analyses

t <- poly(min(hls_gc_subset$binGC):max(hls_gc_subset$binGC), 3)
hls_gc_subset[, paste('ot', 1:3, sep = "")] <- t[hls_gc_subset$binGC, 1:3]

# add wm 
hls_gc_subset$wm <- wm_vector[as.character(hls_gc_subset$participant)]
hls_gc_subset$wm <- as.numeric(hls_gc_subset$wm)


# glimpse(hls_gc_subset)






####################################################
# - Question 1: Are the groups different from      #
#   each other in when they begin to fixate        #
#   on the target?                                 #
#     - test 3 groups at each level of 'condition' #
#     - hypothesis: SS has steeper slope for both  #
#       conditions                                 #
####################################################


# Set SS as reference level 
hls_gc_subset$group <- factor(hls_gc_subset$group, levels = c("ss", "hs", "la"))

# Set condition as factor and set conding to contrast 
hls_gc_subset$condition <- as.factor(hls_gc_subset$condition)

# contrasts(hls_gc_subset$condition)
# contrasts(hls_gc_subset$group)

hls_gc_subset$conditionSum <- C(hls_gc_subset$condition, sum)
# contrasts(hls_gc_subset$conditionSum)

hls_gc_subset$groupSum <- C(hls_gc_subset$group, sum)
# contrasts(hls_gc_subset$groupSum)

## @knitr ignore2

# load the models:
gc_mod_base    <- readRDS('./mySources/models/stress/s2_adv_her_nat/eye_track/gc_mod_base.rds')
gc_mod_group_0 <- readRDS('./mySources/models/stress/s2_adv_her_nat/eye_track/gc_mod_group_0.rds')
gc_mod_group_1 <- readRDS('./mySources/models/stress/s2_adv_her_nat/eye_track/gc_mod_group_1.rds')
gc_mod_group_2 <- readRDS('./mySources/models/stress/s2_adv_her_nat/eye_track/gc_mod_group_2.rds')
gc_mod_cond_0  <- readRDS('./mySources/models/stress/s2_adv_her_nat/eye_track/gc_mod_cond_0.rds')
gc_mod_cond_1  <- readRDS('./mySources/models/stress/s2_adv_her_nat/eye_track/gc_mod_cond_1.rds')
gc_mod_cond_2  <- readRDS('./mySources/models/stress/s2_adv_her_nat/eye_track/gc_mod_cond_2.rds')
gc_mod_full    <- readRDS('./mySources/models/stress/s2_adv_her_nat/eye_track/gc_mod_full.rds')
gc_mod_wm      <- readRDS('./mySources/models/stress/s2_adv_her_nat/eye_track/gc_mod_wm.rds')


# Base model 
if(F){
  gc_mod_base <- lmer(targetProp ~ (ot1+ot2) + 
                 ((ot1+ot2) | participant) + 
                 ((ot1+ot2) | target),
                 control = lmerControl(optimizer = 'bobyqa'), 
                 data = hls_gc_subset, REML = F)
  saveRDS(gc_mod_base, file = "./mySources/models/stress/s2_adv_her_nat/eye_track/gc_mod_base.rds", compress = 'xz')
}


# Add group effect on intercept 
if(F){
  gc_mod_group_0 <- lmer(targetProp ~ (ot1+ot2) + group + 
                    ((ot1+ot2) | participant) + 
                    ((ot1+ot2) | target), 
                    control = lmerControl(optimizer = 'bobyqa'), 
                    data = hls_gc_subset, REML = F)
  saveRDS(gc_mod_group_0, file = "./mySources/models/stress/s2_adv_her_nat/eye_track/gc_mod_group_0.rds", compress = 'xz')
}

# Add group effect on slope
if(F){
 gc_mod_group_1 <- lmer(targetProp ~ (ot1+ot2) + group + 
                   ot1:group + 
                   ((ot1+ot2) | participant) + 
                   ((ot1+ot2) | target),
                   control = lmerControl(optimizer = 'bobyqa'), 
                   data = hls_gc_subset, REML = F)
  saveRDS(gc_mod_group_1, file = "./mySources/models/stress/s2_adv_her_nat/eye_track/gc_mod_group_1.rds", compress = 'xz')
}

# Add group effect on quadratic poly 
if(F){
  gc_mod_group_2 <- lmer(targetProp ~ (ot1+ot2) + group + 
                    ot1:group + ot2:group + 
                    ((ot1+ot2) | participant) + 
                    ((ot1+ot2) | target),
                    control = lmerControl(optimizer = 'bobyqa'), 
                    data = hls_gc_subset, REML = F)
  saveRDS(gc_mod_group_2, file = "./mySources/models/stress/s2_adv_her_nat/eye_track/gc_mod_group_2.rds", compress = 'xz')
}


# Add condition effect on intercept 
if(F){
  gc_mod_cond_0 <- lmer(targetProp ~ (ot1+ot2) * group + condition + 
                   ((ot1+ot2) | participant) + 
                   ((ot1+ot2) | target),
                   control = lmerControl(optimizer = 'bobyqa'), 
                   data = hls_gc_subset, REML = F)
  saveRDS(gc_mod_cond_0, file = "./mySources/models/stress/s2_adv_her_nat/eye_track/gc_mod_cond_0.rds", compress = 'xz')
}

# Add condition effect on slope 
if(F){
  gc_mod_cond_1 <- lmer(targetProp ~ (ot1+ot2) * group + condition + 
                   ot1:condition + 
                   ((ot1+ot2) | participant) + 
                   ((ot1+ot2) | target),
                   control = lmerControl(optimizer = 'bobyqa'), 
                   data = hls_gc_subset, REML = F)
  saveRDS(gc_mod_cond_1, file = "./mySources/models/stress/s2_adv_her_nat/eye_track/gc_mod_cond_1.rds", compress = 'xz')
}

# Add condition effect on quadratic poly 
if(F){
  gc_mod_cond_2 <- lmer(targetProp ~ (ot1+ot2) * group + condition + 
                   ot1:condition + ot2:condition + 
                   ((ot1+ot2) | participant) + 
                   ((ot1+ot2) | target),
                   control = lmerControl(optimizer = 'bobyqa'), 
                   data = hls_gc_subset, REML = F)
  saveRDS(gc_mod_cond_2, file = "./mySources/models/stress/s2_adv_her_nat/eye_track/gc_mod_cond_2.rds", compress = 'xz')
}


# Include all interactions
if(F){
gc_mod_full <- lmer(targetProp ~ (ot1+ot2) * group * condition + 
               ((ot1+ot2) | participant) + 
               ((ot1+ot2) | target),
               control = lmerControl(optimizer = 'bobyqa'), 
               data = hls_gc_subset, REML = F)
  saveRDS(gc_mod_full, file = "./mySources/models/stress/s2_adv_her_nat/eye_track/gc_mod_full.rds", compress = 'xz')
}

if(F){
gc_mod_wm <- lmer(targetProp ~ (ot1+ot2) * group * condition + wm + 
               ((ot1+ot2) | participant) + 
               ((ot1+ot2) | target),
               control = lmerControl(optimizer = 'bobyqa'), 
               data = hls_gc_subset, REML = F)
  saveRDS(gc_mod_full, file = "./mySources/models/stress/s2_adv_her_nat/eye_track/gc_mod_wm.rds", compress = 'xz')
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
      gc_mod_cond_0, 
      gc_mod_cond_1, 
      gc_mod_cond_2, 
      gc_mod_full,
      gc_mod_wm,  test = 'Chisq')

#        Df    AIC    BIC logLik deviance    Chisq Chi Df Pr(>Chisq)     
# object 16 119906 120056 -59937   119874                                # base model 
# ..1    18 119909 120079 -59937   119873   0.3100      2  0.8564336     # add group effect on intercept 
# ..2    20 119902 120090 -59931   119862  11.0322      2  0.0040214 **  # add group effect on slope
# ..3    22 119905 120112 -59931   119861   1.1633      2  0.5589720     # add group effect on quadratic poly 
# ..4    23 119907 120123 -59931   119861   0.0392      1  0.8430199     # add cond effect on intercept
# ..5    24 119898 120124 -59925   119850  11.2161      1  0.0008109 *** # add cond effect on slope 
# ..6    25 119895 120130 -59922   119845   5.2366      1  0.0221171 *   # add cond effect on quadratic poly
# ..7    31 119788 120079 -59863   119726 118.9164      6  < 2.2e-16 *** # full model 
# ..8    32 119787 120088 -59862   119723   2.5104      1  0.1130945     # add wm

# summary(gc_mod_full)
# Fixed effects:
#                                   Estimate Std. Error         df t value  Pr(>|t|)    
# (Intercept)                      5.120e-01  2.511e-02  8.900e+01  20.390   < 2e-16 ***
# ot1                              6.967e-01  1.119e-01  9.800e+01   6.227  1.18e-08 ***
# ot2                              5.378e-01  7.306e-02  1.170e+02   7.362  2.74e-11 ***
# grouphs                          1.903e-03  2.523e-02  7.800e+01   0.075  0.940071    
# groupla                          2.942e-03  2.413e-02  7.800e+01   0.122  0.903257    
# conditionunstressed              5.430e-02  2.569e-02  3.500e+01   2.114  0.041764 *  
# ot1:grouphs                     -3.025e-01  1.186e-01  8.300e+01  -2.550  0.012617 *  
# ot1:groupla                     -5.618e-01  1.138e-01  8.500e+01  -4.935  3.96e-06 ***
# ot2:grouphs                     -1.336e-01  8.637e-02  9.300e+01  -1.546  0.125373    
# ot2:groupla                     -1.238e-01  8.330e-02  9.700e+01  -1.486  0.140404    
# ot1:conditionunstressed          2.832e-01  1.141e-01  4.300e+01   2.483  0.017017 *  
# ot2:conditionunstressed         -2.668e-01  7.251e-02  7.400e+01  -3.680  0.000439 ***
# grouphs:conditionunstressed     -1.119e-02  8.241e-03  8.930e+04  -1.358  0.174495    
# groupla:conditionunstressed     -2.804e-02  8.186e-03  8.929e+04  -3.426  0.000614 ***
# ot1:grouphs:conditionunstressed -7.199e-02  6.276e-02  8.930e+04  -1.147  0.251314    
# ot1:groupla:conditionunstressed  4.191e-01  6.234e-02  8.930e+04   6.723  1.79e-11 ***
# ot2:grouphs:conditionunstressed  2.410e-01  6.275e-02  8.932e+04   3.841  0.000123 ***
# ot2:groupla:conditionunstressed  9.447e-02  6.234e-02  8.931e+04   1.515  0.129673    




# create new df including the fitted model 
data.comp <- data.frame(na.omit(hls_gc_subset), 
                        GCA_Full = fitted(gc_mod_full))
# glimpse(data.comp)


suffix_area <- data.frame(x = 53:58, y = Inf)

condition_namesModHLS <- c(
                    `stressed` = "Paroxytone", 
                    `unstressed` = "Oxytone"
                    )

data.comp %>% 
  ggplot(., aes(x = binGC, y = targetProp, color = group)) + 
  facet_grid(. ~ condition, labeller = as_labeller(condition_namesModHLS)) + 
  geom_area(data = suffix_area, aes(x = x, y = y), inherit = FALSE, alpha = 0.3, fill = 'lightcyan2') +
  stat_summary(fun.data = mean_se, geom = 'errorbar', 
               show.legend = FALSE, size = 0.1) +
  stat_summary(fun.y = mean, geom = 'point', size = 0.2) + 
  stat_summary(aes(y = GCA_Full, color = group), fun.y = mean, geom = 'line', size = 0.4) + 
  xlab("Adjusted time course") +
  ylab("Target fixations") +
  coord_cartesian(ylim = c(0.0, 1.0)) + 
  scale_x_continuous(breaks = c(1, 53), labels = c("Approx.\ntarget\nonset", "Target\nsyllable\nonset")) + 
  scale_color_brewer(palette = "Set1", name = "", labels = c("SS", "HS", "LA")) + 
  theme_bw(base_size = 16, base_family = "Times New Roman") -> stressGCAfullMod

# ggsave('stressGCAfullMod.png', plot = stressGCAfullMod, dpi = 600, device = "png", path = "./mySources/figs/stress/s2_adv_her_nat/eye_track")















# RUN GCA WITHOUT SS GROUP 

hls_gc_hs_la <- hls_gc_subset[hls_gc_subset$group != 'ss', ]

hls_gc_hs_la$group <- droplevels(hls_gc_hs_la$group)

hls_gc_hs_la$groupSum <- C(hls_gc_hs_la$group, sum)

contrasts(hls_gc_hs_la$groupSum)
contrasts(hls_gc_hs_la$conditionSum)



# load the models:
# gc_mod2_base    <- readRDS('./mySources/models/stress/s2_adv_her_nat/eye_track/gc_mod2_base.rds')
# gc_mod2_group_0 <- readRDS('./mySources/models/stress/s2_adv_her_nat/eye_track/gc_mod2_group_0.rds')
# gc_mod2_group_1 <- readRDS('./mySources/models/stress/s2_adv_her_nat/eye_track/gc_mod2_group_1.rds')
# gc_mod2_group_2 <- readRDS('./mySources/models/stress/s2_adv_her_nat/eye_track/gc_mod2_group_2.rds')
# gc_mod2_cond_0  <- readRDS('./mySources/models/stress/s2_adv_her_nat/eye_track/gc_mod2_cond_0.rds')
# gc_mod2_cond_1  <- readRDS('./mySources/models/stress/s2_adv_her_nat/eye_track/gc_mod2_cond_1.rds')
# gc_mod2_cond_2  <- readRDS('./mySources/models/stress/s2_adv_her_nat/eye_track/gc_mod2_cond_2.rds')
# gc_mod2_full    <- readRDS('./mySources/models/stress/s2_adv_her_nat/eye_track/gc_mod2_full.rds')
# gc_mod2_wm      <- readRDS('./mySources/models/stress/s2_adv_her_nat/eye_track/gc_mod2_wm.rds')


# Base model 
if(F){
  gc_mod2_base <- lmer(targetProp ~ (ot1+ot2) + 
                 ((ot1+ot2) | participant) + 
                 ((ot1+ot2) | target),
                 control = lmerControl(optimizer = 'bobyqa'), 
                 data = hls_gc_hs_la, REML = F)
  saveRDS(gc_mod2_base, file = "./mySources/models/stress/s2_adv_her_nat/eye_track/gc_mod2_base.rds", compress = 'xz')
}


# Add group effect on intercept 
if(F){
  gc_mod2_group_0 <- lmer(targetProp ~ (ot1+ot2) + group + 
                    ((ot1+ot2) | participant) + 
                    ((ot1+ot2) | target), 
                    control = lmerControl(optimizer = 'bobyqa'), 
                    data = hls_gc_hs_la, REML = F)
  saveRDS(gc_mod2_group_0, file = "./mySources/models/stress/s2_adv_her_nat/eye_track/gc_mod2_group_0.rds", compress = 'xz')
}

# Add group effect on slope
if(F){
 gc_mod2_group_1 <- lmer(targetProp ~ (ot1+ot2) + group + 
                   ot1:group + 
                   ((ot1+ot2) | participant) + 
                   ((ot1+ot2) | target),
                   control = lmerControl(optimizer = 'bobyqa'), 
                   data = hls_gc_hs_la, REML = F)
  saveRDS(gc_mod2_group_1, file = "./mySources/models/stress/s2_adv_her_nat/eye_track/gc_mod2_group_1.rds", compress = 'xz')
}

# Add group effect on quadratic poly 
if(F){
  gc_mod2_group_2 <- lmer(targetProp ~ (ot1+ot2) + group + 
                    ot1:group + ot2:group + 
                    ((ot1+ot2) | participant) + 
                    ((ot1+ot2) | target),
                    control = lmerControl(optimizer = 'bobyqa'), 
                    data = hls_gc_hs_la, REML = F)
  saveRDS(gc_mod2_group_2, file = "./mySources/models/stress/s2_adv_her_nat/eye_track/gc_mod2_group_2.rds", compress = 'xz')
}


# Add condition effect on intercept 
if(F){
  gc_mod2_cond_0 <- lmer(targetProp ~ (ot1+ot2) * group + condition + 
                   ((ot1+ot2) | participant) + 
                   ((ot1+ot2) | target),
                   control = lmerControl(optimizer = 'bobyqa'), 
                   data = hls_gc_hs_la, REML = F)
  saveRDS(gc_mod2_cond_0, file = "./mySources/models/stress/s2_adv_her_nat/eye_track/gc_mod2_cond_0.rds", compress = 'xz')
}

# Add condition effect on slope 
if(F){
  gc_mod2_cond_1 <- lmer(targetProp ~ (ot1+ot2) * group + condition + 
                   ot1:condition + 
                   ((ot1+ot2) | participant) + 
                   ((ot1+ot2) | target),
                   control = lmerControl(optimizer = 'bobyqa'), 
                   data = hls_gc_hs_la, REML = F)
  saveRDS(gc_mod2_cond_1, file = "./mySources/models/stress/s2_adv_her_nat/eye_track/gc_mod2_cond_1.rds", compress = 'xz')
}

# Add condition effect on quadratic poly 
if(F){
  gc_mod2_cond_2 <- lmer(targetProp ~ (ot1+ot2) * group + condition + 
                   ot1:condition + ot2:condition + 
                   ((ot1+ot2) | participant) + 
                   ((ot1+ot2) | target),
                   control = lmerControl(optimizer = 'bobyqa'), 
                   data = hls_gc_hs_la, REML = F)
  saveRDS(gc_mod2_cond_2, file = "./mySources/models/stress/s2_adv_her_nat/eye_track/gc_mod2_cond_2.rds", compress = 'xz')
}


# Include all interactions
if(F){
gc_mod2_full <- lmer(targetProp ~ (ot1+ot2) * group * condition + 
               ((ot1+ot2) | participant) + 
               ((ot1+ot2) | target),
               control = lmerControl(optimizer = 'bobyqa'), 
               data = hls_gc_hs_la, REML = F)
  saveRDS(gc_mod2_full, file = "./mySources/models/stress/s2_adv_her_nat/eye_track/gc_mod2_full.rds", compress = 'xz')
}

# add wm 
if(F){
gc_mod2_wm <- lmer(targetProp ~ (ot1+ot2) * group * condition + wm +
               ((ot1+ot2) | participant) + 
               ((ot1+ot2) | target),
               control = lmerControl(optimizer = 'bobyqa'), 
               data = hls_gc_hs_la, REML = F)
  saveRDS(gc_mod2_full, file = "./mySources/models/stress/s2_adv_her_nat/eye_track/gc_mod2_wm.rds", compress = 'xz')
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
#       gc_mod2_full, 
#       gc_mod2_wm, test = 'Chisq')

#        Df   AIC   BIC logLik deviance    Chisq Chi Df Pr(>Chisq)    
# object 16 91711 91857 -45839    91679                               # base model 
# ..1    17 91712 91867 -45839    91678  0.4053      1  0.5243727     # add group effect on intercept 
# ..2    18 91714 91879 -45839    91678  0.0110      1  0.9165961     # add group effect on slope
# ..3    19 91715 91889 -45839    91677  0.6345      1  0.4257200     # add group effect on quadratic poly 
# ..4    20 91717 91900 -45839    91677  0.0072      1  0.9324908     # add cond effect on intercept
# ..5    21 91708 91900 -45833    91666 11.7078      1  0.0006224 *** # add cond effect on slope 
# ..6    22 91708 91909 -45832    91664  1.9380      1  0.1638868     # add cond effect on quadratic poly
# ..7    25 91621 91850 -45786    91571 92.6430      3  < 2.2e-16 *** # full model 
# ..8    26 91621 91859 -45785    91569  2.0763      1  0.1496060     # add wm



# summary(gc_mod2_full)

# Fixed effects:
#                                   Estimate Std. Error         df t value  Pr(>|t|)    
# (Intercept)                      5.135e-01  2.636e-02  7.300e+01  19.477  < 2e-16 ***
# ot1                              3.978e-01  1.194e-01  7.300e+01   3.332  0.00135 ** 
# ot2                              4.063e-01  7.675e-02  7.600e+01   5.293 1.12e-06 ***
# groupla                          1.389e-03  2.411e-02  5.300e+01   0.058  0.95429    
# conditionunstressed              4.338e-02  2.792e-02  3.300e+01   1.554  0.12992    
# ot1:groupla                     -2.592e-01  1.142e-01  5.400e+01  -2.270  0.02718 *  
# ot2:groupla                      9.946e-03  8.340e-02  5.900e+01   0.119  0.90548    
# ot1:conditionunstressed          2.111e-01  1.259e-01  3.500e+01   1.677  0.10247    
# ot2:conditionunstressed         -2.639e-02  7.613e-02  4.200e+01  -0.347  0.73065    
# groupla:conditionunstressed     -1.768e-02  7.202e-03  6.814e+04  -2.455  0.01409 *  
# ot1:groupla:conditionunstressed  4.895e-01  5.485e-02  6.814e+04   8.924  < 2e-16 ***
# ot2:groupla:conditionunstressed -1.455e-01  5.484e-02  6.816e+04  -2.652  0.00799 ** 


