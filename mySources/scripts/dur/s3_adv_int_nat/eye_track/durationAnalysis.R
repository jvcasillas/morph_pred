#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Morphosyntactic predictability                                              #
# Duration analyses                                                           #
# 05/20/2017                                                                  #
# Script 2                                                                    #
# Study 3 (interpreters)                                                      #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# clean working directory
rm(list = ls(all = TRUE))

# Set working directory
# setwd("~/Desktop/morph_pred/")
setwd("~/academia/research/in_progress/morph_pred")

library(tidyverse); library(lme4); library(lmerTest); 
library(gridExtra); library(cowplot); library(broom)

# read data
 df_dur <- read_csv("./mySources/data/clean/durationBIN10Clean.csv")
# df_dur <- read_csv("./mySources/data/clean/durationBIN20Clean.csv")
# df_dur50 <- read_csv("./mySources/data/clean/durationBIN50Clean.csv")

# temp1 <- arrange(df_dur, participant, group, bin, exp, wavID, verb, startsentence, word2_c1v1, word3_c1v1, word4_20msafterv1, word4_c1, word4_c1v1, word5, word6, word7, end_sentence, target, targetside, BIN_DURATION, BIN_END_TIME, BIN_SAMPLE_COUNT, BIN_START_TIME, EYE_TRACKED, IA_1_ID, IA_2_ID, targetCount, distractorCount, targetProp, distractorProp, condition, condToken, targetOffset, binN, binAdj)
# temp2 <- arrange(df_dur2, participant, group, bin, exp, wavID, verb, startsentence, word2_c1v1, word3_c1v1, word4_20msafterv1, word4_c1, word4_c1v1, word5, word6, word7, end_sentence, target, targetside, BIN_DURATION, BIN_END_TIME, BIN_SAMPLE_COUNT, BIN_START_TIME, EYE_TRACKED, IA_1_ID, IA_2_ID, targetCount, distractorCount, targetProp, distractorProp, condition, condToken, targetOffset, binN, binAdj)

# identical(temp1, temp2[,-33])
# all.equal(temp1, temp2[,-33], ignore.row.order = TRUE)

# set variables and check it again
df_dur$targetProp <- gsub(",", ".", paste(df_dur$targetProp))
df_dur$distractorProp <- gsub(",", ".", paste(df_dur$distractorProp))

df_dur50$targetProp <- gsub(",", ".", paste(df_dur50$targetProp))
df_dur50$distractorProp <- gsub(",", ".", paste(df_dur50$distractorProp))


# there target word is word 4, so we can cut all data from word 6 on
# What is the highest target offset? (mes/conseguir = 2786)
df_dur %>%
  group_by(target, verb) %>%
  summarise(., offsets = unique(targetOffset)) %>%
  as.data.frame(.)

# what is the highest w6 start? (meses/conseguir: 3011)
df_dur %>%
  group_by(target, verb) %>%
  summarise(., w6start = unique(word6)) %>%
  as.data.frame(.)
# So we can use meses/conseguir:3011 as the cutoff point. 
# What is the bin number corresponding to 3011
high_w6_bin <- round(3011 / 10) + 1


df_dur <- df_dur %>%
  # filter(., bin <= high_w6_bin) %>%
  mutate(., 
         targetProp = as.numeric(targetProp), 
         distractorProp = as.numeric(distractorProp))

# Glimpse of data structure
glimpse(df_dur)





# Token realignment
# - 1st we calculate the lowest max 
# - 2nd we calculate the highest min 
# These values will become the extremes for our time course 

# calculate lowest adjusted bin max in order to determine 
# lowest possible upper bound for the time course
maxes <- df_dur %>%
  group_by(., target, verb, condToken) %>%
  summarize(max = max(binAdj)) %>% 
  as.data.frame(.)
binAdjMaxMin <- min(maxes$max)

# calculate highest adjusted bin minimum in order to determin 
# highest possible lower bound for the time course
mins <- df_dur %>%
  group_by(., target, verb, condToken) %>%
  summarize(min = min(binAdj)) %>% 
  as.data.frame(.)
binAdjMinMax <- max(mins$min)

# subset data based on new ranges 
# dur_short <- df_dur
dur_short <- df_dur %>% filter(., binAdj <= 100 & binAdj >= binAdjMinMax)

# create new adjusted variable that ranges from 1 to max
dur_short$binREadj <- (dur_short$binAdj - binAdjMinMax) + 1


# check it 
glimpse(dur_short)




## @knitr binAdjustments

# Bin adjustment plots
# - First we calculate the target onset for each item 
# - Next we calculate the adjusted onset 

# Where does target end in the shortened time course starts at word 
# prior to target word)? 
offsets <- dur_short %>%
  group_by(., target, verb, condToken) %>%
  summarize(., offset = unique(targetOffset) / 10) %>%
  mutate(., condToken = as.factor(condToken)) %>%
  as.data.frame(.)

# Center time course so that offset = 0
offsetAdj <- dur_short %>%
  group_by(., target, verb, condToken) %>%
  summarize(., offset = (unique(targetOffset) / 10) - 
                        (unique(targetOffset) / 10)) %>%
  mutate(., condToken = as.factor(condToken)) %>%
  as.data.frame(.)

# Where does the target word begin in the time course?
twOnsets <- dur_short %>%
  group_by(., target, verb, condToken) %>%
  summarize(., twOnset = unique(word4_c1v1) / 10) %>%
  mutate(., condToken = as.factor(condToken)) %>%
  as.data.frame(.)

# Adjust to centered time course 
twOnsetAdj <- cbind(offsets, twOnsets[, 4])
twOnsetAdj <- mutate(twOnsetAdj, twOnset = twOnsets[, 4], diff = offset - twOnset, twOnsetAdj = 0 - diff)


# Plot unadjusted durations 
df_dur %>%
  ggplot(., aes(x = bin, y = target)) + 
  geom_point(alpha = 0.2, size = 0.2) + 
  geom_point(data = offsets, aes(x = offset, y = target), color = 'red') +
  geom_point(data = twOnsets, aes(x = twOnset, y = target), color = 'blue') + 
  facet_grid(condToken ~ .) +
  xlab("Time course") + ylab("Items") + 
  theme_bw(base_size = 12, base_family = "Times") -> onsetP


# Plot adjusted durations 
df_dur %>%
  ggplot(., aes(x = binAdj, y = target)) + 
  geom_point(alpha = 0.2, size = 0.2) +
  geom_point(data = offsetAdj, aes(x = offset, y = target), color = 'red') + 
  geom_point(data = twOnsetAdj, aes(x = twOnsetAdj, y = target), color = 'blue') + 
  geom_vline(xintercept = binAdjMinMax, color = 'red') + 
  geom_vline(xintercept = binAdjMaxMin, color = 'red') + 
  scale_y_discrete(position = "right") +
  facet_grid(condToken ~ .) +
  xlab("Adjusted time course") + ylab("Items") + 
  theme_bw(base_size = 12, base_family = "Times") -> onsetAdjP

durBinPlots <- plot_grid(onsetP, onsetAdjP, ncol = 2)

# Bin adjustments look good, save plots as .png file 
# ggsave('durBinPlots.png', plot = durBinPlots, dpi = 300, device = "png", path = "./mySources/figs/dur/general")









# Time sequence for plots 

# Create vector of times for canta 
colTimes <- c(100, 296, 1178, 1890, 2030, 2333, 2557, 2850)
# Create vector of labels 
colSeq   <- c("el", "cocinero", "escogio", "c", "o", "l", "para", "el")
# add labels as names arg for vector colTimes
names(colTimes) <- colSeq
# adjust bins 
colAdj <- (colTimes / 10) - (colTimes[6] / 10)
# turn in into a dataframe 
col_df_temp <- data.frame(group = rep(c("int", "la", "ss"), each = 8), 
                          step  = 1:8, 
                          condition = "monosyllabic", 
                          binAdj     = colAdj, 
                          binN     = colAdj, 
                          targetProp = 0.05, 
                          text  = names(colAdj))


colesTimes <- c(100, 309, 1270, 1987, 2117, 2347, 2467, 2559)
# Create vector of labels 
colesSeq   <- c("el", "cocinero", "escogio", "c", "o", "l", "e", "s")
# add labels as names arg for vector cantaTimes
names(colesTimes) <- colesSeq
# adjust bins 
colesAdj <- (colesTimes / 10) - (colesTimes[6] / 10)
# turn in into a dataframe 
coles_df_temp <- data.frame(group = rep(c("int", "la", "ss"), each = 8), 
                            step  = 1:8, 
                            condition = "bisyllabic", 
                            binAdj     = colesAdj, 
                            binN     = colesAdj, 
                            targetProp     = 0.05, 
                            text  = names(colesAdj))
colesEx <- rbind(col_df_temp, coles_df_temp)

condition_names <- c(
                    `monosyllabic` = "Monosyllabic",
                    `bisyllabic` = "Bisyllabic", 
                    `la` = "LA", 
                    `int` = "IN", 
                    `ss` = "SS" 
                    )

dur_short_subset <- dur_short %>% filter(., group %in% c('la', 'int', 'ss'))

# Relevel factor
dur_short_subset$condition <- factor(dur_short_subset$condition, levels = c("monosyllabic", "bisyllabic"))
dur_short_subset$group <- factor(dur_short_subset$group, levels = c("la", "int", "ss"))


dur_timecourse_p1 <- dur_short_subset %>%
  na.omit(.) %>% 
  filter(., !(target %in% c('mes', 'meses', 'sol', 'soles')), 
             binAdj <= 80 & binAdj >= -100) %>%
  group_by(., participant, binAdj, group, target) %>%
  summarize(., targetProp = mean(targetProp)) %>%
  ggplot(., aes(x = binAdj, y = targetProp)) + 
  facet_grid(group ~ condition, labeller = as_labeller(condition_names)) + 
  geom_vline(xintercept = 0, color = 'grey60') + 
  geom_hline(yintercept = 0.5, color = 'grey60') + 
  stat_summary(fun.data = mean_cl_normal, geom = 'errorbar', width = 0, 
               size = 0.1, show.legend = FALSE, color = 'darkgrey') +
  stat_summary(fun.y = mean, size = 0.5, geom = 'point', color = 'blue') + 
  stat_summary(data = na.omit(dur_short_subset[dur_short_subset$binAdj <= 80 & 
                                               dur_short_subset$binAdj >= -100, ]), aes(x = binAdj, y = distractorProp), 
               fun.data = mean_cl_normal, geom = 'errorbar', width = 0, 
               size = 0.1, color = 'darkgrey') + 
  stat_summary(data = na.omit(dur_short_subset[dur_short_subset$binAdj <= 80 & 
                                               dur_short_subset$binAdj >= -100, ]), aes(x = binAdj, y = distractorProp), 
               fun.y = mean, geom = 'point', size = 0.5, color = 'darkred') + 
  ylab('Proportion of fixations') + 
  xlab('Adjusted time course') + 
  coord_cartesian(ylim = c(0, 1)) + 
  geom_text(data = colesEx[colesEx$step >= 3 & colesEx$step <= 9, ], 
            aes(label = text), color = 'black') + 
  theme_bw(base_size = 16, base_family = "Times") 

# ggsave('dur_timecourse_p1.png', plot = dur_timecourse_p1, dpi = 600, device = "png", path = "./mySources/figs/dur/s3_adv_int_nat/eye_track", width = 10, height = 7, units = "in")
















#######################################
# t-tests                             #
# - Question: can they predict after  #
#   hearing 20ms of target word       #
#   vowel?                            #
# - 'they' refers to each group of    #
#   participants for each type of     #
#   word (mono, di)                   #
# - This analysis does not compare    #
#   groups or conditions              #
#######################################

# We want to analyze proportion of target gaze at target onset 
# so we need to make a subset of the data that only uses the 
# target onset bin (binAdj = 0)

dur_0 <- dur_short_subset %>% filter(., binAdj == 0)

# Quick and dirty mean of target fixations as a function of 
# group and condition (mono, di)

dur_0 %>%
  na.omit() %>%
  filter(., (!target %in% c('mes', 'meses', 'sol', 'soles'))) %>%
  group_by(., group, condition) %>%
  summarise(., meanFix = mean(targetProp), sdFix = sd(targetProp)) 

# We will test this for each group in each condition (mono, di)
# using a one-sided t-test. Specifically, we are testing the 
# hypothesis that the proportion of looks is greater than 
# chance (50%). 
# - H0: u = 0.5
# - Ha: u > 0.5
# The generic code is: t.test(myVector, alternative = "greater", my = 0.50, conf.level = 0.95)

dur_ttest <- dur_0 %>%
  na.omit() %>%
  filter(., (!target %in% c('mes', 'meses', 'sol', 'soles'))) %>%
  group_by(., group, condition, participant) %>%
  summarise(., meanFix = mean(targetProp)) %>% 
  do(tidy(t.test(.$meanFix, alternative = "greater", mu = 0.50, conf.level = 0.95)))

# Convert pvalues from scientific notation
dur_ttest$p.value <- format(dur_ttest$p.value, scientific = F)
dur_ttest$sig <- "N.S."
dur_ttest[dur_ttest$p.value <= 0.05, 'sig'] <- "*"

# Print results
print(as.data.frame(dur_ttest[, c(1:7, 11)]))

# group    condition  estimate  statistic     p.value parameter  conf.low  sig
#    la monosyllabic 0.5163978  0.3829793 0.352217712        30 0.4437270 N.S.
#    la   bisyllabic 0.5962366  2.0264360 0.025845343        30 0.5156327    *
#   int monosyllabic 0.6091667  2.1323141 0.030894311         9 0.5153180    *
#   int   bisyllabic 0.7441667  3.9741222 0.001617207         9 0.6315418    *
#    ss monosyllabic 0.4250000 -1.4497002 0.919375441        22 0.3361638 N.S.
#    ss   bisyllabic 0.5760870  1.3789877 0.090879214        22 0.4813419 N.S.








# We will plot the models 
# This will be almost exactly the same as 
# the previous plot, but it will use the 
# confidence interval from the test we 
# actually conducted
dur_ttest$group <- factor(dur_ttest$group, levels = c("la", "int", "ss"))
dur_ttest$condition <- factor(dur_ttest$condition, levels = c("monosyllabic", "bisyllabic"))


ggplot(dur_ttest, aes(x = group, y = estimate, color = condition, 
                      group = interaction(group, condition), dodge = condition)) +
    geom_hline(yintercept = 0.5, lty = 3) + 
    geom_linerange(aes(ymin = conf.low, ymax = estimate), color = 'grey40',
                     position = position_dodge(width = 0.75), size = 1) +
    geom_point(position = position_dodge(width = 0.75), size = 4) +
    ylim(0, 1.0) + ylab('Target fixations') + xlab('Group') + 
    ggtitle('Mean target fixations and lower-bound 95% confidence interval') + 
    scale_color_brewer(palette = "Set1", name = '', labels = c('Monosyllabic', 'Bisyllabic')) + 
    scale_x_discrete(labels = c('LA', 'INT', 'SS')) + 
    theme_bw(base_size = 16, base_family = 'Times') -> durTargetFixP2

# Looks good, save as .png file. 
# ggsave('durTargetFixP2.png', plot = durTargetFixP2, dpi = 600, device = "png", path = "./mySources/figs/dur/s3_adv_int_nat/eye_track", height = 5, width = 9, unit = "in")





dur_0_prop <- dur_0 %>%
  na.omit() %>%
  filter(., (!target %in% c('mes', 'meses', 'sol', 'soles'))) %>%
  group_by(., group, condition, participant) %>%
  summarise(., meanFix = mean(targetProp))  

# Add working memory covariate
# NEED DATA FOR INT
# read data
my_wm <- read.csv("./mySources/data/raw/my_wm.csv", header = TRUE, quote = "")



prop_0_mod_0 <- lmer(meanFix ~ 1 + 
                    (1 | participant), 
                    data = dur_0_prop, REML = F, control = lmerControl(optimizer = 'bobyqa'))

prop_0_mod_group <- lmer(meanFix ~ 1 + group + 
                        (1 | participant), 
                        data = dur_0_prop, REML = F, control = lmerControl(optimizer = 'bobyqa'))

prop_0_mod_cond <- lmer(meanFix ~ 1 + group + condition + 
                        (1 | participant), 
                        data = dur_0_prop, REML = F, control = lmerControl(optimizer = 'bobyqa'))

prop_0_mod_full <- lmer(meanFix ~ 1 + group * condition + 
                      (1 | participant), 
                       data = dur_0_prop, REML = F, control = lmerControl(optimizer = 'bobyqa'))

anova(prop_0_mod_0, prop_0_mod_group, prop_0_mod_cond, prop_0_mod_full, test = "Chisq")

#        Df     AIC    BIC  logLik deviance  Chisq Chi Df Pr(>Chisq)   
# object  3 16.8573 25.413 -5.4286  10.8573                            
# ..1     5 13.8834 28.144 -1.9417   3.8834 6.9739      2   0.030594 * 
# ..2     6  8.7924 25.905  1.6038  -3.2076 7.0910      1   0.007747 **
# ..3     8 12.1584 34.975  1.9208  -3.8416 0.6340      2   0.728341   

summary(prop_0_mod_cond)

# Fixed effects:
#                      Estimate Std. Error        df t value Pr(>|t|)    
# (Intercept)           0.49929    0.03697 128.00000  13.504  < 2e-16 ***
# groupint              0.12035    0.06145 128.00000   1.959  0.05234 .  
# groupss              -0.05577    0.04650 128.00000  -1.199  0.23258    
# conditionbisyllabic   0.11406    0.04224 128.00000   2.700  0.00787 ** 

# Calculate mean target fixation as a function of group, condition, 
# for each participant. We will plot the mean and calculate the 
# bootstrapped 95% confidence interval and plot it all. 
dur_0_prop %>%
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
    scale_x_discrete(labels = c('LA', 'IN', 'SS')) + 
    ggtitle('Mean target fixations as a function of group and target type') + 
    scale_color_brewer(palette = "Set1", name = '', labels = c('Monosyllabic', 'Bisyllabic')) + 
    theme_bw(base_size = 16, base_family = 'Times') -> durTargetFixP3

# ggsave('durTargetFixP3.png', plot = durTargetFixP3, dpi = 600, device = "png", path = "./mySources/figs/dur/s3_adv_int_nat/eye_track", height = 5, width = 9, unit = "in")















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
#   difference between mono/bisyl                  #
#   items?                                         #
#     - test mono vs. bi for each group            #
#     - hypothesis: steeper slope/earlier break in #
#       mono condition                             #
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

#    - the lowest target word onset is: 'soles/descubrir' @ -35
#    - thus we can use -35 binAdj as the starting point 
#      and be sure we are including the entire target word
#    - we can also do a little higher to lighten the models (-30)

df_dur50 <- read_csv("./mySources/data/clean/durationBIN50Clean.csv")


df_dur50$targetProp <- gsub(",", ".", paste(df_dur50$targetProp))
df_dur50$distractorProp <- gsub(",", ".", paste(df_dur50$distractorProp))

df_dur50 <- df_dur50 %>%
  # filter(., bin <= high_w6_bin) %>%
  mutate(., 
         targetProp = as.numeric(targetProp), 
         distractorProp = as.numeric(distractorProp))


maxes <- df_dur50 %>%
  group_by(., target, verb, condToken) %>%
  summarize(max = max(binAdj)) %>% 
  as.data.frame(.)
binAdjMaxMin <- min(maxes$max)

# calculate highest adjusted bin minimum in order to determin 
# highest possible lower bound for the time course
mins <- df_dur50 %>%
  group_by(., target, verb, condToken) %>%
  summarize(min = min(binAdj)) %>% 
  as.data.frame(.)
binAdjMinMax <- max(mins$min)

# subset data based on new ranges 

dur_short <- df_dur50 %>% filter(., binAdj <= 100 & binAdj >= binAdjMinMax)

# create new adjusted variable that ranges from 1 to max
dur_short$binREadj <- (dur_short$binAdj - binAdjMinMax) + 1


dur_gc_subset <- filter(dur_short, binAdj >= -30 & binAdj <= 30, (!target %in% c('mes', 'meses', 'sol', 'soles')))


# - Readjust time course 
#    - now we will make the time course positive, starting at 1
#    - to do this we will add the lowest 'binAdj' value to each 
#      bin, plus 1 (to avoid starting at 0)

dur_gc_subset$binGC <- dur_gc_subset$binAdj + 31


# - Now we add higher order polynomials for analyses

t <- poly(min(dur_gc_subset$binGC):max(dur_gc_subset$binGC), 3)
dur_gc_subset[, paste('ot', 1:3, sep = "")] <- t[dur_gc_subset$binGC, 1:3]

# glimpse(dur_gc_subset)






####################################################
# - Question 1: Are the groups different from      #
#   each other in when they begin to fixate        #
#   on the target?                                 #
#     - test 3 groups at each level of 'condition' #
#     - hypothesis: SS has steeper slope for both  #
#       conditions                                 #
####################################################


# Set SS as reference level 
dur_gc_subset$group <- factor(dur_gc_subset$group, levels = c("ss", "la", "int"))

# Set condition as factor and set conding to contrast 
dur_gc_subset$condition <- as.factor(dur_gc_subset$condition)

# contrasts(dur_gc_subset$condition)
# contrasts(dur_gc_subset$group)

dur_gc_subset$conditionSum <- C(dur_gc_subset$condition, sum)
# contrasts(dur_gc_subset$conditionSum)

dur_gc_subset$groupSum <- C(dur_gc_subset$group, sum)
# contrasts(dur_gc_subset$groupSum)

## @knitr ignore2

# load the models:
gc_mod_base    <- readRDS('./mySources/models/dur/s3_adv_int_nat/eye_track/gc_mod_base.rds')
gc_mod_group_0 <- readRDS('./mySources/models/dur/s3_adv_int_nat/eye_track/gc_mod_group_0.rds')
gc_mod_group_1 <- readRDS('./mySources/models/dur/s3_adv_int_nat/eye_track/gc_mod_group_1.rds')
gc_mod_group_2 <- readRDS('./mySources/models/dur/s3_adv_int_nat/eye_track/gc_mod_group_2.rds')
gc_mod_cond_0  <- readRDS('./mySources/models/dur/s3_adv_int_nat/eye_track/gc_mod_cond_0.rds')
gc_mod_cond_1  <- readRDS('./mySources/models/dur/s3_adv_int_nat/eye_track/gc_mod_cond_1.rds')
gc_mod_cond_2  <- readRDS('./mySources/models/dur/s3_adv_int_nat/eye_track/gc_mod_cond_2.rds')
gc_mod_full    <- readRDS('./mySources/models/dur/s3_adv_int_nat/eye_track/gc_mod_full.rds')
gc_mod_cube    <- readRDS('./mySources/models/dur/s3_adv_int_nat/eye_track/gc_mod_cube.rds')


# Base model 
if(T){
  gc_mod_base <- lmer(targetProp ~ (ot1+ot2) + 
                 ((ot1+ot2) | participant) + 
                 ((ot1+ot2) | target),
                 control = lmerControl(optimizer = 'bobyqa'), 
                 data = dur_gc_subset, REML = F)
  saveRDS(gc_mod_base, file = "./mySources/models/dur/s3_adv_int_nat/eye_track/gc_mod_base.rds", compress = 'xz')
}


# Add group effect on intercept 
if(T){
  gc_mod_group_0 <- lmer(targetProp ~ (ot1+ot2) + group + 
                    ((ot1+ot2) | participant) + 
                    ((ot1+ot2) | target), 
                    control = lmerControl(optimizer = 'bobyqa'), 
                    data = dur_gc_subset, REML = F)
  saveRDS(gc_mod_group_0, file = "./mySources/models/dur/s3_adv_int_nat/eye_track/gc_mod_group_0.rds", compress = 'xz')
}

# Add group effect on slope
if(T){
 gc_mod_group_1 <- lmer(targetProp ~ (ot1+ot2) + group + 
                   ot1:group + 
                   ((ot1+ot2) | participant) + 
                   ((ot1+ot2) | target),
                   control = lmerControl(optimizer = 'bobyqa'), 
                   data = dur_gc_subset, REML = F)
  saveRDS(gc_mod_group_1, file = "./mySources/models/dur/s3_adv_int_nat/eye_track/gc_mod_group_1.rds", compress = 'xz')
}

# Add group effect on quadratic poly 
if(T){
  gc_mod_group_2 <- lmer(targetProp ~ (ot1+ot2) + group + 
                    ot1:group + ot2:group + 
                    ((ot1+ot2) | participant) + 
                    ((ot1+ot2) | target),
                    control = lmerControl(optimizer = 'bobyqa'), 
                    data = dur_gc_subset, REML = F)
  saveRDS(gc_mod_group_2, file = "./mySources/models/dur/s3_adv_int_nat/eye_track/gc_mod_group_2.rds", compress = 'xz')
}


# Add condition effect on intercept 
if(T){
  gc_mod_cond_0 <- lmer(targetProp ~ (ot1+ot2) * group + condition + 
                   ((ot1+ot2) | participant) + 
                   ((ot1+ot2) | target),
                   control = lmerControl(optimizer = 'bobyqa'), 
                   data = dur_gc_subset, REML = F)
  saveRDS(gc_mod_cond_0, file = "./mySources/models/dur/s3_adv_int_nat/eye_track/gc_mod_cond_0.rds", compress = 'xz')
}

# Add condition effect on slope 
if(T){
  gc_mod_cond_1 <- lmer(targetProp ~ (ot1+ot2) * group + condition + 
                   ot1:condition + 
                   ((ot1+ot2) | participant) + 
                   ((ot1+ot2) | target),
                   control = lmerControl(optimizer = 'bobyqa'), 
                   data = dur_gc_subset, REML = F)
  saveRDS(gc_mod_cond_1, file = "./mySources/models/dur/s3_adv_int_nat/eye_track/gc_mod_cond_1.rds", compress = 'xz')
}

# Add condition effect on quadratic poly 
if(T){
  gc_mod_cond_2 <- lmer(targetProp ~ (ot1+ot2) * group + condition + 
                   ot1:condition + ot2:condition + 
                   ((ot1+ot2) | participant) + 
                   ((ot1+ot2) | target),
                   control = lmerControl(optimizer = 'bobyqa'), 
                   data = dur_gc_subset, REML = F)
  saveRDS(gc_mod_cond_2, file = "./mySources/models/dur/s3_adv_int_nat/eye_track/gc_mod_cond_2.rds", compress = 'xz')
}


# Include all interactions
if(T){
gc_mod_full <- lmer(targetProp ~ (ot1+ot2) * group * condition + 
               ((ot1+ot2) | participant) + 
               ((ot1+ot2) | target),
               control = lmerControl(optimizer = 'bobyqa'), 
               data = dur_gc_subset, REML = F)
  saveRDS(gc_mod_full, file = "./mySources/models/dur/s3_adv_int_nat/eye_track/gc_mod_full.rds", compress = 'xz')
}

if(T){
gc_mod_cube <- lmer(targetProp ~ (ot1+ot2+ot3) * group * condition + 
               ((ot1+ot2+ot3) | participant) + 
               ((ot1+ot2+ot3) | participant:condition) + 
               ((ot1+ot2+ot3) | target),
               control = lmerControl(optimizer = 'bobyqa'), 
               data = dur_gc_subset, REML = F)
  saveRDS(gc_mod_cube, file = "./mySources/models/dur/s3_adv_int_nat/eye_track/gc_mod_cube.rds", compress = 'xz')
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
      gc_mod_full, test = 'Chisq')

#        Df    AIC    BIC logLik deviance   Chisq Chi Df Pr(>Chisq)    
# object 16 117775 117925 -58871   117743                               # base model 
# ..1    18 117773 117943 -58869   117737  5.4043      2    0.06706 .   # add group effect on intercept 
# ..2    20 117774 117963 -58867   117734  2.9591      2    0.22774     # add group effect on slope
# ..3    22 117775 117983 -58866   117731  2.9042      2    0.23408     # add group effect on quadratic poly 
# ..4    23 117775 117991 -58864   117729  2.7054      1    0.10001     # add cond effect on intercept
# ..5    24 117777 118003 -58864   117729  0.0000      1    0.99527     # add cond effect on slope 
# ..6    25 117774 118010 -58862   117724  4.3255      1    0.03754 *   # add cond effect on quadratic poly
# ..7    31 117712 118004 -58825   117650 74.4951      6  4.876e-14 *** # full model 



# summary(gc_mod_full)

# Fixed effects:                                                             
#                                    Estimate Std. Error         df t value  Pr(>|t|)    
# (Intercept)                       5.645e-01  3.153e-02  2.900e+01  17.902   < 2e-16 ***
# ot1                               6.315e-01  1.854e-01  4.500e+01   3.407   0.00139 ** 
# ot2                               7.529e-01  1.563e-01  3.700e+01   4.818  2.44e-05 ***
# groupla                           9.209e-03  2.601e-02  6.500e+01   0.354   0.72447    
# groupint                          7.145e-02  3.578e-02  6.500e+01   1.997   0.05004 .  
# conditionbisyllabic               6.156e-02  3.521e-02  1.200e+01   1.748   0.10532    
# ot1:groupla                       1.908e-01  1.834e-01  7.000e+01   1.041   0.30153    
# ot1:groupint                     -3.216e-02  2.520e-01  6.900e+01  -0.128   0.89885    
# ot2:groupla                      -3.724e-01  1.457e-01  7.300e+01  -2.556   0.01266 *  
# ot2:groupint                     -3.488e-01  2.001e-01  7.200e+01  -1.743   0.08562 .  
# ot1:conditionbisyllabic          -1.887e-01  1.823e-01  1.300e+01  -1.035   0.31890    
# ot2:conditionbisyllabic          -5.675e-01  1.661e-01  1.400e+01  -3.416   0.00429 ** 
# groupla:conditionbisyllabic      -2.260e-02  6.978e-03  8.314e+04  -3.239   0.00120 ** 
# groupint:conditionbisyllabic      1.620e-02  9.507e-03  8.314e+04   1.703   0.08850 .  
# ot1:groupla:conditionbisyllabic   6.608e-02  7.351e-02  8.314e+04   0.899   0.36873    
# ot1:groupint:conditionbisyllabic  6.789e-01  1.002e-01  8.314e+04   6.777  1.23e-11 ***
# ot2:groupla:conditionbisyllabic   2.370e-01  7.351e-02  8.314e+04   3.224   0.00126 ** 
# ot2:groupint:conditionbisyllabic  1.580e-01  1.002e-01  8.314e+04   1.577   0.11471    


# create new df including the fitted model 
data.comp <- data.frame(na.omit(dur_gc_subset), 
                        GCA_Full = fitted(gc_mod_cube))
# glimpse(data.comp)


# I dont remember how I determined this
suffix_area <- data.frame(x = 31:61, y = Inf)

condition_namesMod <- c(
                    `monosyllabic` = "Monosyllabic", 
                    `bisyllabic` = "Bisyllabic" 
                    )

data.comp %>% 
  ggplot(., aes(x = binGC, y = targetProp, color = group, fill = group)) + 
  facet_grid(. ~ condition, labeller = as_labeller(condition_namesMod)) + 
  #geom_area(data = suffix_area, aes(x = x, y = y), inherit = FALSE, alpha = 0.3, fill = 'lightcyan2') +
  #stat_summary(fun.data = mean_se, geom = 'ribbon', 
  #             show.legend = FALSE, alpha = 0.2, color = NA) +
  #stat_summary(fun.y = mean, geom = 'point', size = 0.75) + 
  stat_summary(aes(y = GCA_Full, color = group), fun.y = mean, geom = 'line', size = 0.4) + 
  xlab("Adjusted time course") +
  ylab("Target fixations") +
  #coord_cartesian(ylim = c(0.0, 1.0)) + 
  scale_x_continuous(breaks = c(1, 31), labels = c("Approx.\ntarget\nonset", "Target\nsyllable\noffset")) + 
  scale_color_brewer(palette = "Set1", name = "", labels = c("SS", "LA", "IN")) + 
  theme_bw(base_size = 16, base_family = "Times New Roman") -> durGCAfullMod

ggsave('durGCAfullMod.png', plot = durGCAfullMod, dpi = 600, device = "png", 
         path = "./mySources/figs/dur/s3_int_adv_nat/eye_track", 
         height = 3.25, width = 10.5, unit = "in")

ggsave('durGCAfullMod.png', plot = durGCAfullMod, dpi = 600, device = "png", 
         path = "./mySources/figs/dur/s3_int_adv_nat/eye_track")









time_mat <- poly(sort(unique(data.comp$binGC)), 2) %>%
  polypoly::poly_rescale(1) %>%
  cbind(constant = 1, .)
round(time_mat, 2)

mono_coefs <- fixef(gc_mod_full)[1:3]
bisl_coefs <- mono_coefs + fixef(gc_mod_full)[c(6, 11, 12)]
bisl_coefs

set_colnames <- `colnames<-`

m_mono <- time_mat %*% diag(mono_coefs) %>%
  set_colnames(c("constant", "ot1", "ot2")) 

m_bisl <- time_mat %*% diag(bisl_coefs) %>%
  set_colnames(c("constant", "ot1", "ot2")) 

# Convince ourselves with an example
round(m_bisl, 2)


df_mono <- m_mono %>%
  polypoly::poly_melt() %>%
  tibble::add_column(Condition = "mono")

df_bisl <- m_bisl %>% 
  polypoly::poly_melt() %>%
  tibble::add_column(Condition = "bisl")

df_both <- bind_rows(df_bisl, df_mono) %>% 
  mutate(Condition = factor(Condition, c("mono", "bisl")))

ggplot(df_both) +
  aes(x = observation, y = value, color = Condition) +
  geom_line() + 
  facet_wrap("degree") + 
  theme_bw()