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

library(dplyr); library(tidyr); library(ggplot2); library(plotly)
library(lme4); library(lmerTest); library(gridExtra); library(cowplot)
library(broom)

# read data
# df_lex <- read_csv("./mySources/data/clean/lexicalBIN10Clean.csv")
df_lex <- read_csv("./mySources/data/clean/lexicalBIN10CleanNEW.csv")


# set variables and check it again
df_lex$targetProp <- gsub(",", ".", paste(df_lex$targetProp))
df_lex$distractorProp <- gsub(",", ".", paste(df_lex$distractorProp))

df_lex <- df_lex %>%
  filter(., group %in% c('la', 'int', 'ss') & bin <= 131) %>%
  mutate(., 
         targetProp = as.numeric(targetProp), 
         distractorProp = as.numeric(distractorProp))

# Glimpse of data structure
glimpse(df_lex)






# Token realignment
# - 1st we calculate the lowest max 
# - 2nd we calculate the highest min 
# These values will become the extremes for our time course 

# calculate lowest adjusted bin max in order to determine 
# lowest possible upper bound for the time course
maxes <- df_lex %>%
  group_by(., target, verb, condToken) %>%
  summarize(max = max(binAdj)) %>% 
  as.data.frame(.)
binAdjMaxMin <- min(maxes$max)

# calculate highest adjusted bin minimum in order to determin 
# highest possible lower bound for the time course
mins <- df_lex %>%
  group_by(., target, verb, condToken) %>%
  summarize(min = min(binAdj)) %>% 
  as.data.frame(.)
binAdjMinMax <- max(mins$min)

# subset data based on new ranges 
lex_short <- df_lex %>% filter(., binAdj <= binAdjMaxMin & binAdj >= binAdjMinMax)

# create new adjusted variable that ranges from 1 to max
lex_short$binREadj <- (lex_short$binAdj - binAdjMinMax) + 1


# check it 
glimpse(lex_short)




## @knitr binAdjustments

# Bin adjustment plots
# - First we calculate the target onset for each item 
# - Next we calculate the adjusted onset 

# Where does target suffix begin in the shortened time course 
offsets <- lex_short %>%
  group_by(., target, verb, condToken) %>%
  summarize(., offset = unique(targetOffset) / 10) %>%
#  summarize(., offset = (unique(targetOffset) - unique(word2_c1v1)) / 10) %>%
  mutate(., condToken = as.factor(condToken)) %>%
  as.data.frame(.)

# Center time course so that offset = 0
offsetAdj <- lex_short %>%
  group_by(., target, verb, condToken) %>%
  summarize(., offset = (unique(targetOffset) / 10) - 
                        (unique(targetOffset) / 10)) %>%
#   summarize(., offset = ((unique(targetOffset) - unique(word2_c1v1)) / 10) - 
#                         (unique(targetOffset) - unique(word2_c1v1)) / 10) %>%
  mutate(., condToken = as.factor(condToken)) %>%
  as.data.frame(.)

# Where does the target word begin in the time course?
twOnsets <- lex_short %>%
  group_by(., target, verb, condToken) %>%
  summarize(., twOnset = unique(word2_c1v1) / 10) %>%
  # summarize(., twOnset = (unique(word4_c1v1) - unique(word3_c1v1)) / 10) %>%
  mutate(., condToken = as.factor(condToken)) %>%
  as.data.frame(.)

# Adjust to centered time course 
twOnsetAdj <- cbind(offsets, twOnsets[, 4])
twOnsetAdj <- mutate(twOnsetAdj, twOnset = twOnsets[, 4], diff = offset - twOnset, twOnsetAdj = 0 - diff)


# Plot unadjusted durations 
df_lex %>%
  ggplot(., aes(x = bin, y = target)) + 
  geom_point(alpha = 0.2, size = 0.2) + 
  geom_point(data = offsets, aes(x = offset, y = target), color = 'red') +
  geom_point(data = twOnsets, aes(x = twOnset, y = target), color = 'blue') + 
  facet_grid(condToken ~ .) +
  xlab("Time course") + ylab("Items") + 
  theme_bw(base_size = 12, base_family = "Times") -> onsetP


# Plot adjusted durations 
df_lex %>%
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

lexBinPlots <- plot_grid(onsetP, onsetAdjP, ncol = 2)

# Bin adjustments look good, save plots as .png file 
# ggsave('lexBinPlots.png', plot = lexBinPlots, dpi = 300, device = "png", path = "./mySources/figs/lexical/general")









# Time sequence for plots 

# Create vector of times for canta 
colTimes <- c(100, 296, 1178, 1890, 2030, 2333, 2457, 2767, 2946, 3471)
# Create vector of labels 
colSeq   <- c("el", "cocinero", "escogio", "c", "o", "l", " para", "el", " menu", "end")
# add labels as names arg for vector colTimes
names(colTimes) <- colSeq
# adjust bins 
colAdj <- (colTimes / 10) - (colTimes[6] / 10)
# turn in into a dataframe 
col_df_temp <- data.frame(group = rep(c("int", "la", "ss"), each = 10), 
                          step  = 1:10, 
                          condition = "monosyllabic", 
                          binAdj     = colAdj, 
                          binN     = colAdj, 
                          targetProp = 0.05, 
                          text  = names(colAdj))


colesTimes <- c(100, 309, 1270, 1987, 2117, 2347, 2467, 2559, 2639, 3579)
# Create vector of labels 
colesSeq   <- c("el", "cocinero", "escogio", "c", "o", "l", "e", "s", "  para", "end")
# add labels as names arg for vector cantaTimes
names(colesTimes) <- colesSeq
# adjust bins 
colesAdj <- (colesTimes / 10) - (colesTimes[6] / 10)
# turn in into a dataframe 
coles_df_temp <- data.frame(group = rep(c("int", "la", "ss"), each = 10), 
                            step  = 1:10, 
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


# Relevel factor
lex_short$condition <- factor(lex_short$condition, levels = c("monosyllabic", "bisyllabic"))
lex_short$group <- factor(lex_short$group, levels = c("la", "int", "ss"))


lex_short %>%
  na.omit(.) %>% 
  group_by(., participant, binAdj, group, target, verb) %>%
  summarize(., targetProp = mean(targetProp)) %>%
  ggplot(., aes(x = binAdj, y = targetProp)) + 
  facet_grid(group ~ condition, labeller = as_labeller(condition_names)) + 
  geom_vline(xintercept = 0, color = 'grey60') + 
  geom_hline(yintercept = 0.5, color = 'grey60') + 
  stat_summary(fun.data = mean_cl_boot, geom = 'errorbar', width = 0, size = 0.1,
               show.legend = FALSE, color = 'darkgrey') +
  stat_summary(fun.y = mean, geom = 'point', size = 0.2, color = 'darkgreen') + 
  stat_summary(data = na.omit(lex_short), aes(x = binAdj, y = distractorProp), fun.data = mean_cl_boot, geom = 'errorbar', width = 0, size = 0.1, color = 'darkgrey') + 
  stat_summary(data = na.omit(lex_short), aes(x = binAdj, y = distractorProp), fun.y = mean, geom = 'point', size = 0.2, color = 'red') + 
  ylab('Proportion of fixations') + 
  xlab('Adjusted time course') + 
  ylim(0, 1) + 
  #geom_text(data = colesEx[colesEx$step >= 3 & colesEx$step <= 9, ], aes(label = text), color = 'blue', hjust = "left") + 
  theme_bw(base_size = 16, base_family = "Times") 



















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
  group_by(., group, condition) %>%
  summarise(., meanFix = mean(targetProp))

# We will test this for each group in each condition (mono, di)
# using a one-sided t-test. Specifically, we are testing the 
# hypothesis that the proportion of looks is greater than 
# chance (50%). 
# - H0: u = 0.5
# - Ha: u > 0.5
# The generic code is: t.test(myVector, alternative = "greater", my = 0.50, conf.level = 0.95)

dur_ttest <- dur_0 %>%
  na.omit() %>%
  group_by(., group, condition, participant) %>%
  summarise(., meanFix = mean(targetProp)) %>% 
  do(tidy(t.test(.$meanFix, alternative = "greater", mu = 0.50, conf.level = 0.95)))

# Convert pvalues from scientific notation
dur_ttest$p.value <- format(dur_ttest$p.value, scientific = F)
dur_ttest$sig <- "N.S."
dur_ttest[dur_ttest$p.value <= 0.05, 'sig'] <- "*"

# Print results
print(as.data.frame(dur_ttest[, c(1:7, 11)]))

# group    condition  estimate statistic          p.value parameter  conf.low
#    la monosyllabic 0.7104167  4.675311 0.00003124772865        29 0.6339459
#    la   bisyllabic 0.7552778  7.626087 0.00000001040995        29 0.6984007
#   int monosyllabic 0.6687500  2.510497 0.01664180953451         9 0.5455322
#   int   bisyllabic 0.8462500  9.556776 0.00000260585647         9 0.7798348
#    ss monosyllabic 0.7656056  7.777062 0.00000004703496        22 0.7069609
#    ss   bisyllabic 0.7233696  6.129489 0.00000180339132        22 0.6607937


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
    theme_bw(base_size = 16, base_family = 'Times') -> durTargetFixMODP3

# Looks good, save as .png file. 
# ggsave('durTargetFixMODP3.png', plot = durTargetFixMODP3, dpi = 300, device = "png", path = "./mySources/figs/dur/general")





dur_0_prop <- dur_0 %>%
  na.omit() %>%
  group_by(., group, condition, participant) %>%
  summarise(., meanFix = mean(targetProp))  

# Add working memory covariate
# NEED DATA FOR INT



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

#        Df     AIC     BIC logLik deviance  Chisq Chi Df Pr(>Chisq)  
# object  3 -50.021 -41.512 28.011  -56.021                           
# ..1     5 -46.269 -32.087 28.134  -56.269 0.2478      2    0.88348  
# ..2     6 -45.329 -28.312 28.665  -57.329 1.0605      1    0.30310  
# ..3     8 -46.544 -23.853 31.272  -62.544 5.2142      2    0.07375 .

summary(prop_0_mod_full)

#                               Estimate Std. Error        df t value Pr(>|t|)
# (Intercept)                    0.71042    0.03460 124.18000  20.535   <2e-16
# groupint                      -0.04167    0.06919 124.18000  -0.602    0.548
# groupss                        0.05519    0.05252 124.18000   1.051    0.295
# conditionbisyllabic            0.04486    0.04587  63.00000   0.978    0.332
# groupint:conditionbisyllabic   0.13264    0.09174  63.00000   1.446    0.153
# groupss:conditionbisyllabic   -0.08710    0.06963  63.00000  -1.251    0.216

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
    theme_bw(base_size = 16, base_family = 'Times') -> stressTargetFixLSRLp2
















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

#    - the lowest target word onset is: 'firma' @ -33.75
#    - thus we can use -35 binAdj as the starting point 
#      and be sure we are including the entire target word
#    - we can also do a little higher to lighten the models (-30)
dur_gc_subset <- filter(dur_short_subset, binAdj >= -35 & binAdj <= 5)


# - Readjust time course 
#    - now we will make the time course positive, starting at 1
#    - to do this we will add the lowest 'binAdj' value to each 
#      bin, plus 1 (to avoid starting at 0)

dur_gc_subset$binGC <- dur_gc_subset$binAdj + 36


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


# Base model 
if(F){
  gc_mod_base <- lmer(targetProp ~ (ot1+ot2) + 
                 ((ot1+ot2) | participant) + 
                 ((ot1+ot2) | target),
                 control = lmerControl(optimizer = 'bobyqa'), 
                 data = dur_gc_subset, REML = F)
  saveRDS(gc_mod_base, file = "./mySources/models/dur/s3_adv_int_nat/eye_track/gc_mod_base.rds", compress = 'xz')
}


# Add group effect on intercept 
if(F){
  gc_mod_group_0 <- lmer(targetProp ~ (ot1+ot2) + group + 
                    ((ot1+ot2) | participant) + 
                    ((ot1+ot2) | target), 
                    control = lmerControl(optimizer = 'bobyqa'), 
                    data = dur_gc_subset, REML = F)
  saveRDS(gc_mod_group_0, file = "./mySources/models/dur/s3_adv_int_nat/eye_track/gc_mod_group_0.rds", compress = 'xz')
}

# Add group effect on slope
if(F){
 gc_mod_group_1 <- lmer(targetProp ~ (ot1+ot2) + group + 
                   ot1:group + 
                   ((ot1+ot2) | participant) + 
                   ((ot1+ot2) | target),
                   control = lmerControl(optimizer = 'bobyqa'), 
                   data = dur_gc_subset, REML = F)
  saveRDS(gc_mod_group_1, file = "./mySources/models/dur/s3_adv_int_nat/eye_track/gc_mod_group_1.rds", compress = 'xz')
}

# Add group effect on quadratic poly 
if(F){
  gc_mod_group_2 <- lmer(targetProp ~ (ot1+ot2) + group + 
                    ot1:group + ot2:group + 
                    ((ot1+ot2) | participant) + 
                    ((ot1+ot2) | target),
                    control = lmerControl(optimizer = 'bobyqa'), 
                    data = dur_gc_subset, REML = F)
  saveRDS(gc_mod_group_2, file = "./mySources/models/dur/s3_adv_int_nat/eye_track/gc_mod_group_2.rds", compress = 'xz')
}


# Add condition effect on intercept 
if(F){
  gc_mod_cond_0 <- lmer(targetProp ~ (ot1+ot2) * group + condition + 
                   ((ot1+ot2) | participant) + 
                   ((ot1+ot2) | target),
                   control = lmerControl(optimizer = 'bobyqa'), 
                   data = dur_gc_subset, REML = F)
  saveRDS(gc_mod_cond_0, file = "./mySources/models/dur/s3_adv_int_nat/eye_track/gc_mod_cond_0.rds", compress = 'xz')
}

# Add condition effect on slope 
if(F){
  gc_mod_cond_1 <- lmer(targetProp ~ (ot1+ot2) * group + condition + 
                   ot1:condition + 
                   ((ot1+ot2) | participant) + 
                   ((ot1+ot2) | target),
                   control = lmerControl(optimizer = 'bobyqa'), 
                   data = dur_gc_subset, REML = F)
  saveRDS(gc_mod_cond_1, file = "./mySources/models/dur/s3_adv_int_nat/eye_track/gc_mod_cond_1.rds", compress = 'xz')
}

# Add condition effect on quadratic poly 
if(F){
  gc_mod_cond_2 <- lmer(targetProp ~ (ot1+ot2) * group + condition + 
                   ot1:condition + ot2:condition + 
                   ((ot1+ot2) | participant) + 
                   ((ot1+ot2) | target),
                   control = lmerControl(optimizer = 'bobyqa'), 
                   data = dur_gc_subset, REML = F)
  saveRDS(gc_mod_cond_2, file = "./mySources/models/dur/s3_adv_int_nat/eye_track/gc_mod_cond_2.rds", compress = 'xz')
}


# Include all interactions
if(F){
gc_mod_full <- lmer(targetProp ~ (ot1+ot2) * group * condition + 
               ((ot1+ot2) | participant) + 
               ((ot1+ot2) | target),
               control = lmerControl(optimizer = 'bobyqa'), 
               data = dur_gc_subset, REML = F)
  saveRDS(gc_mod_full, file = "./mySources/models/dur/s3_adv_int_nat/eye_track/gc_mod_full.rds", compress = 'xz')
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

#        Df   AIC   BIC logLik deviance    Chisq Chi Df Pr(>Chisq)    
# object 16 57804 57944 -28886    57772                                # base model 
# ..1    18 57802 57960 -28883    57766   5.7529      2   0.056335 .   # add group effect on intercept 
# ..2    20 57804 57979 -28882    57764   2.6237      2   0.269326     # add group effect on slope
# ..3    22 57806 58000 -28881    57762   1.1935      2   0.550610     # add group effect on quadratic poly 
# ..4    23 57808 58010 -28881    57762   0.5912      1   0.441941     # add cond effect on intercept
# ..5    24 57801 58012 -28877    57753   8.6536      1   0.003264 **  # add cond effect on slope 
# ..6    25 57797 58017 -28874    57747   5.9869      1   0.014412 *   # add cond effect on quadratic poly
# ..7    31 57647 57920 -28793    57585 161.7451      6  < 2.2e-16 *** # full model 



# summary(gc_mod_full)

#Fixed effects:
#                                   Estimate Std. Error         df t value  Pr(>|t|)    
#(Intercept)                       5.862e-01  3.437e-02  3.300e+01  17.055   < 2e-16 ***
#ot1                               8.542e-01  1.101e-01  5.300e+01   7.762  2.62e-10 ***
#ot2                               4.745e-02  5.393e-02  8.300e+01   0.880  0.381585    
#groupla                           1.096e-02  2.621e-02  6.500e+01   0.418  0.677333    
#groupint                          5.275e-03  3.591e-02  6.600e+01   0.147  0.883661    
#conditionbisyllabic               1.095e-01  4.040e-02  1.700e+01   2.711  0.015040 *  
#ot1:groupla                      -2.877e-01  1.074e-01  7.200e+01  -2.678  0.009178 ** 
#ot1:groupint                     -5.303e-01  1.477e-01  7.400e+01  -3.589  0.000593 ***
#ot2:groupla                       6.254e-02  6.281e-02  1.000e+02   0.996  0.321815    
#ot2:groupint                      1.179e-01  8.758e-02  1.100e+02   1.346  0.181004    
#ot1:conditionbisyllabic          -6.864e-01  1.143e-01  2.000e+01  -6.003  7.27e-06 ***
#ot2:conditionbisyllabic          -9.336e-02  5.698e-02  5.000e+01  -1.639  0.107591    
#groupla:conditionbisyllabic      -2.226e-03  8.931e-03  4.790e+04  -0.249  0.803149    
#groupint:conditionbisyllabic      1.120e-01  1.294e-02  4.790e+04   8.653   < 2e-16 ***
#ot1:groupla:conditionbisyllabic   3.242e-01  5.719e-02  4.790e+04   5.668  1.45e-08 ***
#ot1:groupint:conditionbisyllabic  6.195e-01  8.289e-02  4.790e+04   7.474  7.90e-14 ***
#ot2:groupla:conditionbisyllabic  -2.783e-03  5.718e-02  4.792e+04  -0.049  0.961186    
#ot2:groupint:conditionbisyllabic -1.477e-01  8.289e-02  4.790e+04  -1.781  0.074848 .  


# create new df including the fitted model 
data.comp <- data.frame(na.omit(dur_gc_subset), 
                        GCA_Full = fitted(gc_mod_full))
# glimpse(data.comp)


# I dont remember how I determined this
suffix_area <- data.frame(x = 36:41, y = Inf)

condition_namesModLSRL <- c(
                    `ss` = "SS", 
                    `la` = "LA", 
                    `int` = "IN" 
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
  scale_x_continuous(breaks = c(1, 36), labels = c("Approx.\ntarget\nonset", "Target\nsyllable\noffset")) + 
  scale_color_brewer(palette = "Set1", name = "", labels = c("Monosyllabic", "Bisyllabic")) + 
  theme_bw(base_size = 16, base_family = "Times New Roman") -> durGCAfullMod

# ggsave('stressGCAfullMod.png', plot = stressGCAfullMod, dpi = 600, device = "png", path = "./mySources/figs/stress/s1_beg_adv_nat/eye_track")

