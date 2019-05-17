#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Morphosyntactic predictability                                              #
# Duration analyses                                                           #
# 03/30/2017                                                                  #
# Script 2                                                                    #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# clean working directory
rm(list = ls(all = TRUE))

# Set working directory
setwd("~/Desktop/morph_pred/")

library(dplyr); library(tidyr); library(ggplot2); library(plotly)
library(lme4); library(lmerTest); library(gridExtra); library(cowplot)
library(broom)

# read data
df_dur <- read.csv("./mySources/data/clean/durationBIN5Clean.csv", header = TRUE, quote = "", sep = ',')

# Convert condToken variable to factor
df_dur$condToken <- as.factor(df_dur$condToken)

# Glimpse of data structure
glimpse(df_dur)






# Token realignment
# - 1st we calculate the lowest max 
# - 2nd we calculate the highest min 
# These values will become the extremes for our time course 

# calculate lowest max
maxes <- df_dur %>%
  group_by(., condLabel, verb) %>%
  summarize(max = max(binAdj)) %>% 
  as.data.frame(.)
binAdjMaxMin <- min(maxes$max)

# calculate highest low
mins <- df_dur %>%
  group_by(., condLabel, verb) %>%
  summarize(min = min(binAdj)) %>% 
  as.data.frame(.)
binAdjMinMax <- max(mins$min)

# subset data based on new ranges 
dur_short <- df_dur %>% filter(., binAdj <= binAdjMaxMin & binAdj >= binAdjMinMax)

# create new adjusted variable that ranges from 1 to max
dur_short$binREadj <- (dur_short$binAdj - binAdjMinMax) + 1


# add polynomials for growth curve analyses
t <- poly(min(dur_short$binREadj):max(dur_short$binREadj), 3)
dur_short[, paste('ot', 1:3, sep = "")] <- t[dur_short$binREadj, 1:3]

# check it 
glimpse(dur_short)

# Create subsets for individual analyses
hs <- dur_short %>% dplyr::filter(., group == 'hs')
ss <- dur_short %>% dplyr::filter(., group == 'ss')
lb <- dur_short %>% dplyr::filter(., group == 'lb')
la <- dur_short %>% dplyr::filter(., group == 'la')
int <- dur_short %>% dplyr::filter(., group == 'int')


# Check em 
glimpse(hs); glimpse(ss); glimpse(lb); glimpse(la); glimpse(int)








## @knitr binAdjustments

# Bin adjustment plots
# - First we calculate the target onset for each item 
# - Next we calculate the adjusted onset 

onsets <- df_dur %>%
  group_by(., condLabel, verb, condToken) %>%
  summarize(., onset = (unique(word5) - unique(word4_c1v1)) / 5)

onsetAdj <- df_dur %>%
  group_by(., condLabel, verb, condToken) %>%
  summarize(., onset = ((unique(word5) - unique(word4_c1v1)) / 5) - 
                        (unique(word5) - unique(word4_c1v1)) / 5)

# Plot unadjusted durations 
ggplot(onsets, aes(x = condLabel, y = onset, 
                   group = interaction(condToken, condLabel), 
                   color = condToken, dodge = condToken)) + 
  geom_line(data = df_dur, aes(x = condLabel, y = bin), 
            position = position_dodge(width = 1), color = 'black', 
            show.legend = FALSE, size = 0.5) + 
  geom_point(position = position_dodge(width = 1), show.legend = FALSE) + 
  scale_color_brewer(palette = "Set1") + 
  ylab("Time course") + xlab("Items") + 
  coord_flip() + 
  theme_bw(base_size = 12, base_family = "Times") -> onsetP


# Plot adjusted durations 
ggplot(onsetAdj, aes(x = condLabel, y = onset, 
                   group = interaction(condToken, condLabel), 
                   color = condToken, dodge = condToken)) + 
  geom_line(data = df_dur, aes(x = condLabel, y = binAdj), 
            position = position_dodge(width = 1), color = 'black', 
            show.legend = FALSE, size = 0.5) + 
  geom_point(position = position_dodge(width = 1), show.legend = FALSE) + 
  geom_hline(yintercept = binAdjMinMax, color = 'red') + 
  geom_hline(yintercept = binAdjMaxMin, color = 'red') + 
  scale_color_brewer(palette = "Set1") + 
  ylab("Adjusted time course") + xlab("Items") + 
  scale_x_discrete(position = "top") +
  coord_flip() + 
  theme_bw(base_size = 12, base_family = "Times") -> onsetAdjP

durBinPlots <- plot_grid(onsetP, onsetAdjP, ncol = 2)

# Bin adjustments look good, save plots as .png file 
# ggsave('durBinPlots.png', plot = durBinPlots, dpi = 300, device = "png", path = "./mySources/figs/dur/general")














## @knitr durPlotAll

# We make a plot of all the data 
# Does not include models 
# This is just to get an idea of what is going on 

# Relevel factor
dur_short$group <- factor(dur_short$group, levels = c("lb", "la", "hs", "int", "ss"))
dur_short$condition <- factor(dur_short$condition, levels = c("monosyllabic", "disyllabic"))


condition_names <- c(`monosyllabic` = "Monosyllabic", `disyllabic` = "Disyllabic", `lb` = "LB", `la` = "LA", `hs` = "HS", `int` = "IN", `ss` = "SS")

dur_short %>% 
  ggplot(., aes(x = binAdj, y = propT)) + 
  facet_grid(group ~ condition, labeller = as_labeller(condition_names)) + 
  geom_vline(xintercept = 0, color = 'grey60') + 
  geom_hline(yintercept = 0.33, color = 'grey60') + 
  ylim(0, 1) +
  stat_summary(fun.data = mean_se, geom = 'errorbar', width = 0.2, size = 0.1,
               show.legend = FALSE) +
  stat_summary(fun.y = mean, geom = 'point', size = 0.2, color = 'darkgreen') + 
  stat_summary(data = dur_short, aes(x = binAdj, y = propD), fun.data = mean_se, geom = 'errorbar', width = 0.2, size = 0.1) + 
  stat_summary(data = dur_short, aes(x = binAdj, y = propD), fun.y = mean, geom = 'point', size = 0.2, color = 'red') + 
  stat_summary(data = dur_short, aes(x = binAdj, y = propE), fun.data = mean_cl_boot, geom = 'errorbar', width = 0.2, size = 0.1) + 
  stat_summary(data = dur_short, aes(x = binAdj, y = propE), fun.y = mean, geom = 'point', size = 0.2, color = 'blue') + 
  ylab('Proportion of fixations') + 
  xlab('Adjusted time course') + 
  theme_bw(base_size = 12, base_family = "Times") -> p1
  # ggplotly(p1)

# We can run 'ggplotly' for interactive plot 
# Looks good, we save it as a .png file 
# ggsave('durP1.png', plot = p1, dpi = 300, device = "png", path = "./mySources/figs/dur/general")









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

dur_0 <- dur_short %>% filter(., binAdj == 0)

# Quick and dirty mean of target fixations as a function of 
# group and condition (mono, di)

dur_0 %>%
  group_by(., group, condition) %>%
  summarise(., meanFix = mean(propT))

# Calculate mean target fixation as a function of group, condition, 
# for each participant. We will plot the mean and calculate the 
# bootstrapped 95% confidence interval and plot it all. 

dur_0 %>%
  group_by(., group, condition, participant) %>%
  summarise(., meanFix = mean(propT)) %>%
  ggplot(., aes(x = group, y = meanFix, 
                dodge = condition, color = condition,
                group = interaction(group, condition))) +
    geom_hline(yintercept = 0.33, color = "black", size = 0.75, 
               lty = 3) + 
    stat_summary(fun.data = mean_cl_boot, geom = 'errorbar', 
                 position = position_dodge(width = 0.5), 
                 width = 0.35, color = 'grey40') + 
    stat_summary(fun.y = mean, geom = 'point', size = 4,
                 position = position_dodge(width = 0.5)) + 
    ylim(0, 1) + ylab('Target fixations') + xlab('Group') + 
    ggtitle('Mean target fixations as a function of group and target type 20 ms after vowel onset') + 
    scale_color_brewer(palette = "Set1", name = '') + 
    theme_bw(base_size = 12, base_family = 'Times') -> durTargetFixP2

# Looks good, save as .png file
# ggsave('durTargetFixPlot.png', plot = durTargetFixP2, dpi = 300, device = "png", path = "./mySources/figs/dur/general")


# Impressionistically, it looks like maybe the LA
# group predicts above chance for the monosyllabic 
# items, and all groups predict above chance for the 
# disyllabic items

# We will test this for each group in each condition (mono, di)
# using a one-sided t-test. Specifically, we are testing the 
# hypothesis that the proportion of looks is greater than 
# chance (33%). 
# - H0: u = 0.33
# - Ha: u > 0.33
# The generic code is: t.test(myVector, alternative = "greater", my = 0.33, conf.level = 0.95)

dur_ttest <- dur_0 %>%
  group_by(., group, condition, participant) %>%
  summarise(., meanFix = mean(propT)) %>% 
  do(tidy(t.test(.$meanFix, alternative = "greater", mu = 0.33, conf.level = 0.95)))

# Convert pvalues from scientific notation
dur_ttest$p.value <- format(dur_ttest$p.value, scientific = F)
dur_ttest$sig <- "N.S."
dur_ttest[dur_ttest$p.value <= 0.05, 'sig'] <- "*"

# Print results
print(as.data.frame(dur_ttest[, c(1:7, 11)]))

# group   condition      mean         t.value     p.value          df   conf.low     sig
# lb      monosyllabic   0.2799539   -1.3461869   0.905836531644   30   0.2168562   N.S.
# lb      disyllabic     0.4606375    2.7949691   0.004481330107   30   0.3813071      *
# la      monosyllabic   0.3933333    1.6941475   0.050476253655   29   0.3298138   N.S.
# la      disyllabic     0.5443056    5.7666437   0.000001513397   29   0.4811610      *
# hs      monosyllabic   0.4125000    1.4830675   0.074421863856   29   0.3179811   N.S.
# hs      disyllabic     0.4650000    2.6979850   0.005752664612   29   0.3799802      *
# int     monosyllabic   0.3541667    0.4078114   0.345618487952   11   0.2477436   N.S.
# int     disyllabic     0.5833333    2.9894702   0.006154675010   11   0.4311467      *
# ss      monosyllabic   0.3628571    0.9469285   0.176556355937   24   0.3034918   N.S.
# ss      disyllabic     0.5650000    4.6941360   0.000045191428   24   0.4793490      *

# For all groups:
# - significantly above chance for disyllabic targets
# - at or below chance for monosyllabic targets

# We will plot the models 
# This will be almost exactly the same as 
# the previous plot, but it will use the 
# confidence interval from the test we 
# actually conducted
dur_ttest$group <- factor(dur_ttest$group, levels = c("lb", "la", "hs", "int", "ss"))
dur_ttest$condition <- factor(dur_ttest$condition, levels = c("monosyllabic", "disyllabic"))


ggplot(dur_ttest, aes(x = group, y = estimate, color = condition, 
                      group = interaction(group, condition), dodge = condition)) +
    geom_hline(yintercept = 0.33, lty = 3) + 
    geom_linerange(aes(ymin = conf.low, ymax = estimate), color = 'grey40',
                     position = position_dodge(width = 0.75), size = 1) +
    geom_point(position = position_dodge(width = 0.75), size = 4) +
    ylim(0, 1.0) + ylab('Target fixations') + xlab('Group') + 
    scale_color_brewer(palette = "Set1", name = "") + 
    ggtitle('Mean target fixations and lower-bound 95% confidence interval 20 ms after vowel onset') + 
    theme_bw(base_size = 12, base_family = 'Times') -> durTargetFixMODP3

# Looks good, save as .png file. 
# ggsave('durTargetFixMODP3.png', plot = durTargetFixMODP3, dpi = 300, device = "png", path = "./mySources/figs/dur/general")


















################################
# Growth curve analyses        #
# - Compare LA, INT, SS for QP #
# - Looks at target fixation   #
#   slopes                     #
################################





########
# mono #
########

# - Make subset of LA, SS, INT subjects
#   for monosyllabic condition 
# - We limit the area of interest to make 
#   the models faster (start of vowel to 
#   start of following word)

# word4_c1 = vowel onset 
# word5    = start of following word 


word5_onsets <- dur_short %>%
  group_by(., condLabel, verb) %>%
  summarize(max = max(word5)) %>% 
  as.data.frame(.)
word5_onset_max <- max(word5_onsets$max)





mono_dur_qp <- filter(dur_short, group != 'lb' & 
                                 group != 'hs' & 
                                 condition == 'monosyllabic' & 
                                 binAdj >= -5 &
                                 binAdj <= 75)
str(mono_dur_qp)


mono_mod_base  <- lmer(propT ~ (ot1+ot2) + 
                      ((ot1+ot2) | participant) + 
                      ((ot1+ot2) | condLabel), 
                      control = lmerControl(optimizer = 'bobyqa'), 
                      data = mono_dur_qp, REML = F)

mono_mod_group <- lmer(propT ~ (ot1+ot2) + group + 
                      ((ot1+ot2) | participant) + 
                      ((ot1+ot2) | condLabel), 
                      control = lmerControl(optimizer = 'bobyqa'), 
                      data = mono_dur_qp, REML = F)

mono_mod_int <- lmer(propT ~ (ot1+ot2) * group + 
                      ((ot1+ot2) | participant) + 
                      ((ot1+ot2) | condLabel), 
                      control = lmerControl(optimizer = 'bobyqa'), 
                      data = mono_dur_qp, REML = F)

# Is there a main effect of group?
# Is there a difference in the effect
# of the higher order polynomial on the 
# slope of the groups?
anova(mono_mod_base, mono_mod_group, mono_mod_int, test = 'Chisq')
summary(mono_mod_base)

mono_dur_qp %>% 
  ggplot(., aes(x = binAdj, y = propT, color = group)) + 
  geom_vline(xintercept = 0, color = 'grey50') + 
  # stat_summary(fun.y = mean, geom = 'line', size = 0.4) + 
  # stat_summary(fun.y = mean, geom = 'point', size = 0.4) + 
  # stat_summary(fun.data = mean_se, geom = 'pointrange', size = 0.4) + 
  stat_summary(aes(y = fitted(mono_mod_base)), fun.y = mean, 
                   geom = 'line', size = 0.4) + 
  xlab("") +
  ylab("") +
  coord_cartesian(ylim = c(0, 1)) + 
  theme_bw(base_size = 12, base_family = "Times New Roman")








# disyllabic

di_dur_qp <- filter(dur_short, group != 'lb' & 
                                 group != 'hs' & 
                                 condition == 'disyllabic' & 
                                 binAdj >= -5)
str(di_dur_qp)


di_mod_base  <- lmer(propT ~ (ot1+ot2) + 
                      ((ot1+ot2) | participant) + 
                      ((ot1+ot2) | condLabel), 
                      control = lmerControl(optimizer = 'bobyqa'), 
                      data = di_dur_qp, REML = F)

di_mod_group <- lmer(propT ~ (ot1+ot2) + group + 
                      ((ot1+ot2) | participant) + 
                      ((ot1+ot2) | condLabel), 
                      control = lmerControl(optimizer = 'bobyqa'), 
                      data = di_dur_qp, REML = F)

di_mod_int <- lmer(propT ~ (ot1+ot2) * group + 
                      ((ot1+ot2) | participant) + 
                      ((ot1+ot2) | condLabel), 
                      control = lmerControl(optimizer = 'bobyqa'), 
                      data = di_dur_qp, REML = F)

# Is there a main effect of group?
# Is there a difference in the effect
# of the higher order polynomial on the 
# slope of the groups?
anova(di_mod_base, di_mod_group, di_mod_int, test = 'Chisq')
summary(di_mod_group)

di_dur_qp %>% 
  ggplot(., aes(x = binAdj, y = propT, color = group)) + 
  geom_vline(xintercept = 0, color = 'grey50') + 
  # stat_summary(fun.y = mean, geom = 'line', size = 0.4) + 
  # stat_summary(fun.y = mean, geom = 'point', size = 0.4) + 
  # stat_summary(fun.data = mean_se, geom = 'pointrange', size = 0.4) + 
  stat_summary(aes(y = fitted(di_mod_int)), fun.y = mean, 
                   geom = 'line', size = 0.4) + 
  xlab("") +
  ylab("") +
  coord_cartesian(ylim = c(0, 1)) + 
  theme_bw(base_size = 12, base_family = "Times New Roman")
