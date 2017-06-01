# clean working directory
rm(list = ls(all = TRUE))

## @knitr stressLibs

library(tidyverse); library(broom); library(sjPlot)
library(lme4); library(lmerTest); library(gridExtra); library(cowplot)
library(foreign)

## @knitr ignore

# Set working directory
# setwd("~/Desktop/morph_pred/")
setwd("~/academia/research/in_progress/morph_pred")


dataset <- read.spss('./mySources/data/raw/gating.sav', to.data.frame=TRUE)

str(dataset)

dataset <- dataset %>% 
  separate(., col = ID, into = c("keep", "trash"), sep = 1) %>%
  filter(., keep != "P")

df_gt_dur <- dataset %>%
             select(Exp, Group, participant = ExperimentName, Condition, Target, WM, Accuracy) %>% 
             filter(., Exp == 'D', !(Group %in% c("HS", "L"))) %>%
             as.tbl(.)

# remove unwanted characters from column
df_gt_dur$Target <- gsub(" ", "", paste(df_gt_dur$Target))
df_gt_dur$Target <- gsub("\303\263", "o", paste(df_gt_dur$Target))
df_gt_dur$participant <- gsub(" ", "", paste(df_gt_dur$participant))

str(df_gt_dur)


































# Quick and dirty mean of target fixations as a function of 
# group and condition (mono, bi word)
df_gt_dur %>% 
  na.omit(.) %>%
  filter(., !(Target %in% c("mes", "meses", "sol", "soles"))) %>%
  group_by(., Group, Condition) %>%
  summarise(., meanAccuracy = mean(Accuracy), 
               sdAccuracy = sd(Accuracy))


# We will test this for each group in each condition (stressed, untressed)
# using a one-sided t-test. Specifically, we are testing the 
# hypothesis that the proportion of looks is greater than 
# chance (50%). 
# - H0: u = 0.50
# - Ha: u > 0.50
# The generic code is: t.test(myVector, alternative = "greater", my = 0.33, conf.level = 0.95)

duration_gating_ttest <- df_gt_dur %>%
  na.omit(.) %>% 
  filter(., !(Target %in% c("mes", "meses", "sol", "soles"))) %>%
  group_by(., Group, Condition, participant) %>%
  summarise(., meanAccuracy = mean(Accuracy)) %>% 
  do(tidy(t.test(.$meanAccuracy, alternative = "greater", mu = 0.5, conf.level = 0.95)))

# Convert pvalues from scientific notation
duration_gating_ttest$p.value <- format(stress_gating_ttest$p.value, scientific = F)
duration_gating_ttest$sig <- "N.S."
duration_gating_ttest[stress_gating_ttest$p.value <= 0.05, 'sig'] <- "*"

# Print results
print(as.data.frame(duration_gating_ttest[, c(1:7, 11)]))

#   Group Condition  estimate statistic                p.value parameter  conf.low  sig
# 1     S         1 0.8566667 12.999574 0.00000000000117198414        24 0.8097255    *
# 2     S         2 0.3566667 -3.299567 0.99849245852779855426        24 0.2823459 N.S.
# 3    IN         1 0.9166667 10.856203 0.00000016170520155609        11 0.8477397    *
# 4    IN         2 0.2638889 -3.900067 0.99876185453533516956        11 0.1551655 N.S.
# 5    LA         1 0.8878205 14.617621 0.00000000000004692761        25 0.8425018    *
# 6    LA         2 0.4038462 -2.247081 0.98314656205816963741        25 0.3307539 N.S.

# All predict above chance in all conditiona

# We will plot the models 
# This will be almost exactly the same as 
# the previous plot, but it will use the 
# confidence interval from the test we 
# actually conducted
duration_gating_ttest$Group <- factor(duration_gating_ttest$Group, levels = c("LA", "HS", "S"))
duration_gating_ttest$Condition <- factor(duration_gating_ttest$Condition, levels = c("1", "2"))


ggplot(duration_gating_ttest, aes(x = Group, y = estimate, color = Condition, 
                      group = interaction(Group, Condition), dodge = Condition)) +
    geom_hline(yintercept = 0.5, lty = 3) + 
    geom_linerange(aes(ymin = conf.low, ymax = estimate), color = 'grey40',
                     position = position_dodge(width = 0.75), size = 1) +
    geom_point(position = position_dodge(width = 0.75), size = 4) +
    # ylim(0, 1.0) + 
    ylab('% correct') + xlab('Group') + 
    ggtitle('Mean accuracy and lower-bound 95% confidence interval') + 
    scale_color_brewer(palette = "Set1", name = '', labels = c('Paroxytone', 'Oxytone')) + 
    scale_x_discrete(labels = c('LA', 'HS', 'SS')) + 
    theme_bw(base_size = 16, base_family = 'Times') -> stressGatingTT

# Looks good, save as .png file. 
# ggsave('stressGatingTT.png', plot = stressGatingTT, dpi = 600, device = "png", path = "./mySources/figs/stress/s2_adv_her_nat/gating")




## @knitr stressGroupCompare






stress_gating_mod <- df %>%
  na.omit(.) #%>% 
  # group_by(., Group, Condition, participant) %>%
  # summarise(., meanAccuracy = mean(Accuracy)) 

# Add working memory covariate
wm_df <- structure(list(participant = structure(1:72, .Label = c("HS03", 
"HS07", "HS08", "HS09", "HS12", "HS13", "HS14", "HS15", "HS17", 
"HS18", "HS19", "HS20", "HS21", "HS22", "HS23", "HS24", "HS25", 
"HS26", "HS27", "HS28", "HS29", "HS30", "HS31", "LA01", "LA02", 
"LA03", "LA08", "LA09", "LA10", "LA11", "LA12", "LA13", "LA15", 
"LA16", "LA17", "LA18", "LA20", "LA21", "LA22", "LA23", "LA25", 
"LA26", "LA27", "LA28", "LA29", "LA30", "LA31", "SB01", "SB02", 
"SB03", "SB04", "SB05", "SB06", "SB07", "SB08", "SB09", "SB10", 
"SB11", "SB12", "SB13", "SB14", "SB15", "SB16", "SB17", "SB18", 
"SB19", "SB20", "SB21", "SB22", "SB23", "SB24", "SB25"), class = "factor"), 
    wm = c(12L, 11L, 5L, 10L, 9L, 10L, 11L, 11L, 12L, 0L, 9L, 
    11L, 12L, 9L, 8L, 11L, 10L, 7L, 9L, 9L, 11L, 9L, 7L, 10L, 
    13L, 1L, 8L, 10L, 9L, 9L, 8L, 10L, 8L, 10L, 9L, 11L, 7L, 
    8L, 10L, 8L, 8L, 10L, 9L, 9L, 8L, 10L, 11L, 10L, 9L, 6L, 
    10L, 11L, 6L, 11L, 9L, 10L, 7L, 8L, 8L, 12L, 8L, 10L, 11L, 
    10L, 6L, 12L, 8L, 9L, 14L, 12L, 6L, 10L)), .Names = c("participant", 
"wm"), class = "data.frame", row.names = c(NA, -72L))

# remove unwanted participants
wm <- wm_df[!wm_df$participant %in% remove, ]

# set wm dataframe as a named vector 
wm_vector <- setNames(as.character(wm$wm), wm$participant)


# subset the wm vector using the participant column of df as a 
# way to look up the wm values in the vector and assign the values 
# to a new 'wm' column 
stress_gating_mod$wm <- wm_vector[as.character(stress_gating_mod$participant)]
stress_gating_mod$wm <- as.numeric(stress_gating_mod$wm)



gating_mod_0 <- glmer(Accuracy ~ 1 + 
                    (1 | participant), family = "binomial",
                    data = stress_gating_mod, control = glmerControl(optimizer = 'bobyqa'))

gating_mod_group <- glmer(Accuracy ~ 1 + Group + 
                        (1 | participant), family = "binomial",
                        data = stress_gating_mod, control = glmerControl(optimizer = 'bobyqa'))

gating_mod_cond <- glmer(Accuracy ~ 1 + Group + Condition + 
                        (1 | participant), family = "binomial",
                        data = stress_gating_mod, control = glmerControl(optimizer = 'bobyqa'))

gating_mod_full <- glmer(Accuracy ~ 1 + Group * Condition + 
                      (1 | participant), family = "binomial",
                       data = stress_gating_mod, control = glmerControl(optimizer = 'bobyqa'))

gating_mod_wm <- glmer(Accuracy ~ 1 + Group * Condition + wm + 
                      (1 | participant), family = "binomial",
                       data = stress_gating_mod, control = glmerControl(optimizer = 'bobyqa'))

gating_mod_wm_full <- glmer(Accuracy ~ 1 + Group * Condition + wm + wm:Condition +
                      (1 | participant), family = "binomial",
                       data = stress_gating_mod, control = glmerControl(optimizer = 'bobyqa'))


# anova(gating_mod_0, gating_mod_group, gating_mod_cond, gating_mod_full, gating_mod_wm, gating_mod_wm_full, test = "Chisq")

#                    Df    AIC    BIC  logLik deviance  Chisq Chi Df Pr(>Chisq)   
# gating_mod_0        2 2856.7 2868.6 -1426.4   2852.7                             
# gating_mod_group    4 2854.5 2878.2 -1423.3   2846.5 6.1930      2   0.045207  * 
# gating_mod_cond     5 2849.8 2879.4 -1419.9   2839.8 6.6915      1   0.009687  **
# gating_mod_full     7 2849.3 2890.7 -1417.7   2835.3 4.4972      2   0.105548    
# gating_mod_wm       8 2850.8 2898.2 -1417.4   2834.8 0.4941      1   0.482107    
# gating_mod_wm_full  9 2848.3 2901.5 -1415.2   2830.3 4.5360      1   0.033189  * 


# summary(gating_mod_wm_full)

# Fixed effects:
#                   Estimate Std. Error z value Pr(>|z|)    
# (Intercept)        2.34598    0.68960   3.402 0.000669 ***
# GroupHS           -1.07752    0.38018  -2.834 0.004593 ** 
# GroupLA           -0.28208    0.39641  -0.712 0.476728    
# Condition         -0.78643    0.40270  -1.953 0.050835 .  
# wm                -0.10125    0.06684  -1.515 0.129845    
# GroupHS:Condition  0.51122    0.22651   2.257 0.024015 *  
# GroupLA:Condition  0.28202    0.23642   1.193 0.232916    
# Condition:wm       0.08350    0.03917   2.132 0.033040 *  

# Relevel to test hs vs la 
stress_gating_mod$Group <- factor(stress_gating_mod$Group, levels = c("LA", "HS",  "S"))

summary(glmer(Accuracy ~ 1 + Group * Condition + (1 | participant), family = "binomial", data = stress_gating_mod, control = glmerControl(optimizer = 'bobyqa')))

# Fixed effects:
#             Estimate Std. Error z value Pr(>|z|)    
# (Intercept)  1.15348    0.19205   6.006  1.9e-09 ***
# GroupHS     -0.45648    0.18602  -2.454  0.01413 *  
# GroupS      -0.12561    0.18498  -0.679  0.49712    
# Condition    0.24357    0.09393   2.593  0.00951 ** 

# hs and LA different from each other
# la not different from natives

# In the oxytone condition, increased WM is associated with more accuracy



# Calculate mean target fixation as a function of group, condition, 
# for each participant. We will plot the mean and calculate the 
# bootstrapped 95% confidence interval and plot it all. 
df$Group <- factor(df$Group, levels = c("LA", "HS",  "S"))

df %>%
  na.omit(.) %>%
  ggplot(., aes(x = Group, y = Accuracy, 
                dodge = Condition, color = as.factor(Condition),
                group = interaction(Group, Condition))) +
    geom_hline(yintercept = 0.5, color = "black", size = 0.75, lty = 3) +
    stat_summary(fun.data = 'mean_cl_boot', geom = 'errorbar', 
                 position = position_dodge(width = 0.5), 
                 width = 0.35, color = 'grey40') + 
    stat_summary(fun.y = 'mean', geom = 'point', size = 4,
                 position = position_dodge(width = 0.5)) + 
    # coord_cartesian(ylim = c(0, 1)) + 
    ylab('% correct') + xlab('Group') + 
    scale_x_discrete(labels = c('LA', 'HS', 'SS')) + 
    ggtitle('Mean accuracy as a function of group and target type') + 
    scale_color_brewer(palette = "Set1", name = '', labels = c('Paroxytone', 'Oxytone')) + 
    theme_bw(base_size = 16, base_family = 'Times') -> stressGatingP2

# Looks good, save as .png file
# ggsave('stressGatingP2.png', plot = stressGatingP2, dpi = 600, device = "png", path = "./mySources/figs/stress/s2_adv_her_nat/gating")






stress_gating_mod %>%
  na.omit(.) %>%
  group_by(., participant, Condition, wm) %>%
  summarise(., meanAccuracy = mean(Accuracy), 
               sdAccuracy = sd(Accuracy)) %>%
  ggplot(., aes(x = wm, y = meanAccuracy, 
                color = as.factor(Condition))) + 
    geom_jitter() + 
    geom_smooth(method = "lm", se = F) +
    scale_color_brewer(palette = "Set1") + 
    theme_bw()


