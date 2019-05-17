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
duration_gating_ttest$p.value <- format(duration_gating_ttest$p.value, scientific = F)
duration_gating_ttest$sig <- "N.S."
duration_gating_ttest[duration_gating_ttest$p.value <= 0.05, 'sig'] <- "*"

# Print results
print(as.data.frame(duration_gating_ttest[, c(1:7, 11)]))

#   Group Condition  estimate statistic                p.value parameter  conf.low  sig
# 1     S         1 0.8566667 12.999574 0.00000000000117198414        24 0.8097255    *
# 2     S         2 0.3566667 -3.299567 0.99849245852779855426        24 0.2823459 N.S.
# 3    IN         1 0.9166667 10.856203 0.00000016170520155609        11 0.8477397    *
# 4    IN         2 0.2638889 -3.900067 0.99876185453533516956        11 0.1551655 N.S.
# 5    LA         1 0.8878205 14.617621 0.00000000000004692761        25 0.8425018    *
# 6    LA         2 0.4038462 -2.247081 0.98314656205816963741        25 0.3307539 N.S.


# We will plot the models 
# This will be almost exactly the same as 
# the previous plot, but it will use the 
# confidence interval from the test we 
# actually conducted
duration_gating_ttest$Group <- factor(duration_gating_ttest$Group, levels = c("LA", "IN", "S"))
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
    scale_color_brewer(palette = "Set1", name = '', labels = c('Monosyllabic', 'Bisyllabic')) + 
    scale_x_discrete(labels = c('LA', 'IN', 'SS')) + 
    theme_bw(base_size = 16, base_family = 'Times') -> durationGatingTT

# Looks good, save as .png file. 
# ggsave('durationGatingTT.png', plot = durationGatingTT, dpi = 600, device = "png", path = "./mySources/figs/dur/s3_adv_int_nat/gating", height = 5, width = 9, unit = "in")





## @knitr durGroupCompare






duration_gating_mod <- df_gt_dur %>%
  na.omit(.) #%>% 
  # group_by(., Group, Condition, participant) %>%
  # summarise(., meanAccuracy = mean(Accuracy)) 



gating_mod_0 <- glmer(Accuracy ~ 1 + 
                    (1 | participant), family = "binomial",
                    data = duration_gating_mod, control = glmerControl(optimizer = 'bobyqa'))

gating_mod_group <- glmer(Accuracy ~ 1 + Group + 
                        (1 | participant), family = "binomial",
                        data = duration_gating_mod, control = glmerControl(optimizer = 'bobyqa'))

gating_mod_cond <- glmer(Accuracy ~ 1 + Group + Condition + 
                        (1 | participant), family = "binomial",
                        data = duration_gating_mod, control = glmerControl(optimizer = 'bobyqa'))

gating_mod_full <- glmer(Accuracy ~ 1 + Group * Condition + 
                      (1 | participant), family = "binomial",
                       data = duration_gating_mod, control = glmerControl(optimizer = 'bobyqa'))

gating_mod_wm <- glmer(Accuracy ~ 1 + Group * Condition + WM + 
                      (1 | participant), family = "binomial",
                       data = duration_gating_mod, control = glmerControl(optimizer = 'bobyqa'))

gating_mod_wm_full <- glmer(Accuracy ~ 1 + Group * Condition + WM + WM:Condition +
                      (1 | participant), family = "binomial",
                       data = duration_gating_mod, control = glmerControl(optimizer = 'bobyqa'))


# anova(gating_mod_0, gating_mod_group, gating_mod_cond, gating_mod_full, gating_mod_wm, gating_mod_wm_full, test = "Chisq")

#                    Df    AIC    BIC  logLik deviance    Chisq Chi Df Pr(>Chisq)    
# gating_mod_0        2 2619.7 2631.0 -1307.9   2615.7                               
# gating_mod_group    4 2620.8 2643.2 -1306.4   2612.8   2.9911      2    0.22413    
# gating_mod_cond     5 2036.1 2064.1 -1013.0   2026.1 586.6929      1    < 2e-16 ***
# gating_mod_full     7 2035.6 2074.9 -1010.8   2021.6   4.4120      2    0.11014    
# gating_mod_wm       8 2034.2 2079.0 -1009.1   2018.2   3.4891      1    0.06177 .  
# gating_mod_wm_full  9 2036.2 2086.6 -1009.1   2018.2   0.0000      1    0.99514    


# summary(gating_mod_cond)

# Fixed effects:
#             Estimate Std. Error z value Pr(>|z|)    
# (Intercept)   4.7477     0.2622  18.105   <2e-16 ***
# GroupIN      -0.2521     0.2687  -0.938    0.348    
# GroupLA       0.2057     0.2155   0.955    0.340    
# Condition    -2.6302     0.1261 -20.862   <2e-16 ***




# Calculate mean target fixation as a function of group, condition, 
# for each participant. We will plot the mean and calculate the 
# bootstrapped 95% confidence interval and plot it all. 
df_gt_dur$Group <- factor(df_gt_dur$Group, levels = c("LA", "IN", "S"))

df_gt_dur %>%
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
    scale_x_discrete(labels = c('LA', 'IN', 'SS')) + 
    ggtitle('Mean accuracy as a function of group and target type') + 
    scale_color_brewer(palette = "Set1", name = '', labels = c('Monosyllabic', 'Bisyllabic')) + 
    theme_bw(base_size = 16, base_family = 'Times') -> durationGatingP2

# Looks good, save as .png file
# ggsave('durationGatingP2.png', plot = durationGatingP2, dpi = 600, device = "png", path = "./mySources/figs/dur/s3_adv_int_nat/gating", height = 5, width = 9, unit = "in")






duration_gating_mod %>%
  na.omit(.) %>%
  filter(., WM >= 5) %>%
  group_by(., participant, Condition, WM) %>%
  summarise(., meanAccuracy = mean(Accuracy), 
               sdAccuracy = sd(Accuracy)) %>%
  ggplot(., aes(x = WM, y = meanAccuracy, 
                color = as.factor(Condition))) + 
    geom_jitter() + 
    geom_smooth(method = "lm", se = T) +
    scale_color_brewer(palette = "Set1", name = "", labels = c("Monosyllabic", "Bisyllabic")) + 
    theme_bw()


