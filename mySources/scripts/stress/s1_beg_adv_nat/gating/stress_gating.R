# clean working directory
rm(list = ls(all = TRUE))

## @knitr stressLibs

library(plotly); library(tidyverse); library(broom); library(sjPlot)
library(lme4); library(lmerTest); library(gridExtra); library(cowplot)
library(foreign)

## @knitr ignore

# Set working directory
# setwd("~/Desktop/morph_pred/")
setwd("~/academia/research/in_progress/morph_pred")


dataset <- read.spss('./mySources/data/gating.sav', to.data.frame=TRUE)

str(dataset)

df <- dataset %>%
      select(Exp, Group, participant = ExperimentName, Condition, Target, Accuracy) %>% 
      filter(., Exp == 'S', Group != 'HS' & Group != 'IN') %>%
      as.tbl(.)

# remove unwanted characters from column
df$Target <- gsub(" ", "", paste(df$Target))
df$Target <- gsub("\303\263", "o", paste(df$Target))
df$participant <- gsub(" ", "", paste(df$participant))

str(df)
remove <- c("L20", "L21", "L22", "L23", "L30", "L31", "L02", "L05", "L06", "L08", "L10", "L15")
df <- filter(df, !participant %in% remove)



































# Quick and dirty mean of target fixations as a function of 
# group and condition (stressed, unstressed *1st syllable*)
df %>% 
  na.omit(.) %>% 
  filter(., Group == 'L') %>%
  group_by(., Group, participant, Condition) %>%
  summarise(., meanAccuracy = mean(Accuracy))


# We will test this for each group in each condition (stressed, untressed)
# using a one-sided t-test. Specifically, we are testing the 
# hypothesis that the proportion of looks is greater than 
# chance (50%). 
# - H0: u = 0.50
# - Ha: u > 0.50
# The generic code is: t.test(myVector, alternative = "greater", my = 0.33, conf.level = 0.95)

stress_gating_ttest <- df %>%
  na.omit(.) %>% 
  group_by(., Group, Condition, participant) %>%
  summarise(., meanAccuracy = mean(Accuracy)) %>% 
  do(tidy(t.test(.$meanAccuracy, alternative = "greater", mu = 0.5, conf.level = 0.95)))

# Convert pvalues from scientific notation
stress_gating_ttest$p.value <- format(stress_gating_ttest$p.value, scientific = F)
stress_gating_ttest$sig <- "N.S."
stress_gating_ttest[stress_gating_ttest$p.value <= 0.05, 'sig'] <- "*"

# Print results
# print(as.data.frame(stress_gating_ttest[, c(1:7, 11)]))

# Group Condition  estimate statistic              p.value parameter  conf.low   sig
#     S         1 0.7915789  9.978069 0.000000000256540094        24 0.7415836     *
#     S         2 0.7915789  9.212883 0.000000001188597964        24 0.7374312     *
#     L         1 0.7000000  3.289526 0.004692720858268470         9 0.5885485     *
#     L         2 0.5789474  1.310001 0.111321712526483943         9 0.4684746  N.S.
#    LA         1 0.8036437  9.067790 0.000000001111816307        25 0.7464450     *
#    LA         2 0.8198381 11.666952 0.000000000006550283        25 0.7730111     *



# Only LB do not predict above chance in oxytone condition

# We will plot the models 
# This will be almost exactly the same as 
# the previous plot, but it will use the 
# confidence interval from the test we 
# actually conducted
stress_gating_ttest$Group <- factor(stress_gating_ttest$Group, levels = c("L", "LA", "S"))
stress_gating_ttest$Condition <- factor(stress_gating_ttest$Condition, levels = c("1", "2"))


ggplot(stress_gating_ttest, aes(x = Group, y = estimate, color = Condition, 
                      group = interaction(Group, Condition), dodge = Condition)) +
    geom_hline(yintercept = 0.5, lty = 3) + 
    geom_linerange(aes(ymin = conf.low, ymax = estimate), color = 'grey40',
                     position = position_dodge(width = 0.75), size = 1) +
    geom_point(position = position_dodge(width = 0.75), size = 4) +
    ylim(0, 1.0) + ylab('% correct') + xlab('Group') + 
    ggtitle('Mean accuracy and lower-bound 95% confidence interval') + 
    scale_color_brewer(palette = "Set1", name = '', labels = c('Paroxytone', 'Oxytone')) + 
    scale_x_discrete(labels = c('LB', 'LA', 'SS')) + 
    theme_bw(base_size = 16, base_family = 'Times') -> stressGatingTT

# Looks good, save as .png file. 
# ggsave('stressGatingTT.png', plot = stressGatingTT, dpi = 600, device = "png", path = "./mySources/figs/stress/s1_beg_adv_nat/gating")




## @knitr stressGroupCompare






stress_gating_mod <- df %>%
  na.omit(.) %>% 
  group_by(., Group, Condition, participant) %>%
  summarise(., meanAccuracy = mean(Accuracy)) 


gating_mod_0 <- lmer(meanAccuracy ~ 1 + 
                    (1 | participant) , 
                    data = stress_gating_mod, REML = F, control = lmerControl(optimizer = 'bobyqa'))

gating_mod_group <- lmer(meanAccuracy ~ 1 + Group + 
                        (1 | participant), 
                        data = stress_gating_mod, REML = F, control = lmerControl(optimizer = 'bobyqa'))

gating_mod_cond <- lmer(meanAccuracy ~ 1 + Group + Condition + 
                        (1 | participant), 
                        data = stress_gating_mod, REML = F, control = lmerControl(optimizer = 'bobyqa'))

gating_mod_full <- lmer(meanAccuracy ~ 1 + Group * Condition + 
                      (1 | participant), 
                       data = stress_gating_mod, REML = F, control = lmerControl(optimizer = 'bobyqa'))


# anova(gating_mod_0, gating_mod_group, gating_mod_cond, gating_mod_full, test = "Chisq")

#        Df     AIC     BIC logLik deviance   Chisq Chi Df Pr(>Chisq)    
# object  3 -80.399 -71.987 43.200  -86.399                              
# ..1     5 -93.147 -79.127 51.573 -103.147 16.7475      2  0.0002308 ***
# ..2     6 -91.350 -74.526 51.675 -103.350  0.2034      1  0.6519970    
# ..3     8 -90.236 -67.804 53.118 -106.236  2.8856      2  0.2362683     

# summary(gating_mod_group)

# Fixed effects:
#              Estimate Std. Error        df t value Pr(>|t|)    
# (Intercept)   0.79158    0.02242 122.00000  35.303  < 2e-16 ***
# GroupL       -0.15211    0.04195 122.00000  -3.626 0.000421 ***
# GroupLA       0.02016    0.03140 122.00000   0.642 0.522066    

# Relevel to test lb vs la 
stress_gating_mod$Group <- factor(stress_gating_mod$Group, levels = c("L", "LA",  "S"))

summary(lmer(meanAccuracy ~ 1 + Group + (1 | participant), data = stress_gating_mod, REML = F, control = lmerControl(optimizer = 'bobyqa')))
# beginners and LA different from each other




# Calculate mean target fixation as a function of group, condition, 
# for each participant. We will plot the mean and calculate the 
# bootstrapped 95% confidence interval and plot it all. 
df$Group <- factor(df$Group, levels = c("L", "LA",  "S"))

df %>%
  na.omit(.) %>%
  group_by(., Group, Condition, participant) %>%
  summarise(., meanAccuracy = mean(Accuracy)) %>%
  ggplot(., aes(x = Group, y = meanAccuracy, 
                dodge = Condition, color = as.factor(Condition),
                group = interaction(Group, Condition))) +
    geom_hline(yintercept = 0.5, color = "black", size = 0.75, 
               lty = 3) + 
    stat_summary(fun.data = mean_cl_boot, geom = 'errorbar', 
                 position = position_dodge(width = 0.5), 
                 width = 0.35, color = 'grey40') + 
    stat_summary(fun.y = mean, geom = 'point', size = 4,
                 position = position_dodge(width = 0.5)) + 
    coord_cartesian(ylim = c(0, 1)) + ylab('% correct') + xlab('Group') + 
    scale_x_discrete(labels = c('LB', 'LA', 'SS')) + 
    ggtitle('Mean accuracy as a function of group and target type') + 
    scale_color_brewer(palette = "Set1", name = '', labels = c('Paroxytone', 'Oxytone')) + 
    theme_bw(base_size = 16, base_family = 'Times') -> stressGatingLSRLp2

# Looks good, save as .png file
# ggsave('stressGatingLSRLp2.png', plot = stressGatingLSRLp2, dpi = 600, device = "png", path = "./mySources/figs/stress/s1_beg_adv_nat/gating")







