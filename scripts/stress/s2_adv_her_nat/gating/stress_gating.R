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


dataset <- read.spss('./mySources/data/raw/gating.sav', to.data.frame=TRUE)

str(dataset)

df <- dataset %>%
      select(Exp, Group, participant = ExperimentName, WM, Condition, Target, Accuracy) %>% 
      filter(., Exp == 'S', Group != 'IN' & Group != 'L') %>%
      as.tbl(.)

# remove unwanted characters from column
df$Target <- gsub(" ", "", paste(df$Target))
df$Target <- gsub("\303\263", "o", paste(df$Target))
df$participant <- gsub(" ", "", paste(df$participant))

str(df)

# Remove participants according to homogeneity of variance
# test re: working memory
remove <- c("L20", "L21", "L22", "L23", "L30", "L31", "L02", "L05", 
            "L06", "L08", "L10", "L15")
df <- filter(df, !participant %in% remove)

# Remove Targets that are supposed to be in a different exp
remove2 <- c('mus', 'musgo', 'pan', 'panda', 'sur', 'surco')
df <- filter(df, !Target %in% remove2)


scale_this <- function(x) as.vector(scale(x))

df <- df %>% 
  group_by(., Group) %>%
  mutate(., wmScaled = scale_this(WM), 
            coda = ifelse(Target %in% 
              c('bebe', 'bebio', 'llena', 'lleno', 'sube', 'subio', 
                'come', 'comio', 'saca', 'saco', 'lava', 'lavo', 
                'graba', 'grabo'), 
                yes = 0, no = 1),)


str(df)





























# Quick and dirty mean of target fixations as a function of 
# group and condition (stressed, unstressed *1st syllable*)
df %>% 
  na.omit(.) %>% 
  # filter(., Group == 'L') %>%
  group_by(., Group, Condition) %>%
  summarise(., meanAccuracy = mean(Accuracy), sd = sd(Accuracy))


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
  do(tidy(t.test(.$meanAccuracy, alternative = "greater", mu = 0.5, conf.level = 0.99)))

# Convert pvalues from scientific notation
stress_gating_ttest$p.value <- format(stress_gating_ttest$p.value, scientific = F)
stress_gating_ttest$sig <- "N.S."
stress_gating_ttest[stress_gating_ttest$p.value <= 0.05/6, 'sig'] <- "*"

# Print results
print(as.data.frame(stress_gating_ttest[, c(1:7, 11)]))

# Group Condition  estimate statistic          p.value parameter  conf.low sig
#     S         1 0.7600000  7.566555 0.00000004168888        24 0.6743651   *
#     S         2 0.7575000  7.020461 0.00000014683829        24 0.6660913   *
#    HS         1 0.6250000  3.333333 0.00133758284805        25 0.5318085   *
#    HS         2 0.7403846  5.996739 0.00000145456827        25 0.6407669   *
#    LA         1 0.7668269  6.710165 0.00000024742158        25 0.6680077   *
#    LA         2 0.7884615  8.865812 0.00000000171671        25 0.7076051   *

# Only LB do not predict above chance in oxytone condition

# We will plot the models 
# This will be almost exactly the same as 
# the previous plot, but it will use the 
# confidence interval from the test we 
# actually conducted
stress_gating_ttest$Group <- factor(stress_gating_ttest$Group, levels = c("LA", "HS", "S"))
stress_gating_ttest$Condition <- factor(stress_gating_ttest$Condition, levels = c("1", "2"))


ggplot(stress_gating_ttest, aes(x = Group, y = estimate, color = Condition, 
                      group = interaction(Group, Condition), dodge = Condition)) +
    geom_hline(yintercept = 0.5, lty = 3) + 
    geom_linerange(aes(ymin = conf.low, ymax = estimate), color = 'grey40',
                     position = position_dodge(width = 0.75), size = 1) +
    geom_point(position = position_dodge(width = 0.75), size = 4) +
    geom_point(position = position_dodge(width = 0.75), size = 2.75, color = 'grey90') +
    ylim(0.3, 1.0) +
    ggtitle('Mean accuracy and lower-bound 95% confidence interval') + 
    scale_color_brewer(palette = "Set1", name = '', labels = c('Paroxytone', 'Oxytone')) + 
    scale_x_discrete(labels = c('LA', 'HS', 'SS')) + 
    labs(title = 'Response accuracy as a function of stress.', 
         y = 'Target fixations', x = '', caption = 'Mean and lower-bound 99% CI') + 
    theme_minimal(base_size = 16, base_family = 'Times') -> stressGatingTT

# Looks good, save as .png file. 
# ggsave('stressGatingTT.png', plot = stressGatingTT, dpi = 600, device = "png", 
#         path = "./mySources/figs/stress/s2_adv_her_nat/gating",
#         height = 4, width = 7, unit = 'in')






























## @knitr stressGroupCompare






stress_gating_mod <- df %>%
  na.omit(.) %>% 
  filter(., Exp == "S") %>% 
  mutate(., Condition = as.factor(Condition))

stress_gating_mod$Group <- factor(stress_gating_mod$Group, levels = c("S", "LA",  "HS"))

str(stress_gating_mod)

gating_mod_0 <- glmer(Accuracy ~ 1 + 
                     (1 + Condition * coda | participant) + 
                     (1 | Target), 
                     data = stress_gating_mod, 
                     control = glmerControl(optimizer = 'bobyqa'), 
                     family = 'binomial')

gating_mod_group <- glmer(Accuracy ~ 1 + Group + 
                         (1 + Condition * coda | participant) + 
                         (1 | Target), 
                         data = stress_gating_mod, 
                         control = glmerControl(optimizer = 'bobyqa'), 
                         family = 'binomial')

gating_mod_cond <- glmer(Accuracy ~ 1 + Group + Condition + 
                        (1 + Condition * coda | participant) + 
                        (1 | Target), 
                        data = stress_gating_mod, 
                        control = glmerControl(optimizer = 'bobyqa'), 
                        family = 'binomial')

gating_mod_coda <- glmer(Accuracy ~ 1 + Group + Condition + coda + 
                        (1 + Condition * coda | participant) + 
                        (1 | Target), 
                        data = stress_gating_mod, 
                        control = glmerControl(optimizer = 'bobyqa'), 
                        family = 'binomial')

gating_mod_wm  <- glmer(Accuracy ~ 1 + Group + Condition + coda + wmScaled + 
                       (1 + Condition * coda | participant) + 
                       (1 | Target), 
                       data = stress_gating_mod, 
                       control = glmerControl(optimizer = 'bobyqa'), 
                       family = 'binomial')

gating_mod_int1 <- glmer(Accuracy ~ 1 + Group * Condition + coda + wmScaled + 
                        (1 + Condition * coda | participant) + 
                        (1 | Target), 
                        data = stress_gating_mod, 
                        control = glmerControl(optimizer = 'bobyqa'), 
                        family = 'binomial')

gating_mod_int2 <- glmer(Accuracy ~ 1 + Group + Condition * coda + wmScaled + 
                        (1 + Condition * coda | participant) + 
                        (1 | Target), 
                        data = stress_gating_mod, 
                        control = glmerControl(optimizer = 'bobyqa'), 
                        family = 'binomial')

anova(gating_mod_0, gating_mod_group, test = "Chisq") # main effect of group
anova(gating_mod_group, gating_mod_cond, test = "Chisq") # no main effect of stress
anova(gating_mod_cond, gating_mod_coda, test = "Chisq") # no main effect of coda
anova(gating_mod_coda, gating_mod_wm, test = "Chisq") # no effect of WM
anova(gating_mod_wm, gating_mod_int1, test = "Chisq") # no
anova(gating_mod_wm, gating_mod_int2, test = "Chisq") # no


#    AIC    BIC  logLik deviance   Chisq Chi Df Pr(>Chisq)
# 2109.0 2175.9 -1042.5   2085.0                           gating_mod_0        
# 2659.9 2741.3  -1316    2631.9  8.0819      2    0.01758 gating_mod_group *
# 2660.3 2747.5 -1315.2   2630.3  1.6123      1     0.2042 gating_mod_cond     
# 2662.2 2755.2 -1315.1   2630.2  0.131       1     0.7174 gating_mod_coda     
# 2661.7 2760.5 -1313.9   2627.7  2.4976      1      0.114 gating_mod_wm     
# 2664.1 2774.4 -1313.0   2626.1  1.6638      2     0.4352 gating_mod_int1     
# 2663.4 2768.0 -1313.7   2627.4  0.3003      1     0.5837 gating_mod_int2     

gating_mod_final <- glmer(Accuracy ~ 1 + Group  + Group:Condition + 
                        (1 + Condition * coda | participant) + 
                        (1 | Target), 
                        data = stress_gating_mod, 
                        control = glmerControl(optimizer = 'bobyqa'), 
                        family = 'binomial')

summary(gating_mod_final)
confint(gating_mod_final, method = "Wald")

# Fixed effects:
#             Estimate Std. Error z value Pr(>|z|)    
# (Intercept)  1.42635    0.17083   8.350   <2e-16 ***
# GroupLA      0.07896    0.21705   0.364   0.7160    
# GroupHS     -0.48016    0.20915  -2.296   0.0217 *  

# Relevel to test lb vs la 
stress_gating_mod$Group <- factor(stress_gating_mod$Group, levels = c("HS", "LA",  "S"))

summary(glmer(Accuracy ~ 1 + Group  + 
             (1 + Condition * coda | participant) + (1 | Target), 
             data = stress_gating_mod, control = glmerControl(optimizer = 'bobyqa'), 
             family = 'binomial'))
# beginners and LA different from each other
# Fixed effects:
#             Estimate Std. Error z value Pr(>|z|)    
# (Intercept)   0.9462     0.1623   5.829 5.56e-09 ***
# GroupLA       0.5591     0.2078   2.690  0.00714 ** 
# GroupS        0.4802     0.2091   2.296  0.02169 *  


MuMIn::r.squaredGLMM(gating_mod_final)
#       R2m        R2c 
#0.03468807 0.27453053 


# Calculate mean target fixation as a function of group, condition, 
# for each participant. We will plot the mean and calculate the 
# bootstrapped 95% confidence interval and plot it all. 
df$Group <- factor(df$Group, levels = c("LA", "HS",  "S"))

stressGatingCondP2 <- df %>%
  na.omit(.) %>%
  group_by(., Group, Condition, participant) %>%
  summarise(., meanAccuracy = mean(Accuracy)) %>%
  ggplot(., aes(x = Group, y = meanAccuracy, color = as.factor(Condition),
                dodge = Condition, 
                group = interaction(Group, Condition))) +
    geom_hline(yintercept = 0.5, lty = 2) + 
    stat_summary(fun.data = mean_cl_boot, geom = 'pointrange', fun.args = list(conf.int = 0.99),
                 position = position_dodge(width = 0.5), 
                 color = 'black', size = 0.90, show.legend = FALSE) + 
    stat_summary(fun.y = mean, geom = 'point', size = 2.90,
                 position = position_dodge(width = 0.5)) + 
    coord_cartesian(ylim = c(0.4, 1)) + 
    labs(y = '% Correct', x = 'Group', title = 'Response accuracy as a function of group', 
         caption = "Mean +/- 99% CI") + 
    scale_x_discrete(labels = c('LA', 'HS', 'SS')) + 
    scale_color_brewer(name = '', palette = 'Set1', labels = c('Paroxytone', 'Oxytone')) + 
    theme_minimal(base_size = 12, base_family = 'Times') 



# # Looks good, save as .png file
# ggsave('stressGatingP1.png', 
#        plot = stressGatingCondP2, dpi = 600, device = "png", 
#        path = "./mySources/figs/stress/s2_adv_her_nat/gating", 
#        height = 4, width = 7, units = 'in')




stressGatingCodaP2 <- df %>%
  na.omit(.) %>%
  group_by(., Group, coda, participant) %>%
  summarise(., meanAccuracy = mean(Accuracy)) %>%
  ggplot(., aes(x = Group, y = meanAccuracy, shape = as.factor(coda),
                dodge = coda, color = Group,
                group = interaction(Group, coda))) +
    stat_summary(fun.data = mean_cl_boot, geom = 'pointrange', fun.args = list(conf.int = 0.99),
                 position = position_dodge(width = 0.5), 
                 color = 'black', size = 0.90, show.legend = FALSE) + 
    stat_summary(fun.y = mean, geom = 'point', size = 2.90,
                 position = position_dodge(width = 0.5)) + 
    labs(y = '', x = 'Group', caption = '') + 
    scale_x_discrete(labels = c('LB', 'LA', 'SS')) + 
    scale_color_manual(name = '', values = c('grey90', 'grey75', 'grey55'), guide = FALSE) + 
    scale_shape_manual(name = '', values = c(16, 17), labels = c('CV', 'CVC')) + 
    scale_y_continuous(limits = c(0, 1), labels = NULL) + 
    guides(shape = guide_legend(override.aes = list(shape = c(1, 2), color = 'black'))) +
    theme_bw(base_size = 15, base_family = 'Times') +
    theme(legend.position = c(0.15, 0.14), 
          legend.box.just = "left", 
          legend.background = element_rect(fill = "transparent"), 
          legend.key = element_rect(fill = "transparent", colour = "transparent"), 
          legend.key.size = unit(0.75, 'lines'), 
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
          axis.ticks.y=element_blank(), 
          plot.margin = unit(c(7.35, 20, 7.35, 0), "pt")) + 
    annotate("text", x = 0.58, y = 1.0, label = "(b)")

# Looks good, save as .png file
ggsave('stressGatingP2.png', 
       plot = stressGatingCodaP2, dpi = 600, device = "png", 
       path = "./mySources/figs/stress/s1_beg_adv_nat/gating", 
       height = 3.5, width = 5.5, units = 'in')





my_points <- select(df, wmScaled, Accuracy) %>%
  mutate(., Accuracy = recode(Accuracy, `0` = 0.0)) %>% 
  filter(., wmScaled > -3)

# stressGatingWMp3 <- stress_gating_mod %>% 
#   ggplot(., aes(x = wmScaled, y = Accuracy)) + 
#     geom_jitter(width = 0.5, height = 0.03, alpha = 0.03, color = 'black') +
#     geom_smooth(method = 'glm', method.args = list(family = "binomial"), 
#                 se = T, fullrange = TRUE, level = 0.99, color = 'white', 
#                 fill = 'black', size = 0.75) + 
#     scale_x_continuous(expand = c(0, 0), limits = c(-2.5, 2.5)) +
#     coord_cartesian(ylim = c(-0.05, 1.05)) + 
#     labs(y = '% Correct', x = 'Working memory', caption = '+/- 99% CI.') + 
#     theme_minimal(base_size = 12, base_family = 'Times')


# ggsave('stressGatingP3.png', 
#        plot = stressGatingWMp3, dpi = 600, device = "png", 
#        path = "./mySources/figs/stress/s3_adv_int_nat/gating", 
#        height = 4, width = 7, units = 'in')



# arrange all plots together
stress_gating_plot <- grid.arrange(stressGatingLSRLp2, 
                                   stressGatingCodap3, 
                                   stressGatingWMp3, 
                                   ncol = 3)

ggsave('stressGatingAll.png', 
       plot = stress_gating_plot, dpi = 600, device = "png", 
       path = "./mySources/figs/stress/s1_beg_adv_nat/gating", 
       height = 3, width = 8.5, units = 'in')


# Story
# - t-tests
#   - Only NS and LA can predict above chance
# - GLMM
#   - Main effect of group
#   - Natives more accurate than LB, no difference with LA
#   - LA more accurate than LB
#   - marginal main effect of WM





