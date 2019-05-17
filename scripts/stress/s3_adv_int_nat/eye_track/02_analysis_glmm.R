

# read data
df_stress <- read_csv("./mySources/data/clean/stressBIN10iaClean.csv") %>%
  filter(., corr == 1,
            group %in% c('int', 'la', 'ss'),
            !participant %in% c('L01', 'L02', 'L03', 'L04', 'L05',
                                 'L06', 'L07', 'L08', 'L09', 'L10',
                                 'L15', 'L20', 'L21', 'L22', 'L23',
                                 'L26', 'L30', 'L31', 'L33', 'LA04',
                                 'LA06', 'LA07', 'LA14'))







## @knitr stressGroupCompare


# Load wm data and combine with stress_subset_0 (proportion data)
# in order to add working memory as a covariate
wm_df <- read_csv("./mySources/data/raw/wm_all.csv") %>%
         filter(., !(group %in% c("HS", "L")),
                   !(participant %in% c('L01', 'L02', 'L03', 'L04', 'L05',
                                        'L06', 'L07', 'L08', 'L09', 'L10',
                                        'L15', 'L20', 'L21', 'L22', 'L23',
                                        'L26', 'L30', 'L31', 'L33', 'La04',
                                        'LA06', 'LA07', 'LA14')))

scale_this <- function(x) as.vector(scale(x))


# Are groups different from each other?
stress_subset_0_prop <- stress_subset_0 %>%
  na.omit(.) %>%
  select(., group, condition, participant, target, coda,
            targetCount, distractorCount, eLog, wts, targetProp) %>%
  # summarise(., meanFix = mean(targetProp)) %>%
 # ungroup(.) %>%
 # left_join(x = ., y = wm_df[, -1], by = 'participant') %>%
 mutate(., group = factor(group, levels = c("ss", "la", "int")),
           coda = as.factor(coda)) # %>%
 # group_by(., group) %>%
 # mutate(., wmScaled = scale_this(WM))

stress_subset_0_prop$coda <- as.factor(stress_subset_0_prop$coda)
stress_subset_0_prop$codaSum <- C(stress_subset_0_prop$coda, sum)
# contrasts(stress_subset_0_prop$codaSum)

stress_subset_0_prop$condition <- as.factor(stress_subset_0_prop$condition)
stress_subset_0_prop$conditionSum <- C(stress_subset_0_prop$condition, sum)
# contrasts(stress_subset_0_prop$codaSum)



# random effects building

prop_0_ranefA <- glmer(cbind(targetCount, 10 - targetCount) ~ 1 +
                    (1 | participant),
                    data = stress_subset_0_prop, REML = F, family = 'binomial',
                    control = glmerControl(optimizer = 'bobyqa'))

prop_0_ranefB <- glmer(cbind(targetCount, 10 - targetCount) ~ 1 +
                    (1 | participant) + (1 | target),
                    data = stress_subset_0_prop, family = 'binomial',
                    control = glmerControl(optimizer = 'bobyqa'))

anova(prop_0_ranefA, prop_0_ranefB, refit = F) # keep intercept for target

prop_0_ranefC <- glmer(cbind(targetCount, 10 - targetCount) ~ 1 +
                    (1 | participant) + (1 | target) +
                    (0 + conditionSum | participant),
                    data = stress_subset_0_prop, family = 'binomial',
                    control = glmerControl(optimizer = 'bobyqa'))

anova(prop_0_ranefB, prop_0_ranefC, refit = F) # keep uncorrelated slope for condition

prop_0_ranefD <- glmer(cbind(targetCount, 10 - targetCount) ~ 1 +
                    (1 | participant) + (1 | target) +
                    (0 + conditionSum | participant) + (0 + codaSum | participant),
                    data = stress_subset_0_prop, family = 'binomial',
                    control = glmerControl(optimizer = 'bobyqa'))

anova(prop_0_ranefC, prop_0_ranefD, refit = F) # Keep uncorrelated slope for coda

prop_0_ranefE <- glmer(cbind(targetCount, 10 - targetCount) ~ 1 +
                    (1 | participant) + (1 | target) +
                    (0 + conditionSum*codaSum | participant),
                    data = stress_subset_0_prop, family = 'binomial',
                    control = glmerControl(optimizer = 'bobyqa'))

anova(prop_0_ranefD, prop_0_ranefE, refit = F) # Keep condition x coda interaction slope

prop_0_ranefF <- glmer(cbind(targetCount, 10 - targetCount) ~ 1 +
                    (1 | participant) + (1 | target) +
                    (1 + conditionSum*codaSum | participant),
                    data = stress_subset_0_prop, family = 'binomial',
                    control = glmerControl(optimizer = 'bobyqa'))

anova(prop_0_ranefE, prop_0_ranefF, refit = F) # Keep correlated slope



# test fixed effects

prop_0_mod_0 <- glmer(cbind(targetCount, 10 - targetCount) ~ 1 +
                    (1 | participant) + (1 | target) +
                    (1 + conditionSum*codaSum | participant),
                    data = stress_subset_0_prop, family = 'binomial',
                    control = glmerControl(optimizer = 'bobyqa'))

prop_0_mod_group <- glmer(cbind(targetCount, 10 - targetCount) ~ 1 + group +
                    (1 | participant) + (1 | target) +
                    (1 + conditionSum*codaSum | participant),
                    data = stress_subset_0_prop, family = 'binomial',
                    control = glmerControl(optimizer = 'bobyqa'))

prop_0_mod_cond <- glmer(cbind(targetCount, 10 - targetCount) ~ 1 + group + conditionSum +
                    (1 | participant) + (1 | target) +
                    (1 + conditionSum*codaSum | participant),
                    data = stress_subset_0_prop, family = 'binomial',
                    control = glmerControl(optimizer = 'bobyqa'))

prop_0_mod_coda <- glmer(cbind(targetCount, 10 - targetCount) ~ 1 + group + codaSum +
                    (1 | participant) + (1 | target) +
                    (1 + conditionSum*codaSum | participant),
                    data = stress_subset_0_prop, family = 'binomial',
                    control = glmerControl(optimizer = 'bobyqa'))

prop_0_mod_int1 <- glmer(cbind(targetCount, 10 - targetCount) ~ 1 + group + codaSum + group:codaSum +
                    (1 | participant) + (1 | target) +
                    (1 + conditionSum*codaSum | participant),
                    data = stress_subset_0_prop, family = 'binomial',
                    control = glmerControl(optimizer = 'bobyqa'))

prop_0_mod_int2 <- glmer(cbind(targetCount, 10 - targetCount) ~ 1 + group + codaSum + group:conditionSum +
                    (1 | participant) + (1 | target) +
                    (1 + conditionSum*codaSum | participant),
                    data = stress_subset_0_prop, family = 'binomial',
                    control = glmerControl(optimizer = 'bobyqa'))

prop_0_mod_int3 <- glmer(cbind(targetCount, 10 - targetCount) ~ 1 + group + codaSum + conditionSum:codaSum +
                    (1 | participant) + (1 | target) +
                    (1 + conditionSum*codaSum | participant),
                    data = stress_subset_0_prop, family = 'binomial',
                    control = glmerControl(optimizer = 'bobyqa'))

prop_0_mod_full <- glmer(cbind(targetCount, 10 - targetCount) ~ 1 + group * conditionSum * codaSum +
                    (1 | participant) + (1 | target) +
                    (1 + conditionSum*codaSum | participant),
                    data = stress_subset_0_prop, family = 'binomial',
                    control = glmerControl(optimizer = 'bobyqa'))



anova(prop_0_mod_0, prop_0_mod_group) # main effect of group
anova(prop_0_mod_group, prop_0_mod_cond) # no effect of condition
anova(prop_0_mod_group, prop_0_mod_coda) # main effect of coda
anova(prop_0_mod_coda, prop_0_mod_int1) # no group x coda interaction
anova(prop_0_mod_coda, prop_0_mod_int2) # no group condition interaction
anova(prop_0_mod_coda, prop_0_mod_int3) # no condi x coda interaction
anova(prop_0_mod_coda, prop_0_mod_full) # no three way interaction


#    AIC    BIC  logLik deviance  Chisq Chi Df Pr(>Chisq)
# 8419.2 8490.6 -4194.6   8389.2 11.273      2   0.003565  prop_0_mod_group ***
# 8420.6 8496.8 -4194.3   8388.6 0.6178      1     0.4319  prop_0_mod_cond
# 8413.5 8489.7 -4190.7   8381.5 7.7467      1   0.005381  prop_0_mod_coda  **
# 8415.3 8501.0 -4189.7   8379.3 2.1485      2     0.3416  group x coda
# 8412.6 8503.1 -4187.3   8374.6 6.897       3    0.07525  group x condition
# 8412.2 8497.9 -4188.1   8376.2 5.2914      2    0.07096  cond x coda
# 8420.5 8534.8 -4186.3   8372.5 8.9397      8     0.3474  prop_0_mod_full


stress_subset_0_prop$group <- factor(stress_subset_0_prop$group, levels = c("ss", "la",  "int"))

prop_0_mod_final <- glmer(cbind(targetCount, 10 - targetCount) ~ 1 + group + codaSum + conditionSum +
                    (1 | participant) + (1 | target) +
                    (1 + conditionSum*codaSum | participant),
                    data = stress_subset_0_prop, family = 'binomial',
                    control = glmerControl(optimizer = 'bobyqa'))

MuMIn::r.squaredGLMM(prop_0_mod_final)

# summary(prop_0_mod_final)
# confint(prop_0_mod_final, method = "Wald")

# Fixed effects:
#               Estimate Std. Error z value Pr(>|z|)
# (Intercept)     1.6025     0.2886   5.553 2.81e-08 ***
# groupla        -1.1382     0.3274  -3.476 0.000509 ***
# groupint       -0.4926     0.4270  -1.153 0.248711
# codaSum1       -0.8177     0.2258  -3.621 0.000294 ***


# Relevel to test lb vs la
stress_subset_0_prop$group <- factor(stress_subset_0_prop$group, levels = c("int", "la",  "ss"))

summary(glmer(cbind(targetCount, 10 - targetCount) ~ 1 + group + codaSum +
       (1 | participant) + (1 | target) +
       (1 + conditionSum*codaSum | participant),
       data = stress_subset_0_prop, family = 'binomial',
       control = glmerControl(optimizer = 'bobyqa')))

# Fixed effects:
#             Estimate Std. Error z value Pr(>|z|)
# (Intercept)   1.0620     0.3840   2.765  0.00568 **
# groupla      -0.6426     0.4049  -1.587  0.11249
# groupss       0.4879     0.4244   1.150  0.25028
# codaSum1     -0.6111     0.2170  -2.816  0.00486 **




et_ci <- confint(prop_0_mod_final, method = "Wald", level = 0.99) %>%
  as.data.frame(.) %>%
  slice(., 13:16) %>%
  rename(., ciLow = `0.5 %`, ciHi = `99.5 %`)

stressFixModP0 <- broom::tidy(prop_0_mod_final) %>% slice(1:4) %>%
  cbind(., et_ci) %>%
  mutate(., term = recode(term, `(Intercept)` = '(Intercept)',
                                codaSum1 = 'Syllable\nstructure',
                                groupla = 'NIN',
                                groupint = 'IN'),
            term = factor(term, levels = c('IN',
                                           'NIN',
                                           'Syllable\nstructure',
                                           '(Intercept)'))) %>%
  ggplot(., aes(x = estimate, y = term)) +
    geom_vline(xintercept = 0, lty = 3) +
    geom_errorbarh(aes(xmin = ciLow, xmax = ciHi), height = 0.2, size = 0.65) +
    geom_point(size = 3) +
    geom_point(size = 2, color = 'lightgrey') +
    labs(y = 'Term', x = 'Estimate +/- 95% CI') +
    theme_bw(base_size = 15, base_family = 'Times') +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

# ggsave('stressP0.png', plot = stressFixModP0, dpi = 600, device = "png",
#         path = "./mySources/figs/stress/s3_adv_int_nat/eye_track",
#         height = 2.5, width = 7, unit = 'in')





stressFixModsHLS <- stress_subset_0_prop %>%
  mutate(., group = factor(group, levels = c("la", "int",  "ss"))) %>%
  group_by(., group, coda, participant) %>%
  summarise(., meanFix = mean(targetProp)) %>%
  ggplot(., aes(x = group, y = meanFix,
                group = interaction(group, coda), dodge = coda, color = coda)) +
    geom_hline(yintercept = 0.5, lty = 3) +
    stat_summary(fun.data = mean_cl_boot, geom = 'pointrange',
                 position = position_dodge(width = 0.5), color = 'black',
                 size = 0.90, show.legend = FALSE, fun.args = list(conf.int = 0.95)) +
    stat_summary(fun.y = mean, geom = 'point', size = 2.75,
                 position = position_dodge(width = 0.5), show.legend = TRUE) +
    labs(y = 'Target fixations', x = '', caption = 'Mean +/- 99% bootstrap CI.',
         title = "Mean fixations at target syllable offset") +
    coord_cartesian(ylim = c(0.4, 1)) +
    scale_x_discrete(labels = c('LA', 'INT', 'SS')) +
    scale_color_brewer(name = '', palette = 'Set1', labels = c("CV", "CVC")) +
    theme_minimal(base_size = 12, base_family = 'Times')

# ggsave('stressFixModsHLS.png', plot = stressFixModsHLS, dpi = 600, device = "png",
#         path = "./mySources/figs/stress/s3_adv_int_nat/eye_track",
#         height = 4, width = 7, unit = 'in')











# Calculate mean target fixation as a function of group, condition,
# for each participant. We will plot the mean and calculate the
# bootstrapped 95% confidence interval and plot it all.

stressFixModsP1 <- stress_subset_0_prop %>%
  mutate(., group = factor(group, levels = c("la", "int",  "ss"))) %>%
  group_by(., group, condition, participant) %>%
  summarise(., meanFix = mean(targetProp)) %>%
  ggplot(., aes(x = group, y = meanFix, shape = condition,
                group = interaction(group, condition), dodge = condition, color = group)) +
    geom_hline(yintercept = 0.5, lty = 3) +
    stat_summary(fun.data = mean_cl_boot, geom = 'pointrange',
                 position = position_dodge(width = 0.5), color = 'black',
                 size = 0.90, show.legend = FALSE, fun.args = list(conf.int = 0.99)) +
    stat_summary(fun.y = mean, geom = 'point', size = 2.75,
                 position = position_dodge(width = 0.5)) +
    ylim(0, 1) +
    labs(y = '% Correct', x = '', caption = '') +
    scale_x_discrete(labels = c('NIN', 'IN', 'M')) +
    scale_color_manual(name = '', values = c('grey90', 'grey75', 'grey55'), guide = FALSE) +
    scale_shape_manual(name = '', values = c(16, 17), labels = c('Paroxytone', 'Oxytone')) +
    guides(shape = guide_legend(override.aes = list(shape = c(1, 2), color = 'black'))) +
    theme_bw(base_size = 15, base_family = 'Times') +
    theme(legend.position = c(0.26, 0.14),
          legend.box.just = "left",
          legend.background = element_rect(fill = "transparent"),
          legend.key = element_rect(fill = "transparent", colour = "transparent"),
          legend.key.size = unit(0.75, 'lines'),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank())


stressFixModsP2 <- stress_subset_0_prop %>%
  mutate(., group = factor(group, levels = c("la", "int",  "ss"))) %>%
  group_by(., group, coda, participant) %>%
  summarise(., meanFix = mean(targetProp)) %>%
  ggplot(., aes(x = group, y = meanFix, shape = coda,
                group = interaction(group, coda), dodge = coda, color = group)) +
    geom_hline(yintercept = 0.5, lty = 3) +
    stat_summary(fun.data = mean_cl_boot, geom = 'pointrange',
                 position = position_dodge(width = 0.5), color = 'black',
                 size = 0.90, show.legend = FALSE, fun.args = list(conf.int = 0.99)) +
    stat_summary(fun.y = mean, geom = 'point', size = 2.75,
                 position = position_dodge(width = 0.5)) +
    labs(y = '% Correct', x = '', caption = '+/- 99% bootstrap CI.') +
    scale_y_continuous(position = "right", limits = c(0, 1)) +
    scale_x_discrete(labels = c('NIN', 'IN', 'M')) +
    scale_color_manual(name = '', values = c('grey90', 'grey75', 'grey55'), guide = FALSE) +
    scale_shape_manual(name = '', values = c(16, 17), labels = c('CV', 'CVC')) +
    guides(shape = guide_legend(override.aes = list(shape = c(1, 2), color = 'black'))) +
    theme_bw(base_size = 15, base_family = 'Times') +
    theme(legend.position = c(0.16, 0.14),
          legend.box.just = "left",
          legend.background = element_rect(fill = "transparent"),
          legend.key = element_rect(fill = "transparent", colour = "transparent"),
          legend.key.size = unit(0.75, 'lines'),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank())


# arrange all plots together
stress_eyet_plot <- grid.arrange(stressFixModsP1,
                                 stressFixModsP2,
                                   ncol = 2)

# ggsave('stress_eyet_plot.png',
#        plot = stress_eyet_plot, dpi = 600, device = "png",
#       path = "./mySources/figs/stress/s3_adv_int_nat/eye_track",
#        height = 3.5, width = 8.5, units = 'in')















































glimpse(df_stress)
glimpse(df_short_temp)

# calculate proportion of target fixations by sub, by target, by stress cond, by group
# for each landmark in the time course

df_timecourse <- df_short_temp %>%
  select(., participant, group, bin, BIN_START_TIME, BIN_END_TIME, target,
            startsentence:target, condition, coda, targetCount:targetProp, eLog) %>%
  mutate(., bin = bin * 10)

range(df_timecourse$bin)

df_timecourse_startsentence <- df_timecourse %>%
  filter(., startsentence + 1550 >= BIN_START_TIME, startsentence + 1550 <= BIN_END_TIME) %>%
  select(., participant, group, target, condition:eLog) %>%
  mutate(., landmark = 'start_sentence')

df_timecourse_word2_c1v1 <- df_timecourse %>%
  filter(., word2_c1v1 + 1550 >= BIN_START_TIME, word2_c1v1 + 1550 <= BIN_END_TIME) %>%
  select(., participant, group, target, condition:eLog) %>%
  mutate(., landmark = 'word2_c1v1')


df_timecourse_word3_20msafterv1 <- df_timecourse %>%
  filter(., word3_20msafterv1 + 1550 >= BIN_START_TIME, word3_20msafterv1 + 1550 <= BIN_END_TIME) %>%
  select(., participant, group, target, condition:eLog) %>%
  mutate(., landmark = 'word3_20msafterv1')


df_timecourse_word3_c1v1 <- df_timecourse %>%
  filter(., word3_c1v1 + 1550 >= BIN_START_TIME, word3_c1v1 + 1550 <= BIN_END_TIME) %>%
  select(., participant, group, target, condition:eLog) %>%
  mutate(., landmark = 'word3_c1v1')


df_timecourse_word3_c2 <- df_timecourse %>%
  filter(., word3_c2 + 1550 >= BIN_START_TIME, word3_c2 + 1550 <= BIN_END_TIME) %>%
  select(., participant, group, target, condition:eLog) %>%
  mutate(., landmark = 'word3_c2')

df_timecourse_word3_c3 <- df_timecourse %>%
  filter(., word3_c3 + 1550 >= BIN_START_TIME, word3_c3 + 1550 <= BIN_END_TIME)  %>%
  select(., participant, group, target, condition:eLog) %>%
  mutate(., landmark = 'word3_c3')

df_timecourse_word3_suffix <- df_timecourse %>%
  filter(., word3_suffix + 1550 >= BIN_START_TIME, word3_suffix + 1550 <= BIN_END_TIME) %>%
  select(., participant, group, target, condition:eLog) %>%
  mutate(., landmark = 'word3_suffix')

df_timecourse_word4_c1v1 <- df_timecourse %>%
  filter(., word4_c1v1 + 1550 >= BIN_START_TIME, word4_c1v1 + 1550 <= BIN_END_TIME) %>%
  select(., participant, group, target, condition:eLog) %>%
  mutate(., landmark = 'word4_c1v1')

df_timecourse_word5 <- df_timecourse %>%
  filter(., word5 + 1550 >= BIN_START_TIME, word5 + 1550 <= BIN_END_TIME) %>%
  select(., participant, group, target, condition:eLog) %>%
  mutate(., landmark = 'word5')

df_timecourse_end_sentence <- df_timecourse %>%
  filter(., end_sentence + 1550 >= BIN_START_TIME, end_sentence + 1550 <= BIN_END_TIME)  %>%
  select(., participant, group, target, condition:eLog) %>%
  mutate(., landmark = 'end_sentence')

df_landmarks <- do.call("rbind", list(df_timecourse_startsentence,
                      df_timecourse_word2_c1v1,
                      df_timecourse_word3_20msafterv1,
                      df_timecourse_word3_c1v1,
                      df_timecourse_word3_c2,
                      df_timecourse_word3_c3,
                      df_timecourse_word3_suffix,
                      df_timecourse_word4_c1v1,
                      df_timecourse_word5,
                      df_timecourse_end_sentence))

glimpse(df_landmarks)

df_landmarks <- mutate(df_landmarks,
  landmark = factor(landmark, levels = c('start_sentence',
                                         'word2_c1v1',
                                         'word3_c1v1',
                                         'word3_20msafterv1',
                                         'word3_c2',
                                         'word3_c3',
                                         'word3_suffix',
                                         'word4_c1v1',
                                         'word5',
                                         'end_sentence')))




df_landmarks %>%
  na.omit(.) %>%
  filter(., coda == 0,
            landmark %in% c('word3_c1v1', 'word3_20msafterv1', 'word3_suffix', 'word4_c1v1')) %>%
  group_by(., participant, target, group, coda, landmark) %>%
  summarize(., target_fix = mean(targetProp)) %>%
  ggplot(., aes(x = landmark, y = target_fix, shape = group, dodge = group)) +
    geom_hline(yintercept = 0.5, color = 'black') +
    ylim(0, 1) +
    stat_summary(fun.data = mean_cl_boot, geom = 'pointrange', position = position_dodge(width = 0.3)) +
    theme_bw()

df_landmarks %>%
  na.omit(.) %>%
  filter(., coda == 1,
            landmark %in% c('word3_c1v1', 'word3_20msafterv1', 'word3_c2', 'word3_suffix', 'word4_c1v1')) %>%
  group_by(., participant, target, group, coda, landmark) %>%
  summarize(., target_fix = mean(targetProp)) %>%
  ggplot(., aes(x = landmark, y = target_fix, shape = group, dodge = group)) +
    geom_hline(yintercept = 0.5, color = 'black') +
    ylim(0, 1) +
    stat_summary(fun.data = mean_cl_boot, geom = 'pointrange', position = position_dodge(width = 0.3)) +
    theme_bw()

