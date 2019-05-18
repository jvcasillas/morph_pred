

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

# Where does target suffix begin in time course
suffixOnsets <- df_stress %>%
  group_by(., target) %>%
  summarize(., sufOnset = unique(word3_suffix) / 10)

# Center time course so that suffix onset = 0
suffixOnsetAdj <- df_stress %>%
  group_by(., target) %>%
  summarize(., sufOnsetAdj = ((unique(word3_suffix) / 10) -
                              (unique(word3_suffix) / 10)))

# Where does the target word begin in the time course?
twOnsets <- df_stress %>%
  group_by(., target) %>%
  summarize(., twOnset = unique(word3_c1v1) / 10)


# Adjust to centered time course
twOnsetAdj <- cbind(suffixOnsets, twOnsets[, 2])
twOnsetAdj <- mutate(twOnsetAdj, diff = sufOnset - twOnset, twOnsetAdj = 0 - diff)

df_stress %>%
  ggplot(., aes(x = bin, y = target)) +
  geom_path(size = 1) +
  geom_point(data = suffixOnsets, aes(x = sufOnset, y = target), color = 'red') +
  geom_point(data = twOnsets, aes(x = twOnset, y = target), color = 'blue') +
  xlab("Time course") + ylab("Items") +
  theme_bw(base_size = 12, base_family = "Times") -> onsetP

df_stress %>%
  ggplot(., aes(x = binAdj, y = target)) +
  geom_path(size = 1) +
  geom_point(data = suffixOnsetAdj, aes(x = sufOnsetAdj, y = target), color = 'red') +
  geom_point(data = twOnsetAdj, aes(x = twOnsetAdj, y = target), color = 'blue') +
  geom_vline(xintercept = binAdjMinMax, color = 'red') +
  geom_vline(xintercept = binAdjMaxMin, color = 'red') +
  scale_y_discrete(position = "right") +
  xlab("Adjusted time course") + ylab("Items") +
  theme_bw(base_size = 12, base_family = "Times") -> onsetAdjP

stressBinPlots <- plot_grid(onsetP, onsetAdjP, ncol = 2)

# ggsave('stressBinPlots.png', plot = stressBinPlots, dpi = 600, device = "png",
#          path = "./mySources/figs/stress/general", height = 6.5, width = 10, unit = "in")











df_stress_50 %>%
  group_by(group) %>%
  summarize(n = n_distinct(participant))

df_stress_50 <- read_csv("./mySources/data/clean/stressBIN50iaClean_20190513.csv") %>%
  select(., participant, group, target, condition, coda, bin, binAdj,
            binTonsetAlign, binTsuffixAlign, targetCount, distractorCount,
            targetProp, distractorProp, eLog, wts, corr) %>%
  mutate(., group = factor(group, levels = c("ss", "la", "int"))) %>%
  filter(., corr == 1,
            group %in% c('int', 'la', 'ss'),
            !participant %in% c('L01', 'L02', 'L03', 'L04', 'L05',
                                 'L06', 'L07', 'L08', 'L09', 'L10',
                                 'L15', 'L20', 'L21', 'L22', 'L23',
                                 'L26', 'L30', 'L31', 'L33', 'LA04',
                                 'LA06', 'LA07', 'LA14'),
            #binTonsetAlign >= 24 & binTonsetAlign <= 49,
            !target %in% c('cambia', 'cambió'))




condition_names <- c(
                    `stressed` = 'Paroxytone',
                    `unstressed` = 'Oxytone',
                    `0` = 'CV',
                    `1` = 'CVC'
                    )

(df_stress_50 %>%
  na.omit(.) %>%
  filter(., binTonsetAlign >= -2, binTonsetAlign <= 19) %>%
  ggplot(., aes(x = binTonsetAlign, y = targetProp, color = group, shape = group)) +
  facet_grid(. ~ coda, labeller = as_labeller(condition_names)) +
  geom_hline(yintercept = 0.5, color = 'white', size = 2) +
  geom_vline(xintercept = 11, color = 'grey40', lty = 3) +
  geom_vline(xintercept = 5, color = 'grey40', lty = 3) +
  stat_summary(fun.data = mean_cl_boot, geom = 'pointrange',  size = 0.75,
               fun.args = list(conf.int = .95, B = 1000)) +
  stat_summary(fun.y = mean, geom = 'point', color = 'white', alpha = 0.3, size = 1.75) +
  scale_shape_manual(name = "", values = 17:15, labels = c("SS", "LA", "INT")) +
  scale_color_brewer(palette = 'Set1', name = "", labels = c("SS", "LA", "INT")) +
  scale_x_continuous(breaks = c(-4, 1, 6, 11, 16), labels = c("-750", "-500", "-250", "0", "250")) +
  labs(y = 'Proportion of target fixations',
       x = 'Time relative to target syllable offset (ms)',
       caption = "Mean +/- 95% CI") +
  coord_cartesian(ylim = c(0, 1)) +
  annotate("text", x = 11.35, y = 0.02, label = 'Target syllable offset',
           angle = 90, size = 3, hjust = 0) +
  annotate("text", x = 5.35, y = 0.02, label = 'Mean target word onset',
           angle = 90, size = 3, hjust = 0) +
  theme_grey(base_size = 12, base_family = "Times") -> stressP1)

# ggsave('stressP1.png', plot = stressP1, dpi = 600, device = "png",
          # path = "./mySources/figs/stress/s3_adv_int_nat/eye_track",
          # height = 4, width = 9, unit = 'in')






# T-test plots

# We will plot the models
# This will be almost exactly the same as
# the previous plot, but it will use the
# confidence interval from the test we
# actually conducted

(stress_ttest %>%
  ungroup(.) %>%
  mutate(., group = factor(group, levels = c("la", "int", "ss")),
            condition = factor(condition, levels = c("stressed", "unstressed"))) %>%
  ggplot(., aes(x = group, y = estimate, color = condition,
                group = interaction(group, condition), dodge = condition)) +
    geom_hline(yintercept = 0.5, lty = 3) +
    geom_linerange(aes(ymin = conf.low, ymax = estimate), color = 'grey40',
                       position = position_dodge(width = 0.75), size = 1) +
    geom_point(position = position_dodge(width = 0.75), size = 4) +
    geom_point(position = position_dodge(width = 0.75), size = 2.75, color = 'grey90') +
    ylim(0.3, 1.0) +
    scale_color_brewer(palette = "Set1", name = '', labels = c('Oxytone', 'Paroxytone')) +
    scale_x_discrete(labels = c('LA', 'INT', 'SS')) +
    labs(title = 'Mean fixations at target syllable offset',
         y = 'Target fixations', x = '', caption = 'Mean and lower-bound 99% CI') +
    theme_minimal(base_size = 12, base_family = 'Times') -> stressTargetFixMOD1)

# ggsave('stressTargetFixMOD1.png', plot = stressTargetFixMOD1, dpi = 600, device = "png",
#         path = "./mySources/figs/stress/s3_adv_int_nat/eye_track",
#         height = 4, width = 7, unit = 'in')
















# GLMM plots


et_ci <- confint(prop_0_mod_final, method = "Wald", level = 0.99) %>%
  as.data.frame(.) %>%
  slice(., 13:16) %>%
  rename(., ciLow = `0.5 %`, ciHi = `99.5 %`)

stressFixModP0 <- broom::tidy(prop_0_mod_final) %>% slice(1:4) %>%
  cbind(., et_ci) %>%
  mutate(., term = recode(term, `(Intercept)` = '(Intercept)',
                          coda_sum = 'Syllable\nstructure',
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





stressFixModsHLS <- df_stress %>%
  mutate(., group = factor(group, levels = c("la", "int",  "ss")),
         coda = as.factor(coda)) %>%
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

stressFixModsP1 <- df_stress %>%
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
  labs(y = '% Correct', x = NULL, caption = '') +
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


stressFixModsP2 <- df_stress %>%
  mutate(., group = factor(group, levels = c("la", "int",  "ss")),
         coda = as.factor(coda)) %>%
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
  labs(y = '% Correct', x = NULL, caption = '+/- 99% bootstrap CI.') +
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
stress_eyet_plot <- stressFixModsP1 + stressFixModsP2

# ggsave('stress_eyet_plot.png',
#        plot = stress_eyet_plot, dpi = 600, device = "png",
#       path = "./mySources/figs/stress/s3_adv_int_nat/eye_track",
#        height = 3.5, width = 8.5, units = 'in')

