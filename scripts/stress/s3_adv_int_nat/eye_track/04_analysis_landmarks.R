
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

