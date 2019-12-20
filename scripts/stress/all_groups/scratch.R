stress_df %>%
  ggplot(., aes(x = biphon_prob_std, y = phon_prob_std)) +
    geom_point() +
    geom_abline(slope = 1, intercept = 0) +
    geom_smooth(method = "lm")

stress_df %>%
  na.omit() %>%
  distinct(wm_std, pstm_std) %>%
  ggplot(., aes(x = wm_std, y = pstm_std)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0) +
  geom_smooth(method = "lm")


left_join(read_csv(here("data", "raw", "dur_stress_demographics.csv")),
          read_csv(here("data", "raw", "wm_all.csv")))


  mutate(participant = toupper(id),
         wm = as.numeric(wm)) %>%
  filter(., group != "SS") %>%
  select(., participant, pstm, wm)


library(mice)
imp <- mice(mem_dat, m = 5, print = FALSE)
