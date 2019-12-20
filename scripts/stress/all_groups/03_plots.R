


# plot check
full_posterior_averages %>%
  gather(group, estimate, -landmark, -landmark_labs) %>%
  mutate(Estimate = plogis(estimate),
         Landmark = fct_relevel(landmark_labs,
           "Target word\nonset", "V1\nonset", "20 ms\nafter V1",
           "Syllable 1\noffset", "V2\n(suffix)", "Following\nword"),
         group = fct_relevel(group, "lb", "la", "int", "hs", "ss"),
         group = fct_relabel(group, toupper)) %>%
  ggplot(., aes(x = Landmark, y = estimate)) +
  geom_rect(data = tibble(ymin = -0.1, ymax = 0.1), inherit.aes = FALSE,
            aes(ymin = ymin, ymax = ymax, xmin = -Inf, xmax = Inf),
            fill = "lightblue", color = "white", alpha = 0.2) +
  geom_hline(yintercept = 0, lty = 3) +
  stat_gradientinterval(aes(shape = group, fill = group, color = group),
                        position = position_dodge(0.7)) +
  stat_pointinterval(aes(shape = group, fill = group), color = "grey10",
                        position = position_dodge(0.7)) +
  coord_cartesian(ylim = c(-2.25, 6.25)) +
  scale_shape_manual(values = c(21:25), name = NULL) +
  scale_color_brewer(palette = "Dark2", name = NULL) +
  scale_fill_brewer(palette = "Dark2", name = NULL) +
  labs(y = "Estimate", x = "Landmark", title = "Target fixation",
       subtitle = "Log odds of target fixations at diff. segmental landmarks",
       caption = "Posterior means ±95% and 66% CI") +
  landmark_posterior_theme(base_size = 12)

stress_df %>%
  filter(landmark_2 != "tw_coda_start") %>%
  mutate(landmark_2 = fct_relevel(landmark_2,
           "tw_start", "tw_v1_start", "tw_v1_20", "tw_syl1_end",
           "tw_suffix_start", "next_word"),
         group = fct_relevel(group, "lb", "la", "int", "hs", "ss")) %>%
  ggplot(., aes(x = landmark_2, y = targetProp, color = group)) +
  stat_summary(fun.data = mean_cl_boot, geom = "pointrange",
               position = position_dodge(0.5)) +
  scale_color_brewer(palette = "Dark2", name = NULL) +
  coord_cartesian(ylim = c(-0.05, 1.05)) +
  landmark_posterior_theme()
