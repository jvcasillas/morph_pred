
stress_df %>%
  group_by(group) %>%
  distinct(participant) %>%
  count()


# test plot with probs.
stress_df %>%
  filter(landmark_2 != "tw_coda_start",
         !(group %in% c("int", "hs"))) %>%
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



# plot 1
# Log odds estimates avg. over all factors
full_posterior_averages %>%
  pivot_longer(cols = -c("landmark", "landmark_labs"), names_to = "group",
    values_to = "estimate") %>%
  filter(!(group %in% c("int", "hs"))) %>%
  mutate(Estimate = plogis(estimate),
         Landmark = fct_relevel(landmark_labs,
           "Target word\nonset", "V1\nonset", "20 ms\nafter V1",
           "Syllable 1\noffset", "V2\n(suffix)", "Following\nword")) %>%
  ggplot(., aes(x = Landmark, y = estimate)) +
  geom_rect(data = tibble(ymin = -0.1, ymax = 0.1), inherit.aes = FALSE,
            aes(ymin = ymin, ymax = ymax, xmin = -Inf, xmax = Inf),
            fill = "lightblue", color = "white", alpha = 0.2) +
  geom_hline(yintercept = 0, lty = 3) +
  stat_gradientinterval(position = position_dodge(0.7)) +
  stat_pointinterval(color = "grey10", position = position_dodge(0.7),
    shape = 21, point_fill = "white") +
  coord_cartesian(ylim = c(-2.25, 6.25)) +
  labs(y = "Estimate", x = "\nLandmark", title = "Target fixations",
       subtitle = "Log odds of target fixations at diff. segmental landmarks",
       caption = "Posterior means ±95% and 66% CI") +
  landmark_posterior_theme(base_size = 12)

ggsave(
  filename = here("figs", "stress", "s1_beg_adv_nat", "landmarks", "targets.png"),
  plot = last_plot(), height = 3.75, width = 6.5
)


# plot 2
# Log odds estimates by group avg. over all other factors
full_posterior_averages %>%
  gather(group, estimate, -landmark, -landmark_labs) %>%
  filter(!(group %in% c("int", "hs"))) %>%
  mutate(Estimate = plogis(estimate),
         Landmark = fct_relevel(landmark_labs,
           "Target word\nonset", "V1\nonset", "20 ms\nafter V1",
           "Syllable 1\noffset", "V2\n(suffix)", "Following\nword"),
         group = fct_relevel(group, "lb", "la", "ss"),
         group = fct_relabel(group, toupper)) %>%
  ggplot(., aes(x = Landmark, y = estimate)) +
  geom_rect(data = tibble(ymin = -0.1, ymax = 0.1), inherit.aes = FALSE,
            aes(ymin = ymin, ymax = ymax, xmin = -Inf, xmax = Inf),
            fill = "lightblue", color = "white", alpha = 0.2) +
  geom_hline(yintercept = 0, lty = 3) +
  stat_gradientinterval(aes(fill = group), position = position_dodge(0.7),
    show.legend = F) +
  stat_pointinterval(aes(shape = group, fill = group), color = "grey10",
                        position = position_dodge(0.7)) +
  coord_cartesian(ylim = c(-2.25, 6.25)) +
  scale_shape_manual(values = c(21:25), name = NULL) +
  scale_color_brewer(palette = "Dark2", name = NULL) +
  scale_fill_brewer(palette = "Dark2", name = NULL) +
  labs(y = "Estimate", x = "\nLandmark", title = "Target fixations",
       subtitle = "Log odds of target fixations at diff. segmental landmarks",
       caption = "Posterior means ±95% and 66% CI") +
  landmark_posterior_theme(base_size = 12) +
  guides(fill = guide_legend(override.aes = list(size = 6)))

ggsave(
  filename = here("figs", "stress", "s1_beg_adv_nat", "landmarks", "targets_group.png"),
  plot = last_plot(), height = 3.75, width = 6.5
)

# Plot 3
# Log odds estimates by group, syll structure and lexical stress
# avg over phontactic prob and wm
complete_posteriors %>%
  filter(!(group %in% c("int", "hs"))) %>%
  mutate(
    Landmark = fct_relevel(landmark_labs,
      "Target word\nonset", "V1\nonset", "20 ms\nafter V1",
      "Syllable 1\noffset", "V2\n(suffix)", "Following\nword"),
    group = fct_relevel(group, "lb", "la", "ss"),
    group = fct_relabel(group, toupper),
    stress = fct_recode(stress, Oxytone = "oxitone", Paroxytone = "paroxitone"),
    stress = fct_relevel(stress, "Paroxytone"),
    syllable = toupper(syllable)) %>%
  ggplot(., aes(x = Landmark, y = estimate)) +
  facet_grid(syllable ~ stress) +
  geom_rect(data = tibble(ymin = -0.1, ymax = 0.1), inherit.aes = FALSE,
            aes(ymin = ymin, ymax = ymax, xmin = -Inf, xmax = Inf),
            fill = "lightblue", color = "white", alpha = 0.2) +
  geom_hline(yintercept = 0, lty = 3) +
  stat_gradientinterval(aes(fill = group),
                        position = position_dodge(0.7), show.legend = F) +
  stat_pointinterval(aes(shape = group, fill = group), color = "grey10",
                     position = position_dodge(0.7)) +
  coord_cartesian(ylim = c(-5.00, 12.00)) +
  scale_shape_manual(values = c(21:25), name = NULL) +
  scale_color_brewer(palette = "Dark2", name = NULL) +
  scale_fill_brewer(palette = "Dark2", name = NULL) +
  labs(y = "Estimate", x = "\nLandmark", title = "Target fixations",
       subtitle = "Log odds of target fixations at diff. segmental landmarks",
       caption = "Posterior means ±95% and 66% CI") +
  landmark_posterior_theme(base_size = 12) +
  guides(fill = guide_legend(override.aes = list(size = 6))) +
  NULL

ggsave(
  filename = here("figs", "stress", "s1_beg_adv_nat", "landmarks", "targets_group_stress_structure.png"),
  plot = last_plot(), height = 3.75, width = 6.5
)

# Remove irrelevant landmarks
complete_posteriors %>%
  filter(!(group %in% c("int", "hs")),
    landmark_labs %in% c("Syllable 1\noffset", "V2\n(suffix)")) %>%
  mutate(
    Landmark = fct_relevel(landmark_labs, "Syllable 1\noffset", "V2\n(suffix)"),
    group = fct_relevel(group, "lb", "la", "ss"),
    group = fct_relabel(group, toupper),
    stress = fct_recode(stress, Oxytone = "oxitone", Paroxytone = "paroxitone"),
    stress = fct_relevel(stress, "Paroxytone"),
    syllable = toupper(syllable)) %>%
  ggplot(., aes(x = Landmark, y = estimate)) +
  facet_grid(syllable ~ stress) +
  geom_rect(data = tibble(ymin = -0.1, ymax = 0.1), inherit.aes = FALSE,
            aes(ymin = ymin, ymax = ymax, xmin = -Inf, xmax = Inf),
            fill = "lightblue", color = "white", alpha = 0.2) +
  geom_hline(yintercept = 0, lty = 3) +
  stat_gradientinterval(aes(fill = group),
                        position = position_dodge(0.7), show.legend = F) +
  stat_pointinterval(aes(shape = group, fill = group), color = "grey10",
                     position = position_dodge(0.7)) +
  coord_cartesian(ylim = c(-3.00, 7.00)) +
  scale_shape_manual(values = c(21:25), name = NULL) +
  scale_color_brewer(palette = "Dark2", name = NULL) +
  scale_fill_brewer(palette = "Dark2", name = NULL) +
  labs(y = "Estimate", x = "\nLandmark", title = "Target fixations",
       subtitle = "Log odds of target fixations at diff. segmental landmarks",
       caption = "Posterior means ±95% and 66% CI") +
  landmark_posterior_theme(base_size = 8) +
  guides(fill = guide_legend(override.aes = list(size = 6))) +
  NULL

ggsave(
  filename = here("figs", "stress", "s1_beg_adv_nat", "landmarks", "targets_group_stress_structure_reduced.png"),
  plot = last_plot(), height = 3.75, width = 6.5
)

# Plot 4
# target fixations as a function of wm
# overall
# by group
# by group and stress/structure combos


conditions <- data.frame(
  stress_sum = c(-1, 1),
  syllable_sum = c(-1, 1)
  )

int_conditions <- list(
  stress_sum = setNames(c(-1, 1), c("unstressed", "stressed"))
  )

plot(
  conditional_effects(
    mod_ts_suffix_start_full,
    effects = "wm_std:group",
    #conditions = conditions,
    #int_conditions = int_conditions,
    spaghetti = T,
    nsample = 100,
    re_formula = NA
    )
  )
