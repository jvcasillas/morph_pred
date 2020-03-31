# 2020/01/14: Quick and dirty for abstract that nuria needed
# To do:
# Refit simpler model with reduced groups and only WM


complete_posteriors %>%
  filter(group %in% c("ss", "la", "lb")) %>%
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
  #stat_gradientinterval(aes(shape = group, fill = group, color = group),
  #                      position = position_dodge(0.7), show.legend = F) +
  stat_pointinterval(aes(shape = group, fill = group), color = "grey10",
                     position = position_dodge(0.7)) +
  coord_cartesian(ylim = c(-5.00, 12.00)) +
  scale_shape_manual(values = c(21:25), name = NULL) +
  scale_color_brewer(palette = "Dark2", name = NULL) +
  scale_fill_brewer(palette = "Dark2", name = NULL) +
  labs(y = "Estimate", x = "Landmark", title = "Target fixation",
       subtitle = "Log odds of target fixations at diff. segmental landmarks",
       caption = "Posterior means ±95% and 66% CI") +
  landmark_posterior_theme(base_size = 12) +
  guides(fill = guide_legend(override.aes = list(size = 6))) +
  NULL



mod_ts_syl1_end_full
mod_ts_suffix_start_full

df_population <- stress_df %>%
  distinct(wm_std, n_rep, group, stress_sum, syllable_sum,
           phon_prob_std, freq_std) %>%
  mutate(participant = "fake") %>%
  add_fitted_draws(mod_ts_syl1_end_full, allow_new_levels = TRUE)

df_population %>%
  filter(group %in% c("ss", "la", "lb")) %>%
  ggplot(., aes(x = wm_std, y = .value)) +
  facet_grid(. ~ group) +
  stat_lineribbon(
    # new part
    color = "#11111100",
    .width = c(.1, .25, .5, .6, .7, .8, .9, .95)
  ) +
  #scale_x_continuous(
  #  "Days",
  #  breaks = seq(0, 9, by = 2),
  #  minor_breaks = NULL
  #) +
  scale_color_viridis_d(aesthetics = "fill") +
  guides(fill = FALSE) +
  ggtitle("Model estimates the population of participants") +
  theme_grey(base_size = 14)








# Choose settings for interaction
int_conditions_trip <- list(
  wm_std = setNames(c(-1, 0, 1), c("-1", "0", "1")),
  syllable_sum = setNames(c(-1, 1), c("CVC", "CV"))
)

# Select conditions for context_dev
conditions_trip <- make_conditions(mod_ts_syl1_end_full,
                                   vars = c("group")) %>%
  filter(group %in% c("ss", "la", "lb")) %>%
  mutate(., cond__ = fct_recode(cond__, `SS` = "group = ss",
                                `LA` = "group = la",
                                `LB` = "group = lb"))

# Get predictions
# mod_ts_syl1_end_full
# mod_ts_suffix_start_full

three_way_trip <- marginal_effects(
  mod_ts_syl1_end_full, "wm_std:syllable_sum",
  int_conditions = int_conditions_trip,
  conditions = conditions_trip,
  spaghetti = T,
  nsamples = 500
)

# Generate plot
plot_2afc_3way_trip <- plot(
  three_way_trip, plot = F, line_args = list(size = 1.5))[[1]] +
  labs(y = "%Target fixation", x = "Working memory\n(z-score)",
       title = "Model estimates",
       subtitle = "Target fixations at 1st syllable offset") +
  my_2afc_theme() +
  minimal_adj()







my_2afc_theme <- function(...) {
  list(
    scale_color_manual(name = NULL,
                       values = I(alpha(c("darkred", "darkblue"), 0.05))),
    scale_y_continuous(
      breaks = c(0, .25, 0.5, 0.75, 1.0),
      labels = c("0.0%", "25.0%", "50.0%", "75.0%", "100.0%")),
    theme_grey(base_size = 12),
    guides(colour = guide_legend(
      override.aes = list(size = 1, alpha = 1, fill = NA)))
  )
}

minimal_adj <- function(...) {
  list(
    theme(
      plot.title = element_text(size = rel(1.5), face = "bold"),
      plot.subtitle = element_text(size = rel(1.1)),
      plot.caption = element_text(color = "#777777", vjust = 0),
      axis.title = element_text(size = rel(.9), hjust = 0.95, face = "italic"),
      panel.grid.major = element_line(size = 0.15),
      panel.grid.minor = element_line(size = 0.15))
  )
}


legend_adj <- function(...) {
  list(
    theme(
      legend.position = c(0.07, 0.85),
      legend.key = element_blank(),
      legend.background = element_blank(),
      strip.background = element_blank(),
      strip.text.y = element_blank(),
      legend.text.align = 1,
      axis.title.y = element_text(size = rel(.9), hjust = 0.95),
      legend.key.size = unit(0.75, 'lines'),
      legend.title = element_text(size = 8))
  )
}
