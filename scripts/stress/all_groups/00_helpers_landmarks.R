
create_baseline_landmark_posterior <- function() {

  full_posterior <- bind_rows(
    posterior_samples(mod_tw_start_full) %>%
      select(., b_Intercept:b_grouplb) %>%
      transmute(landmark = "tw_start_posterior",
                landmark_labs = "Target word\nonset",
                ss  = b_Intercept,
                hs  = b_Intercept + b_grouphs,
                int = b_Intercept + b_groupint,
                lb  = b_Intercept + b_grouplb,
                la  = b_Intercept + b_groupla),
    posterior_samples(mod_ts_v1_start_full) %>%
      select(., b_Intercept:b_grouplb) %>%
      transmute(landmark = "ts_v1_start_posterior",
                landmark_labs = "V1\nonset",
                ss  = b_Intercept,
                hs  = b_Intercept + b_grouphs,
                int = b_Intercept + b_groupint,
                lb  = b_Intercept + b_grouplb,
                la  = b_Intercept + b_groupla),
    posterior_samples(mod_ts_v1_20_full) %>%
      select(., b_Intercept:b_grouplb) %>%
      transmute(landmark = "ts_v1_20_posterior",
                landmark_labs = "20 ms\nafter V1",
                ss  = b_Intercept,
                hs  = b_Intercept + b_grouphs,
                int = b_Intercept + b_groupint,
                lb  = b_Intercept + b_grouplb,
                la  = b_Intercept + b_groupla),
    posterior_samples(mod_ts_syl1_end_full) %>%
      select(., b_Intercept:b_grouplb) %>%
      transmute(landmark = "ts_syl1_end_posterior",
                landmark_labs = "Syllable 1\noffset",
                ss  = b_Intercept,
                hs  = b_Intercept + b_grouphs,
                int = b_Intercept + b_groupint,
                lb  = b_Intercept + b_grouplb,
                la  = b_Intercept + b_groupla),
    posterior_samples(mod_ts_suffix_start_full) %>%
      select(., b_Intercept:b_grouplb) %>%
      transmute(landmark = "ts_suffix_start_posterior",
                landmark_labs = "V2\n(suffix)",
                ss  = b_Intercept,
                hs  = b_Intercept + b_grouphs,
                int = b_Intercept + b_groupint,
                lb  = b_Intercept + b_grouplb,
                la  = b_Intercept + b_groupla),
    posterior_samples(mod_ts_next_word_full) %>%
      select(., b_Intercept:b_grouplb) %>%
      transmute(landmark = "ts_next_word_posterior",
                landmark_labs = "Following\nword",
                ss  = b_Intercept,
                hs  = b_Intercept + b_grouphs,
                int = b_Intercept + b_groupint,
                lb  = b_Intercept + b_grouplb,
                la  = b_Intercept + b_groupla))

  return(full_posterior)

}


landmark_posterior_theme <- function(base_size = 16) {
  theme_minimal(base_size = base_size) +
    theme(
      plot.title = element_text(size = rel(1.5), face = "bold"),
      plot.subtitle = element_text(size = rel(1.1)),
      plot.caption = element_text(color = "#777777", vjust = 0),
      axis.title = element_text(size = rel(.9), hjust = 0.95, face = "italic"),
      panel.grid.major = element_line(size = rel(.1), color = "#000000"),
      panel.grid.minor = element_line(size = rel(.05), color = "#000000")
    )
}
