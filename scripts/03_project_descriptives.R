# Project descriptives --------------------------------------------------------
#
# Description: get basic project descriptives for sanity checks
# Last update: 5/17/2019
#
# -----------------------------------------------------------------------------


# Source libs -----------------------------------------------------------------

source(here::here("scripts", "00_load_libs.R"))

# -----------------------------------------------------------------------------


stress <- read_csv(here("data", "clean", "stress_50ms_clean.csv"))

stress %>%
  group_by(group) %>%
  summarize(n_participant = n_distinct(participant))




# Age descriptives ------------------------------------------------------------

demo_data <- read_csv(here("data", "raw", "dur_stress_demographics.csv"))


demo_data %>%
  #filter(age < 60) %>%
  group_by(group) %>%
  summarize(mean_age = mean(age), sd_age = sd(age), median_age = median(age),
            min_age = min(age), max_age = max(age))

demo_data %>%
  #filter(age < 60) %>%
  ggplot(., aes(x = group, y = age)) +
    ggbeeswarm::geom_beeswarm(alpha = 0.4) +
    stat_summary(fun.data = mean_sdl, geom = "pointrange", pch = 21, size = 1.5,
                 position = position_nudge(x = 0.1), fill = "white") +
    coord_cartesian(ylim = c(0, max(demo_data$age) + 10)) +
    labs(y = "Age", x = "Group",
         title = "Average age as a function of group",
         caption = "Mean +/- SD") +
    theme_grey(base_family = "Times", base_size = 20)

demo_data %>%
  ggplot(., aes(x = age)) +
    facet_grid(. ~ group) +
    geom_histogram(binwidth = 1, fill = "grey60", color = "black")

# -----------------------------------------------------------------------------
