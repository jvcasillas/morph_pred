# Project descriptives --------------------------------------------------------
#
# Description: get basic project descriptives for sanity checks
# Last update: 6/03/2019
#
# -----------------------------------------------------------------------------


# Source libs -----------------------------------------------------------------

source(here::here("scripts", "00_load_libs.R"))

# -----------------------------------------------------------------------------


stress <- read_csv(here("data", "clean", "stress_50ms_clean.csv"))
View(stress)
stress %>%
  group_by(group) %>%
  summarize(n_participant = n_distinct(participant))




# Age descriptives ------------------------------------------------------------

demo_data <- read_csv(here("data", "raw", "dur_stress_demographics.csv"))

# Remove participant IN17 because we lost eye-tracking data
# (file was corrupted)
# We remove LA09 and LA15 to make the groups homogenous in L2
# proficiency (DELE)

demo_data <- demo_data %>%
  filter(., id != "LA09" & id != "LA15" & id != "IN17")

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


# language experience (years of exposure to the L2)
# we calculate it by substracting age of acquistion (aoa_l2) from their actual age (age)

demo_data <- demo_data %>%
  mutate(., lang_exp = age - aoa_l2)

demo_data %>%
  group_by(., group) %>%
  summarise(., mean_lang_exp = mean(lang_exp),
            sd_lang_exp = sd(lang_exp))


# -----------------------------------------------------------------------------
# time abroad (both learner groups together)

demo_data %>%
  summarise(., max = max(months_abroad), min = min(months_abroad),
           mean = round(mean(months_abroad),2),sd = round(sd(months_abroad),2))

  #     max   min  mean    sd
  # 1   418     0  22.7    60.8

# -----------------------------------------------------------------------------
# years of work for interpreters

# change , to .

demo_data$years_work_int <- as.numeric(demo_data$years_work_int)

demo_data %>%
  filter(., group == "in") %>%
  summarise(., min = min(years_work_int), max = max(years_work_int),
            mean = mean(years_work_int),sd = round(sd(years_work_int),2))

  # min   max  mean    sd
  # 2    35   14.2   9.23

# -----------------------------------------------------------------------------
# Table for BLC manuscript with age and DELE
#
glimpse(demo_data)
demo_data$wm <- as.numeric(demo_data$wm)
demo_data %>%
  group_by(., group) %>%
  summarise(.,  n = n_distinct(id),
            age_mean = round(mean(age),2), age_sd = round(sd(age),2),
            dele_mean = round(mean(dele),2), dele_sd = round(sd(dele),2),
            wm_mn = round(mean(wm, na.rm = TRUE),2), wm_sd = round(sd(wm, na.rm = TRUE),2)) %>% knitr::kable()



