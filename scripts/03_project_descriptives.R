# Project descriptives --------------------------------------------------------
#
# Description: get basic project descriptives for sanity checks
# Last update: 5/17/2019
#
# -----------------------------------------------------------------------------


# Source libs -----------------------------------------------------------------

source(here::here("scripts", "0_load_libs.R"))

# -----------------------------------------------------------------------------


stress <- read_csv(here("data", "clean", "stress_50ms_clean.csv"))

stress %>%
  group_by(group) %>%
  summarize(n_participant = n_distinct(participant))

