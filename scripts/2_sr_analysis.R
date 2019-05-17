# Speech rate analysis --------------------------------------------------------


# Source libs -----------------------------------------------------------------

source(here::here("scripts", "0_load_libs.R"))

# -----------------------------------------------------------------------------


sr_df <- read_csv(here("stimuli", "stress", "output", "stress_speech_rate.csv")) %>%
  rename(., dur = `dur (s)`,
            phonTime = `phonationtime (s)`,
            speechRate = `speechrate (nsyll/dur)`,
            artRate = `articulation rate (nsyll / phonationtime)`,
            asd = `ASD (speakingtime/nsyll)`)

glimpse(sr_df)

mean(sr_df$speechRate)
sd(sr_df$speechRate)

mean(sr_df$dur)
sd(sr_df$dur)

ggplot(sr_df, aes(x = speechRate, y = asd)) +
  geom_point() +
  geom_smooth()
