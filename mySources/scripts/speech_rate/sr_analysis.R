# clean working directory
rm(list = ls(all = TRUE))

library(tidyverse)

setwd("~/academia/research/in_progress/morph_pred")

sr_df <- read_csv("./mySources/stimuli/stress/output/stress_speech_rate.csv") %>% 
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