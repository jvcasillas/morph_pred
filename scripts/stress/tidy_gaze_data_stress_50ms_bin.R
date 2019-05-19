# Morphosyntactic predictability: tidy stress data ----------------------------
#
# This script will load and tidy the raw eye tracking data
# with 50 ms bins and save the output to data/clean
#
# Last update: 05/16/2019
#
# -----------------------------------------------------------------------------




# Source libs -----------------------------------------------------------------

source(here::here("scripts", "0_load_libs.R"))

# -----------------------------------------------------------------------------





# Tidy data -------------------------------------------------------------------

# Read data
stress50 <- read_tsv(here("data", "raw", "stress_50ms.txt")) %>%

  # Filter out unwanted experiments
  filter(., exp == "stress") %>%

  # Remove unnecessary columns
  select(., -TRIAL_LABEL,                   -TRIAL_INDEX,
            -identifier,                    -sentencewav,
            -word1,                         -word2_20msafterv1,
            -word2_c2,                      -word2_c3,
            -word2_v2,                      -word2_c3,
            -word4_20msafterv1,             -word4_c1,
            -word4_suffix,                  -word6,
            -word7,                         -EYE_TRACKED,
            -IA_0_ID,                       -IA_3_ID,
            -IA_4_ID,                       -AVERAGE_IA_3_SAMPLE_COUNT,
            -AVERAGE_IA_4_SAMPLE_COUNT,     -AVERAGE_IA_0_SAMPLE_COUNT,
            -`AVERAGE_IA_0_SAMPLE_COUNT_%`, -`AVERAGE_IA_3_SAMPLE_COUNT_%`,
            -`AVERAGE_IA_4_SAMPLE_COUNT_%`) %>%

  # Rename some columns
  rename(., participant = RECORDING_SESSION_LABEL,
            bin = BIN_INDEX,
            wavID = id,
            targetCount = AVERAGE_IA_1_SAMPLE_COUNT,
            distractorCount = AVERAGE_IA_2_SAMPLE_COUNT,
            targetProp = `AVERAGE_IA_1_SAMPLE_COUNT_%`,
            distractorProp = `AVERAGE_IA_2_SAMPLE_COUNT_%`) %>%

  # Remove weird chars from target list
  # Create condition variable (unstressed/stressed)
  # Create coda variable (0 = no coda, 1 = coda)
  # Create corr variable (0 = incorrect, 1 = correct)
  # Change ',' to '.' in proportion columns
  # Create eLog variable and respective wts
  mutate(., target = as.factor(target),
         target = gsub("\u0097", "o", paste(.$target)),
         condition = ifelse(
           target %in% c('bebió', 'cambió', 'cantó', 'comió', 'compró',
                         'firmó', 'grabó', 'guardó', 'lanzó', 'lavó',
                         'llenó', 'mandó', 'pintó', 'rompió', 'sacó',
                         'subió'),
           yes = "unstressed", no = "stressed"),
         coda = ifelse(
           target %in% c('bebe', 'bebió', 'llena', 'llenó', 'sube', 'subió',
                         'come', 'comió', 'saca', 'sacó', 'lava', 'lavó',
                         'graba', 'grabó'),
           yes = 0, no = 1),
         corr = ifelse(correctresponse == "Lshift" & KEY_PRESSED == "Lshift" |
                         correctresponse == "Rshift" & KEY_PRESSED == "Rshift",
                       yes = 1, no = 0),
         targetProp = as.numeric(gsub(",", ".", paste(.$targetProp))),
         distractorProp = as.numeric(gsub(",", ".", paste(.$distractorProp))),
         eLog = log((targetCount + 0.5) / (50 - targetCount + 0.5)),
         wts = 1 / (targetCount + 0.5) + 1 / (50 - targetCount + 0.5)) %>%

  # Create 'group' column and new 'id'
  # in order to match participant ids with WM df
  separate(., col = participant, into = c("group", "id"), sep = -2, remove = TRUE) %>%

  # Recode groups that have random labels
  # create new participant id labels
  mutate(., group = recode(group, l = 'lb', s = 'ss', `in` = 'int'),
         group_temp = recode(group, lb = 'L', ss = "SB", ss = "SBs",
                             int = "IN", hs = "HS", la = "LA"),
         target = recode(target, gaba = 'graba')) %>%

  # Combine group_temp and id to
  # complete recode of participant IDs (now consistent with WM)
  unite(., col = participant, group_temp, id, sep = "", remove = TRUE) %>%

  # Filter out incorrect responses
  filter(corr == 1) %>%

  # Select necessary columns
  # Gather data to prepare for bin adjustment
  # Get suffix onset label and center at 0 for each
  # participant for each item
  dplyr::select(participant, group, verb, target, condition, coda, target, bin,
                targetCount, targetProp, eLog, wts, word3_c2, word3_c3) %>%
  gather(., landmark, lm_bin, -c(participant:wts)) %>%
  mutate(., lm_bin = (lm_bin / 50) %>% ceiling(.),
         t_onset = if_else(bin == lm_bin, TRUE, FALSE)) %>%
  # Must remove word3_c2 for coda words (its the start of the coda, unnecessary)
  # and remove word3_c3 for non coda words (they dont have c3, lm_bin = 0)
  filter(coda == 1 & landmark == "word3_c3" | coda == 0 & lm_bin != 0) %>%
  group_by(., participant, target) %>%
  mutate(., time_zero = onset_pupil(bin, t_onset, event = c("TRUE"))) %>%
  ungroup(.) %>%
  write_csv(here("data", "clean", "stress_50ms_clean.csv"))

# -----------------------------------------------------------------------------



# Test plot
stress50 %>%
  filter(time_zero > -10) %>%
  ggplot(., aes(x = time_zero, y = targetProp, color = group)) +
  facet_grid(coda ~ condition) +
  geom_vline(xintercept = 4, lty = 3) +
  geom_hline(yintercept = 0.5, color = "white", size = 3) +
  stat_summary(fun.y = mean, geom = "line")



