# Morphosyntactic predictability: tidy landmark stress data -------------------
#
# This script will load and tidy the raw eye tracking data
# with 10 ms bins for each landmark and save the output to data/clean
#
# Last update: 12/14/2019
#
# -----------------------------------------------------------------------------




# Source libs -----------------------------------------------------------------

source(here::here("scripts", "00_load_libs.R"))
word3v1 <- readRDS(here("data", "raw", "stress_word3v1_lookup.rds")) # see 110

# -----------------------------------------------------------------------------





# Tidy data -------------------------------------------------------------------

# Read data
lm10 <- read_tsv(here("data", "raw", "stress_10ms.txt")) %>%

  # Filter out unwanted experiments
  filter(., exp == "stress") %>%

  # Remove unnecessary columns
  dplyr::select(., -TRIAL_LABEL,                   -TRIAL_INDEX,
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
                            'gaba', 'grabó'),
              yes = 0, no = 1),
         corr = ifelse(correctresponse == "Lshift" & KEY_PRESSED == "Lshift" |
                       correctresponse == "Rshift" & KEY_PRESSED == "Rshift",
                       yes = 1, no = 0),
         targetProp = as.numeric(gsub(",", ".", paste(.$targetProp))),
         distractorProp = as.numeric(gsub(",", ".", paste(.$distractorProp))),
         eLog = log((targetCount + 0.5) / (10 - targetCount + 0.5)),
         wts = 1 / (targetCount + 0.5) + 1 / (10 - targetCount + 0.5)) %>%

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

  # Filter only landmarks of interest:
  # NO CODA                             CODA
  # word3_c1v1: target word onset       word3_c1v1: target word onset
  # word3_v1: onset of first vowel      word3_v1: onset of first vowel
  # word3_20msafterv1: 20 ms after v1   word3_20msafterv1: 20 ms after v1
  # word3_c2: syllable 1 offset         word3_c2: coda onset
  # word3_suffix: V2 (suffix)           word3_c3: syllable 1 offset
  # word4_c1v1: Following word          word3_suffix: V2 (suffix)
  #                                     word4_c1v1: Following word
  # NOTE: the original word3_20msafterv1 was wrong.
  #       to fix this word3_v1 was measured by hand and is added below using
  #       a lookup table. word3_20msafterv1 is then calculated manually
  mutate(., word3_v1 = word3v1[target],
            word3_20msafterv1 = word3_v1 + 20) %>%

  # Gather landmarks into single column (lm_bin) and add a 'landmark'
  # column to identify them
  # Divide lm_bin by 10 to convert time in ms to 10ms bins
  # Add 20 to lm_bin to adjust for time necessary to launch an eye movement
  gather(., landmark, lm_bin, c(word3_c1v1, word3_v1, word3_20msafterv1,
                                word3_c2, word3_c3, word3_suffix,
                                word4_c1v1)) %>%
  mutate(., lm_bin = (lm_bin / 10) %>% ceiling() + 20) %>%

  # Rename landmarks so that they make more sense
  mutate(., landmark_2 = case_when(
    landmark == "word3_c1v1" ~ "tw_start",
    landmark == "word3_v1" ~ "tw_v1_start",
    landmark == "word3_20msafterv1" ~ "tw_v1_20",
    landmark == "word3_c2" & coda == 1 ~ "tw_coda_start",
    landmark == "word3_c2" & coda == 0 ~ "tw_syl1_end",
    landmark == "word3_c3" & coda == 1 ~ "tw_syl1_end",
    landmark == "word3_suffix" ~ "tw_suffix_start",
    landmark == "word4_c1v1" ~ "next_word"),
    index = if_else(is.na(landmark_2), 1, 0)) %>%
  filter(index == 0) %>%

  # Group by participant and target, add index column to mark cases
  # where bin is equal to the landmark bin (bin == lm_bin) and then
  # filter to keep only those cases
  group_by(participant, target) %>%
  mutate(index = if_else(bin == lm_bin, 1, 0)) %>%
  filter(index == 1) %>%

  # Add covariates
  # missing frequency for some verbs
  left_join(.,
            read_csv(here("data", "raw", "phonotactic_frequency.csv")) %>%
              select(target, phon_prob, biphon_prob, freq),
            by = "target") %>%
  left_join(., read_delim(here("data", "raw", "wm_pstm_all.tsv"), delim = "\t"),
            by = "participant") %>%
  select(group:targetProp, target, verb, condition, coda, eLog:lm_bin,
         landmark_2, phon_prob:freq, pstm, wm) %>%
  write_csv(here("data", "clean", "stress_landmark_clean.csv"))

# -----------------------------------------------------------------------------



# Test plots
lm10 %>%
  filter(coda == 0, landmark_2 != "tw_v1_start",
         !(verb %in% c("grabar", "gritar"))) %>%
  mutate(., landmark_2 = fct_relevel(landmark_2,
          "tw_start", "tw_v1_20", "tw_syl1_end", "tw_suffix_start")) %>%
  ggplot(., aes(x = landmark_2, y = (targetCount / 10), color = group)) +
    facet_grid(. ~ condition) +
    geom_hline(yintercept = 0.5, color = "white", size = 3) +
    stat_summary(fun.data = mean_se, geom = "pointrange")

lm10 %>%
  filter(coda == 1, landmark_2 != "tw_v1_start",
         !(verb %in% c("grabar", "gritar"))) %>%
  mutate(., landmark_2 = fct_relevel(landmark_2,
          "tw_start", "tw_v1_20", "tw_coda_start", "tw_syl1_end", "tw_suffix_start")) %>%
  ggplot(., aes(x = landmark_2, y = (targetCount / 10), color = group)) +
    facet_grid(. ~ condition) +
    geom_hline(yintercept = 0.5, color = "white", size = 3) +
    stat_summary(fun.data = mean_se, geom = "pointrange")
