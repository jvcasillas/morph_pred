# Bayesian analysis -----------------------------------------------------------
#
# DV: response (0/1)
# IV:
#    - condition (stress: paroxytone, oxyton)
#    - coda (syllable: cv, cvc)
#    - group (ss, in, la, lb, hs)
#    - freq (lexical frequency)
#    - phon_prob (phonotactic probability)
#    - biphon_prob (biphon probability)
#    - wm (working memory)
#
# Notes: phon_prob and biphon_prob are highly correlated so we will only
#        use phon_prob.
#
# fit models at each landmark
# response ~ condition_sum + coda_sum + wm + word_frequency + phonotactic_freq
#
# Null model will test the hypothesis that participants fixate on targets
# above chance at given landmark
#  - intercept is test of independence that B != 0)
#  - this is relevant because logit == 0 == prob 0.5
# We will add predictors and assess models using LOO
#  - categorical predictors are sum coded (-1, 1) and continuous predictors
#    are standardized
#  - intercept is probability of fixating on target at average value of all
#    predictors
#
# -----------------------------------------------------------------------------


# Load data -------------------------------------------------------------------

source(here::here("scripts", "02_load_data.R"))

# -----------------------------------------------------------------------------




# Data prep -------------------------------------------------------------------

# Set reference levels and standardize continuous variables
stress_df <- stress_lm %>%
  filter(!is.na(freq), !is.na(wm)) %>%
  mutate(
    n_rep = 10,
    group = fct_relevel(group, "ss"),
    stress = if_else(condition == "stressed", "paroxytone", "oxytone"),
    stress = fct_relevel(stress, "paroxytone"),
    stress_sum = if_else(stress == "paroxytone", 1, -1),
    syllable = if_else(coda == 0, "CV", "CVC"),
    syllable = fct_relevel(syllable, "CV"),
    syllable_sum = if_else(syllable == "CV", 1, -1),
    freq_log = log(freq),
    freq_std = (freq_log - mean(freq_log)) / sd(freq_log),
    phon_prob_std = (phon_prob - mean(phon_prob)) / sd(phon_prob),
    biphon_prob_log = log(biphon_prob),
    biphon_prob_std = (biphon_prob_log - mean(biphon_prob_log)) / sd(biphon_prob_log),
    wm_std = (wm - mean(wm)) / sd(wm),
    pstm_std = (pstm - mean(pstm)) / sd(pstm)
  )

# -----------------------------------------------------------------------------




# Model setup -----------------------------------------------------------------

# Regularizing, weakly informative priors
priors <- c(
  set_prior("normal(0, 5)", class = "b")
)

priors <- c(prior(normal(0, 3), class = Intercept),
            prior(normal(0, 2), class = b),
            prior(normal(0, 1), class = sd))

# Model formula
full_model <- bf(
  targetCount | trials(n_rep) ~ group + stress_sum + syllable_sum +
  phon_prob_std + freq_std + wm_std +
  group:stress_sum + group:syllable_sum + stress_sum:syllable_sum +
  group:stress_sum:syllable_sum +
  (1 + stress_sum + syllable_sum + wm_std | participant) +
  (1 + phon_prob_std + freq_std | target)
  )

# Separate formula for model that excludes coda factor (only targets w/ codas)
coda_onset_model <- bf(
  targetCount | trials(n_rep) ~ group + stress_sum +
    phon_prob_std + freq_std + wm_std +
    group:stress_sum +
    (1 + stress_sum + wm_std | participant) +
    (1 + phon_prob_std + freq_std | target)
)

# Model formula with log biphon
full_model2 <- bf(
  targetCount | trials(n_rep) ~ group + stress_sum + syllable_sum +
  biphon_prob_std + freq_std + wm_std +
  group:stress_sum + group:syllable_sum + stress_sum:syllable_sum +
  group:biphon_prob_std +
  group:stress_sum:syllable_sum +
  (1 + stress_sum + syllable_sum + wm_std | participant) +
  (1 + phon_prob_std + freq_std | target)
  )

# -----------------------------------------------------------------------------




# Fit models ------------------------------------------------------------------

# BRM at target word onset
mod_tw_start_full <- brm(
  formula = full_model,
  prior = priors,
  warmup = 4000, iter = 8000, chains = 4,
  family = binomial(link = "logit"),
  cores = parallel::detectCores(),
  control = list(max_treedepth = 15, adapt_delta = 0.999),
  sample_prior = T,
  data = filter(stress_df, landmark_2 == "tw_start"),
  file = here("models", "stress", "all_groups", "mod_tw_start_full")
)

# BRM at tw_v1_start
mod_ts_v1_start_full <- brm(
  formula = full_model,
  prior = priors,
  warmup = 4000, iter = 8000, chains = 4,
  family = binomial(link = "logit"),
  cores = parallel::detectCores(),
  control = list(max_treedepth = 15, adapt_delta = 0.999),
  sample_prior = T,
  data = filter(stress_df, landmark_2 == "tw_v1_start"),
  file = here("models", "stress", "all_groups", "mod_ts_v1_start_full")
)

# BRM at tw_v1_20
mod_ts_v1_20_full <- brm(
  formula = full_model,
  prior = priors,
  warmup = 4000, iter = 8000, chains = 4,
  family = binomial(link = "logit"),
  cores = parallel::detectCores(),
  control = list(max_treedepth = 15, adapt_delta = 0.999),
  sample_prior = T,
  data = filter(stress_df, landmark_2 == "tw_v1_20"),
  file = here("models", "stress", "all_groups", "mod_ts_v1_20_full")
)

# BRM at tw_syl1_end
mod_ts_syl1_end_full <- brm(
  formula = full_model,
  prior = priors,
  warmup = 4000, iter = 8000, chains = 4,
  family = binomial(link = "logit"),
  cores = parallel::detectCores(),
  control = list(max_treedepth = 15, adapt_delta = 0.999),
  sample_prior = T,
  data = filter(stress_df, landmark_2 == "tw_syl1_end"),
  file = here("models", "stress", "all_groups", "mod_ts_syl1_end_full")
)

# BRM at tw_syl1_end (biphon)
mod_ts_syl1_end_full2 <- brm(
  formula = full_model2,
  prior = priors,
  warmup = 4000, iter = 8000, chains = 4,
  family = binomial(link = "logit"),
  cores = parallel::detectCores(),
  control = list(max_treedepth = 15, adapt_delta = 0.999),
  sample_prior = T,
  data = filter(stress_df, landmark_2 == "tw_syl1_end"),
  file = here("models", "stress", "all_groups", "mod_ts_syl1_end_full2")
)

# BRM at 1st syllable suffix start
mod_ts_suffix_start_full <- brm(
  formula = full_model,
  prior = priors,
  warmup = 4000, iter = 8000, chains = 4,
  family = binomial(link = "logit"),
  cores = parallel::detectCores(),
  control = list(max_treedepth = 15, adapt_delta = 0.999),
  sample_prior = T,
  data = filter(stress_df, landmark_2 == "tw_suffix_start"),
  file = here("models", "stress", "all_groups", "mod_ts_suffix_start_full")
)

# BRM at 1st syllable suffix start (w/ biphon)
mod_ts_suffix_start_full2 <- brm(
  formula = full_model2,
  prior = priors,
  warmup = 4000, iter = 8000, chains = 4,
  family = binomial(link = "logit"),
  cores = parallel::detectCores(),
  control = list(max_treedepth = 15, adapt_delta = 0.999),
  sample_prior = T,
  data = filter(stress_df, landmark_2 == "tw_suffix_start"),
  file = here("models", "stress", "all_groups", "mod_ts_suffix_start_full2")
)

# BRM at next_word
mod_ts_next_word_full <- brm(
  formula = full_model,
  prior = priors,
  warmup = 4000, iter = 8000, chains = 4,
  family = binomial(link = "logit"),
  cores = parallel::detectCores(),
  control = list(max_treedepth = 15, adapt_delta = 0.999),
  sample_prior = T,
  data = filter(stress_df, landmark_2 == "next_word"),
  file = here("models", "stress", "all_groups", "mod_ts_next_word_full")
)

###########################
# This subset has a CODA  #
# THE MODEL DOES NOT HAVE #
# CODA PREDICTOR          #
###########################

# BRM at tw_coda_start
mod_ts_coda_start_full <- brm(
  formula = coda_onset_model,
  prior = priors,
  warmup = 4000, iter = 8000, chains = 4,
  family = binomial(link = "logit"),
  cores = parallel::detectCores(),
  control = list(max_treedepth = 15, adapt_delta = 0.999),
  sample_prior = T,
  data = filter(stress_df, landmark_2 == "tw_coda_start"),
  file = here("models", "stress", "all_groups", "mod_ts_coda_start_full")
)

# -----------------------------------------------------------------------------
