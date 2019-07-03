# Bayesian analysis -----------------------------------------------------------

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


# Load data -------------------------------------------------------------------

source(here::here("scripts", "02_load_data.R"))

# -----------------------------------------------------------------------------




# Data prep -------------------------------------------------------------------

# Set reference levels and standardize continuous variables
stress_df <- stress_lm %>%
  mutate(n_rep = 10,
         group = fct_relevel(group, "ss"),
         stress = if_else(condition == "stressed", "paroxytone", "oxytone"),
         stress = fct_relevel(stress, "paroxytone"),
         syllable = if_else(coda == 0, "CV", "CVC"),
         syllable = fct_relevel(syllable, "CV"),
         phon_prob_std = (phon_prob - mean(phon_prob)) / sd(phon_prob),
         biphon_prob_std = (biphon_prob - mean(biphon_prob)) / sd(biphon_prob),
         freq_std = (freq - mean(freq)) / sd(freq))

# -----------------------------------------------------------------------------






# Model setup -----------------------------------------------------------------

# Detect cores for model fitting
options(mc.cores = parallel::detectCores())

# Regularizing, weakly informative priors
priors <- c(
  set_prior("normal(0, 10)", class = "b")
)

# Model formula
full_model <- bf(
  targetCount | trials(n_rep) ~ group + stress + syllable +
    phon_prob_std + biphon_prob_std +
    group:stress + group:syllable + stress:syllable + group:stress:syllable +
    (1 + stress * syllable | participant) +
    (1 + group | target)
  )

# -----------------------------------------------------------------------------







# BRM at target word onset ----------------------------------------------------

mod_tw_start_full <- brm(
  formula = full_model,
  prior = priors,
  warmup = 1000, iter = 2000, chains = 4,
  family = binomial(link = "logit"),
  cores = 16,
  data = filter(stress_df, landmark_2 == "tw_start"),
  file = here("models", "stress", "s1_beg_adv_nat", "eye_track",
              "mod_tw_start_full")
)

# -----------------------------------------------------------------------------





# BRM at 1st syllable offset --------------------------------------------------

mod_ts_end_full <- update(
  object = mod_tw_start_full,
  newdata = filter(stress_df, landmark_2 == "tw_syl1_end"),
  file = here("models", "stress", "s1_beg_adv_nat", "eye_track",
              "mod_ts_end_full")
)

# -----------------------------------------------------------------------------







# Diagnostics -----------------------------------------------------------------


fixef(mod_tw_start_full) %>% inv_logit_scaled()

pp_check(mod_tw_start_full, nsamples = 200)
two_full_post  <- posterior_samples(mod_tw_start_full)
tse_full_post  <- posterior_samples(mod_ts_end_full)

hdi(two_full_post)

two_full_post %>%
  select(., b_Intercept:b_grouplb) %>%
  mutate_all(., inv_logit_scaled) %>%
  gather(parameter, estimate) %>%
  ggplot(., aes(x = estimate, y = parameter)) +
    ggridges::geom_density_ridges2(color = "white", fill = "grey80") +
    theme_minimal()

tse_full_post %>%
  select(., b_Intercept:b_grouplb) %>%
  mutate_all(., inv_logit_scaled) %>%
  gather(parameter, estimate) %>%
  ggplot(., aes(x = estimate, y = parameter)) +
    ggridges::geom_density_ridges2(color = "white", fill = "grey80") +
    theme_minimal()


# -----------------------------------------------------------------------------


full_posterior <- bind_rows(
  posterior_samples(mod_tw_start_full) %>%
    mutate(landmark = "tw_start"),
  posterior_samples(mod_ts_end_full) %>%
    mutate(landmark = "ts_end"))

full_posterior %>%
  select(., b_Intercept:b_grouplb, landmark) %>%
  mutate_if(., is.numeric, inv_logit_scaled) %>%
  gather(parameter, estimate, -landmark) %>%
  ggplot(., aes(x = estimate, y = parameter)) +
    facet_grid(. ~ landmark) +
    ggridges::geom_density_ridges2(color = "white", fill = "grey80") +
    theme_minimal()
