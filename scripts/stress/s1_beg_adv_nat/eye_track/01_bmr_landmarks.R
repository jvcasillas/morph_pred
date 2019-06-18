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






# BRM at target word onset ----------------------------------------------------

# Detect cores for model fitting
options(mc.cores = parallel::detectCores())

# Regularizing, weakly informative priors
priors <- c(
  set_prior("normal(0, 10)", class = "b")
)

mod_tw_start_grp <- brm(
  formula = targetCount | trials(n_rep) ~ 0 + group + stress + syllable +
    group:stress + group:syllable + stress:syllable + group:stress:syllable +
    (1 + stress * syllable | participant) +
    (1 | target),
  prior = priors,
  warmup = 1000, iter = 2000, chains = 4,
  family = binomial(link = "logit"),
  cores = 16,
  data = filter(stress_df, landmark_2 == "tw_start"),
  file = here("models", "stress", "s1_beg_adv_nat", "eye_track",
              "mod_tw_start_grp")
)


fixef(mod_tw_start_grp) %>% inv_logit_scaled()



pp_check(mod_tw_start_grp, nsamples = 200)
two_grp_post  <- posterior_samples(mod_tw_start_grp)
hdi(two_grp_post)

two_grp_post %>%
  select(., b_groupss:b_grouplb) %>%
  mutate_all(., inv_logit_scaled) %>%
  gather(parameter, estimate) %>%
  ggplot(., aes(x = estimate, y = parameter)) +
    ggridges::geom_density_ridges2(color = "white", fill = "grey80") +
    theme_minimal()

# -----------------------------------------------------------------------------





# BRM at 1st syllable offset --------------------------------------------------

mod_ts_end_null <- brm(
  formula = targetCount | trials(n_rep) ~ 1 +
    (1 | participant) +
    (1 | target),
  prior = priors,
  warmup = 1000, iter = 2000, chains = 4,
  family = binomial(link = "logit"),
  cores = 16,
  data = filter(stress_df, landmark_2 == "tw_syl1_end"),
  file = here("models", "stress", "s1_beg_adv_nat", "eye_track",
              "mod_ts_end_null")
)

mod_ts_end_grp <- update(
  object = mod_ts_end_null,
  formula = targetCount | trials(n_rep) ~ 1 + group +
    (1 | participant) +
    (1 | target),
  newdata = filter(stress_df, landmark_2 == "tw_syl1_end"),
  file = here("models", "stress", "s1_beg_adv_nat", "eye_track",
              "mod_ts_end_grp")
)

mod_ts_end_full <- update(
  object = mod_ts_end_grp,
  formula = targetCount | trials(n_rep) ~ 1 +
    group + stress + syllable + phon_prob_std + biphon_prob_std +
    group:stress + group:syllable + group:stress:syllable +
    (1 + stress + syllable | participant) +
    (1 | target),
  newdata = filter(stress_df, landmark_2 == "tw_syl1_end"),
  file = here("models", "stress", "s1_beg_adv_nat", "eye_track",
              "mod_ts_end_full")
)

pp_check(mod_tso_full, nsamples = 200)

plot(mod_tso_full)
