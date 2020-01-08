# Tidy posterior distributions ------------------------------------------------
#
# Notes:
# The intercept of each model is the posterior probability of native speaking
# monolinguals of Spanish (ss) fixating on targets at a specific landmark,
# controlling for syllable structure (CV, CVC), lexical stress (oxytone,
# paroxytone), phontactic probability, lexical frequency, and working memory.
#
# This script:
# This script will tidy the posterior distributions of the aforementioned
# models so that we have the probability of fixating on targets for each group
# in each condition controlling for phonotactic probability, lexical
# frequency, and working memory, i.e., when they are at their average values.
#
# Last update: Mon Dec 16 10:30:23 2019
# -----------------------------------------------------------------------------




# Tidy posterior --------------------------------------------------------------
#
# Get posterior distribution of each model
# Calculate prob. for each group in each condition
full_posterior_averages <- create_baseline_landmark_posterior()
ts_coda_start_posterior <- posterior_samples(mod_ts_coda_start_full)




tw_start_posterior <- posterior_samples(mod_tw_start_full) %>%
  select(., starts_with("b_")) %>%
  transmute(
    ss_paroxitone_cv   = b_Intercept + b_stress_sum + b_syllable_sum + `b_stress_sum:syllable_sum`,
    ss_oxitone_cv      = b_Intercept - b_stress_sum + b_syllable_sum - `b_stress_sum:syllable_sum`,
    ss_paroxitone_cvc  = b_Intercept + b_stress_sum - b_syllable_sum - `b_stress_sum:syllable_sum`,
    ss_oxitone_cvc     = b_Intercept - b_stress_sum - b_syllable_sum + `b_stress_sum:syllable_sum`,

    la_paroxitone_cv   = b_Intercept + b_stress_sum + b_syllable_sum + `b_stress_sum:syllable_sum` + b_groupla + `b_groupla:stress_sum` + `b_groupla:syllable_sum` + `b_groupla:stress_sum:syllable_sum`,
    la_oxitone_cv      = b_Intercept - b_stress_sum + b_syllable_sum - `b_stress_sum:syllable_sum` + b_groupla - `b_groupla:stress_sum` + `b_groupla:syllable_sum` - `b_groupla:stress_sum:syllable_sum`,
    la_paroxitone_cvc  = b_Intercept + b_stress_sum - b_syllable_sum - `b_stress_sum:syllable_sum` + b_groupla + `b_groupla:stress_sum` - `b_groupla:syllable_sum` - `b_groupla:stress_sum:syllable_sum`,
    la_oxitone_cvc     = b_Intercept - b_stress_sum - b_syllable_sum + `b_stress_sum:syllable_sum` + b_groupla - `b_groupla:stress_sum` - `b_groupla:syllable_sum` + `b_groupla:stress_sum:syllable_sum`,

    lb_paroxitone_cv   = b_Intercept + b_stress_sum + b_syllable_sum + `b_stress_sum:syllable_sum` + b_grouplb + `b_grouplb:stress_sum` + `b_grouplb:syllable_sum` + `b_grouplb:stress_sum:syllable_sum`,
    lb_oxitone_cv      = b_Intercept - b_stress_sum + b_syllable_sum - `b_stress_sum:syllable_sum` + b_grouplb - `b_grouplb:stress_sum` + `b_grouplb:syllable_sum` - `b_grouplb:stress_sum:syllable_sum`,
    lb_paroxitone_cvc  = b_Intercept + b_stress_sum - b_syllable_sum - `b_stress_sum:syllable_sum` + b_grouplb + `b_grouplb:stress_sum` - `b_grouplb:syllable_sum` - `b_grouplb:stress_sum:syllable_sum`,
    lb_oxitone_cvc     = b_Intercept - b_stress_sum - b_syllable_sum + `b_stress_sum:syllable_sum` + b_grouplb - `b_grouplb:stress_sum` - `b_grouplb:syllable_sum` + `b_grouplb:stress_sum:syllable_sum`,

    hs_paroxitone_cv   = b_Intercept + b_stress_sum + b_syllable_sum + `b_stress_sum:syllable_sum` + b_grouphs + `b_grouphs:stress_sum` + `b_grouphs:syllable_sum` + `b_grouphs:stress_sum:syllable_sum`,
    hs_oxitone_cv      = b_Intercept - b_stress_sum + b_syllable_sum - `b_stress_sum:syllable_sum` + b_grouphs - `b_grouphs:stress_sum` + `b_grouphs:syllable_sum` - `b_grouphs:stress_sum:syllable_sum`,
    hs_paroxitone_cvc  = b_Intercept + b_stress_sum - b_syllable_sum - `b_stress_sum:syllable_sum` + b_grouphs + `b_grouphs:stress_sum` - `b_grouphs:syllable_sum` - `b_grouphs:stress_sum:syllable_sum`,
    hs_oxitone_cvc     = b_Intercept - b_stress_sum - b_syllable_sum + `b_stress_sum:syllable_sum` + b_grouphs - `b_grouphs:stress_sum` - `b_grouphs:syllable_sum` + `b_grouphs:stress_sum:syllable_sum`,

    int_paroxitone_cv  = b_Intercept + b_stress_sum + b_syllable_sum + `b_stress_sum:syllable_sum` + b_groupint + `b_groupint:stress_sum` + `b_groupint:syllable_sum` + `b_groupint:stress_sum:syllable_sum`,
    int_oxitone_cv     = b_Intercept - b_stress_sum + b_syllable_sum - `b_stress_sum:syllable_sum` + b_groupint - `b_groupint:stress_sum` + `b_groupint:syllable_sum` - `b_groupint:stress_sum:syllable_sum`,
    int_paroxitone_cvc = b_Intercept + b_stress_sum - b_syllable_sum - `b_stress_sum:syllable_sum` + b_groupint + `b_groupint:stress_sum` - `b_groupint:syllable_sum` - `b_groupint:stress_sum:syllable_sum`,
    int_oxitone_cvc    = b_Intercept - b_stress_sum - b_syllable_sum + `b_stress_sum:syllable_sum` + b_groupint - `b_groupint:stress_sum` - `b_groupint:syllable_sum` + `b_groupint:stress_sum:syllable_sum`
    ) %>%
  gather(group, estimate) %>%
  separate(group, into = c("group", "stress", "syllable")) %>%
  mutate(landmark = "tw_start_posterior",
         landmark_labs = "Target word\nonset")


ts_v1_start_posterior <- posterior_samples(mod_ts_v1_start_full) %>%
  select(., starts_with("b_")) %>%
  transmute(
    ss_paroxitone_cv   = b_Intercept + b_stress_sum + b_syllable_sum + `b_stress_sum:syllable_sum`,
    ss_oxitone_cv      = b_Intercept - b_stress_sum + b_syllable_sum - `b_stress_sum:syllable_sum`,
    ss_paroxitone_cvc  = b_Intercept + b_stress_sum - b_syllable_sum - `b_stress_sum:syllable_sum`,
    ss_oxitone_cvc     = b_Intercept - b_stress_sum - b_syllable_sum + `b_stress_sum:syllable_sum`,

    la_paroxitone_cv   = b_Intercept + b_stress_sum + b_syllable_sum + `b_stress_sum:syllable_sum` + b_groupla + `b_groupla:stress_sum` + `b_groupla:syllable_sum` + `b_groupla:stress_sum:syllable_sum`,
    la_oxitone_cv      = b_Intercept - b_stress_sum + b_syllable_sum - `b_stress_sum:syllable_sum` + b_groupla - `b_groupla:stress_sum` + `b_groupla:syllable_sum` - `b_groupla:stress_sum:syllable_sum`,
    la_paroxitone_cvc  = b_Intercept + b_stress_sum - b_syllable_sum - `b_stress_sum:syllable_sum` + b_groupla + `b_groupla:stress_sum` - `b_groupla:syllable_sum` - `b_groupla:stress_sum:syllable_sum`,
    la_oxitone_cvc     = b_Intercept - b_stress_sum - b_syllable_sum + `b_stress_sum:syllable_sum` + b_groupla - `b_groupla:stress_sum` - `b_groupla:syllable_sum` + `b_groupla:stress_sum:syllable_sum`,

    lb_paroxitone_cv   = b_Intercept + b_stress_sum + b_syllable_sum + `b_stress_sum:syllable_sum` + b_grouplb + `b_grouplb:stress_sum` + `b_grouplb:syllable_sum` + `b_grouplb:stress_sum:syllable_sum`,
    lb_oxitone_cv      = b_Intercept - b_stress_sum + b_syllable_sum - `b_stress_sum:syllable_sum` + b_grouplb - `b_grouplb:stress_sum` + `b_grouplb:syllable_sum` - `b_grouplb:stress_sum:syllable_sum`,
    lb_paroxitone_cvc  = b_Intercept + b_stress_sum - b_syllable_sum - `b_stress_sum:syllable_sum` + b_grouplb + `b_grouplb:stress_sum` - `b_grouplb:syllable_sum` - `b_grouplb:stress_sum:syllable_sum`,
    lb_oxitone_cvc     = b_Intercept - b_stress_sum - b_syllable_sum + `b_stress_sum:syllable_sum` + b_grouplb - `b_grouplb:stress_sum` - `b_grouplb:syllable_sum` + `b_grouplb:stress_sum:syllable_sum`,

    hs_paroxitone_cv   = b_Intercept + b_stress_sum + b_syllable_sum + `b_stress_sum:syllable_sum` + b_grouphs + `b_grouphs:stress_sum` + `b_grouphs:syllable_sum` + `b_grouphs:stress_sum:syllable_sum`,
    hs_oxitone_cv      = b_Intercept - b_stress_sum + b_syllable_sum - `b_stress_sum:syllable_sum` + b_grouphs - `b_grouphs:stress_sum` + `b_grouphs:syllable_sum` - `b_grouphs:stress_sum:syllable_sum`,
    hs_paroxitone_cvc  = b_Intercept + b_stress_sum - b_syllable_sum - `b_stress_sum:syllable_sum` + b_grouphs + `b_grouphs:stress_sum` - `b_grouphs:syllable_sum` - `b_grouphs:stress_sum:syllable_sum`,
    hs_oxitone_cvc     = b_Intercept - b_stress_sum - b_syllable_sum + `b_stress_sum:syllable_sum` + b_grouphs - `b_grouphs:stress_sum` - `b_grouphs:syllable_sum` + `b_grouphs:stress_sum:syllable_sum`,

    int_paroxitone_cv  = b_Intercept + b_stress_sum + b_syllable_sum + `b_stress_sum:syllable_sum` + b_groupint + `b_groupint:stress_sum` + `b_groupint:syllable_sum` + `b_groupint:stress_sum:syllable_sum`,
    int_oxitone_cv     = b_Intercept - b_stress_sum + b_syllable_sum - `b_stress_sum:syllable_sum` + b_groupint - `b_groupint:stress_sum` + `b_groupint:syllable_sum` - `b_groupint:stress_sum:syllable_sum`,
    int_paroxitone_cvc = b_Intercept + b_stress_sum - b_syllable_sum - `b_stress_sum:syllable_sum` + b_groupint + `b_groupint:stress_sum` - `b_groupint:syllable_sum` - `b_groupint:stress_sum:syllable_sum`,
    int_oxitone_cvc    = b_Intercept - b_stress_sum - b_syllable_sum + `b_stress_sum:syllable_sum` + b_groupint - `b_groupint:stress_sum` - `b_groupint:syllable_sum` + `b_groupint:stress_sum:syllable_sum`
  ) %>%
  gather(group, estimate) %>%
  separate(group, into = c("group", "stress", "syllable")) %>%
  mutate(landmark = "ts_v1_start_posterior",
         landmark_labs = "V1\nonset")


ts_v1_20_posterior <- posterior_samples(mod_ts_v1_20_full) %>%
  select(., starts_with("b_")) %>%
  transmute(
    ss_paroxitone_cv   = b_Intercept + b_stress_sum + b_syllable_sum + `b_stress_sum:syllable_sum`,
    ss_oxitone_cv      = b_Intercept - b_stress_sum + b_syllable_sum - `b_stress_sum:syllable_sum`,
    ss_paroxitone_cvc  = b_Intercept + b_stress_sum - b_syllable_sum - `b_stress_sum:syllable_sum`,
    ss_oxitone_cvc     = b_Intercept - b_stress_sum - b_syllable_sum + `b_stress_sum:syllable_sum`,

    la_paroxitone_cv   = b_Intercept + b_stress_sum + b_syllable_sum + `b_stress_sum:syllable_sum` + b_groupla + `b_groupla:stress_sum` + `b_groupla:syllable_sum` + `b_groupla:stress_sum:syllable_sum`,
    la_oxitone_cv      = b_Intercept - b_stress_sum + b_syllable_sum - `b_stress_sum:syllable_sum` + b_groupla - `b_groupla:stress_sum` + `b_groupla:syllable_sum` - `b_groupla:stress_sum:syllable_sum`,
    la_paroxitone_cvc  = b_Intercept + b_stress_sum - b_syllable_sum - `b_stress_sum:syllable_sum` + b_groupla + `b_groupla:stress_sum` - `b_groupla:syllable_sum` - `b_groupla:stress_sum:syllable_sum`,
    la_oxitone_cvc     = b_Intercept - b_stress_sum - b_syllable_sum + `b_stress_sum:syllable_sum` + b_groupla - `b_groupla:stress_sum` - `b_groupla:syllable_sum` + `b_groupla:stress_sum:syllable_sum`,

    lb_paroxitone_cv   = b_Intercept + b_stress_sum + b_syllable_sum + `b_stress_sum:syllable_sum` + b_grouplb + `b_grouplb:stress_sum` + `b_grouplb:syllable_sum` + `b_grouplb:stress_sum:syllable_sum`,
    lb_oxitone_cv      = b_Intercept - b_stress_sum + b_syllable_sum - `b_stress_sum:syllable_sum` + b_grouplb - `b_grouplb:stress_sum` + `b_grouplb:syllable_sum` - `b_grouplb:stress_sum:syllable_sum`,
    lb_paroxitone_cvc  = b_Intercept + b_stress_sum - b_syllable_sum - `b_stress_sum:syllable_sum` + b_grouplb + `b_grouplb:stress_sum` - `b_grouplb:syllable_sum` - `b_grouplb:stress_sum:syllable_sum`,
    lb_oxitone_cvc     = b_Intercept - b_stress_sum - b_syllable_sum + `b_stress_sum:syllable_sum` + b_grouplb - `b_grouplb:stress_sum` - `b_grouplb:syllable_sum` + `b_grouplb:stress_sum:syllable_sum`,

    hs_paroxitone_cv   = b_Intercept + b_stress_sum + b_syllable_sum + `b_stress_sum:syllable_sum` + b_grouphs + `b_grouphs:stress_sum` + `b_grouphs:syllable_sum` + `b_grouphs:stress_sum:syllable_sum`,
    hs_oxitone_cv      = b_Intercept - b_stress_sum + b_syllable_sum - `b_stress_sum:syllable_sum` + b_grouphs - `b_grouphs:stress_sum` + `b_grouphs:syllable_sum` - `b_grouphs:stress_sum:syllable_sum`,
    hs_paroxitone_cvc  = b_Intercept + b_stress_sum - b_syllable_sum - `b_stress_sum:syllable_sum` + b_grouphs + `b_grouphs:stress_sum` - `b_grouphs:syllable_sum` - `b_grouphs:stress_sum:syllable_sum`,
    hs_oxitone_cvc     = b_Intercept - b_stress_sum - b_syllable_sum + `b_stress_sum:syllable_sum` + b_grouphs - `b_grouphs:stress_sum` - `b_grouphs:syllable_sum` + `b_grouphs:stress_sum:syllable_sum`,

    int_paroxitone_cv  = b_Intercept + b_stress_sum + b_syllable_sum + `b_stress_sum:syllable_sum` + b_groupint + `b_groupint:stress_sum` + `b_groupint:syllable_sum` + `b_groupint:stress_sum:syllable_sum`,
    int_oxitone_cv     = b_Intercept - b_stress_sum + b_syllable_sum - `b_stress_sum:syllable_sum` + b_groupint - `b_groupint:stress_sum` + `b_groupint:syllable_sum` - `b_groupint:stress_sum:syllable_sum`,
    int_paroxitone_cvc = b_Intercept + b_stress_sum - b_syllable_sum - `b_stress_sum:syllable_sum` + b_groupint + `b_groupint:stress_sum` - `b_groupint:syllable_sum` - `b_groupint:stress_sum:syllable_sum`,
    int_oxitone_cvc    = b_Intercept - b_stress_sum - b_syllable_sum + `b_stress_sum:syllable_sum` + b_groupint - `b_groupint:stress_sum` - `b_groupint:syllable_sum` + `b_groupint:stress_sum:syllable_sum`
  ) %>%
  gather(group, estimate) %>%
  separate(group, into = c("group", "stress", "syllable")) %>%
  mutate(landmark = "ts_v1_20_posterior",
         landmark_labs = "20 ms\nafter V1")


ts_syl1_end_posterior <- posterior_samples(mod_ts_syl1_end_full) %>%
  select(., starts_with("b_")) %>%
  transmute(
    ss_paroxitone_cv   = b_Intercept + b_stress_sum + b_syllable_sum + `b_stress_sum:syllable_sum`,
    ss_oxitone_cv      = b_Intercept - b_stress_sum + b_syllable_sum - `b_stress_sum:syllable_sum`,
    ss_paroxitone_cvc  = b_Intercept + b_stress_sum - b_syllable_sum - `b_stress_sum:syllable_sum`,
    ss_oxitone_cvc     = b_Intercept - b_stress_sum - b_syllable_sum + `b_stress_sum:syllable_sum`,

    la_paroxitone_cv   = b_Intercept + b_stress_sum + b_syllable_sum + `b_stress_sum:syllable_sum` + b_groupla + `b_groupla:stress_sum` + `b_groupla:syllable_sum` + `b_groupla:stress_sum:syllable_sum`,
    la_oxitone_cv      = b_Intercept - b_stress_sum + b_syllable_sum - `b_stress_sum:syllable_sum` + b_groupla - `b_groupla:stress_sum` + `b_groupla:syllable_sum` - `b_groupla:stress_sum:syllable_sum`,
    la_paroxitone_cvc  = b_Intercept + b_stress_sum - b_syllable_sum - `b_stress_sum:syllable_sum` + b_groupla + `b_groupla:stress_sum` - `b_groupla:syllable_sum` - `b_groupla:stress_sum:syllable_sum`,
    la_oxitone_cvc     = b_Intercept - b_stress_sum - b_syllable_sum + `b_stress_sum:syllable_sum` + b_groupla - `b_groupla:stress_sum` - `b_groupla:syllable_sum` + `b_groupla:stress_sum:syllable_sum`,

    lb_paroxitone_cv   = b_Intercept + b_stress_sum + b_syllable_sum + `b_stress_sum:syllable_sum` + b_grouplb + `b_grouplb:stress_sum` + `b_grouplb:syllable_sum` + `b_grouplb:stress_sum:syllable_sum`,
    lb_oxitone_cv      = b_Intercept - b_stress_sum + b_syllable_sum - `b_stress_sum:syllable_sum` + b_grouplb - `b_grouplb:stress_sum` + `b_grouplb:syllable_sum` - `b_grouplb:stress_sum:syllable_sum`,
    lb_paroxitone_cvc  = b_Intercept + b_stress_sum - b_syllable_sum - `b_stress_sum:syllable_sum` + b_grouplb + `b_grouplb:stress_sum` - `b_grouplb:syllable_sum` - `b_grouplb:stress_sum:syllable_sum`,
    lb_oxitone_cvc     = b_Intercept - b_stress_sum - b_syllable_sum + `b_stress_sum:syllable_sum` + b_grouplb - `b_grouplb:stress_sum` - `b_grouplb:syllable_sum` + `b_grouplb:stress_sum:syllable_sum`,

    hs_paroxitone_cv   = b_Intercept + b_stress_sum + b_syllable_sum + `b_stress_sum:syllable_sum` + b_grouphs + `b_grouphs:stress_sum` + `b_grouphs:syllable_sum` + `b_grouphs:stress_sum:syllable_sum`,
    hs_oxitone_cv      = b_Intercept - b_stress_sum + b_syllable_sum - `b_stress_sum:syllable_sum` + b_grouphs - `b_grouphs:stress_sum` + `b_grouphs:syllable_sum` - `b_grouphs:stress_sum:syllable_sum`,
    hs_paroxitone_cvc  = b_Intercept + b_stress_sum - b_syllable_sum - `b_stress_sum:syllable_sum` + b_grouphs + `b_grouphs:stress_sum` - `b_grouphs:syllable_sum` - `b_grouphs:stress_sum:syllable_sum`,
    hs_oxitone_cvc     = b_Intercept - b_stress_sum - b_syllable_sum + `b_stress_sum:syllable_sum` + b_grouphs - `b_grouphs:stress_sum` - `b_grouphs:syllable_sum` + `b_grouphs:stress_sum:syllable_sum`,

    int_paroxitone_cv  = b_Intercept + b_stress_sum + b_syllable_sum + `b_stress_sum:syllable_sum` + b_groupint + `b_groupint:stress_sum` + `b_groupint:syllable_sum` + `b_groupint:stress_sum:syllable_sum`,
    int_oxitone_cv     = b_Intercept - b_stress_sum + b_syllable_sum - `b_stress_sum:syllable_sum` + b_groupint - `b_groupint:stress_sum` + `b_groupint:syllable_sum` - `b_groupint:stress_sum:syllable_sum`,
    int_paroxitone_cvc = b_Intercept + b_stress_sum - b_syllable_sum - `b_stress_sum:syllable_sum` + b_groupint + `b_groupint:stress_sum` - `b_groupint:syllable_sum` - `b_groupint:stress_sum:syllable_sum`,
    int_oxitone_cvc    = b_Intercept - b_stress_sum - b_syllable_sum + `b_stress_sum:syllable_sum` + b_groupint - `b_groupint:stress_sum` - `b_groupint:syllable_sum` + `b_groupint:stress_sum:syllable_sum`
  ) %>%
  gather(group, estimate) %>%
  separate(group, into = c("group", "stress", "syllable")) %>%
  mutate(landmark = "ts_syl1_end_posterior",
         landmark_labs = "Syllable 1\noffset")


ts_suffix_start_posterior <- posterior_samples(mod_ts_suffix_start_full) %>%
  select(., starts_with("b_")) %>%
  transmute(
    ss_paroxitone_cv   = b_Intercept + b_stress_sum + b_syllable_sum + `b_stress_sum:syllable_sum`,
    ss_oxitone_cv      = b_Intercept - b_stress_sum + b_syllable_sum - `b_stress_sum:syllable_sum`,
    ss_paroxitone_cvc  = b_Intercept + b_stress_sum - b_syllable_sum - `b_stress_sum:syllable_sum`,
    ss_oxitone_cvc     = b_Intercept - b_stress_sum - b_syllable_sum + `b_stress_sum:syllable_sum`,

    la_paroxitone_cv   = b_Intercept + b_stress_sum + b_syllable_sum + `b_stress_sum:syllable_sum` + b_groupla + `b_groupla:stress_sum` + `b_groupla:syllable_sum` + `b_groupla:stress_sum:syllable_sum`,
    la_oxitone_cv      = b_Intercept - b_stress_sum + b_syllable_sum - `b_stress_sum:syllable_sum` + b_groupla - `b_groupla:stress_sum` + `b_groupla:syllable_sum` - `b_groupla:stress_sum:syllable_sum`,
    la_paroxitone_cvc  = b_Intercept + b_stress_sum - b_syllable_sum - `b_stress_sum:syllable_sum` + b_groupla + `b_groupla:stress_sum` - `b_groupla:syllable_sum` - `b_groupla:stress_sum:syllable_sum`,
    la_oxitone_cvc     = b_Intercept - b_stress_sum - b_syllable_sum + `b_stress_sum:syllable_sum` + b_groupla - `b_groupla:stress_sum` - `b_groupla:syllable_sum` + `b_groupla:stress_sum:syllable_sum`,

    lb_paroxitone_cv   = b_Intercept + b_stress_sum + b_syllable_sum + `b_stress_sum:syllable_sum` + b_grouplb + `b_grouplb:stress_sum` + `b_grouplb:syllable_sum` + `b_grouplb:stress_sum:syllable_sum`,
    lb_oxitone_cv      = b_Intercept - b_stress_sum + b_syllable_sum - `b_stress_sum:syllable_sum` + b_grouplb - `b_grouplb:stress_sum` + `b_grouplb:syllable_sum` - `b_grouplb:stress_sum:syllable_sum`,
    lb_paroxitone_cvc  = b_Intercept + b_stress_sum - b_syllable_sum - `b_stress_sum:syllable_sum` + b_grouplb + `b_grouplb:stress_sum` - `b_grouplb:syllable_sum` - `b_grouplb:stress_sum:syllable_sum`,
    lb_oxitone_cvc     = b_Intercept - b_stress_sum - b_syllable_sum + `b_stress_sum:syllable_sum` + b_grouplb - `b_grouplb:stress_sum` - `b_grouplb:syllable_sum` + `b_grouplb:stress_sum:syllable_sum`,

    hs_paroxitone_cv   = b_Intercept + b_stress_sum + b_syllable_sum + `b_stress_sum:syllable_sum` + b_grouphs + `b_grouphs:stress_sum` + `b_grouphs:syllable_sum` + `b_grouphs:stress_sum:syllable_sum`,
    hs_oxitone_cv      = b_Intercept - b_stress_sum + b_syllable_sum - `b_stress_sum:syllable_sum` + b_grouphs - `b_grouphs:stress_sum` + `b_grouphs:syllable_sum` - `b_grouphs:stress_sum:syllable_sum`,
    hs_paroxitone_cvc  = b_Intercept + b_stress_sum - b_syllable_sum - `b_stress_sum:syllable_sum` + b_grouphs + `b_grouphs:stress_sum` - `b_grouphs:syllable_sum` - `b_grouphs:stress_sum:syllable_sum`,
    hs_oxitone_cvc     = b_Intercept - b_stress_sum - b_syllable_sum + `b_stress_sum:syllable_sum` + b_grouphs - `b_grouphs:stress_sum` - `b_grouphs:syllable_sum` + `b_grouphs:stress_sum:syllable_sum`,

    int_paroxitone_cv  = b_Intercept + b_stress_sum + b_syllable_sum + `b_stress_sum:syllable_sum` + b_groupint + `b_groupint:stress_sum` + `b_groupint:syllable_sum` + `b_groupint:stress_sum:syllable_sum`,
    int_oxitone_cv     = b_Intercept - b_stress_sum + b_syllable_sum - `b_stress_sum:syllable_sum` + b_groupint - `b_groupint:stress_sum` + `b_groupint:syllable_sum` - `b_groupint:stress_sum:syllable_sum`,
    int_paroxitone_cvc = b_Intercept + b_stress_sum - b_syllable_sum - `b_stress_sum:syllable_sum` + b_groupint + `b_groupint:stress_sum` - `b_groupint:syllable_sum` - `b_groupint:stress_sum:syllable_sum`,
    int_oxitone_cvc    = b_Intercept - b_stress_sum - b_syllable_sum + `b_stress_sum:syllable_sum` + b_groupint - `b_groupint:stress_sum` - `b_groupint:syllable_sum` + `b_groupint:stress_sum:syllable_sum`
  ) %>%
  gather(group, estimate) %>%
  separate(group, into = c("group", "stress", "syllable")) %>%
  mutate(landmark = "ts_suffix_start_posterior",
         landmark_labs = "V2\n(suffix)")


ts_next_word_posterior <- posterior_samples(mod_ts_next_word_full) %>%
  select(., starts_with("b_")) %>%
  transmute(
    ss_paroxitone_cv   = b_Intercept + b_stress_sum + b_syllable_sum + `b_stress_sum:syllable_sum`,
    ss_oxitone_cv      = b_Intercept - b_stress_sum + b_syllable_sum - `b_stress_sum:syllable_sum`,
    ss_paroxitone_cvc  = b_Intercept + b_stress_sum - b_syllable_sum - `b_stress_sum:syllable_sum`,
    ss_oxitone_cvc     = b_Intercept - b_stress_sum - b_syllable_sum + `b_stress_sum:syllable_sum`,

    la_paroxitone_cv   = b_Intercept + b_stress_sum + b_syllable_sum + `b_stress_sum:syllable_sum` + b_groupla + `b_groupla:stress_sum` + `b_groupla:syllable_sum` + `b_groupla:stress_sum:syllable_sum`,
    la_oxitone_cv      = b_Intercept - b_stress_sum + b_syllable_sum - `b_stress_sum:syllable_sum` + b_groupla - `b_groupla:stress_sum` + `b_groupla:syllable_sum` - `b_groupla:stress_sum:syllable_sum`,
    la_paroxitone_cvc  = b_Intercept + b_stress_sum - b_syllable_sum - `b_stress_sum:syllable_sum` + b_groupla + `b_groupla:stress_sum` - `b_groupla:syllable_sum` - `b_groupla:stress_sum:syllable_sum`,
    la_oxitone_cvc     = b_Intercept - b_stress_sum - b_syllable_sum + `b_stress_sum:syllable_sum` + b_groupla - `b_groupla:stress_sum` - `b_groupla:syllable_sum` + `b_groupla:stress_sum:syllable_sum`,

    lb_paroxitone_cv   = b_Intercept + b_stress_sum + b_syllable_sum + `b_stress_sum:syllable_sum` + b_grouplb + `b_grouplb:stress_sum` + `b_grouplb:syllable_sum` + `b_grouplb:stress_sum:syllable_sum`,
    lb_oxitone_cv      = b_Intercept - b_stress_sum + b_syllable_sum - `b_stress_sum:syllable_sum` + b_grouplb - `b_grouplb:stress_sum` + `b_grouplb:syllable_sum` - `b_grouplb:stress_sum:syllable_sum`,
    lb_paroxitone_cvc  = b_Intercept + b_stress_sum - b_syllable_sum - `b_stress_sum:syllable_sum` + b_grouplb + `b_grouplb:stress_sum` - `b_grouplb:syllable_sum` - `b_grouplb:stress_sum:syllable_sum`,
    lb_oxitone_cvc     = b_Intercept - b_stress_sum - b_syllable_sum + `b_stress_sum:syllable_sum` + b_grouplb - `b_grouplb:stress_sum` - `b_grouplb:syllable_sum` + `b_grouplb:stress_sum:syllable_sum`,

    hs_paroxitone_cv   = b_Intercept + b_stress_sum + b_syllable_sum + `b_stress_sum:syllable_sum` + b_grouphs + `b_grouphs:stress_sum` + `b_grouphs:syllable_sum` + `b_grouphs:stress_sum:syllable_sum`,
    hs_oxitone_cv      = b_Intercept - b_stress_sum + b_syllable_sum - `b_stress_sum:syllable_sum` + b_grouphs - `b_grouphs:stress_sum` + `b_grouphs:syllable_sum` - `b_grouphs:stress_sum:syllable_sum`,
    hs_paroxitone_cvc  = b_Intercept + b_stress_sum - b_syllable_sum - `b_stress_sum:syllable_sum` + b_grouphs + `b_grouphs:stress_sum` - `b_grouphs:syllable_sum` - `b_grouphs:stress_sum:syllable_sum`,
    hs_oxitone_cvc     = b_Intercept - b_stress_sum - b_syllable_sum + `b_stress_sum:syllable_sum` + b_grouphs - `b_grouphs:stress_sum` - `b_grouphs:syllable_sum` + `b_grouphs:stress_sum:syllable_sum`,

    int_paroxitone_cv  = b_Intercept + b_stress_sum + b_syllable_sum + `b_stress_sum:syllable_sum` + b_groupint + `b_groupint:stress_sum` + `b_groupint:syllable_sum` + `b_groupint:stress_sum:syllable_sum`,
    int_oxitone_cv     = b_Intercept - b_stress_sum + b_syllable_sum - `b_stress_sum:syllable_sum` + b_groupint - `b_groupint:stress_sum` + `b_groupint:syllable_sum` - `b_groupint:stress_sum:syllable_sum`,
    int_paroxitone_cvc = b_Intercept + b_stress_sum - b_syllable_sum - `b_stress_sum:syllable_sum` + b_groupint + `b_groupint:stress_sum` - `b_groupint:syllable_sum` - `b_groupint:stress_sum:syllable_sum`,
    int_oxitone_cvc    = b_Intercept - b_stress_sum - b_syllable_sum + `b_stress_sum:syllable_sum` + b_groupint - `b_groupint:stress_sum` - `b_groupint:syllable_sum` + `b_groupint:stress_sum:syllable_sum`
  ) %>%
  gather(group, estimate) %>%
  separate(group, into = c("group", "stress", "syllable")) %>%
  mutate(landmark = "ts_next_word_posterior",
         landmark_labs = "Following\nword")


complete_posteriors <- bind_rows(
  tw_start_posterior,
  ts_v1_start_posterior,
  ts_v1_20_posterior,
  ts_syl1_end_posterior,
  ts_suffix_start_posterior,
  ts_next_word_posterior
)







complete_posteriors %>%
  ggplot(., aes(x = landmark, y = estimate, color = group)) +
    facet_grid(syllable ~ stress) +
    geom_hline(yintercept = 0, lty = 3) +
    stat_pointinterval(position = position_dodge(0.5)) +
    landmark_posterior_theme(base_size = 12)

