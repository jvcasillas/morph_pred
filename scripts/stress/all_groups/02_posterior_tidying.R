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




posterior_samples(mod_tw_start_full) %>%
  select(., b_Intercept:b_syllable_sum) %>%
  transmute(
    ss_cv_oxitone = b_Intercept - b_stress_sum + b_syllable_sum,
    ss_cv_paroxitone = b_Intercept + b_stress_sum + b_syllable_sum,
    ss_cvc_oxitone = b_Intercept - b_stress_sum - b_syllable_sum,
    ss_cvc_paroxitone = b_Intercept + b_stress_sum - b_syllable_sum,
    la_cv_oxitone = b_Intercept + b_groupla - b_stress_sum + b_syllable_sum,
    la_cv_paroxitone = b_Intercept + b_groupla + b_stress_sum + b_syllable_sum,
    la_cvc_oxitone = b_Intercept + b_groupla - b_stress_sum - b_syllable_sum,
    la_cvc_paroxitone =b_Intercept + b_groupla + b_stress_sum - b_syllable_sum,
    lb_cv_oxitone = b_Intercept + b_grouplb - b_stress_sum + b_syllable_sum,
    lb_cv_paroxitone = b_Intercept + b_grouplb + b_stress_sum + b_syllable_sum,
    lb_cvc_oxitone = b_Intercept + b_grouplb - b_stress_sum - b_syllable_sum,
    lb_cvc_paroxitone = b_Intercept + b_grouplb + b_stress_sum - b_syllable_sum,
    hs_cv_oxitone = b_Intercept + b_grouphs - b_stress_sum + b_syllable_sum,
    hs_cv_paroxitone = b_Intercept + b_grouphs + b_stress_sum + b_syllable_sum,
    hs_cvc_oxitone = b_Intercept + b_grouphs - b_stress_sum - b_syllable_sum,
    hs_cvc_paroxitone = b_Intercept + b_grouphs + b_stress_sum - b_syllable_sum,
    int_cv_oxitone = b_Intercept + b_groupint - b_stress_sum + b_syllable_sum,
    int_cv_paroxitone = b_Intercept + b_groupint + b_stress_sum + b_syllable_sum,
    int_cvc_oxitone = b_Intercept + b_groupint - b_stress_sum - b_syllable_sum,
    int_cvc_paroxitone = b_Intercept + b_groupint + b_stress_sum - b_syllable_sum) %>%
  gather(group, estimate) %>%
  separate(group, into = c("group", "syllable", "stress")) %>%
  mutate(landmark = "start") %>%
  ggplot(., aes(x = landmark, y = estimate, color = group)) +
    facet_grid(syllable ~ stress) +
    geom_hline(yintercept = 0, lty = 3) +
    stat_pointinterval(position = position_dodge(0.5)) +
    landmark_posterior_theme(base_size = 12)








