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

