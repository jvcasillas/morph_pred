# T-tests ---------------------------------------------------------------------
#
# - Question: can participants predict morphology before onset of target
#   suffix?
# - We test this for each group of participants for each type of word
#   (paroxytone, oxytone)
# - This analysis does not compare groups or conditions
#
# -----------------------------------------------------------------------------


# Load data -------------------------------------------------------------------

source(here::here("scripts", "01_load_data.R"))

# -----------------------------------------------------------------------------



# Create stress subset that incldues interpreters, late advanced and natives
# Must check on the participants that are removed
# We want to analyze proportion of target gaze at target onset
# so we need to make a subset of the data that only uses the
# target onset bin (adjusted 200ms for VWP)
# We are using 10ms bins so we want time_zero bin 20 (200/10 = 20)

df_short <- stress10 %>%
  filter(., group %in% c('int', 'la', 'ss'),
            !participant %in% c('L01', 'L02', 'L03', 'L04', 'L05',
                                'L06', 'L07', 'L08', 'L09', 'L10',
                                'L15', 'L20', 'L21', 'L22', 'L23',
                                'L26', 'L30', 'L31', 'L33', 'LA04',
                                'LA06', 'LA07', 'LA14'),
            time_zero == 20)



# Quick and dirty mean of target fixations as a function of
# group and condition (stressed, unstressed *1st syllable*)

df_short %>%
  na.omit(.) %>%
  group_by(., group, condition) %>%
  summarise(., meanFix = mean(targetProp))


# We will test this for each group in each condition (stressed, untressed)
# using a one-sided t-test. Specifically, we are testing the
# hypothesis that the proportion of looks is greater than
# chance (50%).
# - H0: u = 0.50
# - Ha: u > 0.50
# The generic code is:
# t.test(myVector, alternative = "greater", my = 0.33, conf.level = 0.95)

stress_ttest <- df_short %>%
  na.omit(.) %>%
  group_by(., group, condition, participant) %>%
  summarise(., meanFix = mean(targetProp)) %>%
  do(tidy(t.test(.$meanFix, alternative = "greater", mu = 0.5, conf.level = 0.99)))

# Convert pvalues from scientific notation
stress_ttest$p.value <- format(stress_ttest$p.value, scientific = F)
stress_ttest$sig <- "N.S."
stress_ttest[stress_ttest$p.value <= 0.05/6, 'sig'] <- "*"

# Print results
print(as.data.frame(stress_ttest[, c(1:7, 11)]))

saveRDS(stress_ttest, "./reports/stress/stress_int/mods/stress_ttest.rds", compress = "xz")

# group  condition  estimate  statistic       p.value parameter  conf.low  sig
#   int   stressed 0.6344643 3.46115554 0.00357388323         9 0.5248527    *
#   int unstressed 0.6753571 3.37992533 0.00406396949         9 0.5289754    *
#    la   stressed 0.5022928 0.05122216 0.47976999874        26 0.3913462 N.S.
#    la unstressed 0.5841931 2.12785440 0.02149710485        26 0.4861208 N.S.
#    ss   stressed 0.6322078 2.79896749 0.00537726129        21 0.5132880    *
#    ss unstressed 0.7235390 4.76754875 0.00005197742        21 0.6054925    *



