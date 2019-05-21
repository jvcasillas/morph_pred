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






# Data prep -------------------------------------------------------------------

# Create stress subset that incldues interpreters, late advanced and natives
# Must check on the participants that are removed
# We want to analyze proportion of target gaze at target onset
# so we need to make a subset of the data that only uses the
# target onset bin (adjusted 200ms for VWP)
# We are using 10ms bins so we want time_zero bin 20 (200/10 = 20)

df_short <- stress10 %>%
  filter(., group %in% c('int', 'la', 'ss'),
            !participant %in% c("L01", "L02", "L03", "L04", "L05",
                                "L06", "L07", "L08", "L09", "L10",
                                "L15", "L20", "L21", "L22", "L23",
                                "L26", "L30", "L31", "L33", "LA04",
                                "LA06", "LA09", "LA14", "LA15", "LA19"),
            time_zero == 20)

# Old vector of removed participants
# c('L01', 'L02', 'L03', 'L04', 'L05',
#   'L06', 'L07', 'L08', 'L09', 'L10',
#   'L15', 'L20', 'L21', 'L22', 'L23',
#   'L26', 'L30', 'L31', 'L33', 'LA04',
#   'LA06', 'LA07', 'LA14')

# New vector of participants to include
# c("LA01", "LA02", "LA03", "LA05", "LA07",
#   "LA08", "LA10", "LA11", "LA12", "LA13",
#   "LA16", "LA17", "LA18", "LA20", "LA21",
#   "LA22", "LA23", "LA24", "LA25", "LA26",
#   "LA27", "LA28", "LA29", "LA30", "LA31")

# -----------------------------------------------------------------------------






# T-tests ---------------------------------------------------------------------

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
  do(tidy(t.test(.$meanFix, alternative = "greater", mu = 0.5, conf.level = 0.95)))

# Convert pvalues from scientific notation
stress_ttest$p.value <- format(stress_ttest$p.value, scientific = F)
stress_ttest$sig <- "N.S."
stress_ttest[stress_ttest$p.value <= 0.05/6, 'sig'] <- "*"

# Print results
print(as.data.frame(stress_ttest[, c(1:7, 11)]))

saveRDS(stress_ttest, "./reports/stress/stress_int/mods/stress_ttest.rds", compress = "xz")

#   group  condition  estimate  statistic      p.value parameter  conf.low  sig
# 1   int   stressed 0.4892045 -0.2644177 0.6029820843        21 0.4189513 N.S.
# 2   int unstressed 0.5738636  1.8925669 0.0361411122        21 0.5067060 N.S.
# 3    la   stressed 0.5250000  0.6767972 0.2525030221        24 0.4618023 N.S.
# 4    la unstressed 0.5994286  2.1611355 0.0204408034        24 0.5207151 N.S.
# 5    ss   stressed 0.6334821  2.7498264 0.0057043367        23 0.5502873    *
# 6    ss unstressed 0.6763765  4.1900453 0.0001753115        23 0.6042325    *

# -----------------------------------------------------------------------------



