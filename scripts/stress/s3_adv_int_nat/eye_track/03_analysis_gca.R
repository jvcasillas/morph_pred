# Growth curve analyisis ------------------------------------------------------
#
# - Question 1: Are the groups different from each other in when they begin
#   to fixate on the target?
#     - test 3 groups at each level of 'condition'
#     - hypothesis: SS has steeper slope for both conditions
# - Question 2: W/in groups, is the there a difference between
#   oxytone/paroxytone items?
#     - test oxytone vs. paroxytone for each group
#     - hypothesis: steeper slope/earlier break in oxytone condition
#
# -----------------------------------------------------------------------------





# Load data and models --------------------------------------------------------

source(here::here("scripts", "01_load_data.R"))

gca_mods_path  <- here("models", "stress", "s3_adv_int_nat", "eye_track", "gca")
gca_mods       <- list.files(gca_mods_path, pattern = ".rds")
all_gca_mods   <- paste0(gca_mods_path, "/", gca_mods)
all_rds        <- lapply(all_gca_mods, readRDS)
names(all_rds) <- gsub(".rds", "",
                       list.files(gca_mods_path, pattern = ".rds",
                                  full.names = FALSE), fixed = TRUE)

# -----------------------------------------------------------------------------







# Data prep -------------------------------------------------------------------

# - subset using time course
#    - We need to reduce the time course to a relevant time window that
#      that includes enough of the trajectory from before and after the
#      target syllable onset
#    - Importantly, we need to make sure that the adjusted time course
#      is centered at 200ms after the offset of the first syllable
#    - This is because the orthogonal polynomials center the time course,
#      thus the parameter estimates on the intercept and the linear slope
#      are calculated for the midpoint (0).
#    - This has an added bonus of assessing group differences at the mid
#      point (200ms after target syllable offset), which will corroborate
#      the results from the GLMMs.
#    - We can select the appropriate time course subset by selecting the
#      target syllable offset, bin 4 (200ms / 50 = 4), and keeping an
#      equal number of bins on each side:
#                     8 7 6 5 4 3 2 1 X 1 2 3 4 5 6 7 8
#                                     ^
#                     center of time course (bin 4)
#
#
# Number of bins:     1  2  3  4 5 6 7 8 9 10 11 12 13 14 15 16 17
# Actual bin number: -4 -3 -2 -1 0 1 2 3 4  5  6  7  8  9 10 11 12

stress_gc_subset <- stress50 %>%
  filter(., group %in% c('int', 'la', 'ss'),
           !participant %in% c("L01", "L02", "L03", "L04", "L05",
                               "L06", "L07", "L08", "L09", "L10",
                               "L15", "L20", "L21", "L22", "L23",
                               "L26", "L30", "L31", "L33", "LA04",
                               "LA06", "LA09", "LA14", "LA15", "LA19"),
            time_zero >= -4 & time_zero <= 12) %>%
  mutate(., group = fct_relevel(group, "ss", "la", "int"),
            condition_sum = if_else(condition == "stressed", 1, -1),
            coda_sum = if_else(coda == 1, 1, -1)) %>%
  poly_add_columns(., time_zero, degree = 3, prefix = "ot")

# -----------------------------------------------------------------------------
















# Random effects structure ----------------------------------------------------

if(F){
mod_ot1 <-
  lmer(eLog ~ 1 + ot1 +
      (1 + coda_sum + condition_sum + ot1 | participant),
       control = lmerControl(optimizer = 'bobyqa'),
       data = stress_gc_subset, weights = 1/wts, REML = F)

mod_ot2 <-
  update(mod_ot1, . ~ . -(1 + coda_sum + condition_sum + ot1 | participant) +
         ot2 + (1 + coda_sum + condition_sum + ot1 + ot2 | participant))

mod_ot3 <-
  update(mod_ot2, . ~ . -(1 + coda_sum + condition_sum + ot1 + ot2 | participant) +
         ot3 + (1 + coda_sum + condition_sum + ot1 + ot2 + ot3 | participant))

anova(mod_ot1,
      mod_ot2,
      mod_ot3,
      test = "Chisq")
}

# -----------------------------------------------------------------------------







# Test fixed effects ----------------------------------------------------------

# Base model
if(F){
gca_mod_base <-
  lmer(eLog ~ 1 + (ot1 + ot2 + ot3) +
      (1 + coda_sum + condition_sum + (ot1 + ot2 + ot3) | participant),
       control = lmerControl(optimizer = 'bobyqa'),
       data = stress_gc_subset, weights = 1/wts, REML = F)

# Add group effect on intercept
# Add group effect on linear slope
# Add group effect on quadratic poly
# Add group effect of cubic poly
gca_mod_group_0 <- update(gca_mod_base,    . ~ . + group)
gca_mod_group_1 <- update(gca_mod_group_0, . ~ . + ot1:group)
gca_mod_group_2 <- update(gca_mod_group_1, . ~ . + ot2:group)
gca_mod_group_3 <- update(gca_mod_group_2, . ~ . + ot3:group)

anova(gca_mod_base,
      gca_mod_group_0,
      gca_mod_group_1,
      gca_mod_group_2,
      gca_mod_group_3,
      test = 'Chisq')

# Add coda_sum to intercept
# Add coda_sum to linear slope
# Add coda_sum to quadratic effect
# Add coda_sum to cubic poly
gca_mod_coda_0 <- update(gca_mod_group_3, . ~ . + coda_sum)
gca_mod_coda_1 <- update(gca_mod_coda_0, . ~ . + ot1:coda_sum)
gca_mod_coda_2 <- update(gca_mod_coda_1, . ~ . + ot2:coda_sum)
gca_mod_coda_3 <- update(gca_mod_coda_2, . ~ . + ot3:coda_sum)

anova(gca_mod_group_3,
      gca_mod_coda_0,
      gca_mod_coda_1,
      gca_mod_coda_2,
      gca_mod_coda_3,
      test = "Chisq")

# Add 3-way interactions (time x group x coda_sum)
gca_mod_int_1 <- update(gca_mod_coda_3, . ~ . + ot1:group:coda_sum)
gca_mod_int_2 <- update(gca_mod_int_1,  . ~ . + ot2:group:coda_sum)
gca_mod_int_3 <- update(gca_mod_int_2,  . ~ . + ot3:group:coda_sum)

anova(gca_mod_coda_3,
      gca_mod_int_1,
      gca_mod_int_2,
      gca_mod_int_3,
      test = "Chisq")

# Add condition_sum to intercept
# Add condition_sum to linear slope
# Add condition_sum to quadratic poly
# Add condition_sum to cubic poly
gca_mod_cond_0 <- update(gca_mod_int_3,  . ~ . + condition_sum)
gca_mod_cond_1 <- update(gca_mod_cond_0, . ~ . + ot1:condition_sum)
gca_mod_cond_2 <- update(gca_mod_cond_1, . ~ . + ot2:condition_sum)
gca_mod_cond_3 <- update(gca_mod_cond_2, . ~ . + ot3:condition_sum)

anova(gca_mod_int_3,
      gca_mod_cond_0,
      gca_mod_cond_1,
      gca_mod_cond_2,
      gca_mod_cond_3,
      test = "Chisq")

# Add 3-way interactions (time x group x condition_sum)
gca_mod_int_4 <- update(gca_mod_cond_3, . ~ . + ot1:group:condition_sum)
gca_mod_int_5 <- update(gca_mod_int_4,  . ~ . + ot2:group:condition_sum)
gca_mod_int_6 <- update(gca_mod_int_5,  . ~ . + ot3:group:condition_sum)

anova(gca_mod_cond_3,
      gca_mod_int_4,
      gca_mod_int_5,
      gca_mod_int_6,
      test = "Chisq")

}

# -----------------------------------------------------------------------------








# Save models -----------------------------------------------------------------

# Base model
saveRDS(gca_mod_base, compress = 'xz',
        file = here("models", "stress", "s3_adv_int_nat", "eye_track",
                    "gca", "gca_mod_base.rds"))

# Group effects
saveRDS(gca_mod_group_0, compress = 'xz',
        file = here("models", "stress", "s3_adv_int_nat", "eye_track",
                    "gca", "gca_mod_group_0.rds"))
saveRDS(gca_mod_group_1, compress = 'xz',
        file = here("models", "stress", "s3_adv_int_nat", "eye_track",
                    "gca", "gca_mod_group_1.rds"))
saveRDS(gca_mod_group_2, compress = 'xz',
        file = here("models", "stress", "s3_adv_int_nat", "eye_track",
                    "gca", "gca_mod_group_2.rds"))
saveRDS(gca_mod_group_3, compress = 'xz',
        file = here("models", "stress", "s3_adv_int_nat", "eye_track",
                    "gca", "gca_mod_group_3.rds"))

# Coda effects
saveRDS(gca_mod_coda_0, compress = 'xz',
        file = here("models", "stress", "s3_adv_int_nat", "eye_track",
                    "gca", "gca_mod_coda_0.rds"))
saveRDS(gca_mod_coda_1, compress = 'xz',
        file = here("models", "stress", "s3_adv_int_nat", "eye_track",
                    "gca", "gca_mod_coda_1.rds"))
saveRDS(gca_mod_coda_2, compress = 'xz',
        file = here("models", "stress", "s3_adv_int_nat", "eye_track",
                    "gca", "gca_mod_coda_2.rds"))
saveRDS(gca_mod_coda_3, compress = 'xz',
        file = here("models", "stress", "s3_adv_int_nat", "eye_track",
                    "gca", "gca_mod_coda_3.rds"))

# Time coda group interactions
saveRDS(gca_mod_int_1, compress = 'xz',
        file = here("models", "stress", "s3_adv_int_nat", "eye_track",
                    "gca", "gca_mod_int_1.rds"))
saveRDS(gca_mod_int_2, compress = 'xz',
        file = here("models", "stress", "s3_adv_int_nat", "eye_track",
                    "gca", "gca_mod_int_2.rds"))
saveRDS(gca_mod_int_3, compress = 'xz',
        file = here("models", "stress", "s3_adv_int_nat", "eye_track",
                    "gca", "gca_mod_int_3.rds"))

# Condition effects
saveRDS(gca_mod_cond_0, compress = 'xz',
        file = here("models", "stress", "s3_adv_int_nat", "eye_track",
                    "gca", "gca_mod_cond_0.rds"))
saveRDS(gca_mod_cond_1, compress = 'xz',
        file = here("models", "stress", "s3_adv_int_nat", "eye_track",
                    "gca", "gca_mod_cond_1.rds"))
saveRDS(gca_mod_cond_2, compress = 'xz',
        file = here("models", "stress", "s3_adv_int_nat", "eye_track",
                    "gca", "gca_mod_cond_2.rds"))
saveRDS(gca_mod_cond_3, compress = 'xz',
        file = here("models", "stress", "s3_adv_int_nat", "eye_track",
                    "gca", "gca_mod_cond_3.rds"))

# Time group condition interactions
saveRDS(gca_mod_int_4, compress = 'xz',
        file = here("models", "stress", "s3_adv_int_nat", "eye_track",
                    "gca", "gca_mod_int_4.rds"))
saveRDS(gca_mod_int_5, compress = 'xz',
        file = here("models", "stress", "s3_adv_int_nat", "eye_track",
                    "gca", "gca_mod_int_5.rds"))
saveRDS(gca_mod_int_6, compress = 'xz',
        file = here("models", "stress", "s3_adv_int_nat", "eye_track",
                    "gca", "gca_mod_int_6.rds"))

# -----------------------------------------------------------------------------








# Model descriptives ----------------------------------------------------------

# options(scipen = 999)
# summary(all_rds$gca_mod_int_6)
# summary(all_rds$gca_mod_int_3)
# options(scipen = 0)


# Random effects:
# Groups       Name        Variance  Std.Dev.   Corr
# participant (Intercept)    0.4366   0.6608
#              coda_sum      0.1333   0.3650    0.23
#              condition_sum 0.2950   0.5431   -0.13 -0.46
#              ot1           6.7137   2.5911    0.62  0.26  0.04
#              ot2           1.4992   1.2244    0.36 -0.22  0.24  0.60
#              ot3           0.6216   0.7884   -0.27 -0.01  0.01 -0.45  0.08
# Residual                   8.5308   2.9207
# Number of obs: 18981, groups:  participant, 71

# Fixed effects:
#                          Estimate  Std. Error          df t value        Pr(>|t|)
# (Intercept)               0.79921     0.13870    71.26867   5.762 0.0000001958262 ***
# ot1                       4.15427     0.54344    73.00429   7.644 0.0000000000654 ***
# ot2                      -0.60658     0.29032    70.31116  -2.089        0.040298 *
# ot3                      -0.94360     0.22801    65.07726  -4.138        0.000103 ***

# groupla                  -0.11868     0.19331    70.62651  -0.614        0.541213
# ot1:groupla              -0.12142     0.75197    70.91238  -0.161        0.872183
# ot2:groupla               1.36275     0.40570    70.20954   3.359        0.001267 **
# ot3:groupla               0.31861     0.32101    67.54678   0.993        0.324482

# groupint                 -0.28657     0.19933    70.14989  -1.438        0.154973
# ot1:groupint              0.04179     0.77470    70.20114   0.054        0.957137
# ot2:groupint              1.11186     0.41584    68.09693   2.674        0.009382 **
# ot3:groupint             -0.09549     0.32694    63.88240  -0.292        0.771189

# coda_sum                  0.08276     0.04481    66.25607   1.847        0.069250 .
# ot1:coda_sum             -0.67620     0.16385 18033.38549  -4.127 0.0000369334856 ***
# ot2:coda_sum             -0.27515     0.15752 16894.27825  -1.747        0.080708 .
# ot3:coda_sum              0.13886     0.15841 17422.92000   0.877        0.380728

# ot1:groupla:coda_sum      0.80938     0.23244 17898.18535   3.482        0.000499 ***
# ot2:groupla:coda_sum     -0.07077     0.22255 17262.48382  -0.318        0.750495
# ot3:groupla:coda_sum     -0.58210     0.22579 17562.20544  -2.578        0.009945 **

# ot1:groupint:coda_sum    -0.20717     0.23322 17868.00958  -0.888        0.374391
# ot2:groupint:coda_sum    -0.07305     0.22311 17359.80389  -0.327        0.743348
# ot3:groupint:coda_sum    -0.16848     0.22600 17295.77333  -0.745        0.455984

# -----------------------------------------------------------------------------
