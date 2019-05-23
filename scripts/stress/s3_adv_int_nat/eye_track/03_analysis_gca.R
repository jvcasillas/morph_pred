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
list2env(all_rds, globalenv())

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
            coda_sum = if_else(coda == 0, 1, -1)) %>%
  poly_add_columns(., time_zero, degree = 3, prefix = "ot")

# -----------------------------------------------------------------------------
















# Random effects structure ----------------------------------------------------

# Build up random effects to test time terms
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

mod_ot4 <- update(mod_ot3, . ~ . + (1 | target))

anova(mod_ot1, mod_ot2, mod_ot3, mod_ot4)

}

# -----------------------------------------------------------------------------





# Individual models -----------------------------------------------------------

#
# only ss
#

gca_mod_ss_base <-
  lmer(eLog ~ 1 + (ot1 + ot2 + ot3) +
         (1 + coda_sum + condition_sum + (ot1 + ot2 + ot3) | participant) +
         (1 + ot1 + ot2 + ot3 | target),
       control = lmerControl(optimizer = 'bobyqa'), REML = F,
       data = filter(stress_gc_subset, group == "ss"))

# add coda effect to intercept
# add coda effect to linear slope
# add coda effect to quadratic time term
# add coda effect to cubic time term
gca_mod_ss_coda_0 <- update(gca_mod_ss_base,   . ~ . + coda_sum)
gca_mod_ss_coda_1 <- update(gca_mod_ss_coda_0, . ~ . + ot1:coda_sum)
gca_mod_ss_coda_2 <- update(gca_mod_ss_coda_1, . ~ . + ot2:coda_sum)
gca_mod_ss_coda_3 <- update(gca_mod_ss_coda_2, . ~ . + ot3:coda_sum)

anova(gca_mod_ss_base, gca_mod_ss_coda_0, gca_mod_ss_coda_1,
      gca_mod_ss_coda_2, gca_mod_ss_coda_3)

# add condition effect to intercept
# add condition effect to linear slope
# add condition effect to quadratic time term
# add condition effect to cubic time term
gca_mod_ss_cond_0 <- update(gca_mod_ss_coda_3,   . ~ . + condition_sum)
gca_mod_ss_cond_1 <- update(gca_mod_ss_cond_0,   . ~ . + ot1:condition_sum)
gca_mod_ss_cond_2 <- update(gca_mod_ss_cond_1,   . ~ . + ot2:condition_sum)
gca_mod_ss_cond_3 <- update(gca_mod_ss_cond_2,   . ~ . + ot3:condition_sum)

anova(gca_mod_ss_coda_3, gca_mod_ss_cond_0, gca_mod_ss_cond_1,
      gca_mod_ss_cond_2, gca_mod_ss_cond_3)

# add coda x cond int to intercept
# add coda x cond int to linear slope
# add coda x cond int to quadratic time term
# add coda x cond int to cubic time term
gca_mod_ss_int_0 <- update(gca_mod_ss_cond_3, . ~ . + coda_sum:condition_sum)
gca_mod_ss_int_1 <- update(gca_mod_ss_int_0,  . ~ . + ot1:coda_sum:condition_sum)
gca_mod_ss_int_2 <- update(gca_mod_ss_int_1,  . ~ . + ot2:coda_sum:condition_sum)
gca_mod_ss_int_3 <- update(gca_mod_ss_int_2,  . ~ . + ot3:coda_sum:condition_sum)

anova(gca_mod_ss_cond_3, gca_mod_ss_int_0, gca_mod_ss_int_1,
      gca_mod_ss_int_2, gca_mod_ss_int_3)



#
# only la
#

gca_mod_la_base <-
  lmer(eLog ~ 1 + (ot1 + ot2 + ot3) +
         (1 + coda_sum + condition_sum + (ot1 + ot2 + ot3) | participant) +
         (1 + ot1 + ot2 + ot3 | target),
       control = lmerControl(optimizer = 'bobyqa'), REML = F,
       data = filter(stress_gc_subset, group == "la"))

# add coda effect to intercept
# add coda effect to linear slope
# add coda effect to quadratic time term
# add coda effect to cubic time term
gca_mod_la_coda_0 <- update(gca_mod_la_base,   . ~ . + coda_sum)
gca_mod_la_coda_1 <- update(gca_mod_la_coda_0, . ~ . + ot1:coda_sum)
gca_mod_la_coda_2 <- update(gca_mod_la_coda_1, . ~ . + ot2:coda_sum)
gca_mod_la_coda_3 <- update(gca_mod_la_coda_2, . ~ . + ot3:coda_sum)

anova(gca_mod_la_base, gca_mod_la_coda_0, gca_mod_la_coda_1,
      gca_mod_la_coda_2, gca_mod_la_coda_3)

# add condition effect to intercept
# add condition effect to linear slope
# add condition effect to quadratic time term
# add condition effect to cubic time term
gca_mod_la_cond_0 <- update(gca_mod_la_coda_3,   . ~ . + condition_sum)
gca_mod_la_cond_1 <- update(gca_mod_la_cond_0,   . ~ . + ot1:condition_sum)
gca_mod_la_cond_2 <- update(gca_mod_la_cond_1,   . ~ . + ot2:condition_sum)
gca_mod_la_cond_3 <- update(gca_mod_la_cond_2,   . ~ . + ot3:condition_sum)

anova(gca_mod_la_coda_3, gca_mod_la_cond_0, gca_mod_la_cond_1,
      gca_mod_la_cond_2, gca_mod_la_cond_3)

# add coda x cond int to intercept
# add coda x cond int to linear slope
# add coda x cond int to quadratic time term
# add coda x cond int to cubic time term
gca_mod_la_int_0 <- update(gca_mod_la_cond_3, . ~ . + coda_sum:condition_sum)
gca_mod_la_int_1 <- update(gca_mod_la_int_0,  . ~ . + ot1:coda_sum:condition_sum)
gca_mod_la_int_2 <- update(gca_mod_la_int_1,  . ~ . + ot2:coda_sum:condition_sum)
gca_mod_la_int_3 <- update(gca_mod_la_int_2,  . ~ . + ot3:coda_sum:condition_sum)

anova(gca_mod_la_cond_3, gca_mod_la_int_0, gca_mod_la_int_1,
      gca_mod_la_int_2, gca_mod_la_int_3)


#
# only int
#

gca_mod_int_base <-
  lmer(eLog ~ 1 + (ot1 + ot2 + ot3) +
         (1 + coda_sum + condition_sum + (ot1 + ot2 + ot3) | participant) +
         (1 + ot1 + ot2 + ot3 | target),
       control = lmerControl(optimizer = 'bobyqa'), REML = F,
       data = filter(stress_gc_subset, group == "int"))

# add coda effect to intercept
# add coda effect to linear slope
# add coda effect to quadratic time term
# add coda effect to cubic time term
gca_mod_int_coda_0 <- update(gca_mod_int_base,   . ~ . + coda_sum)
gca_mod_int_coda_1 <- update(gca_mod_int_coda_0, . ~ . + ot1:coda_sum)
gca_mod_int_coda_2 <- update(gca_mod_int_coda_1, . ~ . + ot2:coda_sum)
gca_mod_int_coda_3 <- update(gca_mod_int_coda_2, . ~ . + ot3:coda_sum)

anova(gca_mod_int_base, gca_mod_int_coda_0, gca_mod_int_coda_1,
      gca_mod_int_coda_2, gca_mod_int_coda_3)

# add condition effect to intercept
# add condition effect to linear slope
# add condition effect to quadratic time term
# add condition effect to cubic time term
gca_mod_int_cond_0 <- update(gca_mod_int_coda_3,   . ~ . + condition_sum)
gca_mod_int_cond_1 <- update(gca_mod_int_cond_0,   . ~ . + ot1:condition_sum)
gca_mod_int_cond_2 <- update(gca_mod_int_cond_1,   . ~ . + ot2:condition_sum)
gca_mod_int_cond_3 <- update(gca_mod_int_cond_2,   . ~ . + ot3:condition_sum)

anova(gca_mod_int_coda_3, gca_mod_int_cond_0, gca_mod_int_cond_1,
      gca_mod_int_cond_2, gca_mod_int_cond_3)

# add coda x cond int to intercept
# add coda x cond int to linear slope
# add coda x cond int to quadratic time term
# add coda x cond int to cubic time term
gca_mod_int_int_0 <- update(gca_mod_int_cond_3, . ~ . + coda_sum:condition_sum)
gca_mod_int_int_1 <- update(gca_mod_int_int_0,  . ~ . + ot1:coda_sum:condition_sum)
gca_mod_int_int_2 <- update(gca_mod_int_int_1,  . ~ . + ot2:coda_sum:condition_sum)
gca_mod_int_int_3 <- update(gca_mod_int_int_2,  . ~ . + ot3:coda_sum:condition_sum)

anova(gca_mod_int_cond_3, gca_mod_int_int_0, gca_mod_int_int_1,
      gca_mod_int_int_2, gca_mod_int_int_3)

# -----------------------------------------------------------------------------








# Full model ------------------------------------------------------------------

if(F){
# Base model
gca_full_mod_base <-
  lmer(eLog ~ 1 + (ot1 + ot2 + ot3) * coda_sum * condition_sum +
         (1 + coda_sum + condition_sum + (ot1 + ot2 + ot3) | participant) +
         (1 + ot1 + ot2 + ot3 | target),
       control = lmerControl(optimizer = 'bobyqa',
                             optCtrl = list(maxfun = 2e4)),
       data = stress_gc_subset, REML = F)

# Add group effect on intercept
# Add group effect on linear slope
# Add group effect on quadratic poly
# Add group effect of cubic poly
gca_full_mod_group_0 <- update(gca_full_mod_base,    . ~ . + group)
gca_full_mod_group_1 <- update(gca_full_mod_group_0, . ~ . + ot1:group)
gca_full_mod_group_2 <- update(gca_full_mod_group_1, . ~ . + ot2:group)
gca_full_mod_group_3 <- update(gca_full_mod_group_2, . ~ . + ot3:group)

anova(gca_full_mod_base,
      gca_full_mod_group_0,
      gca_full_mod_group_1,
      gca_full_mod_group_2,
      gca_full_mod_group_3)

# Add (group x coda_sum x condition_sum) interaction to intercept
# Add (group x coda_sum x condition_sum) interaction to linear slope
# Add (group x coda_sum x condition_sum) interaction to quadratic time term
# Add (group x coda_sum x condition_sum) interaction to cubic time term
gca_full_mod_int_0 <- update(gca_full_mod_group_3, . ~ . + coda_sum:condition_sum:group)
gca_full_mod_int_1 <- update(gca_full_mod_int_0,   . ~ . + ot1:coda_sum:condition_sum:group)
gca_full_mod_int_2 <- update(gca_full_mod_int_1,   . ~ . + ot2:coda_sum:condition_sum:group)
gca_full_mod_int_3 <- update(gca_full_mod_int_2,   . ~ . + ot3:coda_sum:condition_sum:group)

anova(gca_full_mod_group_3,
      gca_full_mod_int_0,
      gca_full_mod_int_1,
      gca_full_mod_int_2,
      gca_full_mod_int_3)

}

# -----------------------------------------------------------------------------







# Comparison model ------------------------------------------------------------

# Remove SS
gca_int_la <- stress_gc_subset %>%
  filter(group != "ss") %>%
  mutate(group_sum = if_else(group == "int", 1, -1))

gca_mod_int_la <-
  lmer(eLog ~ 1 + (ot1 + ot2 + ot3) + group_sum + coda_sum + condition_sum +
         ot1:group_sum + ot2:group_sum + ot3:group_sum +
         ot1:coda_sum + ot2:coda_sum + ot3:coda_sum +
         ot1:condition_sum + ot2:condition_sum + ot3:condition_sum +
         ot1:group_sum:coda_sum + ot2:group_sum:coda_sum + ot3:group_sum:coda_sum +
         ot1:group_sum:condition_sum + ot2:group_sum:condition_sum + ot3:group_sum:condition_sum +
         (1 + coda_sum + condition_sum + (ot1 + ot2 + ot3) | participant),
       control = lmerControl(optimizer = 'bobyqa'),
       data = gca_int_la, weights = 1/wts, REML = F)

# -----------------------------------------------------------------------------








# Save models -----------------------------------------------------------------

if(F) {
# Build model names programatically
mod_type <- "gca_mod_"
mod_spec <- c("_base", "_coda_0", "_coda_1", "_coda_2", "_coda_3", "_cond_0",
              "_cond_1", "_cond_2", "_cond_3", "_int_0", "_int_1", "_int_2",
              "_int_3")

# Store ind models in list
ind_mods <- mget(c(paste0(mod_type, "ss", mod_spec),
                   paste0(mod_type, "la", mod_spec),
                   paste0(mod_type, "int", mod_spec)))

# Save models as .rds iterating over list
ind_mods %>%
  names(.) %>%
  map(~ saveRDS(ind_mods[[.]], compress = 'xz',
                file = paste0(here("models", "stress", "s3_adv_int_nat",
                                   "eye_track", "gca"), "/", ., ".rds")))


# Store full (ot1, ot2, ot3, group, coda, cond) models in list
full_mods <- mget(c(
  "gca_full_mod_base", "gca_full_mod_group_0", "gca_full_mod_group_1",
  "gca_full_mod_group_2", "gca_full_mod_group_3", "gca_full_mod_int_0",
  "gca_full_mod_int_1", "gca_full_mod_int_2", "gca_full_mod_int_3"))

# Save models as .rds iterating over list
full_mods %>%
  names(.) %>%
  map(~ saveRDS(full_mods[[.]], compress = 'xz',
                file = paste0(here("models", "stress", "s3_adv_int_nat",
                                   "eye_track", "gca"), "/", ., ".rds")))

}

# -----------------------------------------------------------------------------

