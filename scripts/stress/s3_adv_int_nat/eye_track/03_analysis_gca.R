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





# Load data -------------------------------------------------------------------

source(here::here("scripts", "01_load_data.R"))

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
#                     6  5 4 3 2 1 X 1 2 3 4 5 6
#                                  ^
#                     center of time course (bin 4)
#
#
# Number of bins:     1  2 3 4 5 6 7 8 9 10 11 12 13
# Actual bin number: -2 -1 0 1 2 3 4 5 6  7  8  9 10

stress_gc_subset <- stress50 %>%
  filter(., group %in% c('int', 'la', 'ss'),
           !participant %in% c("L01", "L02", "L03", "L04", "L05",
                               "L06", "L07", "L08", "L09", "L10",
                               "L15", "L20", "L21", "L22", "L23",
                               "L26", "L30", "L31", "L33", "LA04",
                               "LA06", "LA09", "LA14", "LA15", "LA19"),
            time_zero >= -2 & time_zero <= 10) %>%
  mutate(., group = fct_relevel(group, "ss", "la", "int"),
            condition_sum = if_else(condition == "stressed", 1, -1),
            coda_sum = if_else(coda == 1, 1, -1)) %>%
  poly_add_columns(., time_zero, degree = 3, prefix = "ot")

# -----------------------------------------------------------------------------












####################################################
# - Question 1: Are the groups different from      #
#   each other in when they begin to fixate        #
#   on the target?                                 #
#     - test 3 groups at each level of 'condition' #
#     - hypothesis: SS has steeper slope for both  #
#       conditions                                 #
####################################################






## @knitr ignore2

# load the models:
gc_mod_base    <- readRDS('./mySources/models/stress/s3_adv_int_nat/eye_track/gc_mod_base.rds')
gc_mod_group_0 <- readRDS('./mySources/models/stress/s3_adv_int_nat/eye_track/gc_mod_group_0.rds')
gc_mod_group_1 <- readRDS('./mySources/models/stress/s3_adv_int_nat/eye_track/gc_mod_group_1.rds')
gc_mod_group_2 <- readRDS('./mySources/models/stress/s3_adv_int_nat/eye_track/gc_mod_group_2.rds')
gc_mod_group_3 <- readRDS('./mySources/models/stress/s3_adv_int_nat/eye_track/gc_mod_group_3.rds')
gc_mod_cond_0  <- readRDS('./mySources/models/stress/s3_adv_int_nat/eye_track/gc_mod_cond_0.rds')
gc_mod_cond_1  <- readRDS('./mySources/models/stress/s3_adv_int_nat/eye_track/gc_mod_cond_1.rds')
gc_mod_cond_2  <- readRDS('./mySources/models/stress/s3_adv_int_nat/eye_track/gc_mod_cond_2.rds')
gc_mod_cond_3  <- readRDS('./mySources/models/stress/s3_adv_int_nat/eye_track/gc_mod_cond_3.rds')
gc_mod_full    <- readRDS('./mySources/models/stress/s3_adv_int_nat/eye_track/gc_mod_full_0.rds')






# Random effects structure ----------------------------------------------------

mod_ot1 <- glmer(cbind(targetCount, 50 - targetCount) ~ 1 + ot1 +
                (1 + coda_sum + condition_sum + ot1 | participant),
                 control = glmerControl(optimizer = 'bobyqa'),
                 data = stress_gc_subset,family = "binomial")
mod_ot2 <- glmer(cbind(targetCount, 50 - targetCount) ~ 1 + ot1 + ot2 +
                (1 + coda_sum + condition_sum + (ot1 + ot2) | participant),
                 control = glmerControl(optimizer = "bobyqa"),
                 data = stress_gc_subset,family = "binomial")
mod_ot3 <- glmer(cbind(targetCount, 50 - targetCount) ~ 1 + ot1 + ot2 + ot3 +
                (1 + coda_sum + condition_sum + (ot1 + ot2) | participant),
                 control = glmerControl(optimizer = "bobyqa"),
                 data = stress_gc_subset,family = "binomial")

anova(mod_ot1, mod_ot2, mod_ot3, test = "Chisq")

# mod_ot3 gives singular fit
isSingular(mod_ot3)

# So we go with cubic model with no ot3 slope

# -----------------------------------------------------------------------------



# Test fixed effects ----------------------------------------------------------

# Base model
if(F){
  gc_mod_base <-
    lmer(eLog ~ 1 + (ot1 + ot2 + ot3) +
        (1 + coda_sum + condition_sum + (ot1 + ot2) | participant),
         control = lmerControl(optimizer = 'bobyqa'),
         data = stress_gc_subset, weights = 1/wts, REML = F)
  saveRDS(gc_mod_base, compress = 'xz',
          file = here("models", "stress", "s3_adv_int_nat", "eye_track",
                      "gca", "gca_mod_base.rds"))
}


# Add group effect on intercept
if(F){
  gc_mod_group_0 <-
    lmer(eLog ~ (ot1 + ot2 + ot3) + group +
        (1 + coda_sum + condition_sum + (ot1 + ot2) | participant),
         control = lmerControl(optimizer = 'bobyqa'),
         data = stress_gc_subset, weights = 1/wts, REML = F)
  saveRDS(gc_mod_group_0, compress = 'xz',
          file = here("models", "stress", "s3_adv_int_nat", "eye_track",
                      "gca", "gca_mod_group_0.rds"))
}

# Add group effect on slope
if(F){
 gc_mod_group_1 <-
   lmer(eLog ~ (ot1 + ot2 + ot3) + group +
        ot1:group +
       (1 + coda_sum + condition_sum + (ot1 + ot2) | participant),
        control = lmerControl(optimizer = 'bobyqa'),
        data = stress_gc_subset, weights = 1/wts, REML = F)
  saveRDS(gc_mod_group_1, compress = 'xz',
          file = here("models", "stress", "s3_adv_int_nat", "eye_track",
                      "gca", "gca_mod_group_1.rds"))
}

# Add group effect on quadratic poly
if(F){
  gc_mod_group_2 <-
    lmer(eLog ~ (ot1 + ot2 + ot3) + group +
         ot1:group + ot2:group +
        (1 + coda_sum + condition_sum + (ot1 + ot2) | participant),
         control = lmerControl(optimizer = 'bobyqa'),
         data = stress_gc_subset, weights = 1/wts, REML = F)
  saveRDS(gc_mod_group_2, compress = 'xz',
          file = here("models", "stress", "s3_adv_int_nat", "eye_track",
                      "gca", "gca_mod_group_2.rds"))
}

anova(gc_mod_base,
      gc_mod_group_0,
      gc_mod_group_1,
      gc_mod_group_2, test = 'Chisq')

# Add coda_sum to intercept

if(F){
gc_mod_coda_0 <-
  lmer(eLog ~ (ot1 + ot2+ ot3) + group + coda_sum +
       ot1:group + ot2:group +
      (1 + coda_sum + condition_sum + (ot1 + ot2) | participant),
       control = lmerControl(optimizer = 'bobyqa'),
       data = stress_gc_subset, weights = 1/wts, REML = F)
  saveRDS(gc_mod_coda_0, compress = 'xz',
          file = here("models", "stress", "s3_adv_int_nat", "eye_track",
                      "gca", "gca_mod_coda_0.rds"))
}

# Add coda_sum to linear slope
if(F){
gc_mod_coda_1 <-
  lmer(eLog ~ (ot1 + ot2 + ot3) + group + coda_sum +
       ot1:group + ot2:group +
       ot1:coda_sum +
      (1 + coda_sum + condition_sum + (ot1 + ot2) | participant),
       control = lmerControl(optimizer = 'bobyqa'),
       data = stress_gc_subset, weights = 1/wts, REML = F)
  saveRDS(gc_mod_coda_1, compress = 'xz',
          file = here("models", "stress", "s3_adv_int_nat", "eye_track",
                      "gca", "gca_mod_coda_1.rds"))
}

# Add coda_sum to quadratic effect
if(F){
  gc_mod_coda_2 <-
    lmer(eLog ~ (ot1 + ot2 + ot3) + group + coda_sum +
         ot1:group + ot2:group +
         ot1:coda_sum + ot2:coda_sum +
        (1 + coda_sum + condition_sum + (ot1 + ot2) | participant),
         control = lmerControl(optimizer = 'bobyqa'),
         data = stress_gc_subset, weights = 1/wts, REML = F)
  saveRDS(gc_mod_coda_2, compress = 'xz',
          file = here("models", "stress", "s3_adv_int_nat", "eye_track",
                      "gca", "gca_mod_coda_2.rds"))
}

# Add coda_sum interactions
if(F){
  gc_mod_coda_full <-
    lmer(eLog ~ 1 +
        (ot1 + ot2 + ot3) * group * coda_sum +
        (1 + coda_sum + condition_sum + (ot1 + ot2) | participant),
         control = lmerControl(optimizer = 'bobyqa',
                               optCtrl = list(maxfun = 2e5)),
         data = stress_gc_subset, weights = 1/wts, REML = F)
  saveRDS(gc_mod_coda_full, compress = 'xz',
          file = here("models", "stress", "s3_adv_int_nat", "eye_track",
                      "gca", "gca_mod_coda_full.rds"))
}

anova(gc_mod_coda_0, gc_mod_coda_1,
      gc_mod_coda_2, gc_mod_coda_full, test = "Chisq")


# Add condition_sum
if(F){
  gc_mod_cond_full <-
    lmer(eLog ~ 1 +
        (ot1 + ot2 + ot3) * group * coda_sum + condition_sum +
        (1 + coda_sum + condition_sum + (ot1 + ot2) | participant),
         control = lmerControl(optimizer = 'bobyqa',
                               optCtrl = list(maxfun = 1e5)),
         data = stress_gc_subset, weights = 1/wts, REML = F)
  saveRDS(gc_mod_cond_full, compress = 'xz',
          file = here("models", "stress", "s3_adv_int_nat", "eye_track",
                      "gca", "gca_mod_cond_full.rds"))
}

anova(gc_mod_coda_full, gc_mod_cond_full) # condition_sum doesnt improve fit
# 80662 80973 -40290    80580 2.5177      1     0.1126

# -----------------------------------------------------------------------------



# Model descriptives ----------------------------------------------------------

# options(scipen = 999)
# summary(gc_mod_coda_full)
# summary(gc_mod_cond_full)
# options(scipen = 0)


# Random effects:
#   Groups      Name          Variance Std.Dev. Corr
# participant (Intercept)   0.4144   0.6437
# coda_sum      0.1449   0.3806    0.35
# condition_sum 0.3300   0.5745   -0.23 -0.59
# ot1           4.3789   2.0926    0.55  0.22  0.00
# ot2           1.1133   1.0551    0.31 -0.09  0.08  0.63
# Residual                  8.5962   2.9319
# Number of obs: 14553, groups:  participant, 71

# Fixed effects:
#                           Estimate   Std. Error           df t value      Pr(>|t|)
# (Intercept)               0.830947     0.135819    69.044490   6.118 0.00000005045 ***
# ot1                       3.229903     0.459626    70.355245   7.027 0.00000000109 ***
# ot2                      -0.330291     0.269470    72.308417  -1.226      0.224286
# ot3                      -0.434760     0.160618 13565.526980  -2.707      0.006802 **

# groupla                  -0.197099     0.189878    68.703313  -1.038      0.302894
# ot1:groupla               0.128481     0.642572    70.020141   0.200      0.842101
# ot2:groupla               0.699797     0.374829    70.915831   1.867      0.066038 .
# ot3:groupla               0.240684     0.223831 13703.554781   1.075      0.282263

# groupint                 -0.354659     0.195708    68.141049  -1.812      0.074363 .
# ot1:groupint              0.382499     0.663009    69.736854   0.577      0.565856
# ot2:groupint              0.925203     0.386362    70.338245   2.395      0.019304 *
# ot3:groupint              0.275606     0.228059 13859.547204   1.208      0.226883

# coda_sum                  0.001020     0.077409    64.139979   0.013      0.989531
# ot1:coda_sum             -0.602962     0.162503 14284.932726  -3.710      0.000208 ***
# ot2:coda_sum             -0.100409     0.157232 14310.611836  -0.639      0.523090
# ot3:coda_sum              0.006918     0.158177 14373.445152   0.044      0.965117

# groupla:coda_sum          0.268366     0.108531    64.614123   2.473      0.016050 *
# ot1:groupla:coda_sum      1.018209     0.229069 14196.198866   4.445 0.00000885601 ***
# ot2:groupla:coda_sum      0.205036     0.219438 14114.633188   0.934      0.350130
# ot3:groupla:coda_sum     -0.131910     0.221266 14104.647841  -0.596      0.551076

# groupint:coda_sum         0.103884     0.111112    62.571776   0.935      0.353410
# ot1:groupint:coda_sum     0.024563     0.232649 14140.803135   0.106      0.915918
# ot2:groupint:coda_sum     0.258691     0.225141 14026.390655   1.149      0.250568
# ot3:groupint:coda_sum     0.008573     0.226816 14211.220139   0.038      0.969849

# -----------------------------------------------------------------------------




# Plot GCA --------------------------------------------------------------------

# create new df including the fitted model
data.comp <- data.frame(
  stress_gc_subset,GCA_Full = fitted(gc_mod_coda_full))



condition_namesGCAMod <- c(
  `stressed` = "Paroxytone",
  `unstressed` = "Oxytone",
  `0` = "CV",
  `1` = "CVC"
  )

(gca_full <- data.comp %>%
  ggplot(., aes(x = time_zero, y = eLog, color = group, shape = group)) +
  facet_grid(. ~ coda, labeller = as_labeller(condition_namesGCAMod)) +
  geom_smooth(method = 'gam', formula = y ~ poly(x, 3), se = F,
              show.legend = FALSE) +
  #stat_summary(aes(y = GCA_Full, color = group), fun.y = mean,
  #             geom = 'line', size = 1.4) +
  stat_summary(fun.data = mean_cl_boot, geom = 'pointrange',  size = 0.75,
               fun.args = list(conf.int = .95, B = 1000)) +
  stat_summary(fun.y = mean, geom = 'point', color = 'white',
               alpha = 0.3, size = 1.75) +
  scale_shape_manual(name = "", values = 17:15,
                     labels = c("SS", "LA", "INT")) +
  labs(x = "Time (ms) relative to target syllable offset",
       y = "Fixation empirical logit",
       caption = "Mean +/- 95% CI") +
  scale_color_brewer(palette = "Set1", name = "", guide = 'legend',
                     labels = c("SS", "LA", "INT")) +
  scale_x_continuous(breaks = c(0, 4, 8), labels = c("0", "200", "400")) +
  theme_grey(base_size = 15, base_family = "Times New Roman"))

# ggsave('stressP2.png', plot = gca_full, dpi = 600, device = "png",
#           path = "./mySources/figs/stress/s3_adv_int_nat/eye_track",
#           height = 4, width = 9, unit = "in")

# -----------------------------------------------------------------------------





















# Plot raw data ---------------------------------------------------------------

df_stress_50 <- stress50 %>%
  filter(.,  group %in% c('int', 'la', 'ss'),
             !participant %in% c("L01", "L02", "L03", "L04", "L05",
                                 "L06", "L07", "L08", "L09", "L10",
                                 "L15", "L20", "L21", "L22", "L23",
                                 "L26", "L30", "L31", "L33", "LA04",
                                 "LA06", "LA09", "LA14", "LA15", "LA19"))

df_stress_50 %>%
  group_by(group) %>%
  summarize(n = n_distinct(participant))



condition_names <- c(
  `stressed` = 'Paroxytone',
  `unstressed` = 'Oxytone',
  `0` = 'CV',
  `1` = 'CVC'
)

(df_stress_50 %>%
    na.omit(.) %>%
    filter(., time_zero >= -10, time_zero <= 20) %>%
    mutate(., group = fct_relevel(group, "ss", "la", "int")) %>%
    ggplot(., aes(x = time_zero, y = targetProp, color = group, shape = group)) +
    facet_grid(condition ~ coda, labeller = as_labeller(condition_names)) +
    geom_hline(yintercept = 0.5, color = 'white', size = 3) +
    geom_vline(xintercept = 0, color = 'grey40', lty = 3) +
    geom_vline(xintercept = 4, color = 'grey40', lty = 3) +
    stat_summary(fun.data = mean_cl_boot, geom = 'pointrange',  size = 0.75,
                 fun.args = list(conf.int = .95, B = 1000)) +
    stat_summary(fun.y = mean, geom = 'point', color = 'white',
                 alpha = 0.3, size = 1.75) +
    scale_shape_manual(name = "", values = 17:15,
                       labels = c("SS", "LA", "INT")) +
    scale_color_brewer(palette = 'Set1', name = "",
                       labels = c("SS", "LA", "INT")) +
    #scale_x_continuous(breaks = c(-4, 1, 6, 11, 16),
    #                   labels = c("-750", "-500", "-250", "0", "250")) +
    labs(y = 'Proportion of target fixations',
         x = 'Time relative to target syllable offset (ms)',
         caption = "Mean +/- 95% CI") +
    coord_cartesian(ylim = c(0, 1)) +
    annotate("text", x = -0.7, y = 0.02, label = 'Target syllable offset',
             angle = 90, size = 3, hjust = 0) +
    annotate("text", x = 3.3, y = 0.02, label = '200ms after target offset',
             angle = 90, size = 3, hjust = 0) +
    theme_grey(base_size = 16, base_family = "Times") -> stressP1)

# ggsave('stressP1.png', plot = stressP1, dpi = 600, device = "png",
# path = "./mySources/figs/stress/s3_adv_int_nat/eye_track",
# height = 4, width = 9, unit = 'in')

