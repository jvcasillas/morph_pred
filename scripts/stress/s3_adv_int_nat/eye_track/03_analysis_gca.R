

# read data
df_stress <- read_csv("./mySources/data/clean/stressBIN10iaClean.csv") %>%
  filter(., corr == 1,
            group %in% c('int', 'la', 'ss'),
            !participant %in% c('L01', 'L02', 'L03', 'L04', 'L05',
                                 'L06', 'L07', 'L08', 'L09', 'L10',
                                 'L15', 'L20', 'L21', 'L22', 'L23',
                                 'L26', 'L30', 'L31', 'L33', 'LA04',
                                 'LA06', 'LA07', 'LA14'))


# Load wm data and combine with stress_subset_0 (proportion data)
# in order to add working memory as a covariate
wm_df <- read_csv("./mySources/data/raw/wm_all.csv") %>%
         filter(., !(group %in% c("HS", "L")),
                   !(participant %in% c('L01', 'L02', 'L03', 'L04', 'L05',
                                        'L06', 'L07', 'L08', 'L09', 'L10',
                                        'L15', 'L20', 'L21', 'L22', 'L23',
                                        'L26', 'L30', 'L31', 'L33', 'La04',
                                        'LA06', 'LA07', 'LA14')))

scale_this <- function(x) as.vector(scale(x))





## @knitr stressGCA


####################################################
# GROWTH CURVE ANALYSIS                            #
# - Question 1: Are the groups different from      #
#   each other in when they begin to fixate        #
#   on the target?                                 #
#     - test 3 groups at each level of 'condition' #
#     - hypothesis: SS has steeper slope for both  #
#       conditions                                 #
# - Question 2: W/in groups, is the there a        #
#   difference between oxytone/paroxyton           #
#   items?                                         #
#     - test oxytone vs. paroxytone for each group #
#     - hypothesis: steeper slope/earlier break in #
#       oxytone condition                          #
####################################################

# Prep
# - subset using time course
#    - we want to find the earliest target word onset
#    - we will substract 'a little more' and use that
#      as the starting point
#    - we already have this information in the 'twOnsetAdj'
#      dataframe

# print(twOnsetAdj)
# min(twOnsetAdj$twOnsetAdj)

#    - the lowest target word onset is: 'firma' @ -51.3
#    - thus we can use -60 binAdj as the starting point
#      and be sure we are including the entire target word
#    - we can also do a little higher to lighten the models (-50)
stress_gc_subset <- filter(df_stress_50, binTonsetAlign >= 33 &
                                         binTonsetAlign <= 43) %>% as.data.frame


# - Readjust time course
#    - now we will make the time course positive, starting at 1
#    - to do this we will add the lowest 'binAdj' value to each
#      bin, plus 1 (to avoid starting at 0)

# lsrl_gc_subset$binGC <- lsrl_gc_subset$binAdj + 51
stress_gc_subset$binGC <- stress_gc_subset$binTonsetAlign - 32


# - Now we add higher order polynomials for analyses

t <- poly(min(stress_gc_subset$binGC):max(stress_gc_subset$binGC), 4)
stress_gc_subset[, paste('ot', 1:4, sep = "")] <- t[stress_gc_subset$binGC, 1:4]

# glimpse(lsrl_gc_subset)






####################################################
# - Question 1: Are the groups different from      #
#   each other in when they begin to fixate        #
#   on the target?                                 #
#     - test 3 groups at each level of 'condition' #
#     - hypothesis: SS has steeper slope for both  #
#       conditions                                 #
####################################################


# Set SS as reference level
stress_gc_subset$group <- factor(stress_gc_subset$group, levels = c("ss", "la", "int"))

# Set condition as factor and set conding to contrast
stress_gc_subset$condition <- as.factor(stress_gc_subset$condition)

# contrasts(stress_gc_subset$condition)
# contrasts(stress_gc_subset$group)

stress_gc_subset$conditionSum <- C(stress_gc_subset$condition, sum)
# contrasts(stress_gc_subset$conditionSum)

stress_gc_subset$coda <- as.factor(stress_gc_subset$coda)
stress_gc_subset$codaSum <- C(stress_gc_subset$coda, sum)
# contrasts(stress_gc_subset$codaSum)


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


# random effects structure
mod_ot1 <- lmer(eLog ~ ot1 + (1 + ot1 | participant) +
                (1 + ot1 | participant:codaSum),
                control = lmerControl(optimizer = 'bobyqa'),
                data = stress_gc_subset, weights = 1/wts, REML = F)
mod_ot2 <- lmer(eLog ~ ot1 + ot2 + (1 + (ot1+ot2) | participant) +
                (1 + (ot1+ot2) | participant:codaSum),
                control = lmerControl(optimizer = 'bobyqa'),
                data = stress_gc_subset, weights = 1/wts, REML = F)
mod_ot3 <- lmer(eLog ~ ot1 + ot2 + ot3 + (1 + (ot1+ot2+ot3) | participant) +
                (1 + (ot1+ot2+ot3) | participant:codaSum),
                control = lmerControl(optimizer = 'bobyqa'),
                data = stress_gc_subset, weights = 1/wts, REML = F)

anova(mod_ot1, mod_ot2, mod_ot3) # best fit with quadratic poly

mod_ot2_simp <- lmer(eLog ~ ot1 + ot2 + (1 + (ot1+ot2) | participant),
                control = lmerControl(optimizer = 'bobyqa'),
                data = stress_gc_subset, weights = 1/wts, REML = F)

mod_ot2_cond <- lmer(eLog ~ ot1 + ot2 + (1 + (ot1+ot2) | participant) +
                (1 + (ot1+ot2) | participant:conditionSum),
                control = lmerControl(optimizer = 'bobyqa'),
                data = stress_gc_subset, weights = 1/wts, REML = F)

anova(mod_ot2_simp, mod_ot2_cond) # cond is better than just part

mod_ot2_coda <- lmer(eLog ~ ot1 + ot2 + (1 + (ot1+ot2) | participant) +
                (1 + (ot1+ot2) | participant:codaSum),
                control = lmerControl(optimizer = 'bobyqa'),
                data = stress_gc_subset, weights = 1/wts, REML = F)

anova(mod_ot2_simp, mod_ot2_coda) # coda is better than simp

mod_ot2_max <- lmer(eLog ~ ot1 + ot2 + (1 + (ot1+ot2) | participant) +
                (1 + (ot1+ot2) | participant:codaSum:conditionSum),
                control = lmerControl(optimizer = 'bobyqa'),
                data = stress_gc_subset, weights = 1/wts, REML = F)

anova(mod_ot2_coda, mod_ot2_max) # coda is better than simp








# Base model
if(T){
  gc_mod_base <- lmer(eLog ~ (ot1+ot2) +
                 (1 + (ot1+ot2) | participant) +
                 (1 + (ot1+ot2) | participant:codaSum),
                 control = lmerControl(optimizer = 'bobyqa'),
                 data = stress_gc_subset, weights = 1/wts, REML = F)
  saveRDS(gc_mod_base, file = "./mySources/models/stress/s3_adv_int_nat/eye_track/gc_mod_base.rds", compress = 'xz')
}




# Add group effect on intercept
if(T){
  gc_mod_group_0 <- lmer(eLog ~ (ot1+ot2) + group +
                    (1 + (ot1+ot2) | participant) +
                    (1 + (ot1+ot2) | participant:codaSum),
                    control = lmerControl(optimizer = 'bobyqa'),
                    data = stress_gc_subset, weights = 1/wts, REML = F)
  saveRDS(gc_mod_group_0, file = "./mySources/models/stress/s3_adv_int_nat/eye_track/gc_mod_group_0.rds", compress = 'xz')
}

# Add group effect on slope
if(T){
 gc_mod_group_1 <- lmer(eLog ~ (ot1+ot2) + group +
                   ot1:group +
                   (1 + (ot1+ot2) | participant) +
                   (1 + (ot1+ot2) | participant:codaSum),
                   control = lmerControl(optimizer = 'bobyqa'),
                   data = stress_gc_subset, weights = 1/wts, REML = F)
  saveRDS(gc_mod_group_1, file = "./mySources/models/stress/s3_adv_int_nat/eye_track/gc_mod_group_1.rds", compress = 'xz')
}

# Add group effect on quadratic poly
if(T){
  gc_mod_group_2 <- lmer(eLog ~ (ot1+ot2) + group +
                    ot1:group + ot2:group +
                    (1 + (ot1+ot2) | participant) +
                    (1 + (ot1+ot2) | participant:codaSum),
                    control = lmerControl(optimizer = 'bobyqa'),
                    data = stress_gc_subset, weights = 1/wts, REML = F)
  saveRDS(gc_mod_group_2, file = "./mySources/models/stress/s3_adv_int_nat/eye_track/gc_mod_group_2.rds", compress = 'xz')
}

anova(gc_mod_base,
      gc_mod_group_0,
      gc_mod_group_1,
      gc_mod_group_2, test = 'Chisq')




# full mod
if(T){
gc_mod_full_0 <- lmer(eLog ~ (ot1+ot2+ot3) * group * codaSum +
               (1 + (ot1+ot2+ot3) | participant) +
               # (1 + (ot1+ot2) | participant:codaSum) +
               (1 + (ot1+ot2+ot3) | participant:conditionSum),
               control = lmerControl(optimizer = 'bobyqa'),
               data = stress_gc_subset, weights = 1/wts, REML = F)
  saveRDS(gc_mod_full_0, file = "./mySources/models/stress/s3_adv_int_nat/eye_track/gc_mod_full_0.rds", compress = 'xz')
}

if(T){
gc_mod_full_1 <- lmer(eLog ~ (ot1+ot2) * group * codaSum + conditionSum +
               (1 + (ot1+ot2) | participant) +
               # (1 + (ot1+ot2) | participant:codaSum) +
               (1 + (ot1+ot2) | participant:conditionSum),
               control = lmerControl(optimizer = 'bobyqa'),
               data = stress_gc_subset, weights = 1/wts, REML = F)
  saveRDS(gc_mod_full_1, file = "./mySources/models/stress/s3_adv_int_nat/eye_track/gc_mod_full_1.rds", compress = 'xz')
}

anova(gc_mod_full_0, gc_mod_full_1) # no main effect of cond
# ..1    32 60048 60281 -29992    59984 3.386      1    0.06575 .



# summary(gc_mod_full_0)



# Fixed effects:
#                         Estimate Std. Error         df t value Pr(>|t|)
# (Intercept)              0.90951    0.14791   68.09840   6.149 4.63e-08 ***
# ot1                      0.31167    0.33240 8438.55185   0.938 0.348457
# ot2                     -0.16466    0.23044   54.49578  -0.715 0.477930

# groupla                 -0.38784    0.19936   68.23686  -1.945 0.055843 .
# ot1:groupla             -0.80774    0.50262   67.96785  -1.607 0.112680
# ot2:groupla              1.04344    0.31158   55.32223   3.349 0.001466 **

# groupint                -0.09539    0.26477   68.37041  -0.360 0.719754
# ot1:groupint             0.07018    0.67222   70.37943   0.104 0.917153
# ot2:groupint            -0.20198    0.41851   59.32022  -0.483 0.631149

# codaSum1                 0.07093    0.05074 9317.63184   1.398 0.162166
# ot1:codaSum1             0.27608    0.17746 7863.15846   1.556 0.119827
# ot2:codaSum1             0.06745    0.17106 6075.19780   0.394 0.693379

# groupla:codaSum1        -0.17768    0.06837 9279.03361  -2.599 0.009371 **
# ot1:groupla:codaSum1    -0.33628    0.23924 7949.21744  -1.406 0.159880
# ot2:groupla:codaSum1     0.07746    0.23258 6122.17549   0.333 0.739124

# groupint:codaSum1       -0.31439    0.09270 9243.38468  -3.392 0.000698 ***
# ot1:groupint:codaSum1    3.30410    0.37224   67.43673   8.876 6.05e-13 ***
# ot2:groupint:codaSum1    0.04064    0.31389 7094.47346   0.129 0.896999




# create new df including the fitted model
data.comp <- data.frame(na.omit(stress_gc_subset),
                        GCA_Full = fitted(gc_mod_full_0))
# glimpse(data.comp)


condition_namesGCAMod <- c(
                    `0` = "CV",
                    `1` = "CVC"
                    )

(gca_full <- data.comp %>%
  ggplot(., aes(x = binTonsetAlign, y = eLog, color = group, shape = group)) +
  facet_grid(. ~ coda, labeller = as_labeller(condition_namesGCAMod)) +
  geom_smooth(method = 'gam', formula = y ~ poly(x, 3), se = F, show.legend = FALSE) +
  # stat_summary(aes(y = GCA_Full, color = group), fun.y = mean, geom = 'line', size = 0.4) +
  stat_summary(fun.data = mean_cl_boot, geom = 'pointrange',  size = 0.75,
               fun.args = list(conf.int = .95, B = 1000)) +
  stat_summary(fun.y = mean, geom = 'point', color = 'white', alpha = 0.3, size = 1.75) +
  scale_shape_manual(name = "", values = 17:15, labels = c("M", "NIN", "IN")) +
  labs(x = "Time relative to target syllable offset (ms)", y = "Fixation empirical logit",
       caption = "Mean +/- 95% CI") +
  scale_color_brewer(palette = "Set1", name = "", guide = 'legend',
                     labels = c("M", "NIN", "IN")) +
  scale_x_continuous(breaks = c(33, 38, 43), labels = c("-250", "0", "250")) +
  theme_grey(base_size = 15, base_family = "Times New Roman"))

# ggsave('stressP2.png', plot = gca_full, dpi = 600, device = "png",
#           path = "./mySources/figs/stress/s3_adv_int_nat/eye_track",
#           height = 4, width = 9, unit = "in")













# OLD OUTPUT

# Fixed effects:
#                         Estimate Std. Error         df t value Pr(>|t|)
# (Intercept)              0.89485    0.14747   68.00000   6.068 6.41e-08 ***
# ot1                      3.33309    0.37630   67.00000   8.858 7.25e-13 ***
# ot2                     -0.16281    0.23645   55.00000  -0.689  0.49401
#
# groupla                 -0.37843    0.19871   68.00000  -1.904  0.06107 .
# ot1:groupla             -0.84030    0.50822   67.00000  -1.653  0.10291
# ot2:groupla              1.04784    0.31978   55.00000   3.277  0.00182 **
#
# groupint                -0.11347    0.26369   68.00000  -0.430  0.66832
# ot1:groupint             0.07456    0.67747   69.00000   0.110  0.91268
# ot2:groupint            -0.18355    0.42840   59.00000  -0.428  0.66989
#
# codaSum1                 0.07384    0.05076 9377.00000   1.455  0.14578
# ot1:codaSum1             0.25070    0.17745 8565.00000   1.413  0.15777
# ot2:codaSum1             0.05033    0.17117 6255.00000   0.294  0.76874
#
# groupla:codaSum1        -0.17708    0.06821 9430.00000  -2.596  0.00945 **
# ot1:groupla:codaSum1    -0.33634    0.23883 8657.00000  -1.408  0.15909
# ot2:groupla:codaSum1     0.09856    0.23274 6338.00000   0.423  0.67196
#
# groupint:codaSum1       -0.29203    0.09172 9378.00000  -3.184  0.00146 **
# ot1:groupint:codaSum1    0.32001    0.32918 8746.00000   0.972  0.33102
# ot2:groupint:codaSum1    0.07729    0.31381 7054.00000   0.246  0.80545


# ot2:groupla              1.04784    0.31978   55.00000   3.277  0.00182 **
# groupla:codaSum1        -0.17708    0.06821 9430.00000  -2.596  0.00945 **
# groupint:codaSum1       -0.29203    0.09172 9378.00000  -3.184  0.00146 **












df_stress_50 <- stress50 %>%
  filter(.,  group %in% c('int', 'la', 'ss'),
             !participant %in% c('L01', 'L02', 'L03', 'L04', 'L05',
                                 'L06', 'L07', 'L08', 'L09', 'L10',
                                 'L15', 'L20', 'L21', 'L22', 'L23',
                                 'L26', 'L30', 'L31', 'L33', 'LA04',
                                 'LA06', 'LA07', 'LA14'))

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
    facet_grid(. ~ coda, labeller = as_labeller(condition_names)) +
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

