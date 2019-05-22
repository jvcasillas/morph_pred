# Time course plots -----------------------------------------------------------
#
# - This script plots the raw eye tracking time course data and the
#   model estimates from a growth curve analysis
#
# -----------------------------------------------------------------------------



# Source scripts and load models ----------------------------------------------

source(here::here("scripts", "01_load_data.R"))

gca_mod_int_3 <- readRDS(here("models", "stress", "s3_adv_int_nat",
                              "eye_track", "gca", "gca_mod_int_6.rds"))

# -----------------------------------------------------------------------------






# Plot raw data ---------------------------------------------------------------

df_stress_50 <- stress50 %>%
  filter(.,  group %in% c('int', 'la', 'ss'),
            !participant %in% c("L01", "L02", "L03", "L04", "L05",
                                "L06", "L07", "L08", "L09", "L10",
                                "L15", "L20", "L21", "L22", "L23",
                                "L26", "L30", "L31", "L33", "LA04",
                                "LA06", "LA09", "LA14", "LA15", "LA19"))

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

# -----------------------------------------------------------------------------




# Plot GCA --------------------------------------------------------------------

# create new df including the fitted model
stress_gca_plot_subset <- df_stress_50 %>%
  filter(time_zero >= -4 & time_zero <= 12) %>%
  mutate(., group = fct_relevel(group, "ss", "la", "int"),
            condition_sum = if_else(condition == "stressed", 1, -1),
            coda_sum = if_else(coda == 1, 1, -1)) %>%
  poly_add_columns(., time_zero, degree = 3, prefix = "ot")

data.comp <- data.frame(
  stress_gca_plot_subset, GCA_Full = fitted(gca_mod_full))



condition_namesGCAMod <- c(
  `stressed` = "Paroxytone",
  `unstressed` = "Oxytone",
  `0` = "CV",
  `1` = "CVC"
)

(gca_full <- data.comp %>%
    ggplot(., aes(x = time_zero, y = eLog, color = group, shape = group)) +
    facet_grid(coda ~ condition,
               labeller = as_labeller(condition_namesGCAMod)) +
    geom_hline(yintercept = 0, color = "white", size = 3) +
    geom_vline(xintercept = 4, color = "white", size = 3) +
    geom_smooth(method = 'gam', formula = y ~ poly(x, 3), se = F,
                show.legend = FALSE) +
    #stat_summary(aes(y = GCA_Full, color = group), fun.y = mean,
    #             geom = 'line', size = 1.4) +
    stat_summary(fun.data = mean_se, geom = 'pointrange',  size = 0.75) +
    stat_summary(fun.y = mean, geom = 'point', size = 1.75, color = "white",
                 alpha = 0.3) +
    scale_shape_manual(name = "", values = 17:15,
                       labels = c("SS", "LA", "INT")) +
    labs(x = "Time (ms) relative to target syllable offset",
         y = "Fixation empirical logit",
         caption = "Mean +/- 95% CI") +
    scale_color_brewer(palette = "Set1", name = "", guide = 'legend',
                       labels = c("SS", "LA", "INT")) +
    scale_x_continuous(breaks = c(-4, 0, 4, 8, 12),
                       labels = c("-200", "0", "200", "400", "600")) +
    theme_grey(base_size = 15, base_family = "Times New Roman"))

# ggsave('stressP2.png', plot = gca_full, dpi = 600, device = "png",
#           path = "./mySources/figs/stress/s3_adv_int_nat/eye_track",
#           height = 4, width = 9, unit = "in")

# -----------------------------------------------------------------------------




# Individual plots ------------------------------------------------------------

data.comp %>%
  ggplot(., aes(x = time_zero, y = eLog, color = factor(coda),
                shape = factor(coda))) +
  facet_grid(group ~ condition) +
  geom_hline(yintercept = 0, color = "white", size = 3) +
  geom_vline(xintercept = 4, color = "white", size = 3) +
  geom_smooth(method = 'gam', formula = y ~ poly(x, 3), se = F,
              show.legend = FALSE) +
  stat_summary(fun.data = mean_se, geom = 'pointrange',  size = 0.75) +
  stat_summary(fun.y = mean, geom = 'point', size = 1.75, color = "white",
               alpha = 0.3) +
  scale_shape_manual(name = "", values = 17:16,
                     labels = c("CV", "CVC")) +
  labs(x = "Time (ms) relative to target syllable offset",
       y = "Fixation empirical logit",
       caption = "Mean +/- 95% CI") +
  scale_color_brewer(palette = "Set1", name = "", guide = 'legend',
                     labels = c("CV", "CVC")) +
  scale_x_continuous(breaks = c(-4, 0, 4, 8, 12),
                     labels = c("-200", "0", "200", "400", "600")) +
  theme_grey(base_size = 15, base_family = "Times New Roman")
