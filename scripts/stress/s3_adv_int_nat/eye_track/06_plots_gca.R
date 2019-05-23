# Time course plots -----------------------------------------------------------
#
# - This script plots the raw eye tracking time course data and the
#   model estimates from a growth curve analysis
#
# -----------------------------------------------------------------------------



# Source scripts and load models ----------------------------------------------

source(here::here("scripts", "01_load_data.R"))

gca_full_mod_full <- readRDS(here("models", "stress", "s3_adv_int_nat",
                              "eye_track", "gca", "gca_full_mod_full.rds"))

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
    ggplot(., aes(x = time_zero, y = targetProp, fill = group, shape = group)) +
    facet_grid(condition ~ coda, labeller = as_labeller(condition_names)) +
    geom_hline(yintercept = 0.5, color = 'white', size = 3) +
    geom_vline(xintercept = 0, color = 'grey40', lty = 3) +
    geom_vline(xintercept = 4, color = 'grey40', lty = 3) +
    stat_summary(fun.data = mean_cl_boot, geom = 'pointrange', size = 0.5,
                 stroke = 0.5, pch = 21) +
    scale_fill_brewer(palette = 'Set1', name = "",
                       labels = c("SS", "LA", "INT")) +
    scale_x_continuous(breaks = c(-10, 0, 10, 20),
                       labels = c("-500", "0", "500", "1000")) +
    labs(y = 'Proportion of target fixations',
         x = 'Time relative to target syllable offset (ms)',
         caption = "Mean +/- 95% CI") +
    annotate("text", x = 3.3, y = 0.02, label = '200ms',
             angle = 90, size = 3, hjust = 0) +
    theme_grey(base_size = 12, base_family = "Times") -> stressP1)

ggsave('stressP1.png', plot = stressP1, dpi = 600, device = "png",
path = here("figs", "stress", "s3_adv_int_nat", "eye_track"),
height = 6, width = 9, unit = 'in')

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
  stress_gca_plot_subset, GCA_Full = fitted(gca_full_mod_int_3))



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
    stat_summary(aes(y = GCA_Full, color = group), fun.y = mean,
                 geom = 'line', size = 1.4) +
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

new_dat_all <- stress_gc_subset %>%
  dplyr::select(group, time_zero, ot1:ot3, coda_sum, condition_sum) %>%
  distinct

# Get models predictions and SE
fits_all <- predictSE(gca_full_mod_int_3, new_dat_all) %>%
  as_tibble %>%
  bind_cols(new_dat_all) %>%
  rename(se = se.fit) %>%
  mutate(ymin = fit - se, ymax = fit + se)

filter(fits_all, time_zero == 4) %>%
  select(group, coda = coda_sum, cond = condition_sum,
         elog = fit, elog_lb = ymin, elog_ub = ymax) %>%
  mutate(prob = plogis(elog),
         prob_lb = plogis(elog_lb),
         prob_ub = plogis(elog_ub)) %>%
  arrange(group)


# Within group differences
fits_all %>%
  mutate(coda = if_else(coda_sum == 1, "CV", "CVC"),
         condition = if_else(condition_sum == 1, "Paroxytone", "Oxytone"),
         condition = fct_relevel(condition, "Paroxytone")) %>%
  ggplot(., aes(x = time_zero, y = fit, ymax = ymax, ymin = ymin,
                fill = coda, color = coda)) +
  facet_grid(group ~ condition) +
  geom_hline(yintercept = 0, lty = 3) +
  geom_vline(xintercept = 4, lty = 3) +
  geom_ribbon(alpha = 0.2, color = NA, show.legend = F) +
  geom_point(aes(shape = coda), size = 3.5, show.legend = F) +
  geom_line(size = 0.9) +
  scale_x_continuous(breaks = c(-4, 0, 4, 8, 12),
                     labels = c("-200", "0", "200", "400", "600")) +
  scale_color_brewer(palette = "Set1", name = "Lexical stress") +
  labs(x = "Time (ms) relative to target syllable offset",
       y = "Empirical logit of looks to target") +
  theme_big

# Comparisons with natives
fits_all %>%
  mutate(coda = if_else(coda_sum == 1, "CV", "CVC"),
         condition = if_else(condition_sum == 1, "Paroxytone", "Oxytone"),
         condition = fct_relevel(condition, "Paroxytone")) %>%
  ggplot(., aes(x = time_zero, y = fit, ymax = ymax, ymin = ymin,
                fill = group, color = group)) +
  facet_grid(coda ~ condition) +
  geom_hline(yintercept = 0, lty = 3) +
  geom_vline(xintercept = 4, lty = 3) +
  geom_ribbon(alpha = 0.2, color = NA, show.legend = F) +
  geom_point(aes(shape = group), size = 3.5, show.legend = F) +
  geom_line(size = 0.9) +
  scale_x_continuous(breaks = c(-4, 0, 4, 8, 12),
                     labels = c("-200", "0", "200", "400", "600")) +
  scale_color_brewer(palette = "Set1", name = "Group") +
  labs(x = "Time (ms) relative to target syllable offset",
       y = "Empirical logit of looks to target") +
  theme_big




# -----------------------------------------------------------------------------








# Effect size plot ------------------------------------------------------------

bind_rows(
  data.comp %>%
    group_by(participant, group, time_zero, coda) %>%
    summarize(competition = mean(eLog)) %>%
    spread(coda, competition) %>%
    mutate(Effect = `1` - `0`,
           condition_type = "Syllable structure",
           fit_type = "raw") %>%
    select(-`0`, -`1`),
  data.comp %>%
    group_by(participant, group, time_zero, condition) %>%
    summarize(competition = mean(eLog)) %>%
    spread(condition, competition) %>%
    mutate(Effect = unstressed - stressed,
           condition_type = "Lexical stress",
           fit_type = "raw") %>%
    select(-unstressed, -stressed),
  data.comp %>%
    group_by(participant, group, time_zero, coda) %>%
    summarize(competition = mean(GCA_Full)) %>%
    spread(coda, competition) %>%
    mutate(Effect = `1` - `0`,
           condition_type = "Syllable structure",
           fit_type = "model") %>%
    select(-`0`, -`1`),
  data.comp %>%
    group_by(participant, group, time_zero, condition) %>%
    summarize(competition = mean(GCA_Full)) %>%
    spread(condition, competition) %>%
    mutate(Effect = unstressed - stressed,
           condition_type = "Lexical stress",
           fit_type = "model") %>%
    select(-unstressed, -stressed)) %>%
  spread(fit_type, Effect) %>%
  #filter(time_zero < 10) %>%
  ggplot(., aes(x = time_zero, y = raw, color = group)) +
    facet_grid(. ~ condition_type) +
    geom_hline(yintercept = 0, color = "white", size = 3) +
    geom_vline(xintercept = 4, color = "white", size = 3) +
    stat_summary(fun.data = mean_se, geom = "pointrange") +
    #stat_summary(aes(y = model, color = group), fun.y = mean,
    #             geom = 'line', size = 1.4) +
    geom_smooth(method = 'gam', formula = y ~ poly(x, 3), se = F,
                show.legend = FALSE) +
    labs(x = "Time (ms) relative to target syllable offset",
         y = "Effect", caption = "Mean +/- 95% CI") +
    scale_color_brewer(palette = "Set1", name = "", guide = 'legend',
                       labels = c("SS", "LA", "INT")) +
    scale_x_continuous(breaks = c(-4, 0, 4, 8, 12),
                       labels = c("-200", "0", "200", "400", "600")) +
    theme_grey(base_size = 15, base_family = "Times New Roman")

# -----------------------------------------------------------------------------











