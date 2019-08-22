# Time course plots -----------------------------------------------------------
#
# - This script plots the raw eye tracking time course data and the
#   model estimates from a growth curve analysis
#
# -----------------------------------------------------------------------------



# Source scripts and load models ----------------------------------------------

source(here::here("scripts", "02_load_data.R"))

# Get path to saved models
gca_mods_path  <- here("models", "stress", "s3_adv_int_nat", "eye_track", "gca")

# Load models as list and store full mod to global env
load(paste0(gca_mods_path, "/full_mods_lang_learn.Rdata"))
load(paste0(gca_mods_path, "/ind_mods_phon.Rdata"))
load(paste0(gca_mods_path, "/nested_model_comparisons_phon.Rdata"))
load(paste0(gca_mods_path, "/model_preds_phon.Rdata"))



# Set path for saving figs
figs_path <- here("figs", "stress", "s3_adv_int_nat", "eye_track", "lang_learn")

# -----------------------------------------------------------------------------











# Plot GCA --------------------------------------------------------------------

# Within group differences
stress_phon_p1 <- model_preds_phon$fits_all_phon %>%
  mutate(`Phonotactic Frequency` = as.factor(phon_std),
         coda = if_else(coda_sum == 1, "CV", "CVC"),
         condition = if_else(condition_sum == 1, "Paroxytone", "Oxytone"),
         condition = fct_relevel(condition, "Paroxytone")) %>%
  ggplot(., aes(x = time_zero, y = fit, ymax = ymax, ymin = ymin,
                fill = coda, color = coda, lty = `Phonotactic Frequency`)) +
  facet_grid(group ~ condition) +
  geom_hline(yintercept = 0, size = 3, color = "white") +
  geom_vline(xintercept = 4, size = 3, color = "white") +
  geom_ribbon(alpha = 0.2, color = NA, show.legend = F) +
  geom_line(size = 0.35) +
  scale_x_continuous(breaks = c(-4, 0, 4, 8, 12),
                     labels = c("-200", "0", "200", "400", "600")) +
  scale_color_brewer(palette = "Set1", name = "Syllable structure") +
  labs(x = "Time (ms) relative to target syllable offset",
       y = "Empirical logit of looks to target") +
  theme_grey(base_size = 10, base_family = "Times") + legend_adj_3


# Comparisons with natives
stress_phon_p2 <- model_preds_phon$fits_all_phon %>%
  mutate(coda = if_else(coda_sum == 1, "CV", "CVC"),
         condition = if_else(condition_sum == 1, "Paroxytone", "Oxytone"),
         condition = fct_relevel(condition, "Paroxytone")) %>%
  filter(phon_std == 0) %>%
  ggplot(., aes(x = time_zero, y = fit, ymax = ymax, ymin = ymin,
                fill = group, color = group)) +
  facet_grid(coda ~ condition) +
  geom_hline(yintercept = 0, size = 3, color = "white") +
  geom_vline(xintercept = 4, size = 3, color = "white") +
  geom_ribbon(alpha = 0.2, color = NA, show.legend = F) +
  geom_line(size = 0.75) +
  geom_point(aes(shape = group), color = "black", size = 1.3, show.legend = F) +
  geom_point(aes(shape = group), size = 0.85, show.legend = F) +
  scale_x_continuous(breaks = c(-4, 0, 4, 8, 12),
                     labels = c("-200", "0", "200", "400", "600")) +
  scale_color_brewer(palette = "Dark2", name = "Group") +
  scale_fill_brewer(palette = "Dark2", name = "Group") +
  labs(x = "Time (ms) relative to target syllable offset",
       y = "Empirical logit of looks to target") +
  theme_grey(base_size = 10, base_family = "Times") + legend_adj_3

ggsave(paste0(figs_path, "/stress_phon_p1.png"), stress_phon_p1, width = 150,
       height = 120, units = "mm", dpi = 600)
ggsave(paste0(figs_path, "/stress_phon_p1.eps"), stress_phon_p1, width = 150,
       height = 120, units = "mm", dpi = 600, device = cairo_ps)
ggsave(paste0(figs_path, "/stress_phon_p2.png"), stress_phon_p2, width = 150,
       height = 120, units = "mm", dpi = 600)
ggsave(paste0(figs_path, "/stress_phon_p2.eps"), stress_phon_p2, width = 150,
       height = 120, units = "mm", dpi = 600, device = cairo_ps)

# -----------------------------------------------------------------------------
