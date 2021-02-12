
# Graphs to help interpret results of the GCA for FIRMA/FIRMÓ + wm
# First graph is for stress, you will need to change the group to get all 3 groups (ss, la, int)
# Second graph is for syllabic structure, same thing with the groups


# Graphs for stress (condition_sum)

test_mod <-
  lmer(eLog ~ 1 + (ot1 + ot2 + ot3) * condition_sum * wm_std +
         (1 + condition_sum + (ot1 + ot2 + ot3) | participant) +
         (1 + ot1 + ot2 + ot3 | target),
       control = lmerControl(optimizer = 'bobyqa'), REML = F,
       data = filter(stress_gc_subset, group %in% c("ss")))

group.labs <- c("non-interpreters", "interpreters", "monolinguals")
names(group.labs) <- c("la", "int", "ss")



# Create design dataframe for predictions
new_dat_wm_minus1 <- stress_gc_subset %>%
  dplyr::select(group, time_zero, ot1:ot3, condition_sum, wm_std) %>%
  distinct(.) %>%
  mutate(wm_std = -1) %>% filter(group %in% c("ss"))
new_dat_wm_0 <- stress_gc_subset %>%
  dplyr::select(group, time_zero, ot1:ot3, condition_sum, wm_std) %>%
  distinct(.) %>%
  mutate(wm_std = 0) %>% filter(group %in% c("ss"))
new_dat_wm_1 <- stress_gc_subset %>%
  dplyr::select(group, time_zero, ot1:ot3, condition_sum, wm_std) %>%
  distinct(.) %>%
  mutate(wm_std = 1) %>% filter(group %in% c("ss"))

# Get model predictions and SE
test_fits <-
  bind_rows(
    predictSE(test_mod, new_dat_wm_minus1) %>%
      as_tibble %>%
      bind_cols(new_dat_wm_minus1),
    predictSE(test_mod, new_dat_wm_0) %>%
      as_tibble %>%
      bind_cols(new_dat_wm_0),
    predictSE(test_mod, new_dat_wm_1) %>%
      as_tibble %>%
      bind_cols(new_dat_wm_1)
  ) %>%
  rename(se = se.fit) %>%
  mutate(ymin = fit - se, ymax = fit + se)





test_fits %>%
  ggplot(., aes(x = time_zero, y = fit, ymax = ymax, ymin = ymin,
                fill = factor(condition_sum), color = factor(condition_sum), lty = factor(wm_std))) +
  facet_grid(. ~ group, labeller = labeller(group = group.labs)) +
  geom_hline(yintercept = 0, size = 3, color = "white") +
  geom_vline(xintercept = 4, size = 3, color = "white") +
  geom_line(size = 1) +
  scale_x_continuous(breaks = c(-4, 0, 4, 8, 12),
                     labels = c("-200", "0", "200", "400", "600")) +
  scale_color_brewer(palette = "Set1", name = "Stress", labels = c("oxytone", "paroxytone")) +
  coord_cartesian(ylim = c(-2, 4)) +
  labs(x = "Time (ms) relative to target syllable offset",
       y = "Empirical logit of looks to target",
       linetype = "Working Memory",
       axis.title = element_text(size = 17)) +
  theme_grey(base_size = 10, base_family = "Times") + legend_adj_3 +
  theme(text = element_text(size = 17),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 12))





# ----------------------------------------------------------------------------------------------

# Graphs for coda

test_mod_coda <-
  lmer(eLog ~ 1 + (ot1 + ot2 + ot3) * coda_sum * wm_std +
         (1 + coda_sum + (ot1 + ot2 + ot3) | participant) +
         (1 + ot1 + ot2 + ot3 | target),
       control = lmerControl(optimizer = 'bobyqa'), REML = F,
       data = filter(stress_gc_subset, group %in% c("ss")))



# Create design dataframe for predictions
new_dat_wm_minus1 <- stress_gc_subset %>%
  dplyr::select(group, time_zero, ot1:ot3, coda_sum, wm_std) %>%
  distinct(.) %>%
  mutate(wm_std = -1) %>% filter(group %in% c("ss"))
new_dat_wm_0 <- stress_gc_subset %>%
  dplyr::select(group, time_zero, ot1:ot3, coda_sum, wm_std) %>%
  distinct(.) %>%
  mutate(wm_std = 0) %>% filter(group %in% c("ss"))
new_dat_wm_1 <- stress_gc_subset %>%
  dplyr::select(group, time_zero, ot1:ot3, coda_sum, wm_std) %>%
  distinct(.) %>%
  mutate(wm_std = 1) %>% filter(group %in% c("ss"))

# Get model predictions and SE
test_fits_coda <-
  bind_rows(
    predictSE(test_mod_coda, new_dat_wm_minus1) %>%
      as_tibble %>%
      bind_cols(new_dat_wm_minus1),
    predictSE(test_mod_coda, new_dat_wm_0) %>%
      as_tibble %>%
      bind_cols(new_dat_wm_0),
    predictSE(test_mod_coda, new_dat_wm_1) %>%
      as_tibble %>%
      bind_cols(new_dat_wm_1)
  ) %>%
  rename(se = se.fit) %>%
  mutate(ymin = fit - se, ymax = fit + se)




test_fits_coda %>%
  ggplot(., aes(x = time_zero, y = fit, ymax = ymax, ymin = ymin,
                fill = factor(coda_sum), color = factor(coda_sum), lty = factor(wm_std))) +
  facet_grid(. ~ group,
             labeller = labeller(group = group.labs))  +
  geom_hline(yintercept = 0, size = 3, color = "white") +
  geom_vline(xintercept = 4, size = 3, color = "white") +
  geom_line(size = 1) +
  scale_x_continuous(breaks = c(-4, 0, 4, 8, 12),
                     labels = c("-200", "0", "200", "400", "600")) +
  scale_color_brewer(palette = "Set1", name = "Coda", labels = c("CVC", "CV")) +
  coord_cartesian(ylim = c(-2, 4)) +
  labs(x = "Time (ms) relative to target syllable offset",
       y = "Empirical logit of looks to target",
       linetype = "Working Memory") +
  theme_grey(base_size = 10, base_family = "Times") + legend_adj_3 +
  theme(text = element_text(size = 17),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 12))


