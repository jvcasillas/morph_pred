source_part <- function(file, start, end, ...) {
  file.lines <- scan(file, what = character(),
                     skip = start-1, nlines = end - start + 1, sep = '\n')
  file.lines.collapsed <- paste(file.lines, collapse = '\n')
  source(textConnection(file.lines.collapsed), ...)
}

source_part(here::here(
  "scripts", "stress", "s3_adv_int_nat", "eye_track", "07_analysis_gca_wm.R"),
  start = 18, end = 111)


# Load model
full_mod_lang_learn_wm <- readRDS(here::here(
  "models", "stress", "s3_adv_int_nat", "eye_track", "gca",
  "full_mod_lang_learn_wm.Rds"))

# Save random effects by participant and by target
random_eff <- ranef(full_mod_lang_learn_wm)


# Separate dataframe for participants
random_part <- random_eff$participant

ran_part  <- setNames(cbind(rownames(random_part), random_part, row.names = NULL),
         c("participant", "intercept", "coda_sum",
           "condition_sum", "wm_std", "ot1", "ot2", "ot3"))

# add wm data
mem_data <- read_csv(here("data", "raw", "dur_stress_demographics.csv"))
mem_data$id <- toupper(mem_data$id)
mem_data <- mem_data %>%
  rename(participant = "id") %>%
  filter(!participant %in% c("LA04",
                             "LA06", "LA09", "LA14", "LA15", "LA19",
                             "IN17", "LA07")) %>%
  select(., participant, wm, group)

random_df <- ran_part %>%
left_join(mem_data, by = "participant")

random_df$wm <- as.numeric(random_df$wm)







########


gca_update <- lmer(
  eLog ~ ot1 + ot2 + ot3 + coda_sum + condition_sum +
  group + wm_std + ot1:coda_sum + ot2:coda_sum +
  ot3:coda_sum + ot1:condition_sum + ot2:condition_sum + ot3:condition_sum +
  coda_sum:condition_sum + ot1:group + ot2:group + ot3:group +
  ot1:coda_sum:condition_sum + ot2:coda_sum:condition_sum +
  ot3:coda_sum:condition_sum + coda_sum:condition_sum:group +
  ot1:coda_sum:condition_sum:group + ot2:coda_sum:condition_sum:group +
  ot3:coda_sum:condition_sum:group +
  (1 + coda_sum + condition_sum + wm_std + (ot1 + ot2 + ot3) | participant) +
  (1 + ot1 + ot2 + ot3 | target) ,
  data = stress_gc_subset, REML = F,
  lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 20000))
)

saveRDS(gca_update, here("models", "stress", "s3_adv_int_nat", "eye_track",
                         "gca", "full_mod_lang_learn_wm.Rds"))

summary(gca_update)

point_colors <- c("#5596bd", "#5F6A72", "#CC0033")

mod_re <- ranef(gca_update)$participant[c("wm_std", "ot1", "ot2", "ot3")] %>%
  as_tibble() %>%
  mutate(id = rownames(ranef(gca_update)$participant)) %>%
  separate(col = id, into = c("group", "id"), sep = 2, remove = F) %>%
  pivot_longer(cols = c("ot1", "ot2", "ot3"), names_to = "poly",
               values_to = "estimate") %>%
  mutate(poly = case_when(
    poly == "ot1" ~ "Linear",
    poly == "ot2" ~ "Quadratic",
    TRUE ~ "Cubic"),
    poly = fct_relevel(poly, "Linear", "Quadratic")) %>%
  rename(WM = wm_std)


corr_tidy <- mod_re %>%
  pivot_wider(id_cols = c("group", "id", "WM"),
              names_from = poly, values_from = estimate) %>%
  select(-group, -id) %>%
  cor() %>%
  as.data.frame() %>%
  tibble::rownames_to_column(var = "predictor1") %>%
  pivot_longer(cols = -predictor1, names_to = "predictor2", values_to = "corr")

p_corr_matrix <- corr_tidy %>%
  ggplot(., aes(x = predictor1, y = predictor2, fill = corr)) +
  geom_tile() +
  scale_fill_gradient2(name = "r", limits = c(-1, 1),
                       low = point_colors[3],
                       high = point_colors[1]) +
  labs(y = NULL, x = NULL) +
  theme_bw() +
  #theme(aspect.ratio = 3/4) +
  NULL


corr_pairs <- corr_tidy %>%
  filter(predictor1 != predictor2) %>%
  distinct(predictors = paste(pmax(predictor1, predictor2),
                              pmin(predictor1, predictor2), sep = " ~ "),
           .keep_all = TRUE)

p_corr_bars <- corr_pairs %>%
  ggplot(aes(x = fct_reorder(predictors, corr), y = corr,
             fill = abs(corr) > 0.5)) +
  geom_col(width = 0.75, color = "grey40") +
  coord_flip() +
  ylim(-1, 1) +
  labs(x = NULL, y = "Pearson product-moment correlation (r)") +
  scale_fill_manual(name = "r > abs(0.5)", breaks = c(TRUE, FALSE),
                    labels = c("Yes", "No"), values = point_colors[c(3, 2)]) +
  theme_bw(base_size = 12, base_family = "Times") +
  #theme(aspect.ratio = 3/4) +
  NULL


text_rs <- corr_pairs %>% filter(predictor1 == "WM") %>%
  select(poly = predictor2, corr) %>%
  mutate(text = paste0("r = ", round(corr, digits = 2)),
         poly = fct_relevel(poly, "Linear", "Cubic"))


p_corr_scatter <- mod_re %>%
  ggplot(., aes(x = WM, y = estimate)) +
  facet_wrap(~ poly, scales = "free_y", ncol = 3) +
  geom_vline(xintercept = 0, lty = 2, size = 0.3) +
  geom_hline(yintercept = 0, lty = 2, size = 0.3) +
  geom_smooth(method = lm, se = T, color = "black", fullrange = T,
              formula = "y ~ x", fill = "grey80") +
  geom_point(aes(fill = group, shape = group), size = 1.5) +
  coord_cartesian(expand = T) +
  scale_fill_manual(name = NULL, values = point_colors) +
  scale_shape_manual(name = NULL, values = 21:23) +
  geom_text(data = text_rs, aes(x = -Inf, y = -Inf, label = text),
            hjust   = -4, vjust   = -1) +
  labs(y = "Estimate", x = "β-WM") +
  theme_bw(base_size = 12, base_family = "Times") +
  NULL

View(random_df)

random_df %>%
  ggplot(aes(wm_std, wm)) +
  geom_point(aes(color = group))


p_corr_matrix
p_corr_bars
p_corr_scatter


ggsave("stress_re_p1.png", plot = p_corr_matrix,
       path = here("figs", "stress", "s3_adv_int_nat", "eye_track", "lang_learn"),
       height = 3.5, width = 10)

ggsave("stress_re_p2.png", plot = p_corr_bars,
       path = here("figs", "stress", "s3_adv_int_nat", "eye_track", "lang_learn"),
       height = 3.5, width = 10)

ggsave("stress_re_p3.png", plot = p_corr_scatter,
       path = here("figs", "stress", "s3_adv_int_nat", "eye_track", "lang_learn"),
       heigh = 3.5, width = 10)



corr_pairs %>%
  select(Term = predictors, r = corr) %>%
  mutate_if(is.numeric, round, digits = 2) %>%
  arrange(desc(r)) %>%
  knitr::kable(format = "pandoc")





# Take parameter estimates for wm and and non-linear time terms (represent
# deviations from fixed effects estimates) and plot them.
# There is a correlation between b wm and b ot2.
# Subjs with higher b wm also have higher b ot2
# Higher wm = more curved line
# more curved line equals faster rate of target fixation
