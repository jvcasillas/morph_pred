
# Graph for HLS 2019, showing probability of looks towards the target
# averaging over all conditions at the offset of the target syllable

stress_prob <- model_preds$target_offset_preds %>%
  mutate(coda = if_else(coda == 1, "CV", "CVC"),
         cond = if_else(cond == 1, "paroxytone", "oxytone")) %>%
  ggplot(aes(x = group, y = prob, color = cond, shape = coda)) +
  geom_point(position = position_dodge(0.5), size = 3) +
  geom_hline(yintercept = 0.5, color = 'white', size = 3) +
  geom_errorbar(aes(ymin = prob_lb, ymax = prob_ub), width = 0.2,
                position = position_dodge(0.5)) +
  labs(y = "prob looks to target", color = "condition")

figs_path <- here("figs", "stress", "s3_adv_int_nat", "eye_track")

ggsave(paste0(figs_path, "/stress_prob.png"), stress_prob, width = 150,
       height = 120, units = "mm", dpi = 600)
