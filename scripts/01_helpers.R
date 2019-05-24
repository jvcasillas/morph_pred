# Helper functions ------------------------------------------------------------

source(here::here("scripts", "00_load_libs.R"))

# -----------------------------------------------------------------------------







# Model helpers ---------------------------------------------------------------

# Workflow adapted from: http://rpubs.com/tjmahr/prettytables_2015

# normal approximation of p-values
normal_approx <- function(ts) {
  2 * (1 - pnorm(abs(ts)))
}

# Create a data-frame with random effect variances and correlations
tidy_ranef_summary <- function(model) {
  vars <- tidy_lme4_variances(model)
  cors <- tidy_lme4_covariances(model) %>% select(-vcov)

  # Create some 1s for the diagonal of the correlation matrix
  self_cor <- vars %>%
    select(-vcov) %>%
    mutate(var2 = var1, sdcor = 1.0) %>%
    na.omit

  # Spread out long-from correlations into a matrix
  cor_matrix <- bind_rows(cors, self_cor) %>%
    mutate(sdcor = fixed_digits(sdcor, 2)) %>%
    spread(var1, sdcor) %>%
    rename(var1 = var2)

  left_join(vars, cor_matrix, by = c("grp", "var1"))
}


tidy_lme4_variances <- . %>%
  VarCorr %>%
  as.data.frame %>%
  filter(is.na(var2)) %>%
  select(-var2)


tidy_lme4_covariances <- . %>%
  VarCorr %>%
  as.data.frame %>%
  filter(!is.na(var2))


# Make a data-frame from an anova object
tidy_anova <- . %>%
  tibble::rownames_to_column("Model") %>%
  tbl_df %>%
  rename(p = `Pr(>Chisq)`, Chi_Df = `Chi Df`)


# Extract fixed effects into a nice data-frame with naive p-values
tidy_lme4 <- . %>%
  tidy("fixed") %>%
  mutate(p = normal_approx(statistic)) %>%
  rename(Parameter = term, Estimate = estimate,
         SE = std.error, t = statistic)


# Rename variables from the main model using markdown (COME BACK TO THIS)
fix_param_names <- . %>%
  str_replace("coda_sum", "Syllable structure") %>%
  str_replace("condition_sum", "Lexical stress") %>%
  str_replace("group", "Group") %>%
  str_replace("Groupint", "Group IN") %>%
  str_replace("Groupla", "Group NIN") %>%
  str_replace("Syllable structure:", "Syllable structure &times; ") %>%
  str_replace("Lexical stress:", "Lexical stress &times; ") %>%
  str_replace("ot(\\d)", "Time^\\1^")  %>%
  str_replace(".Intercept.", "Intercept") %>%
  str_replace(":", " &times; ") %>%
  str_replace("participant", "Participant") %>%
  str_replace("target", "Item")

# Sort random effects groups, and make sure residual comes last
sort_ranef_grps <- function(df) {
  residual <- filter(df, grp == "Residual")
  df %>% filter(grp != "Residual") %>%
    arrange(grp) %>%
    bind_rows(residual)
}

# -----------------------------------------------------------------------------










# Model formatters ------------------------------------------------------------


## Formatters

format_cor <- . %>%
  remove_leading_zero %>%
  leading_minus_sign %>%
  blank_nas


format_fixef_num <- . %>%
  fixed_digits(3) %>%
  leading_minus_sign


# Print with n digits of precision
fixed_digits <- function(xs, n = 2) {
  formatC(xs, digits = n, format = "f")
}


# Print three digits of a p-value, but use the "< .001" notation on tiny values.
format_pval <- function(ps, html = FALSE) {
  small_form <- ifelse(html, "&lt;&nbsp;.001", "< .001")
  ps_chr <- ps %>% fixed_digits(3) %>%
    remove_leading_zero
  ps_chr[ps < 0.001] <- small_form
  ps_chr
}


# Don't print leading zero on bounded numbers.
remove_leading_zero <- function(xs) {
  # Problem if any value is greater than 1.0
  digit_matters <- xs %>% as.numeric %>% abs %>% is_greater_than(1) %>% na.omit
  if (any(digit_matters)) {
    warning("Non-zero leading digit")
  }
  str_replace(xs, "^(-?)0", "\\1")
}

remove_trailing_zero <- function(xs) {
  # Chunk into: optional minus, leading, decimal, trailing, trailing zeros
  str_replace(xs, "(-?)(\\d*)([.])(\\d*)0+$", "\\1\\2\\3\\4")
}


# Use minus instead of hyphen
leading_minus_sign <- . %>% str_replace("^-", "&minus;")

blank_nas <- function(xs) ifelse(is.na(xs), "", xs)

blank_same_as_last <- function(xs) ifelse(is_same_as_last(xs), "", xs)

# Is x[n] the same as x[n-1]
is_same_as_last <- function(xs) {
  same_as_last <- xs == lag(xs)
  # Overwrite NA (first lag) from lag(xs)
  same_as_last[1] <- FALSE
  same_as_last
}

# Break a string into characters
str_tokenize <- . %>% strsplit(NULL) %>% unlist
paste_onto <- function(xs, ys) paste0(ys, xs)

# html tags
parenthesize <- function(xs, parens = c("(", ")")) {
  paste0(parens[1], xs, parens[2])
}

make_html_tagger <- function(tag) {
  template <- paste0("<", tag, ">%s</", tag, ">")
  function(xs) sprintf(template, xs)
}
emphasize <- make_html_tagger("em")
subscript <- make_html_tagger("sub")

# -----------------------------------------------------------------------------



# Printers --------------------------------------------------------------------


## Pretty Printers

# Given an ANOVA model comparison and a model name, give the model comparison
# results in the format "chi-square(x) = y, p = z" (as pretty markdown)
pretty_chi_result <- function(anova_results, model_name) {
  comparison <- anova_results %>%
    tidy_anova %>%
    filter(Model == model_name)

  chi_part <- pretty_chi_df(comparison$Chi_Df, comparison$Chisq)
  p_part <- pretty_p(comparison$p)
  paste0(chi_part, ", ", p_part)
}


# Pretty print chi-square test (i.e., chi^2(dfs) = value)
pretty_chi_df <- function(dfs, chi, html = TRUE) {
  eq <- ifelse(html, "&nbsp;=&nbsp;", " = ")
  sprintf("_&chi;_^2^(%s)%s%.2g", dfs, eq, chi)
}


# Pretty print p-value (for in-line reporting)
pretty_p <- function(p, html = TRUE) {
  space <- ifelse(html, "&nbsp;", " ")
  p_char <- format_pval(p, html) %>% remove_trailing_zero
  p <- ifelse(0.001 <= p, round(p, 3), p)
  # Don't use " = " if formatted p-value would return " < .001"
  p_sep <- ifelse(p < 0.001, space, paste0(space, "=", space))
  sprintf("_p_%s%s", p_sep, p_char)
}

pretty_t <- function(t, html = TRUE) {
  t_formatted <- t %>% fixed_digits(2) %>% leading_minus_sign
  pretty_eq("_t_", t_formatted, html)
}

pretty_se <- function(se, html = TRUE) {
  pretty_eq("SE", fixed_digits(se, 2), html)
}

pretty_eq <- function(lhs, rhs, html = TRUE) {
  eq <- ifelse(html, "&nbsp;=&nbsp;", " = ")
  paste0(lhs, eq, rhs)
}

# Convert a row of fixed effects data-frame to an inline equation sequence:
# Bij = b; SE = se; t = t; p = p
report_fixef_row <- function(df, row) {
  values <- df[row, ]

  # Only B needs special care. Prefix it with an italic gamma and include its
  # subscript.
  values$B %<>% fixed_digits(2) %>% leading_minus_sign
  gamma <- paste0(emphasize("&gamma;"), "~", values$subscript, "~")
  pretty_eq(gamma, values$B)

  equations <- c(pretty_eq(gamma, values$B), pretty_se(values$SE),
                 pretty_t(values$t), pretty_p(values$p))

  # Combine with semicolons
  paste0(equations, collapse = "; ")
}

# logits to proportions
inv_logit <- gtools::inv.logit


# Printing function for reporting models estimates at target offset
target_offset_prob <- function(group, coda, cond) {
  frame <- model_preds$target_offset_preds
  vals <- frame[frame$group == group & frame$coda == coda & frame$cond == cond,
                c("prob", "prob_lb", "prob_ub")]
  vals <- mutate_all(vals, round, digits = 3)
  output <- paste0("Probability = ", vals$prob, "; LB = ",
                   vals$prob_lb, "; UB = ", vals$prob_ub)
  return(output)
}
























# Plots -----------------------------------------------------------------------

# Legends inset
inset_legend <- theme(
  legend.position = c(0.015, 0.98),
  legend.justification = c(0, 1),
  legend.background = element_rect(fill = "white", color = "grey"),
  plot.margin = unit(rep(2, 4), "mm")
)

# Adjustments to legend
legend_adj <- theme(
  legend.position = c(0.09, 0.92),
  legend.key = element_blank(),
  legend.background = element_blank(),
  strip.background = element_blank(),
  axis.title.y = element_text(size = rel(.9), hjust = 0.95),
  axis.title.x = element_text(size = rel(.9)),
  legend.key.size = unit(0.75, 'lines'),
  legend.text = element_text(size = 6),
  legend.title = element_text(size = 7),
  plot.margin = unit(rep(2, 4), "mm"),
  panel.grid.major = element_line(colour = 'grey90', size = 0.15),
  panel.grid.minor = element_line(colour = 'grey90', size = 0.15)
)

legend_adj_2 <- theme(
  legend.position = c(0.05, 0.90),
  legend.key = element_blank(),
  legend.background = element_blank(),
  strip.background = element_blank(),
  axis.title.y = element_text(size = rel(.9), hjust = 0.95),
  axis.title.x = element_text(size = rel(.9)),
  legend.key.size = unit(0.75, 'lines'),
  legend.text = element_text(size = 7),
  legend.title = element_text(size = 8),
  plot.margin = unit(rep(2, 4), "mm"),
  panel.grid.major = element_line(colour = 'grey90', size = 0.15),
  panel.grid.minor = element_line(colour = 'grey90', size = 0.15)
)

# Theme
theme_big <- theme_bw(base_size = 10, base_family = "Times")

# -----------------------------------------------------------------------------
