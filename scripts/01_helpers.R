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
  str_replace("Conditionneutral", "Neutral Cond.") %>%
  str_replace("Conditionfacilitating", "Facilitating Cond.") %>%
  str_replace("Conditionfiller", "Filler Cond.") %>%
  str_replace("ot(\\d)", "Time^\\1^")  %>%
  str_replace("Time.1.", "Time")  %>%
  str_replace(":", " &times; ") %>%
  str_replace(".Intercept.", "Intercept") %>%
  str_replace("Subj", "Child")

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
  legend.position = c(0.09, 0.94),
  legend.key = element_blank(),
  legend.background = element_blank(),
  strip.background = element_blank(),
  axis.title.y = element_text(size = rel(.9), hjust = 0.95),
  legend.key.size = unit(0.75, 'lines'),
  legend.title = element_text(size = 10),
  plot.margin = unit(rep(2, 4), "mm"),
  panel.grid.major = element_line(colour = 'grey90', size = 0.15),
  panel.grid.minor = element_line(colour = 'grey90', size = 0.15)
)

legend_adj_2 <- theme(
  legend.position = c(0.07, 0.90),
  legend.key = element_blank(),
  legend.background = element_blank(),
  strip.background = element_blank(),
  axis.title.y = element_text(size = rel(.9), hjust = 0.95),
  legend.key.size = unit(0.75, 'lines'),
  legend.title = element_text(size = 10),
  plot.margin = unit(rep(2, 4), "mm"),
  panel.grid.major = element_line(colour = 'grey90', size = 0.15),
  panel.grid.minor = element_line(colour = 'grey90', size = 0.15)
)

# Theme
theme_big <- theme_bw(base_size = 12, base_family = "Times")

# -----------------------------------------------------------------------------
