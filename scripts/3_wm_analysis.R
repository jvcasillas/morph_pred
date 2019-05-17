# Working memory analysis -----------------------------------------------------




# Source libs -----------------------------------------------------------------

source(here::here("scripts", "0_load_libs.R"))

# -----------------------------------------------------------------------------




wm_dataset <- read.spss(here("data", "raw", "gating.sav"), to.data.frame = TRUE)

str(wm_dataset)

wm_df <- wm_dataset %>%
      select(Exp, Group, participant = ExperimentName, WM, Condition) %>%
      filter(., Exp == 'S', Group != 'HS' & Group != 'IN') %>%
      as.tbl(.)

# remove unwanted characters from column
wm_df$participant <- gsub(" ", "", paste(wm_df$participant))

str(wm_df)


# Remove participants according to homogeneity of variance
# test re: working memory
remove <- c("L20", "L21", "L22", "L23", "L30", "L31", "L02", "L05",
            "L06", "L08", "L10", "L15")
wm_df2 <- filter(wm_df, !participant %in% remove)


str(wm_df)
str(wm_df2)


wm_df %>%
  ggplot(., aes(x = Group, y = WM)) +
    geom_boxplot()

wm_df2 %>%
  ggplot(., aes(x = Group, y = WM)) +
    geom_boxplot()

unique(wm_df[wm_df$Group == 'L', 'WM'])
unique(wm_df2[wm_df2$Group == 'L', 'WM'])
