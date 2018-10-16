# clean working directory
rm(list = ls(all = TRUE))

## @knitr stressLibs

library(plotly); library(tidyverse); library(broom); library(sjPlot)
library(lme4); library(lmerTest); library(gridExtra); library(cowplot)
library(foreign)

## @knitr ignore

# Set working directory
# setwd("~/Desktop/morph_pred/")
setwd("~/academia/research/in_progress/morph_pred")


wm_dataset <- read.spss('./mySources/data/raw/gating.sav', to.data.frame=TRUE)

str(wm_dataset)

wm_df <- wm_dataset %>%
      select(Exp, Group, participant = ExperimentName, WM, Condition) %>% 
      filter(., Exp == 'S', Group != 'HS' & Group != 'IN') %>%
      as.tbl(.)

# remove unwanted characters from column
wm_df$Target <- gsub(" ", "", paste(wm_df$Target))
wm_df$Target <- gsub("\303\263", "o", paste(wm_df$Target))
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