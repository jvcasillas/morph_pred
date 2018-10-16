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


dele_dataset <- read.spss('./mySources/data/raw/gating.sav', to.data.frame=TRUE)

str(dele_dataset)

dele_df <- dele_dataset %>%
      select(Exp, Group, participant = ExperimentName, dele = DELE) %>% 
      filter(., Exp == 'S', Group != 'HS' & Group != 'IN') %>%
      as.tbl(.)

# remove unwanted characters from column
dele_df$participant <- gsub(" ", "", paste(dele_df$participant))

str(dele_df)


# Remove participants according to homogeneity of variance
# test re: working memory
remove <- c("L20", "L21", "L22", "L23", "L30", "L31", "L02", "L05", 
            "L06", "L08", "L10", "L15")
dele_df2 <- filter(dele_df, !participant %in% remove)


str(dele_df)
str(dele_df2)


dele_df %>%
  ggplot(., aes(x = Group, y = dele)) + 
    geom_boxplot()

dele_df2 %>% 
  group_by(., participant, Group) %>% 
  summarize(., dele = unique(dele)) %>%
  ggplot(., aes(x = Group, y = dele)) + 
    geom_boxplot()

unique(dele_df[dele_df$Group == 'L', 'dele'])
unique(dele_df2[dele_df2$Group == 'LA', 'dele'])


dele_df2 %>% 
  group_by(., participant, Group) %>% 
  summarize(., dele = unique(dele)) %>% 
  mutate(., groupSum = ifelse(Group == "L", yes = -0.5, no = 0.5)) %>%
  lm(dele ~ groupSum, data = .) %>% 
  summary 


dele_df2 %>% 
  group_by(., participant, Group) %>% 
  summarize(., dele = unique(dele)) %>% 
  group_by(., Group) %>% 
  summarize(., mean = mean(dele), sd = sd(dele))

#   Group     mean       sd
#  <fctr>    <dbl>    <dbl>
#       L 16.10000 3.956710
#      LA 46.26923 4.094462

# We administered a revised version of the Diploma de Español como 
# Lengua Extranjera (DELE) test (Eng. Diploma of Spanish as a Foreign 
# Language) to all learner participants. The learners had an 
# average DELE score of 31.18 points; however, the distribution was 
# bimodal. Thus we divided the learners into two groups: late beginning 
# learners (LB) and late advanced (LA). Specifically, the LBs’ scores 
# were 15.08 ± 0.76 standard errors (se) lower than the grand mean 
# (x̅ = 16.10, SD = 3.96). The LA group scores were significantly 
# higher by 30.17 points ± 1.51 se (t = 19.98, p < 0.001, x̅ = 46.27, 
# SD = 4.09).