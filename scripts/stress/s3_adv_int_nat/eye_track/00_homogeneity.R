library(tidyverse); library("TOSTER"); library(foreign)



# Load data
dem_all <- read_csv("./data/raw/dur_stress_demographics.csv")
# to get wm from mon group
dataset <- read.spss('./data/raw/gating.sav', to.data.frame=TRUE)



# Remove participant IN17 because we lost eye-tracking data
# (file was corrupted)
# We remove LA09 and LA15 to make the groups homogenous in L2
# proficiency (DELE)

dem_all <- dem_all %>%
  filter(., id != "IN17") %>%
  filter(., id != "LA09" & id != "LA15")

unique(dem_all$id)

dem_all %>%
  group_by(group) %>%
  summarize(n = n_distinct(id))



# Check info: id, group, dele, pstm, wm, age, aoa_l2, months_abroad,
# l1_use_week, l2_use_week, years_word_int

glimpse(dem_all)

# We're missing wm for LA07 (after computer crashed and we lost data)

# Set wm as numeric

dem_all$wm <- as.numeric(dem_all$wm)

# Create a table with mean + sd for: wm, pstm, dele, aoa_l2, months_abroad,
# l1_use_week, l2_use_week, years_word_int

dem_all %>%
  group_by(., group) %>%
  summarise(., wm = round(mean(wm, na.rm=TRUE),2),
               wm_sd = round(sd(wm, na.rm=TRUE),2),
               pstm = round(mean(pstm),2),
               pstm_sd = round(sd(pstm),2),
               dele = round(mean(dele),2),
               dele_sd = round(sd(dele),2),
               n =  n_distinct(id)) %>% knitr::kable()
unique(dem_all$id)

dem_all %>%
  group_by(., group) %>%
  summarise(aoa_l2 = round(mean(aoa_l2),2),
            aoa_l2_sd = round(sd(aoa_l2),2),
            abroad = round(mean(months_abroad),2),
            abroad_sd = round(sd(months_abroad),2),
            l1 = round(mean(l1_use_week),2),
            l1_sd = round(sd(l1_use_week),2),
            l2 = round(mean(l2_use_week),2),
            l2_sd = round(sd(l2_use_week),2),
            n = n_distinct(id)) %>% knitr::kable()


# Some SD values don't work (not sure why), but they work with this function
aggregate(wm ~ group, data = dem_all, FUN = sd, na.rm=TRUE)
aggregate(pstm ~ group, data = dem_all, FUN = sd)
aggregate(dele ~ group, data = dem_all, FUN = sd)
aggregate(aoa_l2 ~ group, data = dem_all, FUN = sd)

  # |group |    wm| wm_sd| pstm| pstm_sd|  dele| dele_sd|  n|
  # |:-----|-----:|-----:|----:|-------:|-----:|-------:|--:|
  # |in    | 10.27|  2.98| 8.59|    1.47| 48.86|    4.32| 22|
  # |la    |  9.00|  2.15| 7.96|    1.24| 46.60|    3.83| 25|

  # |group | aoa_l2| aoa_l2_sd| abroad| abroad_sd|    l1| l1_sd|    l2| l2_sd|  n|
  # |:-----|------:|---------:|------:|---------:|-----:|-----:|-----:|-----:|--:|
  # |in    |  15.14|      4.46|  34.14|     86.99| 64.09| 11.51| 31.59| 14.59| 22|
  # |la    |  13.72|      3.20|  12.68|     15.13| 72.76| 12.91| 27.24| 12.91| 25|


# Get mean and sd for monolingual WM
dataset %>%
  filter(., Group == "S") %>%
  summarise(., wm_mean = mean(WM),
            wm_sd = round(sd(WM),2),
            n = n_distinct(ExperimentName))

#   wm_mean   wm_sd  n
# 1    9.16   1.89  25

## Homogeneity of variances tests
## All seem ok
bartlett.test(wm ~ group, data = dem_all)
bartlett.test(pstm ~ group, data = dem_all)
bartlett.test(dele ~ group, data = dem_all)
bartlett.test(aoa_l2 ~ group, data = dem_all)
bartlett.test(months_abroad ~ group, data = dem_all) # significant, we can't fix this one
bartlett.test(l1_use_week ~ group, data = dem_all)
bartlett.test(l2_use_week ~ group, data = dem_all)

# wm toast la vs in
# all good
TOSTtwo(m1 = 10.27, sd1 = 2.98, n1 = 22, # in
        m2 = 9.00, sd2 = 2.15, n2 = 25, # la
        low_eqbound_d = -0.3,
        high_eqbound_d = 0.3,
        alpha = 0.05)

# t(37.72) = 0.639, p = 0.737

# wm toast mon vs in
# all good
TOSTtwo(m1 = 10.27, sd1 = 2.98, n1 = 22, # in
        m2 = 9.16, sd2 = 1.89, n2 = 25, # mon
        low_eqbound_d = -0.3,
        high_eqbound_d = 0.3,
        alpha = 0.05)

#  t(34.69) = 0.489, p = 0.686

# wm toast mon vs la
# all good
TOSTtwo(m1 = 9.00, sd1 = 2.15, n1 = 25, # la
        m2 = 9.16, sd2 = 1.89, n2 = 25, # mon
        low_eqbound_d = -0.3,
        high_eqbound_d = 0.3,
        alpha = 0.05)

# t(47.22) = 0.781, p = 0.219

###########################################

# pstm toast la vs in
# all good
TOSTtwo(m1 = 8.59, sd1 = 1.47, n1 = 22, # in
        m2 = 7.96, sd2 = 1.24, n2 = 25, # la
        low_eqbound_d = -0.3,
        high_eqbound_d = 0.3,
        alpha = 0.05)

# dele toast la vs in
# all good now that we removed LA09 and LA15
TOSTtwo(m1 = 48.86, sd1 = 4.32, n1 = 22, # in
        m2 = 46.60, sd2 = 3.83, n2 = 25, # la
        low_eqbound_d = -0.3,
        high_eqbound_d = 0.3,
        alpha = 0.05)

# t(42.36) = 1.887, p = 0.0661

# age of acquistion toast la vs in
# all good
TOSTtwo(m1 = 15.14, sd1 = 4.46, n1 = 22, # in
        m2 = 13.72, sd2 = 3.20, n2 = 25, # la
        low_eqbound_d = -0.3,
        high_eqbound_d = 0.3,
        alpha = 0.05)

# time abroad toast la vs in
# all good, but Barlett is not ok, there's a lot of variance in the int group, what do we do here?
TOSTtwo(m1 = 34.14, sd1 = 86.99, n1 = 22, # in
        m2 = 12.68, sd2 = 15.13, n2 = 25, # la
        low_eqbound_d = -0.3,
        high_eqbound_d = 0.3,
        alpha = 0.05)

# l1 use in a normal week toast la vs in
# They are different, Interpreters use less their L1
TOSTtwo(m1 = 64.09, sd1 = 11.51, n1 = 22, # in
        m2 = 72.76, sd2 = 12.91, n2 = 25, # la
        low_eqbound_d = -0.3,
        high_eqbound_d = 0.3,
        alpha = 0.05)

# l2 use in a normal week toast la vs in
# all good, this might be because some interpreters also use other L2s
# for example, the United Nations interpreters speak French and use it more often than Spanish
TOSTtwo(m1 = 31.59, sd1 = 14.59, n1 = 22, # in
        m2 = 27.24, sd2 = 12.91, n2 = 25, # la
        low_eqbound_d = -0.3,
        high_eqbound_d = 0.3,
        alpha = 0.05)
