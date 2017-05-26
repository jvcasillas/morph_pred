#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Morphosyntactic predictability                                              #
# Stress analyses                                                             #
# 11/09/2016                                                                  #
# Script 1                                                                    #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#


# clean working directory
rm(list = ls(all = TRUE))

# Set working directory
# setwd("~/Desktop/morph_pred/")
setwd("~/academia/research/in_progress/morph_pred")


library(tidyverse); library(plotly)
library(lme4); library(lmerTest)


# Read data
# It is labelled as .xls, but the format is comma separated
# stress1 <- read_tsv("./mySources/data/stressOutput_trialBIN5.xls")
# stress2 <- read_tsv("./mySources/data/stressOutput2_trialBIN5.xls")


# Rerun python scripts to make time course longer
stress3 <- read_tsv("./mySources/data/raw/stressBinOutput.txt")
stress4 <- read_tsv("./mySources/data/raw/stressBinOutput2.txt")

# Check structure of each data frame to make sure they 
# can be combined 
# str(stress1); str(stress2)
str(stress3); str(stress4)


# They can, so we combine them 
stress <- rbind(stress3, stress4)


# rename variables
names(stress)[names(stress)=="RECORDING_SESSION_LABEL"] <- "participant"
names(stress)[names(stress)=="TRIAL_LABEL"] <- "trial"
names(stress)[names(stress)=="BIN_INDEX"] <- "bin"
names(stress)[names(stress)=="id"] <- "wavID"
names(stress)[names(stress)=="RIGHT_IA_1_SAMPLE_COUNT"] <- "targetCount"
names(stress)[names(stress)=="RIGHT_IA_2_SAMPLE_COUNT"] <- "distractorCount"
names(stress)[names(stress)=="RIGHT_IA_1_SAMPLE_COUNT_%"] <- "targetProp"
names(stress)[names(stress)=="RIGHT_IA_2_SAMPLE_COUNT_%"] <- "distractorProp"

# remove unnecessary columns
stress <- select(stress, 
                 -trial, 
                 -identifier, 
                 -sentencewav, 
                 -word1, 
                 -word2_20msafterv1, 
                 -word2_c2, 
                 -word2_c3, 
                 -word2_v2, 
                 -word2_c3, 
                 -word4_20msafterv1, 
                 -word4_c1, 
                 -word4_suffix, 
                 -word6, 
                 -word7
                 -`RIGHT_IA_0_SAMPLE_COUNT_%`,
                 -`RIGHT_IA_0_SAMPLE_COUNT_%`)


stress$target <- as.factor(stress$target)

# Set condition variable (1st syll stressed or unstressed)
stress$condition <- 'stressed'
stress[stress$target == 'bebio' | 
       stress$target == 'cambio' | 
       stress$target == 'canto' | 
       stress$target == 'comio' | 
       stress$target == 'compro' | 
       stress$target == 'firmo' | 
       stress$target == 'grabo' | 
       stress$target == 'guardo' | 
       stress$target == 'lanzo' | 
       stress$target == 'lavo' | 
       stress$target == 'lleno' | 
       stress$target == 'mando' | 
       stress$target == 'pinto' | 
       stress$target == 'rompio' | 
       stress$target == 'saco' | 
       stress$target == 'subio', 'condition'] <- 'unstressed'
stress$condition <- as.factor(stress$condition)

# check number of bins for each participant
stress %>%
  select(., participant, bin) %>%
  aggregate(bin ~ participant, FUN = max, data = .) %>% 
  arrange(., bin)

# The minimum bin = 621/373/186/

# Convert proportions (this creates NAs)
stress$targetProp <- gsub(",", ".", paste(stress$targetProp))
stress$targetProp <- as.numeric(stress$targetProp)

stress$distractorProp <- gsub(",", ".", paste(stress$distractorProp))
stress$distractorProp <- as.numeric(stress$distractorProp)






# Create grouping variable 
stress %>%
  separate(., col = participant, into = c('group', 'id'), sep = -3, remove = FALSE) %>% 
  select(., -id) -> stressReduced

stressReduced[stressReduced$group == 'l', 'group'] <- 'lb'
stressReduced[stressReduced$group == 's', 'group'] <- 'ss'

# Change 'in' level of 'group' to 'int'
# because 'in' is a function in R
stressReduced[stressReduced$group == 'in', 'group'] <- 'int'

stressReduced$group <- as.factor(stressReduced$group)


glimpse(stressReduced)
summary(stressReduced$group)

# Add 1 to bin so it starts at 1
df_stress <- stressReduced
df_stress$bin <- df_stress$bin + 1


# remove distractors 
# stressTargets <- c('bebe', 'bebio', 'cambia', 'cambio', 'canta', 'canto', 'come', 'comio', 'compra', 'compro', 'firma', 'firmo', 'gaba', 'grabo', 'guarda', 'guardo', 'lanza', 'lanzo', 'lava', 'lavo', 'llena', 'lleno', 'manda', 'mando', 'pinta', 'pinto', 'rompe', 'rompio', 'saca', 'saco', 'sube', 'subio')
# df_stress <- df_stress[df_stress$target %in% stressTargets, ]
# df_stress <- droplevels(df_stress)


glimpse(df_stress)
summary(df_stress$group)


# Item-specific time normalization
# - create subset of data based on item (for both stressed and unstressed)
# - use item-specific target onset as 0 (in binN column)
#    - this refers to the start of the word prior to the target for each item 
#    - this 'length' is used to insure all of the target word 'info' is included 
#      for analysis (as opposed to truncating some of the critical vowel)
# - save item-specific bin min and max as variables 
# - count forward from target onset to max 
# - count backward from target onset to min
# - 'move' data 200ms for VWP (duplicate binN column + 20)
# - find common min/max between all items (lowest max, highest min)

summary(df_stress$target)

# bebe UNSTRESSED 
bebe <- filter(df_stress, target == 'bebe')
bebe$binN <- 999
tOSbebe <- round((unique(bebe$word3_suffix) - unique(bebe$word2_c1v1)) / 10)
bebe[bebe$bin == tOSbebe, 'binN'] <- 0
bebe$binN <- bebe$bin - tOSbebe

# bebio stressed
bebio <- filter(df_stress, target == 'bebio')
bebio$binN <- 999
tOSbebio <- round((unique(bebio$word3_suffix) - unique(bebio$word2_c1v1)) / 10)
bebio[bebio$bin == tOSbebio, 'binN'] <- 0
bebio$binN <- bebio$bin - tOSbebio

# cambia stressed 
cambia <- filter(df_stress, target == 'cambia')
cambia$binN <- 999
tOScambia <- round((unique(cambia$word3_suffix) - unique(cambia$word2_c1v1)) / 10)
cambia[cambia$bin == tOScambia, 'binN'] <- 0
cambia$binN <- cambia$bin - tOScambia

# cambio  
cambio <- filter(df_stress, target == 'cambio')
cambio$binN <- 999
tOScambio <- round((unique(cambio$word3_suffix) - unique(cambio$word2_c1v1)) / 10)
cambio[cambio$bin == tOScambio, 'binN'] <- 0
cambio$binN <- cambio$bin - tOScambio

# canta  
canta <- filter(df_stress, target == 'canta')
canta$binN <- 999
tOScanta <- round((unique(canta$word3_suffix) - unique(canta$word2_c1v1)) / 10)
canta[canta$bin == tOScanta, 'binN'] <- 0
canta$binN <- canta$bin - tOScanta

# canto   
canto <- filter(df_stress, target == 'canto')
canto$binN <- 999
tOScanto <- round((unique(canto$word3_suffix) - unique(canto$word2_c1v1)) / 10)
canto[canto$bin == tOScanto, 'binN'] <- 0
canto$binN <- canto$bin - tOScanto

# come  
come <- filter(df_stress, target == 'come')
come$binN <- 999
tOScome <- round((unique(come$word3_suffix) - unique(come$word2_c1v1)) / 10)
come[come$bin == tOScome, 'binN'] <- 0
come$binN <- come$bin - tOScome

# comio 
comio <- filter(df_stress, target == 'comio')
comio$binN <- 999
tOScomio <- round((unique(comio$word3_suffix) - unique(comio$word2_c1v1)) / 10)
comio[comio$bin == tOScomio, 'binN'] <- 0
comio$binN <- comio$bin - tOScomio

# compra 
compra <- filter(df_stress, target == 'compra')
compra$binN <- 999
tOScompra <- round((unique(compra$word3_suffix) - unique(compra$word2_c1v1)) / 10)
compra[compra$bin == tOScompra, 'binN'] <- 0
compra$binN <- compra$bin - tOScompra

# compro  
compro <- filter(df_stress, target == 'compro')
compro$binN <- 999
tOScompro <- round((unique(compro$word3_suffix) - unique(compro$word2_c1v1)) / 10)
compro[compro$bin == tOScompro, 'binN'] <- 0
compro$binN <- compro$bin - tOScompro

# firma 
firma <- filter(df_stress, target == 'firma')
firma$binN <- 999
tOSfirma <- round((unique(firma$word3_suffix) - unique(firma$word2_c1v1)) / 10)
firma[firma$bin == tOSfirma, 'binN'] <- 0
firma$binN <- firma$bin - tOSfirma

# firmo   
firmo <- filter(df_stress, target == 'firmo')
firmo$binN <- 999
tOSfirmo <- round((unique(firmo$word3_suffix) - unique(firmo$word2_c1v1)) / 10)
firmo[firmo$bin == tOSfirmo, 'binN'] <- 0
firmo$binN <- firmo$bin - tOSfirmo

# gaba  
gaba <- filter(df_stress, target == 'gaba')
gaba$binN <- 999
tOSgaba <- round((unique(gaba$word3_suffix) - unique(gaba$word2_c1v1)) / 10)
gaba[gaba$bin == tOSgaba, 'binN'] <- 0
gaba$binN <- gaba$bin - tOSgaba

# grabo 
grabo <- filter(df_stress, target == 'grabo')
grabo$binN <- 999
tOSgrabo <- round((unique(grabo$word3_suffix) - unique(grabo$word2_c1v1)) / 10)
grabo[grabo$bin == tOSgrabo, 'binN'] <- 0
grabo$binN <- grabo$bin - tOSgrabo

# guarda 
guarda <- filter(df_stress, target == 'guarda')
guarda$binN <- 999
tOSguarda <- round((unique(guarda$word3_suffix) - unique(guarda$word2_c1v1)) / 10)
guarda[guarda$bin == tOSguarda, 'binN'] <- 0
guarda$binN <- guarda$bin - tOSguarda

# guardo  
guardo <- filter(df_stress, target == 'guardo')
guardo$binN <- 999
tOSguardo <- round((unique(guardo$word3_suffix) - unique(guardo$word2_c1v1)) / 10)
guardo[guardo$bin == tOSguardo, 'binN'] <- 0
guardo$binN <- guardo$bin - tOSguardo

# lanza  
lanza <- filter(df_stress, target == 'lanza')
lanza$binN <- 999
tOSlanza <- round((unique(lanza$word3_suffix) - unique(lanza$word2_c1v1)) / 10)
lanza[lanza$bin == tOSlanza, 'binN'] <- 0
lanza$binN <- lanza$bin - tOSlanza

# lanzo   
lanzo <- filter(df_stress, target == 'lanzo')
lanzo$binN <- 999
tOSlanzo <- round((unique(lanzo$word3_suffix) - unique(lanzo$word2_c1v1)) / 10)
lanzo[lanzo$bin == tOSlanzo, 'binN'] <- 0
lanzo$binN <- lanzo$bin - tOSlanzo

# lava   
lava <- filter(df_stress, target == 'lava')
lava$binN <- 999
tOSlava <- round((unique(lava$word3_suffix) - unique(lava$word2_c1v1)) / 10)
lava[lava$bin == tOSlava, 'binN'] <- 0
lava$binN <- lava$bin - tOSlava

# lavo  
lavo <- filter(df_stress, target == 'lavo')
lavo$binN <- 999
tOSlavo <- round((unique(lavo$word3_suffix) - unique(lavo$word2_c1v1)) / 10)
lavo[lavo$bin == tOSlavo, 'binN'] <- 0
lavo$binN <- lavo$bin - tOSlavo

# llena  
llena <- filter(df_stress, target == 'llena')
llena$binN <- 999
tOSllena <- round((unique(llena$word3_suffix) - unique(llena$word2_c1v1)) / 10)
llena[llena$bin == tOSllena, 'binN'] <- 0
llena$binN <- llena$bin - tOSllena

# lleno 
lleno <- filter(df_stress, target == 'lleno')
lleno$binN <- 999
tOSlleno <- round((unique(lleno$word3_suffix) - unique(lleno$word2_c1v1)) / 10)
lleno[lleno$bin == tOSlleno, 'binN'] <- 0
lleno$binN <- lleno$bin - tOSlleno

# manda  
manda <- filter(df_stress, target == 'manda')
manda$binN <- 999
tOSmanda <- round((unique(manda$word3_suffix) - unique(manda$word2_c1v1)) / 10)
manda[manda$bin == tOSmanda, 'binN'] <- 0
manda$binN <- manda$bin - tOSmanda

# mando  
mando <- filter(df_stress, target == 'mando')
mando$binN <- 999
tOSmando <- round((unique(mando$word3_suffix) - unique(mando$word2_c1v1)) / 10)
mando[mando$bin == tOSmando, 'binN'] <- 0
mando$binN <- mando$bin - tOSmando

# pinta  
pinta <- filter(df_stress, target == 'pinta')
pinta$binN <- 999
tOSpinta <- round((unique(pinta$word3_suffix) - unique(pinta$word2_c1v1)) / 10)
pinta[pinta$bin == tOSpinta, 'binN'] <- 0
pinta$binN <- pinta$bin - tOSpinta

# pinto  
pinto <- filter(df_stress, target == 'pinto')
pinto$binN <- 999
tOSpinto <- round((unique(pinto$word3_suffix) - unique(pinto$word2_c1v1)) / 10)
pinto[pinto$bin == tOSpinto, 'binN'] <- 0
pinto$binN <- pinto$bin - tOSpinto

# rompe 
rompe <- filter(df_stress, target == 'rompe')
rompe$binN <- 999
tOSrompe <- round((unique(rompe$word3_suffix) - unique(rompe$word2_c1v1)) / 10)
rompe[rompe$bin == tOSrompe, 'binN'] <- 0
rompe$binN <- rompe$bin - tOSrompe

# rompio 
rompio <- filter(df_stress, target == 'rompio')
rompio$binN <- 999
tOSrompio <- round((unique(rompio$word3_suffix) - unique(rompio$word2_c1v1)) / 10)
rompio[rompio$bin == tOSrompio, 'binN'] <- 0
rompio$binN <- rompio$bin - tOSrompio

# saca
saca <- filter(df_stress, target == 'saca')
saca$binN <- 999
tOSsaca <- round((unique(saca$word3_suffix) - unique(saca$word2_c1v1)) / 10)
saca[saca$bin == tOSsaca, 'binN'] <- 0
saca$binN <- saca$bin - tOSsaca

# saco
saco <- filter(df_stress, target == 'saco')
saco$binN <- 999
tOSsaco <- round((unique(saco$word3_suffix) - unique(saco$word2_c1v1)) / 10)
saco[saco$bin == tOSsaco, 'binN'] <- 0
saco$binN <- saco$bin - tOSsaco

# sube
sube <- filter(df_stress, target == 'sube')
sube$binN <- 999
tOSsube <- round((unique(sube$word3_suffix) - unique(sube$word2_c1v1)) / 10)
sube[sube$bin == tOSsube, 'binN'] <- 0
sube$binN <- sube$bin - tOSsube

# subio 
subio <- filter(df_stress, target == 'subio')
subio$binN <- 999
tOSsubio <- round((unique(subio$word3_suffix) - unique(subio$word2_c1v1)) / 10)
subio[subio$bin == tOSsubio, 'binN'] <- 0
subio$binN <- subio$bin - tOSsubio

df_adj <- do.call("rbind", list(bebe,
                                bebio,
                                cambia,
                                cambio,
                                canta,
                                canto,
                                come,
                                comio,
                                compra,
                                compro,
                                firma,
                                firmo,
                                gaba,
                                grabo,
                                guarda,
                                guardo,
                                lanza,
                                lanzo,
                                lava,
                                lavo,
                                llena,
                                lleno,
                                manda,
                                mando,
                                pinta,
                                pinto,
                                rompe,
                                rompio,
                                saca,
                                saco,
                                sube,
                                subio))

levels(df_adj$target)[levels(df_adj$target) == "gaba"] <- "graba"

df_adj$binAdj <- df_adj$binN - 20
df_adj <- as.data.frame(df_adj)


# Add ot1 ot2 ot3 to df
# source('./mySources/scripts/stressCreatePoly.R')
# source('./mySources/scripts/stressCreatePoly5.R')
# source('./mySources/scripts/stressCreatePoly10.R')

# check it 
glimpse(df_adj)

# write table
write.table(df_adj, "./mySources/data/clean/stressBIN10Clean.csv", row.names = F, quote = F, sep = ",")


