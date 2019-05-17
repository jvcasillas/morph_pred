#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Morphosyntactic predictability                                              #
# Stress analyses                                                             #
# 11/09/2016                                                                  #
# Script 1                                                                    #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#


# clean working directory
rm(list = ls(all = TRUE))



library(tidyverse); library(plotly);
library(lme4); library(lmerTest)


# Read data and set encoding
# CRIS: I removed the endoding because it was messing up the accents on target words
stress10 <- read_tsv("./mySources/data/raw/stress_10ms.txt") %>%

#stress20 <- read_tsv("./mySources/data/raw/allExpOutput20_ia.txt",
#                     locale = locale(encoding = "LATIN1")) %>%

#stress50 <- read_tsv("./mySources/data/raw/stress_50ms.txt") %>%

  # Filter out unwanted experiments
  filter(., exp == "stress") %>%

  # Remove unnecessary columns
  # CRIS: I changed the RIGHT_IA_... to AVERAGE_IA_...
  # (new data frame uses cyclopean mode from data viewer)
  select(., -TRIAL_LABEL,                 -TRIAL_INDEX,
            -identifier,                  -sentencewav,
            -word1,                       -word2_20msafterv1,
            -word2_c2,                    -word2_c3,
            -word2_v2,                    -word2_c3,
            -word4_20msafterv1,           -word4_c1,
            -word4_suffix,                -word6,
            -word7,                       -EYE_TRACKED,
            -IA_0_ID,                     -IA_3_ID,
            -IA_4_ID,                     -AVERAGE_IA_3_SAMPLE_COUNT,
            -AVERAGE_IA_4_SAMPLE_COUNT,     -AVERAGE_IA_0_SAMPLE_COUNT,
            -`AVERAGE_IA_0_SAMPLE_COUNT_%`, -`AVERAGE_IA_3_SAMPLE_COUNT_%`,
            -`AVERAGE_IA_4_SAMPLE_COUNT_%`) %>%

  # Rename some columns
  rename(., participant = RECORDING_SESSION_LABEL,
            bin = BIN_INDEX,
            wavID = id,
            targetCount = AVERAGE_IA_1_SAMPLE_COUNT,
            distractorCount = AVERAGE_IA_2_SAMPLE_COUNT,
            targetProp = `AVERAGE_IA_1_SAMPLE_COUNT_%`,
            distractorProp = `AVERAGE_IA_2_SAMPLE_COUNT_%`) %>%

  # Remove weird chars from target list
  # Create condition variable (unstressed/stressed)
  # Create coda variable (0 = no coda, 1 = coda)
  # Create corr variable (0 = incorrect, 1 = correct)
  # Change ',' to '.' in proportion columns
  # Create eLog variable and respective wts
  mutate(., target = as.factor(target),
            target = gsub("\u0097", "o", paste(.$target)),
            condition = ifelse(target %in%
              c('bebió', 'cambió', 'cantó', 'comió', 'compró',
                'firmó', 'grabó', 'guardó', 'lanzó', 'lavó',
                'llenó', 'mandó', 'pintó', 'rompió', 'sacó',
                'subió'),
                yes = "unstressed", no = "stressed"),
            coda = ifelse(target %in%
              c('bebe', 'bebió', 'llena', 'llenó', 'sube', 'subió',
                'come', 'comió', 'saca', 'sacó', 'lava', 'lavó',
                'graba', 'grabó'),
                yes = 0, no = 1),
            corr = ifelse(correctresponse == "Lshift" & KEY_PRESSED == "Lshift" |
                          correctresponse == "Rshift" & KEY_PRESSED == "Rshift",
                          yes = 1, no = 0),
            targetProp = as.numeric(gsub(",", ".", paste(.$targetProp))),
            distractorProp = as.numeric(gsub(",", ".", paste(.$distractorProp))),
            eLog = log((targetCount + 0.5) / (50 - targetCount + 0.5)),
            wts = 1 / (targetCount + 0.5) + 1 / (50 - targetCount + 0.5)) %>%

  # Create 'group' column and new 'id'
  # in order to match participant ids with WM df
  separate(., col = participant, into = c("group", "id"), sep = -2, remove = TRUE) %>%

  # Recode groups that have random labels
  # create new participant id labels
  # Realign bins to start from 1
  mutate(., group = recode(group, l = 'lb', s = 'ss', `in` = 'int'),
            group_temp = recode(group, lb = 'L', ss = "SB",ss = "SBs", int = "IN", hs = "HS", la = "LA"),
            target = recode(target, gaba = 'graba'),
            bin = bin + 1) %>%

  # Combine group_temp and id to
  # complete recode of participant IDs (now consistent with WM)
  unite(., col = participant, group_temp, id, sep = "", remove = TRUE)


glimpse(stress10)
#unique(stress10$group)
#glimpse(stress20)
#glimpse(stress50)
















# Item-specific time normalization
# - create subset of data based on item (for both stressed and unstressed)
# - use item-specific target onset as 0 (in binN column)
#    - this refers to the start of the target for each item
#    - this 'length' is used to insure all of the target word 'info' is included
#      for analysis (as opposed to truncating some of the critical vowel)
# - save item-specific bin min and max as variables


 df_stress <- stress10
# df_stress <- stress20
# df_stress <- stress50

summary(as.factor(df_stress$target))

# bebe UNSTRESSED
bebe <- filter(df_stress, target == 'bebe')
bebe$binTonsetAlign <- 999
tOSbebe <- round(unique(bebe$word3_c1v1) / 10)
bebe[bebe$bin == tOSbebe, 'binTonsetAlign'] <- 0
bebe$binTonsetAlign <- bebe$bin - tOSbebe
bebe$binTsuffixAlign <- 999
tSSbebe <- round(unique(bebe$word3_suffix) / 10)
bebe[bebe$bin == tSSbebe, 'binTsuffixAlign'] <- 0
bebe$binTsuffixAlign <- bebe$bin - tSSbebe

# bebio stressed
bebio <- filter(df_stress, target == 'bebió')
bebio$binTonsetAlign <- 999
tOSbebio <- round(unique(bebio$word3_c1v1) / 10)
bebio[bebio$bin == tOSbebio, 'binTonsetAlign'] <- 0
bebio$binTonsetAlign <- bebio$bin - tOSbebio
bebio$binTsuffixAlign <- 999
tSSbebio <- round(unique(bebio$word3_suffix) / 10)
bebio[bebio$bin == tSSbebio, 'binTsuffixAlign'] <- 0
bebio$binTsuffixAlign <- bebio$bin - tSSbebio

# cambia stressed
cambia <- filter(df_stress, target == 'cambia')
cambia$binTonsetAlign <- 999
tOScambia <- round(unique(cambia$word3_c1v1) / 10)
cambia[cambia$bin == tOScambia, 'binTonsetAlign'] <- 0
cambia$binTonsetAlign <- cambia$bin - tOScambia
cambia$binTsuffixAlign <- 999
tSScambia <- round(unique(cambia$word3_suffix) / 10)
cambia[cambia$bin == tSScambia, 'binTsuffixAlign'] <- 0
cambia$binTsuffixAlign <- cambia$bin - tSScambia

# cambio
cambio <- filter(df_stress, target == 'cambió')
cambio$binTonsetAlign <- 999
tOScambio <- round(unique(cambio$word3_c1v1) / 10)
cambio[cambio$bin == tOScambio, 'binTonsetAlign'] <- 0
cambio$binTonsetAlign <- cambio$bin - tOScambio
cambio$binTsuffixAlign <- 999
tSScambio <- round(unique(cambio$word3_suffix) / 10)
cambio[cambio$bin == tSScambio, 'binTsuffixAlign'] <- 0
cambio$binTsuffixAlign <- cambio$bin - tSScambio

# canta
canta <- filter(df_stress, target == 'canta')
canta$binTonsetAlign <- 999
tOScanta <- round(unique(canta$word3_c1v1) / 10)
canta[canta$bin == tOScanta, 'binTonsetAlign'] <- 0
canta$binTonsetAlign <- canta$bin - tOScanta
canta$binTsuffixAlign <- 999
tSScanta <- round(unique(canta$word3_suffix) / 10)
canta[canta$bin == tSScanta, 'binTsuffixAlign'] <- 0
canta$binTsuffixAlign <- canta$bin - tSScanta

# canto
canto <- filter(df_stress, target == 'cantó')
canto$binTonsetAlign <- 999
tOScanto <- round(unique(canto$word3_c1v1) / 10)
canto[canto$bin == tOScanto, 'binTonsetAlign'] <- 0
canto$binTonsetAlign <- canto$bin - tOScanto
canto$binTsuffixAlign <- 999
tSScanto <- round(unique(canto$word3_suffix) / 10)
canto[canto$bin == tSScanto, 'binTsuffixAlign'] <- 0
canto$binTsuffixAlign <- canto$bin - tSScanto

# come
come <- filter(df_stress, target == 'come')
come$binTonsetAlign <- 999
tOScome <- round(unique(come$word3_c1v1) / 10)
come[come$bin == tOScome, 'binTonsetAlign'] <- 0
come$binTonsetAlign <- come$bin - tOScome
come$binTsuffixAlign <- 999
tSScome <- round(unique(come$word3_suffix) / 10)
come[come$bin == tSScome, 'binTsuffixAlign'] <- 0
come$binTsuffixAlign <- come$bin - tSScome

# comio
comio <- filter(df_stress, target == 'comió')
comio$binTonsetAlign <- 999
tOScomio <- round(unique(comio$word3_c1v1) / 10)
comio[comio$bin == tOScomio, 'binTonsetAlign'] <- 0
comio$binTonsetAlign <- comio$bin - tOScomio
comio$binTsuffixAlign <- 999
tSScomio <- round(unique(comio$word3_suffix) / 10)
comio[comio$bin == tSScomio, 'binTsuffixAlign'] <- 0
comio$binTsuffixAlign <- comio$bin - tSScomio

# compra
compra <- filter(df_stress, target == 'compra')
compra$binTonsetAlign <- 999
tOScompra <- round(unique(compra$word3_c1v1) / 10)
compra[compra$bin == tOScompra, 'binTonsetAlign'] <- 0
compra$binTonsetAlign <- compra$bin - tOScompra
compra$binTsuffixAlign <- 999
tSScompra <- round(unique(compra$word3_suffix) / 10)
compra[compra$bin == tSScompra, 'binTsuffixAlign'] <- 0
compra$binTsuffixAlign <- compra$bin - tSScompra

# compro
compro <- filter(df_stress, target == 'compró')
compro$binTonsetAlign <- 999
tOScompro <- round(unique(compro$word3_c1v1) / 10)
compro[compro$bin == tOScompro, 'binTonsetAlign'] <- 0
compro$binTonsetAlign <- compro$bin - tOScompro
compro$binTsuffixAlign <- 999
tSScompro <- round(unique(compro$word3_suffix) / 10)
compro[compro$bin == tSScompro, 'binTsuffixAlign'] <- 0
compro$binTsuffixAlign <- compro$bin - tSScompro

# firma
firma <- filter(df_stress, target == 'firma')
firma$binTonsetAlign <- 999
tOSfirma <- round(unique(firma$word3_c1v1) / 10)
firma[firma$bin == tOSfirma, 'binTonsetAlign'] <- 0
firma$binTonsetAlign <- firma$bin - tOSfirma
firma$binTsuffixAlign <- 999
tSSfirma <- round(unique(firma$word3_suffix) / 10)
firma[firma$bin == tSSfirma, 'binTsuffixAlign'] <- 0
firma$binTsuffixAlign <- firma$bin - tSSfirma

# firmo
firmo <- filter(df_stress, target == 'firmó')
firmo$binTonsetAlign <- 999
tOSfirmo <- round(unique(firmo$word3_c1v1) / 10)
firmo[firmo$bin == tOSfirmo, 'binTonsetAlign'] <- 0
firmo$binTonsetAlign <- firmo$bin - tOSfirmo
firmo$binTsuffixAlign <- 999
tSSfirmo <- round(unique(firmo$word3_suffix) / 10)
firmo[firmo$bin == tSSfirmo, 'binTsuffixAlign'] <- 0
firmo$binTsuffixAlign <- firmo$bin - tSSfirmo

# graba
graba <- filter(df_stress, target == 'graba')
graba$binTonsetAlign <- 999
tOSgraba <- round(unique(graba$word3_c1v1) / 10)
graba[graba$bin == tOSgraba, 'binTonsetAlign'] <- 0
graba$binTonsetAlign <- graba$bin - tOSgraba
graba$binTsuffixAlign <- 999
tSSgraba <- round(unique(graba$word3_suffix) / 10)
graba[graba$bin == tSSgraba, 'binTsuffixAlign'] <- 0
graba$binTsuffixAlign <- graba$bin - tSSgraba

# grabo
grabo <- filter(df_stress, target == 'grabó')
grabo$binTonsetAlign <- 999
tOSgrabo <- round(unique(grabo$word3_c1v1) / 10)
grabo[grabo$bin == tOSgrabo, 'binTonsetAlign'] <- 0
grabo$binTonsetAlign <- grabo$bin - tOSgrabo
grabo$binTsuffixAlign <- 999
tSSgrabo <- round(unique(grabo$word3_suffix) / 10)
grabo[grabo$bin == tSSgrabo, 'binTsuffixAlign'] <- 0
grabo$binTsuffixAlign <- grabo$bin - tSSgrabo

# guarda
guarda <- filter(df_stress, target == 'guarda')
guarda$binTonsetAlign <- 999
tOSguarda <- round(unique(guarda$word3_c1v1) / 10)
guarda[guarda$bin == tOSguarda, 'binTonsetAlign'] <- 0
guarda$binTonsetAlign <- guarda$bin - tOSguarda
guarda$binTsuffixAlign <- 999
tSSguarda <- round(unique(guarda$word3_suffix) / 10)
guarda[guarda$bin == tSSguarda, 'binTsuffixAlign'] <- 0
guarda$binTsuffixAlign <- guarda$bin - tSSguarda

# guardo
guardo <- filter(df_stress, target == 'guardó')
guardo$binTonsetAlign <- 999
tOSguardo <- round(unique(guardo$word3_c1v1) / 10)
guardo[guardo$bin == tOSguardo, 'binTonsetAlign'] <- 0
guardo$binTonsetAlign <- guardo$bin - tOSguardo
guardo$binTsuffixAlign <- 999
tSSguardo <- round(unique(guardo$word3_suffix) / 10)
guardo[guardo$bin == tSSguardo, 'binTsuffixAlign'] <- 0
guardo$binTsuffixAlign <- guardo$bin - tSSguardo

# lanza
lanza <- filter(df_stress, target == 'lanza')
lanza$binTonsetAlign <- 999
tOSlanza <- round(unique(lanza$word3_c1v1) / 10)
lanza[lanza$bin == tOSlanza, 'binTonsetAlign'] <- 0
lanza$binTonsetAlign <- lanza$bin - tOSlanza
lanza$binTsuffixAlign <- 999
tSSlanza <- round(unique(lanza$word3_suffix) / 10)
lanza[lanza$bin == tSSlanza, 'binTsuffixAlign'] <- 0
lanza$binTsuffixAlign <- lanza$bin - tSSlanza

# lanzo
lanzo <- filter(df_stress, target == 'lanzó')
lanzo$binTonsetAlign <- 999
tOSlanzo <- round(unique(lanzo$word3_c1v1) / 10)
lanzo[lanzo$bin == tOSlanzo, 'binTonsetAlign'] <- 0
lanzo$binTonsetAlign <- lanzo$bin - tOSlanzo
lanzo$binTsuffixAlign <- 999
tSSlanzo <- round(unique(lanzo$word3_suffix) / 10)
lanzo[lanzo$bin == tSSlanzo, 'binTsuffixAlign'] <- 0
lanzo$binTsuffixAlign <- lanzo$bin - tSSlanzo

# lava
lava <- filter(df_stress, target == 'lava')
lava$binTonsetAlign <- 999
tOSlava <- round(unique(lava$word3_c1v1) / 10)
lava[lava$bin == tOSlava, 'binTonsetAlign'] <- 0
lava$binTonsetAlign <- lava$bin - tOSlava
lava$binTsuffixAlign <- 999
tSSlava <- round(unique(lava$word3_suffix) / 10)
lava[lava$bin == tSSlava, 'binTsuffixAlign'] <- 0
lava$binTsuffixAlign <- lava$bin - tSSlava

# lavo
lavo <- filter(df_stress, target == 'lavó')
lavo$binTonsetAlign <- 999
tOSlavo <- round(unique(lavo$word3_c1v1) / 10)
lavo[lavo$bin == tOSlavo, 'binTonsetAlign'] <- 0
lavo$binTonsetAlign <- lavo$bin - tOSlavo
lavo$binTsuffixAlign <- 999
tSSlavo <- round(unique(lavo$word3_suffix) / 10)
lavo[lavo$bin == tSSlavo, 'binTsuffixAlign'] <- 0
lavo$binTsuffixAlign <- lavo$bin - tSSlavo

# llena
llena <- filter(df_stress, target == 'llena')
llena$binTonsetAlign <- 999
tOSllena <- round(unique(llena$word3_c1v1) / 10)
llena[llena$bin == tOSllena, 'binTonsetAlign'] <- 0
llena$binTonsetAlign <- llena$bin - tOSllena
llena$binTsuffixAlign <- 999
tSSllena <- round(unique(llena$word3_suffix) / 10)
llena[llena$bin == tSSllena, 'binTsuffixAlign'] <- 0
llena$binTsuffixAlign <- llena$bin - tSSllena

# lleno
lleno <- filter(df_stress, target == 'llenó')
lleno$binTonsetAlign <- 999
tOSlleno <- round(unique(lleno$word3_c1v1) / 10)
lleno[lleno$bin == tOSlleno, 'binTonsetAlign'] <- 0
lleno$binTonsetAlign <- lleno$bin - tOSlleno
lleno$binTsuffixAlign <- 999
tSSlleno <- round(unique(lleno$word3_suffix) / 10)
lleno[lleno$bin == tSSlleno, 'binTsuffixAlign'] <- 0
lleno$binTsuffixAlign <- lleno$bin - tSSlleno

# manda
manda <- filter(df_stress, target == 'manda')
manda$binTonsetAlign <- 999
tOSmanda <- round(unique(manda$word3_c1v1) / 10)
manda[manda$bin == tOSmanda, 'binTonsetAlign'] <- 0
manda$binTonsetAlign <- manda$bin - tOSmanda
manda$binTsuffixAlign <- 999
tSSmanda <- round(unique(manda$word3_suffix) / 10)
manda[manda$bin == tSSmanda, 'binTsuffixAlign'] <- 0
manda$binTsuffixAlign <- manda$bin - tSSmanda

# mando
mando <- filter(df_stress, target == 'mandó')
mando$binTonsetAlign <- 999
tOSmando <- round(unique(mando$word3_c1v1) / 10)
mando[mando$bin == tOSmando, 'binTonsetAlign'] <- 0
mando$binTonsetAlign <- mando$bin - tOSmando
mando$binTsuffixAlign <- 999
tSSmando <- round(unique(mando$word3_suffix) / 10)
mando[mando$bin == tSSmando, 'binTsuffixAlign'] <- 0
mando$binTsuffixAlign <- mando$bin - tSSmando

# pinta
pinta <- filter(df_stress, target == 'pinta')
pinta$binTonsetAlign <- 999
tOSpinta <- round(unique(pinta$word3_c1v1) / 10)
pinta[pinta$bin == tOSpinta, 'binTonsetAlign'] <- 0
pinta$binTonsetAlign <- pinta$bin - tOSpinta
pinta$binTsuffixAlign <- 999
tSSpinta <- round(unique(pinta$word3_suffix) / 10)
pinta[pinta$bin == tSSpinta, 'binTsuffixAlign'] <- 0
pinta$binTsuffixAlign <- pinta$bin - tSSpinta

# pinto
pinto <- filter(df_stress, target == 'pintó')
pinto$binTonsetAlign <- 999
tOSpinto <- round(unique(pinto$word3_c1v1) / 10)
pinto[pinto$bin == tOSpinto, 'binTonsetAlign'] <- 0
pinto$binTonsetAlign <- pinto$bin - tOSpinto
pinto$binTsuffixAlign <- 999
tSSpinto <- round(unique(pinto$word3_suffix) / 10)
pinto[pinto$bin == tSSpinto, 'binTsuffixAlign'] <- 0
pinto$binTsuffixAlign <- pinto$bin - tSSpinto

# rompe
rompe <- filter(df_stress, target == 'rompe')
rompe$binTonsetAlign <- 999
tOSrompe <- round(unique(rompe$word3_c1v1) / 10)
rompe[rompe$bin == tOSrompe, 'binTonsetAlign'] <- 0
rompe$binTonsetAlign <- rompe$bin - tOSrompe
rompe$binTsuffixAlign <- 999
tSSrompe <- round(unique(rompe$word3_suffix) / 10)
rompe[rompe$bin == tSSrompe, 'binTsuffixAlign'] <- 0
rompe$binTsuffixAlign <- rompe$bin - tSSrompe

# rompio
rompio <- filter(df_stress, target == 'rompió')
rompio$binTonsetAlign <- 999
tOSrompio <- round(unique(rompio$word3_c1v1) / 10)
rompio[rompio$bin == tOSrompio, 'binTonsetAlign'] <- 0
rompio$binTonsetAlign <- rompio$bin - tOSrompio
rompio$binTsuffixAlign <- 999
tSSrompio <- round(unique(rompio$word3_suffix) / 10)
rompio[rompio$bin == tSSrompio, 'binTsuffixAlign'] <- 0
rompio$binTsuffixAlign <- rompio$bin - tSSrompio

# saca
saca <- filter(df_stress, target == 'saca')
saca$binTonsetAlign <- 999
tOSsaca <- round(unique(saca$word3_c1v1) / 10)
saca[saca$bin == tOSsaca, 'binTonsetAlign'] <- 0
saca$binTonsetAlign <- saca$bin - tOSsaca
saca$binTsuffixAlign <- 999
tSSsaca <- round(unique(saca$word3_suffix) / 10)
saca[saca$bin == tSSsaca, 'binTsuffixAlign'] <- 0
saca$binTsuffixAlign <- saca$bin - tSSsaca

# saco
saco <- filter(df_stress, target == 'sacó')
saco$binTonsetAlign <- 999
tOSsaco <- round(unique(saco$word3_c1v1) / 10)
saco[saco$bin == tOSsaco, 'binTonsetAlign'] <- 0
saco$binTonsetAlign <- saco$bin - tOSsaco
saco$binTsuffixAlign <- 999
tSSsaco <- round(unique(saco$word3_suffix) / 10)
saco[saco$bin == tSSsaco, 'binTsuffixAlign'] <- 0
saco$binTsuffixAlign <- saco$bin - tSSsaco

# sube
sube <- filter(df_stress, target == 'sube')
sube$binTonsetAlign <- 999
tOSsube <- round(unique(sube$word3_c1v1) / 10)
sube[sube$bin == tOSsube, 'binTonsetAlign'] <- 0
sube$binTonsetAlign <- sube$bin - tOSsube
sube$binTsuffixAlign <- 999
tSSsube <- round(unique(sube$word3_suffix) / 10)
sube[sube$bin == tSSsube, 'binTsuffixAlign'] <- 0
sube$binTsuffixAlign <- sube$bin - tSSsube

# subio
subio <- filter(df_stress, target == 'subió')
subio$binTonsetAlign <- 999
tOSsubio <- round(unique(subio$word3_c1v1) / 10)
subio[subio$bin == tOSsubio, 'binTonsetAlign'] <- 0
subio$binTonsetAlign <- subio$bin - tOSsubio
subio$binTsuffixAlign <- 999
tSSsubio <- round(unique(subio$word3_suffix) / 10)
subio[subio$bin == tSSsubio, 'binTsuffixAlign'] <- 0
subio$binTsuffixAlign <- subio$bin - tSSsubio

df_stress_adj <- do.call("rbind", list(bebe,
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
                                graba,
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

df_stress_adj$binAdj <- df_stress_adj$binTonsetAlign - 20
#df_stress_adj$binAdj <- df_stress_adj$binTonsetAlign - 10
#df_stress_adj$binAdj <- df_stress_adj$binTonsetAlign - 4

df_stress_adj <- as.data.frame(df_stress_adj)



# check it
glimpse(df_stress_adj)

# write table
 write.table(df_stress_adj, "./mySources/data/clean/stressBIN10iaClean_20190513.csv", row.names = F, quote = F, sep = ",")
# write.table(df_stress_adj, "./mySources/data/clean/stressBIN20iaClean.csv", row.names = F, quote = F, sep = ",")
#write.table(df_stress_adj, "./mySources/data/clean/stressBIN50iaClean_20190513.csv", row.names = F, quote = F, sep = ",")


