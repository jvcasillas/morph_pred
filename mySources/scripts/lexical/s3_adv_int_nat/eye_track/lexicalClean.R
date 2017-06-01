#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Morphosyntactic predictability                                              #
# Duration clean                                                              #
# 05/26/2017                                                                  #
# Study 3 (interpreters)                                                      #
# Script 1                                                                    #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#


# clean working directory
rm(list = ls(all = TRUE))

# Set working directory
# setwd("~/Desktop/morph_pred/")
setwd("~/academia/research/in_progress/morph_pred")


library(tidyverse); library(lme4); library(lmerTest)


# Read data
#lex1_temp <- read_tsv("./mySources/data/raw/lexicalBinOutput.txt")
#lex2_temp <- read_tsv("./mySources/data/raw/lexicalBinOutput2.txt")

lex3_temp <- read_tsv("./mySources/data/raw/allExpOutput.txt")


# Check structure of each data frame to make sure they 
# can be combined 
# str(lex1_temp); str(lex2_temp)

glimpse(lex3_temp)




# They can, so we combine them 
#lex_temp <- rbind(lex1_temp, lex2_temp)

lex_temp <- lex3_temp %>% filter(., exp == "lexical")

glimpse(lex_temp)


# rename variables
names(lex_temp)[names(lex_temp)=="RECORDING_SESSION_LABEL"] <- "participant"
names(lex_temp)[names(lex_temp)=="TRIAL_LABEL"] <- "trial"
names(lex_temp)[names(lex_temp)=="BIN_INDEX"] <- "bin"
names(lex_temp)[names(lex_temp)=="id"] <- "wavID"
names(lex_temp)[names(lex_temp)=="RIGHT_IA_1_SAMPLE_COUNT"] <- "targetCount"
names(lex_temp)[names(lex_temp)=="RIGHT_IA_2_SAMPLE_COUNT"] <- "distractorCount"
names(lex_temp)[names(lex_temp)=="RIGHT_IA_1_SAMPLE_COUNT_%"] <- "targetProp"
names(lex_temp)[names(lex_temp)=="RIGHT_IA_2_SAMPLE_COUNT_%"] <- "distractorProp"

# remove unnecessary columns
lex_temp <- select(lex_temp, 
                   -trial, 
                   -identifier, 
                   -sentencewav, 
                   -word2_v2, 
                   -word3_20msafterv1,
                   -word3_c2, 
                   -word3_c3, 
                   -word3_suffix, 
                   -word4_20msafterv1,
                   -word4_c1,
                   -word4_suffix, 
                   -word6, 
                   -word7,
                   -IA_3_ID,
                   -IA_4_ID,
                   -IA_0_ID,
                   -RIGHT_IA_3_SAMPLE_COUNT,
                   -RIGHT_IA_4_SAMPLE_COUNT,
                   -RIGHT_IA_0_SAMPLE_COUNT,
                   -`RIGHT_IA_3_SAMPLE_COUNT_%`,
                   -`RIGHT_IA_4_SAMPLE_COUNT_%`,
                   -`RIGHT_IA_0_SAMPLE_COUNT_%`)

# Create condition variable (monosyllabic or bisyllabic)
lex_temp$target <- as.factor(lex_temp$target)
lex_bi_syl <- c('marco', 'salsa', 'golpe', 'busto', 'barco', 'gasto', 'parque', 'tosta')
lex_temp$condition <- 'monosyllabic'
lex_temp[lex_temp$target %in% lex_bi_syl, 'condition'] <- 'bisyllabic'


# check number of bins for each participant
# so we can remove fluff
lex_temp %>%
  select(., participant, bin) %>%
  aggregate(bin ~ participant, FUN = max, data = .) %>% 
  arrange(., bin)
# There are a lot with huge bin numbers



# Create grouping variable 
# we use the 'participant' column. All the ids have 2 numbers, but they dont 
# all have 2 letters, so we count from the right
lexReduced <- lex_temp %>%
  separate(., col = participant, into = c('group', 'id'), sep = -3, remove = FALSE) %>% 
  select(., -id)




# Some levels are erroneously labelled 's' or 'l'
lexReduced[lexReduced$group == 'l', 'group'] <- 'lb'
lexReduced[lexReduced$group == 's', 'group'] <- 'ss'

# Change 'in' level of 'group' to 'int'
# because 'in' is a function in R
lexReduced[lexReduced$group == 'in', 'group'] <- 'int'

lexReduced$group <- as.factor(lexReduced$group)

# Add 1 to bin so it starts at 1
lexReduced$bin <- lexReduced$bin + 1


summary(lexReduced$group)
glimpse(lexReduced)





# Add new column to indicate the iteration of the target (1 or 2)
lexReduced$condToken <- 999
lexReduced[lexReduced$target == "mar"    & lexReduced$verb == "tener", 'condToken'] <- 1
lexReduced[lexReduced$target == "marco"  & lexReduced$verb == "tener", 'condToken'] <- 1
lexReduced[lexReduced$target == "sal"    & lexReduced$verb == "costar", 'condToken'] <- 1
lexReduced[lexReduced$target == "salsa"  & lexReduced$verb == "costar", 'condToken'] <- 1
lexReduced[lexReduced$target == "gol"    & lexReduced$verb == "salir", 'condToken'] <- 1
lexReduced[lexReduced$target == "golpe"  & lexReduced$verb == "salir", 'condToken'] <- 1
lexReduced[lexReduced$target == "bus"    & lexReduced$verb == "volver", 'condToken'] <- 1
lexReduced[lexReduced$target == "busto"  & lexReduced$verb == "volver", 'condToken'] <- 1
lexReduced[lexReduced$target == "barco"  & lexReduced$verb == "cerrar", 'condToken'] <- 1
lexReduced[lexReduced$target == "bar"    & lexReduced$verb == "cerrar", 'condToken'] <- 1
lexReduced[lexReduced$target == "gasto"  & lexReduced$verb == "resultar", 'condToken'] <- 1
lexReduced[lexReduced$target == "gas"    & lexReduced$verb == "resultar", 'condToken'] <- 1
lexReduced[lexReduced$target == "parque" & lexReduced$verb == "producir", 'condToken'] <- 1
lexReduced[lexReduced$target == "par"    & lexReduced$verb == "producir", 'condToken'] <- 1
lexReduced[lexReduced$target == "tosta"  & lexReduced$verb == "causar", 'condToken'] <- 1
lexReduced[lexReduced$target == "tos"    & lexReduced$verb == "causar", 'condToken'] <- 1

lexReduced[lexReduced$target == "mar"    & lexReduced$verb == "resultar", 'condToken'] <- 2 # deberia ser parecer
lexReduced[lexReduced$target == "marco"  & lexReduced$verb == "resultar", 'condToken'] <- 2 # deberia ser parecer
lexReduced[lexReduced$target == "sal"    & lexReduced$verb == "proporcionar", 'condToken'] <- 2
lexReduced[lexReduced$target == "salsa"  & lexReduced$verb == "proporcionar", 'condToken'] <- 2
lexReduced[lexReduced$target == "gol"    & lexReduced$verb == "conllevar", 'condToken'] <- 2
lexReduced[lexReduced$target == "golpe"  & lexReduced$verb == "conllevar", 'condToken'] <- 2
lexReduced[lexReduced$target == "bus"    & lexReduced$verb == "salir", 'condToken'] <- 2 # deberia ser venir
lexReduced[lexReduced$target == "busto"  & lexReduced$verb == "salir", 'condToken'] <- 2 # deberia ser venir
lexReduced[lexReduced$target == "barco"  & lexReduced$verb == "vender", 'condToken'] <- 2
lexReduced[lexReduced$target == "bar"    & lexReduced$verb == "vender", 'condToken'] <- 2
lexReduced[lexReduced$target == "gasto"  & lexReduced$verb == "dar", 'condToken'] <- 2
lexReduced[lexReduced$target == "gas"    & lexReduced$verb == "dar", 'condToken'] <- 2
lexReduced[lexReduced$target == "parque" & lexReduced$verb == "vender", 'condToken'] <- 2
lexReduced[lexReduced$target == "par"    & lexReduced$verb == "vender", 'condToken'] <- 2
lexReduced[lexReduced$target == "tosta"  & lexReduced$verb == "generar", 'condToken'] <- 2
lexReduced[lexReduced$target == "tos"    & lexReduced$verb == "generar", 'condToken'] <- 2

# Some of the verbs in the dataframe are not the verbs the 
# used in the stimuli. We change them here
lexReduced[lexReduced$verb == 'resultar' & 
           lexReduced$condToken == 2 & 
           lexReduced$target == "mar", 'verb'] <- "parecer"
lexReduced[lexReduced$verb == 'resultar' & 
           lexReduced$condToken == 2 & 
           lexReduced$target == "marco", 'verb'] <- "parecer"

lexReduced[lexReduced$verb == 'salir' & 
           lexReduced$condToken == 2 & 
           lexReduced$target == "bus", 'verb'] <- "venir"
lexReduced[lexReduced$verb == 'salir' & 
           lexReduced$condToken == 2 & 
           lexReduced$target == "busto", 'verb'] <- "venir"

xtabs(~ target + verb, lexReduced)



# We create a new column, 'targetoffset' where we put the time point in which 
# each target syllable ends
# We have the start and end of the syllable... at some point we will get the start
# and end of the vowel
# word2_c1v1 = m (start)
#            =   (start vowel)
#            =   (end vowel)
# word3_c1v1 = r (end)


# if the taget is mono then the end of the 1st syllable is word3_c1v1
# if the target is bi  then the end of the 1st syllable is word2_c3

lexReduced$targetOffset <- 99999999
lexReduced[lexReduced$condition == 'monosyllabic', 'targetOffset'] <- lexReduced[lexReduced$condition == 'monosyllabic', 'word3_c1v1']
lexReduced[lexReduced$condition == 'bisyllabic', 'targetOffset'] <- lexReduced[lexReduced$condition == 'bisyllabic', 'word2_c3']



df_lex <- lexReduced

df_lex <- droplevels(df_lex)

glimpse(df_lex)


# Item-specific time normalization
# - create subset of data based on item (for both monosyllabic and disyllabic)
# - use item-specific target onset as 0 (in binN column)
#    - this refers to the start of the word prior to the target for each item 
#    - this 'length' is used to insure all of the target word 'info' is included 
#      for analysis (as opposed to truncating some of the critical vowel)
# - save item-specific bin min and max as variables 
# - count forward from target onset to max 
# - count backward from target onset to min
# - 'move' data 200ms for VWP (duplicate binN column + 20)
# - find common min/max between all items (lowest max, highest min)
# - must distinguish between iterations of same words (different preceeding verb)

xtabs(~ target + verb, df_lex)

df_lex %>%
  select(., target, verb, word2_c1v1, targetOffset, word4_c1v1) %>%
  group_by(., target, verb) %>%
  summarise(., start = round(unique(word2_c1v1)),  
              tOS = round(unique(targetOffset)), 
              word4_start = round(unique(word4_c1v1))) %>%
  as.data.frame(.)

# Change target offset for 'gas/resultar'
df_lex[df_lex$target == "gas" & df_lex$verb == "resultar", 'targetOffset'] <- 793

# bar : mono : Cerrar
barCerrar <- filter(df_lex, target == 'bar' & verb == 'cerrar')
barCerrar$binN <- 999
tOSbarCerrar <- round(unique(barCerrar$word2_c1v1) / 10)
barCerrar[barCerrar$bin == tOSbarCerrar, 'binN'] <- 0
barCerrar$binN <- barCerrar$bin - tOSbarCerrar

# bar : mono : Vender
barVender <- filter(df_lex, target == 'bar' & verb == 'vender')
barVender$binN <- 999
tOSbarVender <- round(unique(barVender$word2_c1v1) / 10)
barVender[barVender$bin == tOSbarVender, 'binN'] <- 0
barVender$binN <- barVender$bin - tOSbarVender

# barco : di : Cerrar
barcoCerrar <- filter(df_lex, target == 'barco' & verb == 'cerrar')
barcoCerrar$binN <- 999
tOSbarcoCerrar <- round(unique(barcoCerrar$word2_c1v1) / 10)
barcoCerrar[barcoCerrar$bin == tOSbarcoCerrar, 'binN'] <- 0
barcoCerrar$binN <- barcoCerrar$bin - tOSbarcoCerrar

# barco : di : Vender
barcoVender <- filter(df_lex, target == 'barco' & verb == 'vender')
barcoVender$binN <- 999
tOSbarcoVender <- round(unique(barcoVender$word2_c1v1) / 10)
barcoVender[barcoVender$bin == tOSbarcoVender, 'binN'] <- 0
barcoVender$binN <- barcoVender$bin - tOSbarcoVender




# bus : mono : Volver
busVolver <- filter(df_lex, target == 'bus' & verb == 'volver')
busVolver$binN <- 999
tOSbusVolver <- round(unique(busVolver$word2_c1v1) / 10)
busVolver[busVolver$bin == tOSbusVolver, 'binN'] <- 0
busVolver$binN <- busVolver$bin - tOSbusVolver

# bus : mono : Venir
busVenir <- filter(df_lex, target == 'bus' & verb == 'venir')
busVenir$binN <- 999
tOSbusVenir <- round(unique(busVenir$word2_c1v1) / 10)
busVenir[busVenir$bin == tOSbusVenir, 'binN'] <- 0
busVenir$binN <- busVenir$bin - tOSbusVenir

# busto : di : Volver
bustoVolver <- filter(df_lex, target == 'busto' & verb == 'volver')
bustoVolver$binN <- 999
tOSbustoVolver <- round(unique(bustoVolver$word2_c1v1) / 10)
bustoVolver[bustoVolver$bin == tOSbustoVolver, 'binN'] <- 0
bustoVolver$binN <- bustoVolver$bin - tOSbustoVolver

# busto : di : Venir
bustoVenir <- filter(df_lex, target == 'busto' & verb == 'venir')
bustoVenir$binN <- 999
tOSbustoVenir <- round(unique(bustoVenir$word2_c1v1) / 10)
bustoVenir[bustoVenir$bin == tOSbustoVenir, 'binN'] <- 0
bustoVenir$binN <- bustoVenir$bin - tOSbustoVenir




# gas : mono : Resultar
gasResultar <- filter(df_lex, target == 'gas' & verb == 'resultar')
gasResultar$binN <- 999
tOSgasResultar <- round(unique(gasResultar$word2_c1v1) / 10)
gasResultar[gasResultar$bin == tOSgasResultar, 'binN'] <- 0
gasResultar$binN <- gasResultar$bin - tOSgasResultar

# gas : mono : Dar
gasDar <- filter(df_lex, target == 'gas' & verb == 'dar')
gasDar$binN <- 999
tOSgasDar <- round(unique(gasDar$word2_c1v1) / 10)
gasDar[gasDar$bin == tOSgasDar, 'binN'] <- 0
gasDar$binN <- gasDar$bin - tOSgasDar

# gasto : di : Resultar
gastoResultar <- filter(df_lex, target == 'gasto' & verb == 'resultar')
gastoResultar$binN <- 999
tOSgastoResultar <- round(unique(gastoResultar$word2_c1v1) / 10)
gastoResultar[gastoResultar$bin == tOSgastoResultar, 'binN'] <- 0
gastoResultar$binN <- gastoResultar$bin - tOSgastoResultar

# gasto : di : Dar
gastoDar <- filter(df_lex, target == 'gasto' & verb == 'dar')
gastoDar$binN <- 999
tOSgastoDar <- round(unique(gastoDar$word2_c1v1) / 10)
gastoDar[gastoDar$bin == tOSgastoDar, 'binN'] <- 0
gastoDar$binN <- gastoDar$bin - tOSgastoDar





# gol : mono : Salir
golSalir <- filter(df_lex, target == 'gol' & verb == 'salir')
golSalir$binN <- 999
tOSgolSalir <- round(unique(golSalir$word2_c1v1) / 10)
golSalir[golSalir$bin == tOSgolSalir, 'binN'] <- 0
golSalir$binN <- golSalir$bin - tOSgolSalir

# gol : mono : Conllevar
golConllevar <- filter(df_lex, target == 'gol' & verb == 'conllevar')
golConllevar$binN <- 999
tOSgolConllevar <- round(unique(golConllevar$word2_c1v1) / 10)
golConllevar[golConllevar$bin == tOSgolConllevar, 'binN'] <- 0
golConllevar$binN <- golConllevar$bin - tOSgolConllevar

# golpe : di : Salir
golpeSalir <- filter(df_lex, target == 'golpe' & verb == 'salir')
golpeSalir$binN <- 999
tOSgolpeSalir <- round(unique(golpeSalir$word2_c1v1) / 10)
golpeSalir[golpeSalir$bin == tOSgolpeSalir, 'binN'] <- 0
golpeSalir$binN <- golpeSalir$bin - tOSgolpeSalir

# golpe : di : Conllevar
golpeConllevar <- filter(df_lex, target == 'golpe' & verb == 'conllevar')
golpeConllevar$binN <- 999
tOSgolpeConllevar <- round(unique(golpeConllevar$word2_c1v1) / 10)
golpeConllevar[golpeConllevar$bin == tOSgolpeConllevar, 'binN'] <- 0
golpeConllevar$binN <- golpeConllevar$bin - tOSgolpeConllevar




# mar : mono : Tener
marTener <- filter(df_lex, target == 'mar' & verb == 'tener')
marTener$binN <- 999
tOSmarTener <- round(unique(marTener$word2_c1v1) / 10)
marTener[marTener$bin == tOSmarTener, 'binN'] <- 0
marTener$binN <- marTener$bin - tOSmarTener

# mar : mono : Parecer
marParecer <- filter(df_lex, target == 'mar' & verb == 'parecer')
marParecer$binN <- 999
tOSmarParecer <- round(unique(marParecer$word2_c1v1) / 10)
marParecer[marParecer$bin == tOSmarParecer, 'binN'] <- 0
marParecer$binN <- marParecer$bin - tOSmarParecer

# marco : di : Tener
marcoTener <- filter(df_lex, target == 'marco' & verb == 'tener')
marcoTener$binN <- 999
tOSmarcoTener <- round(unique(marcoTener$word2_c1v1) / 10)
marcoTener[marcoTener$bin == tOSmarcoTener, 'binN'] <- 0
marcoTener$binN <- marcoTener$bin - tOSmarcoTener

# marco : di : Parecer
marcoParecer <- filter(df_lex, target == 'marco' & verb == 'parecer')
marcoParecer$binN <- 999
tOSmarcoParecer <- round(unique(marcoParecer$word2_c1v1) / 10)
marcoParecer[marcoParecer$bin == tOSmarcoParecer, 'binN'] <- 0
marcoParecer$binN <- marcoParecer$bin - tOSmarcoParecer





# par : mono : Producir
parProducir <- filter(df_lex, target == 'par' & verb == 'producir')
parProducir$binN <- 999
tOSparProducir <- round(unique(parProducir$word2_c1v1) / 10)
parProducir[parProducir$bin == tOSparProducir, 'binN'] <- 0
parProducir$binN <- parProducir$bin - tOSparProducir

# par : mono : Vender
parVender <- filter(df_lex, target == 'par' & verb == 'vender')
parVender$binN <- 999
tOSparVender <- round(unique(parVender$word2_c1v1) / 10)
parVender[parVender$bin == tOSparVender, 'binN'] <- 0
parVender$binN <- parVender$bin - tOSparVender

# parque : di : Producir
parqueProducir <- filter(df_lex, target == 'parque' & verb == 'producir')
parqueProducir$binN <- 999
tOSparqueProducir <- round(unique(parqueProducir$word2_c1v1) / 10)
parqueProducir[parqueProducir$bin == tOSparqueProducir, 'binN'] <- 0
parqueProducir$binN <- parqueProducir$bin - tOSparqueProducir

# parque : di : Vender
parqueVender <- filter(df_lex, target == 'parque' & verb == 'vender')
parqueVender$binN <- 999
tOSparqueVender <- round(unique(parqueVender$word2_c1v1) / 10)
parqueVender[parqueVender$bin == tOSparqueVender, 'binN'] <- 0
parqueVender$binN <- parqueVender$bin - tOSparqueVender






# sal : mono : Costar
salCostar <- filter(df_lex, target == 'sal' & verb == 'costar')
salCostar$binN <- 999
tOSsalCostar <- round(unique(salCostar$word2_c1v1) / 10)
salCostar[salCostar$bin == tOSsalCostar, 'binN'] <- 0
salCostar$binN <- salCostar$bin - tOSsalCostar

# sal : mono : Proporcionar
salProporcionar <- filter(df_lex, target == 'sal' & verb == 'proporcionar')
salProporcionar$binN <- 999
tOSsalProporcionar <- round(unique(salProporcionar$word2_c1v1) / 10)
salProporcionar[salProporcionar$bin == tOSsalProporcionar, 'binN'] <- 0
salProporcionar$binN <- salProporcionar$bin - tOSsalProporcionar

# salsa : di : Costar
salsaCostar <- filter(df_lex, target == 'salsa' & verb == 'costar')
salsaCostar$binN <- 999
tOSsalsaCostar <- round(unique(salsaCostar$word2_c1v1) / 10)
salsaCostar[salsaCostar$bin == tOSsalsaCostar, 'binN'] <- 0
salsaCostar$binN <- salsaCostar$bin - tOSsalsaCostar

# salsa : di : Proporcionar
salsaProporcionar <- filter(df_lex, target == 'salsa' & verb == 'proporcionar')
salsaProporcionar$binN <- 999
tOSsalsaProporcionar <- round(unique(salsaProporcionar$word2_c1v1) / 10)
salsaProporcionar[salsaProporcionar$bin == tOSsalsaProporcionar, 'binN'] <- 0
salsaProporcionar$binN <- salsaProporcionar$bin - tOSsalsaProporcionar




# tos : mono : Causar
tosCausar <- filter(df_lex, target == 'tos' & verb == 'causar')
tosCausar$binN <- 999
tOStosCausar <- round(unique(tosCausar$word2_c1v1) / 10)
tosCausar[tosCausar$bin == tOStosCausar, 'binN'] <- 0
tosCausar$binN <- tosCausar$bin - tOStosCausar

# tos : mono : Generar
tosGenerar <- filter(df_lex, target == 'tos' & verb == 'generar')
tosGenerar$binN <- 999
tOStosGenerar <- round(unique(tosGenerar$word2_c1v1) / 10)
tosGenerar[tosGenerar$bin == tOStosGenerar, 'binN'] <- 0
tosGenerar$binN <- tosGenerar$bin - tOStosGenerar

# tosta : di : Causar
tostaCausar <- filter(df_lex, target == 'tosta' & verb == 'causar')
tostaCausar$binN <- 999
tOStostaCausar <- round(unique(tostaCausar$word2_c1v1) / 10)
tostaCausar[tostaCausar$bin == tOStostaCausar, 'binN'] <- 0
tostaCausar$binN <- tostaCausar$bin - tOStostaCausar

# tosta : di : Generar
tostaGenerar <- filter(df_lex, target == 'tosta' & verb == 'generar')
tostaGenerar$binN <- 999
tOStostaGenerar <- round(unique(tostaGenerar$word2_c1v1) / 10)
tostaGenerar[tostaGenerar$bin == tOStostaGenerar, 'binN'] <- 0
tostaGenerar$binN <- tostaGenerar$bin - tOStostaGenerar









# 'rbind' all the dataframes together into 1
df_adj <- do.call("rbind", list(barCerrar, 
                                barVender, 
                                barcoCerrar, 
                                barcoVender, 
                                busVolver, 
                                busVenir, 
                                bustoVolver, 
                                bustoVenir, 
                                gasResultar, 
                                gasDar, 
                                gastoDar, 
                                gastoResultar, 
                                golSalir, 
                                golConllevar, 
                                golpeSalir, 
                                golpeConllevar, 
                                marTener, 
                                marParecer, 
                                marcoTener, 
                                marcoParecer, 
                                parProducir, 
                                parVender, 
                                parqueProducir, 
                                parqueVender, 
                                salCostar, 
                                salProporcionar, 
                                salsaCostar, 
                                salsaProporcionar, 
                                tosCausar, 
                                tosGenerar, 
                                tostaCausar, 
                                tostaGenerar))




# 200ms bin adjustment for VWP
df_adj$binAdj <- df_adj$binN - 20
df_adj <- as.data.frame(df_adj)



# check it 
glimpse(df_adj)

# write table
write.table(df_adj, "./mySources/data/clean/lexicalBIN10CleanNEW.csv", row.names = F, quote = T, sep = ",")

