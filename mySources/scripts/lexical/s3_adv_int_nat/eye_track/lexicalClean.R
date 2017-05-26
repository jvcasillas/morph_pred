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
lex1_temp <- read_tsv("./mySources/data/raw/lexicalBinOutput.txt")
lex2_temp <- read_tsv("./mySources/data/raw/lexicalBinOutput2.txt")

# Check structure of each data frame to make sure they 
# can be combined 
str(lex1_temp); str(lex2_temp)


# They can, so we combine them 
lex_temp <- rbind(lex1_temp, lex2_temp)


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


# bar : mono : cerrar
barCerrar <- filter(df_lex, target == 'bar' & verb == 'cerrar')
barCerrar$binN <- 999
tOSbarCerrar <- round(unique(barCerrar$targetOffset) / 20)
# tOSbarCerrar <- round((unique(barCerrar$targetOffset) - unique(barCerrar$word3_c1v1)) / 20)
barCerrar[barCerrar$bin == tOSbarCerrar, 'binN'] <- 0
barCerrar$binN <- barCerrar$bin - tOSbarCerrar

# bar : mono : Vender
barVender <- filter(df_lex, target == 'bar' & verb == 'vender')
barVender$binN <- 999
tOSbarVender <- round(unique(barVender$targetOffset) / 20)
# tOSbarVender <- round((unique(barVender$targetOffset) - unique(barVender$word3_c1v1)) / 20)
barVender[barVender$bin == tOSbarVender, 'binN'] <- 0
barVender$binN <- barVender$bin - tOSbarVender

# barco : di : Cerrar
barcoCerrar <- filter(df_lex, target == 'barco' & verb == 'cerrar')
barcoCerrar$binN <- 999
tOSbarcoCerrar <- round(unique(barcoCerrar$targetOffset) / 20)
# tOSbarcoCerrar <- round((unique(barcoCerrar$targetOffset) - unique(barcoCerrar$word3_c1v1)) / 20)
barcoCerrar[barcoCerrar$bin == tOSbarcoCerrar, 'binN'] <- 0
barcoCerrar$binN <- barcoCerrar$bin - tOSbarcoCerrar

# barco : di : Vender
barcoVender <- filter(df_lex, target == 'barco' & verb == 'vender')
barcoVender$binN <- 999
tOSbarcoVender <- round(unique(barcoVender$targetOffset) / 20)
# tOSbarcoVender <- round((unique(barcoVender$targetOffset) - unique(barcoVender$word3_c1v1)) / 20)
barcoVender[barcoVender$bin == tOSbarcoVender, 'binN'] <- 0
barcoVender$binN <- barcoVender$bin - tOSbarcoVender




# bus : mono : Venir
busVenir <- filter(df_lex, target == 'bus' & verb == 'venir')
busVenir$binN <- 999
# tOSbusVenir <- round(unique(busVenir$targetOffset) / 20)
tOSbusVenir <- round((unique(busVenir$targetOffset) - unique(busVenir$word3_c1v1)) / 20)
busVenir[busVenir$bin == tOSbusVenir, 'binN'] <- 0
busVenir$binN <- busVenir$bin - tOSbusVenir

# bus : mono : Volver
busVolver <- filter(df_lex, target == 'bus' & verb == 'volver')
busVolver$binN <- 999
tOSbusVolver <- round(unique(busVolver$targetOffset) / 20)
# tOSbusVolver <- round((unique(busVolver$targetOffset) - unique(busVolver$word3_c1v1)) / 20)
busVolver[busVolver$bin == tOSbusVolver, 'binN'] <- 0
busVolver$binN <- busVolver$bin - tOSbusVolver

# busto : di : Venir
bustoVenir <- filter(df_lex, target == 'busto' & verb == 'venir')
bustoVenir$binN <- 999
tOSbustoVenir <- round(unique(bustoVenir$targetOffset) / 20)
# tOSbustoVenir <- round((unique(bustoVenir$targetOffset) - unique(bustoVenir$word3_c1v1)) / 20)
bustoVenir[bustoVenir$bin == tOSbustoVenir, 'binN'] <- 0
bustoVenir$binN <- bustoVenir$bin - tOSbustoVenir

# busto : di : Volver
bustoVolver <- filter(df_lex, target == 'busto' & verb == 'volver')
bustoVolver$binN <- 999
tOSbustoVolver <- round(unique(bustoVolver$targetOffset) / 20)
# tOSbustoVolver <- round((unique(bustoVolver$targetOffset) - unique(bustoVolver$word3_c1v1)) / 20)
bustoVolver[bustoVolver$bin == tOSbustoVolver, 'binN'] <- 0
bustoVolver$binN <- bustoVolver$bin - tOSbustoVolver






# gas : mono : Dar
gasDar <- filter(df_lex, target == 'gas' & verb == 'dar')
gasDar$binN <- 999
tOSgasDar <- round(unique(gasDar$targetOffset) / 20)
# tOSgasDar <- round((unique(gasDar$targetOffset) - unique(gasDar$word3_c1v1)) / 20)
gasDar[gasDar$bin == tOSgasDar, 'binN'] <- 0
gasDar$binN <- gasDar$bin - tOSgasDar

# gas : mono : Resultar
gasResultar <- filter(df_lex, target == 'gas' & verb == 'resultar')
gasResultar$binN <- 999
tOSgasResultar <- round(unique(gasResultar$targetOffset) / 20)
# tOSgasResultar <- round((unique(gasResultar$targetOffset) - unique(gasResultar$word3_c1v1)) / 20)
gasResultar[gasResultar$bin == tOSgasResultar, 'binN'] <- 0
gasResultar$binN <- gasResultar$bin - tOSgasResultar

# gasto : di : Dar
gastoDar <- filter(df_lex, target == 'gasto' & verb == 'dar')
gastoDar$binN <- 999
tOSgastoDar <- round(unique(gastoDar$targetOffset) / 20)
# tOSgastoDar <- round((unique(gastoDar$targetOffset) - unique(gastoDar$word3_c1v1)) / 20)
gastoDar[gastoDar$bin == tOSgastoDar, 'binN'] <- 0
gastoDar$binN <- gastoDar$bin - tOSgastoDar

# gasto : di : Resultar
gastoResultar <- filter(df_lex, target == 'gasto' & verb == 'resultar')
gastoResultar$binN <- 999
tOSgastoResultar <- round(unique(gastoResultar$targetOffset) / 20)
# tOSgastoResultar <- round((unique(gastoResultar$targetOffset) - unique(gastoResultar$word3_c1v1)) / 20)
gastoResultar[gastoResultar$bin == tOSgastoResultar, 'binN'] <- 0
gastoResultar$binN <- gastoResultar$bin - tOSgastoResultar


# STOPPED HERE



# mes : mono : Exigir
mesExigir <- filter(df_lex, target == 'mes' & verb == 'exigir')
mesExigir$binN <- 999
# tOSmesExigir <- round(unique(mesExigir$targetOffset) / 20)
tOSmesExigir <- round((unique(mesExigir$targetOffset) - unique(mesExigir$word3_c1v1)) / 20)
mesExigir[mesExigir$bin == tOSmesExigir, 'binN'] <- 0
mesExigir$binN <- mesExigir$bin - tOSmesExigir
mesExigir$condToken <- 1

# mes : mono : Conseguir
mesConseguir <- filter(df_lex, target == 'mes' & verb == 'conseguir')
mesConseguir$binN <- 999
# tOSmesConseguir <- round(unique(mesConseguir$targetOffset) / 20)
tOSmesConseguir <- round((unique(mesConseguir$targetOffset) - unique(mesConseguir$word3_c1v1)) / 20)
mesConseguir[mesConseguir$bin == tOSmesConseguir, 'binN'] <- 0
mesConseguir$binN <- mesConseguir$bin - tOSmesConseguir
mesConseguir$condToken <- 2

# meses : di : Exigir
mesesExigir <- filter(df_lex, target == 'meses' & verb == 'exigir')
mesesExigir$binN <- 999
# tOSmesesExigir <- round(unique(mesesExigir$targetOffset) / 20)
tOSmesesExigir <- round((unique(mesesExigir$targetOffset) - unique(mesesExigir$word3_c1v1)) / 20)
mesesExigir[mesesExigir$bin == tOSmesesExigir, 'binN'] <- 0
mesesExigir$binN <- mesesExigir$bin - tOSmesesExigir
mesesExigir$condToken <- 1

# meses : di : Conseguir
mesesConseguir <- filter(df_lex, target == 'meses' & verb == 'conseguir')
mesesConseguir$binN <- 999
# tOSmesesConseguir <- round(unique(mesesConseguir$targetOffset) / 20)
tOSmesesConseguir <- round((unique(mesesConseguir$targetOffset) - unique(mesesConseguir$word3_c1v1)) / 20)
mesesConseguir[mesesConseguir$bin == tOSmesesConseguir, 'binN'] <- 0
mesesConseguir$binN <- mesesConseguir$bin - tOSmesesConseguir
mesesConseguir$condToken <- 2


# miel : mono : Hacer
mielHacer <- filter(df_lex, target == 'miel' & verb == 'hacer')
mielHacer$binN <- 999
# tOSmielHacer <- round(unique(mielHacer$targetOffset) / 20)
tOSmielHacer <- round((unique(mielHacer$targetOffset) - unique(mielHacer$word3_c1v1)) / 20)
mielHacer[mielHacer$bin == tOSmielHacer, 'binN'] <- 0
mielHacer$binN <- mielHacer$bin - tOSmielHacer
mielHacer$condToken <- 1

# miel : mono : Transportar
mielTransportar <- filter(df_lex, target == 'miel' & verb == 'transportar')
mielTransportar$binN <- 999
# tOSmielTransportar <- round(unique(mielTransportar$targetOffset) / 20)
tOSmielTransportar <- round((unique(mielTransportar$targetOffset) - unique(mielTransportar$word3_c1v1)) / 20)
mielTransportar[mielTransportar$bin == tOSmielTransportar, 'binN'] <- 0
mielTransportar$binN <- mielTransportar$bin - tOSmielTransportar
mielTransportar$condToken <- 2

# mieles : di : Hacer
mielesHacer <- filter(df_lex, target == 'mieles' & verb == 'hacer')
mielesHacer$binN <- 999
# tOSmielesHacer <- round(unique(mielesHacer$targetOffset) / 20)
tOSmielesHacer <- round((unique(mielesHacer$targetOffset) - unique(mielesHacer$word3_c1v1)) / 20)
mielesHacer[mielesHacer$bin == tOSmielesHacer, 'binN'] <- 0
mielesHacer$binN <- mielesHacer$bin - tOSmielesHacer
mielesHacer$condToken <- 1

# mieles : di : Transportar
mielesTransportar <- filter(df_lex, target == 'mieles' & verb == 'transportar')
mielesTransportar$binN <- 999
# tOSmielesTransportar <- round(unique(mielesTransportar$targetOffset) / 20)
tOSmielesTransportar <- round((unique(mielesTransportar$targetOffset) - unique(mielesTransportar$word3_c1v1)) / 20)
mielesTransportar[mielesTransportar$bin == tOSmielesTransportar, 'binN'] <- 0
mielesTransportar$binN <- mielesTransportar$bin - tOSmielesTransportar
mielesTransportar$condToken <- 2



# rol : mono : Estrenar
rolEstrenar <- filter(df_lex, target == 'rol' & verb == 'estrenar')
rolEstrenar$binN <- 999
# tOSrolEstrenar <- round(unique(rolEstrenar$targetOffset) / 20)
tOSrolEstrenar <- round((unique(rolEstrenar$targetOffset) - unique(rolEstrenar$word3_c1v1)) / 20)
rolEstrenar[rolEstrenar$bin == tOSrolEstrenar, 'binN'] <- 0
rolEstrenar$binN <- rolEstrenar$bin - tOSrolEstrenar
rolEstrenar$condToken <- 1

# rol : mono : Aceptar
rolAceptar <- filter(df_lex, target == 'rol' & verb == 'aceptar')
rolAceptar$binN <- 999
# tOSrolAceptar <- round(unique(rolAceptar$targetOffset) / 20)
tOSrolAceptar <- round((unique(rolAceptar$targetOffset) - unique(rolAceptar$word3_c1v1)) / 20)
rolAceptar[rolAceptar$bin == tOSrolAceptar, 'binN'] <- 0
rolAceptar$binN <- rolAceptar$bin - tOSrolAceptar
rolAceptar$condToken <- 2

# roles : di : Estrenar
rolesEstrenar <- filter(df_lex, target == 'roles' & verb == 'estrenar')
rolesEstrenar$binN <- 999
# tOSrolesEstrenar <- round(unique(rolesEstrenar$targetOffset) / 20)
tOSrolesEstrenar <- round((unique(rolesEstrenar$targetOffset) - unique(rolEstrenar$word3_c1v1)) / 20)
rolesEstrenar[rolesEstrenar$bin == tOSrolesEstrenar, 'binN'] <- 0
rolesEstrenar$binN <- rolesEstrenar$bin - tOSrolesEstrenar
rolesEstrenar$condToken <- 1

# roles : di : Aceptar
rolesAceptar <- filter(df_lex, target == 'roles' & verb == 'aceptar')
rolesAceptar$binN <- 999
# tOSrolesAceptar <- round(unique(rolesAceptar$targetOffset) / 20)
tOSrolesAceptar <- round((unique(rolesAceptar$targetOffset) - unique(rolesAceptar$word3_c1v1)) / 20)
rolesAceptar[rolesAceptar$bin == tOSrolesAceptar, 'binN'] <- 0
rolesAceptar$binN <- rolesAceptar$bin - tOSrolesAceptar
rolesAceptar$condToken <- 2



# sol : mono : Descubrir
solDescubrir <- filter(df_lex, target == 'sol' & verb == 'descubrir')
solDescubrir$binN <- 999
# tOSsolDescubrir <- round(unique(solDescubrir$targetOffset) / 20)
tOSsolDescubrir <- round((unique(solDescubrir$targetOffset) - unique(solDescubrir$word3_c1v1)) / 20)
solDescubrir[solDescubrir$bin == tOSsolDescubrir, 'binN'] <- 0
solDescubrir$binN <- solDescubrir$bin - tOSsolDescubrir
solDescubrir$condToken <- 1

# sol : mono : Escribir
solEscribir <- filter(df_lex, target == 'sol' & verb == 'escribir')
solEscribir$binN <- 999
# tOSsolEscribir <- round(unique(solEscribir$targetOffset) / 20)
tOSsolEscribir <- round((unique(solEscribir$targetOffset) - unique(solEscribir$word3_c1v1)) / 20)
solEscribir[solEscribir$bin == tOSsolEscribir, 'binN'] <- 0
solEscribir$binN <- solEscribir$bin - tOSsolEscribir
solEscribir$condToken <- 2

# soles : di : Descubrir
solesDescubrir <- filter(df_lex, target == 'soles' & verb == 'descubrir')
solesDescubrir$binN <- 999
# tOSsolesDescubrir <- round(unique(solesDescubrir$targetOffset) / 20)
tOSsolesDescubrir <- round((unique(solesDescubrir$targetOffset) - unique(solesDescubrir$word3_c1v1)) / 20)
solesDescubrir[solesDescubrir$bin == tOSsolesDescubrir, 'binN'] <- 0
solesDescubrir$binN <- solesDescubrir$bin - tOSsolesDescubrir
solesDescubrir$condToken <- 1

# soles : di : Escribir
solesEscribir <- filter(df_lex, target == 'soles' & verb == 'escribir')
solesEscribir$binN <- 999
# tOSsolesEscribir <- round(unique(solesEscribir$targetOffset) / 20)
tOSsolesEscribir <- round((unique(solesEscribir$targetOffset) - unique(solesEscribir$word3_c1v1)) / 20)
solesEscribir[solesEscribir$bin == tOSsolesEscribir, 'binN'] <- 0
solesEscribir$binN <- solesEscribir$bin - tOSsolesEscribir
solesEscribir$condToken <- 2



# tul : mono : Coser
tulCoser <- filter(df_lex, target == 'tul' & verb == 'coser')
tulCoser$binN <- 999
# tOStulCoser <- round(unique(tulCoser$targetOffset) / 20)
tOStulCoser <- round((unique(tulCoser$targetOffset) - unique(tulCoser$word3_c1v1)) / 20)
tulCoser[tulCoser$bin == tOStulCoser, 'binN'] <- 0
tulCoser$binN <- tulCoser$bin - tOStulCoser
tulCoser$condToken <- 1

# tul : mono : Pedir
tulPedir <- filter(df_lex, target == 'tul' & verb == 'pedir')
tulPedir$binN <- 999
# tOStulPedir <- round(unique(tulPedir$targetOffset) / 20)
tOStulPedir <- round((unique(tulPedir$targetOffset) - unique(tulPedir$word3_c1v1)) / 20)
tulPedir[tulPedir$bin == tOStulPedir, 'binN'] <- 0
tulPedir$binN <- tulPedir$bin - tOStulPedir
tulPedir$condToken <- 2

# tules : di : Coser
tulesCoser <- filter(df_lex, target == 'tules' & verb == 'coser')
tulesCoser$binN <- 999
# tOStulesCoser <- round(unique(tulesCoser$targetOffset) / 20)
tOStulesCoser <- round((unique(tulesCoser$targetOffset) - unique(tulesCoser$word3_c1v1)) / 20)
tulesCoser[tulesCoser$bin == tOStulesCoser, 'binN'] <- 0
tulesCoser$binN <- tulesCoser$bin - tOStulesCoser
tulesCoser$condToken <- 1

# tules : di : Pedir
tulesPedir <- filter(df_lex, target == 'tules' & verb == 'pedir')
tulesPedir$binN <- 999
# tOStulesPedir <- round(unique(tulesPedir$targetOffset) / 20)
tOStulesPedir <- round((unique(tulesPedir$targetOffset) - unique(tulesPedir$word3_c1v1)) / 20)
tulesPedir[tulesPedir$bin == tOStulesPedir, 'binN'] <- 0
tulesPedir$binN <- tulesPedir$bin - tOStulesPedir
tulesPedir$condToken <- 2



# 'rbind' all the dataframes together into 1
df_adj <- do.call("rbind", list(chalLucir,
                                chalesLucir,
                                chalRegalar,
                                chalesRegalar,
                                colEscoger,
                                colCenar,
                                colesEscoger,
                                colesCenar,
                                gelEmpacar,
                                gelUsar,
                                gelesEmpacar,
                                gelesUsar,
                                mesExigir,
                                mesConseguir,
                                mesesExigir,
                                mesesConseguir,
                                mielHacer,
                                mielTransportar,
                                mielesHacer,
                                mielesTransportar,
                                rolEstrenar,
                                rolAceptar,
                                rolesEstrenar,
                                rolesAceptar,
                                solDescubrir,
                                solEscribir,
                                solesDescubrir,
                                solesEscribir,
                                tulPedir,
                                tulCoser,
                                tulesCoser,
                                tulesPedir))

# 200ms bin adjustment for VWP
df_adj$binAdj <- df_adj$binN - 10
df_adj <- as.data.frame(df_adj)



# check it 
glimpse(df_adj)

# write table
write.table(df_adj, "./mySources/data/clean/durationBIN10Clean.csv", row.names = F, quote = T, sep = ",")

