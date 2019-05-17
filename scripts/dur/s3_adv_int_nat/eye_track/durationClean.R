#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Morphosyntactic predictability                                              #
# Duration clean                                                              #
# 05/22/2017                                                                  #
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
#dur_temp <- read_tsv("./mySources/data/raw/allExpOutput10.txt")
#dur_temp <- read_tsv("./mySources/data/raw/allExpOutput20.txt")
dur_temp <- read_tsv("./mySources/data/raw/allExpOutput50.txt")

# Check structure of each data frame to make sure they 
# can be combined 
#str(dur1_temp); str(dur2_temp)

glimpse(dur_temp)

# dur2 doesn't have 'sentence', which we don't need, so we remove it
# dur1 <- dur1_temp[, -which(colnames(dur1_temp) == "sentence")]

# They can, so we combine them 
#dur_temp <- rbind(dur1, dur2_temp)

dur_temp <- dur_temp %>% filter(., exp == "duracion")

glimpse(dur_temp)

# rename variables
names(dur_temp)[names(dur_temp)=="RECORDING_SESSION_LABEL"] <- "participant"
names(dur_temp)[names(dur_temp)=="TRIAL_LABEL"] <- "trial"
names(dur_temp)[names(dur_temp)=="BIN_INDEX"] <- "bin"
names(dur_temp)[names(dur_temp)=="id"] <- "wavID"
names(dur_temp)[names(dur_temp)=="RIGHT_IA_1_SAMPLE_COUNT"] <- "targetCount"
names(dur_temp)[names(dur_temp)=="RIGHT_IA_2_SAMPLE_COUNT"] <- "distractorCount"
names(dur_temp)[names(dur_temp)=="RIGHT_IA_1_SAMPLE_COUNT_%"] <- "targetProp"
names(dur_temp)[names(dur_temp)=="RIGHT_IA_2_SAMPLE_COUNT_%"] <- "distractorProp"

# remove unnecessary columns
dur_temp <- select(dur_temp, 
                   -trial, 
                   -identifier, 
                   -sentencewav, 
                   -word1, 
                   -word2_20msafterv1, 
                   -word2_c2, 
                   -word2_c3, 
                   -word2_v2, 
                   -word3_20msafterv1,
                   -word3_c2, 
                   -word3_c3, 
                   -word3_suffix, 
                   -word4_suffix, 
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
dur_temp$target <- as.factor(dur_temp$target)

dur_bi_syl <- c('coles', 'chales', 'geles', 'soles', 'tules', 'meses', 'roles', 'mieles')
dur_temp$condition <- 'monosyllabic'
dur_temp[dur_temp$target %in% dur_bi_syl, 'condition'] <- 'bisyllabic'



# check number of bins for each participant
# so we can remove fluff
dur_temp %>%
  select(., participant, bin) %>%
  aggregate(bin ~ participant, FUN = max, data = .) %>% 
  arrange(., bin)
# Everything looks fine



# Create grouping variable 
# we use the 'participant' column. All the ids have 2 numbers, but they dont 
# all have 2 letters, so we count from the right
durReduced <- dur_temp %>%
  separate(., col = participant, into = c('group', 'id'), sep = -3, remove = FALSE) %>% 
  select(., -id)




# Some levels are erroneously labelled 's' or 'l'
durReduced[durReduced$group == 'l', 'group'] <- 'lb'
durReduced[durReduced$group == 's', 'group'] <- 'ss'

# Change 'in' level of 'group' to 'int'
# because 'in' is a function in R
durReduced[durReduced$group == 'in', 'group'] <- 'int'

durReduced$group <- as.factor(durReduced$group)

summary(durReduced$group)
glimpse(durReduced)


durReduced$condToken <- 1
durReduced[durReduced$target == "chal" & durReduced$verb == "regalar", 'condToken'] <- 2
durReduced[durReduced$target == "chales" & durReduced$verb == "regalar", 'condToken'] <- 2
durReduced[durReduced$target == "col" & durReduced$verb == "escoger", 'condToken'] <- 2
durReduced[durReduced$target == "coles" & durReduced$verb == "escoger", 'condToken'] <- 2
durReduced[durReduced$target == "gel" & durReduced$verb == "usar", 'condToken'] <- 2
durReduced[durReduced$target == "geles" & durReduced$verb == "usar", 'condToken'] <- 2
durReduced[durReduced$target == "mes" & durReduced$verb == "conseguir", 'condToken'] <- 2
durReduced[durReduced$target == "meses" & durReduced$verb == "conseguir", 'condToken'] <- 2
durReduced[durReduced$target == "miel" & durReduced$verb == "transportar", 'condToken'] <- 2
durReduced[durReduced$target == "mieles" & durReduced$verb == "transportar", 'condToken'] <- 2
durReduced[durReduced$target == "rol" & durReduced$verb == "aceptar", 'condToken'] <- 2
durReduced[durReduced$target == "roles" & durReduced$verb == "aceptar", 'condToken'] <- 2
durReduced[durReduced$target == "sol" & durReduced$verb == "escribir", 'condToken'] <- 2
durReduced[durReduced$target == "soles" & durReduced$verb == "escribir", 'condToken'] <- 2
durReduced[durReduced$target == "tul" & durReduced$verb == "pedir", 'condToken'] <- 2
durReduced[durReduced$target == "tules" & durReduced$verb == "pedir", 'condToken'] <- 2



# - for the monosyllabic targets, the end of the word/coda/first syllable is: word5
# - for bisyllabic targets, the end of the first consonant of the second syllable is: not available 
# - we have to get the time point of the end of the first consonant in the second syllable
#   in praat, by hand. 
# - We save these values to one of two named vectors (1 for each verb)
# - we add a new column called 'targetOffset' 
# - initialy it is just a copy of 'word5'
# - then we add the vectors of hand calculated times

durReduced$targetOffset <- durReduced$word5

# condToken 1 times
offset_cond1 <- c('coles'  = '2347', # escoger
                  'chales' = '2213', # lucir
                  'geles'  = '2188', # empacar
                  'soles'  = '2131', # descubir
                  'tules'  = '1862', # coser
                  'meses'  = '2566', # exigir
                  'roles'  = '2319', # aceptar
                  'mieles' = '2121'  # hacer
  )

durReduced[durReduced$target == "coles"  & durReduced$verb == "escoger", 'targetOffset']  <- 2347
durReduced[durReduced$target == "chales" & durReduced$verb == "lucir", 'targetOffset']    <- 2213
durReduced[durReduced$target == "geles"  & durReduced$verb == "empacar", 'targetOffset']  <- 2188
durReduced[durReduced$target == "soles"  & durReduced$verb == "descubir", 'targetOffset'] <- 2131
durReduced[durReduced$target == "tules"  & durReduced$verb == "coser", 'targetOffset']    <- 1862
durReduced[durReduced$target == "meses"  & durReduced$verb == "exigir", 'targetOffset']   <- 2566
durReduced[durReduced$target == "roles"  & durReduced$verb == "aceptar", 'targetOffset']  <- 2319
durReduced[durReduced$target == "mieles" & durReduced$verb == "hacer", 'targetOffset']    <- 2121

# condToken 2 times
offset_cond2 <- c('coles'  = '2177', # cenar
                  'chales' = '2289', # regalar
                  'geles'  = '1960', # usar
                  'soles'  = '2109', # escribir
                  'tules'  = '1918', # pedir
                  'meses'  = '2749', # conseguir
                  'roles'  = '2116', # estrenar
                  'mieles' = '2380'  # transportar
  )

durReduced[durReduced$target == "coles"  & durReduced$verb == "cenar", 'targetOffset']       <- 2177
durReduced[durReduced$target == "chales" & durReduced$verb == "regalar", 'targetOffset']     <- 2289
durReduced[durReduced$target == "geles"  & durReduced$verb == "usar", 'targetOffset']        <- 1960
durReduced[durReduced$target == "soles"  & durReduced$verb == "escribir", 'targetOffset']    <- 2109
durReduced[durReduced$target == "tules"  & durReduced$verb == "pedir", 'targetOffset']       <- 1918
durReduced[durReduced$target == "meses"  & durReduced$verb == "conseguir", 'targetOffset']   <- 2749
durReduced[durReduced$target == "roles"  & durReduced$verb == "estrenar", 'targetOffset']    <- 2116
durReduced[durReduced$target == "mieles" & durReduced$verb == "transportar", 'targetOffset'] <- 2380


df_dur <- durReduced

df_dur <- droplevels(df_dur)

glimpse(df_dur)


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


xtabs(~ target + verb, df_dur)


# chal : mono : lucir
chalLucir <- filter(df_dur, target == 'chal' & verb == 'lucir')
chalLucir$binN <- 999
tOSchalLucir <- round(unique(chalLucir$targetOffset) / 50)
#tOSchalLucir <- round((unique(chalLucir$targetOffset) - unique(chalLucir$word3_c1v1)) / 20)
chalLucir[chalLucir$bin == tOSchalLucir, 'binN'] <- 0
chalLucir$binN <- chalLucir$bin - tOSchalLucir
chalLucir$condToken <- 1

# chal : mono : regalar
chalRegalar <- filter(df_dur, target == 'chal' & verb == 'regalar')
chalRegalar$binN <- 999
tOSchalRegalar <- round(unique(chalRegalar$targetOffset) / 50)
#tOSchalRegalar <- round((unique(chalRegalar$targetOffset) - unique(chalRegalar$word3_c1v1)) / 20)
chalRegalar[chalRegalar$bin == tOSchalRegalar, 'binN'] <- 0
chalRegalar$binN <- chalRegalar$bin - tOSchalRegalar
chalRegalar$condToken <- 2

# chales : di : Lucir
chalesLucir <- filter(df_dur, target == 'chales' & verb == 'lucir')
chalesLucir$binN <- 999
tOSchalesLucir <- round(unique(chalesLucir$targetOffset) / 50)
#tOSchalesLucir <- round((unique(chalesLucir$targetOffset) - unique(chalesLucir$word3_c1v1)) / 20)
chalesLucir[chalesLucir$bin == tOSchalesLucir, 'binN'] <- 0
chalesLucir$binN <- chalesLucir$bin - tOSchalesLucir
chalesLucir$condToken <- 1

# chales : di : Regalar
chalesRegalar <- filter(df_dur, target == 'chales' & verb == 'regalar')
chalesRegalar$binN <- 999
tOSchalesRegalar <- round(unique(chalesRegalar$targetOffset) / 50)
#tOSchalesRegalar <- round((unique(chalesRegalar$targetOffset) - unique(chalesRegalar$word3_c1v1)) / 20)
chalesRegalar[chalesRegalar$bin == tOSchalesRegalar, 'binN'] <- 0
chalesRegalar$binN <- chalesRegalar$bin - tOSchalesRegalar
chalesRegalar$condToken <- 2




# col : mono : Cenar
colCenar <- filter(df_dur, target == 'col' & verb == 'cenar')
colCenar$binN <- 999
tOScolCenar <- round(unique(colCenar$targetOffset) / 50)
#tOScolCenar <- round((unique(colCenar$targetOffset) - unique(colCenar$word3_c1v1)) / 20)
colCenar[colCenar$bin == tOScolCenar, 'binN'] <- 0
colCenar$binN <- colCenar$bin - tOScolCenar
colCenar$condToken <- 1

# col : mono : Escoger
colEscoger <- filter(df_dur, target == 'col' & verb == 'escoger')
colEscoger$binN <- 999
tOScolEscoger <- round(unique(colEscoger$targetOffset) / 50)
#tOScolEscoger <- round((unique(colEscoger$targetOffset) - unique(colEscoger$word3_c1v1)) / 20)
colEscoger[colEscoger$bin == tOScolEscoger, 'binN'] <- 0
colEscoger$binN <- colEscoger$bin - tOScolEscoger
colEscoger$condToken <- 2

# coles : di : Cenar
colesCenar <- filter(df_dur, target == 'coles' & verb == 'cenar')
colesCenar$binN <- 999
tOScolesCenar <- round(unique(colesCenar$targetOffset) / 50)
#tOScolesCenar <- round((unique(colesCenar$targetOffset) - unique(colesCenar$word3_c1v1)) / 20)
colesCenar[colesCenar$bin == tOScolesCenar, 'binN'] <- 0
colesCenar$binN <- colesCenar$bin - tOScolesCenar
colesCenar$condToken <- 1

# coles : di : Escoger
colesEscoger <- filter(df_dur, target == 'coles' & verb == 'escoger')
colesEscoger$binN <- 999
tOScolesEscoger <- round(unique(colesEscoger$targetOffset) / 50)
#tOScolesEscoger <- round((unique(colesEscoger$targetOffset) - unique(colesEscoger$word3_c1v1)) / 20)
colesEscoger[colesEscoger$bin == tOScolesEscoger, 'binN'] <- 0
colesEscoger$binN <- colesEscoger$bin - tOScolesEscoger
colesEscoger$condToken <- 2



# gel : mono : Empacar
gelEmpacar <- filter(df_dur, target == 'gel' & verb == 'empacar')
gelEmpacar$binN <- 999
tOSgelEmpacar <- round(unique(gelEmpacar$targetOffset) / 50)
#tOSgelEmpacar <- round((unique(gelEmpacar$targetOffset) - unique(gelEmpacar$word3_c1v1)) / 20)
gelEmpacar[gelEmpacar$bin == tOSgelEmpacar, 'binN'] <- 0
gelEmpacar$binN <- gelEmpacar$bin - tOSgelEmpacar
gelEmpacar$condToken <- 1

# gel : mono : Usar
gelUsar <- filter(df_dur, target == 'gel' & verb == 'usar')
gelUsar$binN <- 999
tOSgelUsar <- round(unique(gelUsar$targetOffset) / 50)
#tOSgelUsar <- round((unique(gelUsar$targetOffset) - unique(gelUsar$word3_c1v1)) / 20)
gelUsar[gelUsar$bin == tOSgelUsar, 'binN'] <- 0
gelUsar$binN <- gelUsar$bin - tOSgelUsar
gelUsar$condToken <- 2

# geles : di : Empacar
gelesEmpacar <- filter(df_dur, target == 'geles' & verb == 'empacar')
gelesEmpacar$binN <- 999
tOSgelesEmpacar <- round(unique(gelesEmpacar$targetOffset) / 50)
#tOSgelesEmpacar <- round((unique(gelesEmpacar$targetOffset) - unique(gelesEmpacar$word3_c1v1)) / 20)
gelesEmpacar[gelesEmpacar$bin == tOSgelesEmpacar, 'binN'] <- 0
gelesEmpacar$binN <- gelesEmpacar$bin - tOSgelesEmpacar
gelesEmpacar$condToken <- 1

# geles : di : Usar
gelesUsar <- filter(df_dur, target == 'geles' & verb == 'usar')
gelesUsar$binN <- 999
tOSgelesUsar <- round(unique(gelesUsar$targetOffset) / 50)
#tOSgelesUsar <- round((unique(gelesUsar$targetOffset) - unique(gelesUsar$word3_c1v1)) / 20)
gelesUsar[gelesUsar$bin == tOSgelesUsar, 'binN'] <- 0
gelesUsar$binN <- gelesUsar$bin - tOSgelesUsar
gelesUsar$condToken <- 2



# mes : mono : Exigir
mesExigir <- filter(df_dur, target == 'mes' & verb == 'exigir')
mesExigir$binN <- 999
tOSmesExigir <- round(unique(mesExigir$targetOffset) / 50)
#tOSmesExigir <- round((unique(mesExigir$targetOffset) - unique(mesExigir$word3_c1v1)) / 20)
mesExigir[mesExigir$bin == tOSmesExigir, 'binN'] <- 0
mesExigir$binN <- mesExigir$bin - tOSmesExigir
mesExigir$condToken <- 1

# mes : mono : Conseguir
mesConseguir <- filter(df_dur, target == 'mes' & verb == 'conseguir')
mesConseguir$binN <- 999
tOSmesConseguir <- round(unique(mesConseguir$targetOffset) / 50)
#tOSmesConseguir <- round((unique(mesConseguir$targetOffset) - unique(mesConseguir$word3_c1v1)) / 20)
mesConseguir[mesConseguir$bin == tOSmesConseguir, 'binN'] <- 0
mesConseguir$binN <- mesConseguir$bin - tOSmesConseguir
mesConseguir$condToken <- 2

# meses : di : Exigir
mesesExigir <- filter(df_dur, target == 'meses' & verb == 'exigir')
mesesExigir$binN <- 999
tOSmesesExigir <- round(unique(mesesExigir$targetOffset) / 50)
#tOSmesesExigir <- round((unique(mesesExigir$targetOffset) - unique(mesesExigir$word3_c1v1)) / 20)
mesesExigir[mesesExigir$bin == tOSmesesExigir, 'binN'] <- 0
mesesExigir$binN <- mesesExigir$bin - tOSmesesExigir
mesesExigir$condToken <- 1

# meses : di : Conseguir
mesesConseguir <- filter(df_dur, target == 'meses' & verb == 'conseguir')
mesesConseguir$binN <- 999
tOSmesesConseguir <- round(unique(mesesConseguir$targetOffset) / 50)
#tOSmesesConseguir <- round((unique(mesesConseguir$targetOffset) - unique(mesesConseguir$word3_c1v1)) / 20)
mesesConseguir[mesesConseguir$bin == tOSmesesConseguir, 'binN'] <- 0
mesesConseguir$binN <- mesesConseguir$bin - tOSmesesConseguir
mesesConseguir$condToken <- 2


# miel : mono : Hacer
mielHacer <- filter(df_dur, target == 'miel' & verb == 'hacer')
mielHacer$binN <- 999
tOSmielHacer <- round(unique(mielHacer$targetOffset) / 50)
#tOSmielHacer <- round((unique(mielHacer$targetOffset) - unique(mielHacer$word3_c1v1)) / 20)
mielHacer[mielHacer$bin == tOSmielHacer, 'binN'] <- 0
mielHacer$binN <- mielHacer$bin - tOSmielHacer
mielHacer$condToken <- 1

# miel : mono : Transportar
mielTransportar <- filter(df_dur, target == 'miel' & verb == 'transportar')
mielTransportar$binN <- 999
tOSmielTransportar <- round(unique(mielTransportar$targetOffset) / 50)
#tOSmielTransportar <- round((unique(mielTransportar$targetOffset) - unique(mielTransportar$word3_c1v1)) / 20)
mielTransportar[mielTransportar$bin == tOSmielTransportar, 'binN'] <- 0
mielTransportar$binN <- mielTransportar$bin - tOSmielTransportar
mielTransportar$condToken <- 2

# mieles : di : Hacer
mielesHacer <- filter(df_dur, target == 'mieles' & verb == 'hacer')
mielesHacer$binN <- 999
tOSmielesHacer <- round(unique(mielesHacer$targetOffset) / 50)
#tOSmielesHacer <- round((unique(mielesHacer$targetOffset) - unique(mielesHacer$word3_c1v1)) / 20)
mielesHacer[mielesHacer$bin == tOSmielesHacer, 'binN'] <- 0
mielesHacer$binN <- mielesHacer$bin - tOSmielesHacer
mielesHacer$condToken <- 1

# mieles : di : Transportar
mielesTransportar <- filter(df_dur, target == 'mieles' & verb == 'transportar')
mielesTransportar$binN <- 999
tOSmielesTransportar <- round(unique(mielesTransportar$targetOffset) / 50)
#tOSmielesTransportar <- round((unique(mielesTransportar$targetOffset) - unique(mielesTransportar$word3_c1v1)) / 20)
mielesTransportar[mielesTransportar$bin == tOSmielesTransportar, 'binN'] <- 0
mielesTransportar$binN <- mielesTransportar$bin - tOSmielesTransportar
mielesTransportar$condToken <- 2



# rol : mono : Estrenar
rolEstrenar <- filter(df_dur, target == 'rol' & verb == 'estrenar')
rolEstrenar$binN <- 999
tOSrolEstrenar <- round(unique(rolEstrenar$targetOffset) / 50)
#tOSrolEstrenar <- round((unique(rolEstrenar$targetOffset) - unique(rolEstrenar$word3_c1v1)) / 20)
rolEstrenar[rolEstrenar$bin == tOSrolEstrenar, 'binN'] <- 0
rolEstrenar$binN <- rolEstrenar$bin - tOSrolEstrenar
rolEstrenar$condToken <- 1

# rol : mono : Aceptar
rolAceptar <- filter(df_dur, target == 'rol' & verb == 'aceptar')
rolAceptar$binN <- 999
tOSrolAceptar <- round(unique(rolAceptar$targetOffset) / 50)
#tOSrolAceptar <- round((unique(rolAceptar$targetOffset) - unique(rolAceptar$word3_c1v1)) / 20)
rolAceptar[rolAceptar$bin == tOSrolAceptar, 'binN'] <- 0
rolAceptar$binN <- rolAceptar$bin - tOSrolAceptar
rolAceptar$condToken <- 2

# roles : di : Estrenar
rolesEstrenar <- filter(df_dur, target == 'roles' & verb == 'estrenar')
rolesEstrenar$binN <- 999
tOSrolesEstrenar <- round(unique(rolesEstrenar$targetOffset) / 50)
#tOSrolesEstrenar <- round((unique(rolesEstrenar$targetOffset) - unique(rolEstrenar$word3_c1v1)) / 20)
rolesEstrenar[rolesEstrenar$bin == tOSrolesEstrenar, 'binN'] <- 0
rolesEstrenar$binN <- rolesEstrenar$bin - tOSrolesEstrenar
rolesEstrenar$condToken <- 1

# roles : di : Aceptar
rolesAceptar <- filter(df_dur, target == 'roles' & verb == 'aceptar')
rolesAceptar$binN <- 999
tOSrolesAceptar <- round(unique(rolesAceptar$targetOffset) / 50)
#tOSrolesAceptar <- round((unique(rolesAceptar$targetOffset) - unique(rolesAceptar$word3_c1v1)) / 20)
rolesAceptar[rolesAceptar$bin == tOSrolesAceptar, 'binN'] <- 0
rolesAceptar$binN <- rolesAceptar$bin - tOSrolesAceptar
rolesAceptar$condToken <- 2



# sol : mono : Descubrir
solDescubrir <- filter(df_dur, target == 'sol' & verb == 'descubrir')
solDescubrir$binN <- 999
tOSsolDescubrir <- round(unique(solDescubrir$targetOffset) / 50)
#tOSsolDescubrir <- round((unique(solDescubrir$targetOffset) - unique(solDescubrir$word3_c1v1)) / 20)
solDescubrir[solDescubrir$bin == tOSsolDescubrir, 'binN'] <- 0
solDescubrir$binN <- solDescubrir$bin - tOSsolDescubrir
solDescubrir$condToken <- 1

# sol : mono : Escribir
solEscribir <- filter(df_dur, target == 'sol' & verb == 'escribir')
solEscribir$binN <- 999
tOSsolEscribir <- round(unique(solEscribir$targetOffset) / 50)
#tOSsolEscribir <- round((unique(solEscribir$targetOffset) - unique(solEscribir$word3_c1v1)) / 20)
solEscribir[solEscribir$bin == tOSsolEscribir, 'binN'] <- 0
solEscribir$binN <- solEscribir$bin - tOSsolEscribir
solEscribir$condToken <- 2

# soles : di : Descubrir
solesDescubrir <- filter(df_dur, target == 'soles' & verb == 'descubrir')
solesDescubrir$binN <- 999
tOSsolesDescubrir <- round(unique(solesDescubrir$targetOffset) / 50)
#tOSsolesDescubrir <- round((unique(solesDescubrir$targetOffset) - unique(solesDescubrir$word3_c1v1)) / 20)
solesDescubrir[solesDescubrir$bin == tOSsolesDescubrir, 'binN'] <- 0
solesDescubrir$binN <- solesDescubrir$bin - tOSsolesDescubrir
solesDescubrir$condToken <- 1

# soles : di : Escribir
solesEscribir <- filter(df_dur, target == 'soles' & verb == 'escribir')
solesEscribir$binN <- 999
tOSsolesEscribir <- round(unique(solesEscribir$targetOffset) / 50)
#tOSsolesEscribir <- round((unique(solesEscribir$targetOffset) - unique(solesEscribir$word3_c1v1)) / 20)
solesEscribir[solesEscribir$bin == tOSsolesEscribir, 'binN'] <- 0
solesEscribir$binN <- solesEscribir$bin - tOSsolesEscribir
solesEscribir$condToken <- 2



# tul : mono : Coser
tulCoser <- filter(df_dur, target == 'tul' & verb == 'coser')
tulCoser$binN <- 999
tOStulCoser <- round(unique(tulCoser$targetOffset) / 50)
#tOStulCoser <- round((unique(tulCoser$targetOffset) - unique(tulCoser$word3_c1v1)) / 20)
tulCoser[tulCoser$bin == tOStulCoser, 'binN'] <- 0
tulCoser$binN <- tulCoser$bin - tOStulCoser
tulCoser$condToken <- 1

# tul : mono : Pedir
tulPedir <- filter(df_dur, target == 'tul' & verb == 'pedir')
tulPedir$binN <- 999
tOStulPedir <- round(unique(tulPedir$targetOffset) / 50)
#tOStulPedir <- round((unique(tulPedir$targetOffset) - unique(tulPedir$word3_c1v1)) / 20)
tulPedir[tulPedir$bin == tOStulPedir, 'binN'] <- 0
tulPedir$binN <- tulPedir$bin - tOStulPedir
tulPedir$condToken <- 2

# tules : di : Coser
tulesCoser <- filter(df_dur, target == 'tules' & verb == 'coser')
tulesCoser$binN <- 999
tOStulesCoser <- round(unique(tulesCoser$targetOffset) / 50)
#tOStulesCoser <- round((unique(tulesCoser$targetOffset) - unique(tulesCoser$word3_c1v1)) / 20)
tulesCoser[tulesCoser$bin == tOStulesCoser, 'binN'] <- 0
tulesCoser$binN <- tulesCoser$bin - tOStulesCoser
tulesCoser$condToken <- 1

# tules : di : Pedir
tulesPedir <- filter(df_dur, target == 'tules' & verb == 'pedir')
tulesPedir$binN <- 999
tOStulesPedir <- round(unique(tulesPedir$targetOffset) / 50)
#tOStulesPedir <- round((unique(tulesPedir$targetOffset) - unique(tulesPedir$word3_c1v1)) / 20)
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








glimpse(df_dur)




# 200ms bin adjustment for VWP
df_adj$binAdj <- df_adj$binN - 22
df_adj <- as.data.frame(df_adj)



# check it 
glimpse(df_adj)

# write table
#write.table(df_adj, "./mySources/data/clean/durationBIN10Clean.csv", row.names = F, quote = T, sep = ",")
#write.table(df_adj, "./mySources/data/clean/durationBIN20Clean.csv", row.names = F, quote = T, sep = ",")
write.table(df_adj, "./mySources/data/clean/durationBIN50Clean.csv", row.names = F, quote = T, sep = ",")
