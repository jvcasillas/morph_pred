#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Morphosyntactic predictability                                              #
# Duration clean                                                              #
# 03/29/2017                                                                  #
# Script 1                                                                    #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#


# clean working directory
rm(list = ls(all = TRUE))

# Set working directory
# setwd("~/Desktop/morph_pred/")
setwd("~/academia/research/in_progress/morph_pred")


library(dplyr); library(tidyr); library(ggplot2); library(plotly)
library(lme4); library(lmerTest)


# Read data
dur1 <- read.csv("./mySources/data/raw/durationOutput_trialBIN5.csv", header = TRUE, quote = "", sep = '\t')
dur2 <- read.csv("./mySources/data/raw/durationOutput2_trialBIN5.csv", header = TRUE, quote = "", sep = '\t')

# Check structure of each data frame to make sure they 
# can be combined 
str(dur1); str(dur2)

# They can, so we combine them 
dur <- rbind(dur1, dur2)


# rename variables
names(dur)[names(dur)=="RECORDING_SESSION_LABEL"] <- "participant"
names(dur)[names(dur)=="TRIAL_LABEL"] <- "trial"
names(dur)[names(dur)=="TRIAL_CONDITION"] <- "condLabel"
names(dur)[names(dur)=="CURRENT_BIN"] <- "bin"
names(dur)[names(dur)=="RIGHT_0_C"] <- "lookElsewhere"
names(dur)[names(dur)=="RIGHT_1_C"] <- "lookTarget"
names(dur)[names(dur)=="RIGHT_2_C"] <- "lookDistractor"

# remove unwanted characters from column (accent marks)
dur$condLabel <- gsub("\227", "o", paste(dur$condLabel))
dur$condLabel <- as.factor(dur$condLabel)


# remove distractors 
durTargets <- c('col', 'coles', 'chal', 'chales', 'gel', 'geles', 'sol', 'soles', 'tul', 'tules', 'mes', 'meses', 'rol', 'roles', 'miel', 'mieles')
dur <- dur[dur$condLabel %in% durTargets, ]
# Check to see if it worked
summary(dur$condLabel)




# Set condition variable (monosyllabic or disyllabic?)
dur$condition <- 'disyllabic'
dur[dur$condLabel == 'col' | 
    dur$condLabel == 'chal' | 
    dur$condLabel == 'gel' | 
    dur$condLabel == 'sol' | 
    dur$condLabel == 'tul' | 
    dur$condLabel == 'mes' | 
    dur$condLabel == 'rol' | 
    dur$condLabel == 'miel', 'condition'] <- 'monosyllabic'
dur$condition <- as.factor(dur$condition)


# check number of bins for each participant
# so we can remove fluff
dur %>%
  select(., participant, bin) %>%
  aggregate(bin ~ participant, FUN = max, data = .) %>% 
  arrange(., bin)

# Participants 'l14' and 'test' have no data so we remove them 
dur <- filter(dur, bin != 0)

# The minimum bin = 425



# Calculate proportions
durReduced <- dur %>%
  select(., participant, trial, condLabel, IA_LABEL_2, TRIAL_START_TIME, 
            BIN_START_TIME, BIN_END_TIME, verb, condition, bin, lookTarget, 
            lookDistractor, lookElsewhere, word3_c1v1:word7) %>%
  group_by(., participant, trial, condLabel, bin) %>% 
  mutate(., word3_c1v1 = as.numeric(word3_c1v1), 
            word4_c1v1 = as.numeric(word4_c1v1),
            word4_v1_20ms = as.numeric(word4_c1) + 20,
            word4_c1 = as.numeric(word4_c1), 
            word4_suffix = as.numeric(word4_suffix), 
            word5 = as.numeric(word5),
            word6 = as.numeric(word6),
            word7 = as.numeric(word7),
            binSize = 5,
            propT = (lookTarget / binSize), 
            propD = (lookDistractor / binSize), 
            propE = (lookElsewhere / binSize)) 

# check number of bins for each participant
durReduced %>%
  select(., participant, bin) %>%
  aggregate(bin ~ participant, FUN = max, data = .) %>% 
  arrange(., bin)


# Create grouping variable 
durReduced %>%
  separate(., col = participant, into = c('group', 'id'), sep = -3, remove = FALSE) %>% 
  select(., -id) -> durReduced

# Some levels are erroneously labelled 's' or 'l'
durReduced[durReduced$group == 'l', 'group'] <- 'lb'
durReduced[durReduced$group == 's', 'group'] <- 'ss'

# Change 'in' level of 'group' to 'int'
# because 'in' is a function in R
durReduced[durReduced$group == 'in', 'group'] <- 'int'

durReduced$group <- as.factor(durReduced$group)

summary(durReduced$group)
glimpse(durReduced)





df_dur <- durReduced

df_dur <- droplevels(df_dur)

glimpse(df_dur)


# Item-specific time normalization
# - create subset of data based on item (for both monosyllabic and disyllabic)
# - use item-specific target onset as 0 (in binN column)
# - save item-specific bin min and max as variables 
# - count forward from target onset to max 
# - count backward from target onset to min
# - 'move' data 200ms for VWP (duplicate binN column + 20)
# - find common min/max between all items (lowest max, highest min)
# - must distinguish between iterations of same words (different preceeding verb)


xtabs(~ condLabel + verb, df_dur)

# chal : mono : lucir
chalLucir <- filter(df_dur, condLabel == 'chal' & verb == 'lucir')
chalLucir$binN <- 999
tOSchalLucir <- round((unique(chalLucir$word4_v1_20ms) - unique(chalLucir$word4_c1v1)) / 5)
chalLucir[chalLucir$bin == tOSchalLucir, 'binN'] <- 0
chalLucir$binN <- chalLucir$bin - tOSchalLucir
chalLucir$condToken <- 1

# chal : mono : regalar
chalRegalar <- filter(df_dur, condLabel == 'chal' & verb == 'regalar')
chalRegalar$binN <- 999
tOSchalRegalar <- round((unique(chalRegalar$word4_v1_20ms) - unique(chalRegalar$word4_c1v1)) / 5)
chalRegalar[chalRegalar$bin == tOSchalRegalar, 'binN'] <- 0
chalRegalar$binN <- chalRegalar$bin - tOSchalRegalar
chalRegalar$condToken <- 2

# chales : di : Lucir
chalesLucir <- filter(df_dur, condLabel == 'chales' & verb == 'lucir')
chalesLucir$binN <- 999
tOSchalesLucir <- round((unique(chalesLucir$word4_v1_20ms) - unique(chalesLucir$word4_c1v1)) / 5)
chalesLucir[chalesLucir$bin == tOSchalesLucir, 'binN'] <- 0
chalesLucir$binN <- chalesLucir$bin - tOSchalesLucir
chalesLucir$condToken <- 1

# chales : di : Regalar
chalesRegalar <- filter(df_dur, condLabel == 'chales' & verb == 'regalar')
chalesRegalar$binN <- 999
tOSchalesRegalar <- round((unique(chalesRegalar$word4_v1_20ms) - unique(chalesRegalar$word4_c1v1)) / 5)
chalesRegalar[chalesRegalar$bin == tOSchalesRegalar, 'binN'] <- 0
chalesRegalar$binN <- chalesRegalar$bin - tOSchalesRegalar
chalesRegalar$condToken <- 2




# col : mono : Cenar
colCenar <- filter(df_dur, condLabel == 'col' & verb == 'cenar')
colCenar$binN <- 999
tOScolCenar <- round((unique(colCenar$word4_v1_20ms) - unique(colCenar$word4_c1v1)) / 5)
colCenar[colCenar$bin == tOScolCenar, 'binN'] <- 0
colCenar$binN <- colCenar$bin - tOScolCenar
colCenar$condToken <- 1

# col : mono : Escoger
colEscoger <- filter(df_dur, condLabel == 'col' & verb == 'escoger')
colEscoger$binN <- 999
tOScolEscoger <- round((unique(colEscoger$word4_v1_20ms) - unique(colEscoger$word4_c1v1)) / 5)
colEscoger[colEscoger$bin == tOScolEscoger, 'binN'] <- 0
colEscoger$binN <- colEscoger$bin - tOScolEscoger
colEscoger$condToken <- 2

# coles : di : Cenar
colesCenar <- filter(df_dur, condLabel == 'coles' & verb == 'cenar')
colesCenar$binN <- 999
tOScolesCenar <- round((unique(colesCenar$word4_v1_20ms) - unique(colesCenar$word4_c1v1)) / 5)
colesCenar[colesCenar$bin == tOScolesCenar, 'binN'] <- 0
colesCenar$binN <- colesCenar$bin - tOScolesCenar
colesCenar$condToken <- 1

# coles : di : Escoger
colesEscoger <- filter(df_dur, condLabel == 'coles' & verb == 'escoger')
colesEscoger$binN <- 999
tOScolesEscoger <- round((unique(colesEscoger$word4_v1_20ms) - unique(colesEscoger$word4_c1v1)) / 5)
colesEscoger[colesEscoger$bin == tOScolesEscoger, 'binN'] <- 0
colesEscoger$binN <- colesEscoger$bin - tOScolesEscoger
colesEscoger$condToken <- 2



# gel : mono : Empacar
gelEmpacar <- filter(df_dur, condLabel == 'gel' & verb == 'empacar')
gelEmpacar$binN <- 999
tOSgelEmpacar <- round((unique(gelEmpacar$word4_v1_20ms) - unique(gelEmpacar$word4_c1v1)) / 5)
gelEmpacar[gelEmpacar$bin == tOSgelEmpacar, 'binN'] <- 0
gelEmpacar$binN <- gelEmpacar$bin - tOSgelEmpacar
gelEmpacar$condToken <- 1

# gel : mono : Usar
gelUsar <- filter(df_dur, condLabel == 'gel' & verb == 'usar')
gelUsar$binN <- 999
tOSgelUsar <- round((unique(gelUsar$word4_v1_20ms) - unique(gelUsar$word4_c1v1)) / 5)
gelUsar[gelUsar$bin == tOSgelUsar, 'binN'] <- 0
gelUsar$binN <- gelUsar$bin - tOSgelUsar
gelUsar$condToken <- 2

# geles : di : Empacar
gelesEmpacar <- filter(df_dur, condLabel == 'geles' & verb == 'empacar')
gelesEmpacar$binN <- 999
tOSgelesEmpacar <- round((unique(gelesEmpacar$word4_v1_20ms) - unique(gelesEmpacar$word4_c1v1)) / 5)
gelesEmpacar[gelesEmpacar$bin == tOSgelesEmpacar, 'binN'] <- 0
gelesEmpacar$binN <- gelesEmpacar$bin - tOSgelesEmpacar
gelesEmpacar$condToken <- 1

# geles : di : Usar
gelesUsar <- filter(df_dur, condLabel == 'geles' & verb == 'usar')
gelesUsar$binN <- 999
tOSgelesUsar <- round((unique(gelesUsar$word4_v1_20ms) - unique(gelesUsar$word4_c1v1)) / 5)
gelesUsar[gelesUsar$bin == tOSgelesUsar, 'binN'] <- 0
gelesUsar$binN <- gelesUsar$bin - tOSgelesUsar
gelesUsar$condToken <- 2



# mes : mono : Exigir
mesExigir <- filter(df_dur, condLabel == 'mes' & verb == 'exigir')
mesExigir$binN <- 999
tOSmesExigir <- round((unique(mesExigir$word4_v1_20ms) - unique(mesExigir$word4_c1v1)) / 5)
mesExigir[mesExigir$bin == tOSmesExigir, 'binN'] <- 0
mesExigir$binN <- mesExigir$bin - tOSmesExigir
mesExigir$condToken <- 1

# mes : mono : Conseguir
mesConseguir <- filter(df_dur, condLabel == 'mes' & verb == 'conseguir')
mesConseguir$binN <- 999
tOSmesConseguir <- round((unique(mesConseguir$word4_v1_20ms) - unique(mesConseguir$word4_c1v1)) / 5)
mesConseguir[mesConseguir$bin == tOSmesConseguir, 'binN'] <- 0
mesConseguir$binN <- mesConseguir$bin - tOSmesConseguir
mesConseguir$condToken <- 2

# meses : di : Exigir
mesesExigir <- filter(df_dur, condLabel == 'meses' & verb == 'exigir')
mesesExigir$binN <- 999
tOSmesesExigir <- round((unique(mesesExigir$word4_v1_20ms) - unique(mesesExigir$word4_c1v1)) / 5)
mesesExigir[mesesExigir$bin == tOSmesesExigir, 'binN'] <- 0
mesesExigir$binN <- mesesExigir$bin - tOSmesesExigir
mesesExigir$condToken <- 1

# meses : di : Conseguir
mesesConseguir <- filter(df_dur, condLabel == 'meses' & verb == 'conseguir')
mesesConseguir$binN <- 999
tOSmesesConseguir <- round((unique(mesesConseguir$word4_v1_20ms) - unique(mesesConseguir$word4_c1v1)) / 5)
mesesConseguir[mesesConseguir$bin == tOSmesesConseguir, 'binN'] <- 0
mesesConseguir$binN <- mesesConseguir$bin - tOSmesesConseguir
mesesConseguir$condToken <- 2


# miel : mono : Hacer
mielHacer <- filter(df_dur, condLabel == 'miel' & verb == 'hacer')
mielHacer$binN <- 999
tOSmielHacer <- round((unique(mielHacer$word4_v1_20ms) - unique(mielHacer$word4_c1v1)) / 5)
mielHacer[mielHacer$bin == tOSmielHacer, 'binN'] <- 0
mielHacer$binN <- mielHacer$bin - tOSmielHacer
mielHacer$condToken <- 1

# miel : mono : Transportar
mielTransportar <- filter(df_dur, condLabel == 'miel' & verb == 'transportar')
mielTransportar$binN <- 999
tOSmielTransportar <- round((unique(mielTransportar$word4_v1_20ms) - unique(mielTransportar$word4_c1v1)) / 5)
mielTransportar[mielTransportar$bin == tOSmielTransportar, 'binN'] <- 0
mielTransportar$binN <- mielTransportar$bin - tOSmielTransportar
mielTransportar$condToken <- 2

# mieles : di : Hacer
mielesHacer <- filter(df_dur, condLabel == 'mieles' & verb == 'hacer')
mielesHacer$binN <- 999
tOSmielesHacer <- round((unique(mielesHacer$word4_v1_20ms) - unique(mielesHacer$word4_c1v1)) / 5)
mielesHacer[mielesHacer$bin == tOSmielesHacer, 'binN'] <- 0
mielesHacer$binN <- mielesHacer$bin - tOSmielesHacer
mielesHacer$condToken <- 1

# mieles : di : Transportar
mielesTransportar <- filter(df_dur, condLabel == 'mieles' & verb == 'transportar')
mielesTransportar$binN <- 999
tOSmielesTransportar <- round((unique(mielesTransportar$word4_v1_20ms) - unique(mielesTransportar$word4_c1v1)) / 5)
mielesTransportar[mielesTransportar$bin == tOSmielesTransportar, 'binN'] <- 0
mielesTransportar$binN <- mielesTransportar$bin - tOSmielesTransportar
mielesTransportar$condToken <- 2



# rol : mono : Estrenar
rolEstrenar <- filter(df_dur, condLabel == 'rol' & verb == 'estrenar')
rolEstrenar$binN <- 999
tOSrolEstrenar <- round((unique(rolEstrenar$word4_v1_20ms) - unique(rolEstrenar$word4_c1v1)) / 5)
rolEstrenar[rolEstrenar$bin == tOSrolEstrenar, 'binN'] <- 0
rolEstrenar$binN <- rolEstrenar$bin - tOSrolEstrenar
rolEstrenar$condToken <- 1

# rol : mono : Aceptar
rolAceptar <- filter(df_dur, condLabel == 'rol' & verb == 'aceptar')
rolAceptar$binN <- 999
tOSrolAceptar <- round((unique(rolAceptar$word4_v1_20ms) - unique(rolAceptar$word4_c1v1)) / 5)
rolAceptar[rolAceptar$bin == tOSrolAceptar, 'binN'] <- 0
rolAceptar$binN <- rolAceptar$bin - tOSrolAceptar
rolAceptar$condToken <- 2

# roles : di : Estrenar
rolesEstrenar <- filter(df_dur, condLabel == 'roles' & verb == 'estrenar')
rolesEstrenar$binN <- 999
tOSrolesEstrenar <- round((unique(rolesEstrenar$word4_v1_20ms) - unique(rolesEstrenar$word4_c1v1)) / 5)
rolesEstrenar[rolesEstrenar$bin == tOSrolesEstrenar, 'binN'] <- 0
rolesEstrenar$binN <- rolesEstrenar$bin - tOSrolesEstrenar
rolesEstrenar$condToken <- 1

# roles : di : Aceptar
rolesAceptar <- filter(df_dur, condLabel == 'roles' & verb == 'aceptar')
rolesAceptar$binN <- 999
tOSrolesAceptar <- round((unique(rolesAceptar$word4_v1_20ms) - unique(rolesAceptar$word4_c1v1)) / 5)
rolesAceptar[rolesAceptar$bin == tOSrolesAceptar, 'binN'] <- 0
rolesAceptar$binN <- rolesAceptar$bin - tOSrolesAceptar
rolesAceptar$condToken <- 2



# sol : mono : Descubrir
solDescubrir <- filter(df_dur, condLabel == 'sol' & verb == 'descubrir')
solDescubrir$binN <- 999
tOSsolDescubrir <- round((unique(solDescubrir$word4_v1_20ms) - unique(solDescubrir$word4_c1v1)) / 5)
solDescubrir[solDescubrir$bin == tOSsolDescubrir, 'binN'] <- 0
solDescubrir$binN <- solDescubrir$bin - tOSsolDescubrir
solDescubrir$condToken <- 1

# sol : mono : Escribir
solEscribir <- filter(df_dur, condLabel == 'sol' & verb == 'escribir')
solEscribir$binN <- 999
tOSsolEscribir <- round((unique(solEscribir$word4_v1_20ms) - unique(solEscribir$word4_c1v1)) / 5)
solEscribir[solEscribir$bin == tOSsolEscribir, 'binN'] <- 0
solEscribir$binN <- solEscribir$bin - tOSsolEscribir
solEscribir$condToken <- 2

# soles : di : Descubrir
solesDescubrir <- filter(df_dur, condLabel == 'soles' & verb == 'descubrir')
solesDescubrir$binN <- 999
tOSsolesDescubrir <- round((unique(solesDescubrir$word4_v1_20ms) - unique(solesDescubrir$word4_c1v1)) / 5)
solesDescubrir[solesDescubrir$bin == tOSsolesDescubrir, 'binN'] <- 0
solesDescubrir$binN <- solesDescubrir$bin - tOSsolesDescubrir
solesDescubrir$condToken <- 1

# soles : di : Escribir
solesEscribir <- filter(df_dur, condLabel == 'soles' & verb == 'escribir')
solesEscribir$binN <- 999
tOSsolesEscribir <- round((unique(solesEscribir$word4_v1_20ms) - unique(solesEscribir$word4_c1v1)) / 5)
solesEscribir[solesEscribir$bin == tOSsolesEscribir, 'binN'] <- 0
solesEscribir$binN <- solesEscribir$bin - tOSsolesEscribir
solesEscribir$condToken <- 2



# tul : mono : Coser
tulCoser <- filter(df_dur, condLabel == 'tul' & verb == 'coser')
tulCoser$binN <- 999
tOStulCoser <- round((unique(tulCoser$word4_v1_20ms) - unique(tulCoser$word4_c1v1)) / 5)
tulCoser[tulCoser$bin == tOStulCoser, 'binN'] <- 0
tulCoser$binN <- tulCoser$bin - tOStulCoser
tulCoser$condToken <- 1

# tul : mono : Pedir
tulPedir <- filter(df_dur, condLabel == 'tul' & verb == 'pedir')
tulPedir$binN <- 999
tOStulPedir <- round((unique(tulPedir$word4_v1_20ms) - unique(tulPedir$word4_c1v1)) / 5)
tulPedir[tulPedir$bin == tOStulPedir, 'binN'] <- 0
tulPedir$binN <- tulPedir$bin - tOStulPedir
tulPedir$condToken <- 2

# tules : di : Coser
tulesCoser <- filter(df_dur, condLabel == 'tules' & verb == 'coser')
tulesCoser$binN <- 999
tOStulesCoser <- round((unique(tulesCoser$word4_v1_20ms) - unique(tulesCoser$word4_c1v1)) / 5)
tulesCoser[tulesCoser$bin == tOStulesCoser, 'binN'] <- 0
tulesCoser$binN <- tulesCoser$bin - tOStulesCoser
tulesCoser$condToken <- 1

# tules : di : Pedir
tulesPedir <- filter(df_dur, condLabel == 'tules' & verb == 'pedir')
tulesPedir$binN <- 999
tOStulesPedir <- round((unique(tulesPedir$word4_v1_20ms) - unique(tulesPedir$word4_c1v1)) / 5)
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
df_adj$binAdj <- df_adj$binN - 20
df_adj <- as.data.frame(df_adj)



# check it 
glimpse(df_adj)

# write table
write.table(df_adj, "./mySources/data/clean/durationBIN5Clean.csv", row.names = F, quote = F, sep = ",")

