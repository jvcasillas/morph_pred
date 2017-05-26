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



df_lex$TOadj <- round(df_lex$targetOffset / 20)
df_lex$binN <- 999
df_lex[df_lex$bin %in% unique(df_lex$TOadj), 'binN'] <- 0
df_lex$binN <- df_lex$bin - df_lex$TOadj


# 200ms bin adjustment for VWP
df_lex$binAdj <- df_lex$binN - 20
df_lex <- as.data.frame(df_lex)



# check it 
glimpse(df_lex)

# write table
write.table(df_lex, "./mySources/data/clean/lexicalBIN10Clean.csv", row.names = F, quote = T, sep = ",")

