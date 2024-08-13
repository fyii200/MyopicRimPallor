######################################################
##                   Fabian Yii                     ##
##               fabian.yii@ed.ac.uk                ##
##    Link SG's dataframe ID to FY's dataframe ID   ##
######################################################

rm(list=ls())
library(dplyr)
library(tidyr)
library(readxl)

########## Part 1: identify FY's participant ID corresponding to each SG's participant ID ##########
## Set up
# # Move up one directory to the and set it as working directory 
# setwd('..') 
# Read data
dSG <- read_excel(paste0('data', .Platform$file.sep, 'rawDataSamGibbon.xlsx'))
dFY <- read_excel(paste0('data', .Platform$file.sep, 'rawDataFabianYii.xlsx'))
# Replace NA with 0
dSG <- dSG %>% replace_na(list(`Townsend deprivation index at recruitment`=0,
                               `Intra-ocular pressure, Goldmann-correlated (left) | Instance 0`=0,
                               `Intra-ocular pressure, Goldmann-correlated (right) | Instance 0`=0,
                               `QC - Image quality (left) | Instance 0`=0,
                               `Spherical power (right) | Instance 0 | Array 0`=0,
                               `Spherical power (left) | Instance 0 | Array 0`=0,
                               `Cylindrical power (right) | Instance 0 | Array 0`=0,
                               `Cylindrical power (left) | Instance 0 | Array 0`=0))
dFY <- dFY %>% replace_na(list(`Townsend deprivation index at recruitment`=0,
                               `Intra-ocular pressure, Goldmann-correlated (left) | Instance 0`=0,
                               `Intra-ocular pressure, Goldmann-correlated (right) | Instance 0`=0,
                               `QC - Image quality (left) | Instance 0`=0,
                               `Spherical power (right) | Instance 0 | Array 0`=0,
                               `Spherical power (left) | Instance 0 | Array 0`=0,
                               `Cylindrical power (right) | Instance 0 | Array 0`=0,
                               `Cylindrical power (left) | Instance 0 | Array 0`=0))

## Link SG's participant ID to FY's participant ID
dFY$`Participant ID SG` <- NA
dFY <- dFY %>% relocate('Participant ID SG', .after='Participant ID')
for(i in 1:nrow(dFY)){
  # Get FY's participant ID in the current loop
  FYid <- dFY$`Participant ID`[i]
  # Identify the row in SG's dataframe that contains the values as the variables in the current row of FY's dataframe
  ind <- which(dSG$`Age at recruitment` == dFY$`Age at recruitment`[i] &
               dSG$`Townsend deprivation index at recruitment` ==  dFY$`Townsend deprivation index at recruitment`[i] &
               dSG$`Intra-ocular pressure, Goldmann-correlated (left) | Instance 0` == dFY$`Intra-ocular pressure, Goldmann-correlated (left) | Instance 0`[i] &
               dSG$`Intra-ocular pressure, Goldmann-correlated (right) | Instance 0` == dFY$`Intra-ocular pressure, Goldmann-correlated (right) | Instance 0`[i] &
               dSG$`QC - Image quality (left) | Instance 0` == dFY$`QC - Image quality (left) | Instance 0`[i] &
               dSG$`Spherical power (right) | Instance 0 | Array 0` == dFY$`Spherical power (right) | Instance 0 | Array 0`[i] &
               dSG$`Spherical power (left) | Instance 0 | Array 0` == dFY$`Spherical power (left) | Instance 0 | Array 0`[i] &
               dSG$`Cylindrical power (right) | Instance 0 | Array 0` == dFY$`Cylindrical power (right) | Instance 0 | Array 0`[i] &
               dSG$`Cylindrical power (left) | Instance 0 | Array 0` == dFY$`Cylindrical power (left) | Instance 0 | Array 0`[i] )
  if(length(ind)==1){
    # Extract SG's participant ID
    SGid <- dSG$subjectID[ind]
    # Save SG's participant ID to the row corresponding to the same participant in FY's dataframe
    dFY$`Participant ID SG`[i] <- SGid
    } else { print(paste('FY id:', FYid, 'has', length(ind), 'match(s)')) }
}
## Save linkage dataframe as csv
names(dFY)[1] <- 'FYid'
names(dFY)[2] <- 'SGid'
write.csv(dFY[,1:2], 'data/link.csv', row.names=FALSE)
######################################## Part 1: COMPLETE ########################################



################# Part 2: link SG's pallor data to FY's main (cleaned) dataframe #################
## Set up
# Read data
d        <- read.csv(paste0('data', .Platform$file.sep, 'cleaned_data_long_all.csv')) # main dataframe
link     <- read.csv(paste0('data', .Platform$file.sep, 'link.csv'))                  # linkage dataframe
pallorRE <- read_excel(paste0('data', .Platform$file.sep, 'SGpallorRE.xlsx'))         # right eye pallor data
pallorLE <- read_excel(paste0('data', .Platform$file.sep, 'SGpallorLE.xlsx'))         # left eye pallor data
# Create new columns in the main dataframe 'd' to store pallor-related metrics
d$palProcError   <- NA
d$palT           <- NA
d$palTI          <- NA
d$palNI          <- NA
d$palN           <- NA
d$palNS          <- NA
d$palTS          <- NA
d$palPMB         <- NA
d$palG           <- NA
d$palNTratio     <- NA
d$palDisc        <- NA
d$palCrowdedness <- NA
d$palDD          <- NA
d$palDM          <- NA
d$palWobbliness  <- NA
d$palTorsion     <- NA
d$palContInt     <- NA
d$palTorsion     <- NA
d$palRejectThres <- NA
# Add pallor data to the main dataframe 'd'
for(i in 1:nrow(d)){
  # Only proceed if fundus is present
  if(d$fundus[i]!=""){
    # Extract id and laterality
    FYid  <- d$id[i]
    FYeye <- d$eye[i]
    # Start adding
    SGid  <- link[which(link$FYid==FYid), ]$SGid
    if(FYeye=="RE"){
      d[i, 79:96] <- pallorRE[which(pallorRE$subjectID==SGid), c(2, 4:13, 18:20, 22:23, 26, 28)]
    } else {
      d[i, 79:96] <- pallorLE[which(pallorLE$subjectID==SGid), c(2, 4:13, 18:20, 22:23, 26, 28)] }
  }
}
## Save linked, cleaned dataframe as csv
write.csv(d, 'data/linked_cleaned_data_long_all.csv', row.names=FALSE)
######################################## Part 2: COMPLETE ########################################
















