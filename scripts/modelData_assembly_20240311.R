## Model data assembly script
# Merge necessary datasets (water checm and lake characteristics)
# AJR: 2024-03-11

## Load libraries and functions
library(tidyverse)
library(here)
library(lubridate)

# Function to convert DMS to DD
convert_dms_to_dd <- function(degrees, minutes, seconds, direction) {
  # Calculate decimal degrees
  dd <- degrees + minutes/60 + seconds/3600
  # Adjust for negative direction
  if (direction %in% c("S", "W")) {
    dd <- -dd
  }
  return(dd)
}


## BsM ##
## lake productivity and chemistry ##
#-----------------------------------------#
bsmChemRaw <- read_csv(here("data","waterChem_bsm.csv"))

## lake variables (class, size, etc.) ##
bsmLakes <- read_csv(here("data","waterbodyAndBathymetry_bsm_20220421.csv")) %>% 
  select(waterbodyID,surfaceArea,lakeVolume,maxDepth,meanDepth,lat,long) %>% 
  distinct() 

## Data manipulation; merging data ####
bsmDat <- bsmChemRaw %>% 
  mutate(yearSample = as.numeric(year(dateSample)),
         monthSample = as.numeric(month(dateSample)),
         daySample = as.numeric(day(dateSample)),
         FMZ = as.factor(FMZ)) %>% 
  rename(predictorVariable = chemistryVariable, predictorValue= chemistryValue)

bsmMorphoDat <- bsmChemRaw %>% 
  select(waterbodyID,secchiDepth) %>% 
  group_by(waterbodyID) %>% 
  summarise(secchiDepth = mean(secchiDepth)) %>% 
  full_join(bsmLakes, by="waterbodyID") # %>% 
  #pivot_longer(-c(waterbodyID,secchiDepth), names_to = "predictorVariable", values_to = "predictorValue") 

bsmNutRatio <- bsmChemRaw %>% 
  select(waterbodyID, dateSample,secchiDepth,chemistryVariable,chemistryValue,BsM_Cycle,FMZ,lakeName) %>% 
  filter(chemistryVariable %in% c("TDP","NNTKUR (ug/L)")) %>% 
  pivot_wider(values_from = chemistryValue,names_from = chemistryVariable, values_fn = mean) %>% 
  mutate(NtoPratio = log(`NNTKUR (ug/L)`+0.1)/log(TDP+0.1),
         FMZ = as.factor(FMZ),
         yearSample = as.numeric(year(dateSample))) %>% 
  filter_if(is.numeric, ~. > 0 & . < Inf) %>%  #filter out instances where 0 or infinity are calculated because numerator/denominator missing
  select(-`NNTKUR (ug/L)`, -TDP) %>% 
  pivot_longer(NtoPratio, names_to = "predictorVariable", values_to = "predictorValue") 

cor.test(bsmNutRatio$secchiDepth,bsmNutRatio$predictorValue) #highly significant cor
modNut <- lm(predictorValue~secchiDepth, data=bsmNutRatio)
summary(modNut) #singificant but low R2

allBsm <- bind_rows(bsmDat,bsmNutRatio) %>% 
  left_join(select(bsmMorphoDat,-secchiDepth),by="waterbodyID")

write_csv(allBsm, here("data","BsM_modelData_20240311.csv"))

## ARU database ##
aruDat <- read_csv(here("data","ontarioAquaticHabitatInventory_oldData.csv")) %>% 
  mutate(latD = as.numeric(substr(LATITUDE,1,2)),
         latM = as.numeric(substr(LATITUDE,3,4)),
         latS = as.numeric(substr(LATITUDE,5,6)),
         latDirection = "N",
         lonD = as.numeric(substr(LONGITUDE,1,2)),
         lonM = as.numeric(substr(LONGITUDE,3,4)),
         lonS = as.numeric(substr(LONGITUDE,5,6)),
         lonDirection = "W") %>% 
  mutate(latARU = convert_dms_to_dd(latD, latM, latS, latDirection),
         longARU = convert_dms_to_dd(lonD, lonM, lonS, lonDirection)) %>%  #Convert DMS to dd.dd
  select(-c(latD,latM,latS,lonD,lonM,lonS))

write_csv(aruDat,here("data","ARU_modelData_20240311.csv"))

## Transform data for stats ####
statDat_prim <- allBsm %>% 
  select(BsM_Cycle,FMZ,waterbodyID,lakeName,predictorVariable,predictorValue,dateSample,yearSample) #Select only necessary variables
statDat_long <- allBsm %>% #Need to join secchi depth
  select(BsM_Cycle,FMZ,waterbodyID,lakeName,dateSample,yearSample,secchiDepth,maxDepth,meanDepth,lat,long) %>% 
  distinct() %>% 
  pivot_longer(c(secchiDepth,maxDepth,meanDepth,lat,long), names_to = "predictorVariable", values_to = "predictorValue") %>% 
  bind_rows(statDat_prim) 

#Turn long data to wide data
statDat_wide <-  statDat_long %>% 
  pivot_wider(names_from = predictorVariable, values_from = predictorValue,
              values_fn = mean) %>%  #"values_fn = mean" only necessary because certain observations contain NAs, 
  #though these NAs make sense (for lake lat/longs can't have a BsM cycle 
  #or dateSample, yet those columns are necessary elsewhere). If this command
  # is not included dataframe becomes a list, which we don't want
  rename(TKN = `NNTKUR (ug/L)`) %>% #Rename, make easier to understand
  filter(!is.na(DOC) & 
           !is.na(secchiDepth) &
           !is.na(TDP) &           #Remove rows where NAs exist. NAs fuck up models
           !is.na(maxDepth) & 
           !is.na(DOC) & 
           !is.na(TKN),
         !is.na(NtoPratio)) # Need to remove NAs from predictor variables

write_csv(statDat_wide, here("data","BsM_WideSampleModelData_20240311.csv"))
