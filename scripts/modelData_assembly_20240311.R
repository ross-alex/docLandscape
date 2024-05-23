## Model data assembly script
# Merge necessary datasets (water checm and lake characteristics)
# AJR: 2024-03-11

## Load libraries and functions
library(tidyverse)
library(here)
library(lubridate)

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

#write_csv(allBsm, here("data","BsM_modelData_20240311.csv"))

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
  filter(!is.na(LATITUDE) | !is.na(LATITUDE)) %>% 
  mutate(latARU = latD + (latM/60) + (latS/60),          #N; DMS to dd.dd
         longAur = -lonD + (lonM/60) + (lonS/60)) %>%    #W; DMS ot dd.dd (need neg for W)
  select(-c(latD,latM,latS,lonD,lonM,lonS))

#write_csv(aruDat,here("data","ARU_modelData_20240311.csv"))

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

#write_csv(statDat_wide, here("data","BsM_WideSampleModelData_20240311.csv"))

## UPDATED BSM DATA ####
# March 19, 2024
# Blair Waslynko provided BsM Cycle 3 data. join to the datafile BsM_modelData_20240311.csv

# Read data
allBsM_update <- read_csv(here("data","BsM_modelData_20240311.csv"))
cyc3Dat <- read_csv(here("data","waterChem_bsm_cycle3.csv"))
cyc3Lakes <- read_csv(here("data","cyc3_lakesBathy.csv")) %>% 
  select(waterbodyID,surfaceArea,lakeVolume,maxDepth,meanDepth,lat,long,yearSample,monthSample,daySample) %>% 
  distinct() 

# Filter cyc3Dat to lakes which already exist in allBsM_update
cyc3Dat_filtered <- cyc3Dat %>% 
  filter(waterbodyID %in% unique(allBsM_update$waterbodyID)) %>% 
  left_join(cyc3Lakes, by="waterbodyID") %>% 
  rename(BsM_Cycle=Cycle,lakeName=`Lake Name`,WaterChemLat=Latitude,WaterChemLong=Longitude) %>% 
  relocate(`COLTR (TCU)`:TDP, .after = last_col()) %>% 
  mutate(FMZ = as.factor(FMZ))
cyc3_long <- cyc3Dat_filtered %>% 
  pivot_longer(cols = `COLTR (TCU)`:TDP, names_to = "predictorVariable", values_to = "predictorValue")

# Merge with old data
allBsM_new <- bind_rows(allBsM_update,cyc3_long)
allBsM_new

#write_csv(allBsM_new, here("data","allBsM_updatedCyc3_20240520.csv"))

## Now updated statDat_wide
statDat_wideCyc3 <- cyc3Dat_filtered %>% 
  rename(TKN = `NNTKUR (ug/L)`) %>% 
  select(-WaterChemLat,-WaterChemLong,-monthSample,-daySample) %>% 
  bind_rows(statDat_wide)

#write_csv(statDat_wideCyc3, here("data","statDatWide_updatedCyc3_20240520.csv"))
