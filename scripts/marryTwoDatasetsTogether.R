## Marrying ARU and BsM databases together
# 2024-03-05; AJR; KB, WD

# To compare contemporary vs. historic DOC concentrations, we must first determine what lakes are found in both databases.
# While ARU_LID and BsM waterbodyID appear to link up, there are many lakes with common names (i.e., Long L., Clear L., etc.)
# where it is worth making sure that we are pairing the right lakes together. First, lakeName character strings were cleaned to 
# remove prefixes or suffixs, and then taken into Excel to manually clean/observe lakeNames to be a consistent format between the two datasets.
# Then, Wade Dombroski and Kennedy Bucci importated both datasets into GIS software and conducted a spatial join that was anchored on
# the BsM dataset. Lakes were then given a common ID to pair contemporary/historic data together. In total, 869 lakes matched between the two 
# datasets. It should be noted that certain BsM 'trend' lakes have repeated samples across lake-years, as the purpose of these lakes is 
# assessing repeatability and trends in sampled data.

## Load libraries and functions ####
library(tidyverse)
library(here)
library(lubridate)
library(lme4)
library(lmerTest)
library(MuMIn)

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

## Read Data ####

## BsM ##
## lake productivity and chemistry ##
#-----------------------------------------#
bsmChemRaw <- read.csv("/Users/alexross/Documents/Lakehead/lakeTrout_lifeHistory/bsmData/fromBlairMNRF/waterChem_bsm.csv",header=TRUE,stringsAsFactors=FALSE)

## lake variables (class, size, etc.) ##
bsmLakes <- read_csv(here("data","waterbodyAndBathymetry_bsm_20220421.csv")) %>% 
  select(waterbodyID,surfaceArea,lakeVolume,maxDepth,meanDepth,lat,long) %>% 
  distinct() 

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
         longARU = convert_dms_to_dd(lonD, lonM, lonS, lonDirection)) #Convert DMS to dd.dd

# Find harmonized dataset for ARU and BsM data ####

bsmDatIDs <- read_csv(here("data","waterbodyAndBathymetry_bsm_20220421.csv")) %>% 
  select(waterbodyID, lakeName, lat, long) %>% 
  mutate(latD = substr(lat,1,2),
         latM1 = substr(lat,4,4),
         latM2 = substr(lat, 5,5),
         longD = substr(long,2,3),
         longM1 = substr(long,5,5),
         longM2 = substr(long, 6,6),
         lakeName = str_replace_all(lakeName, "Lake", ""),
         lakeName = str_trim(lakeName, side = "left")) %>% #remove isntances where first character string is " "
  distinct()

# Some lakes have more than one word associated with them; extract longest word as matching name
bsMDatID_simplified <- bsmDatIDs %>% 
  mutate(waterbodyAndLake = paste(waterbodyID,lakeName,sep=" ")) %>% 
  ungroup() %>% 
  separate_wider_delim(waterbodyAndLake, delim = " ", 
                       names = c("waterbodyID_OG","lakeName_a","lakeName_b"), 
                       too_many = "debug",
                       too_few = "debug")
#write_csv(bsMDatID_simplified, here("data","bsmLakeID_simplify.csv")) #export, manually change wonky multiple word lakes, bring back in

aruDatIDs <- aruDat %>% 
  select(ARU_NM, ARU_LID, LATITUDE, LONGITUDE) %>% 
  rename(lakeName = ARU_NM) %>% 
  mutate(latD = substr(LATITUDE,1,2),
         latM1 = substr(LATITUDE,3,3),
         latM2 = substr(LATITUDE,4,4),
         longD = substr(LONGITUDE,1,2),
         longM1 = substr(LONGITUDE,3,3),
         longM2 = substr(LONGITUDE,4,4),
         lakeName = str_replace_all(lakeName, " L\\.", ""),
         lakeName = str_trim(lakeName, side = "left"))

aruDatID_simplified <- aruDatIDs %>% 
  mutate(waterbodyAndLake = paste(lakeName,sep=" ")) %>% 
  ungroup() %>% 
  separate_wider_delim(waterbodyAndLake, delim = " ", 
                       names = c("lakeName_a","lakeName_b"), 
                       too_many = "debug",
                       too_few = "debug")
#write_csv(aruDatID_simplified, here("data","aruLakeID_simplify.csv")) #export, manually change wonky multiple word lakes, bring back in

## What's been done in Excel with aruDatID_simplified and bsmDatID_simplified
# 1) 570 entries deleted for incomplete lake names(eg., names starting w numbers that will never get matched up)
# 2) 100 entries deleated because numbers at end of lakeName (indicate replicate samples, but not sure how should be treated)
# 3) 128 entries deleted for strange lake suffixs (i.e B., C. I., P., gas station...)
# All lakeNames for both the aru and bsm data have then been manually scanned to try and catch any odd things
# At this point, lakeNames are as best harmanized as possible between the two datasets. 
# Bring them back into R, and merge to find what lakes have been repeated

## Read ARU and BsM lakes back, merge
aruMergeDat <- read_csv(here("data","aruLakeID_simplifiedCleaned_20240304.csv")) %>% 
  select(-lakeName_a,-lakeName_b,-waterbodyAndLake,-waterbodyAndLake_ok,-waterbodyAndLake_pieces,-waterbodyAndLake_remainder) %>% 
  rename(latARU=LATITUDE,longARU=LONGITUDE,
         latD_aru=latD,latM1_aru=latM1,latM2_aru=latM2,longD_aru=longD,longM1_aru=longM1,longM2_aru=longM2)
bsmMergeDat <- read_csv(here("data","bsmLakeID_simplifiedCleaned_20240304.csv")) %>% 
  select(-waterbodyID_OG,-lakeName_a,-lakeName_b,-waterbodyAndLake,-waterbodyAndLake_ok,-waterbodyAndLake_pieces,-waterbodyAndLake_remainder) %>% 
  rename(bsmWaterbodyID=waterbodyID,latBsM=lat,longBsM=long,
         latD_BsM=latD,latM1_BsM=latM1,latM2_BsM=latM2,longD_BsM=longD,longM1_BsM=longM1,longM2_BsM=longM2)

mergeLakeIDs <- aruMergeDat %>% 
  inner_join(bsmMergeDat, by="lakeName", multiple = "all") %>% 
  arrange(lakeName,latD_aru) %>% 
  relocate(lakeName,latARU,latBsM,longARU,longBsM,
           latD_aru,latM1_aru,latM2_aru,latD_BsM,latM1_BsM,latM2_BsM,
           longD_aru,longM1_aru,longM2_aru,longD_BsM,longM1_BsM,longM2_BsM,
           ARU_LID,bsmWaterbodyID)

#write_csv(mergeLakeIDs, here("data","lakeIDsToMerge_doc_20240305.csv"))

## 2024-03-05; made mistake that ARU and BsM coordinate systems were different. To fix, read mergeLakeIDs file and then change
# ARU coordinates with conversion below. Lots of manual cleaning was done from above data objects on lakeNames 
# so don't want to go through that entire process again (ARU was in DMS, needs to be dd.dd)
mergeLakeID_convLL <- read_csv(here("data","lakeIDsToMerge_doc_20240305.csv")) %>% 
  mutate(latD_aru = as.numeric(substr(latARU,1,2)),
         latM_aru = as.numeric(substr(latARU,3,4)),
         latS_aru = as.numeric(substr(latARU,5,6)),
         latDirection = "N",
         lonD_aru = as.numeric(substr(longARU,1,2)),
         lonM_aru = as.numeric(substr(longARU,3,4)),
         lonS_aru = as.numeric(substr(longARU,5,6)),
         lonDirection = "W") %>% 
  mutate(latARU_conv = convert_dms_to_dd(latD_aru, latM_aru, latS_aru, latDirection),
         longARU_conv = convert_dms_to_dd(lonD_aru, lonM_aru, lonS_aru, lonDirection)) %>%  #Convert DMS to dd.dd
  relocate(lakeName,latARU_conv,latBsM, longARU_conv,longBsM,
           latD_aru,latM1_aru,latM2_aru,latD_BsM,latM1_BsM,latM2_BsM,
           longD_aru,longM1_aru,longM2_aru,longD_BsM,longM1_BsM,longM2_BsM,
           ARU_LID,bsmWaterbodyID)
#write_csv(mergeLakeID_convLL, here("data","mergeLakeIDdata_properlyConverted_20240305.csv"))

## Read in GIS spatial join dataset ####
# Read data, and create a common lakeID for both datasets
matchedData <- read_csv(here("data","BsM_ARU_matches_20240306.csv")) %>% 
  group_by(BsM_lakeName) %>% 
  mutate(commonID = paste(BsM_lakeName,row_number(),sep = "_")) %>% 
  ungroup()

# Now that ARU lat/lons have been put into dd.dd the LID actually looks like it has the common format of BsM. See overlap, good double check
commonLakeIDs <- matchedData %>% 
  mutate(duplicateLakeID = if_else(bsmWaterbodyID == ARU_LID, "Y", "N"))
# ~100 don't because 1 or 2 digits off, therefore can't use ARU_LID and bsmWaterbodyID as commonID, use the 'commonID' created in code
# chunk above

#write_csv(matchedData, here("data","commonDataIDs_aruToBsm_20240307.csv"))


