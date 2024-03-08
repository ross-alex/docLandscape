## Is there a relationship between secchi and DOC (and potentially other variables)
# AJR: 2024-02-17

library(tidyverse)
library(here)
library(lubridate)

## Read Data ####

## BsM ##
## lake productivity and chemistry ##
#-----------------------------------------#
bsmChemRaw <- read_csv(here("data","waterChem_bsm.csv"))
  
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
         longARU = convert_dms_to_dd(lonD, lonM, lonS, lonDirection)) %>%  #Convert DMS to dd.dd
  select(-c(latD,latM,latS,lonD,lonM,lonS))
  
## Data manipulation + primary visualizations ####
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
  inner_join(bsmLakes, by="waterbodyID") %>% 
  pivot_longer(-c(waterbodyID,secchiDepth), names_to = "predictorVariable", values_to = "predictorValue") 

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

allBsm <- bind_rows(bsmDat,bsmMorphoDat,bsmNutRatio) 

## Exploratory plots - what correlates with secchi depth
pRough_allVars <- ggplot(allBsm) +
  geom_point(aes(x=secchiDepth,y=predictorValue), alpha = 0.3)+
  facet_wrap(~predictorVariable, scales = "free")
pRough_allVars

pRough_allVars_log <- ggplot(allBsm) +
  geom_point(aes(x=secchiDepth,y=log(predictorValue), alpha = 0.3))+
  facet_wrap(~predictorVariable, scales = "free")
pRough_allVars_log #Log DOC, TDP and N look associated with secchi depth
# It looks like secchi depth does a good job predicting DOC. Other possibly useful 
# correlations w. secchi include Colour, max/meanDepth, TDP, NNTKUR

## Now let's look at DOC as a response variable along with all other variables
docRespDat <- allBsm %>% 
  filter(predictorVariable == "DOC") %>% 
  rename(DOC = predictorValue) %>% 
  select(waterbodyID, dateSampleStart, DOC) %>% 
  right_join(allBsm, by=c("waterbodyID"))

pDOCresp <- ggplot(docRespDat) +
  geom_point(aes(y=DOC,x=(predictorValue)), alpha = 0.3)+
  facet_wrap(~predictorVariable, scales = "free")
pDOCresp

pDOCresp_log <- ggplot(docRespDat) +
  geom_point(aes(y=log(DOC),x=log(predictorValue)), alpha = 0.3)+
  facet_wrap(~predictorVariable, scales = "free")
pDOCresp_log
# Again, DOC, TDP, NNTKUR, depth and Colour look like important independent predictors of DOC

## Goal 1: What is the best predictive DOC model using contemporary data ####
## Transform data for stats ##
statDat_prim <- allBsm %>% 
  select(BsM_Cycle,FMZ,waterbodyID,lakeName,predictorVariable,predictorValue,dateSample,yearSample)
statDat_long <- bsmChemRaw_expl %>% #Need to join secchi depth
  select(BsM_Cycle,FMZ,waterbodyID,lakeName,dateSample,yearSample,secchiDepth) %>% 
  distinct() %>% 
  pivot_longer(secchiDepth, names_to = "predictorVariable", values_to = "predictorValue") %>% 
  bind_rows(statDat_prim)  
statDat_wide <-  statDat_long %>% 
  pivot_wider(names_from = predictorVariable, values_from = predictorValue,
              values_fn = mean) %>%  #in cases where duplicate chem values exist, take mean
  left_join(bsmLakes, by="waterbodyID") %>% 
  mutate(lat = if_else(is.na(lat.x), lat.y, lat.x),
         long = if_else(is.na(long.x), long.y, long.x),
         maxDepth = if_else(is.na(maxDepth.x), maxDepth.y, maxDepth.x),
         meanDepth = if_else(is.na(meanDepth.x), meanDepth.y, meanDepth.x),
         surfaceArea = if_else(is.na(surfaceArea.x), surfaceArea.y, surfaceArea.x),
         lakeVolume = if_else(is.na(lakeVolume.x), lakeVolume.y, lakeVolume.x)) %>% 
  select(-lat.x, -lat.y, -long.x, -long.y, -maxDepth.y, -maxDepth.x, -meanDepth.y, -meanDepth.x,
         -surfaceArea.y, -surfaceArea.x, -lakeVolume.y, -lakeVolume.x) %>% 
  rename(TKN = `NNTKUR (ug/L)`) %>% 
  filter(!is.na(DOC) & 
         !is.na(secchiDepth) &
         !is.na(TDP) & 
         !is.na(maxDepth) & 
         !is.na(DOC) & 
         !is.na(TKN),
         !is.na(NtoPratio)) # Need to remove NAs from predictor variables
statDat_wide1 <- statDat_wide[-c(1480,102,257),] #From initial modeling, these points represent outliers

# What duplicates exist in waterbodyID?
duplicates <- statDat_wide %>%
  filter(duplicated(waterbodyID) | duplicated(waterbodyID, fromLast = TRUE))

#Make stat data but with means across different BsM cycles
statDat_means <- statDat_wide1 %>%
  group_by(waterbodyID,) %>% 
  select(-BsM_Cycle, -dateSample, -yearSample) %>%
  summarise_all(mean, na.rm = TRUE)

##NOTE: Things to think about for analysis
# 1) Choose whether looking at all lakes or only non-eutrophic lakes. What is 'non-eutrophic'? Most/all guidance is for TP, not TDP
# 2) Remove outliers: Do this before modelling, look at plots and remove data accordingly 
# 3) Remove outliers: In diagnostic plots, remove data as needed
# 4) Come up with modelling plan. Full model w. 4-way interaction or more biologically informed? 4-way is significant...that's hard
# 4) Modelling plan: Leave FMZ out of prediction model; both as fixed and random effect (so that we can predict differences across FMZs later)
# 5) Choose a subset of data? Training vs. testing?

## Create a predictive models across all BsM Lakes
## Model selection using dredge (IC)
modFull_4 <- lm(log(DOC)~secchiDepth*log(TDP+0.1)*maxDepth*log(TKN+0.1),
                data = statDat_means, na.action="na.fail")
modSel_4 <-dredge(modFull_4)

modFull_3 <- lm(log(DOC)~secchiDepth*NtoPratio*maxDepth,
              data = statDat_means, na.action="na.fail")
modSel_3 <- dredge(modFull_3)


#Top model - 4 vars
fullMod_4 <- lm(log(DOC)~secchiDepth*log(TDP+0.1)*maxDepth*log(TKN+0.1),
                data = statDat_means)
plot(fullMod_4)
summary(fullMod_4)

#Top model - 3 vars
fullMod_3 <- lm(log(DOC)~secchiDepth*NtoPratio*maxDepth,
                data = statDat_means)
plot(fullMod_3)
summary(fullMod_3)
car::vif(fullMod_3)

fullMod_3_noThreeWay <- lm(log(DOC)~secchiDepth+NtoPratio+maxDepth+
                  secchiDepth:NtoPratio+secchiDepth:maxDepth+NtoPratio:maxDepth,
                data = statDat_means)
plot(fullMod_3_noThreeWay)
summary(fullMod_3_noThreeWay)
car::vif(fullMod_3_noThreeWay) # high vifs

fullMod_3_intReduce <- lm(log(DOC)~secchiDepth+NtoPratio+maxDepth+
                             secchiDepth:NtoPratio+secchiDepth:maxDepth,
                           data = statDat_means)
plot(fullMod_3_intReduce)
summary(fullMod_3_intReduce)
car::vif(fullMod_3_intReduce) #high vifs

fullMod_3_intReduce2 <- lm(log(DOC)~secchiDepth+NtoPratio+maxDepth+
                            secchiDepth:maxDepth,
                          data = statDat_means)
plot(fullMod_3_intReduce2)
summary(fullMod_3_intReduce2)
car::vif(fullMod_3_intReduce2) # interaction vif still > 5 (8.22); strictly additive model best

# Purely additive model
addMod_4 <- lm(log(DOC)~secchiDepth+log(TDP+0.1)+log(TKN+0.1)+maxDepth, data=statDat_means)
plot(addMod_4) 
summary(addMod_4)
car::vif(addMod_4)

addMod_3 <- lm(log(DOC)~secchiDepth+NtoPratio+maxDepth, data=statDat_means)
plot(addMod_3) 
summary(addMod_3)
car::vif(addMod_3)

# What about just the ARU variable set, secchi and maxDepth
aruVarsMod <- lm(log(DOC)~secchiDepth*maxDepth, data=statDat_means, na.action = "na.fail")
aruVarsSel <- dredge(aruVarsMod) #interaction is important

aruVars_finalMod <- lm(log(DOC)~secchiDepth*maxDepth, data=statDat_means)
plot(aruVars_finalMod)
summary(aruVars_finalMod)
car::vif(aruVars_finalMod) #interaction still kind of high (8.14)


## Goal 2: How have DOC conc. changed over time ####

# Visualize ARU data
aruPlotDat <- aruDat %>% 
  select(ARU_NM,LATITUDE,LONGITUDE,SRFAREA,PERIM,DEPMAX,DEPMN, SECCHI, TDS, MEI) %>% 
  pivot_longer(-ARU_NM, names_to = "predictorVariable", values_to = "predictorValue")
  
aruP <- ggplot(aruPlotDat)+
  geom_histogram(aes(x=predictorValue)) +
  facet_wrap(~predictorVariable, scales = "free")
aruP

# With ARU data, estimate TP from TDS and DOC
aruDOC_coefs <- coef(aruVars_finalMod)

aruDat_estimated <- aruDat %>% 
  filter(TDS < 300) %>% #Different linear relationship for saline lakes in Chow-Fraser (1991)
  mutate(logTP = 0.95*log(TDS)-0.669, #This formula is from Chow-Fraser (1991)
         TP = exp(logTP),
         lakeTrophicStatus = case_when(TP < 10 ~ "oligotrophic",
                                       TP > 10 & TP < 20 ~ "mesotrophic",
                                       TP > 20 ~ "eutrophic"),
         logDOC = aruDOC_coefs[1] +
           aruDOC_coefs[2]*SECCHI +
           aruDOC_coefs[3]*DEPMAX +
           aruDOC_coefs[4]*SECCHI*DEPMAX, #Use coefs from reduced BsM data model to estimate ARU DOC
         estDOC = exp(logDOC))

## Create dataset that compares contemporary and historic data
compDat <- 

## Data visualization

## Data filtering - do we remove eutrophic lakes (see Mike email, Cooney email chain)
