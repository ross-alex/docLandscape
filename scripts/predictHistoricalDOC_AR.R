## Is there a relationship between secchi and DOC (and potentially other variables)
# AJR: 2024-02-17

library(tidyverse)
library(here)
library(lubridate)
library(MuMIn)
library(lme4)
library(lmerTest)

## Read Data ####
allBsm <- read_csv(here("data","BsM_modelData_20240311.csv"))
allARU <- read_csv(here("data","ARU_modelData_20240311.csv"))
statDat_wide <- read_csv(here("data","BsM_WideSampleModelData_20240311.csv")) %>% #This dataset is slightly reduced
                  filter(log(TKN) >4) #outliers <4                            #to remove all rows that had at 
                                                                              #at least one NA. See 'modelData_assembly_20240311.R'
                                                                              # script for details (bottom chunk)
combIDs <- read_csv(here("data","commonDataIDs_aruToBsm_20240307.csv")) #has common IDs between BsM and ARU

#Estimate TP from TDS and apply a trophic status
aruDat <- allARU %>% 
  left_join(select(combIDs,ARU_LID,commonID)) %>% 
  mutate(logTP = 0.95*log(TDS)-0.669, #This formula is from Chow-Fraser (1991)
         TP = exp(logTP),
         lakeTrophicStatus = case_when(TP < 10 ~ "oligotrophic",
                                       TP > 10 & TP < 20 ~ "mesotrophic",
                                       TP > 20 ~ "eutrophic")) %>% 
  filter(!is.na(commonID)) %>% 
  rename(aruSecchi = SECCHI)

#Read all secchi depths
secDat <- read_csv(here("data","summerSecchis_bsm.csv")) %>% 
  select(WbyLID,`Water Chem Sample Date`, `Water Chem Secchi Depth (Spring)`,`Limno Sample Date`, `Limno Secchi Depth (Summer)`) %>% 
  rename(waterbodyID = WbyLID, 
         waterChem_dateSample = `Water Chem Sample Date`, 
         limno_dateSample = `Limno Sample Date`,
         springSecchi = `Water Chem Secchi Depth (Spring)`, 
         summerSecchi = `Limno Secchi Depth (Summer)`) %>% 
  relocate(limno_dateSample, .after = waterChem_dateSample) %>% 
  left_join(select(combIDs,bsmWaterbodyID, commonID), by = c("waterbodyID" = "bsmWaterbodyID")) 

secDat_means <- secDat %>% 
  group_by(waterbodyID,commonID) %>% 
  summarise(meanSummerSecchi = mean(summerSecchi),
            meanSpringSecchi = mean(springSecchi))
## Exploratory plots ####
# what correlates with secchi depth
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

## Goal 1: Checking for differences in water clarity; Predicting historical DOC; model selection ####
## The ARU database has secchi and depth, which can be used to predict DOC. 
# A regression can be applied to BsM data where DOC is known, and the parameters from this model can
# then be applied to the ARU database. DOC distribution seems somewhat gamma distributed. Check models
# for linear, log-linear, and gamma dist. for response variable (DOC)
# Additionally, check if removing eutrophic lakes is necessary. 

## Make dataframe using mean values across different BsM cycles
statDat_means <- statDat_wide %>%
  group_by(waterbodyID) %>% 
  select(-BsM_Cycle, -dateSample, -yearSample) %>%
  summarise_all(mean, na.rm = TRUE) %>% 
  left_join(select(combIDs, bsmWaterbodyID,commonID), by = c("waterbodyID" = "bsmWaterbodyID"))

## Step 1: Are there differences in water clarity between the two times
waterClarDiff <- secDat_means %>% 
  select(waterbodyID,commonID,meanSummerSecchi) %>% 
  left_join(select(aruDat, commonID, aruSecchi, lakeTrophicStatus)) %>% 
  filter_all(all_vars(!is.na(.))) %>% 
  mutate(secchiDiff = meanSummerSecchi-aruSecchi) %>% 
  mutate(lakeTrophicStatus = factor(lakeTrophicStatus, levels = c("oligotrophic","mesotrophic","eutrophic")))

waterClar_addedVars <- waterClarDiff %>% 
  left_join(select(aruDat, commonID, latARU,longARU, DEPMAX))

# T test
(clarTtest <- t.test(waterClarDiff$aruSecchi,waterClarDiff$meanSummerSecchi)) #Significant. ARU has deeper secchis = clearer lakes

# Visualize it
boxplot(waterClarDiff$aruSecchi,waterClarDiff$meanSummerSecchi, names = c("Historic","Contemporary"),
        ylab="Secchi Depth (m)", las = 1)
pClar <- ggplot(waterClarDiff) +
  geom_density(aes(x=secchiDiff),fill="grey",alpha=0.5) +
  geom_vline(xintercept = 0, linetype = "dotted", color = "red") +
  ggtitle("Negative values indicate a lake is darker now than historically") +
  theme_bw()
pClar #Mean > 0, but majority data <0. These are using all data. Let's look at trophic status

# What happens in each trophic group? Oligotrophic most negative...
pClar_trophicStatus <- pClar +
  facet_wrap(~lakeTrophicStatus)
pClar_trophicStatus  
# Oligotrophic lakes getting darker, meso also; eutrophic getting lighter

# Let's model it, see if there are variables driving differences
clarMod <- lm(secchiDiff~lakeTrophicStatus+latARU+longARU+DEPMAX, data=waterClar_addedVars)
summary(clarMod)
anova(clarMod) # Trophic status quite important; Northerly lakes are getting darker

## OK - this tells us that water clarity is changing in ways we predict. Now let's model start the DOC modelling process

## Water clarity only increasing in eutrophic lakes, choose oligo and mesotrophic lakes
## Make dataframe using mean values across different BsM cycles
statDat_means_noEutro <- statDat_means %>% 
  filter(TDP < 10) # Wetzel says 20ugL TP and less is meso/oligo; also stats that TDP is ~50% of TP

## Step 2: Create a model from BsM data
## Look at distribution of DOC observations; what distribution is most appropriate
# Linear
plot_DOCdensity <- ggplot(data = statDat_means_noEutro) +
  geom_density(aes(x=DOC), fill="grey", alpha = 0.7) +
  theme_bw()
plot_DOCdensity #right skew

# log-linear
plot_DOCdensity_log <- plot_DOCdensity +
  scale_x_continuous(trans = "log")
plot_DOCdensity_log #left skew

# gamma distribution
# Fit gamma distribution to the observed data
fit <- MASS::fitdistr(statDat_means_noEutro$DOC, densfun = "gamma")

# Visualize the fit
ggplot(statDat_means_noEutro, aes(x = DOC)) +
  geom_histogram(aes(y = ..density..), bins = 30, fill = "skyblue", color = "black", alpha = 0.7) +
  stat_function(fun = dgamma, args = list(shape = fit$estimate["shape"], rate = fit$estimate["rate"]), color = "red", size = 1) +
  labs(x = "DOC", y = "Density") +
  ggtitle("Gamma Distribution Fit to Data") 
# gamma fits well with shape = 4.91, rate = 0.79

## Model ARU data using GLM, gamma distribution
# Fit the initial model
ARUmod_gamma <- glm(DOC ~ secchiDepth * maxDepth, data = statDat_means_noEutro, family = Gamma(link = "log"), na.action = "na.fail")
ARUmodSel <- dredge(ARUmod_gamma) #interaction is best model

ARUmod_gammaBest <- glm(DOC ~ secchiDepth + maxDepth, data = statDat_means, family = Gamma(link = "log"))

#  Model diagnostics
plot(ARUmod_gammaBest)
# Diagnositcs not great but probably fine

# model summaries
summary(ARUmod_gammaBest)

## What about  straight-up linear models
aruVarsMod <- lm(log(DOC)~secchiDepth*maxDepth, data=statDat_means_noEutro, na.action = "na.fail")
aruVarsSel <- dredge(aruVarsMod) #interaction is not important important

aruVars_finalMod <- lm(log(DOC)~secchiDepth+maxDepth, data=statDat_means)
plot(aruVars_finalMod) #diagnostics slightly better
summary(aruVars_finalMod)
anova(aruVars_finalMod)
car::vif(aruVars_finalMod) # not used...only 2 vars

## Use only non-eutrophic data
plot(DOC~TDP, statDat_means_noEutro)
plot(log(DOC)~log(TDP), statDat_means_noEutro)
plot(log(secchiDepth)~log(TDP), statDat_means_noEutro)
plot(DOC~secchiDepth, statDat_means_noEutro)


## Goal 1: Predicting historical DOC; apply BsM eqn coefs to ARU data ####
# With ARU data, estimate TP from TDS and DOC
aruDOC_coefs <- coef(aruVars_finalMod)

aruDat_estimated <- aruDat %>% 
  filter(TDS < 300) %>% #Different linear relationship for saline lakes in Chow-Fraser (1991)
  mutate(logDOC = aruDOC_coefs[1] +
         aruDOC_coefs[2]*aruSecchi +
         aruDOC_coefs[3]*DEPMAX,
         estDOC = exp(logDOC)) %>% 
  select(ARU_NM,ARU_LID,latARU,longARU,logTP,TP,lakeTrophicStatus,logDOC,estDOC,aruSecchi,DEPMAX) %>% 
  rename(maxDepth = DEPMAX)

## Goal 1: Predicting historical DOC; compare ARU to BsM ####
compDat <- aruDat_estimated %>% 
  left_join(select(combIDs, ARU_LID,commonID), by = "ARU_LID") %>% 
  inner_join(select(statDat_means_noEutro, commonID,DOC), by="commonID") %>% 
  filter(!is.na(commonID) & 
           lakeTrophicStatus!="eutrophic")  
  
compTtest <- t.test(compDat$estDOC,compDat$DOC)
compTtest #Mean increase of 0.6 mgL between datasets
compLM <- lm(log(DOC)~log(estDOC), compDat)
plot(compLM)
summary(compLM)
plot(log(DOC)~log(estDOC), compDat, xlab="log Historic DOC Estimate (mg/L)", ylab = "log Contemporary DOC (mg/L)")
abline(a = 0, b = 1, col = "red") 
abline(a = 0.5, b = 0.77) #simple linear mod between variables; higher increases at low DOC lakes

## BsM DOC and secchi depth over time
# Subset lakes that have 2 or more sample years
allBsM_mult <- allBsm %>% 
  group_by(waterbodyID) %>% 
  filter(n_distinct(yearSample) >= 2)

statDatWide_mult <- statDat_wide %>% 
  group_by(waterbodyID) %>% 
  filter(n_distinct(yearSample) >= 2)


#Secchi
BsM_secchiMod <- lmer(log(secchiDepth+0.1)~yearSample + (1|waterbodyID), allBsM_mult) #log transformation necessary
plot(BsM_secchiMod)
summary(BsM_secchiMod)

plot(log(secchiDepth+0.1)~yearSample, allBsM_mult)
abline(lm(log(secchiDepth+0.1)~yearSample, allBsM_mult))

#DOC
BsM_DOCMod <- lmer(DOC~yearSample + (1|waterbodyID), statDatWide_mult)
plot(BsM_DOCMod)
summary(BsM_DOCMod)

plot(DOC~yearSample, statDatWide_mult)
abline(lm(DOC~yearSample, statDatWide_mult))



## Side Analysis - what controls water clarity in Ontario lakes? ####

hist(statDat_means$DOC)
hist(statDat_means$TDP)

clarMod <- lm(secchiDepth~log(DOC)*log(TDP), statDat_means)
summary(clarMod)
clarMod_meLog <- lm(secchiDepth~log(DOC)+log(TDP), statDat_means)
summary(clarMod_meLog)
car::vif(clarMod_meLog) #low because only two variables...
clarMod_me <- lm(secchiDepth~(DOC)+(TDP), statDat_means)
summary(clarMod_me)

par(mfrow = c(1, 3))
plot(secchiDepth~log(DOC), statDat_means)
plot(secchiDepth~log(TDP), statDat_means)
plot(log(TDP)~log(DOC), statDat_means)
par(mfrow = c(1, 1))
cor.test(log(statDat_means$DOC), log(statDat_means$TDP))

#What about on raw scale. This is what Olson et al. 2020 use
plot(TDP~DOC, statDat_means)
abline(lm(TDP~DOC, statDat_means))
cor.test(statDat_means$DOC, statDat_means$TDP)
modPtoCrat <- lm(TDP~DOC, statDat_means)
summary(modPtoCrat)

controlDat <- statDat_means %>% 
  mutate(PtoCRatio = TDP/DOC,
         logPtoCRatio = log(TDP)/log(DOC))
plot(secchiDepth~PtoCRatio, controlDat)


