## Is there a relationship between secchi and DOC (and potentially other variables)
# AJR: 2024-02-17

library(tidyverse)
library(here)
library(lubridate)
library(MuMIn)

###***ALEX, should I Z-score DOC? Doesn't make a ton of sense when comparing the contemporary/historic...I think

## Read Data ####
allBsm <- read_csv(here("data","BsM_modelData_20240311.csv"))
allARU <- read_csv(here("data","ARU_modelData_20240311.csv"))
statDat_wide <- read_csv(here("data","BsM_WideSampleModelData_20240311.csv")) %>% #This dataset is slightly reduced
                  filter(log(TKN) >4) #outliers <4                            #to remove all rows that had at 
                                                                              #at least one NA. See 'modelData_assembly_20240311.R'
                                                                              # script for details (bottom chunk)
combIDs <- read_csv(here("data","commonDataIDs_aruToBsm_20240307.csv")) #has common IDs between BsM and ARU


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

## Goal 1: Predicting historical DOC; model selection ####
## The ARU database has secchi and depth, which can be used to predict DOC. 
# A regression can be applied to BsM data where DOC is known, and the parameters from this model can
# then be applied to the ARU database. DOC distribution seems somewhat gamma distributed. Check models
# for linear, log-linear, and gamma dist. for response variable (DOC)

# Additionally, check if removing eutrophic lakes is necessary. (it doesn't seem it is, so long as TDP is logged)

## Use all data (including eutrophic lakes)

## Make dataframe using mean values across different BsM cycles
statDat_means <- statDat_wide %>%
  group_by(waterbodyID) %>% 
  select(-BsM_Cycle, -dateSample, -yearSample) %>%
  summarise_all(mean, na.rm = TRUE) %>% 
  left_join(select(combIDs, bsmWaterbodyID,commonID), by = c("waterbodyID" = "bsmWaterbodyID"))

## Look at distribution of DOC observations; what distribution is most appropriate
# Linear
plot_DOCdensity <- ggplot(data = statDat_means) +
  geom_density(aes(x=DOC), fill="grey", alpha = 0.7) +
  theme_bw()
plot_DOCdensity #right skew

# log-linear
plot_DOCdensity_log <- plot_DOCdensity +
  scale_x_continuous(trans = "log")
plot_DOCdensity_log #left skew

# gamma distribution
# Fit gamma distribution to the observed data
fit <- MASS::fitdistr(statDat_means$DOC, densfun = "gamma")

# Visualize the fit
ggplot(statDat_means, aes(x = DOC)) +
  geom_histogram(aes(y = ..density..), bins = 30, fill = "skyblue", color = "black", alpha = 0.7) +
  stat_function(fun = dgamma, args = list(shape = fit$estimate["shape"], rate = fit$estimate["rate"]), color = "red", size = 1) +
  labs(x = "DOC", y = "Density") +
  ggtitle("Gamma Distribution Fit to Data") 
# gamma fits well with shape = 4.52, rate = 0.60

## Model ARU data using GLM, gamma distribution
# Step 1: Fit the initial model
ARUmod_gamma <- glm(DOC ~ secchiDepth * maxDepth, data = statDat_means, family = Gamma(link = "log"), na.action = "na.fail")
ARUmodSel <- dredge(ARUmod_gamma) #interaction is best model

ARUmod_gammaBest <- ARUmod_gamma

# Step 2: Model diagnostics
plot(ARUmod_gammaBest)

# Step 3: Model summaries
summary(ARUmod_gammaBest)

## What about with straight-up linear models
aruVarsMod <- lm(log(DOC)~secchiDepth*maxDepth, data=statDat_means, na.action = "na.fail")
aruVarsSel <- dredge(aruVarsMod) #interaction is important

aruVars_finalMod <- lm(log(DOC)~secchiDepth*maxDepth, data=statDat_means)
plot(aruVars_finalMod)
summary(aruVars_finalMod)
car::vif(aruVars_finalMod) #interaction still kind of high (8.14)

aruVars_mainEffects <- lm(log(DOC)~secchiDepth+maxDepth, data=statDat_means)
plot(aruVars_mainEffects)
summary(aruVars_mainEffects)
car::vif(aruVars_mainEffects)


## Use all only non-eutrophic data
# First off, does it make sense to not include eutrophic lakes?
plot(DOC~TDP, statDat_means)
plot(log(DOC)~log(TDP), statDat_means)
plot(log(secchiDepth)~log(TDP), statDat_means)
plot(DOC~secchiDepth, statDat_means)

## Make dataframe using mean values across different BsM cycles
statDat_means_noEutro <- statDat_wide %>%
  group_by(waterbodyID) %>% 
  select(-BsM_Cycle, -dateSample, -yearSample) %>%
  summarise_all(mean, na.rm = TRUE) %>% 
  filter(TDP < 20)

## Look at distribution of DOC observations; what distribution is most appropriate
# Linear
plot_DOCdensity_noEutro <- ggplot(data = statDat_means_noEutro) +
  geom_density(aes(x=DOC), fill="grey", alpha = 0.7) +
  theme_bw()
plot_DOCdensity_noEutro #right skew

# log-linear
plot_DOCdensity_log_noEutro <- plot_DOCdensity_noEutro +
  scale_x_continuous(trans = "log")
plot_DOCdensity_log_noEutro #left skew

# gamma distribution
# Fit gamma distribution to the observed data
fit <- MASS::fitdistr(statDat_means_noEutro$DOC, densfun = "gamma")

# Visualize the fit
ggplot(statDat_means, aes(x = DOC)) +
  geom_histogram(aes(y = ..density..), bins = 30, fill = "skyblue", color = "black", alpha = 0.7) +
  stat_function(fun = dgamma, args = list(shape = fit$estimate["shape"], rate = fit$estimate["rate"]), color = "red", size = 1) +
  labs(x = "DOC", y = "Density") +
  ggtitle("Gamma Distribution Fit to Data") 
# gamma fits well with shape = 4.52, rate = 0.60

## Model ARU data using GLM, gamma distribution
# Step 1: Fit the initial model
ARUmod_gamma_noEutro <- glm(DOC ~ secchiDepth * maxDepth, data = statDat_means_noEutro, family = Gamma(link = "log"), na.action = "na.fail")
ARUmodSel_noEutro <- dredge(ARUmod_gamma_noEutro) #interaction is best model

ARUmod_gammaBest_noEutro <- ARUmod_gamma_noEutro

# Step 2: Model diagnostics
plot(ARUmod_gammaBest_noEutro)

# Step 3: Model summaries
summary(ARUmod_gammaBest_noEutro)

## What about with straight-up linear models
aruVarsMod_noEutro <- lm(log(DOC)~secchiDepth*maxDepth, data=statDat_means_noEutro, na.action = "na.fail")
aruVarsSel_noEutro <- dredge(aruVarsMod_noEutro) #interaction is important

aruVars_finalMod_noEutro <- lm(log(DOC)~secchiDepth*maxDepth, data=statDat_means_noEutro)
plot(aruVars_finalMod_noEutro)
summary(aruVars_finalMod_noEutro)
car::vif(aruVars_finalMod_noEutro) #interaction still kind of high (10.0)

aruVars_mainEffects_noEutro <- lm(log(DOC)~secchiDepth+maxDepth, data=statDat_means_noEutro)
plot(aruVars_mainEffects_noEutro)
summary(aruVars_mainEffects_noEutro)
car::vif(aruVars_mainEffects_noEutro)

##**Take home, do not need to subset eutrophic data. It makes R2 marginally worse, but in any case, doesn't hinder model
## So long everything is on logDOC~logTP scale

## Goal 1: Predicting historical DOC; apply BsM eqn coefs to ARU data ####
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
         estDOC = exp(logDOC)) %>% 
  select(ARU_NM,ARU_LID,latARU,longARU,logTP,TP,lakeTrophicStatus,logDOC,estDOC,SECCHI) 

## Goal 1: Predicting historical DOC; compare ARU to BsM ####
compDat <- aruDat_estimated %>% 
  left_join(select(combIDs, ARU_LID,commonID), by = "ARU_LID") %>% 
  inner_join(select(statDat_means, commonID,DOC), by="commonID") %>% 
  filter(!is.na(commonID))

compTtest <- lm(DOC~estDOC, compDat)
plot(compTtest)
summary(compTtest)
plot(DOC~estDOC, compDat, xlab="Historic DOC Estimate (mg/L)", ylab = "Contemporary DOC (mg/L)", ylim=c(0,21),xlim=c(0,21))
abline(a = 0, b = 1, col = "red") 
abline(a = 1.07, b = 0.96) #simple linear mod between variables; higher increases at low DOC lakes

## BsM DOC and secchi depth over time

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


