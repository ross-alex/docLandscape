## Is there a relationship between secchi and DOC (and potentially other variables)
# AJR: 2024-02-17

library(tidyverse)
library(here)
library(lubridate)
library(MuMIn)
library(lme4)
library(lmerTest)
library(scales)
library(sf)

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
pDOCresp #To deal with wonky DOC~DOC relationship (and some var. in other variables) need to take mean of multiple observations per waterbodyID. Variation shows change b/w sampling events from the same lake

pDOCresp_log <- ggplot(docRespDat) +
  geom_point(aes(y=log(DOC),x=log(predictorValue)), alpha = 0.3)+
  facet_wrap(~predictorVariable, scales = "free")
pDOCresp_log #To deal with wonky DOC~DOC relationship (and some var. in other variables) need to take mean of multiple observations per waterbodyID. Variation shows change b/w sampling events from the same lake
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
  select(-BsM_Cycle, -dateSample) %>%
  summarise_all(mean, na.rm = TRUE) %>% 
  left_join(select(combIDs, bsmWaterbodyID,commonID), by = c("waterbodyID" = "bsmWaterbodyID")) %>% 
  mutate(yearSample = round(yearSample, digits = 0)) # Round year to nearest whole number; when multiple yearSamples exist

# Look at correlations between different variables
cor.test(statDat_means$secchiDepth, statDat_means$TDP)
cor.test(statDat_means$secchiDepth, statDat_means$DOC)
cor.test(statDat_means$DOC, statDat_means$`COLTR (TCU)`)
cor.test(statDat_means$TDP, statDat_means$`COLTR (TCU)`)

## Pre-step... Are there differences in spring + summer secchi's from the same year? Use Bsm
(secDiffTest <- t.test(secDat_means$meanSpringSecchi,waterClarDiff$meanSummerSecchi)) #Significant. 
# Summer secchi's are 0.38 m deeper than spring

boxplot(secDat_means$meanSpringSecchi,waterClarDiff$meanSummerSecchi, names = c("Spring","Summer"),
        ylab="Secchi Depth (m)", las = 1, ylim=c(16,0))

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
        ylab="Secchi Depth (m)", las = 1, ylim=c(16,0))
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
residuals(clarMod)

plot(secchiDiff~latARU,waterClar_addedVars)
abline(lm(secchiDiff~latARU,waterClar_addedVars))
## OK - this tells us that water clarity is changing in ways we predict. Now let's model start the DOC modelling process

## Water clarity only increasing in eutrophic lakes, choose oligo and mesotrophic lakes
## Make dataframe using mean values across different BsM cycles
statDat_means_noEutro <- statDat_means %>% 
  filter(TDP < 10) # Wetzel says 20ugL TP and less is meso/oligo; also stats that TDP is ~50% of TP

## Chose either all lakes (statDat_means) or all non-eutrophic lakes (statDat_means_noEutro)
statDat_dataSet <- statDat_means
statDat_dataSet <- statDat_means_noEutro

## Step 2: Create a model from BsM data
## Look at distribution of DOC observations; what distribution is most appropriate
# Linear
plot_DOCdensity <- ggplot(data = statDat_dataSet) +
  geom_density(aes(x=DOC), fill="grey", alpha = 0.7) +
  theme_bw()
plot_DOCdensity #right skew

# log-linear
plot_DOCdensity_log <- plot_DOCdensity +
  scale_x_continuous(trans = "log")
plot_DOCdensity_log #left skew

# gamma distribution
# Fit gamma distribution to the observed data
fit <- MASS::fitdistr(statDat_dataSet$DOC, densfun = "gamma")

# Visualize the fit
ggplot(statDat_dataSet, aes(x = DOC)) +
  geom_histogram(aes(y = ..density..), bins = 30, fill = "skyblue", color = "black", alpha = 0.7) +
  stat_function(fun = dgamma, args = list(shape = fit$estimate["shape"], rate = fit$estimate["rate"]), color = "red", size = 1) +
  labs(x = "DOC", y = "Density") +
  ggtitle("Gamma Distribution Fit to Data") 
# gamma fits well with shape = 4.91, rate = 0.79

##### FIT MODELS!
## Model ARU data using GLM, gamma distribution - This chunk has outliers, use next
# # Fit the initial model
# ARUmod_gamma <- glm(DOC ~ secchiDepth * maxDepth, data = statDat_dataSet, family = Gamma(link = "log"), na.action = "na.fail")
# ARUmodSel <- dredge(ARUmod_gamma) #interaction is best model
# 
# ## Best models
# # All lakes
# ARUmod_gammaBest <- glm(DOC ~ secchiDepth * maxDepth, data = statDat_dataSet, family = Gamma(link = "log"))
# # All non-eutrophic lakes
# ARUmod_gammaBest <- glm(DOC ~ secchiDepth + maxDepth, data = statDat_dataSet, family = Gamma(link = "log"))
# 
# #  Model diagnostics
# DHARMa::simulateResiduals(ARUmod_gammaBest, plot=T)
# plot(ARUmod_gammaBest) # Diagnositcs not great but probably fine
# car::vif(ARUmod_gammaBest) # VIF of interaction is 8, need to use just additive model
# 
# # model summaries
# summary(ARUmod_gammaBest)

## Gamma dist. has three potential outliers; rows 895, 867, 437. Reduce dataset try again
gammaRedDat <- statDat_dataSet[-c(895,867,437),]

# Fit the initial model
ARUmod_gamma_red <- glm(DOC ~ secchiDepth * maxDepth, data = gammaRedDat, family = Gamma(link = "log"), na.action = "na.fail")
ARUmodSel_red <- dredge(ARUmod_gamma_red) #interaction is best model

## Best models
# All lakes
ARUmod_gammaBest_red <- glm(DOC ~ secchiDepth * maxDepth, data = gammaRedDat, family = Gamma(link = "log"))
# All non-eutrophic lakes
ARUmod_gammaBest_red <- glm(DOC ~ secchiDepth + maxDepth, data = gammaRedDat, family = Gamma(link = "log"))

#  Model diagnostics
DHARMa::simulateResiduals(ARUmod_gammaBest_red, plot=T) #Better
plot(ARUmod_gammaBest_red) # Diagnositcs not great but probably fine
car::vif(ARUmod_gammaBest_red) # VIF of interaction is 8, need to use just additive model

# model summaries
glmSum <- summary(ARUmod_gammaBest_red)

# #### What about  straight-up linear models - outliers in this data chunk, use reduced dataset that follows
# aruVarsMod <- lm(log(DOC)~secchiDepth*maxDepth, data=statDat_dataSet, na.action = "na.fail")
# aruVarsSel <- dredge(aruVarsMod) #interaction is not important important
# 
# ## Best models
# # All lakes
# #aruVars_finalMod <- lm(log(DOC)~secchiDepth*maxDepth, data=statDat_dataSet)
# # All non-eutrophic lakes
# aruVars_finalMod <- lm(log(DOC)~secchiDepth+maxDepth, data=statDat_dataSet) 
# 
# plot(aruVars_finalMod) #diagnostics slightly better
# summary(aruVars_finalMod)
# anova(aruVars_finalMod)
# car::vif(aruVars_finalMod) # VIF of interaction is 8, need to use just additive model

## Some outliers in data set, redo without points 895, 867, 543
statDat_reduce <- statDat_dataSet[-c(895,867,543),]

# Model selection
aruVarsMod_red <- lm(log(DOC)~secchiDepth*maxDepth, data=statDat_reduce, na.action = "na.fail")
aruVarsSel_red <- dredge(aruVarsMod_red) #interaction is not important important

# All lakes
aruVars_finalMod_red <- lm(log(DOC)~secchiDepth*maxDepth, data=statDat_reduce)
# All non-eutrophic lakes
aruVars_finalMod_red <- lm(log(DOC)~secchiDepth+maxDepth, data=statDat_reduce) 

plot(aruVars_finalMod_red) #diagnostics slightly better
summary(aruVars_finalMod_red)
anova(aruVars_finalMod_red)
car::vif(aruVars_finalMod_red) # VIF of interaction is 8, need to use just additive model

## What does this look like?
plot(DOC~TDP, statDat_reduce)
plot(log(DOC)~log(TDP), statDat_reduce)
plot(log(secchiDepth)~log(TDP), statDat_reduce)
plot(log(DOC)~log(secchiDepth), statDat_reduce)


## Goal 1: Predicting historical DOC; apply BsM eqn coefs to ARU data ####
# With ARU data, estimate TP from TDS and DOC
aruDOC_coefs <- coef(ARUmod_gammaBest_red) #GLM
aruDOC_coefs <- coef(aruVars_finalMod) #LM

# Coef DF using GLM
aruDat_estimated <- aruDat %>% 
  rename(secchiDepth = aruSecchi, maxDepth = DEPMAX) %>% 
  filter(TDS < 300) %>% #Different linear relationship for saline lakes in Chow-Fraser (1991)
  mutate(logDOC = predict(ARUmod_gammaBest_red, newdata = .),
         estDOC = exp(logDOC)) %>% 
  select(ARU_NM,ARU_LID,INV_YR,latARU,longARU,logTP,TP,lakeTrophicStatus,estDOC,secchiDepth,maxDepth) %>% 
  rename(aruSecchi = secchiDepth,yearSampleARU = INV_YR)

# Coef DF using LM
aruDat_estimated <- aruDat %>% 
  filter(TDS < 300) %>% #Different linear relationship for saline lakes in Chow-Fraser (1991)
  mutate(logDOC = aruDOC_coefs[1] +
         aruDOC_coefs[2]*aruSecchi +
         aruDOC_coefs[3]*DEPMAX,
         estDOC = exp(logDOC)) %>% 
  select(ARU_NM,ARU_LID,INV_YR,latARU,longARU,logTP,TP,lakeTrophicStatus,logDOC,estDOC,aruSecchi,DEPMAX) %>% 
  rename(maxDepth = DEPMAX,yearSampleARU = INV_YR)

## Goal 1: Predicting historical DOC; compare ARU to BsM ####
# No eutrophic lakes
compDat <- aruDat_estimated %>% 
  left_join(select(combIDs, ARU_LID,commonID), by = "ARU_LID") %>% 
  inner_join(select(statDat_dataSet, yearSample, commonID,DOC), by="commonID") %>% 
  filter(!is.na(commonID) & 
          # lakeTrophicStatus!="eutrophic" &
           !is.na(aruSecchi) &
           !is.na(maxDepth))  %>% 
  mutate(docDiff = DOC-estDOC,
         yearDiff = yearSample - yearSampleARU)
  
compTtest <- t.test(compDat$estDOC,compDat$DOC)
compTtest #Mean increase of 0.58 mgL between datasets when using LM
# Mean increase of 0.43 mgL between datasets when using GLM

# Test rate of change
compLM <- lm(log(DOC)~log(estDOC), compDat)
plot(compLM)
summary(compLM)
residComparison <- residuals(compLM)
compDat_wResid <- cbind(compDat,residComparison)

pDOCcomp <- ggplot(compDat, aes(x=estDOC,y=DOC)) +
  geom_point(aes(colour=latARU)) +
  scale_y_continuous(trans = "log", limits = c(1,18), labels = label_number(digits = 2)) +
  scale_x_continuous(trans = "log", limits = c(1,18), labels = label_number(digits = 2)) +
  geom_abline(linetype="dotted") +
  geom_smooth(method = "lm", colour="black") +
  xlab("Historic DOC estimate (mg/L)") +
  ylab("Contemporary DOC measurement (mg/L)") +
  theme_bw()
pDOCcomp

# Plot historic~contemporary by lake trophic status
pDOCcomp_trophic <- pDOCcomp +
  facet_wrap(~lakeTrophicStatus, nrow = 1)

# Higher increase in historically clearer lakes
ontario_shape <- st_read(here("data","OBM_INDEX.shp")) 

residPLoc <- ggplot(compDat_wResid, aes(x=longARU,y=latARU,colour=residComparison)) +
  geom_point()

residMap <-ggplot() +
  geom_sf(data = ontario_shape, colour="grey", alpha = 0.3) +
  geom_point(data = compDat_wResid, aes(x = longARU, y = latARU, colour = residComparison)) +
  scale_color_gradient2(mid = "white", low = "red", high = "blue", name = "Residuals") +  # Specify the gradient from red to white to blue
  labs(title = "Residuals: red lakes are lighter, blue are darker") +
  theme_minimal()
residMap

## Test/look at differences over time
# Change reference category to oligotrophic
compDat$lakeTrophicStatus <- relevel(factor(compDat$lakeTrophicStatus), ref = "oligotrophic")
compDat$lakeTrophicStatus <- factor(compDat$lakeTrophicStatus, levels = c("oligotrophic","mesotrophic","eutrophic"))

docOverTime <- lm(docDiff~yearDiff+latARU+lakeTrophicStatus+maxDepth, compDat)
plot(docOverTime)
summary(docOverTime)
anova(docOverTime)
TukeyHSD(aov(docDiff~yearDiff+latARU+lakeTrophicStatus+maxDepth, compDat), which = "lakeTrophicStatus")

docDiffMap <-ggplot() +
  geom_sf(data = ontario_shape, colour="grey", alpha = 0.3) +
  geom_point(data = compDat, aes(x = longARU, y = latARU, colour = docDiff,  size = yearDiff)) +
  scale_color_gradient2(mid = "white", low = "blue", high = "brown4", name = "DOC change (mg/L)") +  # Specify the gradient from red to white to blue
  labs(title = "Blue lakes are lighter, brown are darker", size = "Years b/w samples") +
  theme_minimal() +
  facet_wrap(~lakeTrophicStatus, nrow = 1)
docDiffMap

## Summarize how many went up and down, and average for those
docChangeSum <- compDat %>% 
  group_by(lakeTrophicStatus) %>% 
  summarise(n_DOC_increase = sum(docDiff > 0),
            n_DOC_decrease = sum(docDiff < 0),
            percentIncrease = (n_DOC_increase/(n_DOC_increase+n_DOC_decrease))*100,
            avgIncrease_mgL = mean(docDiff > 0),
            avgDecrease_mgL = mean(docDiff < 0),
            medianChange_mgL = median(docDiff),
            maxChange_mgL = max(docDiff), 
            meanYearARU = mean(yearSampleARU),
            meanYearBsM = mean(yearSample))

## BsM DOC and secchi depth over time
# Subset lakes that have 2 or more sample years
allBsM_mult <- allBsm %>% 
  group_by(waterbodyID) %>% 
  filter(n_distinct(yearSample) >= 2) %>% 
  filter(predictorVariable != "TDP" | (predictorVariable == "TDP" & predictorValue < 10.01))


allMult_check <- allBsM_mult %>% 
  select(lakeName,yearSample) %>% 
  distinct()

statDatWide_mult <- statDat_wide %>% 
  group_by(waterbodyID) %>% 
  filter(n_distinct(yearSample) >= 2)  %>% 
  filter(TDP < 10.01) #Remember, this df has slighlty fewer waterbodyIDs because model data with NAs (any) had to be removed


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

# Subset lakes that have 3 or more sample years
allBsM_mult <- allBsm %>% 
  group_by(waterbodyID) %>% 
  filter(n_distinct(yearSample) >= 3)

allMult_check <- allBsM_mult %>% 
  select(lakeName,yearSample) %>% 
  distinct()

statDatWide_mult <- statDat_wide %>% 
  group_by(waterbodyID) %>% 
  filter(n_distinct(yearSample) >= 3)

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

predDOC <- predict(BsM_DOCMod)
BsMDOC_predMultiples <- cbind(statDatWide_mult,predDOC) %>% 
  arrange(waterbodyID,yearSample)
