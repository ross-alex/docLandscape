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
#allBsm <- read_csv(here("data","BsM_modelData_20240311.csv"))
allBsm <- read_csv(here("data","allBsM_updatedCyc3_20240520.csv"))
allARU <- read_csv(here("data","ARU_modelData_20240311.csv"))
#statDat_wide <- read_csv(here("data","BsM_WideSampleModelData_20240311.csv")) %>% 
statDat_wide <- read_csv(here("data","statDatWide_updatedCyc3_20240520.csv")) %>% #This dataset is slightly reduced
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
  summarise(meanSummerSecchi = mean(summerSecchi,na.rm=TRUE),
            meanSpringSecchi = mean(springSecchi,na.rm=TRUE),
            waterClarDiff = meanSummerSecchi-meanSpringSecchi,
            sdSummerSecchi = sd(summerSecchi,na.rm=TRUE),
            sdSpringSecchi = sd(springSecchi,na.rm=TRUE))

## Goal 1: Checking for differences in water clarity; Predicting historical DOC; model selection ####
## The ARU database has secchi and depth, which can be used to predict DOC. 
# A regression can be applied to BsM data where DOC is known, and the parameters from this model can
# then be applied to the ARU database. DOC distribution seems somewhat gamma distributed. Check models
# for linear, log-linear, and gamma dist. for response variable (DOC)
# Additionally, check if removing eutrophic lakes is necessary. 

## Make dataframe using mean values across different BsM cycles; using means as a conservative "contempoary" value.
# Could also take the max (most recent year), though this affects interpretation (some year gaps will be very different; short or long)
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
(secDiffTest <- t.test(secDat_means$meanSpringSecchi,secDat_means$meanSummerSecchi,paired = TRUE)) #Significant. 
# Summer secchi's are 0.31 m deeper than spring. No use accounting for difference - diff. is smaller than mean sd

boxplot(secDat_means$meanSpringSecchi,secDat_means$meanSummerSecchi, names = c("Spring","Summer"),
        ylab="Secchi Depth (m)", las = 1, ylim=c(16,0))

## Step 1: Are there differences in water clarity between the two times
waterClarDiff <- statDat_means %>% 
  select(waterbodyID,commonID,secchiDepth) %>%                          # Using spring BsM secchi; more consistent (and above analysis shows doesn't matter)
  left_join(select(aruDat, commonID, aruSecchi, lakeTrophicStatus)) %>% 
  filter_all(all_vars(!is.na(.))) %>% 
  mutate(secchiDiff = secchiDepth-aruSecchi) %>% 
  mutate(lakeTrophicStatus = factor(lakeTrophicStatus, levels = c("oligotrophic","mesotrophic","eutrophic")))

waterClar_addedVars <- waterClarDiff %>% 
  left_join(select(aruDat, commonID, latARU,longARU, DEPMAX)) %>% 
  mutate(centeredSecchiDiff = scale(secchiDiff,scale=FALSE),
         centeredLat = scale(latARU,scale=FALSE),
         centeredLong = scale(longARU,scale=FALSE),
         centeredMaxDepth = scale(DEPMAX,scale=FALSE))
waterClar_addedVars_clean <- na.omit(waterClar_addedVars)

# T test
(clarTtest <- t.test(waterClarDiff$aruSecchi,waterClarDiff$secchiDepth)) #Significant. ARU has deeper secchis = clearer lakes
# Secchi's are 0.55 m shallower now

# Visualize it
boxplot(waterClarDiff$aruSecchi,waterClarDiff$secchiDepth, names = c("Historic","Contemporary"),
        ylab="Secchi Depth (m)", las = 1, ylim=c(16,0))
pClar <- ggplot(waterClarDiff) +
  geom_density(aes(x=secchiDiff),fill="grey",alpha=0.5) +
  geom_vline(xintercept = 0, linetype = "dotted", color = "red") +
  ggtitle("Negative values indicate a lake is darker now than historically") +
  theme_bw()
pClar #Mean and majority data <0. These are using all data. Let's look at trophic status

# What happens in each trophic group? Oligotrophic most negative...
pClar_trophicStatus <- pClar +
  facet_wrap(~lakeTrophicStatus)
pClar_trophicStatus  
# Oligotrophic lakes getting darker, meso also; eutrophic looks like a normal dist. to me

# Let's model it, see if there are variables driving differences
# *The chunk below likely not going to be used. Not important: secchi differences
# are to know that clarity changes, don't care the variables...but here anyways
clarModDredge <- lm(secchiDiff ~ lakeTrophicStatus*centeredLat*centeredLong*centeredMaxDepth, data = waterClar_addedVars_clean, na.action="na.fail")
(clarDredgeOut <- dredge(clarModDredge))
clarMod <- lm(secchiDiff~lakeTrophicStatus+centeredLat+centeredLong+centeredMaxDepth+
                centeredLat:centeredLong+centeredLong:centeredMaxDepth+centeredLong:lakeTrophicStatus
              , data=waterClar_addedVars_clean)
plot(clarMod) #diagnostics aren't great
summary(clarMod)
anova(clarMod) 

## Step 2: Create a model from BsM data
## BsM dataSet
statDat_dataSet <- statDat_means %>% 
  filter(!is.na(maxDepth)) 

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
# Fit the initial model
ARUmod_gamma <- glm(DOC ~ scale(secchiDepth,scale=FALSE) * scale(maxDepth,scale=FALSE), data = statDat_dataSet, family = Gamma(link = "log"), na.action = "na.fail")
(ARUmodSel <- dredge(ARUmod_gamma)) #interaction has low deltaAIC - use it

## Best models
# All lakes
ARUmod_gammaBest <- glm(DOC ~ scale(secchiDepth,scale=FALSE) * scale(maxDepth,scale=FALSE), data = statDat_dataSet, family = Gamma(link = "log"))

#  Model diagnostics
DHARMa::simulateResiduals(ARUmod_gammaBest, plot=T)
plot(ARUmod_gammaBest) # Diagnositcs not great but probably fine
car::vif(ARUmod_gammaBest) # VIF of interaction is 8, need to use just additive model

# model summaries
summary(ARUmod_gammaBest)

## Gamma dist. has three potential outliers; rows 901. Reduce dataset try again
gammaRedDat <- statDat_dataSet[-901,]

# Fit the initial model
ARUmod_gamma_red <- glm(DOC ~ scale(secchiDepth,scale=FALSE) * scale(maxDepth,scale=FALSE), data = gammaRedDat, family = Gamma(link = "log"), na.action = "na.fail")
(ARUmodSel_red <- dredge(ARUmod_gamma_red)) #interaction is best

## Best models
# All lakes
ARUmod_gammaBest_red <- glm(DOC ~ scale(secchiDepth,scale=FALSE) * scale(maxDepth,scale=FALSE), data = gammaRedDat, family = Gamma(link = "log"))

#  Model diagnostics
DHARMa::simulateResiduals(ARUmod_gammaBest_red, plot=T) #Better
plot(ARUmod_gammaBest_red) # Diagnositcs not great but fine
car::vif(ARUmod_gammaBest_red) # VIFs all good (data needs to be centered, otherwise VIFs too high)

# model summaries
(glmSum <- summary(ARUmod_gammaBest_red))
# Deviance-based R2 - seems legit..
deviance_model <- deviance(glmSum)
deviance_null <- deviance(update(glmSum, . ~ 1))
(R2_deviance <- 1 - (deviance_model / deviance_null)) # R2 0.61

# #### What about  straight-up linear models - outliers in this data chunk, use reduced dataset that follows
aruVarsMod <- lm(log(DOC)~scale(secchiDepth,scale=FALSE) * scale(maxDepth,scale=FALSE), data=statDat_dataSet, na.action = "na.fail")
(aruVarsSel <- dredge(aruVarsMod)) #interaction AICc nearly identical to non-interaction; use interaction 

## Best models
# All lakes
aruVars_finalMod <- lm(log(DOC)~scale(secchiDepth,scale=FALSE) * scale(maxDepth,scale=FALSE), data=statDat_dataSet)

plot(aruVars_finalMod) #diagnostics slightly better not great; outlier
summary(aruVars_finalMod)
anova(aruVars_finalMod)
car::vif(aruVars_finalMod) # VIF of interaction is 7.6, would need to use just additive model but diagnostics whack anyways

## Some outliers in data set, redo without points 901
statDat_reduce <- statDat_dataSet[-c(901),]

# Model selection
aruVarsMod_red <- lm(log(DOC)~scale(secchiDepth,scale=FALSE) * scale(maxDepth,scale=FALSE), data=statDat_reduce, na.action = "na.fail")
(aruVarsSel_red <- dredge(aruVarsMod_red)) #interaction is very important

# All lakes
aruVars_finalMod_red <- lm(log(DOC)~scale(secchiDepth,scale=FALSE) * scale(maxDepth,scale=FALSE), data=statDat_reduce)

plot(aruVars_finalMod_red) #diagnostics slightly better
summary(aruVars_finalMod_red) #R2 = 0.63; very similar to GLM of 0.61
anova(aruVars_finalMod_red)
car::vif(aruVars_finalMod_red) # VIFs all good, with centered data

## Goal 1: Predicting historical DOC; apply BsM eqn coefs to ARU data ####
# With ARU data, estimate TP from TDS and DOC
(aruDOC_coefs_GLM <- coef(ARUmod_gammaBest_redFinal)) #GLM coefs - will be using these, but LM included anyways
(aruDOC_coefs_LM <- coef(aruVars_finalMod_redNoInt)) #LM coefs

# Coef DF using GLM
aruDat_estimated_GLM <- aruDat %>% 
  rename(secchiDepth = aruSecchi, maxDepth = DEPMAX) %>% 
  filter(TDS < 300) %>% #Different linear relationship for saline lakes in Chow-Fraser (1991)
  mutate(logDOC = predict(ARUmod_gammaBest_red, newdata = .),
         estDOC = exp(logDOC)) %>% 
  select(ARU_NM,ARU_LID,INV_YR,latARU,longARU,logTP,TP,lakeTrophicStatus,estDOC,secchiDepth,maxDepth) %>% 
  rename(aruSecchi = secchiDepth,yearSampleARU = INV_YR)

# Coef DF using LM
aruDat_estimated_LM <- aruDat %>% 
  filter(TDS < 300) %>% #Different linear relationship for saline lakes in Chow-Fraser (1991)
  mutate(logDOC = aruDOC_coefs_LM[1] +
         aruDOC_coefs_LM[2]*aruSecchi +
         aruDOC_coefs_LM[3]*DEPMAX,
         estDOC = exp(logDOC)) %>% 
  select(ARU_NM,ARU_LID,INV_YR,latARU,longARU,logTP,TP,lakeTrophicStatus,logDOC,estDOC,aruSecchi,DEPMAX) %>% 
  rename(maxDepth = DEPMAX,yearSampleARU = INV_YR)

## Goal 1: Predicting historical DOC; compare ARU (historic) to BsM (contemp.) ####
compDat_GLM <- aruDat_estimated_GLM %>% 
  distinct() %>%                                       # necessary to eliminate multiple entries for a given ARU_LID (e.g., multiple Shebandowan lakes (seperate units))
  left_join(select(combIDs, ARU_LID,commonID), by = "ARU_LID") %>% 
  inner_join(select(statDat_dataSet, yearSample, commonID,DOC), by="commonID") %>% 
  filter(!is.na(commonID) & 
           !is.na(aruSecchi) &
           !is.na(maxDepth))  %>% 
  mutate(docDiff = DOC-estDOC,
         yearDiff = yearSample - yearSampleARU)

compDat_LM <- aruDat_estimated_LM %>% 
  distinct() %>%                                       # necessary to eliminate multiple entries for a given ARU_LID (e.g., multiple Shebandowan lakes (seperate units))
  left_join(select(combIDs, ARU_LID,commonID), by = "ARU_LID") %>% 
  inner_join(select(statDat_dataSet, yearSample, commonID,DOC), by="commonID") %>% 
  filter(!is.na(commonID) & 
           !is.na(aruSecchi) &
           !is.na(maxDepth))  %>% 
  mutate(docDiff = DOC-estDOC,
         yearDiff = yearSample - yearSampleARU)

# Simple statistical test: 2-sample t-test  
(compTtest_GLM <- t.test(compDat_GLM$estDOC,compDat_GLM$DOC)) #Mean increase of 0.56 mgL between datasets when using GLM
(compTtest_LM <- t.test(compDat_LM$estDOC,compDat_LM$DOC)) #Mean increase of 0.85 mgL between datasets when using LM

# Test rate of change - GLM
compGLM <- lm(log(DOC)~log(estDOC), compDat_GLM)
plot(compGLM) # Good enough
summary(compGLM) # significant increase, R2 0.54; summary not all that useful apart from saying difference b/w time periods
residComparison <- residuals(compGLM)
compDat_wResid <- cbind(compDat_GLM,residComparison)

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
(pDOCcomp_trophic <- pDOCcomp +
  facet_wrap(~lakeTrophicStatus, nrow = 1))

# Higher increase in historically clearer lakes
ontario_shape <- st_read(here("data","OBM_INDEX.shp")) 

# Look at residuals
residPLoc <- ggplot(compDat_wResid, aes(x=longARU,y=latARU,colour=residComparison)) +
  geom_point()

residMap <-ggplot() +
  geom_sf(data = ontario_shape, colour="grey", alpha = 0.3) +
  geom_point(data = compDat_wResid, aes(x = longARU, y = latARU, colour = residComparison)) +
  scale_color_gradient2(mid = "white", low = "red", high = "blue", name = "Residuals") +  # Specify the gradient from red to white to blue
  labs(title = "Residuals: red lakes are lighter, blue are darker") +
  theme_minimal()
residMap

## Now let's model what/where [DOC] is changing 
# Variable justification. We've seen differences in lake trophic status. Seems to be a S/N latitudinal gradient independent of
# other variables, and possible additive effect of having longer gaps between sample events
DOCtempFull <- lm(docDiff~lakeTrophicStatus+scale(latARU,scale=FALSE)+
                    scale(yearDiff,scale=FALSE)+
                    scale(maxDepth,scale=FALSE), compDat_GLM)
plot(DOCtempFull)
car::vif(DOCtempFull)
summary(DOCtempFull)
anova(DOCtempFull) #Differences driven mostly by latitude, max depth, and lake trophic status also significant
TukeyHSD(aov(docDiff~lakeTrophicStatus+scale(latARU,scale=FALSE)+
               scale(yearDiff,scale=FALSE)+
               scale(maxDepth,scale=FALSE), compDat_GLM), which = "lakeTrophicStatus") #Sign. model difference; no groups diff.


## Test/look at differences over time
# Change reference category to oligotrophic
compDat_GLM$lakeTrophicStatus <- relevel(factor(compDat$lakeTrophicStatus), ref = "oligotrophic")
compDat_GLM$lakeTrophicStatus <- factor(compDat$lakeTrophicStatus, levels = c("oligotrophic","mesotrophic","eutrophic"))

docDiffMap <-ggplot() +
  geom_sf(data = ontario_shape, colour="grey", alpha = 0.3) +
  geom_point(data = compDat_GLM, aes(x = longARU, y = latARU, colour = docDiff,  size = yearDiff)) +
  scale_color_gradient2(mid = "white", low = "blue", high = "brown4", name = "DOC change (mg/L)") +  # Specify the gradient from red to white to blue
  labs(title = "Blue lakes are lighter, brown are darker", size = "Years b/w samples") +
  theme_minimal() +
  facet_wrap(~lakeTrophicStatus, nrow = 1)
docDiffMap

## Summarize how many went up and down, and average for those
(docChangeSum <- compDat %>% 
  group_by(lakeTrophicStatus) %>% 
  summarise(n_DOC_increase = sum(docDiff > 0),
            n_DOC_decrease = sum(docDiff < 0),
            percentIncrease = (n_DOC_increase/(n_DOC_increase+n_DOC_decrease))*100,
            avgIncrease_mgL = mean(docDiff > 0),
            avgDecrease_mgL = mean(docDiff < 0),
            medianChange_mgL = median(docDiff),
            maxChange_mgL = max(docDiff), 
            meanYearARU = mean(yearSampleARU),
            meanYearBsM = mean(yearSample)))

## BsM DOC and secchi depth over time
# Subset lakes that have 2 or more sample years
allBsM_mult <- allBsm %>% 
  group_by(waterbodyID) %>% 
  filter(n_distinct(yearSample) >= 2) %>% 
  filter(predictorVariable != "TDP" | (predictorVariable == "TDP" & predictorValue < 10.01))

allMult_check <- allBsM_mult %>% 
  select(lakeName,yearSample) %>% 
  distinct()

statDat_wide_multTrophicCat <- statDat_wide %>% 
  group_by(waterbodyID) %>% 
  filter(n_distinct(yearSample) >= 2)  %>% 
  mutate(lakeTrophicStatus =  case_when(TDP < 5 ~ "oligotrophic",
                                        TDP > 5 & TDP < 10 ~ "mesotrophic",
                                        TDP > 10 ~ "eutrophic")) %>% 
  filter(!is.na(lakeTrophicStatus) & !is.na(yearSample)) %>% 
  arrange(waterbodyID,yearSample)

statDat_wide_multTrophicCat_recent <- statDat_wide_multTrophicCat %>% 
  group_by(waterbodyID) %>% 
  filter(yearSample == max(yearSample))

## DOC - USE THIS; same model formulation as historic comparison

docTemporal <- lm(log(DOC)~ lakeTrophicStatus+
  scale(lat,scale=FALSE)+
  scale(maxDepth,scale=FALSE), data = statDat_wide_multTrophicCat)
plot(docTemporal)
summary(docTemporal)
anova(docTemporal)
#worthwhile to plot partial dependence plots

pTime_nuts <- ggplot(filter(statDat_wide_multTrophicCat,yearSample < 2020)) +
  geom_jitter(aes(x=yearSample,y=DOC, colour=lakeTrophicStatus, size=maxDepth,alpha=lat),width=0.2) +
  geom_smooth(aes(x=yearSample,y=DOC, colour=lakeTrophicStatus), method = "lm")
pTime_nuts

## General map of current DOC concentrations across Ontario
docContempMap <-ggplot() +
  geom_sf(data = ontario_shape, colour="grey", alpha = 0.3) +
  geom_point(data = statDat_wide_multTrophicCat_recent, aes(x = long, y = lat, colour = DOC,  size = maxDepth)) +
  scale_color_gradient2(mid = "white", low = "blue", high = "brown4", name = "DOC (mg/L)") +  # Specify the gradient from red to white to blue
  labs(size = "Max. lake depth (m)") +
  theme_minimal() +
  facet_wrap(~lakeTrophicStatus, nrow = 1)
docContempMap
# need to change gradient of lake depth point size  









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

#Secchi
BsM_secchiMod <- lmer(log(secchiDepth+0.1)~yearSample + (1|waterbodyID), statDat_wide_multTrophicCat) #log transformation necessary
plot(BsM_secchiMod)
summary(BsM_secchiMod)
anova(BsM_secchiMod)

plot(log(secchiDepth+0.1)~yearSample, allBsM_mult)
abline(lm(log(secchiDepth+0.1)~yearSample, allBsM_mult))

#DOC
BsM_DOCMod <- lmer(log(DOC)~yearSample + (1|waterbodyID), statDat_wide_multTrophicCat) #log transformation necessary
plot(BsM_DOCMod)
summary(BsM_DOCMod)
anova(BsM_DOCMod) #No significant difference, but trophic status not included. This has been shown to be important

#DOC
## - the issue here is that all the significant findings were on data w. eutrophic data removed
BsM_DOCMod <- lmer(log(DOC)~yearSample*lakeTrophicStatus+ + (1|waterbodyID), statDat_wide_multTrophicCat) #log transformation necessary
plot(BsM_DOCMod)
summary(BsM_DOCMod)
anova(BsM_DOCMod)

plot(DOC~yearSample, statDatWide_mult)
abline(lm(DOC~yearSample, statDatWide_mult))

