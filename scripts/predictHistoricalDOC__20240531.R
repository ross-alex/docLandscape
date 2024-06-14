## Water clarity differences amongst Ontario lakes: 1960s to present
# This script is starting clean with assessing temporal changes in a continous 
# rather than discrete manner (i.e, 1960s -> 2023 rather than historic vs. contemporary)

## Analysis goals
# 1) Is water clarity (i.e., Secchi) associated w. DOC
# 2) Predict DOC using water clarity
# 3) Has DOC increased over time
# 4) Given changes in DOC, what lakescape predictors best describe variation across Ontario lakes
# 5) In appendix, confirm complimentary changes in #2 and #3 using water clarity (did it decrease?)

library(tidyverse)
library(lubridate)
library(MuMIn)
library(lme4)
library(lmerTest)
library(scales)
library(sf)

## Read Data ####
combIDs <- read_csv(here("data/commonDataIDs_aruToBsm_20240307.csv")) #has common IDs between BsM and ARU
allBsm <- read_csv("data/allBsM_updatedCyc3_20240520.csv") %>%
  left_join(select(combIDs, bsmWaterbodyID, commonID), by = c("waterbodyID" = "bsmWaterbodyID"))
allARU <- read_csv(here("data/ARU_modelData_20240311.csv")) %>% 
  left_join(select(combIDs, ARU_LID, commonID), by = "ARU_LID")
statDat_wide <- read_csv("data/statDatWide_updatedCyc3_20240520.csv") %>% 
  filter(log(TKN) >4) # Clear outliers < 4 (only two observations)                       

# Estimate TP from TDS and apply a trophic status
aruDat <- allARU %>% 
  left_join(select(combIDs,ARU_LID,commonID)) %>% 
  mutate(logTP = 0.95*log(TDS)-0.669,                                         #This formula is from Chow-Fraser (1991)
         TP = exp(logTP),
         lakeTrophicStatus = case_when(TP < 10 ~ "oligotrophic",
                                       TP >= 10 & TP < 20 ~ "mesotrophic",
                                       TP >= 20 ~ "eutrophic")) %>% 
  filter(!is.na(commonID)) %>% 
  rename(aruSecchi = SECCHI,yearSample = INV_YR)

# Create BsM dataframe with only DOC as predictor variable
bsmDOC <- filter(allBsm, predictorVariable == "DOC") %>% 
  filter(!is.na(secchiDepth) & !is.na(maxDepth) & !is.na(predictorValue)) %>%  # Lost 319 values (many Cycle 3 entries incomplete)
  rename(DOC = predictorValue)


## Goal 1 - Is water clarity associated w. DOC? ####
# Was determined that DOC is gamma distributed; use glm

# Fit the initial model - note: This chunk has outliers, use next
clarDOCmod_gamma <- glm(DOC ~ secchiDepth * maxDepth + lat, data = bsmDOC, family = Gamma(link = "log"), na.action = "na.fail")   # added lat in prediction because predictions w. only secchi and maxDepth under-predicted high DOC lakes. Only included as additive effect as there's a S/N gradient in DOC across the province
(clarDOCmod_gamma_sel <- dredge(clarDOCmod_gamma)) #additive model best

# Best model
clarDOCmod_gammaBest <- glm(DOC ~ scale(secchiDepth) * scale(maxDepth) + scale(lat), data = bsmDOC, family = Gamma(link = "log"))

#  Model diagnostics
DHARMa::simulateResiduals(clarDOCmod_gammaBest, plot=T)
plot(clarDOCmod_gammaBest) # Diagnositcs not great; outliers evident
car::vif(clarDOCmod_gammaBest) # Only 2 vars, ok
# model summaries
summary(clarDOCmod_gammaBest)

## Gamma dist. has 2 clear outlier; rows 629 and 1509. Reduce dataset try again (629 becomes even more apparent if only 1509 is removed)
bsmDOCGammaRedDat <- bsmDOC[c(-629, -1509),]

# Fit the initial model
clarDOCmod_gamma_red <- glm(DOC ~ secchiDepth * maxDepth + lat, data = bsmDOCGammaRedDat, family = Gamma(link = "log"), na.action = "na.fail")
(clarDOCmod_gamma_sel_red <- dredge(clarDOCmod_gamma_red)) #interaction is now best

# Best reduced model
clarDOCmod_gammaBest_red <- glm(DOC ~ scale(secchiDepth) * scale(maxDepth) + lat, data = bsmDOCGammaRedDat, family = Gamma(link = "log"))

#  Model diagnostics
DHARMa::simulateResiduals(clarDOCmod_gammaBest_red, plot=T) #Better
plot(clarDOCmod_gammaBest_red) # Diagnositcs not great but fine
car::vif(clarDOCmod_gammaBest_red) # VIFs fine <10

# model summaries
(glmSum <- summary(clarDOCmod_gammaBest_red))
# Deviance-based R2 - seems legit..
deviance_model <- deviance(glmSum)
deviance_null <- deviance(update(glmSum, . ~ 1))
(R2_deviance <- 1 - (deviance_model / deviance_null)) # R2 0.63
# Add predicted values to df
predictedVals <- predict(clarDOCmod_gammaBest_red, type = "response")
bsmDOCGammaRedDat$predictedDOC <- predictedVals

## This shows a significant association of DOC w. water clarity. Use coefficients to predict DOC
pObsFit <- ggplot(bsmDOCGammaRedDat) +
  geom_point(aes(x=DOC,y=predictedDOC, colour = lat)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  theme_bw()
pObsFit

## Goal 2 - Predict historical ARU DOC with BsM eqn coefs ####
# With ARU data, estimate TP from TDS and DOC
(aruDOC_coefs_GLM <- coef(clarDOCmod_gammaBest_red)) #GLM coefs 

# Coef DF using GLM
aruDat_DOCestimated_GLM <- aruDat %>% 
  rename(secchiDepth = aruSecchi, maxDepth = DEPMAX, lat = latARU) %>% 
  filter(TDS < 300) %>% #Different linear relationship for saline lakes in Chow-Fraser (1991)
  mutate(DOC = predict(clarDOCmod_gammaBest_red, newdata = .,type="response")) %>% 
  select(commonID,yearSample,DOC,secchiDepth,maxDepth,lat,longARU,lakeTrophicStatus) %>% 
  rename(long=longARU)

## Goal 3 - Has DOC changed over time? ####

## Mixed-model approach
# Create a common dataframe w. ARU and BsM data. Need:
# yearSampled, secchi, DOC, maxDepth, trophicStatus, commonID/waterbodyID, lat/lon
# *DOC between 1962-1986 is from ARU; estimated using secchi and maxDepth. BsM is 2008-2023
mainDF <- bsmDOC %>% 
  select(commonID,yearSample,DOC,secchiDepth,maxDepth,lat,long) %>% 
  left_join(select(aruDat_DOCestimated_GLM,commonID,lakeTrophicStatus)) %>% 
  bind_rows(aruDat_DOCestimated_GLM) %>% 
  filter(complete.cases(.))

# Distribution of response variable
(docDist_raw <- hist(mainDF$DOC))
(docDist_log <- hist(log(mainDF$DOC)))
(docGamma <- MASS::fitdistr(mainDF$DOC, densfun = "gamma"))
# Visualize gamma fit
ggplot(mainDF, aes(x = DOC)) +
  geom_histogram(aes(y = ..density..), bins = 30, fill = "skyblue", color = "black", alpha = 0.7) +
  stat_function(fun = dgamma, args = list(shape = docGamma$estimate["shape"], rate = docGamma$estimate["rate"]), color = "red", size = 1) +
  labs(x = "DOC", y = "Density") +
  ggtitle("Gamma Distribution Fit to Data") 
# gamma fits well with shape = 5.23, rate = 0.74

# What do the data look like?
plot(DOC~yearSample, mainDF)

# Model w. gamma distribution
docChangeMod_noSlope <- glmer(DOC ~ scale(yearSample) + (1|commonID), data = mainDF, 
                      family = Gamma(link = "log"))
docChangeMod <- glmer(DOC ~ scale(yearSample) + (yearSample|commonID), data = mainDF, 
                      family = Gamma(link = "log"))
docChangeModNull <- glmer(DOC ~ 1 + (yearSample|commonID), data = mainDF, 
                          family = Gamma(link = "log"))
AIC(docChangeMod_noSlope,docChangeMod,docChangeModNull)                         # Null has best AIC by a mile

DHARMa::simulateResiduals(docChangeMod, plot=T) # Absolutely horrible; can't be used
plot(docChangeMod) # Super bad; pattern
summary(docChangeMod)

# Model w. linear distribution
docChangeMod_lm_noSl <- lmer(DOC ~ scale(yearSample) + (1|commonID), data = mainDF)
docChangeMod_lm <- lmer(DOC ~ scale(yearSample) + (yearSample|commonID), data = mainDF)
docChangeMod_lmNull <- lmer(DOC ~ 1 + (yearSample|commonID), data = mainDF)
AIC(docChangeMod_lm, docChangeMod_lm_noSl, docChangeMod_lmNull)   # Model w. year as random slope best
 
DHARMa::simulateResiduals(docChangeMod_lm, plot=T) #Better
plot(docChangeMod_lm)    # some exploding variance but much better than; logging doesn't help, just flips the problem

# View the summary of the model
summary(docChangeMod_lm)            # Significant. Shows a 0.12 mg/L increase per decade
r.squaredGLMM(docChangeMod_lm)      # R2 quite telling. Marginal R2 = 0.006, Conditional R2 = 0.76. Nearly all variance 
                                    # captured at the lake level. Might not be much of an increase OR recent BsM period 
                                    # has no change. See chunk ~30 lines below to look at that (bsmModDat)

plot(DOC~yearSample,mainDF)
abline(lm(DOC~yearSample,mainDF))   # Confirms temporal relationship (or lack thereof)

## Goal 4 - If not changes over time, what explains DOC across Ontario lakes? ####
# Earlier analyses and literature indicates DOC may be affected by other easily accessed variables.
# When included in models, do we improve our understanding of DOC variation across the province?
docChangeMod_allVars_lmm <- lmer(DOC ~ yearSample*lakeTrophicStatus*maxDepth + (yearSample|commonID),         # latitude has been removed as it is now in predictions of historic DOC (can't also be a predictor...)
                                 data = mainDF, na.action = na.fail)
(allVars_modSel <- dredge(docChangeMod_allVars_lmm))   # Best model now just maxDepth and yearSample
                                                       

# docChangeMod_redVars_lmm <- lmer(DOC ~ yearSample+lakeTrophicStatus+maxDepth+lat+             # This model was best prior to removing
#                                    lat:yearSample + (yearSample|commonID),data = mainDF)      # latitude from the model

docChangeMod_redVars_lmm <- lmer(DOC ~ yearSample + maxDepth + (yearSample|commonID),data = mainDF)

DHARMa::simulateResiduals(docChangeMod_redVars_lmm, plot=T) #Ok
plot(docChangeMod_redVars_lmm) 

summary(docChangeMod_redVars_lmm)
r.squaredGLMM(docChangeMod_redVars_lmm)      # R2 better. Marginal R2 = 0.16, Conditional R2 = 0.82. 38% improvement 
                                             # from yearSampled only model

pTime <- ggplot(mainDF) +
  geom_point(aes(x=yearSample,y=DOC, size=lat,colour=maxDepth),alpha=0.7) +
  geom_smooth(aes(x=yearSample,y=DOC), method = "lm") +
  scale_color_viridis_c(option = "plasma") +
  theme_bw()
pTime


## What happens when we just look at BsM data
bsmModDat <- filter(mainDF, yearSample > 1990)  # No ARU data sampled in 90s

bsmModYear <- lmer(DOC~scale(yearSample) + (1|commonID), data=bsmModDat)
bsmModYearNull <- lmer(DOC~ 1 + (1|commonID), data=bsmModDat)
AIC(bsmModYear,bsmModYearNull)                                          # Null is best model; for shits, look at model below
summary(bsmModYear) # No effect when considering just year (and random slope + intercept model isn't working - not enough data; can't fit slope to 1 point)
# No effect in contemporary data suggests that maybe we just take a mean of all the data

# What about when we include more in the model?
docChangeMod_BsM_allVars_lmm <- lmer(DOC ~ scale(yearSample)*lakeTrophicStatus*scale(maxDepth)*scale(lat) + (1|commonID),
                                        data = bsmModDat, na.action = na.fail)    # Couldn't use slopes + intercepts because number of rand. effects > num fixed observations
dredge(docChangeMod_BsM_allVars_lmm)                                              # same best model as all data
# dredge above shows yearSample not even important and should be dropped from model. Clear that there's no effect of year

#***This below was based on
# docChangeMod_BsM_redVars_lmm <- lmer(DOC ~ yearSample+lakeTrophicStatus+maxDepth+lat+
#                                    lat:yearSample + (1|commonID),data = bsmModDat)
# DHARMa::simulateResiduals(docChangeMod_BsM_redVars_lmm, plot=T) #Ok... but not great
# plot(docChangeMod_BsM_redVars_lmm) 
# 
# summary(docChangeMod_BsM_redVars_lmm)
# r.squaredGLMM(docChangeMod_BsM_redVars_lmm)      # R2 better. Marginal R2 = 0.4, Conditional R2 = 0.94. 40% improvement from 
#                                                  # from yearSampled only model; lots of in-lake variance

pTime_BsM <- ggplot(bsmModDat) +
  geom_point(aes(x=yearSample,y=DOC, colour=lakeTrophicStatus, size=lat,alpha=maxDepth)) +
  geom_smooth(aes(x=yearSample,y=DOC, colour=lakeTrophicStatus), method = "lm") +
  theme_bw()
pTime_BsM




