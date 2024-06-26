
## New setup: 
# 1) Change in water clarity within the two sampling periods?
# 2) Can we predict DOC w. lake variables? Secchi, depth, lat?
# 3) Estimate DOC in ARU data
# 4) Are there differences in DOC over time?
# 5) Using BsM data, what predicts DOC across Ontario? Calculate Lake Trophic Status using TDP
# 6) Map + show correlations w. current lake variables and DOC

library(tidyverse)
library(lme4)
library(lmerTest)
library(glmmTMB)
library(MuMIn)

## Read Data ####
combIDs <- read_csv(here("data/commonDataIDs_aruToBsm_20240307.csv")) #has common IDs between BsM and ARU
allBsm <- read_csv("data/allBsM_updatedCyc3_20240520.csv") %>%
  left_join(select(combIDs, bsmWaterbodyID, commonID), by = c("waterbodyID" = "bsmWaterbodyID"))
allARU <- read_csv(here("data/ARU_modelData_20240311.csv")) %>% 
  left_join(select(combIDs, ARU_LID, commonID), by = "ARU_LID") %>% 
  filter(TDS < 300) #Different linear relationship for saline lakes in Chow-Fraser (1991); only ~50/9000 observations being removed
statDat_wide <- read_csv("data/statDatWide_updatedCyc3_20240520.csv") %>% 
  filter(log(TKN) >4) %>%  # Clear outliers < 4 (only two observations)                       
  left_join(select(combIDs, bsmWaterbodyID, commonID), by = c("waterbodyID" = "bsmWaterbodyID"))

# Estimate TP from TDS and apply a trophic status
aruDat <- allARU %>% 
  left_join(select(combIDs,ARU_LID,commonID)) %>% 
  mutate(logTP = 0.95*log(TDS)-0.669,                                         #This formula is from Chow-Fraser (1991)
         TP = exp(logTP),
         lakeTrophicStatus = case_when(TP < 10 ~ "oligotrophic",
                                       TP >= 10 & TP < 20 ~ "mesotrophic",
                                       TP >= 20 ~ "eutrophic")) %>% 
  filter(!is.na(commonID)) %>% 
  rename(secchiDepth = SECCHI,yearSample = INV_YR, maxDepth = DEPMAX, lat=latARU,long=longARU)

# Main dataset combining ARU and BsM
mainDF <- bind_rows(statDat_wide,aruDat) %>% 
  filter(!is.na(yearSample))
todrop <- mainDF %>% group_by(commonID) %>% summarize(count = n()) %>% arrange(count) %>% filter(count == 1)
mainDF <- mainDF[-which(mainDF$commonID %in% todrop$commonID),]
mainDF$scaled_yearSample <- as.vector(scale(mainDF$yearSample, scale = FALSE)) # Create scaledYear. scale=F so that a 1 unit change still represents 1 year 
mainDF$scaled_secchi <- as.vector(scale(mainDF$secchiDepth, scale = FALSE)) # Create scaledYear. scale=F so that a 1 unit change still represents 1 year 
mainDF$scaled_maxDepth <- as.vector(scale(mainDF$maxDepth, scale = FALSE)) # Create scaledYear. scale=F so that a 1 unit change still represents 1 year 
mainDF$samplingProgram <- if_else(mainDF$scaled_yearSample >= 0, "BsM", "ARU")

# 1) Change in water clarity within the two sampling periods? ####
pSecchiTime <- ggplot(data=mainDF,
                      aes(x=yearSample,y=secchiDepth)) +
  geom_point(alpha=0.7) +
  geom_smooth(method = "lm") +
  theme_bw()
pSecchiTime

#What's best model?
secGamma_int <- glmmTMB(secchiDepth ~ scaled_yearSample + (1|commonID),
                        data=mainDF, family=Gamma(link=log))
secGamma_slopeInt <- glmmTMB(secchiDepth ~ scaled_yearSample + (scaled_yearSample|commonID),
                             data=mainDF, family=Gamma(link=log))
secLog_int <- lmer(log(secchiDepth) ~ scaled_yearSample + (1|commonID),
                   data=mainDF)
secLog_slopeInt <- lmer(log(secchiDepth) ~ scaled_yearSample + (scaled_yearSample|commonID),
                   data=mainDF)
AIC(secGamma_int,secGamma_slopeInt,secLog_int,secLog_slopeInt)
# Slopes/Intercept model failed to converge, secLog_int has lowest AIC for reasonable, converging model

# Diagnostics
plot(secLog_int)
qqnorm(residuals(secLog_int))
qqline(residuals(secLog_int))
hist(residuals(secLog_int), breaks = 30, main = "Histogram of Residuals")
lattice::dotplot(ranef(secLog_int, condVar = TRUE), scales = list(relation = "free"))
# Diagnostics look fine

# Let's look at results
(sumSecchi <- summary(secLog_int)) # Highly significant; secchi decreases 0.035 m/year; secchi has
                                   # decreased 2.17m over range of years (62 yrs)
r.squaredGLMM(secLog_int)          # Fixed effects R2 = 0.015; random R2 = 0.74

# 2) Can we predict DOC w. lake variables? Secchi, depth, lat? ####
mainDF_BsM <- filter(mainDF, samplingProgram %in% "BsM")

# Fit the initial model - create full models with log-normal and gamma, look at AIC
clarDOCmod_gamma <- glmmTMB(DOC ~ scaled_secchi * scaled_maxDepth + lat + (1|waterbodyID), 
                            data = mainDF_BsM, family = Gamma(link = "log"), na.action = "na.fail")   # added lat in prediction because predictions w. only secchi and maxDepth under-predicted high DOC lakes. Only included as additive effect as there's a S/N gradient in DOC across the province
(clarDOCmod_gamma_sel <- dredge(clarDOCmod_gamma)) #additive model best
clarDOCmod_LMM<- lmer(DOC ~ scaled_secchi * scaled_maxDepth + lat + (1|waterbodyID), 
                            data = mainDF_BsM,  na.action = "na.fail")    
(clarDOCmod_LMM_sel <- dredge(clarDOCmod_LMM)) #interactive model best
clarDOCmod_logLMM<- lmer(log(DOC) ~ scaled_secchi * scaled_maxDepth + lat + (1|waterbodyID), 
                      data = mainDF_BsM,  na.action = "na.fail")    
(clarDOCmod_logLMM_sel <- dredge(clarDOCmod_logLMM)) #additive model best
# AICc lowest with log-linear model. Do that.

clarDOCmod_linearBest <- lmer(log(DOC) ~ scaled_secchi + scaled_maxDepth + lat + (1|waterbodyID),
                              data = mainDF_BsM)

# Diagnostics
plot(clarDOCmod_linearBest)
qqnorm(residuals(clarDOCmod_linearBest))
qqline(residuals(clarDOCmod_linearBest))
hist(residuals(clarDOCmod_linearBest), breaks = 30, main = "Histogram of Residuals")
lattice::dotplot(ranef(clarDOCmod_linearBest, condVar = TRUE), scales = list(relation = "free"))
# Diagnostics look fine

# Summary
summary(clarDOCmod_linearBest)
r.squaredGLMM(clarDOCmod_linearBest)        # R2 = 0.51 fixed effects, R2 = 0.92 w. random effects 

# Add predicted values to df
predictedVals <- predict(clarDOCmod_linearBest, type = "response")
mainDF_BsM$predictedDOC <- predictedVals

## This shows a significant association of DOC w. water clarity. Use coefficients to predict DOC
pObsFit <- ggplot(mainDF_BsM) +
  geom_point(aes(x=DOC,y=exp(predictedDOC), colour = lat, size=maxDepth), alpha=0.7) +    # exp y-axis to get back on raw scale
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  scale_colour_viridis_c(option="inferno", direction = -1) +
  theme_bw()
pObsFit

# 3) Estimate DOC in ARU data #### 
docFixed <- fixef(clarDOCmod_linearBest)                        # Extract fixed effect coefs
docModCoefs <- as.data.frame(ranef(clarDOCmod_linearBest)) %>%  # Extract random intercepts for each lake
  select(grp,condval) %>%                                       # subset only waterbodyID and random intercept offset
  rename(waterbodyID = grp, interceptR = condval) %>%           # rename variables
  mutate(intercept = docFixed[1] + interceptR,                  # For each lake add the random intercept offset
         scaledSecchi_coef = docFixed[2],
         scaledMaxDepth_coef = docFixed[3],
         lat_coef = docFixed[4]) %>% 
  left_join(select(combIDs,bsmWaterbodyID,commonID),by=c("waterbodyID"="bsmWaterbodyID")) %>%   # add commonID for later sorting
  filter(!is.na(commonID)) 

# Coef DF using GLMM
# **Cody, not sure if the 'mainDF_estimateDOC' is legit? In most cases we'd just be interested in estimating a response
# at the population level (so just coefficients for the fixed effects, as is done in mainDF_estimateDOC) but in this
# case it seems useful to estimate lake-specific intercepts using information from the random intercept. Is estimating
# an intercept for each lake (offset from grand mean) legit? Seems like it should be, but I've always just estimated 
# coefs from the fixed effects. Furthermore, the obs. vs fitted plot above is taking random intercepts into account
mainDF_estimateDOC <- mainDF %>% 
  group_by(commonID) %>% 
  mutate(updatedDOC = case_when(samplingProgram == "BsM" ~ DOC,       # For BsM data, leave as measured DOC
                            samplingProgram == "ARU" ~                # For ARU data, calculate DOC based off model coefs. 
                              exp(docFixed[1] +                       # Exp estimate to get on raw scale
                              docFixed[2] * scaled_secchi +
                              docFixed[3] * scaled_maxDepth +
                              docFixed[4] * lat))) %>% 
  select(lakeName,commonID,waterbodyID,samplingProgram,yearSample,lat,long,maxDepth,meanDepth,surfaceArea,secchiDepth,TDP,DOC,
         TP,lakeTrophicStatus,scaled_secchi,scaled_maxDepth,scaled_yearSample,updatedDOC)

mainDF_estimateDOC_wRandInt <- mainDF %>% 
  left_join(select(docModCoefs,intercept,scaledSecchi_coef,scaledMaxDepth_coef,lat_coef,commonID), by="commonID") %>% 
  group_by(commonID) %>% 
  mutate(updatedDOC = case_when(samplingProgram == "BsM" ~ DOC,      
                            samplingProgram == "ARU" ~ 
                              exp(intercept +
                              scaledSecchi_coef * scaled_secchi +
                              scaledMaxDepth_coef * scaled_maxDepth +
                              lat_coef * lat))) %>% 
  select(lakeName,commonID,waterbodyID,samplingProgram,yearSample,lat,long,maxDepth,meanDepth,surfaceArea,secchiDepth,TDP,DOC,
         TP,lakeTrophicStatus,scaled_secchi,scaled_maxDepth,scaled_yearSample,updatedDOC)

# 4) Are there differences in DOC over time? ####
pDOCTime <- ggplot(data=mainDF_estimateDOC_wRandInt,
                      aes(x=yearSample,y=updatedDOC)) +
  geom_point(alpha=0.7) +
  geom_smooth(method = "lm") +
  theme_bw()
pDOCTime

# Using IC approach, assemble different models to look for changes in DOC over time
docChangeMod_noSlope_glmmTMB <- glmmTMB(updatedDOC ~ scaled_yearSample + (1|commonID),
                                        data = mainDF_estimateDOC_wRandInt,
                                        family = Gamma(link = "log"))
docChangeMod_Slope_glmmTMB <- glmmTMB(updatedDOC ~ scaled_yearSample + (scaled_yearSample|commonID),
                                      data = mainDF_estimateDOC_wRandInt,
                                      family = Gamma(link = "log"))          # MODEL DIDN'T CONVERGE
docChangeMod_lm_noSl <- lmer(updatedDOC ~ scaled_yearSample + (1|commonID), data = mainDF_estimateDOC_wRandInt)
docChangeMod_lm <- lmer(updatedDOC ~ scaled_yearSample + (scaled_yearSample|commonID), data = mainDF_estimateDOC_wRandInt)  # MODEL SINGULAR FIT
docChangeMod_lm_noSl_log <- lmer(log(updatedDOC) ~ scaled_yearSample+ (1|commonID), data = mainDF_estimateDOC_wRandInt)
docChangeMod_lm_log <- lmer(log(updatedDOC) ~ scaled_yearSample + (scaled_yearSample|commonID), data = mainDF_estimateDOC_wRandInt) # MODEL SINGULAR FIT
AIC(docChangeMod_noSlope_glmmTMB,docChangeMod_Slope_glmmTMB,docChangeMod_lm_noSl,docChangeMod_lm,docChangeMod_lm_noSl_log,docChangeMod_lm_log)
# Best model by far is log-normal random intercept model

# Parameterize best model for DOC over time
bestDOC_timeMod <- lmer(log(updatedDOC) ~ scaled_yearSample+ (1|commonID), data = mainDF_estimateDOC_wRandInt)

# Diagnostics
plot(bestDOC_timeMod)
qqnorm(residuals(bestDOC_timeMod))
qqline(residuals(bestDOC_timeMod))
hist(residuals(bestDOC_timeMod), breaks = 30, main = "Histogram of Residuals")
lattice::dotplot(ranef(bestDOC_timeMod, condVar = TRUE), scales = list(relation = "free"))
# Diagnostics look fine

# Model summaries
summary(bestDOC_timeMod)              # Strongly significant
r.squaredGLMM(bestDOC_timeMod)        # R2 = 0.0016 fixed effects, R2 = 0.95 w. random effects 
# CODY - what's the math for sorting out what the raw estimate from the scaled coefficent again?
# 0.001 mg/L per year isn't right, is it? that would mean 0.062 change over 62 years...which based
# on the plot of DOC over time can't be right (even though the mean hasn't changed that much)

# 5) Using BsM data (the best we got), what predicts DOC across Ontario? Calculate Lake Trophic Status using TDP ####
mainDF_onlyBsM <- mainDF_estimateDOC_wRandInt %>% 
  filter(samplingProgram == "BsM") %>% 
  mutate(lakeTrophicStatus = case_when(TDP < 5 ~ "oligotrophic",
                                       TDP >= 5 & TDP < 10 ~ "mesotrophic",
                                       TDP >= 10 ~ "eutrophic"),
         scaledTDP = as.vector(scale(TDP, scale = FALSE)))
mainDF_onlyBsM$lakeTrophicStatus <- factor(mainDF_onlyBsM$lakeTrophicStatus, 
                                           levels = c("oligotrophic","mesotrophic","eutrophic"))  # just changing order for plotting

# Log-normal models perform much better, just start w. full log-normal model for IC
docSpatial <- lmer(log(DOC) ~ scaled_yearSample*scaled_maxDepth*scaledTDP*lat + (1|waterbodyID), 
                   data=mainDF_onlyBsM, na.action = "na.fail")
dredge(docSpatial)

# Also see what trophic status says rather than continous TDP
docSpatial_catTrophicStatus <- lmer(log(DOC) ~ scaled_yearSample*scaled_maxDepth*lakeTrophicStatus*lat + (1|waterbodyID), 
                                    data=mainDF_onlyBsM, na.action = "na.fail")
dredge(docSpatial_catTrophicStatus)

# Overwhelmingly, best model is simple additive model
docSpatialBest <- lmer(log(DOC) ~ lat + scaled_maxDepth + scaledTDP
                       + (1|waterbodyID), data=mainDF_onlyBsM)
summary(docSpatialBest)
r.squaredGLMM(docSpatialBest)

docSpatialBest_cat <- lmer(log(DOC) ~ lat + scaled_maxDepth + lakeTrophicStatus
                       + (1|waterbodyID), data=mainDF_onlyBsM)
summary(docSpatialBest_cat)
r.squaredGLMM(docSpatialBest_cat)

pSpatialDOC<- ggplot(mainDF_onlyBsM, aes(x=yearSample,y=DOC,colour=lat)) +
  geom_point(aes(size=maxDepth), alpha=0.7) +    # exp y-axis to get back on raw scale
  scale_colour_viridis_c(option="inferno", direction = -1) +
  theme_bw()+
  facet_wrap(~lakeTrophicStatus)
pSpatialDOC


# 6) Map + show correlations w. current lake variables and DOC ####



# EXTRA ####

# ARU
plot(aruSecchi~yearSample, aruDat)

aruClarMod <- lm(aruSecchi~yearSample, aruDat)
plot(aruClarMod)
aruClarMod_log <- lm(log(aruSecchi)~yearSample, aruDat)
plot(aruClarMod_log)
aruClarMod_gamma <- glm(aruSecchi~yearSample, family = Gamma(link = "log"), aruDat)
plot(aruClarMod_gamma)
AIC(aruClarMod, aruClarMod_log, aruClarMod_gamma)     # log model definitely best

(sumClarARU <- summary(aruClarMod_log))  # Significant INCREASE in clarity, though R2 basically 0
# So no, no decrease in water clarity

# BsM
plot(meanSecchi ~ yearSample, bsmDOC)
abline(lm(meanSecchi ~ yearSample, bsmDOC))  # Decrease in clarity over time
pSecchiBsM <- ggplot(bsmDOC, aes(x=yearSample,y=meanSecchi)) +
  geom_point(aes(colour=lat, size = maxDepth), alpha=0.6) +
  geom_smooth(method = "lm")
pSecchiBsM

bsmClarMod <- lmer(meanSecchi~yearSample + (1|commonID), bsmDOC)
plot(bsmClarMod)
bsmClarMod_log <- lmer(log(meanSecchi)~yearSample + (1|commonID), bsmDOC)
plot(bsmClarMod_log)
bsmClarMod_gamma <- glmmTMB(meanSecchi~yearSample + (1|commonID), family = Gamma(link = "log"), data=bsmDOC)
plot(bsmClarMod_gamma)
AIC(bsmClarMod, bsmClarMod_log, bsmClarMod_gamma)     # log model definitely best

(sumClarBsM <- summary(bsmClarMod_log))  # Significant decrease in clarity
r.squaredGLMM(bsmClarMod_log) # Marginal R2 < 1%; Conditional is 78% so lots of within-lake variance
# Water clarity is decreasing 1cm per year between 2008-2023 (so 15cm since 2008)

# What about DOC
plot(meanDOC ~ yearSample, bsmDOC)
abline(lm(meanDOC ~ yearSample, bsmDOC))  # Decrease in clarity over time

bsmDOCMod <- lmer(meanDOC~yearSample + (1|commonID), bsmDOC)
plot(bsmDOCMod)
bsmDOCMod_log <- lmer(log(meanDOC)~yearSample + (1|commonID), bsmDOC)
plot(bsmDOCMod_log)
bsmDOCMod_gamma <- glmmTMB(meanDOC~yearSample + (1|commonID), family = Gamma(link = "log"), data=bsmDOC)
plot(bsmDOCMod_gamma)
AIC(bsmDOCMod, bsmDOCMod_log, bsmDOCMod_gamma)     # log model definitely best

(sumClarBsM <- summary(bsmDOCMod_log))  # Significant decrease in clarity, though R2 basically 0
# Water clarity is decreasing 1cm per year between 2008-2023 (so 15cm since 2008)


