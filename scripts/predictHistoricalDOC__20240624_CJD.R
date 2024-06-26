
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
library(visreg)
library(here)

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

#What's best model? #CODY: updated so all use glmmTMB
secGamma_int <- glmmTMB(secchiDepth ~ scaled_yearSample + (1|commonID),
                        data=mainDF, family=Gamma(link="log"))

secGamma_slopeInt <- glmmTMB(secchiDepth ~ scaled_yearSample + (scaled_yearSample|commonID),
                             data=mainDF, family=Gamma(link="log"))

secLog_int <- update(secGamma_int, family = lognormal(link="log"))

secLog_slopeInt <- update(secGamma_slopeInt, family=lognormal(link="log"))

AIC(secGamma_int,secGamma_slopeInt,secLog_int,secLog_slopeInt)

# All converged. Gamma better than lognormal. Random slopes better than just random intercept but I still 
# prefer the random intercept for ease of interpretation - although fine with either

# Diagnostics: CODY rejigged this
bestmodel = secGamma_int
DHARMa::simulateResiduals(bestmodel, plot=T)
hist(residuals(bestmodel), breaks = 30, main = "Histogram of Residuals")
# Diagnostics look fine

# Let's look at results. CODY: I rejjiged this for glmmTMB models
(sumSecchi <- summary(bestmodel)) 
exp(sumSecchi$coefficients$cond[, "Estimate"])#secchi decreases by 0.37% per year (1-0.9963)
visreg(bestmodel, "scaled_yearSample", scale="response")
r.squaredGLMM(bestmodel)          # Fixed effects R2 = 0.017; random R2 = 0.75

####################################################################

# 2) Can we predict DOC w. lake variables? Secchi, depth, lat? ####
mainDF_BsM <- filter(mainDF, samplingProgram %in% "BsM")

# Fit the initial model - create full models with log-normal and gamma, look at AIC
clarDOCmod_gamma <- glmmTMB(DOC ~ scaled_secchi * scaled_maxDepth + lat + (1|waterbodyID), 
                            data = mainDF_BsM, family = Gamma(link = "log"), na.action = "na.fail")   # added lat in prediction because predictions w. only secchi and maxDepth under-predicted high DOC lakes. Only included as additive effect as there's a S/N gradient in DOC across the province
(clarDOCmod_gamma_sel <- dredge(clarDOCmod_gamma)) #additive model best
clarDOCmod_LMM<- glmmTMB(DOC ~ scaled_secchi * scaled_maxDepth + lat + (1|waterbodyID), 
                            data = mainDF_BsM,  na.action = "na.fail")    
(clarDOCmod_LMM_sel <- dredge(clarDOCmod_LMM)) #interactive model best
clarDOCmod_logLMM<- glmmTMB(DOC ~ scaled_secchi * scaled_maxDepth + lat + (1|waterbodyID), 
                      data = mainDF_BsM, family = lognormal(link="log"),  na.action = "na.fail")    
(clarDOCmod_logLMM_sel <- dredge(clarDOCmod_logLMM)) #additive model best
# AICc lowest with gamma. Do that.

bestmodel2 <- glmmTMB(DOC ~ scaled_secchi + scaled_maxDepth + lat + (1|waterbodyID),
                              data = mainDF_BsM, family = Gamma(link = "log"))

# Diagnostics
DHARMa::simulateResiduals(bestmodel2, plot=T)
hist(residuals(bestmodel2), breaks = 30, main = "Histogram of Residuals")
# Diagnostics look fine

# Summary
summary(bestmodel2)
r.squaredGLMM(bestmodel2)        # R2 = 0.51 fixed effects, R2 = 0.92 w. random effects 

# Add predicted values to df
predictedVals <- predict(bestmodel2, type = "response", re.form=NA) #for discussion (should random effects be used or not)
mainDF_BsM$predictedDOC <- predictedVals
cor(mainDF_BsM$DOC, mainDF_BsM$predictedDOC) 

## Relationship between observed and predicted DOC for BSM data
pObsFit <- ggplot(mainDF_BsM) +
  geom_point(aes(x=DOC,y=predictedDOC, colour = lat, size=maxDepth), alpha=0.7) +    # exp y-axis to get back on raw scale
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  scale_colour_viridis_c(option="inferno", direction = -1) +
  theme_bw()
pObsFit


# 3) Estimate DOC in ARU data ####  Cody's way (maybe simpler)
mainDF$rowID = 1:nrow(mainDF)
newdata = mainDF %>% filter(samplingProgram == "ARU")
newdata$ARU_DOC_predicted = predict(bestmodel2, newdata= newdata, type = "response", re.form = NA) #prediction only on fixed effects

updatedDOC = NA
for(i in 1:nrow(mainDF)){
  if(mainDF$samplingProgram[i] == "ARU"){ 
    updatedDOC[i] = newdata$ARU_DOC_predicted[which(newdata$rowID == i)]
  } else {
    updatedDOC[i] = mainDF$DOC[i]
  }
}

mainDF$updatedDOC = updatedDOC

##should we drop lakes that aren't repeat sampled?? 
##probably doesn't matter but might be worth doing just in case

# 4) Are there differences in DOC over time? ####
pDOCTime <- ggplot(data=mainDF,
                      aes(x=yearSample,y=updatedDOC)) +
  geom_point(alpha=0.7) +
  geom_smooth(method = "lm") +
  theme_bw()
pDOCTime

# Using IC approach, assemble different models to look for changes in DOC over time
M1 <- glmmTMB(updatedDOC ~ scaled_yearSample + (1|commonID),
                                        data = mainDF,
                                        family = Gamma(link = "log"))
M2 <- glmmTMB(updatedDOC ~ scaled_yearSample + (scaled_yearSample|commonID),
                                      data = mainDF,
                                      family = Gamma(link = "log"))          # MODEL DIDN'T CONVERGE
M3 = update(M1, family = lognormal(link="log"))
M4 = update(M2, family = lognormal(link="log"))
AIC(M1,M2,M3,M4)

# Gamma models are best. 

# Parameterize best model for DOC over time
bestmodel3 = M2

# Diagnostics
DHARMa::simulateResiduals(bestmodel3, plot=T)
hist(residuals(bestmodel3), breaks = 30, main = "Histogram of Residuals")
# Diagnostics look fine

# Model summaries
(sumDOC <- summary(bestmodel3)) 
exp(sumDOC$coefficients$cond[, "Estimate"]) #DOC increases by 0.06% per year
visreg(bestmodel3, "scaled_yearSample", scale="response")
r.squaredGLMM(bestmodel3)          # pretty much all lake effect



###CODY: I HAVE NOT UPDATED BELOW




# 5) Using BsM data (the best we got), what predicts DOC across Ontario? Calculate Lake Trophic Status using TDP ####
mainDF_onlyBsM <- mainDF %>% 
  filter(samplingProgram == "BsM") %>% 
  mutate(lakeTrophicStatus = case_when(TDP < 5 ~ "oligotrophic",
                                       TDP >= 5 & TDP < 10 ~ "mesotrophic",
                                       TDP >= 10 ~ "eutrophic"),
         scaledTDP = as.vector(scale(TDP, scale = FALSE)))
mainDF_onlyBsM$lakeTrophicStatus <- factor(mainDF_onlyBsM$lakeTrophicStatus, 
                                           levels = c("oligotrophic","mesotrophic","eutrophic"))  # just changing order for plotting

# Log-normal models perform much better, just start w. full log-normal model for IC
docSpatial <- glmmTMB(DOC ~ scaled_yearSample*scaled_maxDepth*scaledTDP*lat + (1|waterbodyID), 
                   data=mainDF_onlyBsM, na.action = "na.fail", family=Gamma(link="log"))
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
