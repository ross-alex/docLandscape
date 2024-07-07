
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
  filter(!is.na(yearSample)) %>% filter(!is.na(commonID))
mainDF = mainDF %>% filter(!is.na(secchiDepth)) %>% filter(!is.na(maxDepth))

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

# All converged. Gamma better than lognormal. Random slopes better than just random intercept 
# 

# Diagnostics: CODY rejigged this
bestmodel = secGamma_slopeInt
DHARMa::simulateResiduals(bestmodel, plot=T)
hist(residuals(bestmodel), breaks = 30, main = "Histogram of Residuals")
# Diagnostics look fine

# Let's look at results. CODY: I rejjiged this for glmmTMB models
(sumSecchi <- summary(bestmodel)) 
exp(sumSecchi$coefficients$cond[, "Estimate"])#secchi decreases by 0.37% per year (1-0.9963)
visreg(bestmodel, "scaled_yearSample", scale="response", partial = TRUE) # I think
r.squaredGLMM(bestmodel)          # Fixed effects R2 = 0.016; random R2 = 0.78

#nicer plot
mainDF$secchi_fit = predict(bestmodel, type = "response", re.form = NA)
ggplot(mainDF)+
  geom_point(aes(x=yearSample, y = secchiDepth), alpha=0.7, size=2.5)+
  geom_smooth(aes(x=yearSample, y = secchi_fit), 
              method="glm",
              formula = y~x,
              method.args = list(family = Gamma(link = 'log')))+
  scale_x_continuous(name="Year")+
  scale_y_continuous(name = "Secchi depth (m)")+
  theme_minimal()


# 2) Can we predict DOC w. lake variables? Secchi, depth, lat? ####
mainDF_BsM <- filter(mainDF, samplingProgram %in% "BsM")

# Fit the initial model - create full models with log-normal and gamma, look at AIC
clarDOCmod_gamma <- glmmTMB(DOC ~ scaled_secchi * scaled_maxDepth + lat + (1|commonID), 
                            data = mainDF_BsM, family = Gamma(link = "log"), na.action = "na.fail")   # added lat in prediction because predictions w. only secchi and maxDepth under-predicted high DOC lakes. Only included as additive effect as there's a S/N gradient in DOC across the province
(clarDOCmod_gamma_sel <- dredge(clarDOCmod_gamma)) #additive model best
# clarDOCmod_LMM<- glmmTMB(DOC ~ scaled_secchi * scaled_maxDepth + lat + (1|commonID),
#                             data = mainDF_BsM,  na.action = "na.fail")
# (clarDOCmod_LMM_sel <- dredge(clarDOCmod_LMM))
# clarDOCmod_logLMM<- glmmTMB(DOC ~ scaled_secchi * scaled_maxDepth + lat + (1|commonID),
#                       data = mainDF_BsM, family = lognormal(link="log"),  na.action = "na.fail")
# (clarDOCmod_logLMM_sel <- dredge(clarDOCmod_logLMM))
# AICc lowest with gamma. Do that.

bestmodel2 <- glmmTMB(DOC ~ scaled_secchi + scaled_maxDepth + lat + (1|commonID),
                              data = mainDF_BsM, family = Gamma(link = "log"))

# Diagnostics
DHARMa::simulateResiduals(bestmodel2, plot=T)
hist(residuals(bestmodel2), breaks = 30, main = "Histogram of Residuals")
# Diagnostics look fine

# Summary
summary(bestmodel2)
r2_secchiMod <- r.squaredGLMM(bestmodel2)        # R2 = 0.49 fixed effects, R2 = 0.92 w. random effects 

# Partial variance explained 
# Fit reduced models
model_no_secchi <- update(bestmodel2, . ~ . - scaled_secchi)
model_no_depth <- update(bestmodel2, . ~ . - scaled_maxDepth)
model_no_lat <- update(bestmodel2, . ~ . - lat)

# Calculate R² for reduced models
r2_no_secchi <- r.squaredGLMM(model_no_secchi)
r2_no_depth <- r.squaredGLMM(model_no_depth)
r2_no_lat <- r.squaredGLMM(model_no_lat)

# Calculate partial R²
partial_R2_secchi <- r2_secchiMod[1] - r2_no_secchi[1]
partial_R2_depth <- r2_secchiMod[1] - r2_no_depth[1]
partial_R2_lat <- r2_secchiMod[1] - r2_no_lat[1]

# Print results
cat("Partial R² for scaled_secchi: ", partial_R2_secchi, "\n")
cat("Partial R² for scaled_maxDepth: ", partial_R2_depth, "\n")
cat("Partial R² for lat: ", partial_R2_lat, "\n")

# Add predicted values to df
predictedVals <- predict(bestmodel2, type = "response", re.form=NA) #for discussion (should random effects be used or not)
mainDF_BsM$predictedDOC <- predictedVals
cor(mainDF_BsM$DOC, mainDF_BsM$predictedDOC)  # obs-fitted have 0.74 correlation - good
cor.test(mainDF_BsM$DOC, mainDF_BsM$predictedDOC)

## Relationship between observed and predicted DOC for BSM data
pObsFit <- ggplot(mainDF_BsM) +
  geom_point(aes(x=DOC,y=predictedDOC, colour = lat, size=maxDepth), alpha=0.7) +    # exp y-axis to get back on raw scale
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  scale_colour_viridis_c(option="inferno", direction = -1) +
  ylab("Model Predicted DOC (mg/L)") +
  xlab("Measured DOC (mg/L)") +
  labs(colour="Latitude", size="Max. depth (m)")+
  theme_minimal()
pObsFit


# 3) Estimate DOC in ARU data #### 
mainDF$rowID = 1:nrow(mainDF)
newdata = mainDF %>% filter(samplingProgram == "ARU")
newdata$ARU_DOC_predicted = predict(bestmodel2, newdata= newdata, type = "response", re.form = NA) #prediction only on fixed effects
#newdata$ARU_DOC_predicted_wRE = predict(bestmodel2, newdata= newdata, type = "response") #

updatedDOC = NA
#updatedDOC_RE = NA
for(i in 1:nrow(mainDF)){
  if(mainDF$samplingProgram[i] == "ARU"){ 
    updatedDOC[i] = newdata$ARU_DOC_predicted[which(newdata$rowID == i)]
 #   updatedDOC_RE[i] = newdata$ARU_DOC_predicted_wRE[which(newdata$rowID == i)]
    
  } else {
    updatedDOC[i] = mainDF$DOC[i]
   # updatedDOC_RE[i] = mainDF$DOC[i]
  }
}

mainDF$updatedDOC = updatedDOC
#mainDF$updatedDOC_RE = updatedDOC_RE


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
                                      family = Gamma(link = "log"))          
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
exp(sumDOC$coefficients$cond[, "Estimate"]) #DOC increases by 0.0008% per year
visreg(bestmodel3, "scaled_yearSample", scale="response")
r.squaredGLMM(bestmodel3)          # pretty much all lake effect

mainDF$DOC_fit = predict( bestmodel3, type = "response", re.form = NA)
mainDF$DOC_fitWRandom = predict( bestmodel3, type = "response", re.form = NULL)  # To plot just predicted lines
ggplot(mainDF)+
  geom_point(aes(x=yearSample, y = updatedDOC), alpha=0.7, size=2.5)+
  geom_smooth(aes(x=yearSample, y = DOC_fit), 
              method="glm",
              formula = y~x,
              method.args = list(family = Gamma(link = 'log')))+
  scale_x_continuous(name="Year")+
  scale_y_continuous(name = "DOC (mg/L)")+
  theme_minimal()

ggplot(mainDF, aes(x = yearSample, y = DOC_fitWRandom, colour=lat, group = commonID)) +
  geom_line(alpha=0.3) +
  theme_minimal() +
  labs(title = "Fitted Relationships Including Random Effects",
       x = "Year",
       y = "DOC (mg/L)")



# 5) Using BsM data (the best we got), what predicts DOC across Ontario? Calculate Lake Trophic Status using TDP ####
mainDF_onlyBsM <- mainDF %>% 
  filter(samplingProgram == "BsM") %>% 
  mutate(lakeTrophicStatus = case_when(TDP < 5 ~ "oligotrophic",
                                       TDP >= 5 & TDP < 10 ~ "mesotrophic",
                                       TDP >= 10 ~ "eutrophic"),
         scaledTDP = as.vector(scale(TDP, scale = FALSE)))

mainDF_onlyBsM$lakeTrophicStatus <- factor(mainDF_onlyBsM$lakeTrophicStatus, 
                                           levels = c("oligotrophic","mesotrophic","eutrophic"))  # just changing order for plotting


# Model it all: full, biologically informed IC approach
docSpatial_catTrophicStatus <- glmmTMB(DOC ~ scaled_yearSample+
                                         scaled_maxDepth+
                                         lakeTrophicStatus+
                                         lat +
                                         scaled_yearSample:lakeTrophicStatus+
                                         scaled_yearSample:lat+
                                         (1|commonID), 
                                       data=mainDF_onlyBsM, family=Gamma(link="log"), na.action = "na.fail")

#dredgeresults = dredge(docSpatial_catTrophicStatus)

bestmodel4=docSpatial_catTrophicStatus

# Investigate the results
summary(bestmodel4)
DHARMa::simulateResiduals(bestmodel4, plot=T)
hist(residuals(bestmodel4), breaks = 30, main = "Histogram of Residuals")
r.squaredGLMM(bestmodel4)         

#plot some stuff

#Reminder: consider taking the depth variable out of the pSpatialDOC plot and 
#if you want to show a depth thing do it as a separate plot (and maybe in the supplement)

pSpatialDOC<- ggplot(mainDF_onlyBsM, aes(x=yearSample,y=DOC,colour=lat)) +
  geom_point(aes(size=maxDepth), alpha=0.7) +
  scale_colour_viridis_c(option="inferno", direction = -1) +
  theme_bw()+
  facet_wrap(~lakeTrophicStatus)
pSpatialDOC

mainDF_onlyBsM$predictedDOC = predict(bestmodel4, type = "response", re.form = NA)

ggplot(mainDF_onlyBsM)+
  geom_point(aes(y=DOC, x=maxDepth), alpha=0.7, size=2.5)+
  geom_smooth(aes(x=maxDepth, y = predictedDOC), 
              method="glm",
              formula = y~x,
              method.args = list(family = Gamma(link = 'log')))+
  scale_x_continuous(name="Maximum depth (m)")+
  scale_y_continuous(name = "DOC (mg/L)")+
  theme_minimal()

ggplot(mainDF_onlyBsM)+
  geom_point(aes(y=DOC, x=yearSample, colour=lat), alpha=0.7, size=2.5)+
  geom_smooth(aes(x=yearSample, y = predictedDOC), 
              method="glm",
              formula = y~x,
              method.args = list(family = Gamma(link = 'log')))+
  scale_colour_viridis_c(option="inferno", direction = -1) +
  facet_wrap(~lakeTrophicStatus, labeller = labeller(lakeTrophicStatus = function(x) str_to_title(x)))+
  scale_x_continuous(name="Year")+
  scale_y_continuous(name = "DOC (mg/L)")+
  labs(colour="Latitude")+
  theme_minimal()
