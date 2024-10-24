
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
library(patchwork)
source("scripts/theme_DOC.R") #load ggplot theme


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
car::Anova(bestmodel,type = "III")

#nicer plot
mainDF$secchi_fit = predict(bestmodel, type = "response", re.form = NA) # Fixed effect fit
mainDF$secchi_fitWRand = predict(bestmodel, type = "response", re.form = NULL) # Random effect fit
secLinearPlot <- ggplot(mainDF)+
  geom_point(aes(x=yearSample, y = secchiDepth, fill=secchiDepth), alpha=0.7, size=2.5, shape=21)+
  geom_smooth(aes(x=yearSample, y = secchi_fit), 
              method="glm",
              formula = y~x,
              method.args = list(family = Gamma(link = 'log')))+
  scale_x_continuous(name="Year")+
  scale_y_continuous(name = "Secchi depth (m)")+
  scale_fill_gradientn(colours = c("#3D2B1F","#704214", "#D2B48C", "#F5DEB3")) +
  theme_DOC()+
  theme(legend.position = "none")
secLinearPlot

# 2) Can we predict DOC w. lake variables? Secchi, depth, lat? ####
mainDF_BsM <- filter(mainDF, samplingProgram %in% "BsM")

# Fit the initial model - create full models with log-normal and gamma, look at AIC
clarDOCmod_gamma <- glmmTMB(DOC ~ scaled_secchi * scaled_maxDepth + lat + (1|commonID), 
                            data = mainDF_BsM, family = Gamma(link = "log"), na.action = "na.fail")   # added lat in prediction because predictions w. only secchi and maxDepth under-predicted high DOC lakes. Only included as additive effect as there's a S/N gradient in DOC across the province
(clarDOCmod_gamma_sel <- dredge(clarDOCmod_gamma)) #full model best
# clarDOCmod_LMM<- glmmTMB(DOC ~ scaled_secchi * scaled_maxDepth + lat + (1|commonID),
#                             data = mainDF_BsM,  na.action = "na.fail")
# (clarDOCmod_LMM_sel <- dredge(clarDOCmod_LMM))
# clarDOCmod_logLMM<- glmmTMB(DOC ~ scaled_secchi * scaled_maxDepth + lat + (1|commonID),
#                       data = mainDF_BsM, family = lognormal(link="log"),  na.action = "na.fail")
# (clarDOCmod_logLMM_sel <- dredge(clarDOCmod_logLMM))
# AICc lowest with gamma. Do that.

bestmodel2 <- glmmTMB(DOC ~ scaled_secchi * scaled_maxDepth + lat + (1|commonID),
                              data = mainDF_BsM, family = Gamma(link = "log"))

# Diagnostics
DHARMa::simulateResiduals(bestmodel2, plot=T)
hist(residuals(bestmodel2), breaks = 30, main = "Histogram of Residuals")
# Diagnostics look fine

# Summary
summary(bestmodel2)
r2_secchiMod <- r.squaredGLMM(bestmodel2)        # R2 = 0.49 fixed effects, R2 = 0.92 w. random effects 
car::Anova(bestmodel2,type = "III")

# Partial variance explained 
# Fit reduced models
model_no_secchi <- update(bestmodel2, . ~ . - scaled_secchi)
model_no_depth <- update(bestmodel2, . ~ . - scaled_maxDepth)
model_no_lat <- update(bestmodel2, . ~ . - lat)
model_no_int <- update(bestmodel2, . ~ . - scaled_secchi:scaled_maxDepth)

# Calculate R² for reduced models
r2_no_secchi <- r.squaredGLMM(model_no_secchi)
r2_no_depth <- r.squaredGLMM(model_no_depth)
r2_no_lat <- r.squaredGLMM(model_no_lat)
r2_no_int <- r.squaredGLMM(model_no_int)

# Calculate partial R²
partial_R2_secchi <- r2_secchiMod[1] - r2_no_secchi[1]
partial_R2_depth <- r2_secchiMod[1] - r2_no_depth[1]
partial_R2_lat <- r2_secchiMod[1] - r2_no_lat[1]
partial_R2_int <- r2_secchiMod[1] - r2_no_int[1]


# Print results
cat("Partial R² for scaled_secchi: ", partial_R2_secchi, "\n")
cat("Partial R² for scaled_maxDepth: ", partial_R2_depth, "\n")
cat("Partial R² for lat: ", partial_R2_lat, "\n")
cat("Partial R² for int: ", partial_R2_int, "\n")

# Add predicted values to df
predictedVals <- predict(bestmodel2, type = "response", re.form=NA) #for discussion (should random effects be used or not)
mainDF_BsM$predictedDOC <- predictedVals
cor(mainDF_BsM$DOC, mainDF_BsM$predictedDOC)  # obs-fitted have 0.75 correlation - good
cor.test(mainDF_BsM$DOC, mainDF_BsM$predictedDOC)

## Relationship between observed and predicted DOC for BSM data
pObsFit <- ggplot(mainDF_BsM) +
  geom_point(aes(x=DOC,y=predictedDOC, colour = lat, size=maxDepth), alpha=0.7) +    # exp y-axis to get back on raw scale
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  scale_colour_viridis_c(option="inferno", direction = -1) +
  ylab("Model Predicted DOC (mg/L)") +
  xlab("Measured DOC (mg/L)") +
  labs(colour="Latitude", size="Max. depth (m)")+
  theme_minimal() +
  theme_DOC()
pObsFit

# ggsave(filename = ("plotsTables/plots/obsFittedDOC_20240713.png"),device = "png", plot = pObsFit, 
#        width = 7, height = 6, units = "in", dpi = 450)

# 2b) If we do the same with secchi, do predictors have the same sign ####
# i.e., do we see water clarity changing in a complimentary way to DOC?
bestmodel2b <- glmmTMB(secchiDepth ~ scaled_maxDepth + lat + (1|commonID),
                      data = mainDF_BsM, family = Gamma(link = "log"))

# Diagnostics
DHARMa::simulateResiduals(bestmodel2b, plot=T)
hist(residuals(bestmodel2b), breaks = 30, main = "Histogram of Residuals")
# Diagnostics look fine

# Summary
summary(bestmodel2b)             # Beauty, response moving in same direction as DOC, and coefficients (i.e., slopes) are similar
r2_secchiMod <- r.squaredGLMM(bestmodel2b)        # R2 = 0.30 fixed effects, R2 = 0.77 w. random effects 
car::Anova(bestmodel2,type = "III")

# Partial variance explained 
# Fit reduced models
model_no_depth <- update(bestmodel2b, . ~ . - scaled_maxDepth)
model_no_lat <- update(bestmodel2b, . ~ . - lat)

# Calculate R² for reduced models
r2_no_depth <- r.squaredGLMM(model_no_depth)
r2_no_lat <- r.squaredGLMM(model_no_lat)

# Calculate partial R²
partial_R2_depth <- r2_secchiMod[1] - r2_no_depth[1]
partial_R2_lat <- r2_secchiMod[1] - r2_no_lat[1]


# Print results
cat("Partial R² for scaled_maxDepth: ", partial_R2_depth, "\n")
cat("Partial R² for lat: ", partial_R2_lat, "\n")

# Add predicted values to df
predictedVals <- predict(bestmodel2b, type = "response", re.form=NA) #for discussion (should random effects be used or not)
mainDF_BsM$predictedSecchi <- predictedVals
cor(mainDF_BsM$secchiDepth, mainDF_BsM$predictedSecchi)  # obs-fitted have 0.48 correlation - good
cor.test(mainDF_BsM$secchiDepth, mainDF_BsM$predictedSecchi)


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
# Using ALL data
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
pDOCLinear <- ggplot(mainDF)+
  geom_point(aes(x=yearSample, y = updatedDOC, fill=updatedDOC), alpha=0.7, size=2.5, shape = 21) +
  geom_smooth(aes(x=yearSample, y = DOC_fit), 
              method="glm",
              formula = y~x,
              method.args = list(family = Gamma(link = 'log')))+
  scale_x_continuous(name="Year")+
  scale_y_continuous(name = "DOC (mg/L)")+
  scale_fill_gradientn(colours = c("#F5DEB3","#D2B48C","#704214","#3D2B1F")) +
  theme_DOC()+
  theme(legend.position="none")
pDOCLinear

ggplot(mainDF, aes(x = yearSample, y = DOC_fitWRandom, colour=lat, group = commonID)) +
  geom_line(alpha=0.3) +
  theme_minimal() +
  labs(title = "Fitted Relationships Including Random Effects",
       x = "Year",
       y = "DOC (mg/L)")

# 4a) Are there differences in Secchi over time, when environmental variables are also considered ####
# i.e., same effect direction and magnitidue between BsM data set AND BsM + AHI datasets?
mainDF_secMod <- mainDF %>% 
  group_by(commonID) %>%
  mutate(lakeTrophicStatus = ifelse(is.na(lakeTrophicStatus), 
                                    unique(lakeTrophicStatus[!is.na(lakeTrophicStatus)]), 
                                    lakeTrophicStatus)) %>% 
  filter(!is.na(scaled_yearSample) & 
         !is.na(lakeTrophicStatus) &
         !is.na(lat) &
         !is.na(secchiDepth))

mainDF_secMod$lakeTrophicStatus <- factor(mainDF_secMod$lakeTrophicStatus, 
                                           levels = c("oligotrophic","mesotrophic","eutrophic"))  # just changing order for plotting

# Use model formulation of final DOC model, sub in secchi depths and the whole dataset
secchiSpatial_catTrophicStatus_allData <- glmmTMB(secchiDepth ~ scaled_yearSample+
                                            scaled_maxDepth+
                                            lakeTrophicStatus+
                                            lat +
                                            scaled_yearSample:lakeTrophicStatus+
                                            scaled_yearSample:lat+
                                            (1|commonID), 
                                          data=mainDF_secMod, family=Gamma(link="log"), na.action = "na.fail")


summary(secchiSpatial_catTrophicStatus_allData)

#plot some stuff

mainDF_secMod$predictedSecchi = predict(secchiSpatial_catTrophicStatus_allData, type = "response", re.form = NA)
mainDF_secMod$predictedSecchi_wRandom = predict(secchiSpatial_catTrophicStatus_allData, type = "response",re.form = NULL)

pMaxDepth_secAll <- ggplot(mainDF_secMod)+
  geom_point(aes(y=secchiDepth, x=maxDepth, fill = secchiDepth), alpha=0.7, size=2.5, shape = 21)+
  geom_smooth(aes(x=maxDepth, y = predictedSecchi), 
              method="glm",
              formula = y~x,
              method.args = list(family = Gamma(link = 'log')))+
  scale_x_continuous(name="Maximum depth (m)")+
  scale_y_continuous(name = "Secchi depth (m)")+
  scale_fill_gradientn(colours = c("#3D2B1F","#704214", "#D2B48C", "#F5DEB3")) +
  theme_DOC()+
  theme(legend.position = "none")
pMaxDepth_secAll

pSpatialSecAll <- ggplot(data=arrange(mainDF_secMod,yearSample))+
  geom_point(aes(y=secchiDepth, x=yearSample, colour=lat), alpha=0.7, size=2.5)+
  # geom_line(aes(y=predictedSecchi, x=as.numeric(yearSample), 
  #              group = interaction(lakeTrophicStatus,lat)))+
  stat_smooth(aes(x=yearSample, y = predictedSecchi),
              method="glm",
              formula = y~x,
              method.args = list(family = Gamma(link = 'log')))+
  scale_colour_viridis_c(option="inferno", direction = -1) +
  facet_wrap(~lakeTrophicStatus, labeller = labeller(lakeTrophicStatus = function(x) str_to_title(x)))+
  scale_x_continuous(name="Year")+
  #scale_y_reverse() +
  #scale_y_reverse(transform = scales::transform_reverse())+
  scale_y_continuous(name = "Secchi depth (m)") +  
  labs(colour="Latitude")+
  theme_DOC() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
pSpatialSecAll

spatioTemp_panelPlotAll <- pSpatialSecAll + pMaxDepth_secAll + plot_layout(ncol = 1)

# Save as png

ggsave(filename = ("plotsTables/plots/spatioTemp_panelPlotAll_20240830.png"),device = "png", plot = spatioTemp_panelPlotAll,
        width = 8, height = 6,bg="white")

# 5) Using BsM data (the best we got), what predicts DOC across Ontario? Calculate Lake Trophic Status using TDP ####
mainDF_onlyBsM <- mainDF %>% 
  filter(samplingProgram == "BsM") %>% 
  mutate(lakeTrophicStatus = case_when(TDP < 5 ~ "oligotrophic",
                                       TDP >= 5 & TDP < 10 ~ "mesotrophic",
                                       TDP >= 10 ~ "eutrophic"),
         scaledTDP = as.vector(scale(TDP, scale = FALSE)))

mainDF_onlyBsM$lakeTrophicStatus <- factor(mainDF_onlyBsM$lakeTrophicStatus, 
                                           levels = c("oligotrophic","mesotrophic","eutrophic"))  # just changing order for plotting

## Let's look just at BsM and Year, AND THEN, do full spatial

                   
# Using IC approach, assemble different models to look for changes in DOC over time
M1a <- glmmTMB(updatedDOC ~ scaled_yearSample + (1|commonID),
  data = mainDF_onlyBsM,
  family = Gamma(link = "log"))
M2a <- glmmTMB(updatedDOC ~ scaled_yearSample + (scaled_yearSample|commonID),
  data = mainDF_onlyBsM,
  family = Gamma(link = "log"))          
M3a = update(M1a, family = lognormal(link="log"))
M4a = update(M2a, family = lognormal(link="log"))
AIC(M1a,M2a,M3a,M4a)
# Gamma models are best; random intercept marginally better than random slopes
                   
# Parameterize best model for DOC over time
bestmodel3a = M1a
                   
# Diagnostics
DHARMa::simulateResiduals(bestmodel3a, plot=T)
hist(residuals(bestmodel3a), breaks = 30, main = "Histogram of Residuals")
# Diagnostics look fine
                   
# Model summaries
(sumDOCa <- summary(bestmodel3a)) 
exp(sumDOCa$coefficients$cond[, "Estimate"]) #DOC increases by 0.0008% per year
visreg(bestmodel3a, "scaled_yearSample", scale="response")
r.squaredGLMM(bestmodel3a)          # pretty much all lake effect
                   
mainDF_onlyBsM$DOC_fit = predict( bestmodel3a, type = "response", re.form = NA)
mainDF_onlyBsM$DOC_fitWRandom = predict( bestmodel3a, type = "response", re.form = NULL)  # To plot just predicted lines

                  
# Using only BsM data plot DOC ~ time
pDOCTime_bsm <- ggplot(data=mainDF_onlyBsM,aes(x=yearSample,y=updatedDOC)) +
  geom_point(aes(x=yearSample, y = updatedDOC, fill=updatedDOC), alpha=0.7, size=2.5, shape = 21) +
  #geom_smooth(method = "lm") +
  scale_fill_gradientn(colours = c("#F5DEB3","#D2B48C","#704214","#3D2B1F")) +
  labs(y="DOC (mg/L)",x="Year") +
  theme_DOC()+
  theme(legend.position = "none")
pDOCTime_bsm

# ggsave(filename = ("plotsTables/plots/DOCTime_bsm_20240829.svg"),device = "svg", plot = pDOCTime_bsm,
#        width = 8, height = 6)

# Using only BsM plot Secchi ~ time
pSecTime_bsm <- ggplot(mainDF_onlyBsM)+
  geom_point(aes(x=yearSample, y = secchi_fitWRand, fill=secchi_fitWRand), alpha=0.7, size=2.5, shape = 21) +
  geom_smooth(aes(x=yearSample, y = secchi_fitWRand),
  method="glm",
  formula = y~x,
  method.args = list(family = Gamma(link = 'log')))+
  scale_x_continuous(name="Year")+
  scale_fill_gradientn(colours = c("#3D2B1F","#704214","#D2B48C","#F5DEB3")) +
  scale_y_continuous(name = "Secchi depth (m)")+
  theme_DOC() +
  theme(legend.position = "none")
pSecTime_bsm

 # ggsave(filename = ("plotsTables/plots/SecTime_bsm_20240829.svg"),device = "svg", plot = pSecTime_bsm,
 #        width = 8, height = 6)
                   
# Model it all: full, biologically informed IC approach
docSpatial_catTrophicStatus <- glmmTMB(DOC ~ scaled_yearSample+
                                         scaled_maxDepth+
                                         lakeTrophicStatus+
                                         lat +
                                         scaled_yearSample:lakeTrophicStatus+
                                         scaled_yearSample:lat+
                                         (1|commonID), 
                                       data=mainDF_onlyBsM, family=Gamma(link="log"), na.action = "na.fail")
docSpatial_catTrophicStatus_SI <- glmmTMB(DOC ~ scaled_yearSample+
                                         scaled_maxDepth+
                                         lakeTrophicStatus+
                                         lat +
                                         scaled_yearSample:lakeTrophicStatus+
                                         scaled_yearSample:lat+
                                         (scaled_yearSample|commonID), 
                                       data=mainDF_onlyBsM, family=Gamma(link="log"), na.action = "na.fail")
docSpatial_catTrophicStatus_logN_I <- glmmTMB(DOC ~ scaled_yearSample+
                                         scaled_maxDepth+
                                         lakeTrophicStatus+
                                         lat +
                                         scaled_yearSample:lakeTrophicStatus+
                                         scaled_yearSample:lat+
                                         (1|commonID), 
                                       data=mainDF_onlyBsM, family=lognormal(link="log"), na.action = "na.fail")
docSpatial_catTrophicStatus_SI_logN_SI <- glmmTMB(DOC ~ scaled_yearSample+
                                            scaled_maxDepth+
                                            lakeTrophicStatus+
                                            lat +
                                            scaled_yearSample:lakeTrophicStatus+
                                            scaled_yearSample:lat+
                                            (scaled_yearSample|commonID), 
                                          data=mainDF_onlyBsM, family=lognormal(link="log"), na.action = "na.fail")
AIC(docSpatial_catTrophicStatus,docSpatial_catTrophicStatus_SI,
    docSpatial_catTrophicStatus_logN_I,docSpatial_catTrophicStatus_SI_logN_SI)

#dredgeresults = dredge(docSpatial_catTrophicStatus)

bestmodel4=docSpatial_catTrophicStatus

# Investigate the results
summary(bestmodel4)
Anova(bestmodel4)
DHARMa::simulateResiduals(bestmodel4, plot=T)
hist(residuals(bestmodel4), breaks = 30, main = "Histogram of Residuals")
r.squaredGLMM(bestmodel4)         

#plot some stuff

#Reminder: consider taking the depth variable out of the pSpatialDOC plot and 
#if you want to show a depth thing do it as a separate plot (and maybe in the supplement)

mainDF_onlyBsM$predictedDOC = predict(bestmodel4, type = "response", re.form = NA)
mainDF_onlyBsM$predictedDOC_wRandom = predict(bestmodel4, type = "response",re.form = NULL)

pMaxDepth <- ggplot(mainDF_onlyBsM)+
  geom_point(aes(y=DOC, x=maxDepth, fill=DOC), alpha=0.7, size=2.5, shape=21)+
  geom_smooth(aes(x=maxDepth, y = predictedDOC), 
              method="glm",
              formula = y~x,
              method.args = list(family = Gamma(link = 'log')))+
  scale_x_continuous(name="Maximum depth (m)")+
  scale_y_continuous(name = "DOC (mg/L)")+
  scale_fill_gradientn(colours = c("#F5DEB3","#D2B48C","#704214","#3D2B1F")) +
  theme_DOC()+
  theme(legend.position="none")
pMaxDepth

pSpatialDOC <- ggplot(mainDF_onlyBsM)+
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
  theme_DOC() 
pSpatialDOC

# 5b) Using all only BsM data, do we see consistent trends between secchi depths across the province? ####
# i.e., do we see the same statistical signs and relative magnitudes between DOC and secchi?
secchiTimeMod <- glmmTMB(secchiDepth ~ scaled_yearSample+ (1|commonID), 
                         data=mainDF_onlyBsM, family=Gamma(link="log"), na.action = "na.fail")
secchiTimeModb <- glmmTMB(secchiDepth ~ scaled_yearSample + (scaled_yearSample|commonID),
               data = mainDF_onlyBsM,
               family = Gamma(link = "log"))          
secchiTimeModc = update(secchiTimeMod, family = lognormal(link="log"))
secchiTimeModd = update(secchiTimeModb, family = lognormal(link="log"))
AIC(secchiTimeMod,secchiTimeModb,secchiTimeModc,secchiTimeModd)

summary(secchiTimeModb)

# Use model formulation of final DOC model, sub in secchi depths and the whole dataset
secchiSpatial_catTrophicStatus <- glmmTMB(secchiDepth ~ scaled_yearSample+
                                         scaled_maxDepth+
                                         lakeTrophicStatus+
                                         lat +
                                         scaled_yearSample:lakeTrophicStatus+
                                         scaled_yearSample:lat+
                                         (1|commonID), 
                                       data=mainDF_onlyBsM, family=Gamma(link="log"), na.action = "na.fail")

#Investigate the results
summary(secchiSpatial_catTrophicStatus)
DHARMa::simulateResiduals(secchiSpatial_catTrophicStatus, plot=T)
hist(residuals(secchiSpatial_catTrophicStatus), breaks = 30, main = "Histogram of Residuals")
r.squaredGLMM(secchiSpatial_catTrophicStatus)         

#plot some stuff

#Reminder: consider taking the depth variable out of the pSpatialDOC plot and 
#if you want to show a depth thing do it as a separate plot (and maybe in the supplement)

mainDF_onlyBsM$predictedSecchi = predict(secchiSpatial_catTrophicStatus, type = "response", re.form = NA)
mainDF_onlyBsM$predictedSecchi_wRandom = predict(secchiSpatial_catTrophicStatus, type = "response",re.form = NULL)

pMaxDepth_sec <- ggplot(mainDF_onlyBsM)+
  geom_point(aes(y=secchiDepth, x=maxDepth, fill = secchiDepth), alpha=0.7, size=2.5, shape = 21)+
  geom_smooth(aes(x=maxDepth, y = predictedSecchi), 
              method="glm",
              formula = y~x,
              method.args = list(family = Gamma(link = 'log')))+
  scale_x_continuous(name="Maximum depth (m)")+
  scale_y_continuous(name = "Secchi depth (m)")+
  scale_fill_gradientn(colours = c("#3D2B1F","#704214", "#D2B48C", "#F5DEB3")) +
  theme_DOC()+
  theme(legend.position = "none")
pMaxDepth_sec

pSpatialSec <- ggplot(data=arrange(mainDF_onlyBsM,yearSample))+
  geom_point(aes(y=secchiDepth, x=yearSample, colour=lat), alpha=0.7, size=2.5)+
  # geom_line(aes(y=predictedSecchi, x=as.numeric(yearSample), 
  #              group = interaction(lakeTrophicStatus,lat)))+
   stat_smooth(aes(x=yearSample, y = predictedSecchi),
            method="glm",
            formula = y~x,
            method.args = list(family = Gamma(link = 'log')))+
  scale_colour_viridis_c(option="inferno", direction = -1) +
  facet_wrap(~lakeTrophicStatus, labeller = labeller(lakeTrophicStatus = function(x) str_to_title(x)))+
  scale_x_continuous(name="Year")+
  #scale_y_reverse() +
  #scale_y_reverse(transform = scales::transform_reverse())+
  scale_y_continuous(name = "Secchi depth (m)") +  
  labs(colour="Latitude")+
  theme_DOC() 
pSpatialSec


# 6) What are differences in DOC for each lake across time periods ####
# Summarise each lakes fitted relationship min/max values

# Ontario shapefile: 
ontario_shapefile <- st_read("data/Province/Province.shp")
ontario_shapefile <-st_transform(ontario_shapefile, crs = 4326)

## Change in DOC and Secchi with only BsM
docBsMTimeDat <- mainDF_onlyBsM %>% 
  group_by(commonID) %>% 
  summarise(initDOC = predictedDOC_wRandom[which.min(yearSample)],
            finalDOC = predictedDOC_wRandom[which.max(yearSample)],
            #initDOC = updatedDOC[which.min(yearSample)],
            #finalDOC = updatedDOC[which.max(yearSample)],
            diffDOC = finalDOC-initDOC,
            yearDiffDOC = max(yearSample) - min(yearSample),
            initSecchi = predictedSecchi_wRandom[which.min(yearSample)],
            finalSecchi = predictedSecchi_wRandom[which.max(yearSample)],
            #initSecchi = predictedSecchi[which.min(yearSample)],
            #finalSecchi = predictedSecchi[which.max(yearSample)],
            diffSecchi = finalSecchi-initSecchi,
            sampleSize = n()) %>% 
  left_join(select(statDat_wide,commonID,lat,long)) %>% 
  distinct() %>% 
  filter(!is.na(lat)& yearDiffDOC >0) %>% 
  mutate(dataset = "all")
docBsMTimeDat

## Change in DOC and Secchi with AHI and BsM
docTimeDat <- mainDF %>% 
  group_by(commonID) %>% 
  summarise(initDOC = DOC_fitWRandom[which.min(yearSample)],
            finalDOC = DOC_fitWRandom[which.max(yearSample)],
            diffDOC = finalDOC-initDOC,
            yearDiffDOC = max(yearSample) - min(yearSample),
            initSecchi = secchi_fitWRand[which.min(yearSample)],
            finalSecchi = secchi_fitWRand[which.max(yearSample)],
            diffSecchi = finalSecchi-initSecchi,
            sampleSize = n()) %>% 
  left_join(select(statDat_wide,commonID,lat,long)) %>% 
  distinct() %>% 
  filter(!is.na(lat)) %>% 
  mutate(dataset = "all")
docTimeDat

## Map of DOC, only BsM
pDOCtimeDiff_BsM <- ggplot() +
  geom_sf(data = ontario_shapefile, fill = NA, color = "black") +
  geom_point(data=docBsMTimeDat, aes(y=lat,x=long,fill=diffDOC,size=yearDiffDOC), 
             shape = 21, stroke = 0.5, colour = "black") +
  scale_fill_gradient2(
    low = "#aa95cd",           
    mid = "white",         
    high = "brown",     
    midpoint = 0)+  
  ylim(c(42.5,55)) +
  labs(y="Latitude",x="Longitude",fill="Δ DOC (mg/L)",size="Sample period (years)") +
  theme_DOC()  +
  theme(legend.position = "top")
pDOCtimeDiff_BsM

# ggsave(filename = ("plotsTables/plots/doc_bsmMap_20240829.svg"),device = "svg", plot = pDOCtimeDiff_BsM,
#        width = 8, height = 6)

## Map of Secchi, only BsM
pSecchitimeDiff_BsM <- ggplot() +
  geom_sf(data = ontario_shapefile, fill = NA, color = "black") +
  geom_point(data=docBsMTimeDat, aes(y=lat,x=long,fill=diffSecchi,size=yearDiffDOC), 
             shape = 21, stroke = 0.5, colour = "black") +
  scale_fill_gradient2(
    high = "#aa95cd",            
    mid = "white",         
    low = "brown",     
    midpoint = 0)+  
  ylim(c(42.5,55)) +
  labs(y="Latitude",x="Longitude",fill="Δ Secchi depth (m)",size="Sample period (years)") +
  theme_DOC() +
  theme(legend.position = "none")
pSecchitimeDiff_BsM

# ggsave(filename = ("plotsTables/plots/waterClar_bsmMap_20240829.svg"),device = "svg", plot = pSecchitimeDiff_BsM,
#        width = 8, height = 6)

## Map of DOC, including estimated DOC from AHI
pDOCtimeDiff <- ggplot() +
  geom_sf(data = ontario_shapefile, fill = NA, color = "black") +
  geom_point(data=docTimeDat, aes(y=lat,x=long,fill=diffDOC,size=yearDiffDOC), 
             shape = 21, stroke = 0.5, colour = "black") +
  scale_fill_gradient2(
    low = "navyblue",          
    mid = "white",         
    high = "brown",     
    midpoint = 0)+  
  ylim(c(42.5,55)) +
  labs(y="Latitude",x="Longitude",fill="Δ DOC (mg/L)",size="Sample period (years)") +
  theme_DOC() +
  theme(legend.position = "none")
pDOCtimeDiff

## Map of Secchi, including AHI
pSecchitimeDiff <- ggplot() +
  geom_sf(data = ontario_shapefile, fill = NA, color = "black") +
  geom_point(data=docTimeDat, aes(y=lat,x=long,fill=diffSecchi,size=yearDiffDOC), 
             shape = 21, stroke = 0.5, colour = "black") +
  scale_fill_gradient2(
    high = "navyblue",          
    mid = "white",         
    low = "brown",     
    midpoint = 0)+  
  ylim(c(42.5,55)) +
  labs(y="Latitude",x="Longitude",fill="Δ Secchi depth (m)",size="Sample period (years)") +
  theme_DOC() +
  theme(legend.position = "none")
pSecchitimeDiff

#### Patchwork - assemble plots ####
library(svglite)
waterClar_panelPlot <- pDOCLinear+secLinearPlot+pDOCtimeDiff+pSecchitimeDiff + plot_layout(heights = c(1.5,2))
spatioTemp_panelPlot <- pSpatialDOC + pMaxDepth + plot_layout(ncol = 1)
spatioTempSecchi_panelPlot <- pSpatialSec + pMaxDepth_sec + plot_layout(ncol = 1)
# Patch

# Save as svg
 # ggsave(filename = ("plotsTables/plots/waterClar_panelPlot_updated_20240830.svg"),device = "svg", plot = waterClar_panelPlot,
 #        width = 12, height = 12)
# 
# ggsave(filename = ("plotsTables/plots/spatioTemp_panelPlot_20240807.svg"),device = "svg", plot = spatioTemp_panelPlot, 
#         width = 8, height = 6,bg="white")
# 
# ggsave(filename = ("plotsTables/plots/spatioTempSec_panelPlot_20240807.svg"),device = "svg", plot = spatioTempSecchi_panelPlot,
#        width = 8, height = 6,bg="white")


## Create Tables ####
lakeLoc <- read_csv("data/lakeLocations.csv") %>% 
  group_by(commonID) %>% 
  distinct()

#write_csv(lakeLoc, "lakeLocations_forWade_20240714.csv")

### Extra Plots ####
docTP <- ggplot(mainDF_onlyBsM) +
  geom_point(aes(y=updatedDOC,x=TDP,colour=lat))+
  #geom_smooth(aes(y=updatedDOC,x=TDP), method = "lm")+
  scale_colour_viridis_c(option="inferno", direction = -1) +
  facet_wrap(~lakeTrophicStatus, scales="free_x")+
  scale_y_continuous(name = "DOC (mg/L)") +
  scale_x_continuous(name = "Total dissolved P (ug/L)") +
  theme_DOC()
docTP

