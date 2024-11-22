## Last script update. This reflects results as they're presented in the MS
# AJR: 2024-09-01

# Load libraries
library(tidyverse)
library(here)
library(lme4)
library(lmerTest)
library(glmmTMB)
library(MuMIn)
library(visreg)
library(patchwork)
library(svglite)
library(sf)
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
# Ontario shapefile: 
ontario_shapefile <- st_read("data/Province/Province.shp")
ontario_shapefile <-st_transform(ontario_shapefile, crs = 4326)

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

#write_csv(mainDF,"data/dataForDryadRepo_landscapeDOC_AR_20241122.csv")

# Can we predict DOC w. lake variables? Secchi, depth, lat? ####
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
  ylim(c(0,21))+
  xlim(c(0,21)) +
  theme_minimal() +
  theme_DOC()
pObsFit

# ggsave(filename = ("plotsTables/plots/obsFittedDOC_20241118.png"),device = "png", plot = pObsFit,
#        width = 7, height = 6, units = "in", dpi = 450)

# If we do the same with secchi, do predictors have the same sign ####
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


# Estimate DOC in ARU data #### 
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


## Change in DOC and Secchi; ONLY BsM ####

# Create BsM only DF
mainDF_onlyBsM <- mainDF %>% 
  filter(samplingProgram == "BsM") %>% 
  mutate(lakeTrophicStatus = case_when(TDP < 5 ~ "oligotrophic",
                                       TDP >= 5 & TDP < 10 ~ "mesotrophic",
                                       TDP >= 10 ~ "eutrophic"),
         scaledTDP = as.vector(scale(TDP, scale = FALSE)))

mainDF_onlyBsM$lakeTrophicStatus <- factor(mainDF_onlyBsM$lakeTrophicStatus, 
                                           levels = c("oligotrophic","mesotrophic","eutrophic"))  # just changing order for plotting

## DOC
# Using IC approach, assemble different models to look for changes in DOC between 2008-2023
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
exp(sumDOCa$coefficients$cond[, "Estimate"]) #DOC decrease by 0.0032% per year; NS
visreg(bestmodel3a, "scaled_yearSample", scale="response")
r.squaredGLMM(bestmodel3a)          # pretty much all lake effect

mainDF_onlyBsM$DOC_fit = predict( bestmodel3a, type = "response", re.form = NA)
mainDF_onlyBsM$DOC_fitWRandom = predict( bestmodel3a, type = "response", re.form = NULL)  # To plot just predicted lines

# Plot it:
pDOCTime_bsm <- ggplot(data=mainDF_onlyBsM,aes(x=yearSample,y=updatedDOC)) +
  geom_jitter(aes(x=yearSample, y = updatedDOC, fill=updatedDOC), alpha=0.7, size=2.5, shape=21, width = 0.25) +
  #geom_smooth(method = "lm") +      # Don't plot line as relationship not significant
  scale_fill_gradientn(colours = c("#F5DEB3","#D2B48C","#704214","#3D2B1F")) +
  labs(y="DOC (mg/L)",x="Year") +
  theme_DOC()+
  theme(legend.position = "none")
pDOCTime_bsm

# ggsave(filename = ("plotsTables/plots/DOCTime_bsm_20241118.svg"),device = "svg", plot = pDOCTime_bsm,
#        width = 8, height = 6)

## Secchi - do we see consistent trends with DOC?
secchiTimeMod <- glmmTMB(secchiDepth ~ scaled_yearSample+ (1|commonID), 
                         data=mainDF_onlyBsM, family=Gamma(link="log"), na.action = "na.fail")
secchiTimeModb <- glmmTMB(secchiDepth ~ scaled_yearSample + (scaled_yearSample|commonID),
                          data = mainDF_onlyBsM,
                          family = Gamma(link = "log"))          
secchiTimeModc = update(secchiTimeMod, family = lognormal(link="log"))
secchiTimeModd = update(secchiTimeModb, family = lognormal(link="log"))
AIC(secchiTimeMod,secchiTimeModb,secchiTimeModc,secchiTimeModd)

secchiBsmMod_year <- secchiTimeModb

# Diagnostics
DHARMa::simulateResiduals(secchiBsmMod_year, plot=T)
hist(residuals(secchiBsmMod_year), breaks = 30, main = "Histogram of Residuals")
# Diagnostics look fine

# Model summaries
(sumSecchiYear <- summary(secchiBsmMod_year)) 
exp(sumSecchiYear$coefficients$cond[, "Estimate"]) #Secchi decrease by 1.27% per year; significant
r.squaredGLMM(secchiBsmMod_year)          # pretty much all lake effect

mainDF_onlyBsM$Secchi_fit = predict( secchiBsmMod_year, type = "response", re.form = NA)
mainDF_onlyBsM$Secchi_fitWRandom = predict( secchiBsmMod_year, type = "response", re.form = NULL)  # To plot just predicted lines

secchiBsmYearPlot <- ggplot(mainDF_onlyBsM)+
  geom_jitter(aes(x=yearSample, y = Secchi_fitWRandom, fill=secchiDepth), alpha=0.7, size=2.5, shape=21, width = 0.25)+
  geom_smooth(aes(x=yearSample, y = Secchi_fitWRandom), 
              method="glm",
              formula = y~x,
              method.args = list(family = Gamma(link = 'log')))+
  scale_x_continuous(name="Year")+
  scale_y_continuous(name = "Secchi depth (m)")+
  scale_fill_gradientn(colours = c("#3D2B1F","#704214", "#D2B48C", "#F5DEB3")) +
  theme_DOC()+
  theme(legend.position = "none")
secchiBsmYearPlot

# ggsave(filename = ("plotsTables/plots/secchiTime_bsm_20241118.svg"),device = "svg", plot = secchiBsmYearPlot,
#        width = 8, height = 6)

# Summarize modelled differences in Secchi over time period for each trophic group
# Extract data from plot above - want the min/max of the fitted lines
plotData_secNoEnv <- ggplot_build(secchiBsmYearPlot)
# Extract the data related to geom_smooth
lineData_secNoEnv <- plotData_secNoEnv$data[[2]]  # The second layer corresponds to geom_smooth

# Group by lakeTrophicStatus and calculate the min/max of the fitted values
sumsecNoEnv_timeBsM <- lineData_secNoEnv %>%
  summarise(minSec = min(y),
            maxSec = max(y),
            diffSec = maxSec - minSec,
            yearDiff = max(x)-min(x),
            diffPerDecade = (diffSec/yearDiff)*10)


## Change in DOC and Secchi with only BsM
docBsMTimeDat <- mainDF_onlyBsM %>% 
  group_by(commonID) %>% 
  summarise(initDOC = updatedDOC[which.min(yearSample)],
            finalDOC = updatedDOC[which.max(yearSample)],
            diffDOC = finalDOC-initDOC,
            yearDiffDOC = max(yearSample) - min(yearSample),
            initSecchi = secchiDepth[which.min(yearSample)],
            finalSecchi = secchiDepth[which.max(yearSample)],
            diffSecchi = finalSecchi-initSecchi,
            sampleSize = n()) %>% 
  left_join(select(statDat_wide,commonID,lat,long)) %>% 
  distinct() %>% 
  filter(!is.na(lat)& yearDiffDOC >0) %>% 
  mutate(dataset = "all")
docBsMTimeDat

## Make maps of deltaDOC and deltaSecchi
## Map of DOC, only BsM
pDOCtimeDiff_BsM <- ggplot() +
  geom_sf(data = ontario_shapefile, fill = NA, color = "black") +
  geom_point(data=docBsMTimeDat, aes(y=lat,x=long,fill=diffDOC,size=yearDiffDOC), 
             shape = 21, stroke = 0.5, colour = "black", alpha=0.7) +
  scale_fill_gradient2(
    low = "#aa95cd",           
    mid = "white",         
    high = "brown",     
    midpoint = 0)+  
  ylim(c(42.5,55)) +
  labs(y="Latitude",x="Longitude",fill="Δ DOC (mg/L)",size="Sample period (years)") +
  theme_DOC()  +
  theme(legend.position = "none")
pDOCtimeDiff_BsM

pDOCTimeDiff_BsM_legTop = pDOCtimeDiff_BsM +
  theme(legend.position = "top")  #This is just to add legend to plot (custom)
pDOCTimeDiff_BsM_legRight = pDOCtimeDiff_BsM +
  theme(legend.position = "right")  #This is just to add legend to plot (custom)

# ggsave(filename = ("plotsTables/plots/doc_bsmMap_20240829.svg"),device = "svg", plot = pDOCtimeDiff_BsM,
#        width = 8, height = 6)
# ggsave(filename = ("plotsTables/plots/doc_bsmMap_legTop_20240829.svg"),device = "svg", plot = pDOCTimeDiff_BsM_legTop,
#        width = 8, height = 6)
# ggsave(filename = ("plotsTables/plots/doc_bsmMap_legRight_20240829.svg"),device = "svg", plot = pDOCTimeDiff_BsM_legRight,
#        width = 8, height = 6)

## Map of Secchi, only BsM
pSecchitimeDiff_BsM <- ggplot() +
  geom_sf(data = ontario_shapefile, fill = NA, color = "black") +
  geom_point(data=docBsMTimeDat, aes(y=lat,x=long,fill=diffSecchi,size=yearDiffDOC), 
             shape = 21, stroke = 0.5, colour = "black", alpha=0.7) +
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

pSecchiTimeDiff_BsM_legTop = pSecchitimeDiff_BsM +
  theme(legend.position = "top")  #This is just to add legend to plot (custom)
# ggsave(filename = ("plotsTables/plots/secchi_bsmMap_legTop_20240901.svg"),device = "svg", plot = pSecchiTimeDiff_BsM_legTop,
#        width = 8, height = 6)

waterClarBsM_panelPlot <- pDOCTime_bsm+secchiBsmYearPlot+pDOCtimeDiff_BsM+pSecchitimeDiff_BsM + plot_layout(heights = c(1.5,2))

# Save as svg
# ggsave(filename = ("plotsTables/plots/waterClarBsM_panelPlot_20241118.svg"),device = "svg", plot = waterClarBsM_panelPlot,
#        width = 12, height = 12)

## Change in DOC and Secchi with environmental variables; ONLY BsM ####
# Model it all: full, biologically informed IC approach
docSpatial_catTrophicStatus <- glmmTMB(updatedDOC ~ scaled_yearSample+
                                         scaled_maxDepth+
                                         lakeTrophicStatus+
                                         lat +
                                         scaled_yearSample:lakeTrophicStatus+
                                         scaled_yearSample:lat+
                                         (1|commonID), 
                                       data=mainDF_onlyBsM, family=Gamma(link="log"), na.action = "na.fail")
docSpatial_catTrophicStatus_SI <- glmmTMB(updatedDOC ~ scaled_yearSample+
                                            scaled_maxDepth+
                                            lakeTrophicStatus+
                                            lat +
                                            scaled_yearSample:lakeTrophicStatus+
                                            scaled_yearSample:lat+
                                            (scaled_yearSample|commonID), 
                                          data=mainDF_onlyBsM, family=Gamma(link="log"), na.action = "na.fail")
docSpatial_catTrophicStatus_logN_I <- glmmTMB(updatedDOC ~ scaled_yearSample+
                                                scaled_maxDepth+
                                                lakeTrophicStatus+
                                                lat +
                                                scaled_yearSample:lakeTrophicStatus+
                                                scaled_yearSample:lat+
                                                (1|commonID), 
                                              data=mainDF_onlyBsM, family=lognormal(link="log"), na.action = "na.fail")
docSpatial_catTrophicStatus_SI_logN_SI <- glmmTMB(updatedDOC ~ scaled_yearSample+
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

bestmodel4=glmmTMB(updatedDOC ~ scaled_yearSample+
                  scaled_maxDepth+
                  lakeTrophicStatus+
                  lat +
                  scaled_yearSample:lakeTrophicStatus+
                  (1|commonID), 
                  data=mainDF_onlyBsM, family=Gamma(link="log"))

# Investigate the results
summary(bestmodel4)
car::Anova(bestmodel4)
DHARMa::simulateResiduals(bestmodel4, plot=T)
hist(residuals(bestmodel4), breaks = 30, main = "Histogram of Residuals")
r.squaredGLMM(bestmodel4) 

# Predict **only use this for maxDepth**; can't parse out mean lakeTrophicStatus
mainDF_onlyBsM$predictedDOC_wEnvVars = predict(bestmodel4, type = "response", re.form = NA)
mainDF_onlyBsM$predictedDOC_wRandomEnvVars = predict(bestmodel4, type = "response",re.form = NULL)

pMaxDepthBsM <- ggplot(mainDF_onlyBsM)+
  geom_point(aes(y=DOC, x=maxDepth, fill=DOC), alpha=0.5, size=2.5, shape=21)+
  geom_smooth(aes(x=maxDepth, y = predictedDOC_wEnvVars), 
              method="glm",
              formula = y~x,
              method.args = list(family = Gamma(link = 'log')))+
  scale_x_continuous(name="Maximum depth (m)")+
  scale_y_continuous(name = "DOC (mg/L)")+
  ylim(c(0,21)) +
  scale_fill_gradientn(colours = c("#F5DEB3","#D2B48C","#704214","#3D2B1F")) +
  theme_DOC()+
  theme(legend.position="none")
pMaxDepthBsM

# Plot partial regressions of interaction between year and trophic status
newdata = mainDF_onlyBsM
newdata$lat2 = newdata$lat
newdata$lat = mean(newdata$lat)
newdata$scaled_maxDepth = mean(newdata$scaled_maxDepth)

newdata$predictedDOC_wEnvVars = predict(bestmodel4,
                                        newdata = newdata,
                                        type = "response", re.form = NA)

pSpatialDOCBsM <- ggplot(newdata)+
  geom_jitter(aes(y=DOC, x=yearSample, colour=lat2), alpha=0.3, size=2, width=0.22)+
  geom_smooth(aes(x=yearSample, y = predictedDOC_wEnvVars),
              method="glm",
              formula = y~x,
              method.args = list(family = Gamma(link = 'log')))+
  scale_colour_viridis_c(option="inferno", direction = -1) +
  facet_wrap(~lakeTrophicStatus, labeller = labeller(lakeTrophicStatus = function(x) str_to_title(x)))+
  ylim(c(0,21)) +
  labs(colour="Latitude", y="DOC (mg/L)", x="Year" )+
  theme_DOC()
pSpatialDOCBsM

#Below is using 'predict' function data where all model variables are influencing fitted relationship. DOES NOT SHOW INTERACTION B/W YEAR & TROPHIC STATUS
# pSpatialDOCBsM <- ggplot(mainDF_onlyBsM)+
#   #geom_point(aes(y=DOC, x=yearSample, colour=lat, size = scaled_maxDepth), alpha=0.7)+
#   geom_jitter(aes(y=DOC, x=yearSample, colour=lat, size = scaled_maxDepth), alpha=0.7, width=0.22)+
#   geom_smooth(aes(x=yearSample, y = predictedDOC_wEnvVars),
#               method="glm",
#               formula = y~x,
#               method.args = list(family = Gamma(link = 'log')))+
#   scale_colour_viridis_c(option="inferno", direction = -1) +
#   facet_wrap(~lakeTrophicStatus, labeller = labeller(lakeTrophicStatus = function(x) str_to_title(x)))+
#   scale_x_continuous(name="Year")+
#   scale_y_continuous(name = "DOC (mg/L)")+
#   labs(colour="Latitude")+
#   theme_DOC()
# pSpatialDOCBsM


# Summarize modelled differences in DOC over time period for each trophic group
# Extract data from plot above - want the min/max of the fitted lines
plotData <- ggplot_build(pSpatialDOCBsM)
# Extract the data related to geom_smooth
lineData <- plotData$data[[2]]  # The second layer corresponds to geom_smooth

# Group by lakeTrophicStatus and calculate the min/max of the fitted values
sumDOC_timeBsM <- lineData %>%
  group_by(PANEL) %>%  # PANEL corresponds to lakeTrophicStatus facet (1, Oligo; 2, Meso; 3, Eutro)
  summarise(minDOC = min(y),
            maxDOC = max(y),
            diffDOC = maxDOC - minDOC,
            yearDiff = max(x)-min(x),
            diffPerDecade = (diffDOC/yearDiff)*10)


## Now look at Secchi
# Use model formulation of final DOC model, sub in secchi depths and the whole dataset
secchiSpatial_catTrophicStatus <- glmmTMB(secchiDepth ~ scaled_yearSample+
                                            scaled_maxDepth+
                                            lakeTrophicStatus+
                                            lat +
                                            scaled_yearSample:lakeTrophicStatus+
                                            (1|commonID), 
                                          data=mainDF_onlyBsM, family=Gamma(link="log"))

#Investigate the results
summary(secchiSpatial_catTrophicStatus)
DHARMa::simulateResiduals(secchiSpatial_catTrophicStatus, plot=T)
hist(residuals(secchiSpatial_catTrophicStatus), breaks = 30, main = "Histogram of Residuals")
r.squaredGLMM(secchiSpatial_catTrophicStatus)      

# Plot some stuff
mainDF_onlyBsM$predictedSecchi_wEnvVars = predict(secchiSpatial_catTrophicStatus, type = "response", re.form = NA)
mainDF_onlyBsM$predictedSecchi_wRandomEnvVars = predict(secchiSpatial_catTrophicStatus, type = "response",re.form = NULL)

pMaxDepth_sec_BsM <- ggplot(mainDF_onlyBsM)+
  geom_point(aes(y=secchiDepth, x=maxDepth, fill = secchiDepth), alpha=0.5, size=2.5, shape = 21)+
  geom_smooth(aes(x=maxDepth, y = predictedSecchi_wEnvVars), 
              method="glm",
              formula = y~x,
              method.args = list(family = Gamma(link = 'log')))+
  scale_x_continuous(name="Maximum depth (m)")+
  scale_y_continuous(name = "Secchi depth (m)")+
  scale_fill_gradientn(colours = c("#3D2B1F","#704214", "#D2B48C", "#F5DEB3")) +
  theme_DOC()+
  #ylim(c(0,21)) +
  theme(legend.position = "none")
pMaxDepth_sec_BsM

## Plot only partial regression containing trophic status * year interaction
newdata_sec = mainDF_onlyBsM
newdata_sec$lat2 = newdata_sec$lat
newdata_sec$lat = mean(newdata_sec$lat)
newdata_sec$scaled_maxDepth = mean(newdata_sec$scaled_maxDepth)

newdata_sec$predictedSecchi_wEnvVars = predict(secchiSpatial_catTrophicStatus,
                                        newdata = newdata_sec,
                                        type = "response", re.form = NA)

pSpatialSecBsM <- ggplot(newdata_sec)+
  geom_jitter(data = mainDF_onlyBsM, aes(y=secchiDepth, x=yearSample, colour=lat), alpha=0.3, size=2, width=0.22)+
  geom_smooth(aes(x=yearSample, y = predictedSecchi_wEnvVars),
              method="glm",
              formula = y~x,
              method.args = list(family = Gamma(link = 'log')))+
  scale_colour_viridis_c(option="inferno", direction = -1) +
  facet_wrap(~lakeTrophicStatus, labeller = labeller(lakeTrophicStatus = function(x) str_to_title(x)))+
  scale_x_continuous(name="Year")+
  scale_y_continuous(name = "Secchi depth (m)") +  
  labs(colour="Latitude")+
 # ylim(c(0,21)) +
  theme_DOC() 
pSpatialSecBsM

# Summarize modelled differences in Secchi over time period for each trophic group
# Extract data from plot above - want the min/max of the fitted lines
plotData_sec <- ggplot_build(pSpatialSecBsM)
# Extract the data related to geom_smooth
lineData_sec <- plotData_sec$data[[2]]  # The second layer corresponds to geom_smooth

# Group by lakeTrophicStatus and calculate the min/max of the fitted values
sumSec_timeBsM <- lineData_sec %>%
  group_by(PANEL) %>%  # PANEL corresponds to lakeTrophicStatus facet (1, Oligo; 2, Meso; 3, Eutro)
  summarise(minSec = min(y),
            maxSec = max(y),
            diffSec = maxSec - minSec,
            yearDiff = max(x)-min(x),
            diffPerDecade = (diffSec/yearDiff)*10)

spatioTempBsM_panelPlot <- pSpatialDOCBsM + pMaxDepthBsM + plot_layout(ncol = 1)
spatioTempSecchiBsM_panelPlot <- pSpatialSecBsM + pMaxDepth_sec_BsM + plot_layout(ncol = 1)

# ggsave(filename = ("plotsTables/plots/spatioTempBsM_panelPlot_20241118.svg"),device = "svg", plot = spatioTempBsM_panelPlot,
#         width = 8, height = 6,bg="white")
# 
# ggsave(filename = ("plotsTables/plots/spatioTempSecchiBsM_panelPlot_20241118.svg"),device = "svg", plot = spatioTempSecchiBsM_panelPlot,
#        width = 8, height = 6,bg="white")


## Change in DOC and Secchi; BsM and AHI ####
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
exp(sumDOC$coefficients$cond[, "Estimate"]) #DOC increases by 0.08% per year; 1.05 mg/L increase over 62 years (4.9% increase over time); 0.3mg/L per decade
r.squaredGLMM(bestmodel3)          # pretty much all lake effect

mainDF$DOC_fit_allData_noEnv = predict( bestmodel3, type = "response", re.form = NA)
mainDF$DOC_fitWRandom_allData_noEnv = predict( bestmodel3, type = "response", re.form = NULL)  # To plot just predicted lines

# Plot it
pDOC_allData_time <- ggplot(mainDF)+
  geom_point(aes(x=yearSample, y = updatedDOC, fill=updatedDOC), alpha=0.7, size=2.5, shape = 21) +
  geom_smooth(aes(x=yearSample, y = DOC_fit_allData_noEnv), 
              method="glm",
              formula = y~x,
              method.args = list(family = Gamma(link = 'log')))+
  scale_x_continuous(name="Year")+
  scale_y_continuous(name = "DOC (mg/L)")+
  scale_fill_gradientn(colours = c("#F5DEB3","#D2B48C","#704214","#3D2B1F")) +
  theme_DOC()+
  theme(legend.position="none")
pDOC_allData_time

# Summarize modelled differences in DOC over time period 
# Extract data from plot above - want the min/max of the fitted lines
plotData_DOC <- ggplot_build(pDOC_allData_time)
# Extract the data related to geom_smooth
lineData_DOC <- plotData_DOC$data[[2]]  # The second layer corresponds to geom_smooth

# Group by lakeTrophicStatus and calculate the min/max of the fitted values
sumDOC_allYear <- lineData_DOC %>%
  summarise(minDOC = min(y),
            maxDOC = max(y),
            diffDOC = maxDOC - minDOC,
            yearDiff = max(x)-min(x),
            diffPerDecade = (diffDOC/yearDiff)*10)

## Now Secchi
M1_sec <- glmmTMB(secchiDepth ~ scaled_yearSample + (1|commonID),
              data = mainDF,
              family = Gamma(link = "log"))
M2_sec <- glmmTMB(secchiDepth ~ scaled_yearSample + (scaled_yearSample|commonID),
              data = mainDF,
              family = Gamma(link = "log"))          
M3_sec = update(M1, family = lognormal(link="log"))
M4_sec = update(M2, family = lognormal(link="log"))
AIC(M1_sec,M2_sec,M3_sec,M4_sec)
# Gamma models are best. 

# Parameterize best model for DOC over time
bestmodel3_sec = M2_sec

# Diagnostics
DHARMa::simulateResiduals(bestmodel3_sec, plot=T)
hist(residuals(bestmodel3_sec), breaks = 30, main = "Histogram of Residuals")
# Diagnostics look fine

# Model summaries
(sumSec_allDat_noEnv <- summary(bestmodel3_sec)) 
exp(sumSec_allDat_noEnv$coefficients$cond[, "Estimate"]) #Secchi decreases by 0.37% per year; (3.61% decrease per decade); 0.18 m decrease per decade
visreg(sumSec_allDat_noEnv, "scaled_yearSample", scale="response")
r.squaredGLMM(sumSec_allDat_noEnv)          # pretty much all lake effect

mainDF$Secchi_fit_allData_noEnv = predict( bestmodel3_sec, type = "response", re.form = NA)
mainDF$Secchi_fitWRandom_allData_noEnv = predict( bestmodel3_sec, type = "response", re.form = NULL)  # To plot just predicted lines

#Plot it
pSecchi_allData_time <- ggplot(mainDF)+
  geom_point(aes(x=yearSample, y = secchiDepth, fill=secchiDepth), alpha=0.7, size=2.5, shape = 21) +
  geom_smooth(aes(x=yearSample, y = Secchi_fit_allData_noEnv), 
              method="glm",
              formula = y~x,
              method.args = list(family = Gamma(link = 'log')))+
  scale_x_continuous(name="Year")+
  scale_y_continuous(name = "Secchi depth (m)")+
  scale_fill_gradientn(colours = c("#3D2B1F","#704214","#D2B48C","#F5DEB3")) +
  theme_DOC()+
  theme(legend.position="none")
pSecchi_allData_time

# Summarize modelled differences in DOC over time period for each trophic group
# Extract data from plot above - want the min/max of the fitted lines
plotData_secAllYear <- ggplot_build(pSecchi_allData_time)
# Extract the data related to geom_smooth
lineData_secAllYear <- plotData_secAllYear$data[[2]]  # The second layer corresponds to geom_smooth

# Group by lakeTrophicStatus and calculate the min/max of the fitted values
sumsecNoEnv_allYear <- lineData_secAllYear %>%
  summarise(minSec = min(y),
            maxSec = max(y),
            diffSec = maxSec - minSec,
            yearDiff = max(x)-min(x),
            diffPerDecade = (diffSec/yearDiff)*10)

## Change in DOC and Secchi with AHI and BsM
docTimeDat <- mainDF %>% 
  group_by(commonID) %>% 
  summarise(initDOC = DOC_fitWRandom_allData_noEnv[which.min(yearSample)],
            finalDOC = DOC_fitWRandom_allData_noEnv[which.max(yearSample)],
            diffDOC = finalDOC-initDOC,
            yearDiffDOC = max(yearSample) - min(yearSample),
            initSecchi = Secchi_fitWRandom_allData_noEnv[which.min(yearSample)],
            finalSecchi = Secchi_fitWRandom_allData_noEnv[which.max(yearSample)],
            diffSecchi = finalSecchi-initSecchi,
            sampleSize = n()) %>% 
  left_join(select(statDat_wide,commonID,lat,long)) %>% 
  distinct() %>% 
  filter(!is.na(lat)) %>% 
  mutate(dataset = "all")
docTimeDat

## Map of DOC, including estimated DOC from AHI
pDOCtimeDiff <- ggplot() +
  geom_sf(data = ontario_shapefile, fill = NA, color = "black") +
  geom_point(data=docTimeDat, aes(y=lat,x=long,fill=diffDOC,size=yearDiffDOC), 
             shape = 21, stroke = 0.5, colour = "black", alpha=0.7) +
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
             shape = 21, stroke = 0.5, colour = "black",alpha=0.7) +
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

waterClar_panelPlot <- pDOC_allData_time+pSecchi_allData_time+pDOCtimeDiff+pSecchitimeDiff + plot_layout(heights = c(1.5,2))

# ggsave(filename = ("plotsTables/plots/waterClar_panelPlot_updated_20241118.svg"),device = "svg", plot = waterClar_panelPlot,
#        width = 12, height = 12)


## Using ALL data, see if BsM enviro results hold w. Secchi ####
# Expect similar direction and magnitude of effects
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
                                                    (1|commonID), 
                                                  data=mainDF_secMod, family=Gamma(link="log"))


summary(secchiSpatial_catTrophicStatus_allData)

#plot some stuff

mainDF_secMod$predictedSecchi = predict(secchiSpatial_catTrophicStatus_allData, type = "response", re.form = NA)
mainDF_secMod$predictedSecchi_wRandom = predict(secchiSpatial_catTrophicStatus_allData, type = "response",re.form = NULL)

# Plot partial regressions of interaction between year and trophic status
newdata_allSec = mainDF_secMod
newdata_allSec$lat2 = newdata_allSec$lat
newdata_allSec$lat = mean(newdata_allSec$lat)
newdata_allSec$scaled_maxDepth = mean(newdata_allSec$scaled_maxDepth)

newdata_allSec$predictedSecchi = predict(secchiSpatial_catTrophicStatus_allData,
                                        newdata = newdata_allSec,
                                        type = "response", re.form = NA)

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
  geom_point(aes(y=secchiDepth, x=yearSample, colour=lat), alpha=0.3, size=2)+
  stat_smooth(data=newdata_allSec, aes(x=yearSample, y = predictedSecchi),
              method="glm",
              formula = y~x,
              method.args = list(family = Gamma(link = 'log')),
              se=TRUE)+
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

# ggsave(filename = ("plotsTables/plots/spatioTemp_panelPlotAll_20241118.svg"),device = "svg", plot = spatioTemp_panelPlotAll,
#        width = 8, height = 6,bg="white")

# Summarize modelled differences in Secchi over time period for each trophic group
# Extract data from plot above - want the min/max of the fitted lines
plotData_secAll <- ggplot_build(pSpatialSecAll)
# Extract the data related to geom_smooth
lineData_secAll <- plotData_secAll$data[[2]]  # The second layer corresponds to geom_smooth

# Group by lakeTrophicStatus and calculate the min/max of the fitted values
sumSec_timeAll <- lineData_secAll %>%
  group_by(PANEL) %>%  # PANEL corresponds to lakeTrophicStatus facet (1, Oligo; 2, Meso; 3, Eutro)
  summarise(minSec = min(y),
            maxSec = max(y),
            diffSec = maxSec - minSec,
            yearDiff = max(x)-min(x),
            diffPerDecade = (diffSec/yearDiff)*10)



