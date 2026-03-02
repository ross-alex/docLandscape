## Last script update. This reflects results as they're presented in the MS
# AJR: 2024-09-01

#**ALEX, note from March 10,2025 - maybe add elevation into predictive model. Easy to get info
#* otherwise, maybe add watershed slope (less easy, but useful predictor)
#* 
#* NOTE FROM June 11, 2025 - the equation and code for estimating DOC in historical AHI data using the formulas fit
#* to BsM data have been removed. Variability was too high in these data for observing an effect, even if one was seen

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
lakeElevationDat <- read_csv("data/BsM_ARU_elevation_20240320.csv")  ## Extract elevation data
statDat_wide <- read_csv("data/statDatWide_updatedCyc3_20240520.csv") %>% 
  filter(log(TKN) >4) %>%  # Clear outliers < 4 (only two observations)                       
  left_join(select(combIDs, bsmWaterbodyID, commonID), by = c("waterbodyID" = "bsmWaterbodyID"))
# Ontario shapefile: 
ontario_shapefile <- st_read("data/Province/Province.shp")
ontario_shapefile <-st_transform(ontario_shapefile, crs = 4326)
ontario_ohnWaterbody_latlon <- read_csv("data/ohnWaterbodies_csvWithCoords.csv") %>% 
  mutate(geoRegion = if_else(Latitude > 47, "North (> 47° N)", "South (< 47° N)"))

mainDF <- read_csv("data/dataForDryadRepo_landscapeDOC_AR_20241122.csv") %>% 
  left_join(select(combIDs, bsmWaterbodyID, commonID), by = "commonID") %>% 
  left_join(select(lakeElevationDat,bsmWaterbodyID,elevation), by="bsmWaterbodyID")

# Create a BsM dataset, set lake trophic status
mainDF_onlyBsM <- mainDF %>% 
  filter(samplingProgram == "BsM") %>% 
  mutate(lakeTrophicStatus = case_when(TDP < 5 ~ "oligotrophic",
                                       TDP >= 5 & TDP < 10 ~ "mesotrophic",
                                       TDP >= 10 ~ "eutrophic"),
         scaledTDP = as.vector(scale(TDP, scale = FALSE)))

mainDF_onlyBsM$lakeTrophicStatus <- factor(mainDF_onlyBsM$lakeTrophicStatus, 
                                           levels = c("oligotrophic","mesotrophic","eutrophic"))  # just changing order for plotting
# Correlation between DOC and Secchi ####
cor.test(log(mainDF_onlyBsM$secchiDepth),log(mainDF_onlyBsM$DOC))
 
pCor_secchiDOC <- ggplot(mainDF_onlyBsM) +
  geom_point(aes(x=secchiDepth,y=DOC), alpha=0.7) +    # exp y-axis to get back on raw scale
  geom_smooth(aes(x=(secchiDepth),y=(DOC)), method = "lm") +
  scale_colour_viridis_c(option="inferno", direction = -1) +
  scale_x_continuous(trans='log', labels = function(x) signif(x, 3)) +
  scale_y_continuous(trans='log', labels = function(x) signif(x, 3)) +
  ylab("DOC (mg/L)") +
  xlab("Secchi depth (m)") +
  # ylim(c(0,21))+
  # xlim(c(0,21)) +
  theme_DOC()
pCor_secchiDOC

# Correlation between DOC and TDP ####
cor.test(log(mainDF_onlyBsM$TDP),log(mainDF_onlyBsM$DOC))

# Also see what a linear model w. trophic status does ####

pCor_TDPDOC <- ggplot(mainDF_onlyBsM) +
  geom_point(aes(x=TDP,y=DOC), alpha=0.7) +    # exp y-axis to get back on raw scale
  geom_smooth(aes(x=(TDP),y=(DOC)), method = "lm") +
  scale_colour_viridis_c(option="inferno", direction = -1) +
  scale_x_continuous(trans='log', labels = function(x) signif(x, 3)) +
  scale_y_continuous(trans='log', labels = function(x) signif(x, 3)) +
  ylab("DOC (mg/L)") +
  xlab("Total dissolved phosphorus (ug/L)") +
  # ylim(c(0,21))+
  # xlim(c(0,21)) +
  theme_DOC() 
pCor_TDPDOC

pCor_TDPDOC <- ggplot(mainDF_onlyBsM) +
  geom_point(aes(x=TDP,y=DOC), alpha=0.7) +    # exp y-axis to get back on raw scale
  geom_smooth(aes(x=(TDP),y=(DOC)), method = "lm") +
  scale_colour_viridis_c(option="inferno", direction = -1) +
  scale_x_continuous(trans='log', labels = function(x) signif(x, 3)) +
  scale_y_continuous(trans='log', labels = function(x) signif(x, 3)) +
  ylab("DOC (mg/L)") +
  xlab("Total dissolved phosphorus (ug/L)") +
  # ylim(c(0,21))+
  # xlim(c(0,21)) +
  theme_DOC() +
  facet_wrap(~lakeTrophicStatus, nrow = 1, scales = "free_x")
pCor_TDPDOC



## Change in DOC and Secchi; ONLY BsM ####

## DOC
# Using IC approach, assemble different models to look for changes in DOC between 2008-2023
M1a <- glmmTMB(DOC ~ scaled_yearSample + (1|commonID),
               data = mainDF_onlyBsM,
               family = Gamma(link = "log"))
M2a <- glmmTMB(DOC ~ scaled_yearSample + (scaled_yearSample|commonID),
               data = mainDF_onlyBsM,
               family = Gamma(link = "log"))          
M3a = update(M1a, family = lognormal(link="log"))
M4a = update(M2a, family = lognormal(link="log"))
AIC(M1a,M2a,M3a,M4a)
# Gamma models are best; random intercept marginally better than random slopes. To understand how within lake 
# variance affects DOC respones across the province, it makes sense to use the random slopes and intercepts model

# Parameterize best model for DOC over time
bestmodel3a = M2a

# Diagnostics
DHARMa::simulateResiduals(bestmodel3a, plot=T)
hist(residuals(bestmodel3a), breaks = 30, main = "Histogram of Residuals")
# Diagnostics look fine

# Model summaries
(sumDOCa <- summary(bestmodel3a)) 
exp(sumDOCa$coefficients$cond[, "Estimate"]) #DOC INCREASE by 0.0021% per year; NS
r.squaredGLMM(bestmodel3a)          # pretty much all lake effect

mainDF_onlyBsM$DOC_fit = predict( bestmodel3a, type = "response", re.form = NA)
mainDF_onlyBsM$DOC_fitWRandom = predict( bestmodel3a, type = "response", re.form = NULL)  # To plot just predicted lines

## Plot it:
# Plot data points - none significant model, so no need for plotting a line
pDOCTime_bsm <- ggplot(data=mainDF_onlyBsM,aes(x=yearSample,y=DOC_fitWRandom)) +
  geom_jitter(aes(x=yearSample, y = DOC_fitWRandom, fill=DOC_fitWRandom), alpha=0.7, size=2.5, shape=21, width = 0.25) +
  #geom_smooth(method = "lm") +      # Don't plot line as relationship not significant
  scale_fill_gradientn(colours = c("#F5DEB3","#D2B48C","#704214","#3D2B1F")) +
  labs(y="DOC (mg/L)",x="Year") +
  theme_DOC()+
  theme(legend.position = "none")
pDOCTime_bsm

# ggsave(filename = ("plotsTables/plots/DOCTime_bsm_20241118.svg"),device = "svg", plot = pDOCTime_bsm,
#        width = 8, height = 6)

# Plot random slope/intercept lines
pDOCTime_bsm_randomEffects <- ggplot(data = mainDF_onlyBsM, aes(x = yearSample, y = DOC_fitWRandom, group = commonID)) +
  geom_line(alpha = 0.4) +
  labs(y = "Predicted DOC (mg/L)", x = "Year",) +
  theme_DOC() +
  theme(legend.position = "none") 
pDOCTime_bsm_randomEffects

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
exp(sumSecchiYear$coefficients$cond[, "Estimate"]) #Secchi decrease by 1.28% per year; significant
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

# Plot random slope/intercept lines
pSecchiTime_bsm_randomEffects <- ggplot(data = mainDF_onlyBsM, aes(x = yearSample, y = Secchi_fitWRandom, group = commonID)) +
  geom_line(alpha = 0.4) +
  labs(y = "Predicted Secchi (m)", x = "Year") +
  theme_DOC() +
  theme(legend.position = "none") 
pSecchiTime_bsm_randomEffects


# Summarize modelled differences in Secchi over time period for each trophic group
# Extract data from plot above - want the min/max of the fitted lines
plotData_secNoEnv <- ggplot_build(secchiBsmYearPlot)
# Extract the data related to geom_smooth
lineData_secNoEnv <- plotData_secNoEnv$data[[2]]  # The second layer corresponds to geom_smooth

# Calculate the min/max of the fitted values
sumsecNoEnv_timeBsM <- lineData_secNoEnv %>%
  summarise(minSec = min(y),
            maxSec = max(y),
            diffSec = maxSec - minSec,
            yearDiff = max(x)-min(x),
            diffPerDecade = (diffSec/yearDiff)*10)


## Change in DOC and Secchi with only BsM
docBsMTimeDat <- mainDF_onlyBsM %>% 
  group_by(commonID) %>% 
  summarise(initDOC = DOC[which.min(yearSample)],
            finalDOC = DOC[which.max(yearSample)],
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
  #geom_point(data = ontario_ohnWaterbody_latlon, 
  #           aes(x=Longitude,y=Latitude),size=1,alpha=0.1) +    #If you only want a map of all Ontario lakes, use this geom_point (disregards DOC viz)
  geom_point(data=docBsMTimeDat, aes(y=lat,x=long,fill=diffDOC,size=yearDiffDOC), 
              shape = 21, stroke = 0.5, colour = "black", alpha=0.7) +
  scale_fill_gradient2(
     low = "#aa95cd",           
     mid = "white",         
     high = "brown",     
     midpoint = 0)+  
  ylim(c(42.5,57)) +
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
  theme(legend.position = "right")  #This is just to add legend to plot (custom)
# ggsave(filename = ("plotsTables/plots/secchi_bsmMap_legTop_20240901.svg"),device = "svg", plot = pSecchiTimeDiff_BsM_legTop,
#        width = 8, height = 6)

waterClarBsM_panelPlot <- pDOCTime_bsm+secchiBsmYearPlot+pDOCtimeDiff_BsM+pSecchitimeDiff_BsM + plot_layout(heights = c(1.5,2))

# Save as svg
# ggsave(filename = ("plotsTables/plots/waterClarBsM_panelPlot_20241118.svg"),device = "svg", plot = waterClarBsM_panelPlot,
#        width = 12, height = 12)

## Change in DOC and Secchi with environmental variables; ONLY BsM ####
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
                                          data=mainDF_onlyBsM, family=Gamma(link="log"), na.action = "na.fail")  #This model won't fit. Too many variables for predicting random 
                                                                                                                # slopes (some cases only 2 data points; ok with only 
                                                                                                                # one predictor, not ok with all these)
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

bestmodel4=glmmTMB(DOC ~ scaled_yearSample+
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

# Need 95%CIs to plot; function to predict from bootstrapped models
bootFun <- function(model) {
  predict(model, newdata = newdata, type = "response", re.form = NA)
}

# Perform bootstrapping (adjust nsim for precision vs. speed)
bootResults <- bootMer(bestmodel4, FUN = bootFun, nsim = 1000, re.form = NA)

# Compute confidence intervals
newdata$predictedDOC_wEnvVars <- predict(bestmodel4, newdata = newdata, type = "response", re.form = NA)
newdata$lwr <- apply(bootResults$t, 2, quantile, probs = 0.025)  # 2.5% CI
newdata$upr <- apply(bootResults$t, 2, quantile, probs = 0.975)  # 97.5% CI

pSpatialDOCBsM <- ggplot(newdata)+
  geom_jitter(aes(y=DOC, x=yearSample, colour=lat2), alpha=0.3, size=2, width=0.22)+
  geom_smooth(aes(x=yearSample, y = predictedDOC_wEnvVars),
              method="glm",
              formula = y~x,
              se = TRUE,
              method.args = list(family = Gamma(link = 'log')))+
  geom_ribbon(aes(x = yearSample, ymin = lwr, ymax = upr), fill = "#C6DBEF", alpha = 0.2) +  # Confidence interval
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

# Need 95%CIs to plot; function to predict from bootstrapped models
bootFun_sec <- function(model) {
  predict(model, newdata = newdata_sec, type = "response", re.form = NA)
}

# Perform bootstrapping (adjust nsim for precision vs. speed)
bootResults_sec <- bootMer(secchiSpatial_catTrophicStatus, FUN = bootFun_sec, nsim = 10, re.form = NA)

# Compute confidence intervals
newdata_sec$predictedSecchi_wEnvVars <- predict(secchiSpatial_catTrophicStatus, newdata = newdata_sec, type = "response", re.form = NA)
newdata_sec$lwr <- apply(bootResults_sec$t, 2, quantile, probs = 0.025)  # 2.5% CI
newdata_sec$upr <- apply(bootResults_sec$t, 2, quantile, probs = 0.975)  # 97.5% CI


pSpatialSecBsM <- ggplot(newdata_sec)+
  geom_jitter(data = mainDF_onlyBsM, aes(y=secchiDepth, x=yearSample, colour=lat), alpha=0.3, size=2, width=0.22)+
  geom_smooth(aes(x=yearSample, y = predictedSecchi_wEnvVars),
              method="glm",
              formula = y~x,
              method.args = list(family = Gamma(link = 'log')))+
  geom_ribbon(aes(x = yearSample, ymin = lwr, ymax = upr), fill = "#C6DBEF", alpha = 0.2) +  # Confidence interval
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
  theme(legend.position = "right")
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
  theme(legend.position = "right")
pSecchitimeDiff

waterClar_panelPlot <- pDOC_allData_time+pSecchi_allData_time+pDOCtimeDiff+pSecchitimeDiff + plot_layout(heights = c(1.5,2))

# ggsave(filename = ("plotsTables/plots/waterClar_panelPlot_updated_20241118.svg"),device = "svg", plot = waterClar_panelPlot,
#        width = 12, height = 12)

## Change in DOC using only oligotrophic lakes; BsM and AHI ####
# Using IC approach, assemble different models to look for changes in DOC over time
mainDF_AHI <- mainDF %>% 
  filter(samplingProgram == "ARU")

mainDF_oligo_BsM <- mainDF %>% 
  mutate(lakeTrophicStatus = case_when(TDP < 5 ~ "oligotrophic",
                                       TDP >= 5 & TDP < 10 ~ "mesotrophic",
                                       TDP >= 10 ~ "eutrophic"),
         scaledTDP = as.vector(scale(TDP, scale = FALSE))) %>% 
  filter(lakeTrophicStatus == "oligotrophic") 

mainDF_oligo_all <- mainDF_oligo %>% 
  bind_rows(filter(mainDF_AHI, commonID %in% unique(mainDF_oligo_BsM$commonID)))
  

M1 <- glmmTMB(updatedDOC ~ scaled_yearSample + (1|commonID),
              data = mainDF_oligo_all,
              family = Gamma(link = "log"))
M2 <- glmmTMB(updatedDOC ~ scaled_yearSample + (scaled_yearSample|commonID),
              data = mainDF_oligo_all,
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

mainDF_oligo_all$DOC_fit_allData_noEnv = predict( bestmodel3, type = "response", re.form = NA)
mainDF_oligo_all$DOC_fitWRandom_allData_noEnv = predict( bestmodel3, type = "response", re.form = NULL)  # To plot just predicted lines

# Plot it
pDOC_allData_time <- ggplot(mainDF_oligo_all)+
  geom_point(aes(x=yearSample, y = updatedDOC, fill=updatedDOC), alpha=0.7, size=2.5, shape = 21) +
  geom_smooth(aes(x=yearSample, y = DOC_fit_allData_noEnv), 
              method="glm",
              formula = y~x,
              method.args = list(family = Gamma(link = 'log')))+
  scale_x_continuous(name="Year")+
  scale_y_continuous(name = "DOC (mg/L)")+
  #scale_fill_gradientn(colours = c("#F5DEB3","#D2B48C","#704214","#3D2B1F")) +
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
  theme(legend.position = "right")
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
  theme(legend.position = "right")
pSecchitimeDiff

waterClar_panelPlot <- pDOC_allData_time+pSecchi_allData_time+pDOCtimeDiff+pSecchitimeDiff + plot_layout(heights = c(1.5,2))

# ggsave(filename = ("plotsTables/plots/waterClar_panelPlot_updated_20241118.svg"),device = "svg", plot = waterClar_panelPlot,
#        width = 12, height = 12)



## Bootstrap Secchi results - 15 yr slope legit or not? ####

## Bootstrap years to see if observed slope of Secchi depth from 15yr time series is stronger
# then by random chance (i.e., if years were shuffled). Randomly draw 15 year segments (i.e.,
# different 'start' times) and create a distribution of slope values

## Filter for time ranges
# 15-yr BsM data range
contempBsM_years <- sort(unique(mainDF_onlyBsM$yearSample))
contempBsM_period_length <- length(unique(contempBsM_years))
# BsM and AHI data range
all_years <- sort(unique(mainDF$yearSample))
all_period_length <- length(unique(all_years))

n_boot <- 999

## Go back to Mikes Lecture code -- figure out how to shuffle years **** ##

# Step 1: Fit model to the actual present period
present_mod <- glmmTMB(secchiDepth ~ scaled_yearSample + (scaled_yearSample | commonID),
                       data = mainDF_onlyBsM,
                       family = Gamma(link = "log"))

observed_slope <- summary(present_mod)$coefficients$cond["scaled_yearSample", "Estimate"]

# Step 2: Bootstrapping: shuffle years
boot_slopes <- numeric(n_boot)

set.seed(42)  # reproducibility
for (i in 1:n_boot) {
  # Randomly choose a starting year such that you can get a full 15-year period
  possible_start_years <- 1962:1972
  #possible_start_years <- all_years[all_years <= max(all_years) - 14]
  start_year <- sample(possible_start_years, 1)
  selected_years <- start_year:(start_year + 14)
  
  # Filter your data to the selected 15-year period
  boot_data <- mainDF %>%
    filter(yearSample %in% selected_years) %>%
    mutate(realYear = yearSample,
           yearSample_fake = scale(yearSample))  # scale the actual years
  
  # Run bootstrap model - has to just be fixed-effect model; not enough data where lakeIDs > 1 because AHI data only had 1 lake visit
  mod <- glmmTMB(secchiDepth ~ yearSample_fake,
                 data = boot_data,
                 family = Gamma(link = "log"))
  
  # collect bootstrap slopes
  boot_slopes[i] <- summary(mod)$coefficients$cond["yearSample_fake", "Estimate"]
}

boot_slopes <- na.omit(boot_slopes)
# Step 3: Plot results
df_boot <- data.frame(slope = boot_slopes)
ggplot(df_boot, aes(x = slope)) +
  geom_histogram(binwidth = 0.01, fill = "skyblue", color = "black") +
  geom_vline(xintercept = observed_slope, color = "red", size = 1.2, linetype = "dashed") +
  labs(title = "Bootstrapped Slopes vs Observed Slope",
       x = "Slope (scaled yearSample)", y = "Frequency") +
  theme_minimal()

# Optional: empirical p-value
empirical_p <- mean(abs(boot_slopes) >= abs(observed_slope))
cat("Empirical p-value:", empirical_p, "\n")


# Step 2b: Bootstrapping: shuffle years randomly (not 15 consecutive years)
boot_slopes <- numeric(n_boot)

set.seed(42)  # reproducibility
for (i in 1:n_boot) {
  # Randomly choose a starting year such that you can get a full 15-year period
  #possible_start_years <- 1962:1972
  possible_start_years <- all_years[all_years <= max(all_years) - 16]
  #selected_years <- sample(possible_start_years,15)  # use this one if don't want any BsM
  selected_years <- sample(all_years,15)   # use this one if want to use all years (including BsM)
  
  # Filter your data to the selected 15-year period
  boot_data <- mainDF %>%
    filter(yearSample %in% selected_years) %>%
    mutate(realYear = yearSample,
           yearSample_fake = scale(yearSample))  # scale the actual years
  
  # Run bootstrap model - has to just be fixed-effect model; not enough data where lakeIDs > 1 because AHI data only had 1 lake visit
  mod <- glmmTMB(secchiDepth ~ yearSample_fake,
                 data = boot_data,
                 family = Gamma(link = "log"))
  
  # collect bootstrap slopes
  boot_slopes[i] <- summary(mod)$coefficients$cond["yearSample_fake", "Estimate"]
}

boot_slopes <- na.omit(boot_slopes)
# Step 3: Plot results
df_boot <- data.frame(slope = boot_slopes)
ggplot(df_boot, aes(x = slope)) +
  geom_histogram(binwidth = 0.01, fill = "skyblue", color = "black") +
  geom_vline(xintercept = observed_slope, color = "red", size = 1.2, linetype = "dashed") +
  labs(title = "Bootstrapped Slopes vs Observed Slope",
       x = "Slope (scaled yearSample)", y = "Frequency") +
  theme_minimal()

# Optional: empirical p-value
empirical_p <- mean(abs(boot_slopes) >= abs(observed_slope))
cat("Empirical p-value:", empirical_p, "\n")

## Testing observed secchi slope; different than chance? ####
# Try adding a "period" term into the model, and test for an interaction
# If significant interaction, we then know something about period differences, if no interaction
# then it's a general trend over time

# Step 1: Fit model to the actual present period
present_mod <- glmmTMB(secchiDepth ~ scaled_yearSample + (scaled_yearSample | commonID),
                       data = mainDF_onlyBsM,
                       family = Gamma(link = "log"))

observed_slope <- summary(present_mod)$coefficients$cond["scaled_yearSample", "Estimate"]

# Step 2: Use all data, include an interaction term
interactionSec_mod <- glmmTMB(secchiDepth ~ scaled_yearSample*samplingProgram + (scaled_yearSample | commonID),
                              data = mainDF,
                              family = Gamma(link = "log"))

summary(interactionSec_mod) ## Significant interaction, time points are different

# Step 3; Compare above model to that of common slope
commonslope_Secmod <- glmmTMB(secchiDepth ~ scaled_yearSample + (scaled_yearSample | commonID),
                            data = mainDF,
                            family = Gamma(link = "log"))

# Predict this common slope out
min_scaled_year <- min(mainDF$scaled_yearSample)
max_scaled_year <- max(mainDF$scaled_yearSample)
predDat_commonSlope <- data.frame(scaled_yearSample = seq(min_scaled_year, max_scaled_year, length.out = 100))

predDat_commonSlope$predicted_secchiDepth <- predict(commonslope_Secmod,
                                                  newdata = predDat_commonSlope,
                                                  re.form = NA, # Ignore random effects for common slope
                                                  type = "response")



# Step 4; plot it all
ggplot(mainDF, aes(x = scaled_yearSample, y= secchiDepth, colour = samplingProgram)) +
  geom_point(aes(fill = samplingProgram), shape = 21, colour = "black", stroke = 0.5, alpha = 0.7) +
  geom_smooth(aes(x=scaled_yearSample, y = secchiDepth), 
              method="glm",
              formula = y~x,
              method.args = list(family = Gamma(link = 'log'))) +
  scale_x_continuous(breaks = c(-40,-20,0, 20), labels = c(1960,1980,2000,2020)) +
  # Add the common slope line
  geom_line(data = predDat_commonSlope,
            aes(x = scaled_yearSample, y = predicted_secchiDepth),
            inherit.aes = FALSE, # Important to avoid inheriting 'colour' mapping
            linetype = "dashed", # Optional: make it dashed to distinguish
            color = "black",     # Optional: make it black for commonality
            size = 1) +          # Optional: adjust line thickness
  scale_fill_manual(values = c("#56B4E9", "#E69F00")) +
  scale_colour_manual(values = c("#56B4E9", "#E69F00")) +
  labs(x= "Year Sample", y = "Secchi depth (m)") +
  guides(fill = "none", colour = "none") +
  theme(legend.position = "none") +
  theme_DOC()


## Using ALL data, see if BsM enviro results hold w. Secchi ####
# Expect similar direction and magnitude of effects
mainDF_secMod <- mainDF %>%
  group_by(commonID) %>%
  mutate(lakeTrophicStatus = case_when(
    TDP < 5 ~ "oligotrophic",
    TDP >= 5 & TDP < 10 ~ "mesotrophic",
    TDP >= 10 ~ "eutrophic")) %>%
  fill(lakeTrophicStatus, .direction = "downup")  # Fill NAs for ARU observations within each commonID

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

# ggsave(filename = ("plotsTables/plots/spatioTemp_panelPlotAll_20241122.svg"),device = "svg", plot = spatioTemp_panelPlotAll,
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



## Histogram of lake size and trophic statuses
histSize <- ggplot(ontario_ohnWaterbody_latlon) +
  #geom_histogram(aes(x=SYSTEM_CALCULATED_AREA))
  geom_histogram(aes(x=(SYSTEM_CALCULATED_AREA/10000))) +
  scale_x_log10(labels = scales::label_number()) +  # Removes scientific notation  theme_minimal()
  xlab("Surface area (ha)") +
  ylab("Frequency") +
  theme_DOC() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  facet_wrap(~geoRegion)
histSize

histSize <- ggplot(ontario_ohnWaterbody_latlon) +
  #geom_histogram(aes(x=SYSTEM_CALCULATED_AREA))
  geom_histogram(aes(x=(SYSTEM_CALCULATED_AREA/10000))) +
  scale_x_log10(labels = scales::label_number()) +  # Removes scientific notation  theme_minimal()
  xlab("Surface area (ha)") +
  ylab("Frequency") +
  theme_DOC() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  facet_wrap(~geoRegion)
histSize

mainDF_onlyBsM_wSA <- mainDF_onlyBsM %>% 
  left_join(select(statDat_wide,commonID,surfaceArea)) %>% 
  distinct()

histTrophicStatus_loc <- ggplot(mainDF_onlyBsM_wSA) +
  geom_histogram(aes(x=surfaceArea)) +
  facet_wrap(~lakeTrophicStatus) +
  scale_x_log10(labels = scales::label_number()) +  # Removes scientific notation  theme_minimal()
  theme_DOC() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) 
histTrophicStatus_loc

